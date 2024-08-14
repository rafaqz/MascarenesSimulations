
namedkeys(nt::NamedTuple{K}) where K = NamedTuple{K}(K)

function agg_aux(masks::NamedTuple, dems, lcs, aggfactor, last_year)
    map(masks, dems, lcs) do args...
        agg_aux(args..., aggfactor, last_year)
    end
end
function agg_aux(mask_orig::AbstractArray, dem_orig, lc_orig, aggfactor, last_year)
    mask = Rasters.aggregate(Rasters.Center(), mask_orig, aggfactor; skipmissingval=true)
    # Clip roughness at a maximum of 500
    # rgh = Float32.(replace_missing(Rasters.aggregate(maximum, min.(roughness(replace_missing(dem_orig, 0)), 500)./ 500, aggfactor), NaN))
    dem = Float32.(replace_missing(Rasters.aggregate(Rasters.Center(), dem_orig, aggfactor; skipmissingval=true), 0)) .* mask
    lc_ag = Rasters.aggregate(mean, rebuild(lc_orig; missingval=nothing), (X(aggfactor), Y(aggfactor)); skipmissingval=true)
    lc_ag1 = Rasters.extend(lc_ag; to=(Ti(Sampled(1500:1:2020; sampling=Intervals(Start())))), missingval=false)
    map(layers(lc_ag1)) do A
        broadcast_dims!(identity, view(A, Ti=1500..1600), view(A, Ti=At(1600)))
    end
    lc = map(layers(lc_ag1)...) do xs...
        Float32.(NV{keys(lc_orig),length(keys(lc_orig))}(xs))
    end .* mask
    (; mask, dem, lc)
end

function gpu_cleanup(A)
    x, y, ti = lookup(A, (X, Y, Ti))
    Rasters.set(A,
        Ti => Sampled(first(ti) - 0.5:last(ti) - 0.5; sampling=Intervals(Start())),
        X => LinRange(first(x), last(x), length(x)),
        Y => LinRange(first(y), last(y), length(y)),
    )
end

function generate_predator_effect!(f, x, pred_pop, pred_suscept)
    ThreadsX.map!(x, pred_pop) do pop
        map(f, predator_effect(pop, pred_suscept))
    end end

function generate_predator_effect(f, pred_pop::Union{AbstractArray{<:Any,2},AbstractArray{<:Any,3}}, pred_suscept)
    xs = ThreadsX.map(pred_pop) do pop
        map(f, predator_effect(pop, pred_suscept))
    end
    rebuild(pred_pop, xs)
end

function predict_timeline(endemic_ruleset::Ruleset, islands, pred_response; kw...)
    map(islands) do island
        _predict_timeline(endemic_ruleset, island, pred_response; kw...)
    end
end
function _predict_timeline(endemic_ruleset::Ruleset, island, pred_response; range_times, kw...)
    output_kw = island.output_kw
    aux = output_kw.aux
    pred_suscept = predator_suceptibility(pred_response, aux.endemic_traits)
    (; pred_pop, pred_effect) = aux
    generate_predator_effect!(tanh, pred_effect, pred_pop, pred_suscept)
    output = TransformedOutput(island.endemic_init; output_kw...) do f, (i, t)
        # Take the sum of each replicates slice
        replicates = eachslice(f.endemic_presence; dims=3)
        presence_sums = parent(ThreadsX.map(sum, replicates))
        ranges = if t in range_times
            copy(f.endemic_presence)
        else
            missing
        end
        (; presence_sums, ranges)
    end
    sim!(output, endemic_ruleset; kw..., proc=CPUGPU())
    return output
end

scale_params(x, cc) = Float32.(x ./ cc .* minimum(cc))

function extinction_objective(x, p; kw...)
    (; last_year, extant_extension, endemic_ruleset, islands, parameters, loss, range_times, pred_carrycap, island_ranges, obs) = p
    (; plot_obs, params_obs, loss_obs, throwit) = obs
    throwit[] && error()
    # Update parameter values scaled for predator carrycap
    parameters[:val] = scale_params(x, pred_carrycap)
    pred_response = stripparams(parameters)
    # Run the simulations
    island_timelines = predict_timeline(endemic_ruleset, islands, pred_response; range_times, kw...)

    # Update plot
    params_obs[] .= x
    notify(params_obs)
    # Calculate loss
    island_losses = map(island_timelines, island_ranges, islands, plot_obs) do timeline, ranges, island, obs
        dates = extinction_dates_from_sim(timeline, island, obs; last_year, extant_extension)
        range_loss = sum_extinction_ranges(timeline, ranges, obs) |> sum
        date_loss = map(loss, island.extinction_dates, dates.mean) |> sum
        (; date_loss, range_loss)
    end
    loss = sum(island_losses) do l
        @show l
        # Just add losses. its not perfect but its stable
        # And these objectives are highly correlated, not opposing
        l.date_loss + 2l.range_loss
    end
    @show loss
    loss_obs[] = "loss: $loss"
    notify(loss_obs)

    return loss
end

function sum_extinction_ranges(timeline, ranges, obs=nothing)
    mapreduce(+, lookup(ranges, Ti)) do t
        known = ranges[At(t)].endemic_presence
        predicted = timeline[At(t)].ranges
        replicates = eachslice(predicted; dims=3)
        if obs isa NamedTuple
            div = obs.diversity_obs
            div[] .= mean(sum, predicted; dims=3)
            div[] .= (x -> x == 0.0 ? NaN : x).(div[])
            notify(div)
        end
        # Count the number of matching presence/absence in each replicate
        ntotal = length(islands.mus.aux.mask)
        x = ThreadsX.map(replicates) do predicted
            mapreduce(.!=, +, predicted, known)
        end |> mean
    end |> mean
end

function extinction_dates_from_sim(timeline, island, obs=nothing;
    last_year, extant_extension,
)
    firstyear = first(island.output_kw.tspan)
    years_present = sum(timeline) do (slice, _)
        map(s -> s .> 0, slice)
    end
    extinction_years = map(years_present) do yp
        map(yp) do y
            y1 = y + firstyear - 1
            # Handle not-going-extinct
            y1 >= last_year ? y1 + extant_extension : y1
        end
    end
    ext_mean = mean(extinction_years)
    ext_std = std(extinction_years)
    if obs isa NamedTuple
        dates = obs.dates_obs
        dates[] .= ext_mean
        notify(dates)
    end
    return (years=extinction_years, mean=ext_mean, std=ext_std)
end

function extinction_dates_from_tables(tables, EndemicNVs, island_keys; last_year, extant_extension)
    # Extract extinction dates
    map(tables, EndemicNVs, island_keys) do table, EndemicNV, key
        map(table[!, Symbol(key, :_extinct)]) do x
            if ismissing(x)
                # use the last year of the simulation + extan_extinction as the "not extinct yet" extinction date
                last_year + extant_extension
            else
                parse(Int, first(split(x, ':')))
            end
        end |> EndemicNV
    end
end

function extinction_forward(x, p; kw...)
    (; endemic_ruleset, extant_extension, islands, parameters, last_year, loss, range_times, pred_carrycap) = p
    # Update
    parameters[:val] = scale_params(x, pred_carrycap)
    pred_response = stripparams(parameters)
    @show pred_response
    timelines = predict_timeline(endemic_ruleset, islands, pred_response; range_times, kw...)
    dates = map(timelines, islands) do t, i
        extinction_dates_from_sim(t, i; last_year, extant_extension)
    end
    return map((dates, timelines) -> (; dates, timelines), dates, timelines)
end


function endemic_traits(tables::NamedTuple, NVs)
    map(tables, NVs) do table, NV
        endemic_traits(table, NV)
    end
end
function endemic_traits(table, NV)
    # hunting_preference = NV(table.Hunting_preference)
    ismammal = NV(table.Group .== "mammal")
    isbird = NV(table.Group .== "bird")
    isreptile = NV(table.Group .== "reptile")
    isgroundnesting = NV(Float32.(table.Ground_nesting))
    flightlessness = NV(Float32.(table.Flightlessness))
    return (; ismammal, isbird, isreptile, isgroundnesting, flightlessness)#, hunting_preference)
end


function predator_response_params(pred_keys)
    # These dummy values will be optimised. Its likely too many for that but we
    # can reduce the numver of predators, removing mouse, macaque etc
    pred_response_raw = (;
        ismammal = (;
            cat =         0.1,
            black_rat =   0.1,
            norway_rat =  0.1,
            mouse =       0.1,
            pig =         0.1,
            wolf_snake =  0.1,
            macaque =     0.1,
        ),
        isbird = (;
            cat =         0.3,
            black_rat =   0.3,
            norway_rat =  0.3,
            mouse =       0.1,
            pig =         0.1,
            wolf_snake =  0.1,
            macaque =     0.1,
        ),
        isreptile = (;
            cat =         0.3,
            black_rat =   0.3,
            norway_rat =  0.3,
            mouse =       0.1,
            pig =         0.1,
            wolf_snake =  0.1,
            macaque =     0.1,
        ),
        # isgroundnesting = (;
        #     cat =         0.2,
        #     black_rat =   0.01,
        #     norway_rat =  0.02,
        #     mouse =       0.01,
        #     pig =         0.05,
        #     wolf_snake =  0.0,
        #     macaque =     0.0,
        # ),
        flightlessness = (;
            cat =         0.4,
            black_rat =   0.2,
            norway_rat =  0.4,
            mouse =       0.1,
            pig =         0.1,
            wolf_snake =  0.1,
            macaque =     0.1,
        ),
    )

    response_keys = NamedTuple{keys(pred_response_raw)}(keys(pred_response_raw))
    pred_response = map(response_keys, pred_response_raw) do k, ps
        map(NamedTuple(ps[pred_keys]), NamedTuple{pred_keys}(pred_keys)) do val, label
            Param(val; label=Symbol(k, :_, label), bounds=(0, 1))
        end
    end
end
