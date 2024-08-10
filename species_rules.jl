using DynamicGrids
using Dispersal
using LandscapeChange
using ModelParameters
using StaticArrays
using Distributions
using Setfield

const DG = DynamicGrids
const DD = DimensionalData
const NV = NamedVector

include("functions.jl")

# Struct rules

struct InteractiveCarryCap{R,W,CC,CS,I} <: Dispersal.GrowthRule{R,W}
    carrycap::CC
    carrycap_scaling::CS
    inputs::I
end
function InteractiveCarryCap{R,W}(; carrycap, carrycap_scaling, inputs=(;),) where {R,W}
    InteractiveCarryCap{R,W}(carrycap, carrycap_scaling, inputs)
end

@inline Base.@assume_effects :foldable function DynamicGrids.applyrule(
    data, rule::InteractiveCarryCap,
    populations::NamedVector{Keys}, I,
) where Keys
    local_inputs = get(data, rule.inputs, I)
    return calc_carrycaps(local_inputs, populations, rule.carrycap, rule.carrycap_scaling)
end

function calc_carrycaps(local_inputs, populations, carrycap, carrycap_scaling)
    relative_pop = NamedTuple(populations ./ carrycap) # combine populations
    params = merge(relative_pop, NamedTuple(local_inputs))
    scaling = map(carrycap_scaling) do val_f
        f = DynamicGrids._unwrap(val_f)
        (oneunit(eltype(populations)) + f(params))
    end |> NamedVector 
    # Carrycap cant be zero for numerical reasons, but make the minimum very small
    absolute_min_carrycap = oneunit(eltype(carrycap)) .* 1f-10
    new_carrycaps = max.(absolute_min_carrycap, carrycap .* scaling)

    return new_carrycaps
end

struct ExtirpationRisks{R,W,F,S,T,PR,PS,PP,PC,E,SE,RR} <: DynamicGrids.NeighborhoodRule{R,W}
    f::F
    stencil::S
    traits::T
    pred_response::PR
    pred_suscept::PS
    pred_pop::PP
    pred_carrycap::PC
    pred_effect::E
    stochastic_extirpation::SE
    recouperation_rates::RR
end
function ExtirpationRisks{R,W}(; f, stencil, traits, pred_response, pred_suscept, pred_pop, pred_carrycap, pred_effect, stochastic_extirpation, recouperation_rates) where {R,W}
    ExtirpationRisks{R,W}(f, stencil, traits, pred_response, pred_suscept, pred_pop, pred_carrycap, pred_effect, stochastic_extirpation, recouperation_rates)
end

# endemic_recouperation_rule = let recouperation_rate_aux=Aux{:recouperation_rate}()
#     Neighbors{:endemic_presence}(Moore(1)) do data, hood, presences, I
#         # any(presences) || return presences
#         recouperation_rate = DG.get(data, recouperation_rate_aux)
#         nbr_sums = foldl(hood; init=Base.reinterpret.(UInt8, zero(first(hood)))) do x, y
#             Base.reinterpret(UInt8, x) + Base.reinterpret(UInt8, y)
#         end
#         map(presences, nbr_sums, recouperation_rate) do p, n_nbrs, rr
#             if p
#                 true
#             elseif n_nbrs > 0
#                 rand(Float32) < (n_nbrs * rr / length(hood))
#             else
#                 false
#             end
#         end
#     end
# end

@inline function DynamicGrids.applyrule(data, rule::ExtirpationRisks, endemic_presences, I)
    recouperation_rates = get(data, rule.recouperation_rates)

    # Get the effect of predators on each endemic
    pred_effect = if isnothing(rule.pred_effect)
        pred_relative_pop = get(data, rule.pred_pop, I) ./ rule.pred_carrycap
        pred_suscept = get(data, rule.pred_suscept)
        map(rule.f, predator_effect(pred_relative_pop, pred_suscept))
    else
        # We have already precalculated the effect
        get(data, rule.pred_effect, I)
    end

    hood = stencil(rule)
    # Count neighbors to UInt8. Otherwise we get Int64 and use a lot of registers
    n_neighbors = foldl(hood; init=Base.reinterpret.(UInt8, zero(first(hood)))) do x, y
        Base.reinterpret(UInt8, x) + Base.reinterpret(UInt8, y)
    end

    updated_presence = map(endemic_presences, pred_effect, n_neighbors, recouperation_rates) do present, effect, n_neighbor, rr
        if present
            rand(typeof(effect)) > (effect * ((length(hood) / 2) / (n_neighbor + 1)) + rule.stochastic_extirpation)
        elseif n > 0
            rand(typeof(effect)) < (n_neighbor * rr / length(hood))
        else
            false
        end
    end
    return updated_presence
end

@inline function DynamicGrids.applyrule(data, rule::ExtirpationRisks{Grids,Grids}, (endemic_presences, causes), I) where Grids<:Tuple
    recouperation_rates = get(data, rule.recouperation_rates)

    hood = stencil(rule)
    # Count neighbors to UInt8. Otherwise we get Int64 and use a lot of registers
    n_neighbors = foldl(hood; init=Base.reinterpret.(UInt8, zero(first(hood)))) do x, y
        Base.reinterpret(UInt8, x) + Base.reinterpret(UInt8, y)
    end

    # Get the effect of predators on each endemic
    pred_relative_pop = get(data, rule.pred_pop, I) ./ rule.pred_carrycap .* 32
    pred_suscept = get(data, rule.pred_suscept)
    pred_effects = predator_effects(pred_relative_pop, pred_suscept)

    results = map(endemic_presences, pred_effects, n_neighbors, recouperation_rates, causes) do present, effects, n_neighbor, rr, cs
        if present
            effect = rule.f(sum(effects))
            extirpation_probability = effect * length(hood) / (n_neighbor + 1) + rule.stochastic_extirpation
            still_present = rand(typeof(first(pred_relative_pop))) > extirpation_probability
            updated_causes = if still_present 
                cs
            else
                effects
            end
        elseif n_neighbor > 0
            still_present = rand(typeof(first(pred_relative_pop))) < (n_neighbor * rr / length(hood))
            updated_causes = if still_present
                zero(cs)
            else
                cs
            end
        else
            still_present = false
            updated_causes = cs
        end
        still_present, updated_causes
    end
    updated_presence = map(first, results)
    updated_causes = map(last, results)

    return updated_presence, updated_causes
end

@inline function DynamicGrids.modifyrule(rule::ExtirpationRisks, data::AbstractSimData)
    if isnothing(rule.pred_effect)
        pred_response = get(data, rule.pred_response)
        # mass_response = get(data, rule.mass_response)
        traits = get(data, rule.traits)
        @set! rule.pred_suscept = predator_suceptibility(pred_response, traits)
    end
    @set! rule.traits = nothing # Simplify for GPU argument size
    @set rule.pred_response = nothing # Simplify for GPU argument size
end

Base.@assume_effects :foldable function predator_effects(pred_pop, pred_suscept::NamedVector{<:Any,<:Any,<:NamedVector{K}}) where K
    is = ntuple(identity, length(K))
    map(is) do i
        map(pred_suscept, pred_pop) do ps, pp
            Float32(ps[i] * pp)
        end
    end |> NamedVector{K,length(K)}
end
Base.@assume_effects :foldable function predator_effect(pred_pop, pred_suscept)
    mapreduce(+, pred_suscept, pred_pop) do ps, pp
         map(xs -> map(Float32 ∘ *, xs, pp), ps)
    end
end

# function predator_suceptibility(mass_response, pred_response, traits)
#     pred_suscept = mapreduce(+, pred_response, traits) do pr, t
#         map(mass_response, pr) do m, p
#             t .* p .* m
#         end
#     end ./ (32 * 8^2)
# end

function predator_suceptibility(pred_response, traits)
    mapreduce(+, pred_response, traits) do pr, t
        map(NamedVector(pr)) do p
            NamedVector(t) .* p
        end |> NamedVector
    end ./ (32 * 8^2)
end


# Define shorthand rules, inits and aux data

function def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    aggscale = aggfactor^2,
    replicates=nothing,
    pred_pops_aux,
    island_keys = NamedTuple{keys(island_endemic_tables)}(keys(island_endemic_tables)),
    first_year, last_year,
    extant_extension,
    pred_keys=(:cat, :black_rat, :norway_rat, :mouse, :pig, :wolf_snake, :macaque),
    pred_response=predator_response_params(pred_keys),
    EndemicNVs = begin
        map(island_endemic_tables) do endemic_table
            ek = Tuple(Symbol.(replace.(endemic_table.Species, Ref(' ' => '_'))))
            NamedVector{ek,length(ek)}
        end
    end,
    island_extinction_dates = extinction_dates_from_tables(island_endemic_tables, EndemicNVs, island_keys, extant_extension),
    mean_prey_mass = (;
        cat =        (41.0, 51.0), # Pearre and Maaas 1998
        black_rat =  (10.0, 10.0), # made up
        norway_rat = (8.0f0, 8.0), # made up
        mouse =      (3.0f0, 5.0), # made up
        pig =        (100.0, 100.0), # made up
        wolf_snake = (9.0f0, 7.0), # Estimated from Fritts 1993
        macaque =    (40.0, 40.0), # made up
   )[pred_keys],
    island_mass_response = map(island_endemic_tables, EndemicNVs) do table, EndemicNV
        endemic_mass = EndemicNV(table.Mass)
        map(mean_prey_mass) do (mean, std)
            dist = Distributions.Normal(mean, 2std) # Doubled because of the skew
            scalar = 1 / pdf(dist, mean)
            map(endemic_mass) do pm
                pdf(dist, pm) * scalar
            end
        end
    end,
    # TODO: This is made up
    island_recouperation_rates = map(EndemicNVs) do EndemicNV
        Float32.(ones(EndemicNV) .* 1.0)
    end,
    # These are taken from the literature in contexts where it seems also applicable to the Mascarenes
    carrycap = Float32.(NV(;
        cat =        0.01,
        black_rat =  30.0, # This may be up to 100/ha? See Harper & Bunbury 2015.
        norway_rat = 15.0, # This one is more of a guess
        mouse =      52.0,
        pig =        0.02, # Tuned to have ~600 total in mauritius in 2009 (actually 714, more significant digits are not justifiable).
        wolf_snake = 10.0, # As is this one
        macaque =    0.5,
    ) .* aggscale)[pred_keys],
    spread_rate = NV(;
        cat =         20.0f0,
        black_rat =   1.0f0,
        norway_rat =  1.0f0,
        mouse =       0.5f0,
        pig =         100.0f0,
        wolf_snake =  1.0f0,
        macaque =     5.0f0,
    )[pred_keys],
    # pred_funcs = (;
    #     cat =        p -> 1.0f0p.black_rat + 0.3f0p.norway_rat + 1.0f0p.mouse + 10f0p.urban + 2f0p.cleared,
    #     black_rat  = p -> -0.2f0p.cat - 0.1f0p.norway_rat - 0.1f0p.mouse + 0.5f0p.native + 0.3f0p.abandoned + 0.3f0p.forestry + 1p.urban,
    #     norway_rat = p -> -0.1f0p.cat - 0.1f0p.black_rat - 0.1f0p.mouse + 1.5f0p.urban - 0.2f0p.native,
    #     mouse =      p -> -0.3f0p.cat - 0.2f0p.black_rat - 0.2f0p.norway_rat + 0.8f0p.cleared + 1.5f0p.urban,
    #     pig =        p -> 0.0f0p.native - 0.0f3p.abandoned - 2f0p.urban - 1.0f0p.cleared,
    #     wolf_snake = p -> -0.2f0p.cat + 0.2f0p.black_rat + 0.3f0p.mouse - 0.5f0p.urban + 0.3f0p.native,
    #     macaque =    p -> 1.0f0p.abandoned + 0.7f0p.forestry + 0.4f0p.native - 1.0f0p.urban - 0.8f0p.cleared
    # )[pred_keys],
    pred_funcs = (;
        cat =        p -> 2.0f0p.black_rat + 0.5f0p.norway_rat + 10f0p.urban + 2f0p.cleared,
        black_rat  = p -> -0.2f0p.cat - 0.1f0p.norway_rat + 0.5f0p.native + 0.3f0p.abandoned + 0.3f0p.forestry + 1p.urban,
        norway_rat = p -> -0.1f0p.cat - 0.1f0p.black_rat + 1.5f0p.urban - 0.2f0p.native,
        mouse =      p -> -0.3f0p.cat - 0.2f0p.black_rat - 0.2f0p.norway_rat + 0.8f0p.cleared + 1.5f0p.urban,
    )[pred_keys]
)
    pred_df = filter(r -> Symbol(r.name) in pred_keys, pred_df)
    moore = Moore{3}()

    #= Assumptions
    1. cats suppress rodents to some extent, black rats more than norway rats (size selection - norway rats are above 250g)
    2. cats live near people with maybe 2 orders of magnitude higher density than far from people
    3. black rats live anywhere including forest (they're efficient and good climbers), but are outcompeted by norway rats in cities and/or in the presence of cats
    4. norway rats could live everywhere but lose in the forest in competition with black rats, because theyre bad climbers and are less efficient in feeding/metabolism overall. They dominate in coastal areas because theyre better around water, and in high resource conditions because they're bigger.
    5. mice are also implicated in bird deaths, just to make things worse... they outcompete rats in farmland but not forests, and coexist in urban areas.
    6. pigs probably don't care about any of these things and just go wherever they like.
    =#

    # Parameters



    # These are parametries from reading qualitative literature and quantitative literature without
    # enough context to really use the number.
    # Maybe we should fit them to some observations in specific points at specific times.
    # Refs: Smucker et al 2000 - Hawaiii cats, rats, mice

    # pred_funcs = (;
    #     cat =        p -> 10f0p.urban + 2f0p.cleared,
    #     black_rat  = p -> 0.5f0p.native + 0.3f0p.abandoned + 0.3f0p.forestry + 1p.urban,
    #     norway_rat = p -> 1spec .5f0p.urban - 0.2f0p.native,
    #     mouse =      p -> 0.8f0p.cleared + 1.5f0p.urban,
    #     pig =        p -> 0.5f0p.native + 0.4f0p.abandoned - 2f0p.urban - 1.0f0p.cleared,
    #     wolf_snake = p -> 0.5f0p.urban + 0.3f0p.native,
    #     macaque =    p -> 1.0f0p.abandoned + 0.7f0p.forestry + 0.4f0p.native - 1.0f0p.urban - 0.8f0p.cleared
    # )[pred_keys]



    # These need to somewhat balance low growth rates. They are almost totally made up.
    # The units are in pixels - it needs fixing to the aggregation size.

    island_endemic_traits = endemic_traits(island_endemic_tables, EndemicNVs)
    # gecko_mass = 8 # estimated mean of multiple species
    # skink_mass = 3 # estimated mean of multiple species
    # mouse_mass = 16.25
    # wolf_snake_pre_mass = round(0.48gecko_mass + 0.30mouse_mass + 0.22skink_mass)

    # How much these species are supported outside of this system
    # populations = NV(cat=0.01, black_rat=25.0, norway_rat=10.0, pig=0.3)

    # suitabilities = (human_dependency => human_intensity, forest_preference => forest_density)
    # scale_carrycap(populations, carrycap, interactions, suitabilities)

    kernel = DispersalKernel(
        stencil=moore,
        formulation=ExponentialKernel(Param(1.0f0, bounds=(0.0f0, 2.0f0))),
        cellsize=1.0f0,
    )
    stencil_masks = map(auxs) do aux
        StencilArray(aux.mask, kernel; padding=Halo{:out}())
    end
    stencil_dems = map(auxs) do aux
        StencilArray(aux.dem, moore; padding=Halo{:out}())
    end

    PredNV = NamedVector{pred_keys,length(pred_keys)}
    pred_names = PredNV(pred_keys)
    pred_indices = PredNV(ntuple(identity, length(pred_keys)))
    island_names = NamedTuple{keys(auxs)}(keys(auxs))
    pred_init_nvs = map(pred_indices) do i
        x = zeros(Float16, length(pred_keys))
        x[i] = 50 * aggfactor
        PredNV(x)
    end

    introductions = map(island_names) do key
        island_df = filter(r -> r.island == string(key) && Symbol(r.species) in pred_keys, introductions_df)
        display(island_df)
        map(eachrow(island_df)) do r
            init = pred_init_nvs[Symbol(r.species)]
            (; year=r.year, geometry=(X=r.lon, Y=r.lat), init)
        end
    end

    pred_rmax = Float32.(PredNV(pred_df.rmax))
    pred_max_density = Float32.(PredNV(pred_df.max_density))
    pred_pops = map(auxs) do aux
        map(_ -> map(_ -> 0.0f0, pred_rmax), aux.mask)
    end
    pred_carrycaps = map(auxs) do aux
        map(_ -> carrycap, aux.mask)
    end

    pred_carrycap_rule = InteractiveCarryCap{:pred_pop,:pred_carrycap}(;
        carrycap,
        carrycap_scaling=map(Val, pred_funcs),
        inputs=Aux{:landcover}(),
    )

    pred_growth_rule = LogisticGrowth{:pred_pop}(;
        rate=pred_rmax,
        carrycap=DynamicGrids.Grid{:pred_carrycap}(),
        timestep=1,
        nsteps_type=Float32,
    )

    # Every species is everywhere initially, in this dumb model
    endemic_presences = map(auxs, island_endemic_traits) do aux, traits
        map(aux.mask) do m
            map(_ -> m, traits.ismammal) # traits.mammal is just as a map source, not used
        end
    end
    island_causes = map(auxs, island_endemic_traits) do aux, traits
        map(aux.mask) do m
            map(_ -> zero(pred_rmax), traits.ismammal) # traits.mammal is just as a map source, not used
        end
    end


    #### Rules ##########################################################333

    introduction_rule = let introductions_aux=Aux{:introductions}()
        SetGrid{:pred_pop}() do data, pred_pop, t
            D = dims(DG.init(data).pred_pop)
            current_year = currenttime(data)
            intros = get(data, introductions_aux)
            foreach(intros) do intro
                if intro.year == current_year
                    p = intro.geometry
                    I = DimensionalData.dims2indices(D, (X(Contains(p.X)), Y(Contains(p.Y))))[1:2]
                    # x = view(pred_pop, I..., :) .+ (intro.init,)
                    # @show x I size(intro.init) size(pred_pop)
                    pred_pop[I..., :] .= view(pred_pop, I..., :) .+ (intro.init,)
                end
            end
        end
    end

    habitat = let native=Aux{:native}()
        Cell{:presences}() do data, presences, I
            hp = get(data, native, I)
            habitat_requirement = DG.aux(data).habitat_requirement
            map(presences, habitat_requirement) do present, hs
                if present
                    rand() < hp * hs
                else
                    false
                end
            end
        end
    end

    clearing_rule = let landcover=Aux{:landcover}()
        Cell{:endemic_presence}() do data, presences, I
            lc = get(data, landcover, I)
            presences .& (lc.native > 0.2f0)
        end
    end

    # clearing_rule = let landcover=Aux{:landcover}(), native_needs=native_needs
    #     Cell{:endemic_presence}() do data, presences, I
    #         lc = get(data, landcover, I)
    #         presences .& map((lc, n) -> lc * n > rand(Float32), lc.native, native_needs)
    #     end
    # end

    pred_kernels = map(spread_rate) do s
        DispersalKernel(
            stencil=moore,
            # formulation=ExponentialKernel(Param(s, bounds=(0.000000000001, 100.0)))
            formulation=ExponentialKernel(s),
            cellsize=1.0f0,
        )
    end
    pred_spread_rule = let demaux=Aux{:dem}(), aggfactor=aggfactor, pred_kernels=pred_kernels, carrycap=carrycap
        SetNeighbors{Tuple{:pred_pop,:pred_carrycap}}(pred_kernels[1]) do data, hood, (Ns, _), I
            Ns === zero(Ns) && return nothing
            dem = get(data, demaux)
            reps = DynamicGrids.replicates(data)
            carrycap_nbrs = if isnothing(reps)
                DG.neighbors(DG.grids(data).pred_carrycap, I)
            else
                DG.neighbors(DG.grids(data).pred_carrycap, (I..., reps))
            end
            elev_nbrs = DG.neighbors(dem, I)
            @inbounds elev_center = dem[I...]

            sum = zero(Ns) # TODO needs cell size here
            cellsize = (100 * aggfactor)

            # Randomise hood starting position to avoid directional artifacts in output
            start = rand(0:length(hood)-1)
            @inbounds for ix in eachindex(hood)
                # Rotate indices in relation to starting point
                i = start + ix
                if i > length(hood)
                    i = i - length(hood)
                end
                sp = carrycap_nbrs[i] ./ carrycap
                any(map(isnan, sp)) && error("sp is NaN")
                Ih = DG.indices(hood, I)[i]
                ks = getindex.(DynamicGrids.kernel.(pred_kernels), i)
                d = DG.distances(hood)[i]
                e = elev_nbrs[i]
                propagules = trunc.(Float32.(Ns .* sp .* rand(typeof(Ns)) .^ 3 .* ks .* 40))
                sum1 = sum + propagules
                # If we run out of propagules
                if any(sum1 .> Ns)
                    propagules = min.(sum1, Ns) .- sum
                    sum = sum + propagules
                else
                    sum = sum1
                end
                add!(data[:pred_pop], propagules, Ih...)
            end
            @inbounds sub!(data[:pred_pop], sum, I...)
            return nothing
        end
    end # let

    risks_rule = ExtirpationRisks{Tuple{:endemic_presence,:causes}}(;
        f=tanh,
        stencil=Moore(1),
        traits=Aux{:endemic_traits}(),
        pred_response,
        pred_suscept=nothing,
        pred_effect=nothing,
        pred_pop=Grid{:pred_pop}(),
        pred_carrycap=carrycap,
        stochastic_extirpation=0.005f0/aggfactor,
        recouperation_rates=Aux{:recouperation_rates}(),
    )

    aux_pred_risks_rule = ExtirpationRisks{:endemic_presence,:causes}(;
        f=tanh,
        stencil=Moore(1),
        traits=Aux{:endemic_traits}(),
        pred_response,
        pred_suscept=nothing,
        pred_effect=Aux{:pred_effect}(),
        pred_pop=nothing,
        pred_carrycap=carrycap,
        stochastic_extirpation=0.005f0/aggfactor,
        recouperation_rates=Aux{:recouperation_rates}(),
    )
    tspans = map(introductions) do intros
        first_year:last_year
    end
    inits = map(pred_pops, pred_carrycaps, endemic_presences, island_causes) do pred_pop, pred_carrycap, endemic_presence, causes
        (; pred_pop, pred_carrycap, endemic_presence, causes)
    end
    pred_inits = map(pred_pops, pred_carrycaps) do pred_pop, pred_carrycap
        (; pred_pop, pred_carrycap)
    end
    endemic_inits = map(endemic_presences) do endemic_presence
        (; endemic_presence)
    end

    pred_ruleset = Ruleset(
        pred_carrycap_rule,
        introduction_rule,
        pred_spread_rule,
        pred_growth_rule;
        boundary=Remove()
    )

    endemic_ruleset = Ruleset(
        Chain(aux_pred_risks_rule, clearing_rule);
        boundary=Remove()
    )

    ruleset = Ruleset(
        DynamicGrids.rules(pred_ruleset)...,
        risks_rule, clearing_rule;
        boundary=Remove()
    )

    rules = (;
        introduction_rule,
        pred_carrycap_rule,
        pred_spread_rule,
        pred_growth_rule,
        risks_rule,
        aux_pred_risks_rule,
        clearing_rule,
    )

    pred_effects = map(pred_pops_aux, island_mass_response, island_endemic_traits) do pred_pop, mass_response, traits
        pr = ModelParameters.stripparams(pred_response)
        pred_suscept = predator_suceptibility(pr, traits)
        isnothing(pred_pop) ? nothing : generate_predator_effect(risks_rule.f, pred_pop, pred_suscept)
    end

    outputs_kw = map(island_names, tspans, auxs, pred_pops_aux, island_endemic_traits, pred_effects, island_recouperation_rates) do island, tspan, aux1, pred_pop, endemic_traits, pred_effect, recouperation_rates
        aux = (;
            introductions=getproperty(introductions, island),
            # habitat_requirement=endemictraits.habitat_requirements,
            dem=getproperty(stencil_dems, island),
            recouperation_rates,
            pred_pop,
            endemic_traits,
            pred_effect,
            landcover=aux1.lc
        )
        (; aux, mask=getproperty(stencil_masks, island), replicates, tspan)
    end
    outputs = map(inits, outputs_kw) do init, kw
        ResultOutput(init; kw...)
    end
    pred_outputs = map(pred_inits, outputs_kw) do init, kw
        if isnothing(replicates)
            trans_output = TransformedOutput(init; kw...) do f
                Array(f.pred_pop)
            end
        else
            ResultOutput(init; kw...)
        end
    end
    endemic_outputs = map(endemic_inits, outputs_kw) do init, kw
        ResultOutput(init; kw...)
    end
    outputs = map(inits, outputs_kw) do init, kw
        ResultOutput(init; kw...)
    end

    islands = map(
        island_keys, inits, endemic_inits, pred_inits, outputs, endemic_outputs, pred_outputs, outputs_kw, auxs, island_mass_response, island_extinction_dates
    ) do key, init, endemic_init, pred_init, output, endemic_output, pred_output, output_kw, aux, mass_response, extinction_dates
        (; key, init, endemic_init, pred_init, output, endemic_output, pred_output, output_kw, aux, mass_response, extinction_dates)
    end

    return (; ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response)
end
