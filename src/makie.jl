Makie.set_theme!(theme_light())

const COLORMAPS = [:reds, :blues, :greens, :magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]

function makie_sim(init, ruleset; maxpops=zero(eltype(init.pred_pop)), landcover=nothing, tspan, kw...)
    MakieOutput(init;
        kw...,
        tspan,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
) do (; layout, frame, time)

        colorrange_obs = map(x -> Observable((zero(x), oneunit(x))), maxpops)
        ax_lc = Axis(layout[1, 1]; title="Landcover")
        if !isnothing(landcover)
            # Landcover
            lc = lift(time) do i
                replace_missing(landcover[Ti(Near(tspan[i]))], NaN32)
            end
            # hidexdecorations!(ax_lc; grid=false)
            # hideydecorations!(ax_lc; grid=false)
            Makie.image!(ax_lc, lc; colormap=:batlow, colorrange=(0, 6), interpolate=false)
        end

        # Predators
        pred_keys = propertynames(frame.pred_pop[][1])
        npreds = length(pred_keys)
        ncols = npreds + 1
        pred_axes = map(2:ncols) do i
            Axis(layout[1, i]; title=_title(pred_keys[i-1]))
        end
        hidexdecorations!.(pred_axes; grid=false)
        hideydecorations!.(pred_axes; grid=false)
        predators = map(1:npreds) do i
            Observable(rebuild(init.pred_pop, (x -> iszero(x) ? NaN : Float64(x)).(getindex.(frame.pred_pop[], i))))
        end
        on(frame.pred_pop) do pred_pop
            foreach(maximum(pred_pop), colorrange_obs) do m, obs 
                obs[] = (obs[][1], max(m, obs[][2]))
                notify(obs)
            end
            foreach(predators, 1:npreds) do pred, i
                pred[] .= (x -> iszero(x) ? NaN : Float64(x)).(getindex.(pred_pop, i))
                notify(pred)
            end
        end
        foreach(2:ncols, pred_axes, predators, pred_keys, COLORMAPS[1:npreds], colorrange_obs) do i, ax, pred, k, colormap, cr
            p = Makie.image!(ax, pred; colormap, colorrange=cr, interpolate=false)
            Colorbar(layout[1, i, Right()], p)
        end

        # Endemics
        extinct_keys = propertynames(frame.endemic_presence[][1])
        n_extinct = length(extinct_keys)
        extinct_strings = collect(string.(extinct_keys))
        menus = map(1:ncols) do i
            Menu(layout[3, i]; default=extinct_strings[i], options=extinct_strings)
        end
        extinct_axes = map(1:ncols, menus) do i, m
            title = lift(m.selection) do s
                _title(s)
            end
            ax = Axis(layout[2, i]; title)
        end
        hidexdecorations!.(extinct_axes; grid=false)
        hideydecorations!.(extinct_axes; grid=false)
        extincts = map(1:ncols) do i
            Observable(replace_missing(rebuild(init.pred_pop, Float32.(getindex.(frame.endemic_presence[], i))), NaN32))
        end
        causes = map(extincts) do e
            Observable((_ -> RGBA(0.0f0, 0.0f0, 0.0f0, 0.0f0)).(DimArray(e[])))
        end

        foreach(extincts, causes, menus) do extinct, cause, menu
            onany(frame.causes, menu.selection) do causes_vecs, selection
                i = findfirst(==(selection), extinct_strings)
                extinct[] .= replace(Float32.(getindex.(frame.endemic_presence[], i)), 0.0f0 => NaN32)
                function to_rgba(x, i)
                    sp = x[i]
                    pred_causes = (sp.cat, sp.norway_rat, sp.black_rat)
                    if any(map(>(0), pred_causes))
                        RGBA(min.(1.0, pred_causes ./ sum(pred_causes) .* 0.7)..., 1.0f0)
                    else
                        RGBA(pred_causes..., 0.0f0)
                    end
                end
                cause[] .= to_rgba.(causes_vecs, i)
                notify(extinct)
                notify(cause)
            end
        end

        endemic_cmaps = map(1:ncols) do i
            cgrad(ColorScheme([RGB{Float64}(0.0, 0.0, 0.0), RGB{Float64}(i, 0.1i, 1/i)]), 2, categorical=true)
        end
        foreach(extinct_axes, extincts, causes, endemic_cmaps) do ax, extinct, cause, colormap
            Makie.image!(ax, extinct; colorrange=(0.0, 1.0), colormap=:solar, interpolate=false)
            Makie.image!(ax, cause; interpolate=false)
        end

        # Link
        linkaxes!(pred_axes..., extinct_axes..., ax_lc)
        return nothing
    end
end

function mk_pred(init, ruleset; maxpops=zero(eltype(init.pred_pop)), landcover=nothing, tspan, kw...)
    MakieOutput(init;
        kw...,
        tspan,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do (; layout, frame, time)

        colorrange_obs = map(x -> Observable((zero(x), oneunit(x))), maxpops)
        # Landcover
        ax_lc = Axis(layout[1, 1]; title="Landcover")
        if !isnothing(landcover)
            lc = lift(time) do i
                replace_missing(landcover[Ti(Near(tspan[i]))], NaN32)
            end
            hidexdecorations!(ax_lc; grid=false)
            hideydecorations!(ax_lc; grid=false)
            Makie.image!(ax_lc, lc; colormap=:batlow, colorrange=(0, 6), interpolate=false)
        end

        # Predators
        pred_keys = propertynames(frame[].pred_pop[1])
        npreds = length(pred_keys)
        ncols = npreds + 1
        pred_axis_inds = CartesianIndices((2, ncols ÷ 2))[2:end]
        pred_sums = map(1:npreds) do _
            Observable("")
        end
        pred_axes = map(enumerate(pred_axis_inds), pred_sums) do (i, I), xlabel
            Axis(layout[Tuple(I)...]; title=_title(pred_keys[i]), xlabel)
        end
        # hidexdecorations!.(pred_axes; grid=false)
        hideydecorations!.(pred_axes; grid=false)
        pred_obs = map(1:npreds) do i
            Observable(rebuild(init.pred_pop, (x -> iszero(x) ? NaN : Float64(x)).(getindex.(frame[].pred_pop, i))))
        end
        on(frame.pred_pop) do pred_pop
            foreach(maximum(pred_pop), colorrange_obs) do m, obs 
                obs[] = (obs[][1], max(m, obs[][2]))
                notify(obs)
            end
            foreach(pred_obs, 1:npreds) do pred, i
                pred[] .= (x -> iszero(x) ? NaN : Float64(x)).(getindex.(pred_pop, i))
                pred_sums[i][] = string(sum(x -> x[i], pred_pop))
                notify(pred)
                notify(pred_sums[i])
            end
        end
        foreach(pred_axis_inds, pred_axes, pred_obs, COLORMAPS[1:npreds], colorrange_obs) do i, ax, pred, colormap, cr
            p = Makie.image!(ax, pred; colormap, colorrange=cr, interpolate=false)
            Colorbar(layout[Tuple(i)..., Right()], p)
        end

        # Link
        linkaxes!(ax_lc, pred_axes...)
        return nothing
    end
end

function mk_endemic(init, ruleset; 
    landcover=nothing, pred_pop, tspan,
    maxpops=map(i -> maximum(getindex.(pred_pops_aux.mus, i)), 1:length(first(pred_pop))),
    kw...
)
    MakieOutput(init;
        kw...,
        tspan,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do (; layout, frame, time)
          
        ax_lc = Axis(layout[1, 1]; title="Landcover")
        # Landcover
        if !isnothing(landcover)
            lc = lift(time) do i
                replace_missing(landcover[Ti(Near(tspan[i]))], NaN32)
            end
            # hidexdecorations!(ax_lc; grid=false)
            # hideydecorations!(ax_lc; grid=false)
            Makie.image!(ax_lc, lc; colormap=:batlow, colorrange=(0, 6), interpolate=false)
        end

        # Predators
        pred_keys = propertynames(pred_pop[1])
        npreds = length(pred_keys)
        ncols = npreds + 1
        pred_axes = map(2:ncols) do i
            Axis(layout[1, i]; title=_title(pred_keys[i-1]))
        end
        hidexdecorations!.(pred_axes; grid=false)
        hideydecorations!.(pred_axes; grid=false)
        pred_pop1 = view(pred_pop, Ti=1)
        predator_obs = map(1:npreds) do i
            Observable(rebuild(pred_pop1, (x -> iszero(x) ? NaN32 : Float32(x)).(getindex.(pred_pop1, i))))
        end
        on(time) do t
            foreach(predator_obs, 1:npreds) do pred, i
                pred_pop_t = view(pred_pop, Ti=t)
                pred[] .= (x -> iszero(x) ? NaN32 : Float32(x)).(getindex.(pred_pop_t, i))
                notify(pred)
            end
        end
        foreach(2:ncols, pred_axes, predator_obs, pred_keys, maxpops) do i, ax, pred, k, maxpop
            p = Makie.image!(ax, pred; colormap=:navia, colorrange=(0, maxpop), interpolate=false)
            Colorbar(layout[1, i, Right()], p)
        end

        # Endemics
        extinct_keys = propertynames(frame[].endemic_presence[1])
        n_extinct = length(extinct_keys)
        extinct_strings = collect(string.(extinct_keys))
        menus = map(1:ncols) do i
            Menu(layout[3, i]; default=extinct_strings[i], options=extinct_strings)
        end
        extinct_axes = map(1:ncols, menus) do i, m
            title = lift(m.selection) do s
                _title(s)
            end
            ax = Axis(layout[2, i]; title)
        end
        hidexdecorations!.(extinct_axes; grid=false)
        hideydecorations!.(extinct_axes; grid=false)
        extincts = map(1:ncols) do i
            Observable(replace_missing(rebuild(init.pred_pop, Float32.(getindex.(frame[].endemic_presence, i))), NaN32))
        end

        foreach(extincts, causes, menus) do extinct, cause, menu
            onany(frame, menu.selection) do f, selection
                i = findfirst(==(selection), extinct_strings)
                extinct[] .= replace(Float32.(getindex.(f.endemic_presence, i)), 0.0f0 => NaN32)
                if haskey(init, :causes)
                    cause[] .= f.causes
                end
                notify(extinct)
            end
        end

        endemic_cmaps = map(1:ncols) do i
            cgrad(ColorScheme([RGB{Float64}(0.0, 0.0, 0.0), RGB{Float64}(i, 0.1i, 1/i)]), 2, categorical=true)
        end
        foreach(extinct_axes, extincts, causes, endemic_cmaps) do ax, extinct, cause, colormap
            Makie.image!(ax, extinct; colorrange=(0.0, 1.0), colormap, interpolate=false)
            Makie.image!(ax, cause; colorrange=(0.0, 1.0), colormap, interpolate=false, alpha=0.5)
        end

        # Link
        linkaxes!(pred_axes..., extinct_axes..., ax_lc)
        return nothing
    end
end

_title(s) = replace(titlecase(string(s)), "_"=>" ")
