function mk(init, ruleset; carrycaps, kw...)
    MakieOutput(init;
        kw...,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do (; layout, frame, time)
        pred_keys = propertynames(frame[].pred_pop[1])
        ncols = length(pred_keys)
        extinct_keys = propertynames(frame[].endemic_presence[1])
        n_extinct = length(extinct_keys)
        extinct_strings = collect(string.(extinct_keys))
        menus = map(1:ncols) do i
            Menu(layout[3, i]; default=extinct_strings[i], options=extinct_strings)
        end
        pred_axes = map(1:ncols) do i
            Axis(layout[1, i]; title=string(pred_keys[i]))
        end
        extinct_axes = map(1:ncols, menus) do i, m
            ax = Axis(layout[2, i]; title=m.selection)
        end
        linkaxes!(pred_axes..., extinct_axes...)
        predators = map(1:ncols) do i
            Observable(Array((x -> iszero(x) ? NaN : Float64(x)).(getindex.(frame[].pred_pop, i))))
        end
        extincts = map(1:ncols) do i
            Observable(Array(getindex.(frame[].endemic_presence, i)))
        end

        on(frame) do f
            foreach(predators, 1:ncols) do  pred, i
                pred[] .= (x -> iszero(x) ? NaN : Float64(x)).(getindex.(frame[].pred_pop, i))
                notify(pred)
            end
        end
        foreach(extincts, menus) do extinct, menu
            onany(frame, menu.selection) do f, selection
                i = findfirst(==(selection), extinct_strings)
                extinct[] .= getindex.(f.endemic_presence, i)
                notify(extinct)
            end
        end
        foreach(pred_axes, predators, pred_keys, carrycaps) do ax, pred, k, cc
            Makie.image!(ax, pred; colorrange=(0.0, cc), colormap=:magma, interpolate=false)
        end
        foreach(extinct_axes, extincts, Iterators.Cycle([:blues, :reds, :greens])) do ax, extinct, colormap
            Makie.image!(ax, extinct; colorrange=(0.0, 1.0), colormap, interpolate=false)
        end
        return nothing
    end
end

function mk_pred(init, ruleset, mask; carrycap, kw...)
    MakieOutput(init;
        kw...,
        fps=100,
        store=false,
        ruleset,
        sim_kw=(; printframe=true),
    ) do fig, frame, time
        # nanmask = map(x -> x ? 1.0f0 : NaN32, mask)
        pred_keys = propertynames(frame[].pred_pop[1])
        ncols = length(pred_keys)
        pred_axes = map(1:ncols) do i
            Axis(fig[1, i]; title=string(pred_keys[i]))
        end
        linkaxes!(pred_axes...)
        predators = map(1:ncols) do i
            Observable(Array(Float32.(getindex.(frame[].pred_pop, i))))
        end

        on(frame) do f
            foreach(predators, 1:ncols) do pred, i
                pred[] .= getindex.(f.pred_pop, i)
                notify(pred)
            end
        end
        foreach(pred_axes, predators, pred_keys, carrycap) do ax, pred, k, cc
            Makie.image!(ax, pred; colorrange=(0, cc), colormap=:viridis, interpolate=false)
        end
        return nothing
    end
end


