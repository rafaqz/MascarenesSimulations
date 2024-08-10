
namedkeys(nt::NamedTuple{K}) where K = NamedTuple{K}(K)

function plot_lc(lc_raster::Raster)
    Plots.plot(lc_raster; 
        color=palette(:Paired_12, 12), 
        clims=(1, 13), size=(2000,2000),
        colorbar_ticks=map(Pair, 0:12, lc_categories),
    )
end

function plot_lc_makie(lc_raster::Raster)
    fig = Makie.Figure()
    ax, hm = Makie.heatmap(fig[1, 1], parent(parent(dims(lc_raster, X))), parent(parent(dims(lc_raster, Y))), parent(read(lc_raster)),
        colormap=cgrad(:cyclic_mygbm_30_95_c78_n256, 13, categorical=true), colorrange=(0, 13)
    )
    # ax.aspect = Makie.AxisAspect(1)
    # Makie.Colorbar(fig[1, 2], hm; 
        # ticks=(0:12, lc_categories),
    # )
    return fig
end

