
# Build auxiliary rasters

# Apply boolean mask corrections passed in as data.
# Each patch is a NamedTuple with:
#   island :: Symbol          — key into the masks NamedTuple
#   dims   :: Tuple           — DimensionalData selectors (e.g. X(Between(a,b)))
# Cells matching all selectors are set to false (excluded from simulation).
function apply_mask_patches!(masks, patches)
    for patch in patches
        view(getproperty(masks, patch.island), patch.dims...) .= false
    end
    return masks
end

function load_aux(; masks, dems, landcover_paths, aggfactor, last_year, mask_patches=())
    sim_setup_file = joinpath(basepath, "cache/sim_setup_$aggfactor.jld2")
    auxs = if isfile(sim_setup_file)
        println("Loading aux data from jld...")
        let
            f = jldopen(sim_setup_file, "r")
            auxs = f["auxs"]
            close(f)
            auxs
        end
    else
        let
            # netcdf has the annoying center locus for time
            lc_predictions = map(landcover_paths) do path
                RasterStack(path) |>
                    x -> maybeshiftlocus(Start(), x) |>
                    x -> DD.set(x, Ti => Int.(lookup(x, Ti))) |>
                    x -> rebuild(Rasters.modify(BitArray, x); missingval=false)
            end
            apply_mask_patches!(masks, mask_patches)
            auxs = agg_aux(masks, dems, lc_predictions, aggfactor, last_year)
            jldsave(sim_setup_file; auxs)
            auxs
        end
    end
end

# Merge landcover to a single integer layer for Makie visualisation.
function graphic_landcover(auxs)
    map(auxs) do aux
        map(aux.lc) do lcs
            sum(map(.*, ntuple(UInt8, length(lcs)), lcs))
        end
    end
end
