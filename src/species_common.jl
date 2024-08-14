
using Dispersal
using JLD2
using StaticArrays
using Statistics
using Stencils

# includet("optimisation.jl")
include("raster_common.jl")
include("species_rules.jl")
include("species_tables.jl")

# landcover paths
lc_predictions_paths = (
    mus="$outputdir/lc_predictions_mus.nc",
    reu="$outputdir/lc_predictions_reu.nc",
    rod="$outputdir/lc_predictions_rod.nc",
)

island_keys = NamedTuple{(:mus, :reu, :rod)}((:mus, :reu, :rod))

# forested = gpu_cleanup(modify(BitArray, Raster("forested.nc")))

# Set aggregation
aggfactor = 16

# And the last year of the simulation
first_year = 1550
last_year = 2018
extant_extension = 0


# Build auxiliary rasters
sim_setup_file = "../cache/sym_setup_$aggfactor.jld2"

auxs = if isfile(sim_setup_file)
    println("Loading aux data from jld...")
    let
        f = jldopen(sim_setup_file, "r")
        # pred_pops_aux = f["pred_pops_aux"];
        auxs = f["auxs"];
        close(f)
        auxs
    end
else
    println("Loading and aggregating aux data manually...")
    let
        # netcdf has the annoying center locus for time
        lc_predictions = map(lc_predictions_paths) do path
            lc_predictions = RasterStack(path) |>
                x -> maybeshiftlocus(Start(), x) |>
                x -> DD.set(x, Ti => Int.(lookup(x, Ti))) |>
                x -> rebuild(Rasters.modify(BitArray, x); missingval=false)
        end
        # # Remove islands of rodrigues
        masks.rod[X=60 .. 63.33, Y=-19.775 .. -19.675] .= false
        masks.rod[X=63 .. 65, Y= -19.8 .. -19.775] .= false
        # And one pixel in Muaritius that creates a bug, false the whole row
        masks.mus[Y = -19.985 .. -19] .= false
        # Makie.plot(masks.mus)
        auxs = agg_aux(masks, dems, lc_predictions, aggfactor, last_year)
        jldsave(sim_setup_file;
            auxs#, pred_pops_aux,
        );
        auxs
    end
end

# Merge landcover to a single layer for makie visualisation
lc_graphic = map(auxs) do aux
    map(aux.lc) do lcs
        sum(map(.*, ntuple(UInt8, length(lcs)), lcs))
    end
end
