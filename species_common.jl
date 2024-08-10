using Statistics
using StaticArrays
using Dispersal
using XLSX
using TerminalPager
using Rasters
using GLMakie
using NCDatasets
using Geomorphometry
using Stencils
using Statistics
using Setfield
using ConstructionBase
using JLD2

island_keys = NamedTuple{(:mus, :reu, :rod)}((:mus, :reu, :rod))

# includet("optimisation.jl")
include("species_rules.jl")
include("species_tables.jl")
include("raster_common.jl")
# forested = gpu_cleanup(modify(BitArray, Raster("forested.nc")))
# Set aggregation
aggfactor = 16
# And the last year of the simulation
first_year = 1550
last_year = 2018
extant_extension = 0

# Build auxiliary rasters
lc_predictions_paths = (
    mus="$outputdir/lc_predictions_mus.nc",
    reu="$outputdir/lc_predictions_reu.nc",
    rod="$outputdir/lc_predictions_rod.nc",
)

sim_setup_file = "sym_setupj$aggfactor.jld2"

auxs = if isfile(sim_setup_file)
    let
        f = jldopen(sim_setup_file, "r")
        # pred_pops_aux = f["pred_pops_aux"];
        auxs = f["auxs"];
        close(f)
        auxs
    end
else
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

# just for makie visualisation, kind of merges the pixel colors...
lc_all = map(auxs) do aux
    map(aux.lc) do lcs
        sum(map(.*, ntuple(UInt8, length(lcs)), lcs))
    end
end
