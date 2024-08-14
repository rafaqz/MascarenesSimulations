module MascarenesSimulations

using ArchGDAL
using ColorSchemes
using Colors
using DimensionalData
using Dispersal
using Extents
using GeoInterface
using GeometryBasics
using GADM
using Makie
using JLD2
using NCDatasets
using Rasters
using RasterDataSources
using StaticArrays
using Statistics
using Stencils
using Tables
using Unitful

using DimensionalData.LookupArrays
using Rasters: Between, trim, Band

export load_tables, load_rasters, load_aux, graphic_landcover

export define_simulations

export makie_sim

include("species_common.jl")
include("species_tables.jl")
include("species_rules.jl")
include("raster_common.jl")
include("makie.jl")

end
