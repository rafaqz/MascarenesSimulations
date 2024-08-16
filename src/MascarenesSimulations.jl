module MascarenesSimulations

using ArchGDAL
using ColorSchemes
using Colors
using DimensionalData
using DynamicGrids
using Dispersal
using Distributions
using Extents
using GeoInterface
using GeometryBasics
using GADM
using JLD2
using LandscapeChange
using Makie
using ModelParameters
using NCDatasets
using Rasters
using RasterDataSources
using Setfield
using StaticArrays
using Statistics
using Stencils
using Tables
using ThreadsX
using Unitful

using DimensionalData.LookupArrays
using Rasters: Between, trim, Band

const DG = DynamicGrids
const DD = DimensionalData
const NV = NamedVector

export load_tables, load_rasters, load_aux, graphic_landcover

export define_simulations

export makie_sim

const basepath = realpath(joinpath(@__DIR__, ".."))

include("common.jl")
include("functions.jl")
include("aux.jl")
include("rules.jl")
include("rasters.jl")
include("tables.jl")
include("makie.jl")

end
