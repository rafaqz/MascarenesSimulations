println("Loading packages...")

using ArchGDAL
using DimensionalData
using GeoInterface
using GeometryBasics
using GADM
using Extents
using Rasters
using RasterDataSources
using Tables
using Unitful

using DimensionalData.LookupArrays
using Rasters: Between, trim, Band

println("Loading functions...")
include("common.jl")

years = 1638, 1773, 1835, 1872, 1935, "present"
lc_years = 1638, 1773, 1835, 1872, 1935, "present"
lc_year_keys = map(y -> "lc_$y", lc_years)

println("Getting borders...")
mus = GADM.get("MUS")
gdal_borders = map((
    mus=GADM.get("MUS"),
    reu=GADM.get("REU"),
    rod=GADM.get("MUS"),
)) do table
    Tables.columntable(table).geom[1]
end
borders = (
   mus=GeoInterface.convert(GeometryBasics, gdal_borders.mus),
   reu=GeoInterface.convert(GeometryBasics, gdal_borders.reu),
   rod=GeoInterface.convert(GeometryBasics, gdal_borders.rod),
)
island_bounds = (
    # mus=((57.1, 57.9), (-20.6, -19.8)), # with islands
    mus=((57.1, 57.9), (-20.6, -19.949)),
    reu=((55.0, 56.0), (-22.0, -20.0)),
    rod =((63.0, 64.0), (-20.0, -19.0)),
)
tiles = getraster(SRTM; bounds=island_bounds.mus)
dem1 = Raster(tiles[1]; name=:DEM)
dem2 = Raster(tiles[2]; name=:DEM)
border_selectors = map(island_bounds) do bb
    (X(Between(bb[1])), Y(Between(bb[2])))
end

println("Getting DEMs...")
# Mauritius is right over the split in the tiles
m1 = view(dem1, border_selectors.mus...)
m2 = view(dem2, border_selectors.mus...)
mus_dem = replace_missing(trim(cat(m1, m2; dims=Y); pad=10))
reu_tile  = getraster(SRTM; bounds=island_bounds.reu)[1]
reu_dem = replace_missing(trim(view(Raster(reu_tile), border_selectors.reu...); pad=10))
rod_tile = getraster(SRTM; bounds=island_bounds.rod)[1]
rod_dem = replace_missing(trim(view(Raster(rod_tile), border_selectors.rod...); pad=10))
dems = map(fix_order, (mus=mus_dem, reu=reu_dem, rod=rod_dem))
elevation = map(d -> d .* u"m", dems)

masks = map(d -> rebuild(boolmask(d); name=:mask), dems)

println("Getting native vegetation...")
mus_native_veg_tif_path = "/home/raf/PhD/Mascarenes/Data/Generated/mus_native_veg.tif"
reu_native_veg_tif_path = "/home/raf/PhD/Mascarenes/Data/Generated/reu_all_natives.tif"
native_veg = (;
    mus=Raster(mus_native_veg_tif_path) ,
    reu=Raster(reu_native_veg_tif_path),
)

mus_veg_path = "/home/raf/PhD/Mascarenes/Data/Selected/Mauritius/Undigitised/page33_mauritius_vegetation_colored.tif"
reu_veg_path = "/home/raf/PhD/Mascarenes/Data/Dominique/Vegetation_Rasters/pastveg3.tif"
original_veg = (;
    mus=reorder(replace_missing(Raster(mus_veg_path), 0), masks.mus),
    reu=reorder(resample(replace_missing(Raster(reu_veg_path), 0); to=masks.reu), masks.reu),
)

mus_rod_pop_density = Raster("/home/raf/PhD/Mascarenes/Data/Population/raster/population_mus_2018-10-01.tif"; lazy=true)
reu_pop_density = Raster("/home/raf/PhD/Mascarenes/Data/Population/raster/population_reu_2018-10-01.tif"; lazy=true)
pop_density = map(fix_order, (
    mus=mus_rod_pop_density[Extents.extent(dems.mus)],
    reu=reu_pop_density[Extents.extent(dems.reu)],
    rod=mus_rod_pop_density[Extents.extent(dems.rod)],
))


# Homiisland Landcover
lc_dir = joinpath(datadir, "Landcover/")
lc_names = (
  :No_Data,
  :Continuous_urban,
  :Discontinuous_urban,
  :Forest,
  :Shrub_vegetation,
  :Herbaceaous_vegetation,
  :Mangrove,
  :Barren_land,
  :Water,
  :Sugarcane,
  :Pasture,
  :UnusedIndex,
  :Other_cropland,
)
lc_2017_category_groups = (
    forest_or_abandoned = (:Forest, :Mangrove, :Shrub_vegetation, :Barren_land),
    urban = (:Continuous_urban, :Discontinuous_urban),
    cleared = (:Sugarcane, :Other_cropland, :Pasture),
    uncertain = (:Herbaceaous_vegetation,),
)

lc_2017_categories = NamedTuple{lc_names}((Int32.(0:12)...,))
