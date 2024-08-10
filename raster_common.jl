println("Loading packages...")
using GeoJSON
using GeoInterface
using GeometryBasics
using GADM
using Shapefile
using RasterDataSources
using Rasters
using Rasters: Between, trim, Band
using Unitful
using Distributions
using DimensionalData
using DimensionalData.LookupArrays

println("Loading functions...")
includet("common.jl")
includet("slope.jl")
includet("functions.jl")
# includet("lost_land_images.jl")

years = 1638, 1773, 1835, 1872, 1935, "present"
lc_years = 1638, 1773, 1835, 1872, 1935, "present"
lc_year_keys = map(y -> "lc_$y", lc_years)

println("Getting borders...")
gdal_borders = (
    mus=GADM.get("MUS").geom[1],
    reu=GADM.get("REU").geom[1],
    rod=GADM.get("MUS").geom[1],
)
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
# Plots.plot(mauritius_dem)
reu_tile  = getraster(SRTM; bounds=island_bounds.reu)[1]
reu_dem = replace_missing(trim(view(Raster(reu_tile), border_selectors.reu...); pad=10))
rod_tile = getraster(SRTM; bounds=island_bounds.rod)[1]
rod_dem = replace_missing(trim(view(Raster(rod_tile), border_selectors.rod...); pad=10))
dems = map(fix_order, (mus=mus_dem, reu=reu_dem, rod=rod_dem))
elevation = map(d -> d .* u"m", dems)

# mauritius_proj_dem = Raster("/home/raf/PhD/Mascarenes/Data/Norder/LS factor/DEM/DEM100x100_Resample.img"; crs=EPSG(3337))
# mauritius_proj_dem = rebuild(mauritius_proj_dem; missingval=minimum(mauritius_proj_dem))

masks = map(d -> rebuild(boolmask(d); name=:mask), dems)
println("Calculating slope...")
slope_stacks = map(dems) do dem
    slopeaspect(dem, FD3Linear(); cellsize=111.0)
end

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

using Rasters, Plots, ArchGDAL, Extents
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
