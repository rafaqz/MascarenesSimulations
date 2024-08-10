# This file precomputes and save a lot of high-cost rasters
# so we can just load them from disk during normal work.

include("raster_common.jl")
include("ports.jl")
include("nearest.jl")

# Distances

# To coasts
distance_to_coasts = map(dems, borders) do dem, border 
    coast = boolmask(border; to=dem, shape=:line)
    dists = nearest_distances(coast)
    rebuild(mask(dists; with=dem); name=:to_coasts)
end
foreach(island_keys, distance_to_coasts) do i, raster
    mkpath(joinpath(distancedir, string(i)))
    write(joinpath(distancedir, string(i), "to_coast.tif"), raster)
end

foreach(island_keys, dems, ports) do i, dem, classes
    foreach(namedkeys(classes), classes) do k, locations
        ports = falses(dims(dem))
        for loc in locations
            selectors = Y(Contains(loc[1])), X(Contains(loc[2]))
            if DimensionalData.hasselection(ports, selectors)
                ports[selectors...] = true
            end
        end
        dists = nearest_distances(ports)
        rast = rebuild(mask(dists; with=dem); name=Symbol("to_$(k)_ports"))
        dir = joinpath(distancedir, string(i))
        write(joinpath(dir, "to_$(k)_ports.tif"), rast)
    end
end

distance_to_water = map(island_keys, dems, watermasks) do i, dem, watermask
    rast = mask(nearest_distances(watermask); with=dem)
    write(joinpath(distancedir, string(i), "to_water.tif"), rast)
    rast
end

# Landcover
lc_dir = joinpath(datadir, "Landcover/")
lc_shapes = (
    mus=(
        shapepath=joinpath(lc_dir, "mus/cla_maurice_fin.shp"), crspath=joinpath(lc_dir, "mus/cla_maurice_fin.prj"),), reu=(
        shapepath = joinpath(lc_dir, "reu/cla_run_2014_fin_2975.shp"),
        crspath = joinpath(lc_dir, "reu/cla_run_2014_fin_2975.prj"),
    ),
    rod=(
        shapepath = joinpath(lc_dir, "rod/cla_rod_fin.shp"),
        crspath = joinpath(lc_dir, "rod/cla_rod_fin.prj"),
    )
)
lc_categories = NamedTuple{lc_names}((Int32.(0:12)...,))

function rasterize_lc(template, shape_file, crs_file; res=nothing, categories)
    lc_shape = Shapefile.Table(shape_file)
    lc_crs = WellKnownText(readlines(crs_file)[1])
    lc_df = DataFrame(lc_shape)
    lc_raster = Raster(similar(template, Int32); missingval=typemin(Int32))
    lc_raster .= typemin(Int32)
    if !isnothing(res)
        lc_raster = read(resample(lc_raster, res; crs=lc_crs))
    end
    display(dims(lc_raster))
    # Order of rasterization matters?... (probably should calculate areas?)
    fillvals = [
        categories.No_Data,
        categories.Water,
        categories.Herbaceaous_vegetation,
        categories.Shrub_vegetation,
        categories.Barren_land,
        categories.Other_cropland,
        categories.Sugarcane,
        categories.Pasture,
        categories.Forest,
        categories.Mangrove,
        categories.Continuous_urban,
        categories.Discontinuous_urban,
    ]
    for fillval in fillvals
        rows = filter(x -> x.ocsol_num == fillval, lc_df)
        if length(rows.geometry) > 0
            fillname = first(eachrow(rows)).ocsol_name
            # This should use `last` but its not optimised currently
            rasterize!(last, lc_raster, rows.geometry; fill=fillval, threaded=false)
        end
    end
    return lc_raster
end

lc_rasterized_high_res = map(dems, lc_shapes) do dem, (; shapepath, crspath)
    rasterize_lc(dem, shapepath, crspath; categories=lc_categories, res=90) |> read
end
lc_rasterized_high_res.rod |> Rasters.rplot
lc_rasterized = map(lc_rasterized_high_res, dems) do A, to
    resample(A; to, method=:near)
end

# Write to tif
foreach(lc_rasterized, island_keys) do lc, island
    write(joinpath(datadir, "Generated", "Landcover", "$(island)_landcover.tif"), lc; force=true)
    rasters = Dict()
end

# Mauritius landcover Shapefiles
tbl = Shapefile.Table("/home/raf/PhD/Mascarenes/Data/Claudia/Black River Gorges and other shapes/npcs.shp") |> DataFrame
tbl = Shapefile.Table("/home/raf/PhD/Mascarenes/Data/Claudia/Demo/GIS WILD LIFE FOUNDATION/SHAPE FILE/land_use_WGS_region.shp") |> DataFrame
union(tbl.LAND_USE)
swamps = filter(:LAND_USE => x -> !ismissing(x) && x == "Marsh or Swamp", tbl)
cleared = filter( :LAND_USE => 
  x -> !ismissing(x) && x in ("Other Plantation", "Sugar Cane", "Other PlanTation", "Tea Plantation", "Sugar cane"), tbl)
# swamps_rast = boolmask(swamps; to=masks.mus)
cleared_rast = boolmask(cleared; to=masks.mus)
other_rast = .!(cleared_rast)
lc = cleared_rast .* 1 .+ other_rast .* 2 .* masks.mus
write("/home/raf/PhD/Mascarenes/Data/Generated/Landcover/mus_wlf_shape.tif", lc)

# ag = Raster("/home/raf/PhD/Mauritius/Data/Deforestation-Mauritius/rasters/agr_suitab/agr_suitab_h/w001001.adf"; crs=EPSG(3337))
# resample(ag; to=dems.mus) |> plot

natparks_paths = (
    mus="/home/raf/PhD/Mascarenes/Data/NationalParks/mus.geojson", 
    reu="/home/raf/PhD/Mascarenes/Data/NationalParks/reu.geojson",
)
map(natparks_paths, island_keys[(:mus, :reu)]) do path, island
    reu_natparks_json = GeoJSON.read(path)
    natparks_raster = UInt8.(boolmask(natparks_json; to=dems.reu))
    write("/home/raf/PhD/Mascarenes/Data/Generated/NationalParks/$island.tif", natparks_raster; force=true)
end
Rasters.rplot(natparks_raster)

reu_veg_rast = resample(Raster("/home/raf/PhD/Mascarenes/Data/Dominique/Vegetation_Rasters/present_veg4.tif"); to=dems.reu)
write("/home/raf/PhD/Mascarenes/Data/Generated/reu_native.tif", reu_veg_rast; force=true)
reu_all_natives = rebuild(UInt8.((x -> x in 1:20).(reu_veg_rast)); missingval=0x00)
write("/home/raf/PhD/Mascarenes/Data/Generated/reu_all_natives.tif", reu_all_natives; force=true)

reu_veg_rast = Raster("/home/raf/PhD/Mascarenes/Data/Generated/reu_native.tif")
