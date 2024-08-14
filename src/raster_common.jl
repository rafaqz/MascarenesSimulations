
function load_rasters()
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

    println("Getting DEMs...")
    tiles = getraster(SRTM; bounds=island_bounds.mus)
    dem1 = Raster(tiles[1]; name=:DEM)
    dem2 = Raster(tiles[2]; name=:DEM)
    border_selectors = map(island_bounds) do bb
        (X(Between(bb[1])), Y(Between(bb[2])))
    end
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

    (; borders, masks, elevation, dems)
end

fix_order(A) = reorder(A, ForwardOrdered)
