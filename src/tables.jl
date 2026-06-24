using CSV
using DataFrames
using ConstructionBase

# parse_year_range and is_present are defined in common.jl

function load_tables()
    tables_path = realpath(joinpath(dirname(pathof(MascarenesSimulations)), "../tables"))
    pred_df          = CSV.read("$tables_path/invasives.csv",    DataFrame)
    introductions_df = CSV.read("$tables_path/introductions.csv", DataFrame)
    mascarene_species_csv = "$tables_path/mascarene_species.csv"

    all_species = CSV.read(mascarene_species_csv, DataFrame) |>
        x -> subset(x, :Species => ByRow(!ismissing); skipmissing=true)

    island_tables = map(island_keys) do key
        df = DataFrame(subset(all_species, key => x -> coalesce.(is_present.(x), false)))
        df.extinct    = parse_year_range.(df[!, "$(key)_extinct"])
        df.introduced = parse_year_range.(df[!, "$(key)_introduced"])
        df
    end
    island_endemic_tables = map(island_tables) do tbl
        # TODO: add missing-mass rows and remove the Mass filter
        DataFrame(subset(tbl, :Origin => ByRow(==("Endemic")), :Mass => ByRow(!ismissing); skipmissing=true))
    end

    get_species_names(table) = Tuple(Symbol.(replace.(skipmissing(table.Species), Ref(" " => "_"))))
    island_names         = NamedTuple{keys(island_tables)}(keys(island_tables))
    island_endemic_names = map(get_species_names, island_endemic_tables)
    all_endemic_names    = union(island_endemic_names...)

    (; pred_df, introductions_df, island_names, island_endemic_names, island_tables, island_endemic_tables)
end
