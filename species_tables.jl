using CSV
using DataFrames
using ConstructionBase

include("common.jl")

# Import tabular data
pred_df = CSV.read("tables/invasives.csv", DataFrame)
introductions_df = CSV.read("tables/introductions.csv", DataFrame)
mascarene_species_csv = "tables/mascarene_species.csv"

all_species = CSV.read(mascarene_species_csv, DataFrame) |> 
    x -> subset(x, :Species => ByRow(!ismissing); skipmissing=true)
endemic_species = CSV.read(mascarene_species_csv, DataFrame) |> 
    x -> subset(x, :Origin => ByRow(==("Endemic")); skipmissing=true)
island_tables = map(island_keys) do key
    df = DataFrame(subset(all_species, key => x -> .!ismissing.(x)))
    df.extinct = map(df[!, "$(key)_extinct"]) do e
        ismissing(e) ? missing : eval(Meta.parse(e))::UnitRange
    end
    df.introduced = map(df[!, "$(key)_introduced"]) do e
        ismissing(e) ? missing : eval(Meta.parse(e))::UnitRange
    end
    df
end
island_endemic_tables = map(island_tables) do tbl
    # TODO add missing mass rows and remove the missing Mass filter
    DataFrame(subset(tbl, :Origin => ByRow(==("Endemic")), :Mass => ByRow(!ismissing); skipmissing=true))
end
get_species_names(table) = Tuple(Symbol.(replace.(skipmissing(table.Species), Ref(" " => "_"))))
island_names = NamedTuple{keys(island_tables)}(keys(island_tables))
island_endemic_names = map(get_species_names, island_endemic_tables)
all_endemic_names = union(island_endemic_names...)

nspecies = length(all_endemic_names)
nislandspecies = map(length, island_endemic_names)
nislandcounts = map(island_names, island_endemic_names) do k, names
    isonisland = map(names) do name
        map(island_names -> name in island_names, island_endemic_names)
    end
    per_island = map(island_names) do k1
        count(x -> x[k1], isonisland)
    end
    onall = count(all, isonisland)
    per_island = map(x -> x - onall, per_island)
    endemic = per_island[k] - sum(per_island[Not(k)])
    out = (; per_island..., onall)
    # Adjust out so current island is only endemics
    ConstructionBase.setproperties(out, NamedTuple{(k,)}((endemic,)))
end

