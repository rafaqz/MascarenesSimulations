
const lc_categories = (native=1, cleared=2, abandoned=3, urban=4, forestry=5, water=6)
const category_names = NamedTuple{keys(lc_categories)}(keys(lc_categories))
const island_keys = (; mus=:mus, reu=:reu, rod=:rod)

# Parse a CSV year-range string ("1600" or "1600:1700") into a UnitRange.
# Defined here so it is available to both tables.jl and functions.jl.
function parse_year_range(s::AbstractString)
    parts = split(s, ':')
    parse(Int, parts[1]):parse(Int, parts[end])
end
parse_year_range(::Missing) = missing

# Normalise island-presence column values to Bool.
# The source CSV uses inconsistent encodings (p, x, y, TRUE, xipi, …) for presence;
# any non-missing value means present.
is_present(::Missing) = missing
is_present(_)         = true
