
const lc_categories = (native=1, cleared=2, abandoned=3, urban=4, forestry=5, water=6)
const category_names = NamedTuple{keys(lc_categories)}(keys(lc_categories))
const island_keys = (; mus=:mus, reu=:reu, rod=:rod)
