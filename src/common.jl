workdir = "/home/raf/PhD/Mascarenes"
datadir = joinpath(workdir, "Data")
outputdir = joinpath(datadir, "Generated")
distancedir = joinpath(outputdir, "Distances")

lc_categories = (native=1, cleared=2, abandoned=3, urban=4, forestry=5, water=6)
category_names = NamedTuple{keys(lc_categories)}(keys(lc_categories))
island_keys = (; mus=:mus, reu=:reu, rod=:rod)

fix_order(A) = reorder(A, ForwardOrdered)

