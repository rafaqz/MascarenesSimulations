# Force remove master environment before activate
# to make sure it works
# deleteat!(Base.LOAD_PATH, 2:3)

using MascarenesSimulations
using DynamicGrids
using GLMakie

# Settings

aggfactor = 16
first_year = 1550
last_year = 2018
extant_extension = 0

# Choose predator subset
pred_keys = (:cat, :black_rat, :norway_rat) # (:cat, :black_rat, :norway_rat, :mouse, :pig, :macaque)

# Choose an island
k = :reu
k = :rod
k = :mus

# Landcover paths
datadir = "/home/raf/PhD/Mascarenes/Data/Generated"

landcover_paths = (
    mus="$datadir/lc_predictions_mus.nc",
    reu="$datadir/lc_predictions_reu.nc",
    rod="$datadir/lc_predictions_rod.nc",
)

# Load data
(; pred_df, introductions_df, island_names, island_endemic_names, island_tables, island_endemic_tables) = load_tables()
(; borders, masks, elevation, dems) = load_rasters()
auxs = load_aux(; masks, landcover_paths, aggfactor, last_year)

# Run full invasive/endemic simulations

# Define rules and outputs
(; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = define_simulations(
    pred_df, 
    introductions_df, 
    island_endemic_tables, 
    auxs, 
    aggfactor; 
    replicates=nothing, 
    first_year, 
    last_year, 
    extant_extension,
    pred_keys,
    pred_pops_aux=map(_ -> nothing, dems),
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];

# Run
@time sim!(output, ruleset; proc=SingleCPU(), printframe=true);

# Makie visual simulations
lc_graphic = graphic_landcover(auxs)
mkoutput = makie_sim(init, ruleset; landcover=lc_graphic[k], output_kw..., ncolumns=4)


# If you need to debug performance
# you should get over 20 frames a second for all pred + endemic rules

# using ProfileView
# @profview 1 + 1 # warmup
# sim!(output, ruleset; proc=SingleCPU(), printframe=true, tspan=1550:1551);
# Then profile a single frame
# @profview sim!(output, ruleset; proc=SingleCPU(), printframe=true, tspan=1550:1551);




# Endemic-only sims

# Simulate and store invasive predator population dynamics
# pred_pop_jld = "../cache/pred_pops_$aggfactor.jld2"
# if isfile(pred_pop_jld)
#     _jld = jldopen(pred_pop_jld, "r")
#     pred_pops_aux = _jld["pred_pops_aux"];
#     close(_jld)
# else
#     (; ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response) = define_simulations(
#         pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
#         replicates=nothing, pred_keys, first_year, last_year, extant_extension,
#         pred_pops_aux = map(_ -> nothing, dems),
#     );
#     pred_pops_aux = map(islands) do island
#         (; pred_output, init) = island
#         @time sim!(pred_output, pred_ruleset; proc=SingleCPU(), printframe=true);
#         A = cat(pred_output...; dims=3)
#         DimArray(A, (dims(init.pred_pop)..., dims(pred_output)...))
#     end

#     # Store so we don't have to run the above
#     jldsave(pred_pop_jld; pred_pops_aux, pred_response);
# end
# sum(getproperty.(pred_pops_aux.rod[Ti=At(2009)], :cat))

# (; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = define_simulations(
#     pred_df, 
#     introductions_df, 
#     island_endemic_tables, 
#     auxs, 
#     aggfactor; 
#     replicates=nothing, 
#     first_year, 
#     last_year, 
#     extant_extension,
#     pred_keys,
#     pred_pops_aux,
# );

# mkoutput = mk_pred(init, pred_ruleset; landcover=lc_all[k], output_kw...)

