deleteat!(Base.LOAD_PATH, 2:3)

include("species_common.jl")
include("species_tables.jl")
include("species_rules.jl")
include("makie.jl")

# pred_keys = (:cat, :black_rat, :norway_rat, :mouse, :pig, :macaque)
pred_keys = (:cat, :black_rat, :norway_rat)

EndemicNVs = map(island_endemic_tables) do endemic_table
    ek = Tuple(Symbol.(replace.(endemic_table.Species, Ref(' ' => '_'))))
    NamedVector{ek,length(ek)}
end

# Simulate and store invasive predator population dynamics
pred_pop_jld = "pred_pops_$aggfactor.jld2"
if isfile(pred_pop_jld)
    _jld = jldopen(pred_pop_jld, "r")
    pred_pops_aux = _jld["pred_pops_aux"];
    close(_jld)
else
    (; ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response) = def_syms(
        pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
        replicates=nothing, pred_keys, first_year, last_year, extant_extension,
        pred_pops_aux = map(_ -> nothing, dems),
    );
    pred_pops_aux = map(islands) do island
        (; pred_output, init) = island
        @time sim!(pred_output, pred_ruleset; proc=SingleCPU(), printframe=true);
        A = cat(pred_output...; dims=3)
        DimArray(A, (dims(init.pred_pop)..., dims(pred_output)...))
    end

    # Store so we don't have to run the above
    jldsave(pred_pop_jld; pred_pops_aux, pred_response);
end
# sum(getproperty.(pred_pops_aux.rod[Ti=At(2009)], :cat))


# Run full invasive/endemic simulations

# Choose and island
k = :reu
k = :rod
k = :mus

# Defin rules and outputs
(; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, first_year, last_year, extant_extension,
    pred_pops_aux=map(_ -> nothing, dems),
    pred_keys,
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];

@time sim!(output, ruleset; proc=SingleCPU(), printframe=true);

# TODO debug performance
# using ProfileView
# @profview sim!(output, ruleset; proc=SingleCPU(), printframe=true, tspan=1550:1552);

mkoutput = mk(init, ruleset; landcover=lc_all[k], output_kw..., ncolumns=4)

# mkoutput = mk_pred(init, pred_ruleset; landcover=lc_all[k], output_kw...)

# @time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);
# @time sim!(output, ruleset; proc=SingleCPU(), printframe=true);

