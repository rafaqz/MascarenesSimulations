include("species_common.jl")

# 600 feral pigs in mauritius 2009 !! Lubisi et al

# Set up and run simulations
k = :reu
k = :rod
k = :mus
# pred_keys = (:cat, :black_rat, :norway_rat, :mouse, :pig, :macaque)
# pred_keys = (:cat, :black_rat, :norway_rat, :pig, :mouse)
pred_keys = (:cat, :black_rat, :norway_rat)

include("makie.jl")
include("species_rules.jl")

EndemicNVs = begin
    map(island_endemic_tables) do endemic_table
        ek = Tuple(Symbol.(replace.(endemic_table.Species, Ref(' ' => '_'))))
        NamedVector{ek,length(ek)}
    end
end
(; ruleset, rules, pred_ruleset, endemic_ruleset, islands, pred_response) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, pred_keys, first_year, last_year, extant_extension,
    pred_pops_aux = map(_ -> nothing, dems),
    # pred_funcs,
    island_recouperation_rates = map(EndemicNVs) do EndemicNV
        Float32.(ones(EndemicNV) .* 0.05)
    end,
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];
k
sim!(output, ruleset; proc=SingleCPU(), printframe=true, tspan=1600:1610);

# mkoutput = mk_pred(init, pred_ruleset; landcover=lc_all[k], output_kw...)
mkoutput = mk(init, ruleset; landcover=lc_all[k], output_kw..., ncolumns=4)
# display(mkoutput)


# Optimisation
pred_pops_aux = map(islands) do island
    (; pred_output, init) = island
    @time sim!(pred_output, pred_ruleset; proc=SingleCPU(), printframe=true);
    A = cat(pred_output...; dims=3)
    DimArray(A, (dims(init.pred_pop)..., dims(pred_output)...))
end
sum(getproperty.(pred_pops_aux.mus[Ti=At(2009)], :pig))
sum(getproperty.(pred_pops_aux.mus[Ti=At(2009)], :macaque))
sum(getproperty.(pred_pops_aux.rod[Ti=At(2009)], :cat))

# Store so we don't have to run the above
jldsave("sym_setup2_$aggfactor.jld2";
    auxs, pred_pops_aux, pred_response
);


(; ruleset, rules, pred_ruleset, endemic_ruleset, islands) = def_syms(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor; 
    replicates=nothing, first_year, last_year, extant_extension,
    pred_pops_aux=map(_ -> nothing, dems),
    pred_keys,
);
(; output, endemic_output, pred_output, init, output_kw) = islands[k];
# @time sim!(endemic_output, endemic_ruleset; proc=SingleCPU(), printframe=true);
# @time sim!(output, ruleset; proc=SingleCPU(), printframe=true);
#
