# Struct rules

# --- Model constants ---------------------------------------------------------
const MIN_CARRYCAP_FACTOR         = 1f-10   # avoids division-by-zero in logistic growth
const CLEARING_NATIVE_THRESHOLD   = 0.2f0   # native cover below which a cell is cleared
const BASE_STOCHASTIC_EXTIRPATION = 0.005f0 # per-year background extinction risk (before aggfactor scaling)
const PRED_POP_PRESSURE_SCALE     = 32f0    # scales pred relative pop before extirpation calculation
const PRED_SUSCEPT_NORM           = 32 * 8^2 # normalization denominator in predator_susceptibility
const PROPAGULE_SPREAD_SCALE      = 40f0    # overall scaling of dispersal propagule counts
const PROPAGULE_SKEW_EXPONENT     = 3       # rand exponent giving sparse (right-skewed) propagule events
# -----------------------------------------------------------------------------

struct InteractiveCarryCap{R,W,CC,CS,I} <: Dispersal.GrowthRule{R,W}
    carrycap::CC
    carrycap_scaling::CS
    inputs::I
end
function InteractiveCarryCap{R,W}(; carrycap, carrycap_scaling, inputs=(;),) where {R,W}
    InteractiveCarryCap{R,W}(carrycap, carrycap_scaling, inputs)
end

@inline Base.@assume_effects :foldable function DynamicGrids.applyrule(
    data, rule::InteractiveCarryCap,
    populations::NamedVector{Keys}, I,
) where Keys
    local_inputs = get(data, rule.inputs, I)
    return calc_carrycaps(local_inputs, populations, rule.carrycap, rule.carrycap_scaling)
end

function calc_carrycaps(local_inputs, populations, carrycap, carrycap_scaling)
    relative_pop = NamedTuple(populations ./ carrycap)
    params = merge(relative_pop, NamedTuple(local_inputs))
    scaling = map(carrycap_scaling) do val_f
        f = DynamicGrids._unwrap(val_f)
        (oneunit(eltype(populations)) + f(params))
    end |> NamedVector
    absolute_min_carrycap = oneunit(eltype(carrycap)) .* MIN_CARRYCAP_FACTOR
    new_carrycaps = max.(absolute_min_carrycap, carrycap .* scaling)
    return new_carrycaps
end

struct ExtirpationRisks{R,W,F,S,T,PR,PS,PP,PC,E,SE,RR} <: DynamicGrids.NeighborhoodRule{R,W}
    f::F
    stencil::S
    traits::T
    pred_response::PR
    pred_suscept::PS
    pred_pop::PP
    pred_carrycap::PC
    pred_effect::E
    stochastic_extirpation::SE
    recuperation_rates::RR
end
function ExtirpationRisks{R,W}(; f, stencil, traits, pred_response, pred_suscept, pred_pop, pred_carrycap, pred_effect, stochastic_extirpation, recuperation_rates) where {R,W}
    ExtirpationRisks{R,W}(f, stencil, traits, pred_response, pred_suscept, pred_pop, pred_carrycap, pred_effect, stochastic_extirpation, recuperation_rates)
end

# Single-output version (presence only, no causes tracking). Intended for future
# use; currently superseded by the tuple version below.
@inline function DynamicGrids.applyrule(data, rule::ExtirpationRisks, endemic_presences, I)
    recuperation_rates = get(data, rule.recuperation_rates)

    pred_effect = if isnothing(rule.pred_effect)
        pred_relative_pop = get(data, rule.pred_pop, I) ./ rule.pred_carrycap
        pred_suscept = get(data, rule.pred_suscept)
        map(rule.f, predator_effect(pred_relative_pop, pred_suscept))
    else
        get(data, rule.pred_effect, I)
    end

    hood = stencil(rule)
    # Count neighbors to UInt8. Otherwise we get Int64 and use a lot of registers
    n_neighbors = foldl(hood; init=Base.reinterpret.(UInt8, zero(first(hood)))) do x, y
        Base.reinterpret(UInt8, x) + Base.reinterpret(UInt8, y)
    end

    updated_presence = map(endemic_presences, pred_effect, n_neighbors, recuperation_rates) do present, effect, n_neighbor, rr
        if present
            rand(typeof(effect)) > (effect * ((length(hood) / 2) / (n_neighbor + 1)) + rule.stochastic_extirpation)
        elseif n_neighbor > 0
            rand(typeof(effect)) < (n_neighbor * rr / length(hood))
        else
            false
        end
    end
    return updated_presence
end

# Causes-tracking version: reads and writes both endemic_presence and causes grids.
@inline function DynamicGrids.applyrule(data, rule::ExtirpationRisks{Grids,Grids}, (endemic_presences, causes), I) where Grids<:Tuple
    recuperation_rates = get(data, rule.recuperation_rates)

    hood = stencil(rule)
    # Count neighbors to UInt8. Otherwise we get Int64 and use a lot of registers
    n_neighbors = foldl(hood; init=Base.reinterpret.(UInt8, zero(first(hood)))) do x, y
        Base.reinterpret(UInt8, x) + Base.reinterpret(UInt8, y)
    end

    pred_relative_pop = get(data, rule.pred_pop, I) ./ rule.pred_carrycap .* PRED_POP_PRESSURE_SCALE
    pred_suscept = get(data, rule.pred_suscept)
    pred_effects = predator_effects(pred_relative_pop, pred_suscept)

    results = map(endemic_presences, pred_effects, n_neighbors, recuperation_rates, causes) do present, effects, n_neighbor, rr, cs
        if present
            effect = rule.f(sum(effects))
            extirpation_probability = effect * length(hood) / (n_neighbor + 1) + rule.stochastic_extirpation
            still_present = rand(typeof(first(pred_relative_pop))) > extirpation_probability
            updated_causes = if still_present
                cs
            else
                effects
            end
        elseif n_neighbor > 0
            still_present = rand(typeof(first(pred_relative_pop))) < (n_neighbor * rr / length(hood))
            updated_causes = if still_present
                zero(cs)
            else
                cs
            end
        else
            still_present = false
            updated_causes = cs
        end
        still_present, updated_causes
    end
    updated_presence = map(first, results)
    updated_causes = map(last, results)

    return updated_presence, updated_causes
end

@inline function DynamicGrids.modifyrule(rule::ExtirpationRisks, data::AbstractSimData)
    if isnothing(rule.pred_effect)
        pred_response = get(data, rule.pred_response)
        traits = get(data, rule.traits)
        @set! rule.pred_suscept = predator_susceptibility(pred_response, traits)
    end
    @set! rule.traits = nothing        # Simplify for GPU argument size
    @set rule.pred_response = nothing  # Simplify for GPU argument size
end

Base.@assume_effects :foldable function predator_effects(pred_pop, pred_suscept::NamedVector{<:Any,<:Any,<:NamedVector{K}}) where K
    is = ntuple(identity, length(K))
    map(is) do i
        map(pred_suscept, pred_pop) do ps, pp
            Float32(ps[i] * pp)
        end
    end |> NamedVector{K,length(K)}
end
Base.@assume_effects :foldable function predator_effect(pred_pop, pred_suscept)
    mapreduce(+, pred_suscept, pred_pop) do ps, pp
         map(xs -> map(Float32 ∘ *, xs, pp), ps)
    end
end

function predator_susceptibility(pred_response, traits)
    mapreduce(+, pred_response, traits) do pr, t
        map(NamedVector(pr)) do p
            NamedVector(t) .* p
        end |> NamedVector
    end ./ PRED_SUSCEPT_NORM
end

# ---------------------------------------------------------------------------
# build_rules — island-independent rules and rulesets
# ---------------------------------------------------------------------------

function build_rules(pred_df, aggfactor;
    aggscale = aggfactor^2,
    pred_keys = (:cat, :black_rat, :norway_rat, :mouse, :pig, :wolf_snake, :macaque),
    pred_response = predator_response_params(pred_keys),
    # Prey-mass data (Pearre & Maaas 1998 for cats; estimates for the rest)
    mean_prey_mass = (;
        cat =        (41.0, 51.0),
        black_rat =  (10.0, 10.0),
        norway_rat = (8.0,  8.0),
        mouse =      (3.0,  5.0),
        pig =        (100.0, 100.0),
        wolf_snake = (9.0,  7.0),   # Estimated from Fritts 1993
        macaque =    (40.0, 40.0),
    )[pred_keys],
    # Carrying capacities from Harper & Bunbury 2015 and literature estimates
    carrycap = Float32.(NV(;
        cat =        0.01,
        black_rat =  30.0,   # may be up to 100/ha — Harper & Bunbury 2015
        norway_rat = 15.0,
        mouse =      52.0,
        pig =        0.02,   # tuned to ~600 in Mauritius 2009
        wolf_snake = 10.0,
        macaque =    0.5,
    ) .* aggscale)[pred_keys],
    spread_rate = NV(;
        cat =        20.0f0,
        black_rat =   1.0f0,
        norway_rat =  1.0f0,
        mouse =       0.5f0,
        pig =        100.0f0,
        wolf_snake =  1.0f0,
        macaque =     5.0f0,
    )[pred_keys],
    #= Predator carrying-capacity scaling as function of neighbour densities and landcover.
    Assumptions (Refs: Smucker et al 2000 — Hawaii cats, rats, mice):
      1. Cats suppress rodents; black rats more than norway rats (size selection >250 g)
      2. Cats live near people — density ~2 orders of magnitude higher in urban cells
      3. Black rats live anywhere but are outcompeted by norway rats in cities / near cats
      4. Norway rats dominate urban/coastal; lose in forest (bad climbers, less efficient)
      5. Mice outcompete rats in farmland, coexist in urban areas
      6. Pigs are largely independent of other predators
    =#
    pred_funcs = (;
        cat =        p -> 2.0f0p.black_rat + 0.5f0p.norway_rat + 10f0p.urban + 2f0p.cleared,
        black_rat  = p -> -0.2f0p.cat - 0.1f0p.norway_rat + 0.5f0p.native + 0.3f0p.abandoned + 0.3f0p.forestry + 1p.urban,
        norway_rat = p -> -0.1f0p.cat - 0.1f0p.black_rat + 1.5f0p.urban - 0.2f0p.native,
        mouse =      p -> -0.3f0p.cat - 0.2f0p.black_rat - 0.2f0p.norway_rat + 0.8f0p.cleared + 1.5f0p.urban,
    )[pred_keys],
)
    pred_df    = filter(r -> Symbol(r.name) in pred_keys, pred_df)
    PredNV     = NamedVector{pred_keys,length(pred_keys)}
    pred_rmax  = Float32.(PredNV(pred_df.rmax))
    pred_indices = PredNV(ntuple(identity, length(pred_keys)))
    pred_init_nvs = map(pred_indices) do i
        x = zeros(Float16, length(pred_keys))
        x[i] = 50 * aggfactor
        PredNV(x)
    end

    moore  = Moore{3}()
    kernel = DispersalKernel(
        stencil=moore,
        formulation=ExponentialKernel(Param(1.0f0, bounds=(0.0f0, 2.0f0))),
        cellsize=1.0f0,
    )
    pred_kernels = map(spread_rate) do s
        DispersalKernel(stencil=moore, formulation=ExponentialKernel(s), cellsize=1.0f0)
    end

    pred_carrycap_rule = InteractiveCarryCap{:pred_pop,:pred_carrycap}(;
        carrycap,
        carrycap_scaling=map(Val, pred_funcs),
        inputs=Aux{:landcover}(),
    )
    pred_growth_rule = LogisticGrowth{:pred_pop}(;
        rate=pred_rmax,
        carrycap=DynamicGrids.Grid{:pred_carrycap}(),
        timestep=1,
        nsteps_type=Float32,
    )
    introduction_rule = let introductions_aux=Aux{:introductions}()
        SetGrid{:pred_pop}() do data, pred_pop, t
            D = dims(DG.init(data).pred_pop)
            current_year = currenttime(data)
            intros = get(data, introductions_aux)
            foreach(intros) do intro
                if intro.year == current_year
                    p = intro.geometry
                    I = DimensionalData.dims2indices(D, (X(Contains(p.X)), Y(Contains(p.Y))))[1:2]
                    pred_pop[I..., :] .= view(pred_pop, I..., :) .+ (intro.init,)
                end
            end
        end
    end
    clearing_rule = let landcover=Aux{:landcover}()
        Cell{:endemic_presence}() do data, presences, I
            lc = get(data, landcover, I)
            presences .& (lc.native > CLEARING_NATIVE_THRESHOLD)
        end
    end
    pred_spread_rule = let demaux=Aux{:dem}(), aggfactor=aggfactor, pred_kernels=pred_kernels, carrycap=carrycap
        SetNeighbors{Tuple{:pred_pop,:pred_carrycap}}(pred_kernels[1]) do data, hood, (Ns, _), I
            Ns === zero(Ns) && return nothing
            dem = get(data, demaux)
            reps = DynamicGrids.replicates(data)
            carrycap_nbrs = if isnothing(reps)
                DG.neighbors(DG.grids(data).pred_carrycap, I)
            else
                DG.neighbors(DG.grids(data).pred_carrycap, (I..., reps))
            end
            elev_nbrs = DG.neighbors(dem, I)
            @inbounds elev_center = dem[I...]

            sum = zero(Ns)
            cellsize = 100 * aggfactor  # cell size in meters at base 100m resolution

            # Randomise hood starting position to avoid directional artifacts in output
            start = rand(0:length(hood)-1)
            @inbounds for ix in eachindex(hood)
                i = start + ix
                if i > length(hood)
                    i = i - length(hood)
                end
                sp = carrycap_nbrs[i] ./ carrycap
                any(map(isnan, sp)) && error("sp is NaN")
                Ih = DG.indices(hood, I)[i]
                ks = getindex.(DynamicGrids.kernel.(pred_kernels), i)
                d = DG.distances(hood)[i]
                e = elev_nbrs[i]
                propagules = trunc.(Float32.(Ns .* sp .* rand(typeof(Ns)) .^ PROPAGULE_SKEW_EXPONENT .* ks .* PROPAGULE_SPREAD_SCALE))
                sum1 = sum + propagules
                if any(sum1 .> Ns)
                    propagules = min.(sum1, Ns) .- sum
                    sum = sum + propagules
                else
                    sum = sum1
                end
                add!(data[:pred_pop], propagules, Ih...)
            end
            @inbounds sub!(data[:pred_pop], sum, I...)
            return nothing
        end
    end
    risks_rule = ExtirpationRisks{Tuple{:endemic_presence,:causes}}(;
        f=tanh,
        stencil=Moore(1),
        traits=Aux{:endemic_traits}(),
        pred_response,
        pred_suscept=nothing,
        pred_effect=nothing,
        pred_pop=Grid{:pred_pop}(),
        pred_carrycap=carrycap,
        stochastic_extirpation=BASE_STOCHASTIC_EXTIRPATION/aggfactor,
        recuperation_rates=Aux{:recuperation_rates}(),
    )
    aux_pred_risks_rule = ExtirpationRisks{:endemic_presence,:causes}(;
        f=tanh,
        stencil=Moore(1),
        traits=Aux{:endemic_traits}(),
        pred_response,
        pred_suscept=nothing,
        pred_effect=Aux{:pred_effect}(),
        pred_pop=nothing,
        pred_carrycap=carrycap,
        stochastic_extirpation=BASE_STOCHASTIC_EXTIRPATION/aggfactor,
        recuperation_rates=Aux{:recuperation_rates}(),
    )

    pred_ruleset = Ruleset(
        pred_carrycap_rule, introduction_rule, pred_spread_rule, pred_growth_rule;
        boundary=Remove()
    )
    endemic_ruleset = Ruleset(Chain(aux_pred_risks_rule, clearing_rule); boundary=Remove())
    ruleset = Ruleset(
        DynamicGrids.rules(pred_ruleset)..., risks_rule, clearing_rule;
        boundary=Remove()
    )
    rules = (;
        introduction_rule, pred_carrycap_rule, pred_spread_rule, pred_growth_rule,
        risks_rule, aux_pred_risks_rule, clearing_rule,
    )

    return (;
        ruleset, pred_ruleset, endemic_ruleset, rules,
        pred_rmax, carrycap, PredNV, pred_keys, pred_init_nvs, pred_kernels,
        pred_response, mean_prey_mass, aggfactor, aggscale, moore, kernel, risks_rule,
    )
end

# ---------------------------------------------------------------------------
# build_island — per-island init arrays, outputs, and aux data
# ---------------------------------------------------------------------------

function build_island(key, endemic_table, aux, introductions_df, shared;
    first_year, last_year, extant_extension, replicates,
    pred_pops_aux = nothing,
    EndemicNV = begin
        ek = Tuple(Symbol.(replace.(endemic_table.Species, Ref(' ' => '_'))))
        NamedVector{ek,length(ek)}
    end,
    # Extinction date defaults: start year of the recorded range, or last_year+extension
    extinction_dates = map(endemic_table[!, Symbol(key, :_extinct)]) do x
        ismissing(x) ? last_year + extant_extension : parse(Int, first(split(string(x), ':')))
    end |> EndemicNV,
    mass_response = begin
        endemic_mass = EndemicNV(endemic_table.Mass)
        map(shared.mean_prey_mass) do (mean, std)
            dist = Distributions.Normal(mean, 2std) # Doubled to account for distribution skew
            scalar = 1 / pdf(dist, mean)
            map(endemic_mass) do pm; pdf(dist, pm) * scalar end
        end
    end,
    recuperation_rates = Float32.(ones(EndemicNV) .* 1.0), # TODO: derive from literature
)
    (; pred_rmax, carrycap, PredNV, pred_keys, pred_init_nvs, pred_kernels,
       pred_response, aggfactor, moore, kernel, risks_rule) = shared

    stencil_mask = StencilArray(aux.mask, kernel; padding=Halo{:out}())
    stencil_dem  = StencilArray(aux.dem, moore; padding=Halo{:out}())

    traits = endemic_traits(endemic_table, EndemicNV)

    island_df = filter(r -> r.island == string(key) && Symbol(r.species) in pred_keys, introductions_df)
    display(island_df)
    island_introductions = map(eachrow(island_df)) do r
        (; year=r.year, geometry=(X=r.lon, Y=r.lat), init=pred_init_nvs[Symbol(r.species)])
    end

    pred_pop      = map(_ -> map(_ -> 0.0f0, pred_rmax), aux.mask)
    pred_carrycap = map(_ -> carrycap, aux.mask)
    # Every endemic species is present everywhere initially
    endemic_presence = map(aux.mask) do m
        map(_ -> m, traits.ismammal) # ismammal used only as a shape template
    end
    causes = map(aux.mask) do m
        map(_ -> zero(pred_rmax), traits.ismammal)
    end

    pred_effect = if isnothing(pred_pops_aux)
        nothing
    else
        pr = ModelParameters.stripparams(pred_response)
        pred_suscept = predator_susceptibility(pr, traits)
        generate_predator_effect(risks_rule.f, pred_pops_aux, pred_suscept)
    end

    init         = (; pred_pop, pred_carrycap, endemic_presence, causes)
    pred_init    = (; pred_pop, pred_carrycap)
    endemic_init = (; endemic_presence)
    tspan        = first_year:last_year

    output_kw = (;
        aux=(;
            introductions = island_introductions,
            dem           = stencil_dem,
            recuperation_rates,
            pred_pop      = pred_pops_aux,
            endemic_traits = traits,
            pred_effect,
            landcover     = aux.lc,
        ),
        mask = stencil_mask,
        replicates,
        tspan,
    )

    output = ResultOutput(init; output_kw...)
    pred_output = if isnothing(replicates)
        TransformedOutput(pred_init; output_kw...) do f
            Array(f.pred_pop)
        end
    else
        ResultOutput(pred_init; output_kw...)
    end
    endemic_output = ResultOutput(endemic_init; output_kw...)

    return (; key, init, endemic_init, pred_init, output, endemic_output, pred_output,
              output_kw, aux, mass_response, extinction_dates)
end

# ---------------------------------------------------------------------------
# define_simulations — thin coordinator: calls build_rules then build_island
# ---------------------------------------------------------------------------

function define_simulations(
    pred_df, introductions_df, island_endemic_tables, auxs, aggfactor;
    aggscale     = aggfactor^2,
    replicates   = nothing,
    pred_pops_aux,
    island_keys  = NamedTuple{keys(island_endemic_tables)}(keys(island_endemic_tables)),
    first_year, last_year, extant_extension,
    pred_keys    = (:cat, :black_rat, :norway_rat, :mouse, :pig, :wolf_snake, :macaque),
    pred_response = predator_response_params(pred_keys),
    mean_prey_mass = (;
        cat =        (41.0, 51.0),
        black_rat =  (10.0, 10.0),
        norway_rat = (8.0,  8.0),
        mouse =      (3.0,  5.0),
        pig =        (100.0, 100.0),
        wolf_snake = (9.0,  7.0),
        macaque =    (40.0, 40.0),
    )[pred_keys],
    carrycap = Float32.(NV(;
        cat =        0.01,
        black_rat =  30.0,
        norway_rat = 15.0,
        mouse =      52.0,
        pig =        0.02,
        wolf_snake = 10.0,
        macaque =    0.5,
    ) .* aggscale)[pred_keys],
    spread_rate = NV(;
        cat =        20.0f0,
        black_rat =   1.0f0,
        norway_rat =  1.0f0,
        mouse =       0.5f0,
        pig =        100.0f0,
        wolf_snake =  1.0f0,
        macaque =     5.0f0,
    )[pred_keys],
    pred_funcs = (;
        cat =        p -> 2.0f0p.black_rat + 0.5f0p.norway_rat + 10f0p.urban + 2f0p.cleared,
        black_rat  = p -> -0.2f0p.cat - 0.1f0p.norway_rat + 0.5f0p.native + 0.3f0p.abandoned + 0.3f0p.forestry + 1p.urban,
        norway_rat = p -> -0.1f0p.cat - 0.1f0p.black_rat + 1.5f0p.urban - 0.2f0p.native,
        mouse =      p -> -0.3f0p.cat - 0.2f0p.black_rat - 0.2f0p.norway_rat + 0.8f0p.cleared + 1.5f0p.urban,
    )[pred_keys],
)
    shared = build_rules(pred_df, aggfactor;
        aggscale, pred_keys, pred_response, mean_prey_mass, carrycap, spread_rate, pred_funcs)

    islands = map(island_keys, island_endemic_tables, auxs, pred_pops_aux) do key, table, aux, ppa
        build_island(key, table, aux, introductions_df, shared;
            first_year, last_year, extant_extension, replicates, pred_pops_aux=ppa)
    end

    return (; shared.ruleset, shared.rules, shared.pred_ruleset, shared.endemic_ruleset,
              islands, shared.pred_response)
end
