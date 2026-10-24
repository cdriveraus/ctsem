using ComponentArrays

export EKFParameters
"""
    UNSET_PARAMETER

Fill value for a parameter slot that nothing has written yet.

The full parameter vector is materialised at *mutable* and at *fixed*
positions only, so a model-matrix cell owned solely by a state-dependent
transform is covered by neither: it holds whatever the buffer was filled
with until its transform group runs. Reading it before then is always a
bug, and the fill value decides whether that bug is visible.

Zero is the worst possible choice, because it is what a well-behaved
model matrix is mostly made of: it annihilates the products it enters and
reads back as a plausible fixed value. The row-1 defect fixed in
`kalman_filters.jl` survived precisely that way -- and one attempt to
reproduce it saw nothing, because the test model happened to fix the
initial latent mean at zero, so multiplying it by the unwritten cell hid
the fault a second time.

This follows the convention ctsem already uses for a missing predictor: a
number no real parameter can take, finite so that ForwardDiff partials
stay clean rather than turning into NaNs that travel silently, and large
enough that any likelihood it reaches is absurd on sight.
"""
const UNSET_PARAMETER = 99999.0

# The templates ctModelWriter writes by default, as Julia source with index 1
# standing for any: DRIFT's diagonal (continuous and discrete time), DIFFUSION,
# standard deviations (MANIFESTVAR, T0VAR), the population block, identity,
# intercepts and means, loadings, thresholds. Every spec carries a group for
# each of these, empty when unused -- see `_group_regular_transforms`. Listed
# from the parameter tables of default models of every kind; the R suite fails
# when a default model emits a single-index template that is not here.
const _DEFAULT_TEMPLATE_SOURCES = (
    "-(1e-06 + 2 * log1p_exp(-(2 * param[1])))",
    "1/(1 + exp(-(2 * param[1])))",
    "10 * log1p_exp(2 * param[1]) + 1e-10",
    "1e-10 + 5 * log1p_exp(2 * param[1])",
    "1e-10 + 1 * log1p_exp(2 * param[1] - 1)",
    "param[1]",
    "10 * param[1]",
    "0.5 + 5 * param[1]",
    "2 * log1p_exp(2 * param[1])",
)
const _DEFAULT_TEMPLATE_TYPES = Any[]
const _DEFAULT_TEMPLATE_LOCK = ReentrantLock()

"""The closure types of the default templates, built once per session."""
function _default_template_types()
    lock(_DEFAULT_TEMPLATE_LOCK) do
        if isempty(_DEFAULT_TEMPLATE_TYPES)
            for src in _DEFAULT_TEMPLATE_SOURCES
                push!(_DEFAULT_TEMPLATE_TYPES, typeof(_regular_transform_closure(src)))
            end
        end
        return _DEFAULT_TEMPLATE_TYPES
    end
end

"""The default templates' cache keys."""
ctsem_default_template_keys() = [_regular_transform_key(s) for s in _DEFAULT_TEMPLATE_SOURCES]

export ctsem_nondefault_templates
"""
    ctsem_nondefault_templates(transforms)

The single-index transforms among `transforms` whose template is not a default
one, joined by `" | "`, or `""` when there are none. The R suite asks this of
default models, so that a template ctModelWriter starts writing by default
fails a test instead of quietly making models of equal dimensions different
types again. A string, not a vector: an empty vector does not cross the bridge.
"""
function ctsem_nondefault_templates(transforms::AbstractVector)
    defaults = Set(ctsem_default_template_keys())
    out = String[]
    for t in transforms
        s = String(t)
        isempty(s) && continue
        _regular_transform_index(s) === nothing && continue
        _regular_transform_key(s) in defaults || push!(out, s)
    end
    return join(unique(out), " | ")
end

"""
    _group_regular_transforms(transforms, mutables)

Group the regular transforms by template type.

Returns `(groups, cells, tf)`: a tuple with one `Vector` of same-typed closures
per distinct template, and for each group the `all_params` position each member
writes and the position it held in the original order.

Groups are ordered by the template's type name so that two models using the same
templates produce the *same* tuple type whatever order their cells appear in --
which is the whole point, since that type is what `EKFParameters` is
parameterised by.

Every default template has a group whether or not the model uses it, so the
tuple's type depends on the templates only when a model brings one of its own.
Without that, two models of the same dimensions differing in which defaults they
use -- one frees a loading, the other does not -- were different types, and the
whole pipeline, Hessian duals included, compiled again for the second: about
40 s on dev2 for a one-latent model. An empty group costs a loop that runs no
times.
"""
function _group_regular_transforms(transforms::Vector{Any}, mutables::BitVector)
    cells = findall(mutables)
    length(cells) == length(transforms) || throw(DimensionMismatch(
        "$(length(transforms)) regular transforms for $(length(cells)) mutable " *
        "positions; they are produced one per mutable cell"))
    types = unique(typeof.(transforms))
    for T in _default_template_types()
        T in types || push!(types, T)
    end
    # Every other template is an expression (`_ctsem_regular_expr`), and that
    # group too is carried empty when unused.
    CTSEMExpr in types || push!(types, CTSEMExpr)
    sort!(types; by = string)
    groups = Any[]
    gcells = Vector{Int}[]
    gtf = Vector{Int}[]
    for T in types
        members = [k for k in eachindex(transforms) if typeof(transforms[k]) === T]
        push!(groups, T[transforms[k] for k in members])
        push!(gcells, Int[cells[k] for k in members])
        push!(gtf, Int[k for k in members])
    end
    return (Tuple(groups), gcells, gtf)
end

"""
    CTSEMHybrid

Random-effect deviations the subject filter integrates itself rather than the
Laplace layer: one per item, an own variance (`popunique`) on a cell that
enters a mean affinely -- a manifest intercept (MANIFESTMEANS, `kind` 1), an
initial state mean (T0MEANS, `kind` 2) or a continuous intercept (CINT,
`kind` 3) -- whose transform is affine. For fixed deviations a linear Gaussian
filter's covariances and gains do not depend on them -- they only shift
means -- so the log likelihood is exactly quadratic in them and their
Gaussian prior integrates out in closed form (de Jong's augmented filter).
The filter carries `X`, the state mean's dependence on a unit of each cell's
value: a manifest intercept's column starts at zero and shifts its item's
observation directly (`E`), an initial mean's starts at `e_l` and propagates
by the discrete drift, a continuous intercept's gains the discrete
intercept's Jacobian in it at every prediction substep (`_hybrid_predict!`).
It accumulates `S = sum V' S_t^-1 V` and `s = sum V' S_t^-1 y_t` with
`V = Lambda X + E` (`_hybrid_update!`); the pass ends by adding
`-1/2 logdet(M) + 1/2 b' M^-1 b`, `M = I + G' S G`, `b = G' s`
(`_hybrid_marginal!`), where the cells' deviations are `G u`, `u ~ N(0, I)`.
Exact for a linear Gaussian model; the R side sends it for no other.

Two sources of `G`. Own variances (`block` false): each item its own sd on a
reduced-rank level (`popunique`), `G` diagonal, the items' sds times their
slopes, the level's shared loadings left to the Laplace layer. A full-rank
level's effects (`block` true): the level's covariance taken in the order
`perm` -- the effects the layer keeps first, the items last -- and factored,
`[[C_LL, 0], [C_FL, C_FF]]`; the layer keeps the leading columns as its
coordinates, so the items' conditional mean `C_FL u_L` is in the member's
values, and the filter integrates the trailing block, `G = slope .* C_FF`, the
conditional covariance's factor. In both, a subject that never observes a
manifest item has that row of `G` dropped and keeps every column.

Mutable and empty by default, because the Laplace objective that knows which
deviations are its own is built after the `EKFParameters` it wraps:
`ctsem_laplace_objective` fills it in place.

  * `kind[h]`: 1 a manifest intercept, 2 an initial state mean, 3 a
    continuous intercept (the discrete one in a discrete-time model), 4 a
    time dependent predictor's effect.
  * `target[h]`: the manifest item `h` shifts (`kind` 1), else the latent.
  * `pred[h]`: the time dependent predictor a `kind` 4 effect multiplies, else
    0.
  * `psi_index[h]`: the raw position of its own variance.
  * `effect[h]`: its index among the subject level's effects, whose `re_index`
    is the raw position the deviation shifts (reporting adds it there).
  * `sd_scale[h]`, `slope[h]`: its sdscale, and the intercept transform's slope
    in the raw value, so an own variance's deviation in manifest units is
    `slope * (log1p_exp(2 raw - 1) sd_scale + 1e-10) * u`, `u ~ N(0, 1)`.
  * `block`, and for it the level's construction (`_level_popcov`'s
    arguments: `lv_sd`, `lv_cor`, `lv_sdfix`, `lv_corfix`, `lv_scale`,
    `lv_code`), the order `perm` (the layer's effects, then the items in item
    order), and `lparams`, the raw positions of the level's free scales and
    correlations, which the closing term's cotangent reaches; `version`,
    counted up each time the fields are set, which the filter's cached factor
    is checked against (`_hybrid_block_cached!`).
"""
mutable struct CTSEMHybrid
    kind::Vector{Int}
    target::Vector{Int}
    psi_index::Vector{Int}
    effect::Vector{Int}
    sd_scale::Vector{Float64}
    slope::Vector{Float64}
    block::Bool
    lv_sd::Vector{Int}
    lv_cor::Vector{Int}
    lv_sdfix::Vector{Float64}
    lv_corfix::Vector{Float64}
    lv_scale::Vector{Float64}
    lv_code::Int
    perm::Vector{Int}
    lparams::Vector{Int}
    version::Int
    pred::Vector{Int}
end
CTSEMHybrid() = CTSEMHybrid(Int[], Int[], Int[], Int[], Float64[], Float64[], false,
    Int[], Int[], Float64[], Float64[], Float64[], 0, Int[], Int[], 0, Int[])
_hybrid_on(h::CTSEMHybrid) = !isempty(h.effect)
"""The model's hybrid item count, the capacity every per-item array is sized by."""
_hybrid_items(h::CTSEMHybrid) = length(h.effect)

"""
    CTSEMStationaryElements

Single initial conditions started at their stationary values: T0MEANS cells,
and blocks of T0VAR cells, that the model names 'stationary'
(`_ctsem_stationary_elements_mean!`, `_ctsem_stationary_elements_cov!`). Every
other initial condition is T0MEANS and T0VAR as written. `stationary = TRUE`
is the whole of both, and the two are not combined.

  * `inclosure` -- per latent, whether its stationary moments are needed: the
    marked latents and every latent that drives them through DRIFT. The set is
    shut under driving, so its dynamics are autonomous and its moments exist
    whatever the latents outside it do.
  * `means` -- the latents whose T0MEANS cell is marked.
  * `blocks` -- the latents whose T0VAR cells are marked, grouped: one latent
    whose diagonal alone is, or several whose cells among themselves all are.

Empty, and off, unless R sent the marks (`ekf_from_columns`).
"""
mutable struct CTSEMStationaryElements
    inclosure::Vector{Bool}
    means::Vector{Int}
    blocks::Vector{Vector{Int}}
end
CTSEMStationaryElements() = CTSEMStationaryElements(Bool[], Int[], Vector{Int}[])
_stationary_elements_on(s::CTSEMStationaryElements) = !isempty(s.inclosure)
# None, for callers that have no model's to pass. Never filled.
const _CTSEM_NO_STATIONARY_ELEMENTS = CTSEMStationaryElements()

"""
    CTSEMAnchors

Parameters whose population value is stated rather than estimated, with
random effects about it: a person-level latent level with mean zero, say.
Such a parameter has no coordinate in the raw vector the optimiser and the
samplers move. Its population value is a constant placed after the `nfree`
free parameters in each subject's parameter vector, which a cell reads like
any other raw parameter (`param[nfree + j]`), and which the random effects
shift like any other (`CTSEMLaplaceSpec` levels name it in `re_index`).

`_materialize_subject_values!` appends the constants to a vector of the free
parameters alone -- a population evaluation -- and takes a vector that
already carries them, a member's shifted parameters (`_laplace_extend`), as
it is. Empty, and off, unless R sent them (`ekf_from_columns`).
"""
mutable struct CTSEMAnchors
    nfree::Int
    values::Vector{Float64}
end
CTSEMAnchors() = CTSEMAnchors(0, Float64[])

"""
    EKFParameters

Container for parameter metadata used by the continuous-time EKF routines.

The fields describe which flattened parameter positions are mutable or fixed,
which transforms should be applied before filtering, and the `ComponentArray`
axis used to view the flattened parameter vector as named model matrices.
"""

# The regular transforms are held twice, and the second is the one the hot
# loops read.
#
# `regular_transforms` is a `Vector{Any}`: it keeps the original order for the
# handful of cold callers that want it, and holds no type information, so it
# does not put one entry per free cell into `EKFParameters`'s own type. That
# tuple was what made two models differing in template order or parameter count
# different types, and a new model type costs seconds of specialisation on its
# first evaluation while the same type costs nothing.
#
# `regular_groups` is a small tuple with one `Vector{T}` per *distinct template
# type*. Indexing one gives a concrete type, so the transform call is a static
# dispatch and nothing is boxed -- which erasing the tuple on its own would have
# lost, at 4% on every subject evaluation. And the tuple's type now depends only
# on the set of templates, which is a handful across every model ctsem writes,
# so essentially all models share one specialisation.
#
# `regular_group_cells` and `regular_group_tf` carry, per group member, where its
# result belongs in `all_params` and which position it held in the original
# order -- the latter because the reverse pass indexes its supports that way.
# ---- The parameter vector's layout, as data ----------------------------------
#
# A ComponentArrays axis spells every matrix's name, shape and offset into the
# type, so two models of different size are different types and every function
# that touches the parameter vector compiled once per model. `CTSEMLayout` keeps
# the same facts in fields, and `CTSEMPars` answers `pars.DRIFT` with a view of
# one concrete type whatever the model. The names the engine reads are fixed
# slots, resolved at compile time when the name is a literal (`_matrix_slot`);
# any other name a parameter table carries is found by search.

const _CTSEM_MATRIX_SLOTS = (:PARS, :T0MEANS, :LAMBDA, :DRIFT, :DIFFUSION,
    :MANIFESTVAR, :MANIFESTMEANS, :CINT, :T0VAR, :TDPREDEFFECT, :JAx, :Jy, :Jtd,
    :THRESHOLDS, :ASYMPTOTES, :RESPONSE, :RAWPOPVAR)

@eval @inline function _matrix_slot(s::Symbol)
    $([:(s === $(QuoteNode(n)) && return $i) for (i, n) in enumerate(_CTSEM_MATRIX_SLOTS)]...)
    return 0
end

struct CTSEMLayout
    names::Vector{Symbol}
    first::Vector{Int}
    nrows::Vector{Int}
    ncols::Vector{Int}
    # For each of `_CTSEM_MATRIX_SLOTS`, its index in `names`, or 0.
    slot::Vector{Int}
    total::Int
end

function CTSEMLayout(names, nrows, ncols, ranges)
    syms = Symbol[Symbol(n) for n in names]
    slot = zeros(Int, length(_CTSEM_MATRIX_SLOTS))
    for (i, s) in enumerate(syms)
        k = _matrix_slot(s)
        k > 0 && (slot[k] = i)
    end
    total = isempty(ranges) ? 0 : maximum(last, ranges)
    return CTSEMLayout(syms, Int[Base.first(r) for r in ranges], Int.(nrows), Int.(ncols),
        slot, total)
end

"""
The layout a ComponentArrays `Axis` describes (`retrieve_axes`): what a caller
building `EKFParameters` by hand -- the engine's tests -- has at hand.
"""
function CTSEMLayout(ax::ComponentArrays.Axis)
    names = Symbol[]
    nr = Int[]
    nc = Int[]
    ranges = UnitRange{Int}[]
    for (name, va) in pairs(ComponentArrays.indexmap(ax))
        r = ComponentArrays.viewindex(va)
        r = r isa Integer ? (r:r) : UnitRange{Int}(r)
        sz = va isa ComponentArrays.ViewAxis ? size(va.ax) : (length(r),)
        push!(names, name)
        push!(nr, sz[1])
        push!(nc, length(sz) > 1 ? sz[2] : 1)
        push!(ranges, r)
    end
    return CTSEMLayout(names, nr, nc, ranges)
end
_ctsem_layout(l::CTSEMLayout) = l
_ctsem_layout(ax::ComponentArrays.Axis) = CTSEMLayout(ax)

@inline function _layout_index(l::CTSEMLayout, s::Symbol)
    k = _matrix_slot(s)
    k > 0 && return @inbounds l.slot[k]
    i = findfirst(==(s), l.names)
    return i === nothing ? 0 : i
end

"""
    CTSEMPars(data, layout)

The flat parameter vector `data` with named, shaped access to its matrices:
`p.DRIFT` is a view. An `AbstractVector`, so indexing, broadcasting and
`copyto!` act on `data`; `getdata(p)` is `data` itself.
"""
struct CTSEMPars{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    layout::CTSEMLayout
end

Base.size(p::CTSEMPars) = size(getfield(p, :data))
Base.IndexStyle(::Type{<:CTSEMPars}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(p::CTSEMPars, i::Int) = getfield(p, :data)[i]
Base.@propagate_inbounds Base.setindex!(p::CTSEMPars, v, i::Int) = (getfield(p, :data)[i] = v)
ComponentArrays.getdata(p::CTSEMPars) = getfield(p, :data)
# Aliasing checks (`unalias`, broadcasting) ask what memory an array holds: the
# data's, not an identity hashed from the whole struct.
Base.dataids(p::CTSEMPars) = Base.dataids(getfield(p, :data))

@noinline _ctsem_no_matrix(s) = throw(ArgumentError("the model has no matrix $(s)"))

@inline function Base.getproperty(p::CTSEMPars, s::Symbol)
    (s === :data || s === :layout) && return getfield(p, s)
    l = getfield(p, :layout)
    i = _layout_index(l, s)
    i == 0 && _ctsem_no_matrix(s)
    @inbounds begin
        f = l.first[i]
        nr = l.nrows[i]
        nc = l.ncols[i]
    end
    return reshape(view(getfield(p, :data), f:(f + nr * nc - 1)), nr, nc)
end

Base.hasproperty(p::CTSEMPars, s::Symbol) = _layout_index(getfield(p, :layout), s) != 0
Base.propertynames(p::CTSEMPars, private::Bool=false) = Tuple(getfield(p, :layout).names)

struct EKFParameters{RG,PT,UT,TT,AX,FV}
    mutables::BitVector
    transform_indices::Vector{Bool}
    predict_transforms_indices::Vector{Bool}
    update_transforms_indices::Vector{Bool}
    td_transforms_indices::Vector{Bool}
    regular_transforms::Vector{Any}
    regular_groups::RG
    regular_group_cells::Vector{Vector{Int}}
    regular_group_tf::Vector{Vector{Int}}
    predict_transforms::PT
    update_transforms::UT
    td_transforms::TT
    parnumber::Vector{Int}
    parameter_axis::AX
    fixed_indices::Vector{Bool}
    fixed_values::FV
    ti_parameter_indices::Vector{Int}
    ti_predictor_indices::Vector{Int}
    ti_coefficient_indices::Vector{Int}
    diffusion_state_indices::Vector{Int}
    # Discrete time is not a different filter, only a different discretization:
    # DRIFT, CINT and DIFFUSION are already the one-step quantities, so the
    # matrix exponential, the Lyapunov solve and the intercept solve that turn
    # continuous parameters into per-interval ones all collapse to identities.
    continuous_time::Bool
    # Which manifest variables are Gaussian (0), binary (1) or ordinal (2).
    # Last because the constructor takes it last: a field declared in one order
    # and supplied in another put a Bool into this slot and failed with a
    # convert error pointing at neither.
    #
    # Empty for the Gaussian-only models that were all this engine handled
    # before, so the filter can skip the branch rather than test a vector of
    # zeros on every row.
    manifesttype::Vector{Int}
    # Category count per manifest variable, used only by the ordinal ones: an
    # ordinal variable with `K` categories reads the first `K-1` columns of
    # THRESHOLDS and ignores the rest, so variables with different numbers of
    # categories can share one rectangular matrix.
    ncategories::Vector{Int}
    # Per manifest variable, the row of ASYMPTOTES it reads -- that matrix has
    # one row per binary indicator -- or 0 for the plain two parameter
    # logistic. A property of the model rather than of the values, so it
    # cannot be recovered by looking at them: a guessing parameter of zero is
    # an ordinary place for an optimizer to be.
    asymptoterow::Vector{Int}
    # Censoring limits per manifest variable, read only by the censored ones.
    # Known constants rather than parameters: a scale's floor and ceiling are
    # properties of the instrument, not things the data can inform. Empty means
    # no variable is censored, which is every model that does not ask for it.
    censormin::Vector{Float64}
    censormax::Vector{Float64}
    # Which covariance construction the model asked for, the same code the stan
    # path reads from `standata$choleskymats`: 0 for the unconstrained
    # correlation square root, 2 for covmattransform='z'. Declared last and
    # supplied last, for the reason the `manifesttype` comment above gives.
    covmatcode::Int
    # Which augmented states the population covariance occupies, in the order
    # of the RAWPOPVAR matrix rows: `population_indices[j]` is the state that
    # population row `j` describes. Those indices are not disjoint from the
    # model own latents -- an individually varying T0MEANS gets no carrier
    # state, so its own latent state is what its population row refers to.
    #
    # Empty means the model has no separate population covariance and T0VAR is
    # the whole of the initial covariance, which is every model that does not
    # use intoverpop.
    population_indices::Vector{Int}
    # Where that matrix sits in the flat parameter vector, as first and last
    # positions. A name lookup would not do: `pars` is a ComponentVector whose
    # axis is a type parameter, so `pars.RAWPOPVAR` has to compile for every
    # model, and for a model without the block it cannot. `0:-1` when absent.
    population_range::UnitRange{Int}
    # Which construction builds the population block, which need not be the one
    # the rest of the model's covariances use. A reduced-rank population
    # covariance is a factor -- code 1 -- while DIFFUSION and T0VAR keep the
    # model's own setting, and forcing all four to agree would change three
    # matrices to fix one. Defaults to `covmatcode`, so a model that has not
    # asked for anything different is exactly as it was.
    population_covmatcode::Int
    # The state-unit conversion for each population row, applied where the
    # block is placed. An appended carrier is 1, because its T0MEANS uses the
    # identity transform and the consuming cell does the scaling; an
    # individually varying T0MEANS carries its cell multiplier*meanscale,
    # because the carrier is the model latent itself and its covariance has to
    # be in that latent's units. Deliberately NOT folded into the block's own
    # transform: the same conversion is applied at the same point on the stan
    # path, and two implementations of one rule about this matrix is what the
    # previous arrangement cost.
    population_scale::Vector{Float64}
    # How many leading states are genuine dynamics rather than the static
    # coordinates a random effect augments the state with.
    #
    # The continuous form's local affine offset needs `JAx` inverted, and `JAx`
    # is exactly singular on a static coordinate -- zero row and zero column --
    # so that solve runs over this leading block. Restricting it is exact rather
    # than an approximation: order the states [d, s] and `JAx` is block
    # triangular with zero static rows while the offset is genuinely zero there,
    # so the sub-solve equals the full one, and the static-to-dynamic coupling
    # still arrives through the full matrix exponential.
    #
    # This is Stan's `1:nlatent` (`ctModelWriter.R`:
    # `mdivide_left(JAx[1:nlatent,1:nlatent], ...)`), NOT its `derrind`.
    # `derrind` additionally drops a latent that merely has no diffusion of its
    # own and no coupling to anything that has some -- whose offset is *not*
    # zero -- and using it took that state's CINT out of the likelihood
    # altogether: gradient exactly zero, so the optimizer left the parameter at
    # its starting value and the summary reported that as an estimate.
    #
    # 0 means the whole state vector. Declared last and supplied last, for the
    # reason the `manifesttype` comment above gives.
    affine_dim::Int
    # Whether the latent processes start from the distribution they settle
    # into rather than from T0MEANS and T0VAR; see `_ctsem_stationary!`.
    # Declared last and supplied last, for the reason the `manifesttype`
    # comment above gives.
    stationary::Bool
    # Parallel to `ti_parameter_indices`: the value of a TI effect fixed to one,
    # whose `ti_coefficient_indices` entry is 0 -- no coefficient to estimate,
    # the value times the predictor (`_ctsem_ti_effect`). Declared last and
    # supplied last, for the reason the `manifesttype` comment above gives.
    ti_fixed_values::Vector{Float64}
    # Where a subject's TI predictor values sit in its materialised parameter
    # vector, after the raw parameters (`_materialize_subject_values!`), or
    # `0:-1` when no regular transform reads one. A T0MEANS or T0VAR cell is a
    # regular transform -- a function of that vector alone, evaluated before
    # any row -- so a predictor it reads has to be in the vector; the row
    # groups read theirs from the row context instead (`CTSEMRowContext`).
    # Declared last and supplied last, for the reason the `manifesttype`
    # comment above gives.
    ti_tail::UnitRange{Int}
    # Deviations the filter integrates itself; see `CTSEMHybrid`. Empty unless a
    # Laplace objective filled it. Declared last and supplied last.
    hybrid::CTSEMHybrid
    # Single initial conditions at their stationary values; see
    # `CTSEMStationaryElements`. Empty unless R sent the marks. Declared last
    # and supplied last.
    stationary_elements::CTSEMStationaryElements
    # Parameters with a stated population value; see `CTSEMAnchors`. Empty
    # unless R sent them. Declared last and supplied last.
    anchors::CTSEMAnchors

    # The constructor ensures that the provided vectors are of the correct types and converts them if necessary.
    function EKFParameters(
        mutables,
        transform_indices,
        predict_transforms_indices,
        update_transforms_indices,
        td_transforms_indices,
        regular_transforms,
        predict_transforms,
        update_transforms,
        td_transforms,
        parnumber,
        parameter_axis,
        fixed_indices,
        fixed_values,
        ti_parameter_indices=Int[],
        ti_predictor_indices=Int[],
        ti_coefficient_indices=Int[],
        diffusion_state_indices=Int[],
        continuous_time::Bool=true,
        manifesttype=Int[],
        ncategories=Int[],
        asymptoterow=Int[],
        censormin=Float64[],
        censormax=Float64[],
        covmatcode::Int=0,
        population_indices=Int[],
        population_range::UnitRange{Int}=0:-1,
        population_covmatcode::Union{Nothing,Integer}=nothing,
        population_scale=Float64[],
        affine_dim::Int=0,
        stationary::Bool=false,
        ti_fixed_values=Float64[],
        ti_tail::UnitRange{Int}=0:-1,
    )
        # Every argument from `ti_parameter_indices` on is optional and
        # positional, so inserting a field shifts every caller that passes the
        # ones after it -- silently, because the types are compatible either
        # way: `Vector{Int}([0.0])` succeeds, so a censoring limit of zero
        # becomes an asymptote row. `asymptoterow` was inserted ahead of
        # these two and that is exactly what happened, and the only symptom was
        # generated data clamped at the wrong limit.
        #
        # R builds both limits as n.manifest length vectors or neither, so
        # unequal lengths mean the arguments did not land where the caller
        # meant them to.
        isempty(ti_fixed_values) ||
            length(ti_fixed_values) == length(ti_parameter_indices) ||
            throw(ArgumentError("ti_fixed_values must be empty or one per TI effect"))
        isempty(ti_fixed_values) && any(iszero, ti_coefficient_indices) &&
            throw(ArgumentError("a TI effect with coefficient 0 is fixed and needs ti_fixed_values"))
        length(censormin) == length(censormax) || throw(ArgumentError(
            "censormin and censormax have lengths $(length(censormin)) and " *
            "$(length(censormax)); they are per-manifest vectors and must " *
            "match. A mismatch usually means the optional positional " *
            "arguments are one slot out."))

        parameter_axis = _ctsem_layout(parameter_axis)
        regular_vector = Vector{Any}(collect(regular_transforms))
        regular_groups, regular_group_cells, regular_group_tf =
            _group_regular_transforms(regular_vector, BitVector(mutables))
        predict_transforms_tuple = _ctsem_transform_container(predict_transforms)
        update_transforms_tuple = _ctsem_transform_container(update_transforms)
        td_transforms_tuple = _ctsem_transform_container(td_transforms)
        fixed_values_vector = _concrete_float_vector(fixed_values)

        # Creation of the EKFParameters struct with the appropriate types for each field
        return new{
            typeof(regular_groups),
            typeof(predict_transforms_tuple),
            typeof(update_transforms_tuple),
            typeof(td_transforms_tuple),
            typeof(parameter_axis),
            typeof(fixed_values_vector),
        }(
            BitVector(mutables),
            Vector{Bool}(transform_indices),
            Vector{Bool}(predict_transforms_indices),
            Vector{Bool}(update_transforms_indices),
            Vector{Bool}(td_transforms_indices),
            regular_vector,
            regular_groups,
            regular_group_cells,
            regular_group_tf,
            predict_transforms_tuple,
            update_transforms_tuple,
            td_transforms_tuple,
            Vector{Int}(parnumber),
            parameter_axis,
            Vector{Bool}(fixed_indices),
            fixed_values_vector,
            Vector{Int}(ti_parameter_indices),
            Vector{Int}(ti_predictor_indices),
            Vector{Int}(ti_coefficient_indices),
            Vector{Int}(diffusion_state_indices),
            continuous_time,
            Vector{Int}(manifesttype),
            Vector{Int}(ncategories),
            Vector{Int}(asymptoterow),
            Vector{Float64}(censormin),
            Vector{Float64}(censormax),
            covmatcode,
            Vector{Int}(population_indices),
            population_range,
            population_covmatcode === nothing ? Int(covmatcode) :
                Int(population_covmatcode),
            isempty(population_scale) ?
                ones(Float64, length(population_indices)) :
                Vector{Float64}(population_scale),
            affine_dim,
            stationary,
            isempty(ti_fixed_values) ? zeros(Float64, length(ti_parameter_indices)) :
                Vector{Float64}(ti_fixed_values),
            ti_tail,
            CTSEMHybrid(),
            CTSEMStationaryElements(),
            CTSEMAnchors(),
        )
    end
end

# Compatibility constructor for existing callers without TD transform metadata.
function EKFParameters(mutables, transform_indices, predict_transforms_indices,
    update_transforms_indices, regular_transforms, predict_transforms,
    update_transforms, parnumber, parameter_axis, fixed_indices, fixed_values)
    td_transforms_indices = falses(length(mutables))
    return EKFParameters(mutables, transform_indices, predict_transforms_indices,
        update_transforms_indices, td_transforms_indices, regular_transforms,
        predict_transforms, update_transforms, Function[], parnumber,
        parameter_axis, fixed_indices, fixed_values)
end

function _concrete_float_vector(values)
    isempty(values) && return Float64[]
    T = mapreduce(typeof, promote_type, values)
    T <: AbstractFloat || throw(ArgumentError("fixed_values must contain floating-point values"))
    return Vector{T}(values)
end
