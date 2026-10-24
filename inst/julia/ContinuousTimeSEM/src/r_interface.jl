using ComponentArrays

# The R-facing constructor for `EKFParameters`.
#
# This takes the parameter table as **plain column vectors**, not as a
# `DataFrame`. That is a deployment decision rather than a stylistic one:
# DataFrames was this package's single most expensive dependency, and it was
# used only here, only as a row container. Dropping it takes the dependency
# closure from 111 packages to 56, a fresh install from 268 MB to 124 MB, and
# `using ContinuousTimeSEM` from ~8.7 s to ~3.8 s -- a cost paid in every R
# session that touches the backend. It remains a test-only dependency, because
# a `DataFrame` literal is still the most readable way to write a table in a
# test; see `test/table_helpers.jl`.
#
# Absent entries are encoded as sentinels rather than `missing`, so every column
# is concretely typed: `parnumber == 0` means "not a free parameter",
# `isnan(value)` means "no fixed value", and `""` means "no transform". The R
# side already has to normalise NA anyway, so this moves no work, and it keeps
# `Union{Missing,T}` out of the hot constructor.

"""
    log1p_exp(x)

Return `log(1 + exp(x))`.

This helper is intended for transform expressions passed through the R
interface.
"""
log1p_exp(x) = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))

"""
    inv_logit(x)
    logit(x)

The logistic function and its inverse.

Symbolic Jacobian cells use `inv_logit` for every logistic they differentiate
(`ctJacobian`, R), so it is written not to overflow: from the side where `exp`
is at most 1. `1 / (1 + exp(-x))` was right in value everywhere, but below
x = -709 `exp(-x)` is Inf, and a dual number's partial through it is Inf / Inf.
"""
inv_logit(x) = x >= 0 ? inv(1 + exp(-x)) : (e = exp(x); e / (1 + e))
logit(x) = log(x / (1 - x))

"""
    generate_transform_string(s)

Wrap a scalar parameter transform expression string as a Julia lambda string.

The returned string has the form `"param -> ..."` and is parsed later by
`ekf_from_columns`.
"""
function generate_transform_string(s)
    return "param -> " * s
end

"""
    generate_complex_transform_string(s)

Wrap a state-dependent transform expression string as a Julia lambda string.

The returned function accepts a `CTSEMRowContext` and is parsed later by
`ekf_from_columns`.
"""
function generate_complex_transform_string(s)
    expression = replace(s, r"\bstate\b" => "ctx.state", r"\bPARS\b" => "ctx.pars.PARS")
    for matrix_name in ("T0MEANS", "LAMBDA", "DRIFT", "DIFFUSION", "MANIFESTVAR",
        "MANIFESTMEANS", "CINT", "T0VAR", "TDPREDEFFECT", "JAx", "Jtd", "Jy")
        expression = replace(expression, Regex("\\b" * matrix_name * "\\s*\\[") =>
            "ctx.pars." * matrix_name * "[")
    end
    return "ctx -> " * expression
end

"""
    _transform_closure(source)

The compiled closure for one transform expression, reused across models.

`Meta.parse |> eval` mints a *new anonymous function type* every time it runs,
and that type is a type parameter of `EKFParameters`. So two models with
identical transform strings are, to the compiler, different types -- and the
whole EKF pipeline, which is generic in that type, specialises again for each
of them. Measured on a one-latent model: the first `ctsem_laplace_evaluate`
took 34.5 s where the second took 0.02 s, and a second model built from the
same strings but a fresh `ekf_from_columns` call paid 28.6 s more. A
replication study fitting the same model to a hundred datasets paid that a
hundred times.

Caching by the *source string* removes it. The same string yields the same
closure object, so `typeof` agrees, so `EKFParameters` agrees, so the compiled
code is reused: the second model of a shape is free.

The lock is not decoration. `ekf_from_columns` is ordinarily called from one
thread, but nothing in its contract says so, and a `Dict` torn by two
concurrent inserts fails in ways that look like a model error.
"""
function _transform_closure(source::AbstractString)
    key = String(source)
    lock(_TRANSFORM_CACHE_LOCK) do
        cached = get(_TRANSFORM_CACHE, key, nothing)
        cached === nothing || return cached
        built = eval(Meta.parse(key))
        _TRANSFORM_CACHE[key] = built
        return built
    end
end

"""
    _regular_transform_closure(source)

The compiled closure for one *regular* transform, with its parameter index
carried as a value rather than baked into the code.

`ctModelWriter` writes the index into the expression -- `param[3]`, `param[7]` --
so every cell of every model compiles to its own anonymous function *type*, and
`EKFParameters` is parameterised by a tuple of them. Two models are then
different types however similar they are, and the whole EKF pipeline specialises
again for each: measured at 20-70 s on the first evaluation of a model shape,
paid again for the next model.

Rewriting `param[3]` to `param[i]` under a closure over `i` collapses that. All
cells sharing a *template* now share one type and differ only in a field, so a
model is typed by which templates it uses, and -- since `EKFParameters` groups
its regular transforms by type rather than holding one tuple entry per cell --
not by how many cells use each or which index they read. The engine's own
transforms come from a handful of templates, so a second model reuses the
first one's compiled code even when its size and layout differ.

Almost every regular transform reads exactly one free parameter, so there is
exactly one index to lift out. The exception is a T0MEANS or T0VAR cell
composed from several PARS parameters, which reads one index per parameter; a
template cannot represent that with a single `i`, so `_regular_transform_index`
returns `nothing` for it and it falls through to `_transform_closure` with its
literal indices, minting its own specialisation. That is the right trade for a
cell that appears once or twice in a model, and it is why the fall-through is a
deliberate branch rather than a guess.
"""
function _regular_transform_closure(expression::AbstractString; compile::Bool=false)
    text = String(expression)
    index = _regular_transform_index(text)
    # Outside the default templates: an expression, so the model's type does
    # not depend on it (transform_expr.jl) -- unless asked to compile.
    if !compile && (index === nothing || !_ctsem_default_template(text))
        e = _ctsem_regular_expr(text)
        e === nothing || return e
    end
    index === nothing && return _transform_closure(generate_transform_string(text))
    # `invokelatest` because the maker may have been `eval`ed into this module a
    # moment ago, and a method compiled after the current frame started is not
    # visible to it. This runs once per model cell, so the dynamic call costs
    # nothing worth measuring -- and the closure it returns is called normally.
    return Base.invokelatest(_transform_closure(_regular_transform_key(text)), index)
end

"""
The single parameter index a regular transform reads, or `nothing` when it
reads none or several -- see `_regular_transform_closure` for what happens to
the several case.
"""
function _regular_transform_index(text::AbstractString)
    indices = Set{Int}()
    for m in eachmatch(r"param\[\s*(\d+)\s*\]", text)
        push!(indices, parse(Int, m.captures[1]))
    end
    return length(indices) == 1 ? first(indices) : nothing
end

"""The cache key for a regular transform: its source with the index lifted out."""
function _regular_transform_key(expression::AbstractString)
    text = String(expression)
    _regular_transform_index(text) === nothing &&
        return generate_transform_string(text)
    template = replace(text, r"param\[\s*\d+\s*\]" => "param[ctsem_index]")
    return "ctsem_index -> (param -> " * template * ")"
end

const _TRANSFORM_CACHE = Dict{String,Function}()
const _TRANSFORM_CACHE_LOCK = ReentrantLock()

export ctsem_transforms_cached
"""
    ctsem_transforms_cached(sources)

Whether every one of these transform expressions already has a closure.

This is how the R side knows, *before* it builds anything, whether a model is a
shape this session has compiled for. A model whose transforms are all cached
reuses the compiled filter and evaluates immediately; one with a new expression
mints new closure types and specialises the whole pipeline again, which takes
tens of seconds and is worth saying out loud rather than looking like a hang.
"""
function ctsem_transforms_cached(regular, complex=String[]; compile::Bool=false)
    lock(_TRANSFORM_CACHE_LOCK) do
        # The cache is keyed by the *lambda* the builder compiles, not by the
        # expression the caller wrote, so the same wrapping has to be applied
        # here -- and it differs between the two kinds. Asking with the raw
        # expression finds nothing, always, which reads as "never compiled".
        for source in regular
            expression = String(source)
            isempty(expression) && continue
            # Regular transforms are cached by *template*, with the parameter
            # index lifted out -- see `_regular_transform_closure` -- so the
            # question is whether this expression's template is known, not
            # whether this exact expression has been seen.
            compile || _ctsem_default_template(expression) ||
                _ctsem_regular_expr(expression) === nothing || continue
            haskey(_TRANSFORM_CACHE, _regular_transform_key(expression)) ||
                return false
        end
        for source in complex
            expression = String(source)
            isempty(expression) && continue
            compile || _ctsem_expr(expression, nothing) === nothing || continue
            haskey(_TRANSFORM_CACHE,
                generate_complex_transform_string(expression)) || return false
        end
        return true
    end
end

"""
    map_fixed_values!(target, val_pos, values)

Write `values` into `target[val_pos]` in place and return `nothing`.
"""
function map_fixed_values!(target, val_pos, values)
    target[val_pos] .= values
    return nothing
end

"""
    map_fixed_values!(target, val_pos::AbstractVector{Bool}, values)

Mask-indexed specialization that writes in place without allocating.

`target[mask] .= values` materializes a temporary for the masked selection on
every call. This runs once per subject per objective evaluation, which made it
the second largest allocation source in the primal filter (6.3 MB per 20
evaluations of a 20-latent model) despite doing almost no arithmetic.
"""
function map_fixed_values!(target, val_pos::AbstractVector{Bool}, values)
    k = 0
    @inbounds for i in eachindex(val_pos)
        if val_pos[i]
            k += 1
            target[i] = values[k]
        end
    end
    return nothing
end

"""
    map_fixed_values(target, val_pos, values)

Write `values` into `target[val_pos]` in place and return `target`.
"""
function map_fixed_values(target, val_pos, values)
    target[val_pos] .= values
    return target
end


export params_to_all_params
"""
    params_to_all_params(params, mutables)

Expand the vector of free parameter values into a full parameter vector.

Positions marked `true` in `mutables` receive entries from `params`; other
positions are left uninitialized for later filling by fixed values.
"""
function params_to_all_params(params, mutables)
    all_params = similar(params, length(mutables))
    all_params[mutables] .= params
    return all_params
end

"""
    retrieve_names_and_dims(matrix, row, col)

Return matrix names and dimensions described by the parameter table's
coordinate columns.

The returned tuple is `(names, nrows, ncols)`, where names are `Symbol`s in
order of first appearance. Dimensions come from the largest row/column index
seen for each matrix, so a table listing cells in any order yields the same
layout.
"""
function retrieve_names_and_dims(matrix, row, col)
    names = Symbol[]
    nrows = Int[]
    ncols = Int[]
    positions = Dict{Symbol,Int}()
    @inbounds for i in eachindex(matrix)
        name = Symbol(matrix[i])
        slot = get(positions, name, 0)
        if slot == 0
            push!(names, name)
            push!(nrows, Int(row[i]))
            push!(ncols, Int(col[i]))
            positions[name] = length(names)
        else
            nrows[slot] = max(nrows[slot], Int(row[i]))
            ncols[slot] = max(ncols[slot], Int(col[i]))
        end
    end
    return (names, nrows, ncols)
end

"""
    retrieve_layout(matrix, row, col)

Where each matrix sits in the flat parameter vector: `(names, nrows, ncols,
ranges)`, with `ranges[j]` the column-major slice belonging to `names[j]`.

One function because two callers need the same answer and one of them cannot
notice being wrong. `retrieve_axes` uses the ranges to build the
ComponentArrays axis; `ekf_from_columns` uses them to find the population
block. That block's range indexes a `view` into the parameter vector, so a
range computed independently and disagreeing by one would not error -- it would
read the wrong slice, build a covariance out of it, and return a number.
"""
function retrieve_layout(matrix, row, col)
    names, nrows, ncols = retrieve_names_and_dims(matrix, row, col)
    endings = cumsum([a * b for (a, b) in zip(nrows, ncols)])
    startings = vcat(1, endings[1:end-1] .+ 1)
    ranges = [starts:ends for (starts, ends) in zip(startings, endings)]
    return (names, nrows, ncols, ranges)
end

"""
    retrieve_axes(matrix, row, col)

Build a `ComponentArrays.Axis` from the matrix layout encoded in the parameter
table's coordinate columns.

Each distinct matrix name becomes a shaped component whose row and column
dimensions are inferred from the `row` and `col` columns.
"""
function retrieve_axes(matrix, row, col)
    names, nrows, ncols, ranges = retrieve_layout(matrix, row, col)
    symbols = tuple(names...)
    axes = ViewAxis.(ranges, ShapedAxis.(zip(nrows, ncols)))
    named_axes = NamedTuple{symbols}(axes)
    return Axis(named_axes)
end

export ekf_from_columns
"""
    ekf_from_columns(matrix, row, col, parnumber, value, transform,
                     predicttransform, updatetransform, tdtransform;
                     ti_parameter, ti_predictor, ti_coefficient,
                     diffusion_state_indices)

Construct an `EKFParameters` object from a parameter specification given as
column vectors, one entry per model-matrix cell.

  * `matrix`, `row`, `col` — the cell's coordinates. Matrix names may be
    `String`s or `Symbol`s.
  * `parnumber` — index of the free parameter occupying the cell, or `0`.
  * `value` — the cell's fixed value, or `NaN` if it has none.
  * `transform` — the free parameter's transform as Julia source, or `""`.
  * `predicttransform`, `updatetransform`, `tdtransform` — state-dependent
    expressions evaluated before the prediction, measurement and TD-impulse
    steps respectively, or `""`.
  * `ti_parameter`, `ti_predictor`, `ti_coefficient` — one entry per
    TI-predictor effect. Keyword arguments, so a model without TI predictors
    simply omits them; that also keeps a caller from having to marshal an empty
    vector across a language boundary, which is not something every bridge
    does well.
  * `diffusion_state_indices` — the states carrying their own diffusion,
    defaulting to all of them.
  * `affine_dim` — how many leading states are genuine dynamics rather than
    the static coordinates a random effect augments the state with. The
    continuous form's intercept solve runs over that leading block, because
    `JAx` is singular on a static coordinate. `0` (the default) means the
    whole state vector. Distinct from `diffusion_state_indices`: a latent
    with no diffusion of its own is still genuine dynamics and still has an
    intercept.
  * `manifesttype` — `0` for a Gaussian manifest variable, `1` for binary,
    `2` for ordinal. Empty (the default) means every variable is Gaussian,
    which lets the filter skip the branch entirely rather than test a vector of
    zeros.
  * `ncategories` — number of categories per manifest variable. Read only for
    the ordinal ones, which take their thresholds from the first
    `ncategories - 1` columns of `THRESHOLDS`.
  * `asymptoterow` — per manifest variable, the row of `ASYMPTOTES` (lower,
    upper; one row per binary indicator) whose asymptotes its response
    probability runs between, or `0` for the plain two parameter logistic.
  * `continuous_time` — `false` for a discrete-time model, whose DRIFT, CINT and
    DIFFUSION are already the one-step quantities.
  * `stationary` — `true` to start the latent processes from the distribution
    DRIFT, CINT and DIFFUSION settle into, in place of T0MEANS and T0VAR. See
    `_ctsem_stationary!` for what it requires of the model.
  * `stationary_mean`, `stationary_block`, `stationary_closure` — single
    initial conditions at their stationary values instead
    (`CTSEMStationaryElements`), each one entry per latent: whether its
    T0MEANS cell is marked; the block its marked T0VAR cells belong to, or 0;
    and whether its moments are needed, the marked latents and every latent
    that drives them.
  * `compile_transforms` — `true` to compile every state-dependent cell and
    every parameter transform outside the default templates into its own
    closure, as the engine did before transforms were data. The model then
    has its own type and compiles the filter for itself on first use; its
    transforms then cost what compiled code costs, which for a model of many
    large state-dependent cells can be well worth the compile.

State-dependent cells, and parameter transforms outside the default
templates, are held as `CTSEMExpr` data (transform_expr.jl), so a model's type
does not depend on them and a model of a new shape reuses the code every other
model compiled. An expression outside that vocabulary, and every one under
`compile_transforms`, is `Meta.parse`d and `eval`ed into a closure, once per
distinct expression (`_transform_closure`).
"""
function ekf_from_columns(matrix, row, col, parnumber, value, transform,
    predicttransform, updatetransform, tdtransform;
    ti_parameter=Int[], ti_predictor=Int[], ti_coefficient=Int[],
    diffusion_state_indices=Int[], continuous_time::Bool=true,
    manifesttype=Int[], ncategories=Int[], asymptoterow=Int[],
    censormin=Float64[],
    censormax=Float64[], covmatcode::Int=0, population_indices=Int[],
    population_covmatcode::Union{Nothing,Integer}=nothing,
    population_scale=Float64[], affine_dim::Int=0, stationary::Bool=false,
    ti_fixed_value=Float64[], ti_tail_offset::Integer=0, n_tipreds::Integer=0,
    compile_transforms::Bool=false, stationary_mean=Int[], stationary_block=Int[],
    stationary_closure=Int[], nfree::Integer=0, anchor_values=Float64[])

    n = length(matrix)
    length(row) == n && length(col) == n ||
        throw(DimensionMismatch("matrix, row and col must have one entry per cell"))
    for (name, column) in (("parnumber", parnumber), ("value", value),
        ("transform", transform), ("predicttransform", predicttransform),
        ("updatetransform", updatetransform), ("tdtransform", tdtransform))
        length(column) == n ||
            throw(DimensionMismatch("column `$name` has $(length(column)) entries, expected $n"))
    end

    axis = retrieve_axes(matrix, row, col)
    values = ComponentVector(fill(NaN, n), axis)

    # A regular transform that reads a TI predictor -- a T0MEANS or T0VAR cell
    # (`.ctJuliaResolveStaticRefs`) -- reads it from the tail of the subject's
    # materialised vector, `ti_tail`, so `tipreds[j]` compiles to
    # `param[ti_tail_offset + j]`. A cell that reads predictors and no
    # parameter arrives with parnumber 0 and takes the first predictor it
    # reads as its representative position, which the gradient's support
    # discovery checks the transform reads (`_ctsem_regular_transform_supports`).
    ti_tail, regular_text, regular_parnumber = _ctsem_ti_tail_rendering(
        transform, parnumber, ti_tail_offset, n_tipreds)

    # Which flattened position each free parameter occupies.
    map_from = ComponentVector(zeros(Int, n), axis)
    @inbounds for i in 1:n
        pn = regular_parnumber[i]
        pn == 0 && continue
        @view(map_from[Symbol(matrix[i])])[Int(row[i]), Int(col[i])] = pn
    end
    par_pos = getdata(map_from) .|> !iszero
    map_from = map_from[par_pos]

    # Regular (population) transforms.
    reg_tfs = ComponentVector(Vector{Any}(undef, n), axis)
    @inbounds for i in 1:n
        pn = regular_parnumber[i]
        pn == 0 && continue
        expression = isempty(regular_text[i]) ? "param[$(pn)]" : regular_text[i]
        @view(reg_tfs[Symbol(matrix[i])])[Int(row[i]), Int(col[i])] =
            _regular_transform_closure(expression; compile=compile_transforms)
    end
    tfs_pos = par_pos
    reg_tfs = getdata(reg_tfs)[tfs_pos]

    # State-dependent transforms, split by the point in the row at which the
    # forward filter needs them to be current.
    # Each held as a `CTSEMExpr` when it is in that vocabulary, so the model's
    # type does not depend on its expressions (transform_expr.jl); otherwise,
    # or when asked, compiled into a closure.
    layout = CTSEMLayout(retrieve_layout(matrix, row, col)...)
    complex_transform(source) = begin
        e = compile_transforms ? nothing : _ctsem_expr(source, layout)
        e === nothing || return e
        _transform_closure(_ctsem_flat_cells(generate_complex_transform_string(source), layout))
    end
    predict_tfs = ComponentVector(Array{Any}(undef, n), axis)
    update_tfs = ComponentVector(Array{Any}(undef, n), axis)
    td_tfs = ComponentVector(Array{Any}(undef, n), axis)
    @inbounds for i in 1:n
        target_row = Int(row[i])
        target_col = Int(col[i])
        name = Symbol(matrix[i])
        if !isempty(predicttransform[i])
            @view(predict_tfs[name])[target_row, target_col] =
                complex_transform(String(predicttransform[i]))
        end
        if !isempty(updatetransform[i])
            @view(update_tfs[name])[target_row, target_col] =
                complex_transform(String(updatetransform[i]))
        end
        if !isempty(tdtransform[i])
            @view(td_tfs[name])[target_row, target_col] =
                complex_transform(String(tdtransform[i]))
        end
    end

    ptfs_raw = getdata(predict_tfs)
    utfs_raw = getdata(update_tfs)
    ttfs_raw = getdata(td_tfs)
    ptf_pos = [isassigned(ptfs_raw, idx) for idx in eachindex(ptfs_raw)]
    utf_pos = [isassigned(utfs_raw, idx) for idx in eachindex(utfs_raw)]
    ttf_pos = [isassigned(ttfs_raw, idx) for idx in eachindex(ttfs_raw)]
    predict_tfs = predict_tfs[ptf_pos]
    update_tfs = update_tfs[utf_pos]
    td_tfs = td_tfs[ttf_pos]

    # Fixed values.
    @inbounds for i in 1:n
        v = value[i]
        isnan(v) && continue
        @view(values[Symbol(matrix[i])])[Int(row[i]), Int(col[i])] = v
    end
    fixed_positions = getdata(values) .|> !isnan
    fixed_values = values[fixed_positions]

    # The population block's slice of the flat vector, read from the same
    # layout the axis is built from rather than derived again. `RAWPOPVAR` is
    # present only when the R side emits it, which it does for a model whose
    # population covariance is its own matrix rather than part of T0VAR.
    population_range = let (bnames, _, _, branges) =
            retrieve_layout(matrix, row, col)
        slot = findfirst(isequal(:RAWPOPVAR), bnames)
        # Parenthesised: `? 0:-1 :` puts the range operator and the
        # conditional's own colon in the same expression, which parses as
        # something else entirely.
        slot === nothing ? (0:-1) : branges[slot]
    end

    sp = EKFParameters(par_pos, tfs_pos, ptf_pos, utf_pos, ttf_pos,
        reg_tfs, predict_tfs, update_tfs, td_tfs, map_from, layout,
        fixed_positions, fixed_values, Int.(ti_parameter),
        Int.(ti_predictor), Int.(ti_coefficient), Int.(diffusion_state_indices),
        continuous_time, Int.(manifesttype), Int.(ncategories),
        Int.(asymptoterow),
        Float64.(censormin), Float64.(censormax), covmatcode,
        Int.(population_indices), population_range,
        population_covmatcode === nothing ? covmatcode :
            Int(population_covmatcode),
        Float64.(population_scale), affine_dim, stationary,
        Float64.(ti_fixed_value), ti_tail)
    _ctsem_set_stationary_elements!(sp.stationary_elements, stationary_mean,
        stationary_block, stationary_closure)
    # Parameters with a stated population value (`CTSEMAnchors`): their cells
    # read `param[nfree + j]`, and the TI predictor tail follows them, which
    # `ti_tail_offset` already says.
    if !isempty(anchor_values)
        nfree > 0 || throw(ArgumentError("anchor_values needs nfree"))
        sp.anchors.nfree = Int(nfree)
        sp.anchors.values = Float64.(anchor_values)
    end
    return sp
end

"""
    _ctsem_set_stationary_elements!(se, mean, block, closure)

Fill `CTSEMStationaryElements` from R's per-latent vectors: nonzero `mean`
entries mark T0MEANS cells, equal nonzero `block` entries put latents in one
T0VAR block, nonzero `closure` entries are the latents whose moments are
needed. Nothing when `closure` is empty.
"""
function _ctsem_set_stationary_elements!(se::CTSEMStationaryElements, mean, block,
        closure)
    isempty(closure) && return se
    se.inclosure = Vector{Bool}(closure .!= 0)
    se.means = findall(!=(0), mean)
    ids = sort(unique(filter(!=(0), Int.(block))))
    se.blocks = [findall(==(id), Int.(block)) for id in ids]
    all(i -> se.inclosure[i], se.means) && all(B -> all(i -> se.inclosure[i], B),
        se.blocks) || throw(ArgumentError(
        "every latent with a stationary element must be in the closure"))
    return se
end

const _CTSEM_TIPRED_REF = r"\btipreds\s*\[\s*(\d+)\s*\]"

const _CTSEM_CELL_REF = r"ctx\.pars\.([A-Za-z_][A-Za-z0-9_]*)\[\s*(\d+)\s*(?:,\s*(\d+)\s*)?\]"

"""
    _ctsem_flat_cells(text, layout)

A compiled transform's source with each `ctx.pars.M[i, j]` read as
`ctx.pars[k]`, `k` the cell's place in the flat parameter vector. The named
read looks the layout up and builds a view every time it runs, where the
ComponentArrays axis this replaced made it a constant offset: a transform
compiled from the named form ran slower than the same cell as a `CTSEMExpr`,
which resolves its cells at build time. A cell the layout does not have is
left as written, so it fails where it always did.
"""
function _ctsem_flat_cells(text::AbstractString, layout::CTSEMLayout)
    return replace(String(text), _CTSEM_CELL_REF => function (ref)
        m = match(_CTSEM_CELL_REF, ref)
        i = _layout_index(layout, Symbol(m.captures[1]))
        i == 0 && return ref
        nr = layout.nrows[i]
        nc = layout.ncols[i]
        r = parse(Int, m.captures[2])
        if m.captures[3] === nothing
            1 <= r <= nr * nc || return ref
            return string("ctx.pars[", layout.first[i] + r - 1, "]")
        end
        c = parse(Int, m.captures[3])
        (1 <= r <= nr && 1 <= c <= nc) || return ref
        return string("ctx.pars[", layout.first[i] + (c - 1) * nr + (r - 1), "]")
    end)
end

"""
    _ctsem_ti_tail_rendering(transform, parnumber, offset, ntipreds)

The regular transforms' texts with each `tipreds[j]` read rewritten to
`param[offset + j]`, the parameter number of each cell (a predictor-only cell's
first read for its 0), and the tail `offset+1:offset+ntipreds` -- `0:-1` when
no transform reads a predictor. See `EKFParameters.ti_tail`.
"""
function _ctsem_ti_tail_rendering(transform, parnumber, offset::Integer,
    ntipreds::Integer)
    texts = String[String(t) for t in transform]
    numbers = Int[Int(p) for p in parnumber]
    reads = [occursin(_CTSEM_TIPRED_REF, t) for t in texts]
    any(reads) || return (0:-1), texts, numbers
    (offset > 0 && ntipreds > 0) || throw(ArgumentError(string(
        "a T0MEANS or T0VAR cell reads a TI predictor, but the parameter ",
        "vector's length (ti_tail_offset = ", offset, ") or the number of TI ",
        "predictors (n_tipreds = ", ntipreds, ") was not given")))
    for i in findall(reads)
        js = [parse(Int, m.captures[1]) for m in eachmatch(_CTSEM_TIPRED_REF, texts[i])]
        all(j -> 1 <= j <= ntipreds, js) || throw(ArgumentError(string(
            "transform `", texts[i], "` reads a TI predictor outside 1:", ntipreds)))
        texts[i] = replace(texts[i], _CTSEM_TIPRED_REF =>
            s -> string("param[", offset + parse(Int, match(r"\d+", s).match), "]"))
        numbers[i] == 0 && (numbers[i] = offset + first(js))
    end
    return (offset + 1):(offset + ntipreds), texts, numbers
end
