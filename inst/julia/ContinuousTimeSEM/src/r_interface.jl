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

ctsem writes these out longhand (`1/(1+exp(-x))`) in the transforms it
generates, so nothing currently needs them -- but the C++ engine's expression
parser accepts them by name, and the two engines should accept the same
vocabulary rather than differing over which spellings of the same function work.
"""
inv_logit(x) = 1 / (1 + exp(-x))
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
    retrieve_axes(matrix, row, col)

Build a `ComponentArrays.Axis` from the matrix layout encoded in the parameter
table's coordinate columns.

Each distinct matrix name becomes a shaped component whose row and column
dimensions are inferred from the `row` and `col` columns.
"""
function retrieve_axes(matrix, row, col)
    names, nrows, ncols = retrieve_names_and_dims(matrix, row, col)
    endings = cumsum([a * b for (a, b) in zip(nrows, ncols)])
    startings = vcat(1, endings[1:end-1] .+ 1)
    ranges = [starts:ends for (starts, ends) in zip(startings, endings)]
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
  * `continuous_time` — `false` for a discrete-time model, whose DRIFT, CINT and
    DIFFUSION are already the one-step quantities.

Expression strings are `Meta.parse`d and `eval`ed into closures here, once per
model. That is what makes this engine model-agnostic without a compile step.
"""
function ekf_from_columns(matrix, row, col, parnumber, value, transform,
    predicttransform, updatetransform, tdtransform;
    ti_parameter=Int[], ti_predictor=Int[], ti_coefficient=Int[],
    diffusion_state_indices=Int[], continuous_time::Bool=true)

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

    # Which flattened position each free parameter occupies.
    map_from = ComponentVector(zeros(Int, n), axis)
    @inbounds for i in 1:n
        pn = Int(parnumber[i])
        pn == 0 && continue
        @view(map_from[Symbol(matrix[i])])[Int(row[i]), Int(col[i])] = pn
    end
    par_pos = getdata(map_from) .|> !iszero
    map_from = map_from[par_pos]

    # Regular (population) transforms.
    reg_tfs = ComponentVector(Vector{Function}(undef, n), axis)
    @inbounds for i in 1:n
        pn = Int(parnumber[i])
        pn == 0 && continue
        expression = isempty(transform[i]) ? "param[$(pn)]" : String(transform[i])
        @view(reg_tfs[Symbol(matrix[i])])[Int(row[i]), Int(col[i])] =
            expression |> generate_transform_string |> Meta.parse |> eval
    end
    tfs_pos = par_pos
    reg_tfs = getdata(reg_tfs)[tfs_pos]

    # State-dependent transforms, split by the point in the row at which the
    # forward filter needs them to be current.
    predict_tfs = ComponentVector(Array{Function}(undef, n), axis)
    update_tfs = ComponentVector(Array{Function}(undef, n), axis)
    td_tfs = ComponentVector(Array{Function}(undef, n), axis)
    @inbounds for i in 1:n
        target_row = Int(row[i])
        target_col = Int(col[i])
        name = Symbol(matrix[i])
        if !isempty(predicttransform[i])
            @view(predict_tfs[name])[target_row, target_col] =
                String(predicttransform[i]) |> generate_complex_transform_string |> Meta.parse |> eval
        end
        if !isempty(updatetransform[i])
            @view(update_tfs[name])[target_row, target_col] =
                String(updatetransform[i]) |> generate_complex_transform_string |> Meta.parse |> eval
        end
        if !isempty(tdtransform[i])
            @view(td_tfs[name])[target_row, target_col] =
                String(tdtransform[i]) |> generate_complex_transform_string |> Meta.parse |> eval
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

    return EKFParameters(par_pos, tfs_pos, ptf_pos, utf_pos, ttf_pos,
        reg_tfs, predict_tfs, update_tfs, td_tfs, map_from, axis,
        fixed_positions, fixed_values, Int.(ti_parameter),
        Int.(ti_predictor), Int.(ti_coefficient), Int.(diffusion_state_indices),
        continuous_time)
end
