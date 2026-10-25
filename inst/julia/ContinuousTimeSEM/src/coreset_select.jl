# A coreset of units: those whose weighted terms best follow the sum of every
# unit's. One forward selection (`_CoresetPath`) serves the optimiser's stages
# (`_ctsem_coreset_path`), tempered SMC and NUTS on a coreset
# (`_smc_coreset_fit`) and the importance sampler (`_is_coreset_choose`); each
# brings the rows it fits and tests a size on rows the selection did not see.
# What a size is worth is each caller's own, because the four are hurt by
# different errors (review/CORESET-one-rule-2026-10-25.md):
#
# - The optimiser's stages need an accuracy: a stage's step, relative to the
#   full objective's, within a tolerance, so that stages contract on the full
#   maximum and the curvature pairs they share come from one function.
# - The samplers correct a coreset's error exactly -- importance weights, the
#   backward kernel, delayed acceptance -- so an error costs efficiency only,
#   and a size is priced: expected full evaluations per unit of what the
#   sampler is after, against what it costs without a coreset.

"""
    _ctsem_nnls(AtA, Atb)

Lawson and Hanson's non-negative least squares, from the normal equations:
the `x >= 0` minimising `|A x - b|^2`, given `A'A` and `A'b`. Each passive set
is solved by a Cholesky factorization of its block, with a ridge of `1e-12` of
the diagonal's mean so that collinear units do not stop it.
"""
function _ctsem_nnls(AtA::AbstractMatrix, Atb::AbstractVector)
    n = length(Atb)
    x = zeros(n)
    n == 0 && return x
    P = falses(n)
    w = collect(Float64, Atb)
    tol = 1e-10 * max(1.0, maximum(abs, w))
    ridge = 1e-12 * max(sum(AtA[i, i] for i in 1:n) / n, eps())
    solve(idx) = begin
        B = AtA[idx, idx]
        for i in eachindex(idx)
            B[i, i] += ridge
        end
        F = cholesky!(Symmetric(B); check=false)
        issuccess(F) ? F \ Atb[idx] : pinv(AtA[idx, idx]) * Atb[idx]
    end
    for _ in 1:(3n)
        j = 0; best = tol
        for i in 1:n
            (!P[i] && w[i] > best) && (best = w[i]; j = i)
        end
        j == 0 && break
        P[j] = true
        z = zeros(n)
        for _ in 1:(3n)
            idx = findall(P)
            z = zeros(n)
            z[idx] = solve(idx)
            # A loop, not `all(... for i in idx)`: the generator would close
            # over `z`, assigned twice, and box it.
            positive = true
            for i in idx
                z[i] > 0 || (positive = false; break)
            end
            positive && break
            alpha = Inf
            for i in idx
                z[i] <= 0 && (alpha = min(alpha, x[i] / (x[i] - z[i])))
            end
            x .+= alpha .* (z .- x)
            for i in idx
                x[i] <= 1e-12 && (P[i] = false; x[i] = 0.0)
            end
        end
        x = z
        w = Atb .- AtA * x
    end
    return x
end

"""The sizes tried, 1.25 times apart, from one unit to `N`."""
function _ctsem_coreset_sizes(N::Integer)
    sizes = Int[]
    m = 1
    while m <= N
        push!(sizes, m)
        m = max(m + 1, round(Int, 1.25m))
    end
    return sizes
end

"""
    _CoresetPath

One forward-selection path over the columns of `A` (rows by units) for `b`,
held so that it grows a size at a time (`_coreset_grow!`). Each step the unit
whose column, orthogonalised against the units in (a rank-one downdate a
unit), has the largest positive projection on what is left of `b`, squared
over its orthogonalised length, joins, while the units' `cost` stays within
`budget`. What is left is `b`'s projection off the units' span (`refit =
false`), or what non-negative weights on them leave (`refit = true`), which
keeps choosing units that help with a positive weight, at an NNLS solve a
step on a Gram matrix grown a column at a time. The `start` units join first,
in order, whatever they cost. The path ends once nothing positive is left to
add, the next unit would pass the budget, or what is left is within
`resid_tol` of `b`'s length (`done`).
"""
mutable struct _CoresetPath
    A::Matrix{Float64}
    b::Vector{Float64}
    cost::Vector{Float64}
    budget::Float64
    refit::Bool
    start::Vector{Int}
    resid_tol::Float64
    Ap::Matrix{Float64}
    nrm::Vector{Float64}
    floor0::Float64
    r::Vector{Float64}
    nb::Float64
    chosen::Vector{Int}
    spent::Float64
    along::Vector{Float64}
    t::Vector{Float64}
    G::Matrix{Float64}
    g::Vector{Float64}
    w::Vector{Float64}
    taken::Int
    done::Bool
end

function _coreset_path(A::AbstractMatrix{Float64}, b::AbstractVector{Float64},
        cost::AbstractVector{<:Real}, budget::Real; refit::Bool=false,
        start::AbstractVector{<:Integer}=Int[], resid_tol::Real=0.0)
    N = Base.size(A, 2)
    Am = A isa Matrix{Float64} ? A : Matrix{Float64}(A)
    Ap = copy(Am)
    nrm = vec(sum(abs2, Ap; dims=1))
    return _CoresetPath(Am, collect(Float64, b), collect(Float64, cost), Float64(budget),
        refit, collect(Int, start), Float64(resid_tol), Ap, nrm,
        1e-14 * maximum(nrm; init=0.0), collect(Float64, b), sqrt(sum(abs2, b)), Int[],
        0.0, zeros(N), zeros(N), zeros(0, 0), Float64[], Float64[], 0, false)
end

"""Grow path `P` to `target` units, or as far as it goes; whether it got there."""
function _coreset_grow!(P::_CoresetPath, target::Integer)
    N = Base.size(P.A, 2)
    target = min(Int(target), N)
    Ap = P.Ap
    nrm = P.nrm
    chosen = P.chosen
    while length(chosen) < target && !P.done
        if P.taken < length(P.start)
            P.taken += 1
            j = P.start[P.taken]
        else
            if P.resid_tol > 0 && sqrt(sum(abs2, P.r)) <= P.resid_tol * P.nb
                P.done = true
                break
            end
            mul!(P.along, transpose(Ap), P.r)
            j = 0; top = 0.0
            @inbounds for U in 1:N
                a = P.along[U]
                if a > 0 && nrm[U] > P.floor0
                    s = a^2 / nrm[U]
                    s > top && (top = s; j = U)
                end
            end
            if j == 0 || P.spent + P.cost[j] > P.budget
                P.done = true
                break
            end
        end
        push!(chosen, j)
        P.spent += P.cost[j]
        q = Ap[:, j]
        nq = norm(q)
        if nq > 0
            q ./= nq
            mul!(P.t, transpose(Ap), q)
            LinearAlgebra.BLAS.ger!(-1.0, q, P.t, Ap)
            # The lengths by downdating, and afresh every eighth unit, so
            # that rounding cannot leave a spent column looking long.
            if length(chosen) % 8 == 0
                sum!(abs2, reshape(nrm, 1, N), Ap)
            else
                @inbounds for U in 1:N
                    nrm[U] = max(nrm[U] - P.t[U]^2, 0.0)
                end
            end
        end
        for U in chosen
            nrm[U] = 0.0
        end
        if P.refit
            P.G, P.g = _coreset_gram_grow(P.G, P.g, P.A, chosen, P.b)
            P.w = _ctsem_nnls(P.G, P.g)
            P.r .= P.b .- view(P.A, :, chosen) * P.w
        elseif nq > 0
            P.r .-= q .* dot(q, P.r)
        end
    end
    return length(chosen) >= target
end

"""The non-negative least-squares weights of path `P`'s units."""
function _coreset_weights(P::_CoresetPath)
    P.refit && return copy(P.w)
    As = P.A[:, P.chosen]
    return _ctsem_nnls(transpose(As) * As, transpose(As) * P.b)
end

"""
    _coreset_select(A, b, cost, budget, sizes, value; refit=false, start=Int[],
        resid_tol=0.0)

One forward-selection path (`_CoresetPath`) over the columns of `A` for `b`,
scored at each size in `sizes`, ascending: there the weights are
non-negative least squares (`_ctsem_nnls`) and `value(units, weights)`
returns `(score, bound, extra)`: the best score above zero is kept, and the
search stops once `bound`, what any larger coreset could score, is no better,
or the path ends. Returns `(units, weights, extra)` for the best, or
`nothing`.
"""
function _coreset_select(A::AbstractMatrix{Float64}, b::AbstractVector{Float64},
        cost::AbstractVector{<:Real}, budget::Real, sizes, value;
        refit::Bool=false, start::AbstractVector{<:Integer}=Int[], resid_tol::Real=0.0)
    P = _coreset_path(A, b, cost, budget; refit=refit, start=start, resid_tol=resid_tol)
    best = nothing
    bestv = 0.0
    for target in sizes
        reached = _coreset_grow!(P, target)
        isempty(P.chosen) && return nothing
        ws = _coreset_weights(P)
        v, bound, extra = value(P.chosen, ws)
        if v > bestv
            best = (units=copy(P.chosen), weights=ws, extra=extra)
            bestv = v
        end
        # No larger coreset can do better, or there is none.
        (bound <= bestv || !reached) && break
    end
    return best
end

"""
    _coreset_cv(A, b, cost, budget, sizes, splits, judge; refit=false,
        resid_tol=0.0, centre=false)

The out-of-sample test the samplers' coresets are sized by: for each split
(the rows of `A` it holds out), a selection path on the other rows, the paths
grown together through `sizes`. At each size `judge(paths, held)` sees every
split's path (`_CoresetPath`, its units `chosen` and `_coreset_weights`) and
held-out rows, and returns `(score, bound, extra)`, as `_coreset_select`'s
`value` does. With `centre`, each split's rows are centred on its own
training rows' means, for a target whose level is free. Returns `(size,
extra)` at the best score above zero, or `nothing`. The search stops once
`bound` is no better than the best, or any path ends.

The units are chosen out of sample, as their weights are: chosen on every
row and only reweighted per split, a coreset's held-out error read the
selection's own fit to the rows held out.
"""
function _coreset_cv(A::AbstractMatrix{Float64}, b::AbstractVector{Float64},
        cost::AbstractVector{<:Real}, budget::Real, sizes, splits, judge;
        refit::Bool=false, resid_tol::Real=0.0, centre::Bool=false)
    nrows = Base.size(A, 1)
    paths = map(splits) do held
        train = setdiff(1:nrows, held)
        At = A[train, :]
        bt = b[train]
        if centre
            At .-= sum(At; dims=1) ./ length(train)
            bt = bt .- sum(bt) / length(bt)
        end
        _coreset_path(At, bt, cost, budget; refit=refit, resid_tol=resid_tol)
    end
    best = nothing
    bestv = 0.0
    for target in sizes
        reached = true
        for P in paths
            reached &= _coreset_grow!(P, target)
        end
        any(P -> isempty(P.chosen), paths) && break
        v, bound, extra = judge(paths, splits)
        if v > bestv
            best = (size=maximum(length(P.chosen) for P in paths), extra=extra)
            bestv = v
        end
        (bound <= bestv || !reached) && break
    end
    return best
end

"""Rows `1:n` split into `k` folds at random (`rng`), a vector of each fold's
rows."""
function _coreset_folds(n::Integer, k::Integer, rng::AbstractRNG)
    perm = randperm(rng, Int(n))
    edges = round.(Int, range(0, n; length=k + 1))
    return [perm[(edges[f] + 1):edges[f + 1]] for f in 1:k]
end

"""`G` and `g`, the Gram matrix of the units in and their columns' products with
`b`, grown by the newest of `chosen`."""
function _coreset_gram_grow(G::Matrix{Float64}, g::Vector{Float64},
        A::AbstractMatrix{Float64}, chosen::Vector{Int}, b::AbstractVector{Float64})
    k = length(chosen)
    Gn = Matrix{Float64}(undef, k, k)
    @inbounds for c in 1:(k - 1), i in 1:(k - 1)
        Gn[i, c] = G[i, c]
    end
    aj = view(A, :, chosen[k])
    for i in 1:(k - 1)
        v = dot(view(A, :, chosen[i]), aj)
        Gn[i, k] = v
        Gn[k, i] = v
    end
    Gn[k, k] = dot(aj, aj)
    gn = Vector{Float64}(undef, k)
    copyto!(gn, 1, g, 1, k - 1)
    gn[k] = dot(aj, b)
    return Gn, gn
end
