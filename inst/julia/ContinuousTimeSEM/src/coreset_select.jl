# A coreset of units: those whose weighted terms best follow the sum of every
# unit's. One forward selection serves the optimiser's stages
# (`_ctsem_coreset_path`), tempered SMC (`_smc_coreset_fit`) and the
# importance sampler; each brings the rows it fits and the rule that judges a
# size (review/CORESET-one-selector-2026-10-24.md).

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
            all(z[i] > 0 for i in idx) && break
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
    _coreset_select(A, b, cost, budget, sizes, value; refit=false, start=Int[],
        resid_tol=0.0)

Forward selection over the columns of `A` (rows by units) for `b`. Each step
the unit whose column, orthogonalised against the units in (a rank-one
downdate a unit), has the largest positive projection on what is left of `b`,
squared over its orthogonalised length, joins, while the units' `cost` stays
within `budget`. What is left is `b`'s projection off the units' span
(`refit = false`), or what non-negative weights on them leave (`refit =
true`), which keeps choosing units that help with a positive weight, at an
NNLS solve a step on a Gram matrix grown a column at a time. The `start` units
join first, in order, whatever they cost. The path stops early once nothing
positive is left to add, or what is left is within `resid_tol` of `b`'s length.

At each size in `sizes`, ascending, the weights are non-negative least
squares (`_ctsem_nnls`) and `value(units, weights)` returns `(score, bound,
extra)`: the best score above zero is kept, and the search stops once `bound`,
what any larger coreset could score, is no better. Returns `(units, weights,
extra)` for the best, or `nothing`.
"""
function _coreset_select(A::AbstractMatrix{Float64}, b::AbstractVector{Float64},
        cost::AbstractVector{<:Real}, budget::Real, sizes, value;
        refit::Bool=false, start::AbstractVector{<:Integer}=Int[], resid_tol::Real=0.0)
    N = Base.size(A, 2)
    Ap = Matrix{Float64}(A)
    nrm = vec(sum(abs2, Ap; dims=1))
    floor0 = 1e-14 * maximum(nrm; init=0.0)
    r = collect(Float64, b)
    nb = sqrt(sum(abs2, b))
    chosen = Int[]
    spent = 0.0
    along = zeros(N)
    t = zeros(N)
    G = zeros(0, 0)
    g = Float64[]
    w = Float64[]
    taken = 0
    best = nothing
    bestv = 0.0
    for target in sizes
        target = min(Int(target), N)
        while length(chosen) < target
            if taken < length(start)
                taken += 1
                j = Int(start[taken])
            else
                resid_tol > 0 && sqrt(sum(abs2, r)) <= resid_tol * nb && break
                mul!(along, transpose(Ap), r)
                j = 0; top = 0.0
                @inbounds for U in 1:N
                    a = along[U]
                    if a > 0 && nrm[U] > floor0
                        s = a^2 / nrm[U]
                        s > top && (top = s; j = U)
                    end
                end
                j == 0 && break
                spent + cost[j] > budget && break
            end
            push!(chosen, j)
            spent += cost[j]
            q = Ap[:, j]
            nq = norm(q)
            if nq > 0
                q ./= nq
                mul!(t, transpose(Ap), q)
                LinearAlgebra.BLAS.ger!(-1.0, q, t, Ap)
                # The lengths by downdating, and afresh every eighth unit, so
                # that rounding cannot leave a spent column looking long.
                if length(chosen) % 8 == 0
                    sum!(abs2, reshape(nrm, 1, N), Ap)
                else
                    @inbounds for U in 1:N
                        nrm[U] = max(nrm[U] - t[U]^2, 0.0)
                    end
                end
            end
            for U in chosen
                nrm[U] = 0.0
            end
            if refit
                G, g = _coreset_gram_grow(G, g, A, chosen, b)
                w = _ctsem_nnls(G, g)
                r .= b .- view(A, :, chosen) * w
            elseif nq > 0
                r .-= q .* dot(q, r)
            end
        end
        isempty(chosen) && return nothing
        ws = if refit
            copy(w)
        else
            As = A[:, chosen]
            _ctsem_nnls(transpose(As) * As, transpose(As) * b)
        end
        v, bound, extra = value(chosen, ws)
        if v > bestv
            best = (units=copy(chosen), weights=ws, extra=extra)
            bestv = v
        end
        # No larger coreset can do better, or there is none.
        (bound <= bestv || length(chosen) < target) && break
    end
    return best
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
