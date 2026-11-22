## Zero covariances under covmattransform 'z'.
##
## 'z' builds a correlation matrix as exp(A) normalised to a unit diagonal, A
## hollow and symmetric from the cells below the diagonal (cov_expm.jl). A cell
## of A held at 0 is not a zero of exp(A) once other cells connect its two
## variables: with three variables and the other two cells at 0.5 it leaves a
## correlation of 0.115. So a cell the model fixes at zero between connected
## variables (`CTSEMZeroCells`) is not held at 0. It is solved, with the
## others like it, so that exp(A) is zero there, which the normalisation
## keeps.
##
## For any values of the other cells there is exactly one solution. exp(A) is
## the gradient of tr exp(A), which is strictly convex on symmetric matrices
## and grows without bound along every hollow direction, so the closed cells
## are its minimiser given the open ones: Newton on a convex function, in as
## many unknowns as there are such cells. The open cells then reach every
## covariance matrix with those zeros, each once, as plain 'z' reaches every
## covariance matrix.
##
## The values are solved in Float64 through the eigendecomposition of A
## (`_ctsem_symeig!`, the engine's Jacobi), which gives the Newton matrix
## dexp(A)_Z / da_Z in closed form: about z d^3 multiplications an iteration
## for z cells of a d x d matrix, and 3 to 7 iterations. Derivatives are not
## carried through those iterations. At a number type that carries them the
## solved cells take their derivative parts from chord steps on the value's
## Newton matrix, a_Z -= H^-1 exp(A)_Z, each of which makes one more order of
## derivative exact for one exponential at that type. The reverse pass folds
## the solved cells' cotangent into the open cells' through the same matrix
## (`_zcov_fold!`). Everything is cached with the exponential (cov_expm.jl),
## so a matrix is solved once for each value of its cells.
##
## Where the correlation matrix is singular to machine precision the Newton
## matrix does not factor and there is no solution to find: the construction
## then returns NaN, an invalid point, as a non-finite matrix does. That is at
## coordinates of about 10 and beyond, far past the bound on correlations.

# How small the largest correlation left at a solved cell must be, and what a
# solve that stops improving above that may still be accepted at.
const _ZCOV_TOL = 1e-13
const _ZCOV_ACCEPT = 1e-8
const _ZCOV_MAXIT = 60

"""Working storage of the zero solve at one dimension, in the task's scratch."""
mutable struct CTSEMZeroSolve
    d::Int
    eig::CTSEMSymEigWork
    a::Matrix{Float64}
    trial::Matrix{Float64}
    ev::Vector{Float64}
    Psi::Matrix{Float64}
    W::Matrix{Float64}
    TW::Matrix{Float64}
    needed::Vector{Bool}
    H::Matrix{Float64}
    g::Vector{Float64}
    gkeep::Vector{Float64}
    step::Vector{Float64}
end

CTSEMZeroSolve(d::Int) = CTSEMZeroSolve(d, CTSEMSymEigWork(d), zeros(d, d), zeros(d, d),
    zeros(d), zeros(d, d), zeros(d, d), zeros(d, d), fill(false, d), zeros(0, 0),
    Float64[], Float64[], Float64[])

struct _CTSEMZeroSolveKey end

function _zcov_work(d::Int, z::Int)
    table = get!(Dict{Int,CTSEMZeroSolve}, _ctsem_scratch(),
        _CTSEMZeroSolveKey())::Dict{Int,CTSEMZeroSolve}
    w = get(table, d, nothing)
    if w === nothing
        w = CTSEMZeroSolve(d)
        table[d] = w
    end
    if size(w.H, 1) != z
        w.H = zeros(z, z)
        resize!(w.g, z)
        resize!(w.gkeep, z)
        resize!(w.step, z)
    end
    return w
end

"""
How many orders of derivative a number type carries, which is how many chord
steps make a solved cell's derivative parts exact.
"""
_zcov_order(::Type) = 0
_zcov_order(::Type{ForwardDiff.Dual{G,V,N}}) where {G,V,N} = 1 + _zcov_order(V)

# A type whose values are finer than Float64 (BigFloat, in the adjoint's
# precision checks) takes further chord steps to carry the Float64 solution to
# its own precision: each gains about twelve digits.
_zcov_steps(::Type{T}) where {T} =
    _zcov_order(T) + (_primal_type(T) === Float64 ? 0 : 6)

# The eigenpairs of `a` into `w.eig`, `exp(value - top)` into `w.ev`, and
# `exp(a)` at the zero cells, on that scale, into `w.g`. Returns `log tr exp(a)`,
# the largest correlation left at a zero cell, and `tr exp(a)` on that scale.
function _zcov_evaluate!(w::CTSEMZeroSolve, a::Matrix{Float64}, cells::CTSEMZeroCells)
    d = w.d
    e = _ctsem_symeig!(w.eig, a)
    vals, V = e.values, e.vectors
    top = vals[d]
    total = 0.0
    @inbounds for k in 1:d
        w.ev[k] = exp(vals[k] - top)
        total += w.ev[k]
    end
    worst = 0.0
    @inbounds for c in eachindex(cells.rows)
        i = cells.rows[c]; j = cells.cols[c]
        yij = 0.0; yii = 0.0; yjj = 0.0
        for k in 1:d
            ek = w.ev[k]
            yij += V[i, k] * ek * V[j, k]
            yii += V[i, k] * ek * V[i, k]
            yjj += V[j, k] * ek * V[j, k]
        end
        w.g[c] = yij
        r = abs(yij) / sqrt(yii * yjj)
        # Not `max`, which would pass a NaN by.
        (r <= worst) || (worst = r)
    end
    return top + log(total), worst, total
end

# The Newton matrix `d exp(a - shift I)_Z / d a_Z` into `w.H`, from the
# eigenpairs `_zcov_evaluate!` left: the Frechet derivative along a cell's
# direction is `V ((V' E V) .* Psi) V'`, `Psi` the exponential's divided
# differences.
function _zcov_newton_matrix!(w::CTSEMZeroSolve, cells::CTSEMZeroCells, shift::Float64)
    d = w.d
    z = length(cells.rows)
    vals, V = w.eig.values, w.eig.vectors
    Psi, W, TW, H = w.Psi, w.W, w.TW, w.H
    @inbounds for b in 1:d, a in 1:d
        la = vals[a] - shift
        lb = vals[b] - shift
        delta = la - lb
        Psi[a, b] = abs(delta) < 1e-3 ? exp((la + lb) / 2) * (1 + delta * delta / 24) :
            (exp(la) - exp(lb)) / delta
    end
    fill!(w.needed, false)
    @inbounds for c in 1:z
        w.needed[cells.rows[c]] = true
    end
    @inbounds for c in 1:z
        k = cells.rows[c]; l = cells.cols[c]
        for b in 1:d, a in 1:d
            W[a, b] = Psi[a, b] * (V[k, a] * V[l, b] + V[l, a] * V[k, b])
        end
        for i in 1:d
            w.needed[i] || continue
            for b in 1:d
                s = 0.0
                for a in 1:d
                    s += V[i, a] * W[a, b]
                end
                TW[i, b] = s
            end
        end
        for cp in c:z
            i = cells.rows[cp]; j = cells.cols[cp]
            s = 0.0
            for b in 1:d
                s += TW[i, b] * V[j, b]
            end
            H[cp, c] = s
            H[c, cp] = s
        end
    end
    return H
end

# `x <- (U'U)^-1 x` for the upper factor `_ctsem_cholesky!` leaves; `x` may be
# of any number type.
function _zcov_factor_solve!(x::AbstractVector, U::Matrix{Float64}, z::Int)
    @inbounds for i in 1:z
        s = x[i]
        for k in 1:(i - 1)
            s -= U[k, i] * x[k]
        end
        x[i] = s / U[i, i]
    end
    @inbounds for i in z:-1:1
        s = x[i]
        for k in (i + 1):z
            s -= U[i, k] * x[k]
        end
        x[i] = s / U[i, i]
    end
    return x
end

"""
    _zcov_solve!(w, cells) -> Bool

Solve the cells `cells` of the hollow symmetric `w.a` so that `exp(w.a)` is
zero there, by Newton's method on `log tr exp`, which is convex, with a
backtracking line search. `true` with the solution in `w.a` and its eigenpairs
in `w.eig`; `false` where there is none to find in double precision.
"""
function _zcov_solve!(w::CTSEMZeroSolve, cells::CTSEMZeroCells)
    d = w.d
    z = length(cells.rows)
    a, trial = w.a, w.trial
    phi, r, total = _zcov_evaluate!(w, a, cells)
    (isfinite(phi) && isfinite(r)) || return false
    for _ in 1:_ZCOV_MAXIT
        r <= _ZCOV_TOL && return true
        H = _zcov_newton_matrix!(w, cells, w.eig.values[d])
        if !_ctsem_cholesky!(H, z)
            return r <= _ZCOV_ACCEPT
        end
        slope = 0.0
        @inbounds for c in 1:z
            w.step[c] = -w.g[c]
            w.gkeep[c] = w.g[c]
        end
        _zcov_factor_solve!(w.step, H, z)
        @inbounds for c in 1:z
            slope += w.gkeep[c] * w.step[c]
        end
        slope = 2 * slope / total
        t = 1.0
        accepted = false
        phi_t = phi; r_t = r; total_t = total
        while t > 1e-10
            copyto!(trial, a)
            @inbounds for c in 1:z
                i = cells.rows[c]; j = cells.cols[c]
                trial[i, j] += t * w.step[c]
                trial[j, i] += t * w.step[c]
            end
            phi_t, r_t, total_t = _zcov_evaluate!(w, trial, cells)
            # The last term lets a step through once the gain is under the
            # value's roundoff, where Newton is at its solution.
            if isfinite(phi_t) && isfinite(r_t) &&
               phi_t <= phi + 1e-4 * t * slope + 8 * eps(Float64) * (1 + abs(phi))
                accepted = true
                break
            end
            t /= 2
        end
        if !accepted
            # `w.eig` is the last trial's; put the point's own back.
            _zcov_evaluate!(w, a, cells)
            return r <= _ZCOV_ACCEPT
        end
        copyto!(a, trial)
        # At roundoff the residual stops halving: that is the solution.
        stalled = r <= _ZCOV_ACCEPT && r_t > r / 2
        phi, r, total = phi_t, r_t, total_t
        stalled && return true
    end
    return r <= _ZCOV_ACCEPT
end

# The diagonal shift `_fill_hollow!` gives a hollow matrix: the largest row sum
# of absolute values, a bound on its largest eigenvalue.
function _zcov_shift(a::Matrix{Float64}, d::Int)
    shift = 0.0
    @inbounds for i in 1:d
        rowsum = 0.0
        for j in 1:d
            j == i || (rowsum += abs(a[i, j]))
        end
        rowsum > shift && (shift = rowsum)
    end
    return shift
end

"""
    _zcov_solved!(S, Hfac, cells, d, Y, W1, eb, gT) -> Bool

`S` is a hollow symmetric matrix of any number type with its solved cells at
zero. Writes their solution into it and shifts its diagonal as `_fill_hollow!`
does. `Hfac` receives the upper Cholesky factor of the Newton matrix at the
solution, for `exp` of the shifted matrix, which `_zcov_fold!` reads. `Y`,
`W1`, `eb` and `gT` are scratch for the chord steps. `false` where the solve
found no solution.
"""
function _zcov_solved!(S::AbstractMatrix{T}, Hfac::Matrix{Float64}, cells::CTSEMZeroCells,
    d::Int, Y::AbstractMatrix{T}, W1::AbstractMatrix{T}, eb, gT::Vector{T}) where {T}
    z = length(cells.rows)
    w = _zcov_work(d, z)
    @inbounds for j in 1:d, i in 1:d
        w.a[i, j] = i == j ? 0.0 : Float64(_primal(S[i, j]))
    end
    _zcov_solve!(w, cells) || return false
    shift = _zcov_shift(w.a, d)
    _zcov_newton_matrix!(w, cells, shift)
    copyto!(Hfac, w.H)
    _ctsem_cholesky!(Hfac, z) || return false
    @inbounds for c in 1:z
        i = cells.rows[c]; j = cells.cols[c]
        S[i, j] = T(w.a[i, j])
        S[j, i] = T(w.a[i, j])
    end
    steps = _zcov_steps(T)
    if steps > 0
        @inbounds for i in 1:d
            S[i, i] = T(-shift)
        end
        for _ in 1:steps
            my_exp!(Y, S, W1, eb, d)
            @inbounds for c in 1:z
                gT[c] = Y[cells.rows[c], cells.cols[c]]
            end
            _zcov_factor_solve!(gT, Hfac, z)
            @inbounds for c in 1:z
                i = cells.rows[c]; j = cells.cols[c]
                S[i, j] -= gT[c]
                S[j, i] -= gT[c]
            end
        end
    end
    _hollow_shift!(S, d, T)
    return true
end

"""
    _zcov_fold!(L, Abar, A, cells, Hfac, d, Y, Lam, fb, q, lam)

The solved cells are functions of the open ones, so their cotangent belongs to
the open cells. `Abar` is the cotangent of the shifted, solved matrix `A` with
every cell taken as free, a cell's own being `Abar[i,j] + Abar[j,i]`. Writes
into `L` the matrix whose entry `(i, j)` is to be subtracted from open cell
`(i, j)`'s cotangent: with `H` the Newton matrix and `q` the solved cells'
cotangents, `L = Frechet(A)[Lam]`, `Lam` holding `H^-1 q` at the solved cells.
At a type that carries derivatives `H^-1 q` is refined on the value's factor,
once for each order, as the solved cells themselves are.
"""
function _zcov_fold!(L::AbstractMatrix{T}, Abar::AbstractMatrix{T}, A::AbstractMatrix{T},
    cells::CTSEMZeroCells, Hfac::Matrix{Float64}, d::Int, Y::AbstractMatrix{T},
    Lam::AbstractMatrix{T}, fb, q::Vector{T}, lam::Vector{T}) where {T}
    z = length(cells.rows)
    @inbounds for c in 1:z
        i = cells.rows[c]; j = cells.cols[c]
        lam[c] = Abar[i, j] + Abar[j, i]
    end
    _zcov_factor_solve!(lam, Hfac, z)
    steps = _zcov_steps(T)
    for s in 0:steps
        @inbounds for j in 1:d, i in 1:d
            Lam[i, j] = zero(T)
        end
        @inbounds for c in 1:z
            i = cells.rows[c]; j = cells.cols[c]
            Lam[i, j] = lam[c]
            Lam[j, i] = lam[c]
        end
        my_exp_frechet!(Y, L, A, Lam, fb)
        s == steps && break
        # What is left of `H lam = q` at this type, solved on the value's factor.
        @inbounds for c in 1:z
            i = cells.rows[c]; j = cells.cols[c]
            q[c] = Abar[i, j] + Abar[j, i] - (L[i, j] + L[j, i]) / 2
        end
        _zcov_factor_solve!(q, Hfac, z)
        @inbounds for c in 1:z
            lam[c] += q[c]
        end
    end
    return L
end
