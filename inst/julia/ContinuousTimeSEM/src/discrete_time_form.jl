using ComponentArrays
using LinearAlgebra
using ForwardDiff

################################################################################
# Discrete-time EKF log-likelihood and gradient function
################################################################################

"""
    DiscretizationCache{T}

Last-value cache for the two expensive, row-invariant pieces of
`_compute_discrete_time_form!`.

Both are guarded by an `O(n^2)` comparison against the inputs that produced
them, so a miss costs almost nothing and a hit skips an `O(n^3)` factorization:

  * **`exp`**: `eJAx = exp(JAx * Δt)` is reused when both `JAx` and `Δt` are
    unchanged. That covers the common panel design where every subject is
    measured on the same wave schedule -- and correctly misses on genuinely
    irregular observation times.
  * **`lyap`**: the asymptotic-diffusion solve `X` depends on `JAx` and the
    diffusion covariance and **not on `Δt` at all**, so for any model whose
    drift and diffusion are not state-dependent it is the same at every row of
    every subject. This is the important one: `schur!` was 38% of the primal
    filter's runtime and 60 of its 68 MB of allocations on a 20-latent model.

A state-dependent model changes `JAx` per row, so it simply misses every time
and pays only the comparison. Nothing here assumes linearity; the guard checks
the actual inputs rather than trusting a model-level flag.
"""
mutable struct DiscretizationCache{T}
    exp_JAx::Matrix{T}
    exp_dt::T
    exp_out::Matrix{T}
    exp_valid::Bool
    lyap_JAx::Matrix{T}
    lyap_Q::Matrix{T}
    lyap_out::Matrix{T}
    lyap_valid::Bool
end

function DiscretizationCache(::Type{T}, n::Int, k::Int) where {T}
    return DiscretizationCache{T}(
        zeros(T, n, n), zero(T), zeros(T, n, n), false,
        zeros(T, k, k), zeros(T, k, k), zeros(T, k, k), false)
end

"""
    _blocks_identical(a, b, d1, d2)

Whether the leading `d1 x d2` blocks of `a` and `b` are identical.

Plain `==` is the right comparison here **including for `ForwardDiff.Dual`
entries**, which is worth stating because it is easy to assume otherwise:
ForwardDiff defines `==` on same-tag duals as
`==(value(x), value(y)) && ==(partials(x), partials(y))` (see `dual.jl`), so
two duals agreeing in value but not in derivative correctly compare unequal.
If that were value-only, these caches would reuse a matrix exponential whose
derivative belongs to a different point -- a silently wrong gradient on the
default ForwardDiff path. `test_adjoint_primitives.jl` pins that ForwardDiff
behaviour, so a change to it fails there rather than here.
"""
@inline function _blocks_identical(a, b, d1::Int, d2::Int)
    @inbounds for j in 1:d2, i in 1:d1
        a[i, j] == b[i, j] || return false
    end
    return true
end

@inline function _copy_block!(dest, src, d1::Int, d2::Int)
    @inbounds for j in 1:d2, i in 1:d1
        dest[i, j] = src[i, j]
    end
    return dest
end

# `nothing` disables caching, for the convenience constructor below and any
# caller that has no workspace to hang a cache on.
@inline _exp_cache_hit(::Nothing, JAx, Δt, n::Int) = false
@inline _exp_cache_store!(::Nothing, JAx, Δt, out, n::Int) = nothing
@inline _lyap_cache_hit(::Nothing, JAx, Q, k::Int) = false
@inline _lyap_cache_store!(::Nothing, JAx, Q, out, k::Int) = nothing

@inline function _exp_cache_hit(cache::DiscretizationCache, JAx, Δt, n::Int)
    return cache.exp_valid && cache.exp_dt == Δt && _blocks_identical(cache.exp_JAx, JAx, n, n)
end

@inline function _exp_cache_store!(cache::DiscretizationCache, JAx, Δt, out, n::Int)
    _copy_block!(cache.exp_JAx, JAx, n, n)
    cache.exp_dt = Δt
    _copy_block!(cache.exp_out, out, n, n)
    cache.exp_valid = true
    return nothing
end

@inline function _lyap_cache_hit(cache::DiscretizationCache, JAx, Q, k::Int)
    return cache.lyap_valid && _blocks_identical(cache.lyap_JAx, JAx, k, k) &&
        _blocks_identical(cache.lyap_Q, Q, k, k)
end

@inline function _lyap_cache_store!(cache::DiscretizationCache, JAx, Q, out, k::Int)
    _copy_block!(cache.lyap_JAx, JAx, k, k)
    _copy_block!(cache.lyap_Q, Q, k, k)
    _copy_block!(cache.lyap_out, out, k, k)
    cache.lyap_valid = true
    return nothing
end

"""
    _compute_discrete_time_form!(discrete_ca, buffer, DIFFUSIONcov, pars, Δt, exp_buffer, lyap_buffer)

Compute continuous-to-discrete dynamics for one EKF time interval.

The function writes the discrete drift, intercept, and process covariance into
`discrete_ca` using the supplied scratch buffers. `DIFFUSIONcov` is the
continuous-time process covariance and `Δt` is the elapsed time.
"""
function _compute_discrete_time_form!(discrete_ca, buffer, DIFFUSIONcov, pars, Δt, exp_buffer, lyap_buffer)
    n = size(DIFFUSIONcov, 1)
    diffusion_state_indices = collect(1:n)
    diffusion_buffer = _make_square_buffer(eltype(DIFFUSIONcov), n)
    discretization_buffer = _make_discretization_buffer(eltype(DIFFUSIONcov), n)
    return _compute_discrete_time_form!(
        discrete_ca,
        buffer,
        DIFFUSIONcov,
        pars,
        Δt,
        exp_buffer,
        lyap_buffer,
        zeros(eltype(DIFFUSIONcov), n),
        diffusion_state_indices,
        diffusion_buffer,
        discretization_buffer,
        exp_buffer.dim,
        nothing,
    )
end

function _compute_discrete_time_form!(discrete_ca, buffer, DIFFUSIONcov, pars, Δt, exp_buffer,
    lyap_buffer, state, diffusion_state_indices, diffusion_buffer,
    discretization_buffer, dim::Val{d}, cache=nothing) where {d}
    # Stan propagates the local affine EKF model, not the raw DRIFT matrix:
    # f(x) = DRIFT * x + CINT, J = JAx, c = f(x) - J * x.
    if _exp_cache_hit(cache, pars.JAx, Δt, d)
        _copy_block!(discrete_ca.eJAx, cache.exp_out, d, d)
    else
        copyto!(exp_buffer.As, pars.JAx)
        rmul!(exp_buffer.As, Δt)
        my_exp!(discrete_ca.eJAx, exp_buffer.As, buffer.intermediate, exp_buffer, dim)
        _exp_cache_store!(cache, pars.JAx, Δt, discrete_ca.eJAx, d)
    end
    copyto!(discrete_ca.dDRIFT, discrete_ca.eJAx)

    # The Lyapunov solve is restricted to the genuine diffusion block. Static
    # augmented coordinates have no diffusion and are excluded, as in Stan's
    # derrind subset; their covariance still propagates through full eJAx.
    fill!(discrete_ca.dDIFFUSION, zero(eltype(discrete_ca.dDIFFUSION)))
    kdim = length(diffusion_state_indices)
    fill!(discrete_ca.dINT, zero(eltype(discrete_ca.dINT)))

    # Stan forms the local affine correction only in the dynamic state block.
    # The full eJAx above still supplies the dynamic-from-static contribution;
    # keeping this solve in the dynamic block avoids inverting the singular
    # augmented Jacobian induced by static random-effect coordinates.
    @inbounds for i in 1:kdim
        ii = diffusion_state_indices[i]
        affine = pars.CINT[ii]
        for j in 1:d
            affine += (pars.DRIFT[ii, j] - pars.JAx[ii, j]) * state[j]
        end
        diffusion_buffer.r[i] = affine
    end
    @inbounds for j in 1:kdim, i in 1:kdim
        ii = diffusion_state_indices[i]
        jj = diffusion_state_indices[j]
        diffusion_buffer.intermediate[i, j] = pars.JAx[ii, jj]
    end
    @inbounds for i in 1:kdim
        ii = diffusion_state_indices[i]
        correction = -diffusion_buffer.r[i]
        for q in 1:kdim
            qq = diffusion_state_indices[q]
            correction += discrete_ca.eJAx[ii, qq] * diffusion_buffer.r[q]
        end
        diffusion_buffer.s[i] = correction
    end
    _solve_square_system!(diffusion_buffer.intermediate, diffusion_buffer.s,
        view(exp_buffer.piv, 1:kdim), Val(kdim))
    @inbounds for i in 1:kdim
        discrete_ca.dINT[diffusion_state_indices[i]] = diffusion_buffer.s[i]
    end

    @inbounds for j in 1:kdim, i in 1:kdim
        ii = diffusion_state_indices[i]
        jj = diffusion_state_indices[j]
        discretization_buffer.input[i, j] = pars.JAx[ii, jj]
        diffusion_buffer.intermediate[i, j] = DIFFUSIONcov[ii, jj]
    end
    dynamic_jacobian = view(discretization_buffer.input, 1:kdim, 1:kdim)
    # `diffusion_buffer.intermediate` holds the diffusion covariance here and
    # is overwritten immediately below, so the cache has to copy it now.
    if _lyap_cache_hit(cache, dynamic_jacobian, diffusion_buffer.intermediate, kdim)
        _copy_block!(diffusion_buffer.out, cache.lyap_out, kdim, kdim)
    else
        my_lyap!(diffusion_buffer.out, dynamic_jacobian, diffusion_buffer.intermediate, lyap_buffer)
        _lyap_cache_store!(cache, dynamic_jacobian, diffusion_buffer.intermediate,
            diffusion_buffer.out, kdim)
    end
    @inbounds for r in 1:kdim, i in 1:kdim
        value = zero(eltype(diffusion_buffer.intermediate))
        for q in 1:kdim
            value += discrete_ca.eJAx[diffusion_state_indices[i], diffusion_state_indices[q]] *
                diffusion_buffer.out[q, r]
        end
        diffusion_buffer.intermediate[i, r] = value
    end
    @inbounds for j in 1:kdim, i in 1:kdim
        value = diffusion_buffer.out[i, j]
        for r in 1:kdim
            value -= diffusion_buffer.intermediate[i, r] *
                discrete_ca.eJAx[diffusion_state_indices[j], diffusion_state_indices[r]]
        end
        discrete_ca.dDIFFUSION[diffusion_state_indices[i], diffusion_state_indices[j]] = value
    end
    return 
end

"""
    _compute_one_step_form!(discrete_ca, DIFFUSIONcov, pars, state,
                            diffusion_state_indices, dim)

The discrete-time case, where there is nothing to discretize.

A discrete-time model's DRIFT, CINT and DIFFUSION are already the one-step
quantities, so the three expensive pieces of the continuous form collapse: the
transition is `JAx` itself, the process noise is the diffusion covariance, and
the intercept is the local affine offset with no solve around it. `Δt` plays no
part -- a discrete model advances one step per row whatever the recorded
interval, which is also what Stan does.

The affine offset is formed over every state rather than only the diffusing
ones: with no solve to keep away from the singular augmented block there is no
reason to restrict it, and it is zero on the static states in any case.
"""
function _compute_one_step_form!(discrete_ca, DIFFUSIONcov, pars, state,
    diffusion_state_indices, dim::Val{d}) where {d}
    copyto!(discrete_ca.eJAx, pars.JAx)
    copyto!(discrete_ca.dDRIFT, pars.JAx)

    @inbounds for i in 1:d
        affine = pars.CINT[i]
        for j in 1:d
            affine += (pars.DRIFT[i, j] - pars.JAx[i, j]) * state[j]
        end
        discrete_ca.dINT[i] = affine
    end

    fill!(discrete_ca.dDIFFUSION, zero(eltype(discrete_ca.dDIFFUSION)))
    kdim = length(diffusion_state_indices)
    @inbounds for j in 1:kdim, i in 1:kdim
        ii = diffusion_state_indices[i]
        jj = diffusion_state_indices[j]
        discrete_ca.dDIFFUSION[ii, jj] = DIFFUSIONcov[ii, jj]
    end
    return discrete_ca
end
