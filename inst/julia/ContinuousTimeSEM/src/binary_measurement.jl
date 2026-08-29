"""
Binary observations, integrated rather than linearised.

# Why not the usual extended-Kalman treatment

ctsem's Stan path handles a binary indicator by moment-matching it to a
Gaussian: predict `p = inv_logit(η)`, take a Jacobian of that, and run the
ordinary Kalman update with `p(1-p)` added to the innovation covariance. It is
the standard extended-filter move and it is biased, because the covariance step
`P -= K H P` uses a linearisation of a saturating function and so understates
how uncertain the state still is. The understatement compounds over time.

Measured on ctsem's Stan path: process noise for a latent seen only through
binary indicators comes back at 0.71 of truth with 5 indicators and 0.46 with
30 -- worse with more data, which is what a systematic bias looks like. Against
a reference fit that samples the states exactly, a drift of -0.935 was reported
as -0.095.

# What this does instead

The Bernoulli likelihood depends on the state only through the scalar linear
predictor `η = λ'x + μ`, and under the predicted Gaussian that is a *scalar*
normal: `η ~ N(λ'x̂ + μ, λ'Pλ)`. So the entire difficulty is one-dimensional,
and Gauss-Hermite integrates it directly:

    Z      = E[ P(y | η) ]                 the observation's marginal likelihood
    E[η|y] = E[ η P(y | η) ] / Z
    V[η|y] = E[ η² P(y | η) ] / Z - E[η|y]²

Those are exact up to quadrature error -- no linearisation anywhere. Lifting
back to the full state needs no approximation either, because `x | η` is
linear-Gaussian:

    c = Pλ,  s² = λ'Pλ
    x̂  ← x̂ + c (E[η|y] - η̂) / s²
    P  ← P - (c c' / s²) (1 - V[η|y] / s²)

The bracketed factor is the whole point: the ordinary update has `1` there,
which is the `V[η|y] → 0` limit. Keeping the term is what stops the filter
believing the state is pinned down.

# Cost

A few Newton steps to locate the scalar mode, then one logistic evaluation per
node and O(n) vector work -- against an O(n²m) matrix update and a Cholesky for
the Gaussian path. Several binary indicators observed at one row are
conditionally independent given the state, so they are applied one at a time,
each an exact scalar update with a Gaussian projection between, which is both
cheaper and more accurate than one joint linearisation.
"""

"""
Nodes for the scalar integration.

Twenty-one, chosen by measurement rather than by taste. Relative error in the
posterior variance against dense numerical integration, with the rule centred
on the posterior mode:

    nodes    s=1      s=2      s=4      s=8
        7   3.4e-5   2.2e-3   2.0e-1   3.3e-1
       21   1.2e-10  3.1e-6   5.5e-3   3.8e-2
       41   1.1e-14  2.8e-9   6.1e-5   1.9e-2

`s` is the predicted standard deviation of the linear predictor. Seven is ample
for the `s <= 2` that ordinary data produces and visibly not enough at `s = 4`,
which happens with a diffuse initial covariance and during optimisation while
the parameters are still poor -- exactly when a bad likelihood does most damage.
Twenty-one costs twenty-one logistic evaluations against a Cholesky for the
Gaussian path, which is not a trade worth economising on.

Large `s` is where this is hardest: the posterior is then dominated by the
likelihood and genuinely skewed, so a mode-centred Gaussian rule has real work
to do. It is also where the answer is least determined, so the error matters
least.
"""
const _CTSEM_BINARY_NODES = Ref(21)

"""Newton steps for the scalar mode. Log-concave, so this converges hard."""
const _CTSEM_BINARY_NEWTON = Ref(6)

"""
    _binary_mode(ηbar, s2, y)

`(mode, curvature)` of `log N(η; ηbar, s²) + log P(y | η)`.

Strictly concave -- the prior contributes `-1/s²` and the Bernoulli
`-p(1-p)`, both negative -- so Newton from the prior mean converges quickly and
cannot diverge. The gradient of the likelihood term is bounded by one, which
bounds the step and keeps this well behaved even when `s` is large and the
observation is nearly deterministic.
"""
@inline function _binary_mode(ηbar::T, s2::T, y::Real) where {T}
    target = T(y > 0.5 ? 1 : 0)
    η = ηbar
    precision = inv(s2)
    curvature = precision
    @inbounds for _ in 1:_CTSEM_BINARY_NEWTON[]
        p = inv(one(T) + exp(-η))
        gradient = -(η - ηbar) * precision + (target - p)
        curvature = precision + p * (one(T) - p)
        η -= -gradient / curvature   # Newton on a concave objective
    end
    return (η, curvature)
end

"""
    _binary_moments(ηbar, s, y, nodes, weights)

`(logZ, mean, variance)` of `η` given one Bernoulli observation.

Adaptive Gauss-Hermite: the rule is centred on the mode of the scalar posterior
and scaled by its curvature, so the integrand it sees is close to the Gaussian
the rule is exact for. A fixed rule centred on the prior instead sees
`inv_logit(ηbar + √2 s t)`, which becomes a step in `t` as `s` grows -- and a
polynomial rule converges slowly on a step, which is precisely where the
measured error was.

`_gauss_hermite` returns physicists' nodes, so `∫f(t)e^{-t²}dt ≈ Σ wᵢf(tᵢ)`;
re-centring turns that into `∫h(η)dη ≈ √2σ̂ Σ wᵢ exp(tᵢ²) h(η̂ + √2σ̂tᵢ)`.
"""
@inline function _binary_moments(ηbar::T, s::T, y::Real, nodes, weights) where {T}
    s2 = s * s
    s2 > zero(T) || return (zero(T), ηbar, zero(T))
    mode, curvature = _binary_mode(ηbar, s2, y)
    scale = sqrt(T(2) / curvature)

    Z = zero(T)
    M1 = zero(T)
    M2 = zero(T)
    observed_one = y > 0.5
    halfprec = inv(T(2) * s2)
    @inbounds for i in eachindex(nodes)
        t = T(nodes[i])
        η = mode + scale * t
        # `p` for a one, `1-p` for a zero, written so neither saturates: at
        # large |η| one underflows to zero and the other to one, and computing
        # the small one directly keeps its logarithm finite.
        likelihood = observed_one ? inv(one(T) + exp(-η)) : inv(one(T) + exp(η))
        deviation = η - ηbar
        # The `exp(t²)` undoes the rule's own kernel; the prior density is then
        # carried explicitly rather than folded into the nodes.
        w = T(weights[i]) * exp(t * t - deviation * deviation * halfprec) *
            likelihood
        Z += w
        M1 += w * η
        M2 += w * η * η
    end
    # A zero means every node put zero probability on the observation, which is
    # a state so far from the data that the row carries no usable information.
    # The caller treats it as an invalid evaluation.
    Z > zero(T) || return (T(-Inf), ηbar, s2)
    mean = M1 / Z
    variance = M2 / Z - mean * mean
    # Z above is √(2π)s times the marginal likelihood: the scale factor and the
    # prior's normalising constant are both outside the sum.
    logZ = log(Z) + log(scale) - log(sqrt(T(2) * T(pi)) * s)
    return (logZ, mean, max(variance, zero(T)))
end

"""
    _ekf_binary_update!(ws, λ, μ, y, n)

One Bernoulli observation, applied exactly in the scalar direction it informs.

Returns the log marginal likelihood of the observation, or `-Inf` if the
predicted state gives it no support.
"""
function _ekf_binary_update!(ws, λ, μ, y::Real, n::Int)
    T = eltype(ws.state)
    nodes, weights = _gauss_hermite(_CTSEM_BINARY_NODES[])

    # c = P λ and s² = λ'Pλ, the predicted mean and variance of η.
    c = view(ws.bufferQ.r, 1:n)
    P = ws.P_predict.data
    s2 = zero(T)
    ηbar = μ
    @inbounds for i in 1:n
        acc = zero(T)
        for j in 1:n
            acc += P[i, j] * λ[j]
        end
        c[i] = acc
        s2 += λ[i] * acc
        ηbar += λ[i] * ws.state[i]
    end
    s2 = max(s2, zero(T))
    s = sqrt(s2)

    logZ, ηpost, vpost = _binary_moments(ηbar, s, y, nodes, weights)
    isfinite(logZ) || return T(-Inf)

    # With no predicted variance in this direction the observation cannot move
    # the state, and the division below would be 0/0. The likelihood still
    # counts.
    s2 > zero(T) || return logZ

    shift = (ηpost - ηbar) / s2
    shrink = (one(T) - vpost / s2) / s2
    @inbounds for i in 1:n
        ws.state[i] += c[i] * shift
    end
    # P ← P - (c c'/s²)(1 - V[η|y]/s²). Written on the lower triangle and
    # mirrored, matching how the rest of the filter keeps its covariances.
    @inbounds for j in 1:n, i in j:n
        ws.P_predict.data[i, j] -= shrink * c[i] * c[j]
    end
    _copy_lower_to_upper!(ws.P_predict.data, ws.state_dim)
    return logZ
end

"""
    _binary_moment_derivatives(ηbar, s2, y)

`(logZ, m, v, dlogZ_da, dlogZ_db, dm_da, dm_db, dv_da, dv_db)` where `a = ηbar`
and `b = s²`.

# Why this differentiates the quadrature rather than the moments

The tilted moments have exact derivatives, obtained by writing the Gaussian's
own derivatives inside the integral: with `d = m - a` and central moments `v`,
`κ₃`, `κ₄`,

    ∂logZ/∂a = d/b            ∂m/∂a = v/b          ∂v/∂a = κ₃/b
    ∂logZ/∂b = (v+d²-b)/(2b²) ∂m/∂b = (κ₃+2dv)/(2b²) ∂v/∂b = (κ₄+2dκ₃-v²)/(2b²)

Those were implemented first and they are correct -- for the *exact* moments.
The forward pass does not return the exact moments, it returns a quadrature
approximation of them, and a reverse pass has to differentiate the function the
forward pass actually computed. Where the quadrature is imperfect the two come
apart: on a one-observation model with a predicted sd near 2.6, the exact-moment
derivative disagreed with a finite difference of the objective by 2%, which is
not an error an optimiser should be asked to work around.

So the rule is differentiated directly, in two dual components over a scalar
loop. That is a small cost -- one extra evaluation of a 21-node logistic sum --
and it is exactly consistent with the forward pass by construction, which is the
property that matters here. Consistency beats elegance: a slightly-wrong
gradient is worse than a slightly-expensive one.
"""
function _binary_moment_derivatives(ηbar::T, s2::T, y::Real) where {T}
    nodes, weights = _gauss_hermite(_CTSEM_BINARY_NODES[])
    if s2 <= zero(T)
        z = zero(T)
        return (z, ηbar, z, z, z, one(T), z, z, one(T))
    end
    triple = function (ab)
        logZ, m, v = _binary_moments(ab[1], sqrt(ab[2]), y, nodes, weights)
        return [logZ, m, v]
    end
    at = [ηbar, s2]
    value = triple(at)
    isfinite(value[1]) || return (T(-Inf), value[2], value[3],
        zero(T), zero(T), one(T), zero(T), zero(T), one(T))
    J = ForwardDiff.jacobian(triple, at)
    return (value[1], value[2], value[3],
        J[1, 1], J[1, 2], J[2, 1], J[2, 2], J[3, 1], J[3, 2])
end

"""
Unconstrained magnitude past which ctsem's transforms are numerically flat.

`log1p_exp(2x)` and friends saturate once `exp(-2|x|)` underflows relative to
one, which is around `|x| = 18`; twenty gives a little room without reaching
into any region a real estimate occupies. Nothing ctsem parameterises has a
meaningful value out there -- a drift of `-1e-9` and a drift of `-1e-15` are
the same model -- so a fit that lands beyond it has stopped for arithmetic
reasons rather than statistical ones.
"""
const _CTSEM_SATURATION = Ref(20.0)
