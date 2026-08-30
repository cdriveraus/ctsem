"""
Categorical observations, integrated rather than linearised.

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

# Ordinal

Nothing above is Bernoulli-specific except `P(y | η)` itself, so an ordinal
variable is the same machinery with the cumulative logit in that slot:

    P(y <= k | η) = inv_logit(τ_k - η)

Binary is its `K = 2` case with a single threshold at zero, and the two agree
to the last bit -- which is why `_category_likelihood` and `_category_score`
carry both and the binary fast path is an optimisation rather than a separate
model. The reverse pass generalises unchanged apart from one extra cotangent,
for the thresholds themselves.

The same substitution would admit a Poisson or a censored observation. What it
will not admit is nonlinear *dynamics*: those are multivariate in the
prediction step, and none of this applies to them.
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

Ordinal is easier at the same node count, measured the same way over the four
categories of a `τ = (-1.0, 0.4, 1.9)` variable:

    nodes    s=0.5    s=1      s=2      s=4      s=8
        7   8.4e-08  2.1e-05  7.7e-04  1.9e-02  1.3e-01
       21   5.5e-14  7.6e-11  7.1e-08  5.3e-06  3.7e-04
       41   5.3e-14  8.2e-14  3.8e-11  9.2e-08  3.8e-05

An interior category is a bump rather than a step, and a mode-centred rule has
less work to do on a bump. The end categories are the binary case again.

Large `s` is where this is hardest: the posterior is then dominated by the
likelihood and genuinely skewed, so a mode-centred Gaussian rule has real work
to do. It is also where the answer is least determined, so the error matters
least.
"""
const _CTSEM_BINARY_NODES = Ref(21)


"""Newton steps for the scalar mode. Log-concave, so this converges hard."""
const _CTSEM_BINARY_NEWTON = Ref(6)

"""
Predicted variance below which the observation is treated as exact.

The scalar update divides by `s²` once for the mean shift and twice for the
covariance shrink, and the reverse pass divides by `s²` cubed. A variance that
is *positive but tiny* passes a `> 0` test and then overflows those: `b³`
underflows to zero at `b = 6e-103` and `1/b²` overflows at `7.5e-155`, so a
denormal variance produces `Inf - Inf` and a NaN cotangent. It arrives by
ordinary means -- one subject's TI-predictor shift on a variance parameter is
quite enough -- and the NaN then lands on that parameter and every threshold
sharing the row.

Below this the exact update is replaced by its `s² -> 0` limit, which is the
Fisher-scoring form `shift = score`, `shrink = information`: both bounded, both
continuous, and *more* accurate here than the exact expression, whose
cancellation costs it about half a percent by `s² = 1e-12`. The error the
substitution introduces is `O(s² h)`, so 4e-11 relative at the boundary.

Used by the forward update and the reverse pass alike. They have to agree on
where the boundary is, or the adjoint differentiates a function the forward
never computed.
"""
const _CTSEM_MIN_VARIANCE = Ref(1e-100)

"""
The observation kinds the scalar quadrature path handles, matching ctsem's
`manifesttype`: 1 binary, 2 ordinal, 3 count.

The kind travels as an integer beside the thresholds rather than being inferred
from them. Inferring worked while there were two kinds -- no thresholds meant
Bernoulli -- and it is exactly the wrong shape for a third, because a count
also has no thresholds and would have been silently treated as binary. Every
call site passes it explicitly for the same reason: a default here would put
that failure back.
"""
const CTSEM_OBS_BINARY = 1
const CTSEM_OBS_ORDINAL = 2
const CTSEM_OBS_COUNT = 3
const CTSEM_OBS_CENSORED = 4

"""
    _norm_logcdf(z)

`log Phi(z)`, accurate into the tail.

A censored observation at its limit contributes exactly this, and the tail is
where it matters: an observation pinned at a floor while the model predicts
well above it has `z` far negative, and `Phi(z)` underflows to zero long before
the log of it stops being a perfectly ordinary number. The engine's own
`_standard_normal_cdf` is an Abramowitz-Stegun approximation good to about
7.5e-08 *absolute*, which is no relative accuracy at all once `Phi` is 1e-20 --
and a likelihood needs the relative kind.

`logerfc` gives it directly and differentiably. Checked against the asymptotic
expansion: at `z = -100` this returns -5005.524209 where the series gives
-5005.524209, and its first three derivatives are finite there, which the
Laplace path needs because it differentiates this twice more.
"""
@inline _norm_logcdf(z::Real) = logerfc(-z / sqrt(oftype(float(z), 2))) -
    log(oftype(float(z), 2))

"""`log phi(z)`, the standard normal log density."""
@inline _norm_logpdf(z::Real) = -z * z / 2 - log(sqrt(2 * oftype(float(z), pi)))

"""
    _mills(z)

`phi(z) / Phi(z)`, the inverse Mills ratio, as the exponential of a difference
of logarithms rather than a quotient.

Formed directly, the quotient is `0/0` in the tail: both parts underflow at
around `z = -38`, where the ratio itself is a perfectly well behaved `38.03`.
Through the logarithms each part stays of moderate size and the exponential is
taken of their difference, which tends to `log(-z)`.
"""
@inline _mills(z::Real) = exp(_norm_logpdf(z) - _norm_logcdf(z))

"""
    _censored_at(y, limit, upper::Bool)

Whether an observation sits at a censoring limit, decided on values alone.

ForwardDiff breaks comparison ties lexicographically on the partials, so a limit
that happens to carry a derivative seed does not compare equal to the number it
equals. Measured: `5.0 >= 5.0` is true, and `5.0 >= Dual(5.0, [0,1,0])` is
false. That is exactly the case here -- the limits share a vector with the
standard deviation, so they get seeded whenever it does -- and the effect was
that an observation *at* the upper limit silently took the interior branch
during differentiation, contributing a Gaussian density where the forward pass
had contributed a tail probability. The gradient for MANIFESTVAR came out 1%
wrong with two censored observations and 263% wrong with two hundred and
eighty-nine.

Censoring is a property of the data and of limits that are constants, so the
comparison has no business seeing derivative information at all.
"""
@inline _primal(x::Real) = x
@inline _primal(x::ForwardDiff.Dual) = _primal(ForwardDiff.value(x))
@inline _censored_at(y::Real, limit, upper::Bool) =
    upper ? y >= _primal(limit) : y <= _primal(limit)

"""
    _censor_limits(extras, ::Type{T})

The `(lower, upper, sd)` a censored row carries, from the same slot the ordinal
thresholds use.

A censored observation is the only non-Gaussian kind with a measurement error
of its own, so its parameters do not fit in `manifesttype` and `ncategories`
alone. They ride in the extras view because that view already flows everywhere
the kernel does -- through the quadrature, through the mode solve, and into the
adjoint's record -- so the standard deviation is differentiated with everything
else rather than needing a second channel.
"""
@inline function _censor_limits(extras, ::Type{T}) where {T}
    length(extras) >= 3 || return (T(-Inf), T(Inf), one(T))
    # Returned as they are rather than converted to `T`. `T` is the linear
    # predictor's type, and converting to it truncates whenever the extras
    # carry derivative information the predictor does not -- which is exactly
    # the case when the standard deviation is being differentiated. Arithmetic
    # against the predictor promotes correctly on its own.
    return (extras[1], extras[2], extras[3])
end

"""
Largest linear predictor a count observation is allowed to reach.

The Poisson rate is `exp(η)`, which overflows to `Inf` at `η = 710`, and an
infinite rate is worse than a large one: the score becomes `y - Inf` and the
information `Inf`, so Newton's step is `-Inf/Inf` and the mode solve returns
NaN rather than walking back. Clamping the exponent keeps both finite and the
step bounded, so a trial point out here is merely bad rather than poisonous.

Two hundred rather than seven hundred so that the rate squares without
overflowing as well. For real data it never binds: a rate of `exp(20)` is
already half a billion events, and a parameter that reaches this is caught by
the saturation guard long before.
"""
const _CTSEM_COUNT_MAX_LOG_RATE = Ref(200.0)

"""
Hard ceiling on the count-generation walk.

Generating a count inverts its marginal distribution one value at a time, so
the work is linear in the value drawn. A rate large enough for this to bind is
already far outside what these models are for, and an unbounded loop on a bad
parameter draw is worse than a capped one.
"""
const _CTSEM_COUNT_GENERATE_MAX = Ref(100000)

"""
    _log_factorial(y)

`log(y!)`, for the Poisson normalising constant.

Constant in the linear predictor, so it changes no mode, no score and no
posterior moment -- but it *is* part of the log likelihood the fit reports, and
leaving it out would make counts incomparable with every other likelihood in
the package. Computed rather than taken from SpecialFunctions, which the engine
does not otherwise depend on; the observation is data, so this never needs a
derivative.
"""
@inline function _log_factorial(y::Real)
    n = Int(round(y))
    n <= 1 && return 0.0
    if n < 16
        acc = 0.0
        for i in 2:n
            acc += log(i)
        end
        return acc
    end
    # Stirling with the first two correction terms: better than 1e-12 relative
    # from n = 16 up, which is far finer than a constant offset needs.
    x = float(n)
    return 0.5 * log(2 * pi * x) + x * log(x) - x +
        inv(12 * x) - inv(360 * x^3)
end

"""
    _category_score(η, y, thresholds, kind)

`(d log P(y|η)/dη, -d² log P(y|η)/dη²)`, the score and observed information of
one categorical observation.

# Why there is no division here

Write the interval probability as a product rather than a difference. With
`a = τ_{k-1} - η` and `b = τ_k - η`, and `F` the logistic CDF,

    P = F(b) - F(a) = F(-a) F(b) (1 - exp(-(b - a)))

and `b - a` is the *gap between adjacent thresholds*, which is what the model
parameterises. Differentiating the logarithm of that product term by term
gives, with no cancellation and no quotient anywhere,

    d log P / dη  =  F(a) - F(-b)
    -d² log P/dη² =  f(a) + f(b)

both bounded, both continuous as the gap closes. The end categories are the
same formulas with the missing tail set to zero, which is why there is no
branch on `k` beyond fetching the thresholds that exist.

The obvious implementation instead computes `(f(a) - f(b)) / (F(b) - F(a))`,
and that is what this was. It is correct and it is unusable: the denominator is
a difference of two nearly equal numbers, so it loses every significant digit
as the gap closes, and the curvature squares the result. A threshold gap is a
free parameter under a positive transform, so an optimizer reaches tiny gaps
routinely -- measured at a raw value of -7.5, giving a gap of 6e-7, the
gradient came back with eight of twenty-one entries NaN while the objective was
a perfectly ordinary -785. Sixty-eight trial points in one fit were thrown away
for it, the line search ran out of room, and the fit stopped 0.6 log units
short with a gradient of 11.

The *value* is large near a closed gap -- `d log P / dΔ` really is about `1/Δ`
-- but that is a true derivative, not an artefact, and the transform's own
Jacobian cancels it exactly: a gap of `2 log1p_exp(2x)` has `dΔ/dx ≈ 2Δ` when
`Δ` is small, so the raw gradient stays of order one.

Binary with a threshold at zero is the `K = 2` case and gives exactly
`(y - p, p(1-p))`, which is what the Bernoulli fast path computes -- they agree
identically, not just numerically, which is why the fast path can stay.
"""
@inline function _category_score(η::T, y::Real, thresholds,
    kind::Int) where {T}
    if kind == CTSEM_OBS_COUNT
        # `y - λ` and `λ`. Log-concave in η like the others, so the Newton
        # solve above it is unchanged; the clamp is what keeps λ finite.
        λ = exp(min(η, T(_CTSEM_COUNT_MAX_LOG_RATE[])))
        return (T(y) - λ, max(λ, floatmin(T)))
    end
    if kind == CTSEM_OBS_CENSORED
        lower, upper, sd = _censor_limits(thresholds, T)
        prec = inv(sd * sd)
        # For a censored value the score is the inverse Mills ratio and the
        # information is `lambda(z)(z + lambda(z))`, which is the standard
        # truncated-normal result and is non-negative for every `z` -- so the
        # scalar posterior stays log-concave and the Newton solve above still
        # cannot diverge. Uncensored, both collapse to the Gaussian forms.
        if _censored_at(y, lower, false)
            z = (lower - η) / sd
            λ = _mills(z)
            return (-λ / sd, max(λ * (z + λ) * prec, floatmin(T)))
        elseif _censored_at(y, upper, true)
            w = (η - upper) / sd
            λ = _mills(w)
            return (λ / sd, max(λ * (w + λ) * prec, floatmin(T)))
        end
        return ((T(y) - η) * prec, prec)
    end
    if kind == CTSEM_OBS_BINARY || isempty(thresholds)
        p = inv(one(T) + exp(-η))
        return (T(y > 0.5 ? 1 : 0) - p, p * (one(T) - p))
    end
    k = Int(y)
    n = length(thresholds)
    # `F(a)`, zero below the first threshold, and `F(-b)`, zero above the last.
    Fa = k <= 1 ? zero(T) : inv(one(T) + exp(η - thresholds[k - 1]))
    Fnb = k > n ? zero(T) : inv(one(T) + exp(thresholds[k] - η))
    information = Fa * (one(T) - Fa) + Fnb * (one(T) - Fnb)
    return (Fa - Fnb, max(information, floatmin(T)))
end

"""
    _binary_mode(ηbar, s2, y, thresholds, kind)

`(mode - ηbar, curvature)` of `log N(η; ηbar, s²) + log P(y | η)`.

The *offset* rather than the mode, because everything downstream wants the
offset and forming it by subtraction afterwards throws away exactly the digits
that matter. When the prior is tight the mode sits a hair from `ηbar`, and
`mode - ηbar` computed from two numbers of order `ηbar` keeps only the digits
`ηbar` has to spare.

Strictly concave -- the prior contributes `-1/s²` and the observation a
non-positive term, since both the Bernoulli likelihood and a difference of
logistic CDFs are log-concave -- so Newton from the prior mean converges
quickly and cannot diverge. The score is bounded by one in absolute value,
which bounds the step and keeps this well behaved even when `s` is large and
the observation is nearly deterministic.
"""
@inline function _binary_mode(ηbar::T, s2::T, y::Real, thresholds,
    kind::Int) where {T}
    offset = zero(T)
    precision = inv(s2)
    curvature = precision
    @inbounds for _ in 1:_CTSEM_BINARY_NEWTON[]
        score, information = _category_score(ηbar + offset, y, thresholds, kind)
        # `-offset * precision`, not `-(η - ηbar) * precision`: the prior's
        # score is exact this way rather than a difference of two numbers of
        # order ηbar.
        gradient = -offset * precision + score
        curvature = precision + information
        offset += gradient / curvature   # Newton on a concave objective
    end
    return (offset, curvature)
end

"""
    _category_likelihood(η, y, thresholds, kind)

`P(y | η)` for an ordinal observation under the cumulative logit model:

    P(y <= k | η) = inv_logit(τ_k - η)

so category `k` has probability `inv_logit(τ_k - η) - inv_logit(τ_{k-1} - η)`,
with `τ_0 = -Inf` and `τ_K = +Inf`. Binary is the two-category case with a
single threshold at zero, and gives `inv_logit(η)` for a one -- which is why
`thresholds` being empty means binary and needs no separate code path.

`_category_loglikelihood` is the same identity with `log F(z) = -log1p_exp(-z)`
substituted, and it is the one the quadrature uses. A category probability
underflows to zero once the linear predictor is a few hundred away from the
threshold that bounds it, which an optimiser reaches while its parameters are
still poor; the probability form then makes the whole observation impossible
and the row invalid, where the log form is merely a large negative number. On
one 25-subject Laplace fit that difference was 38 of 78 trial points whose
inner mode solve had nothing to work with.

Written as a product of tails rather than as a difference of CDFs: at large
`|η|` one CDF rounds to one and the difference to zero, losing the category's
probability entirely, and as the gap between two thresholds closes the
difference loses its significant digits long before the probability itself
stops being representable. The product form has neither problem -- see
`_category_score`, which needs the same identity for its derivatives.
"""
@inline function _category_loglikelihood(η::T, y::Real, thresholds,
    kind::Int) where {T}
    if kind == CTSEM_OBS_COUNT
        # `y η - λ - log y!`. The clamp matches `_category_score`: the two have
        # to describe the same function or the mode solve chases a likelihood
        # the quadrature is not integrating.
        λ = exp(min(η, T(_CTSEM_COUNT_MAX_LOG_RATE[])))
        return T(y) * η - λ - T(_log_factorial(y))
    end
    if kind == CTSEM_OBS_CENSORED
        lower, upper, sd = _censor_limits(thresholds, T)
        # An observation at or beyond a limit is the probability of being there
        # at all; one inside is the ordinary Gaussian density. Compared with
        # `<=` and `>=` so a value sitting exactly on the limit counts as
        # censored, which is how censored data is recorded.
        if _censored_at(y, lower, false)
            return _norm_logcdf((lower - η) / sd)
        elseif _censored_at(y, upper, true)
            return _norm_logcdf((η - upper) / sd)
        end
        z = (T(y) - η) / sd
        return -z * z / 2 - log(sd) - log(sqrt(2 * T(pi)))
    end
    if kind == CTSEM_OBS_BINARY || isempty(thresholds)
        return y > 0.5 ? -log1p_exp(-η) : -log1p_exp(η)
    end
    k = Int(y)
    n = length(thresholds)
    k <= 1 && return -log1p_exp(η - thresholds[1])
    k > n && return -log1p_exp(thresholds[n] - η)
    gap = thresholds[k] - thresholds[k - 1]
    gap > zero(gap) || return T(-Inf)
    return -log1p_exp(thresholds[k - 1] - η) - log1p_exp(η - thresholds[k]) +
        log(-expm1(-gap))
end

@inline function _category_likelihood(η::T, y::Real, thresholds,
    kind::Int) where {T}
    (kind == CTSEM_OBS_COUNT || kind == CTSEM_OBS_CENSORED) &&
        return exp(_category_loglikelihood(η, y, thresholds, kind))
    (kind == CTSEM_OBS_BINARY || isempty(thresholds)) &&
        return y > 0.5 ? inv(one(T) + exp(-η)) : inv(one(T) + exp(η))
    k = Int(y)
    n = length(thresholds)
    # Below the first threshold, or above the last: one tail, no subtraction.
    k <= 1 && return inv(one(T) + exp(η - thresholds[1]))
    k > n && return inv(one(T) + exp(thresholds[n] - η))
    # F(-a) F(b) (1 - exp(-gap)), not F(b) - F(a): see `_category_score` for
    # what the difference costs. `expm1` keeps the last factor accurate for a
    # gap far below the point where `1 - exp(-gap)` would round to zero.
    gap = thresholds[k] - thresholds[k - 1]
    Fna = inv(one(T) + exp(thresholds[k - 1] - η))
    Fb = inv(one(T) + exp(η - thresholds[k]))
    return Fna * Fb * (-expm1(-gap))
end

"""
    _as_scalar_type(T, thresholds)

The thresholds as element type `T`.

The two derivative helpers each differentiate one of `_binary_moments`' two
kinds of argument while holding the other fixed, so they call it with a dual
predictor and plain thresholds, or the reverse. Resolving that by promotion
does not work: under `intoverpop='laplace'` the fixed side is *already* a dual
carrying the Laplace seed tags, ForwardDiff's `promote_rule` has to order two
unrelated tags to combine them, and it throws rather than guess. Converting
instead needs no ordering -- widening a value into a dual is always defined --
and each caller knows which side is which.

The `Tuple{}` method is what keeps a binary observation free of all of it.
"""
@inline _as_scalar_type(::Type{T}, thresholds) where {T} = T.(thresholds)
@inline _as_scalar_type(::Type{T}, thresholds::Tuple{}) where {T} = thresholds

"""
    _binary_moments(ηbar, s, y, nodes, weights, thresholds, kind)

`(logZ, mean - ηbar, variance)` of `η` given one categorical observation.

# Why the mean comes back as an offset, and the variance about the mode

The filter never wants the posterior mean; it wants `(mean - ηbar)/s²`, and it
wants `(1 - variance/s²)/s²`. Both are differences that vanish as `s²` does,
and computing them from a mean and a variance of order one destroys them.

The variance is the worse of the two. Accumulated as `E[η²] - E[η]²` it is a
difference of two numbers near `ηbar²`, so it carries a relative error of about
`eps * ηbar²/v`; when `v` is `1e-8` and `ηbar` is `0.3` that is one part in
`4e-9`, and `1 - v/s²` -- itself of size `v * h` -- comes out with no correct
digits at all. Measured before this change, on a real fit's numbers: at
`s² = 1e-8` the covariance shrink came out as `1.024` where the answer is
`0.439`, and at `1e-10` as **`-2215`**. A negative shrink makes the filter
*widen* the covariance on an observation, which is not an approximation of
anything.

Accumulating the second moment about the mode instead, and the first as an
offset from `ηbar`, removes both cancellations: the quantities summed are the
small ones to begin with. A near-zero `T0VAR` is an ordinary place for an
optimizer to look, so this is not a corner case -- it was reached on eight of
ten starting draws of a 25-subject fit, and the Laplace route, which
differentiates all of this twice more, turned it into NaN gradients and a
line search with nothing left to accept.

Adaptive Gauss-Hermite: the rule is centred on the mode of the scalar posterior
and scaled by its curvature, so the integrand it sees is close to the Gaussian
the rule is exact for. A fixed rule centred on the prior instead sees
`inv_logit(ηbar + √2 s t)`, which becomes a step in `t` as `s` grows -- and a
polynomial rule converges slowly on a step, which is precisely where the
measured error was.

`_gauss_hermite` returns physicists' nodes, so `∫f(t)e^{-t²}dt ≈ Σ wᵢf(tᵢ)`;
re-centring turns that into `∫h(η)dη ≈ √2σ̂ Σ wᵢ exp(tᵢ²) h(η̂ + √2σ̂tᵢ)`.
"""
@inline function _binary_moments(ηbar::T, s::T, y::Real, nodes, weights,
    thresholds, kind::Int) where {T}
    s2 = s * s
    # A degenerate prior in this direction: `η` is known exactly, so the
    # observation contributes its likelihood *at that point* and moves nothing.
    #
    # Returning `logZ = 0` here instead -- as this did -- says the observation
    # was certain, and that is not a harmless edge case. It makes a vanishing
    # predicted variance *pay*: every categorical observation whose prior
    # variance collapses stops costing anything, so on a model whose T0VAR is
    # free the optimizer is rewarded for driving it to zero, and buys about
    # seventy log units of nothing on twenty-five subjects with two indicators
    # at the first occasion. The objective is also discontinuous there, jumping
    # from `log P(y | η̂)` to `0` the moment the variance underflows, which is a
    # cliff for a line search to fall off rather than a region to search.
    #
    # `log P(y | η̂)` is both the right answer and the continuous limit of the
    # integral, so nothing has to know where the boundary is.
    if !(s2 > T(_CTSEM_MIN_VARIANCE[]))
        return (_category_loglikelihood(ηbar, y, thresholds, kind),
            zero(T), zero(T))
    end
    mode_offset, curvature = _binary_mode(ηbar, s2, y, thresholds, kind)
    scale = sqrt(T(2) / curvature)

    # Accumulated relative to the largest weight seen so far, so the sums are
    # of numbers no larger than one and `Z` is never smaller than one.
    #
    # Unnormalised, the weights are `exp(t² - deviation²/2s²)` times a category
    # probability, and both factors range over many orders of magnitude: a wide
    # prior with a likelihood that confines `η` to one category puts the far
    # nodes' weights near the bottom of the floating point range. `logZ`
    # survives that, because it takes a logarithm; `M1/Z` does not, because the
    # derivative of a quotient divides by `Z²`, which underflows to zero while
    # `Z` is still representable. The result is a finite value with NaN
    # partials -- observed at a predicted variance of 23, where nothing looks
    # extreme at all.
    #
    # Rescaling costs one comparison per node and a multiply on the three
    # accumulators each time the maximum moves, which for a mode-centred rule
    # is a handful of times at most. It is exact algebra, so it changes no
    # value and no derivative.
    # Log-sum-exp over the nodes, rescaled to the largest exponent seen so far.
    #
    # Two underflows are being avoided at once. The likelihood itself vanishes
    # when the linear predictor is far from every threshold, and multiplying it
    # in would make the whole node -- and often the whole observation -- zero;
    # carried as a logarithm it is just a large negative number, and the
    # observation stays usable. And the accumulated sums have to stay away from
    # the bottom of the range because `M1/Z` divides by `Z²` under
    # differentiation, which underflows while `Z` is still representable and
    # leaves a finite value with NaN partials.
    #
    # The rescaling is exact algebra, so it changes no value and no derivative,
    # and it costs a comparison per node plus a multiply on three accumulators
    # each time the maximum moves.
    Z = zero(T)
    M1 = zero(T)
    M2 = zero(T)
    emax = T(-Inf)
    halfprec = inv(T(2) * s2)
    @inbounds for i in eachindex(nodes)
        t = T(nodes[i])
        centred = scale * t              # η - mode
        deviation = mode_offset + centred  # η - ηbar
        # The `t²` undoes the rule's own kernel; the prior density is then
        # carried explicitly rather than folded into the nodes.
        e = t * t - deviation * deviation * halfprec +
            _category_loglikelihood(ηbar + deviation, y, thresholds, kind)
        if e > emax
            ratio = isfinite(emax) ? exp(emax - e) : zero(T)
            Z *= ratio
            M1 *= ratio
            M2 *= ratio
            emax = e
        end
        u = T(weights[i]) * exp(e - emax)
        Z += u
        M1 += u * deviation
        M2 += u * centred * centred
    end
    isfinite(emax) || return (T(-Inf), zero(T), s2)
    # A zero means every node put zero probability on the observation, which is
    # a state so far from the data that the row carries no usable information.
    # The caller treats it as an invalid evaluation.
    Z > zero(T) || return (T(-Inf), zero(T), s2)
    offset = M1 / Z                      # posterior mean - ηbar
    spread = offset - mode_offset        # posterior mean - mode
    variance = M2 / Z - spread * spread
    # Z above is √(2π)s times the marginal likelihood: the scale factor and the
    # prior's normalising constant are both outside the sum.
    logZ = log(Z) + emax + log(scale) - log(sqrt(T(2) * T(pi)) * s)
    return (logZ, offset, max(variance, zero(T)))
end

"""
    _ordinal_thresholds!(ws, pars, row)

The cumulated thresholds for manifest variable `row`, as a view into the
workspace scratch. Empty for a Gaussian or binary variable, which is what makes
the binary path fall through to its own two-outcome branch everywhere.

# Why the matrix holds gaps rather than thresholds

Thresholds must increase, and an optimiser handed `K-1` unconstrained cells
will cross them -- at which point the category between the crossed pair has
probability zero, the likelihood is `-Inf`, and there is no gradient pointing
back out. So THRESHOLDS holds `τ₁` in its first column and the *gap to the
previous threshold* in the rest, each gap passed through a positive transform
on the R side, and this function accumulates them. Ordering then holds by
construction and no constraint has to be enforced anywhere.

The alternative -- writing `τ₁ + exp(δ₂) + ...` into each cell's transform
string -- keeps thresholds in the matrix but makes a transform read several
free parameters, which the adjoint's parameter layer does not support and which
would cost every model a ForwardDiff gradient per cell to support. A running
sum here costs a handful of additions on a vector of length `K-1`.
"""
@inline function _ordinal_thresholds!(ws, pars, row::Int)
    types = ws.manifesttype
    # A censored row carries its limits and its own standard deviation in the
    # same slot. The standard deviation is `MANIFESTVAR`'s diagonal entry,
    # which is a standard deviation before `sdcovsqrt2cov` turns the matrix
    # into a covariance -- and taking it from there rather than from the
    # assembled covariance is what lets the reverse pass hand its cotangent
    # straight back to `MANIFESTVAR` without going through that construction.
    # Legitimate because a censored row is updated on its own, sequentially, so
    # it is never correlated with another row anyway.
    if row <= length(types) && types[row] == 4
        length(ws.thresholds) >= 3 || return view(ws.thresholds, 1:0)
        @inbounds begin
            ws.thresholds[1] = row <= length(ws.censormin) ?
                ws.censormin[row] : -Inf
            ws.thresholds[2] = row <= length(ws.censormax) ?
                ws.censormax[row] : Inf
            ws.thresholds[3] = pars.MANIFESTVAR[row, row]
        end
        return view(ws.thresholds, 1:3)
    end
    hasproperty(pars, :THRESHOLDS) || return view(ws.thresholds, 1:0)
    (row <= length(types) && types[row] == 2) ||
        return view(ws.thresholds, 1:0)
    ncat = row <= length(ws.ncategories) ? ws.ncategories[row] : 0
    k = min(max(ncat - 1, 0), size(pars.THRESHOLDS, 2))
    k == 0 && return view(ws.thresholds, 1:0)
    raw = view(pars.THRESHOLDS, row, :)
    @inbounds begin
        ws.thresholds[1] = raw[1]
        for j in 2:k
            ws.thresholds[j] = ws.thresholds[j - 1] + raw[j]
        end
    end
    return view(ws.thresholds, 1:k)
end

"""
    _ekf_binary_update!(ws, λ, μ, y, n, thresholds, kind)

One categorical observation, applied exactly in the scalar direction it informs.

Returns the log marginal likelihood of the observation, or `-Inf` if the
predicted state gives it no support.
"""
function _ekf_binary_update!(ws, λ, μ, y::Real, n::Int, thresholds,
    kind::Int)
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

    logZ, ηoffset, vpost = _binary_moments(ηbar, s, y, nodes, weights,
        thresholds, kind)
    isfinite(logZ) || return T(-Inf)

    # With no predicted variance in this direction the observation cannot move
    # the state, and the division below would be 0/0. The likelihood still
    # counts.
    s2 > T(_CTSEM_MIN_VARIANCE[]) || return logZ

    shift = ηoffset / s2
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
    _binary_moment_derivatives(ηbar, s2, y, thresholds, kind)

`(logZ, m, v, dlogZ_da, dlogZ_db, dm_da, dm_db, dv_da, dv_db)` where `a = ηbar`
and `b = s²`, and `m` is the posterior mean's *offset* from `a` -- see
`_binary_moments` for why nothing here works with the mean itself.

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
function _binary_moment_derivatives(ηbar::T, s2::T, y::Real,
    thresholds, kind::Int) where {T}
    nodes, weights = _gauss_hermite(_CTSEM_BINARY_NODES[])
    if s2 <= zero(T)
        z = zero(T)
        return (z, z, z, z, z, one(T), z, z, one(T))
    end
    triple = function (ab)
        τ = _as_scalar_type(eltype(ab), thresholds)
        logZ, m, v = _binary_moments(ab[1], sqrt(ab[2]), y, nodes, weights,
            τ, kind)
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
    _binary_threshold_derivatives(ηbar, s2, y, thresholds, kind)

`∂(logZ, m, v)/∂τ` as a `3 x length(thresholds)` matrix.

Separate from `_binary_moment_derivatives` rather than folded into it so the
binary and Gaussian paths pay nothing for ordinal support: this is called only
when a variable actually has thresholds. Differentiating the quadrature rather
than the exact moments, for the same reason given there -- the reverse pass has
to differentiate the function the forward pass computed, not the one it
approximates.
"""
function _binary_threshold_derivatives(ηbar::T, s2::T, y::Real,
    thresholds, kind::Int) where {T}
    k = length(thresholds)
    (k == 0 || s2 <= zero(T)) && return zeros(T, 3, k)
    nodes, weights = _gauss_hermite(_CTSEM_BINARY_NODES[])
    s = sqrt(s2)
    triple = function (τ)
        D = eltype(τ)
        logZ, m, v = _binary_moments(D(ηbar), D(s), y, nodes, weights, τ, kind)
        return [logZ, m, v]
    end
    at = collect(T, thresholds)
    isfinite(triple(at)[1]) || return zeros(T, 3, k)
    return ForwardDiff.jacobian(triple, at)
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

"""
    _standard_normal_cdf(z)

Φ(z), by the Zelen & Severo rational approximation (A&S 26.2.17).

Needed to turn one of the standard normals R supplies into the uniform a
Bernoulli draw wants. `Random` is a dependency and `rand()` would be simpler,
but the whole point of taking the base normals from R is that `set.seed()`
governs generation; drawing here would put half the randomness outside the
user's control.

Absolute error below 8e-8, which is nothing against a Bernoulli threshold, and
this is generation rather than likelihood so it is never differentiated.
"""
@inline function _standard_normal_cdf(z::Real)
    t = inv(1 + 0.2316419 * abs(z))
    poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 +
        t * (-1.821255978 + t * 1.330274429))))
    tail = exp(-z * z / 2) / sqrt(2 * pi) * poly
    return z >= 0 ? 1 - tail : tail
end
