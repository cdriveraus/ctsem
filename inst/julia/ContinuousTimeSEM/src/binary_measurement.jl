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


"""
Iterations the scalar mode solve may take before giving up.

One budget for every kind, where there used to be six for a logistic and
fourteen for a count. The two numbers existed because the walk was undamped and
capped, so the budget had to cover the whole distance in fixed-size steps and
the distance differs by kind. With a line search the iteration is globally
convergent on this objective and stops on its own as soon as the step it would
take is worth less than the tolerance, so this is a bound against pathology
rather than a tuning parameter: it is reached only if something is wrong, and
raising it does not buy accuracy in any case that was already converging.
"""
const _CTSEM_MODE_MAXITER = Ref(30)

"""
Halvings the mode solve's line search may try before it concludes the step is
not an improvement at all.

Twenty, as the Laplace inner solve uses, which is the same rule applied to the
same kind of objective: `2^-20` of a Newton step is far below the point where
an accepted step could still matter, so exhausting them means the objective is
not improvable in that direction rather than that the step was too long.
"""
const _CTSEM_MODE_BACKTRACKS = 20

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
    _count_dispersion(extras, ::Type{T})

The log-scale dispersion a count row carries, from the same slot the ordinal
thresholds and the censored limits use. Zero when the row carries nothing,
which is the equidispersed Poisson.

# Why this one does not enter the likelihood

Every other extra is a parameter of `P(y | eta)`: a threshold moves the category
boundaries, a censoring limit moves where the instrument saturates. This one is
not. A count's dispersion is a Gaussian term *added to* the linear predictor, so
`y | x` is a Poisson-lognormal mixture and the mixing is over exactly the scalar
the quadrature already integrates:

    eta = lam'x + mu + e,  e ~ N(0, s2)   =>   eta | x ~ N(lam'xhat + mu, lam'P lam + s2)

So it joins the variance the rule integrates over, and `_category_loglikelihood`
never sees it. That is also why the state update needs no new algebra:
`cov(x, eta)` is still `P lam`, so the projection that lifts the scalar posterior
back to the state holds with `lam'P lam + s2` in place of `lam'P lam` -- see
`_ekf_binary_update!`.

# Why a count gets one and a binary does not

Not a preference: for a binary or an ordinal indicator this parameter is not
identified. With a probit link the mixture is exact,

    E_e Phi(lam'x + mu + e) = Phi((lam'x + mu) / sqrt(1 + s2))

so the dispersion is absorbed into the loadings and thresholds and nothing in
the data separates it from them; the logistic link differs only in that the
absorption is approximate. A count has no such freedom, because the Poisson's
variance is locked to its mean: the dispersion shifts the mean by `s2/2`, which
MANIFESTMEANS absorbs, and multiplies the variance by a factor nothing else can
produce. So it is identified, and it is the only parameter in the model that
moves the variance-to-mean ratio.

Returned raw rather than converted to `T`, for the reason `_censor_limits`
gives: `T` is the predictor's type and converting truncates the derivative
information the dispersion carries and the predictor does not.
"""
@inline function _count_dispersion(extras, ::Type{T}) where {T}
    isempty(extras) && return zero(T)
    return extras[1]
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

It is also what keeps a *generated* count finite when a parameter draw puts the
predictor out where `exp(η)` is not: a saturated draw is then a large finite
count, which plots and summarises as the nonsense it is, rather than an `Inf`
that takes the rest of the row's likelihood with it.
"""
const _CTSEM_COUNT_MAX_LOG_RATE = Ref(200.0)

"""
Hard ceiling on the count-generation walk.

Generating a count inverts its distribution one value at a time, so the work is
linear in the value drawn. `_ctsem_draw_count` hands over to a normal
approximation at `_CTSEM_POISSON_NORMAL_RATE`, so its walk is already bounded
by a rate of five hundred; this is the ceiling for the case that bound does not
cover, a `u` in the last representable sliver below one.
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

Kept in floating point throughout rather than counting in `Int`. A count is
unbounded, so a value past `typemax(Int64)` is reachable -- from generation
under a saturated rate, or from data someone hands us -- and there it is
Stirling's series that is wanted, evaluated on the float. Converting first
turned that into an `InexactError` from inside the likelihood, which named a
number and nothing else.
"""
@inline function _log_factorial(y::Real)
    x = float(round(y))
    x <= 1 && return 0.0
    if x < 16
        acc = 0.0
        for i in 2:Int(x)
            acc += log(i)
        end
        return acc
    end
    # Stirling with the first two correction terms: better than 1e-12 relative
    # from n = 16 up, which is far finer than a constant offset needs.
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
        F = inv(one(T) + exp(-η))
        length(thresholds) >= 2 ||
            return (T(y > 0.5 ? 1 : 0) - F, F * (one(T) - F))
        # With asymptotes the likelihood is `q = flat + (d-c)F` for a one and
        # `flat + (d-c)(1-F)` for a zero. Writing `u` for `dq/dη` in absolute
        # value, the score is `+/-u/q` and the information is
        # `(u^2 - q u')/q^2`. At `c = 0, d = 1` this collapses to the pair
        # above, which is the check worth keeping in mind when reading it.
        #
        # This is only reached where the predicted variance has collapsed, and
        # from the adjoint's matching branch. The mode solve never sees it: the
        # asymptote route integrates through `_asymptote_moments`, which hands
        # the quadrature a plain binary row.
        c, d = _binary_asymptotes(thresholds, T)
        flat = y > 0.5 ? c : one(c) - d
        span = d - c
        q = flat + span * (y > 0.5 ? F : one(T) - F)
        u = span * F * (one(T) - F)
        du = span * F * (one(T) - F) * (one(T) - 2 * F)
        # `dq/dη` is `u` for a one and `-u` for a zero; `d^2q/dη^2` is `du`
        # and `-du`. Both signs cancel in the information.
        score = (y > 0.5 ? u : -u) / q
        # Returned as it is, and it is genuinely negative in places: a three
        # parameter logistic is not log-concave for a correct response, so the
        # curvature of its log likelihood changes sign. The other kinds clamp
        # at `floatmin` because theirs cannot, and clamping here would report
        # a convex region as a flat one -- and hand any Newton step that
        # divided by it something of order 1e308.
        #
        # Nothing takes such a step. `_binary_moments` sends an item with
        # asymptotes to `_asymptote_moments`, which integrates it as a mixture
        # precisely so that no mode solve ever meets this surface, and the
        # only other caller is the adjoint's degenerate-variance branch, which
        # uses the score and discards this.
        information = (u * u - q * (y > 0.5 ? du : -du)) / (q * q)
        return (score, information)
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
    _logaddexp(a, b)

`log(exp(a) + exp(b))` without forming either exponential.

One helper rather than the same three lines at each site: the asymptote
likelihood and the asymptote moments both add a component in log space, and
both have a component that legitimately underflows -- a guessing probability of
zero, or a two parameter logistic term whose mass is far into the tail.
"""
@inline function _logaddexp(a::Real, b::Real)
    isfinite(a) || return b
    isfinite(b) || return a
    m = max(a, b)
    m + log(exp(a - m) + exp(b - m))
end

"""
    _binary_asymptotes(extras, T)

A binary item's lower and upper asymptotes, `(c, d)`, or `(0, 1)` when it has
none.

Presence is decided by the length of the extras, never by comparing the values
to zero and one. A free asymptote sitting at zero is an ordinary place for an
optimizer to be, and `iszero` on a dual number tests the value while the
partials are the thing that would be lost -- the same tie-breaking hazard
`_censored_at` is written to avoid.

The pair arrives already accumulated, as ordinal thresholds do: the matrix
holds `c` and then a gap, each a single free parameter with its own bounded
transform, and `_ordinal_thresholds!` turns the gap into `d = c + (1-c)g`. That
keeps `0 <= c < d <= 1` without any cell's transform having to read another
cell's parameter, which is the constraint the parameter layer imposes.
"""
@inline function _binary_asymptotes(extras, ::Type{T}) where {T}
    length(extras) >= 2 || return (zero(T), one(T))
    return (extras[1], extras[2])
end

"""
    _censored_moments(ηbar, s, y, thresholds)

`(logZ, mean - ηbar, variance)` for a censored observation, in closed form.

A censored row is a Gaussian observation of a Gaussian state, so everything
the filter wants about it can be written down. The other kinds need quadrature
because a logistic or a Poisson likelihood against a Gaussian prior has no
elementary integral; this one does, and integrating it numerically was only
ever costing accuracy. Measured against the closed form, the 21 node rule was
1.1e-02 out at a predicted sd of 5 and 7.4e-02 at 20 -- the worst of the four
kinds, on the one where no approximation was needed at all.

Write `V = sd^2 + s^2` for the variance of the observation `y* = η + ε`, which
is what is actually observed or known to lie beyond a limit.

*Inside the limits* the observation is `y` itself. Its density is
`N(y; ηbar, sqrt(V))`, and the posterior is the usual precision-weighted
combination: the mean moves `s^2/V` of the way from the prior to the
observation and the variance is `s^2 sd^2 / V`.

*At a limit* only `y* <= lower` (or `y* >= upper`) is known. Then `(η, y*)` is
jointly Gaussian with covariance `s^2`, so conditioning on a one-sided event in
`y*` gives the truncated-normal moments through that covariance. With `b` the
standardised distance to the limit, on the side the mass is, and `lambda` the
inverse Mills ratio at `b`:

    logZ      = log Phi(b)
    mean      = ηbar -/+ (s^2 / sqrt(V)) lambda
    variance  = s^2 - (s^4 / V) (b lambda + lambda^2)

`b lambda + lambda^2` is one minus the truncated variance of a standard normal,
and it is the one place here that cancels: both terms grow like `b^2` while
their sum tends to one. At `b = -40` that is three digits of sixteen, which is
where deep censoring lives and is comfortably enough.

Nothing special is needed for a zero measurement standard deviation, which the
formulas reduce to a noiseless observation, or for an infinite limit on one
side, which `_censored_at` never reports as censored.
"""
@inline function _censored_moments(ηbar::T, s::T, y::Real, thresholds) where {T}
    lower, upper, sd = _censor_limits(thresholds, T)
    s2 = s * s
    V = sd * sd + s2
    rootV = sqrt(V)
    atlower = _censored_at(y, lower, false)
    atupper = _censored_at(y, upper, true)
    if !atlower && !atupper
        deviation = y - ηbar
        z = deviation / rootV
        logZ = -z * z / 2 - log(rootV) - log(sqrt(2 * T(pi)))
        return (logZ, deviation * s2 / V, s2 * sd * sd / V)
    end
    # The standardised distance to the limit, signed so that the mass the
    # observation reports always lies below `b`. The upper case is the lower
    # one in `-η`, which flips the mean shift and leaves the variance alone.
    b = atlower ? (lower - ηbar) / rootV : (ηbar - upper) / rootV
    logZ = _norm_logcdf(b)
    lambda = _mills(b)
    shift = (s2 / rootV) * lambda
    offset = atlower ? -shift : shift
    variance = s2 - (s2 * s2 / V) * (b * lambda + lambda * lambda)
    return (logZ, offset, max(variance, zero(variance)))
end

"""
    _mode_start(ηbar, y, thresholds, kind)

Where the mode iteration begins, as an offset from `ηbar`: the prior mean moved
to where this observation's own likelihood puts its mass.

Starting at `ηbar` is starting wherever the prior happens to be, which for an
observation the predicted state makes improbable is arbitrarily far from the
answer. The line search will still get there, but it gets there by halving, and
the distance it has to cover is then set by the prior rather than by the
likelihood. Projecting first does that part in closed form and leaves a
distance the logistic sets, which is a few units whatever `s` is.

Every kind that reaches the iteration, in one place. The count case has had
such a start since its own mode solve was fixed, and ordinal and binary got one
when the same failure was found there. Censored is absent because it does not
come here at all: it is solved in closed form by `_censored_moments`.

Written as branches on the observation rather than as a `clamp` against
infinities, so that no infinity meets a dual number: an infinite partial times
a zero is the NaN that travels silently.
"""
@inline function _mode_start(ηbar::T, y::Real, thresholds, kind::Int) where {T}
    if kind == CTSEM_OBS_COUNT
        # The likelihood's own mode, `log(y + 1/2)`. The half keeps `y = 0`
        # finite, which is the commonest observation in floor-heavy count data
        # and the corner a bare `log y` start sent back to `ηbar`.
        return log(T(y) + T(0.5)) - ηbar
    end
    if kind == CTSEM_OBS_ORDINAL && !isempty(thresholds)
        n = length(thresholds)
        k = Int(y)
        k <= 1 && return min(zero(T), thresholds[1] - ηbar)
        k > n && return max(zero(T), thresholds[n] - ηbar)
        lo = thresholds[k - 1] - ηbar
        hi = thresholds[k] - ηbar
        return min(max(zero(T), lo), hi)
    end
    # Binary, whose cut is at zero: the band is above it for a one and below it
    # for a zero.
    return y > 0.5 ? max(zero(T), -ηbar) : min(zero(T), -ηbar)
end

"""
    _mode_objective(ηbar, precision, offset, y, thresholds, kind)

`log N(η; ηbar, s²) + log P(y | η)` at `η = ηbar + offset`, up to a constant.

What the line search compares. The prior term is written from the offset rather
than from `η - ηbar`, for the reason `_binary_mode` writes its gradient that
way: the two differ by a cancellation when the prior is tight.
"""
@inline _mode_objective(ηbar::T, precision::T, offset::T, y::Real, thresholds,
    kind::Int) where {T} =
    -offset * offset * precision / 2 +
        _category_loglikelihood(ηbar + offset, y, thresholds, kind)

"""
    _newton_gain(gradient, step)

How much objective one Newton step is predicted to gain: `g² / 2M`, which for a
`step` of `g / M` is `g * step / 2`.

The scalar twin of `_laplace_stationary_gain`, deliberately the same quantity
written the same way, so that the two mode solves in this package declare
convergence by one rule rather than by two. The argument for preferring it to a
bound on the gradient is written out at `_laplace_inner_tolerance`: a gradient
bound asks for a number of digits that depends on how large the objective
happens to be, while the predicted gain is in the units of the value and
survives a reparameterisation of the thing being solved for.
"""
@inline function _newton_gain(gradient::T, step::T) where {T}
    gain = gradient * step / 2
    isfinite(gain) ? abs(gain) : T(Inf)
end

"""
    _mode_tolerance(value)

How much predicted gain the mode solve may leave unclaimed.

Relative to the value with an absolute floor, as `_laplace_inner_tolerance` is
and for the same reason. The mode places a quadrature rule rather than being
reported, and the rule's own error at these node counts is far above this when
the mode is right, so asking for more is asking for precision the answer does
not carry.
"""
@inline _mode_tolerance(value::T) where {T} = T(1e-12) * (one(T) + abs(value))

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
logistic CDFs are log-concave. Concavity makes the mode unique; it does not on
its own make undamped Newton find it, and the argument that used to stand here
-- that a score bounded by one bounds the step -- had the implication
backwards. A Newton step is `gradient / curvature`, so a bounded score with a
*vanishing* information is the dangerous combination, not a safe one: it is
what makes the step `±s²`. What bounds it is the line search below, which
takes a step only when it improves the objective; `_mode_start` then shortens
the journey rather than making it safe.
"""
@inline function _binary_mode(ηbar::T, s2::T, y::Real, thresholds,
    kind::Int) where {T}
    precision = inv(s2)
    offset = _mode_start(ηbar, y, thresholds, kind)
    value = _mode_objective(ηbar, precision, offset, y, thresholds, kind)
    score, information = _category_score(ηbar + offset, y, thresholds, kind)
    curvature = precision + information
    @inbounds for _ in 1:_CTSEM_MODE_MAXITER[]
        # `-offset * precision`, not `-(η - ηbar) * precision`: the prior's
        # score is exact this way rather than a difference of two numbers of
        # order ηbar.
        gradient = -offset * precision + score
        step = gradient / curvature
        _newton_gain(gradient, step) <= _mode_tolerance(value) && break
        # Concavity makes the Newton direction an ascent direction; it does not
        # make the Newton *step* an improvement, and that distinction is the
        # whole history of this function. An undamped step is
        # `gradient / curvature`, so wherever the likelihood contributes no
        # curvature -- an over-predicted zero for a count, a category the
        # predicted state makes improbable for an ordinal or a binary -- the
        # only curvature left is the prior's `1/s²` and the step is the
        # gradient times `s²`. That overshoots, the score flips sign, and the
        # iteration oscillates rather than converging. It was patched twice,
        # with a per-kind cap and a per-kind iteration budget, before a line
        # search replaced both: a step is taken when it improves the objective
        # and halved when it does not, which needs nothing chosen per kind and
        # cannot oscillate.
        scale = one(T)
        accepted = false
        for _ in 1:_CTSEM_MODE_BACKTRACKS
            candidate = offset + scale * step
            trial = _mode_objective(ηbar, precision, candidate, y, thresholds,
                kind)
            if isfinite(trial) && trial >= value
                offset = candidate
                value = trial
                score, information = _category_score(ηbar + offset, y,
                    thresholds, kind)
                curvature = precision + information
                accepted = true
                break
            end
            scale /= 2
        end
        # No scale of an ascent direction improved the value, so what is left
        # to gain has fallen below the objective's own roundoff. That is what
        # being at a mode looks like in floating point, and it is the
        # conclusion `_laplace_newton_unit_mode` reaches from the same evidence.
        accepted || break
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
        # Two, three and four parameter logistic in one expression. With
        # `P(y=1) = c + (d-c)F(η)` both responses are a constant plus a
        # multiple of the plain logistic term -- `c` and `d-c` for a one,
        # `1-d` and `d-c` for a zero -- so the same two lines cover all of
        # them, and the plain binary case is `c = 0, d = 1` where the constant
        # drops out. Added in log space because the constant is legitimately
        # zero there and the logistic term legitimately underflows.
        length(thresholds) >= 2 || return y > 0.5 ? -log1p_exp(-η) :
            -log1p_exp(η)
        c, d = _binary_asymptotes(thresholds, T)
        flat = y > 0.5 ? c : one(c) - d
        logF = y > 0.5 ? -log1p_exp(-η) : -log1p_exp(η)
        return _logaddexp(log(flat), log(d - c) + logF)
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
    # observation contributes its likelihood at that point and moves nothing.
    if !(s2 > T(_CTSEM_MIN_VARIANCE[]))
        return (_category_loglikelihood(ηbar, y, thresholds, kind),
            zero(T), zero(T))
    end
    # Three strategies, chosen by what the observation is rather than by how
    # hard it looks. Closed form where one exists, a mixture where the
    # likelihood is a combination of things that have one, and the rule only
    # for what is left.
    kind == CTSEM_OBS_CENSORED &&
        return _censored_moments(ηbar, s, y, thresholds)
    kind == CTSEM_OBS_BINARY && length(thresholds) >= 2 &&
        return _asymptote_moments(ηbar, s, y, nodes, weights, thresholds)
    return _binary_quadrature(ηbar, s, y, nodes, weights, thresholds, kind)
end

"""
    _asymptote_moments(ηbar, s, y, nodes, weights, thresholds)

`(logZ, mean - ηbar, variance)` for a binary item with asymptotes, as a
mixture rather than as one integral.

A three or four parameter logistic likelihood is not log-concave in `η`. For a
correct response `P = c + (1-c)F(η)` is a constant plus an increasing bounded
term, so `log P` is convex and then concave, and against a wide enough prior
the posterior has two modes -- measured at 10 items with `a = 1.7` and
`c = 0.2`, two maxima for any prior mean below about -4.3, separated by as
much as ten units and by as little as 0.01 of log posterior. Everything in the
quadrature path assumes one mode: the solve is claimed globally convergent
because the objective is concave, and the rule is centred on *the* mode. Handed
a bimodal posterior it would centre on whichever mode it found and integrate
the wrong bump, silently, which is the failure this file has already had once.

None of that has to be faced, because the likelihood is a mixture before it is
anything else. `P(y|η) = A + B q(y|η)`, with `q` the plain logistic term and

    A = c,      B = d - c      for a correct response
    A = 1 - d,  B = d - c      for an incorrect one

so the posterior is a two component mixture of the *prior* -- which is known in
closed form and needs no integration at all -- and the *plain binary
posterior*, which is log-concave and is exactly what the rule already handles
well. Integrating the components separately and combining their moments is
algebraically identical to integrating the lump, and every piece of it is
unimodal. Checked against direct integration of the three parameter posterior
on a fine grid, including at prior means where that posterior is bimodal:
agreement to 1e-11, which is the grid's own error.

The moments combine in offset coordinates, relative to `ηbar`, rather than as
absolute means. The prior component contributes an offset of exactly zero and a
variance of `s^2`; forming the same thing from absolute first and second
moments would subtract two numbers of order `ηbar^2`, which is the cancellation
`_binary_quadrature` accumulates about the mode to avoid.
"""
@inline function _asymptote_moments(ηbar::T, s::T, y::Real, nodes, weights,
    thresholds) where {T}
    s2 = s * s
    c, d = _binary_asymptotes(thresholds, T)
    flat = y > 0.5 ? c : one(c) - d
    span = d - c
    # An item with no span carries no information about `η`: every response is
    # the constant, so the posterior is the prior. `<=` rather than `==`
    # because a free asymptote pair can cross before the optimizer is pulled
    # back, and a negative span is not a likelihood.
    if !(span > zero(span))
        return (log(max(flat, zero(flat))), zero(T), s2)
    end
    logZq, offsetq, varq = _binary_quadrature(ηbar, s, y, nodes, weights, (),
        CTSEM_OBS_BINARY)
    logspan = log(span) + logZq
    logflat = log(flat)
    logZ = _logaddexp(logflat, logspan)
    isfinite(logZ) || return (T(-Inf), zero(T), s2)
    # The weight on the logistic component. Taken as a ratio of logarithms so
    # that a vanishing `flat` gives exactly one rather than `0/0`.
    w = exp(logspan - logZ)
    offset = w * offsetq
    second = (one(w) - w) * s2 + w * (varq + offsetq * offsetq)
    return (logZ, offset, max(second - offset * offset, zero(second)))
end

"""
    _binary_quadrature(ηbar, s, y, nodes, weights, thresholds, kind)

The adaptive Gauss-Hermite rule: what is used for the kinds whose likelihood
against a Gaussian prior has no elementary integral, and whose scalar posterior
is log-concave so that the mode the rule is centred on is unique.
"""
@inline function _binary_quadrature(ηbar::T, s::T, y::Real, nodes, weights,
    thresholds, kind::Int) where {T}
    s2 = s * s
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
    # A count row carries its dispersion in the same slot and for the same
    # reason: it is MANIFESTVAR's diagonal entry as a standard deviation, taken
    # before `sdcovsqrt2cov` assembles the matrix so that the reverse pass can
    # hand its cotangent straight back. A count is updated on its own,
    # sequentially, so it is never correlated with another row and the
    # off-diagonal it would otherwise need does not exist.
    if row <= length(types) && types[row] == CTSEM_OBS_COUNT
        isempty(ws.thresholds) && return view(ws.thresholds, 1:0)
        @inbounds ws.thresholds[1] = pars.MANIFESTVAR[row, row]
        return view(ws.thresholds, 1:1)
    end
    hasproperty(pars, :THRESHOLDS) || return view(ws.thresholds, 1:0)
    # A binary row with asymptotes reads two cells of the same matrix: the
    # lower asymptote and then a gap, accumulated here into the upper one so
    # that `0 <= c < d <= 1` holds without a cell's transform having to read
    # another cell's parameter. That is the constraint the parameter layer
    # imposes and the reason the ordinal thresholds below are gaps too.
    #
    # `nasymptotes` rather than the values decides whether the row has them: a
    # free guessing parameter sitting at zero is an ordinary place for an
    # optimizer to be, and inferring from the numbers would lose it there.
    if row <= length(types) && types[row] == CTSEM_OBS_BINARY
        na = row <= length(ws.nasymptotes) ? ws.nasymptotes[row] : 0
        (na >= 1 && size(pars.THRESHOLDS, 2) >= 2 &&
            length(ws.thresholds) >= 2) || return view(ws.thresholds, 1:0)
        @inbounds begin
            c = pars.THRESHOLDS[row, 1]
            ws.thresholds[1] = c
            ws.thresholds[2] = c + (one(c) - c) * pars.THRESHOLDS[row, 2]
        end
        return view(ws.thresholds, 1:2)
    end
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
    # A count's dispersion is additive and Gaussian on this same scalar, so it
    # joins the variance rather than the density -- see `_count_dispersion`.
    # Everything below is unchanged by it: `c` is still the covariance between
    # the state and the linear predictor and `s2` is still that predictor's
    # variance, which is all the projection below uses.
    if kind == CTSEM_OBS_COUNT
        σ = _count_dispersion(thresholds, T)
        s2 += σ * σ
    end
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
    # Nested: this runs inside the adjoint, which `ctsem_hessian` differentiates.
    J = _ctsem_nested_jacobian(triple, at)
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
    return _ctsem_nested_jacobian(triple, at)
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
