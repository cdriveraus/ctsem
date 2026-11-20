# Score tests for relaxing a fitted model

Tests, from one fitted model and without refitting it, whether the model
would fit better if it were relaxed: parameters that differ between
subjects (a missing random effect), parameters that depend on a time
independent covariate the model does not include, or dynamics that
depend on the latent state (a missing nonlinearity). Each test is a
score (Lagrange multiplier) test: the slope of the log likelihood in the
direction of the relaxation, at the fitted model, with the fitted
parameters projected out. Each also gives a one-step estimate of the
size of the relaxation. These estimates are for screening, so refit with
the relaxations that matter before interpreting their sizes.

## Usage

``` r
ctScoreTest(
  fit,
  covariates = NULL,
  tests = c("heterogeneity", "covariates", "statedependence", "volatility", "fixed",
    "level", "instability"),
  form = c("auto", "linear", "spline", "factor", "threshold"),
  df = 3,
  parameters = NULL,
  order = 2L,
  variance = c("robust", "model"),
  nboot = 1999L,
  exact = 5L,
  intoverpop = c("laplace", "quadrature"),
  cores = getOption("mc.cores", 2L),
  verbose = TRUE
)
```

## Arguments

- fit:

  An optimised fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) with
  `backend = 'julia'`.

- covariates:

  Optional data frame of subject level covariates to test, with a column
  named as the fit's subject id (or `id`) and one column per covariate.
  Long data with repeated rows per subject can be used if each covariate
  is constant within subject. Without an id column, give one row per
  subject in the fit's subject order.

- tests:

  Which families to run: any of `'heterogeneity'`, `'covariates'`,
  `'statedependence'`, `'volatility'`, `'fixed'`, `'level'` and
  `'instability'`. By default all but `'instability'`. Covariates are
  tested only when `covariates` is given or `'covariates'` is asked for
  (then the model's own predictors alone), and levels only where the
  model has a random intercept.

- form:

  How a covariate enters: `'auto'` (indicators for a factor, character
  or logical covariate, linear otherwise), `'linear'`, `'spline'`,
  `'factor'` or `'threshold'` (a step at an unknown cut). One value, or
  one per covariate.

- df:

  Degrees of freedom of the spline basis.

- parameters:

  Optional parameter names to restrict the heterogeneity and covariate
  tests to.

- order:

  Highest degree of the state dependence terms. 2 tests each drift
  coefficient varying linearly with each state. 3 adds variation with
  products and squares of states, which is needed to see a symmetric
  dependence, such as regulation that strengthens at both extremes.

- variance:

  `'robust'` (the default) or `'model'`; see the section on variance.

- nboot:

  Bootstrap resamples for the robust p-values.

- exact:

  For `statedependence`, `volatility` and `fixed` with the robust
  variance: how many of each family's terms are corrected exactly for
  the estimation of the fitted parameters, which takes two evaluations
  of the extended model per term. Every term is first tested without
  that correction, which is conservative, and the strongest `exact`
  terms are then corrected (`projected` in the tables). `Inf` corrects
  every term.

- intoverpop:

  How the extended families (`statedependence`, `volatility`, `fixed`,
  `level`) integrate the random effects on a Laplace fit. `'laplace'`
  (the default) uses the Laplace approximation, a fast screen.
  `'quadrature'` integrates every subject's random effects by
  quadrature, as `ctFit(intoverpop = 'quadrature')` does, which removes
  the Laplace error from these families and takes longer. The other
  families use the fit's own objective either way.

- cores:

  Engine cores for the evaluations.

- verbose:

  Report progress.

## Value

An object of class `ctScoreTest`, a list with `heterogeneity` (a data
frame, one row per parameter or fixed `MANIFESTMEANS` cell), `dimension`
(the largest eigenvalue, its p-value and the direction, when two or more
parameters are tested), `blocks` (on a reduced rank level, `dimension`,
a row per block with the one-step `sd` along its direction,
`covariance`, a row per pair of blocks, and `directions`, the direction
found for each tested block), `covariates` (one element per covariate,
each a list of its `joint` test and a data frame `parameters`, including
the model's own predictors, marked by `joint$inmodel`),
`statedependence` (a list of `terms`, a data frame, and `joint`, the
tests per equation and overall), `volatility` (the same, for the
innovation sds), `fixed` (a list of `cells` and their `joint` test),
`instability` (a list of `parameters`, one row per ordering and
parameter, and `joint`, per ordering), and `info` (how the fit
integrated its random effects, the number of subjects and parameters,
the variance used, and `errors`, why any family was not tested). Tables
have a `class` column where relaxations are alike, and the extended
families a `projected` column, whether a term was corrected exactly (see
`exact`).

## Heterogeneity

Tests letting each parameter vary over subjects, as `ctFit` would fit
it, for each free parameter that does not vary already and is not an
initial condition (`T0MEANS` and `T0VAR`). The test covers the
parameter's variance over subjects together with its covariances with
what already varies: each subject's initial state (the free `T0MEANS`
parameters) and the parameters with random effects, listed in `with`.
Even a model without random effects has initial states, so these
covariances are tested there too. They matter because a varying `CINT`
intercept that is unrelated to the initial means makes subjects drift
apart from a common start, while the same intercept correlated with each
subject's start is a random intercept.

The covariances are tested as ordinary score tests, and the variance,
given them, as a one sided test (Stram & Lee, 1994), so `LM` is referred
to an equal mixture of chi-square distributions on `df - 1` and `df`
degrees of freedom. `z` is the one sided test of the variance alone,
with its covariances held at zero, which has more power when the new
effect is independent of what already varies. `p` is the smaller of the
two tests' p-values, calibrated by resampling subjects. A negative `z`
beside a large `LM` means the parameter varies only together with what
already varies, as an intercept does with each subject's start. No
assumption is made about how the initial states relate to the dynamics.

The statistic is the information matrix test, which is the score test
for neglected heterogeneity (Chesher, 1984). Other misspecification can
also move it. A drift that varies over subjects also moves the
diffusion's test, for example, and a missing nonlinearity moves several.
`sd` is the one-step estimate of the subjects' standard deviation on the
raw scale, with the covariances estimated alongside, and `sd_natural`
approximately the same on the parameter's own scale. They are `NA` where
the test comes from the covariances and the variance alone is estimated
as zero, as often happens for an intercept, and then only a refit can
give the size. To act on a result, set `indvarying = TRUE` for the
parameter and the `T0MEANS` parameters, refit, and test again.
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
does this.

Each `MANIFESTMEANS` cell the model fixes is also tested, as a random
intercept on that indicator around its fixed value (rows named
`MANIFESTMEANS[...]`). With one indicator per latent this overlaps with
the covariances of the `CINT` rows, with more indicators it is the
indicator's own intercept, and for a latent with no free intercept it is
the only test of its level. To act on it, give the cell a parameter with
the fixed value as its stated mean (`ctPar(mean=)`) and
`indvarying = TRUE`.

On a reduced rank level (`poprank`), a parameter that varies only
through the shared dimensions, without a variance of its own, is tested
for variation outside those dimensions. Each such test is of the
variance alone. A parameter whose test is already determined by the
fitted model, as for an intercept that loads on the shared dimensions,
is not tested (`NA`) and is left out of `dimension`.

`dimension` tests whether subjects vary along some combination of the
tested parameters, a new dimension of variation, with its covariances
with what already varies held at zero. The statistic is the largest
eigenvalue of the standardised score matrix of the variances, tested
against a bootstrap over subjects, and its eigenvector gives the
direction, as relative loadings on the raw scale. It can detect
variation spread over several parameters that no single test picks up,
and is less sensitive than the best single test when one parameter
carries it. It is not computed over a reduced rank level's parameters,
whose blocks are tested instead.

On a reduced rank level, `blocks` tests what its blocks of shared
dimensions (`popblocks`, or without them the whole level as one block)
leave out, beside the own variances. `blocks$dimension` tests, within
each block with more effects than dimensions, whether the block needs
another shared dimension, with the block's fitted loadings projected
out. The rows for single parameters, which test own variances, cannot
show this. `blocks$covariance` tests, for each pair of blocks, the
covariances between their dimensions, which the model holds at zero,
with a joint two sided score test. An example is whether a response
style shared by an item's thresholds covaries with the dynamic level. A
pair in which one block loads on every effect of the other already has
these covariances in the model, and is listed with a note but not
tested. Both are calibrated by resampling subjects and adjusted together
(`p_adj`). An effect whose population mean is stated is tested through
the model with that mean free and set to its stated value, and only the
fit's own parameters are projected out.
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
acts on both, adding a dimension to a block or letting blocks share
dimensions.

## Covariates

For each covariate, tests the effects it would have on every tested
parameter if it were added as a time independent predictor (on the raw
scale), including effects on the population standard deviations of the
parameters that vary. The `joint` test of all of a covariate's effects
is the screen. For each parameter, `estimate` is the one-step effect per
unit of the covariate, estimated jointly over the parameters with the
other fitted parameters held at their estimates, so that a covariate
that shifts one parameter is not attributed to another that is
correlated with it. Holding the other parameters keeps the estimate
close to what a refit would give. `estimate_natural` gives it
approximately on the parameter's own scale, and `p` tests it. A
covariate enters linearly (centred), as treatment indicators for a
factor, or as a natural spline with `df` degrees of freedom, for which
only `p` is given. Subjects with a missing value contribute at the mean.
To act on a result, add the covariate to the data and to `TIpredNames`,
and refit.

The model's own time independent predictors are tested as well, each on
the parameters it does not yet act on, using the values the fit read.
Their results have `joint$inmodel` set to `TRUE`. A predictor added with
an effect on one parameter may act on another too, and that can show
only once other changes are made. To act on a result, set the
predictor's effect on that parameter (`<predictor>_effect` in
`model$pars`) and refit.
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
does that.

`form = 'threshold'` tests for a step at an unknown cut instead. The
joint test is computed at each cut from the 15th to the 85th percentile,
and the largest (Andrews, 1993), which is the split a score based tree
makes, is compared with the same scan with the covariate permuted over
subjects (Hothorn & Lausen, 2003). This accounts for the search and
needs no refitting. The effects reported are those of being above the
best cut, given as `cut` in the joint row.

On a reduced rank level, a covariate's effect on the spread of a varying
parameter is tested as a scaling of that parameter's loadings and own sd
together, reported as `scale(<parameter>)`: the effect on its log
standard deviation per unit of the covariate, with the correlations held
fixed. With about 100 subjects these rows reject a correct model
somewhat too often (6.5 to 8 per cent at .05), so read them only once
the covariate's joint test has picked it out.

## State dependence

For each latent equation and each pair of latent states `k <= m`, tests
a term `d * eta_k * eta_m` in that equation's drift, that is, a drift
coefficient that varies linearly with a state (Luukkonen, Saikkonen &
Terasvirta, 1988). A term `eta_k * eta_m` in the equation of `eta_i` is
`DRIFT[i,k]` varying with `eta_m`, or when `k != m` equally `DRIFT[i,m]`
varying with `eta_k`. The data cannot tell these apart, and `reads`
gives both. `estimate` is the one-step change in that drift coefficient
per unit of the state. It is `NA` when the fitted model cannot size the
term, which happens when the departure is strong, and the test itself is
still valid. `order = 3` adds the terms of degree three. The terms use
the state predicted from earlier observations, as the filter does. To
act on a result, write the dependence into the drift cell, for example
`DRIFT[1,1] = 'dr11 + d * eta2'` with
`PARS = c('dr11|-log1p_exp(-param)', 'd')`, and refit.

## Volatility

`tests = 'volatility'` tests whether the size of each latent's
innovations depends on a state, by scaling `DIFFUSION[i, i]`, the
innovation sd, by `exp(v * eta_m)` and testing `v = 0`, for every latent
with innovations and every state. Variance that rises with a process's
own level, as when negative affect is more variable when high, is the
diagonal term. As for state dependence, the predicted state is used. To
act on a result, write the scaling into the cell, for example
`DIFFUSION[2,2] = 'sd2 * exp(v * eta2)'` with
`PARS = c('sd2|log1p_exp(param)', 'v')`, and refit.
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
does this. These terms make the model nonlinear, so the test and a refit
use the extended Kalman filter.

## Level

`tests = 'level'` tests whether a parameter of the dynamics or of the
error variances depends on a subject's level, where the level is a
random intercept. For each free `DRIFT`, `DIFFUSION` and `MANIFESTVAR`
parameter that does not vary over subjects, including a `PARS` parameter
that those matrices read in an expression, and each `MANIFESTMEANS` or
`CINT` intercept that does, it tests a term `g * intercept` added to the
parameter, at `g = 0`, wherever the model reads the parameter. When
subjects' levels are in the intercepts, the latent states are deviations
from them, so the state dependence and volatility tests only see
dependence on the deviation. A parameter that depends on a subject's
usual level, such as persistence or innovation size that differs between
people who are usually high and those usually low, otherwise shows only
as heterogeneity, with no direction. These terms test that direction.
The score is taken from the extended model's gradient, averaged over
each subject's posterior for the intercept. `estimate` is the one-step
change in the parameter, on its own scale, per unit of the intercept. To
act on a result, let the parameter vary over subjects, so its
correlation with the intercept is estimated.
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
does this. A model with no random intercept has nothing to test.

## Fixed cells

For each cell the model fixes at zero that could sensibly be freed (a
`DRIFT` cross effect, a `DIFFUSION` or `MANIFESTVAR` covariance between
two that have a variance, a `LAMBDA` cross loading, a `TDPREDEFFECT`),
tests freeing it. This is the modification index (Sorbom, 1989), and
`estimate` is the expected parameter change. `CINT` and `MANIFESTMEANS`
are not offered, since with the other one free, either is identified
only through the initial transient. State dependence, volatility and
fixed cells are tested on one extended model, so asking for several
costs one preparation.

## Stability

`tests = 'instability'`, run only when asked for, tests whether each
parameter is stable across the data, with the data ordered by time, and
by each latent's state as predicted from earlier observations. It uses
each data row's score, whose partial sums behave as a Brownian bridge
under the model (Nyblom, 1989). A parameter that changes over the study,
or differs at high and low levels of a state, bends them. Ordered by a
state, this is a test of state dependence of any form. The p-values come
from a wild bootstrap over rows. It can take a long time on large data.

## Alike relaxations

Two relaxations whose contributions from each subject are nearly
collinear (correlation of 0.9 or more) would explain the same signal,
and the data cannot tell them apart (MacCallum et al., 1993). They share
a letter in `class`, and `NA` means alike to none. Refit with the most
plausible one, then test again.

## Variance and p-values

By default (`variance = 'robust'`) the variance of each test is the
spread over subjects of their projected contributions, and its p-value
comes from a bootstrap over subjects of the studentised statistic
(`nboot` resamples, from a fixed random stream, so a call is
reproducible and leaves the session's random numbers alone). This keeps
each test valid when the model is wrong in other ways.
`variance = 'model'` uses the fitted model's curvature for the variance
and chi-square and normal reference distributions. It is always used
with fewer than 20 subjects.

## Reading the results

The tests are screens. A badly misspecified model shows many significant
tests at once, and the one-step sizes are local approximations,
especially for drift parameters. Relax what stands out most, refit, and
test again. Within each family `p_adj` is the Benjamini-Hochberg
adjustment over the tests shown.

The tests use the fit's own objective, so they are for the `intoverpop`
it was fitted with. They need an optimised fit, and refuse a sampled
one. Heterogeneity and covariate tests need subjects as the independent
units, so a model with a level above the subject can only be tested for
state dependence. If a family cannot be tested on a model, the reason is
given as a warning and in `info$errors`, and the other families' results
stand. With `priors = TRUE`, the new terms are tested on the likelihood,
without a prior.

On a Laplace fit with random dynamics beside random intercepts, the
state dependence and volatility families use the one-step Laplace
approximation for both the fit and its extension (see
`optimcontrol$laplace_onestep` in
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)). This
is a fast screen whose approximation error stays in the scores, and
`info$onestep` names the families concerned. `intoverpop = 'quadrature'`
integrates every subject's random effects for these families instead
(`info$quadrature`).

## References

Boos, D. D. (1992). On generalized score tests. The American
Statistician, 46, 327-333.

Chesher, A. (1984). Testing for neglected heterogeneity. Econometrica,
52, 865-872.

Stram, D. O., & Lee, J. W. (1994). Variance components testing in the
longitudinal mixed effects model. Biometrics, 50, 1171-1177.

Self, S. G., & Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. Journal of the American Statistical Association, 82,
605-610.

Hothorn, T., & Lausen, B. (2003). On the exact distribution of maximally
selected rank statistics. Computational Statistics & Data Analysis, 43,
121-137.

Andrews, D. W. K. (1993). Tests for parameter instability and structural
change with unknown change point. Econometrica, 61, 821-856.

Nyblom, J. (1989). Testing for the constancy of parameters over time.
Journal of the American Statistical Association, 84, 223-230.

Hall, D. B., & Praestgaard, J. T. (2001). Order-restricted score tests
for homogeneity in generalised linear and nonlinear mixed models.
Biometrika, 88, 739-751.

MacCallum, R. C., Wegener, D. T., Uchino, B. N., & Fabrigar, L. R.
(1993). The problem of equivalent models in applications of covariance
structure analysis. Psychological Bulletin, 114, 185-199.

Luukkonen, R., Saikkonen, P., & Terasvirta, T. (1988). Testing linearity
against smooth transition autoregressive models. Biometrika, 75,
491-499.

Sorbom, D. (1989). Model modification. Psychometrika, 54, 371-384.

## See also

[`vignette('score-tests')`](https://cdriveraus.github.io/ctsem/articles/score-tests.md)
for three worked examples.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- ctFit(dat, model, backend = 'julia')
st <- ctScoreTest(fit, covariates = subjectdata)
st
st$heterogeneity
} # }
```
