# New in ctsem 4.x

This page is for people who have used ctsem 3.x. It covers what version 4
changed and added, grouped by topic, with the calls involved. Later 4.x
releases add to it. [NEWS.md](NEWS.md) lists changes one line each.

Most of what is new depends on a second estimation engine, written in Julia,
which is now the default; unless a section says otherwise, what it describes
needs that engine. Models, data and most function calls carry over from 3.x
unchanged, and `ctFit(data, model, backend = 'stan')` fits as 3.x did.

## The Julia engine

```r
ctJuliaInstall()   # once
fit <- ctFit(data, model)
```

`ctJuliaInstall()` installs whatever is missing: Julia itself (a fixed version,
checked against a checksum), the JuliaConnectoR bridge package, and the
engine's Julia packages. It asks before each download. The first fit after
installing or updating ctsem builds a compiled image of the engine, which takes
several minutes; later sessions load it in a few seconds.
`ctJuliaStatus()` reports what was found.

What changes in practice:

- **Speed.** Fits that reach the same maximum take between a third and two
  fifths of Stan's time on one core, and the gap widens with more cores (see
  [Timings](#timings) below).
- **No compiling per model.** With Stan, a nonlinear model is compiled to C++
  before it can be fitted, which took about six minutes for the model in the
  table below, and again in each new R session. The Julia engine runs every
  model on the same compiled code, so a nonlinear model starts fitting as
  quickly as a linear one.
- **Threads.** The likelihood, gradient and Hessian are computed in parallel
  over subjects (or over groups, with nested levels) within one process.
  `cores` sets how many CPU cores a call may use, and defaults to
  `getOption('mc.cores', 2)`. The more work a fit has per subject, the more
  it gains: a panel with random auto-effects fitted 10.7 times faster on
  sixteen cores than on one.
- **Exact gradients.** The engine computes gradients exactly, with a
  hand-written reverse pass through the filter, and the Hessian behind the
  standard errors from those gradients.

### Timings

Seconds for a complete fit with each backend's defaults, standard errors
included, measured on 16 reserved cores of an AMD EPYC 9654 server on
2026-11-16 (ctsem 4.0.0 development version, R 4.6.1, Julia 1.13.3, rstan
2.32.7). Julia times are for a second fit in the same R session; the first fit
of a session also starts Julia, which added about 4 seconds. Stan's 16-core
times include starting its worker processes.

| model | Julia, 1 core | Julia, 4 | Julia, 16 | Stan, 1 core | Stan, 16 |
|---|---|---|---|---|---|
| The README's first example: 53 cows, 12,026 observations, 23 parameters | 18.0 | 9.6 | 7.1 | 56.3 | 27.5 |
| Simulated panel: 400 people, 40 observations each at irregular times, 23 parameters | 18.5 | 9.8 | 7.6 | 47.3 | 28.8 |
| The same design with random auto-effects as well, 34 parameters | 637 | 173 | 60 | | |

In each of the first two rows both backends reached the same log likelihood.
The third row was a first fit in its session and was not run on Stan, whose
filter would treat the random auto-effects differently (see Random effects
below).

A nonlinear version of the first example, with the effect of fat on yield
depending on fat, fitted to ten cows at 4 cores, shows the cost of compiling:
the Julia engine took 64 seconds as the first fit of a new R session. Stan took
402 and 441 seconds in two new sessions, and 76 seconds for a second fit in the
same session, so about six minutes of each first fit was compiling.

## Random effects

**Laplace approximation by default.** In 3.x every random effect was carried
as an extra state of the Kalman filter (`intoverpop = 'augmented'`). That is
exact when an effect shifts a mean and every indicator is Gaussian. For effects
on the dynamics it relies on a linearisation, which shrinks their estimated
variance. Version 4 integrates random effects out with a Laplace approximation
instead (`intoverpop = 'laplace'`, chosen automatically when optimising). In a
simulation with a random auto-effect, the augmented filter recovered about two
thirds of the true standard deviation and the Laplace approximation recovered
all of it. Effects that shift a mean in a linear Gaussian model (manifest
intercepts, initial means, continuous intercepts, time dependent predictor
effects) are still integrated exactly, inside each subject's filter.

**The approximation is checked.** After every Laplace fit, ctsem estimates the
approximation error for each subject and warns where correcting it would move
an estimate by more than about a tenth of its standard error. Two remedies are
available when it does:

- `ctFit(..., intoverpop = 'quadrature')`, or after the fit
  `ctFitUncertainty(fit, intoverpop = 'quadrature')`, refits with the flagged
  subjects integrated by adaptive quadrature.
- `ctFitUncertainty(fit, 'sample')`, or `ctFit(..., optimize = FALSE)`,
  samples the joint posterior of parameters and random effects, which involves
  no approximation of the integral.

**Nested levels.** `id = c('person', 'couple')` (innermost first) fits random
effects at several levels. A parameter's `indvarying` field names the levels it
varies at, for example `'cint||couple'` or `ctPar('cint', indvarying = 'couple')`.
See the [multilevel Laplace vignette](vignettes/laplace-multilevel.qmd).

**Restricted covariances.** `ctFit(poprank = )` gives the random effects
covariance a reduced rank, `popunique` adds a variance of each effect's own
(the factor analytic form), and `popblocks` arranges the shared dimensions in
blocks. These help when many parameters vary.

`intoverpop = 'augmented'` is still available on both backends.

## Measurement models

`manifesttype` in `ctModel()` takes one entry per observed variable, and types
can be mixed in one model:

| code | type | extra arguments |
|---|---|---|
| 0 | continuous (Gaussian) | |
| 1 | binary | `ASYMPTOTES` for guessing and ceiling parameters |
| 2 | ordinal | `ncategories` |
| 3 | count, Poisson-lognormal | the `MANIFESTVAR` diagonal is the dispersion |
| 4 | censored | `censormin`, `censormax`, and `RESPONSE` for rounding and scale |

Binary indicators existed in 3.x on Stan; the others are new, and all of them
work with random effects and nonlinear dynamics. Loadings, intercepts and
residual variances may depend on the latent state. For continuous indicators
such a measurement is now linearised at each observation's posterior mode
rather than at the prediction, which in tests against an exact grid filter
cut the error of the likelihood five- to thirtyfold. See the
[measurement models](vignettes/measurement-models.qmd),
[IRT items](vignettes/irt-items.qmd) and
[ordinal or continuous](vignettes/ordinal-or-continuous.qmd) vignettes.

## Nonlinear dynamics

Where the drift depends on the state more than quadratically (a cubic
restoring force, for example), the filter's prediction step matches moments by
cubature instead of linearising. `nlcontrol = list(nsubsteps = 'auto')`
chooses the number of integration steps for each observation interval from
the model's own linearisation error.

## Predictors and weights

- **Smooth effects** of time independent predictors:
  `ctModel(..., TIpredNames = ~ s(age) + sex)`. The smoothness is estimated
  from the data. `ctTIpredEffects(fit, whichTIpreds = 'age', plot = TRUE)`
  draws the estimated effect.
- **Missing predictor values** are imputed when optimising and sampled when
  sampling.
- **Predictors inside expressions.** A time independent predictor can be read
  by name in any model cell, so a known constant per subject can enter a
  matrix directly.
- **Fixed predictor effects** such as `'cint||||age=0.5'` are fitted. Stan in
  3.x ignored them, and now refuses them.
- **Sampling weights**, one per subject: `ctFit(..., weights = w)`, a weighted
  pseudo-likelihood with robust standard errors. See the
  [sampling weights vignette](vignettes/sampling-weights.qmd).

## Uncertainty and sampling

- An optimised fit reports standard errors from the Hessian, together with a
  check of whether the normal approximation they rest on holds. The
  check's verdict appears in a message and in `summary()`.
- `ctFitUncertainty(fit, method)` recomputes a fit's uncertainty: `'hessian'`,
  `'checked'` (the Hessian's normal, checked and importance sampled where it is
  wrong), `'is'`, `'sample'`, `'sandwich'`, `'opg'`, `'bootstrap'`,
  `'fullbootstrap'` or `'stored'`. It replaces `ctOptimUncertainty()`, which
  remains for old code.
- Sampling, with `ctFit(optimize = FALSE)` or `ctFitUncertainty(fit, 'sample')`,
  runs four chains at once and continues until effective sample size targets
  are reached. It samples the joint posterior of parameters and random
  effects by default. See the [uncertainty vignette](vignettes/uncertainty.qmd).

## Checking a model, and building on it

| function | what it does |
|---|---|
| `ctReport(fit)` | runs the reporting and checking functions and writes their output to a folder, with an index |
| `ctIdentify(data, model)` | which parameters the data can inform, before fitting |
| `ctFitProfile(fit)` | profile likelihood for chosen parameters |
| `ctScoreTest(fit)` | tests, without refitting, whether the model is missing random effects, covariate effects or state dependence ([vignette](vignettes/score-tests.qmd)) |
| `ctFitAdequacy(fit)` | tests of the one-step-ahead prediction errors, by family |
| `ctModelImprove(fit)`, `ctModelProposals(fit)` | step-by-step model search from those tests |
| `ctModelHigherOrder()`, `ctModelSlowComponent()`, `ctModelPeriodicComponent()`, `ctModelUniqueProcess()`, `ctModelLatentLevel()` | add a second order, a slow process, a cycle, an item-specific process, or move levels into the latent processes |
| `ctLOO(fit, method = 'psis')` | cross validation of a Laplace fit without refitting |
| `ctLaplaceCheck(fit)`, `ctParticleLik(fit)` | measure the error of the Laplace approximation and of the Kalman filter |
| `ctTracePlot(fit)` | the optimisation path or the sampling chains |
| `anova()`, `coef()`, `vcov()`, `logLik()` | the usual generics, for fits from either backend |

## Views of the dynamics

`ctNetwork()` and `ctNetworkPlot()` draw the temporal and contemporaneous
networks of a fitted system, `ctPhasePortrait()` its phase portrait,
`ctStateDependencePlot()` how model matrices vary with the state in a
nonlinear model, and `ctVarianceDecomposition()` splits each indicator's
variance into between person, process and measurement parts.

## Specifying models

- Cells take named fields as well as the positional `|` form:
  `'cint, indvarying=TRUE, sdscale=0.5, tipreds=c(age)'`. `ctPar()` writes
  either form, and `ctCov()` writes a covariance matrix as cells.
- Covariance matrices are parameterised by a matrix logarithm, so every
  positive definite matrix can be reached. An off-diagonal parameter is not a
  correlation (for two variables it is the correlation's Fisher z); the
  correlations, with intervals, are `DIFFUSIONcor`, `T0cor` and `MANIFESTcor`
  in `summary(fit)$parmatrices` and `ctSummaryMatrices()`.
- `ctModel(driftScaled = c('DIFFUSION', 'CINT'))` expresses the diffusion and
  intercept as what accumulates over one time unit, which makes uncertainty
  and sampling better behaved when the drift is weakly determined.
- `ctGenerate()` generates from any model the engine can fit, random effects
  at every level and non-Gaussian indicators included. `values =` sets
  population values on the scales `summary()` reports, `datastruct =` fills in
  a given design, and `fromPriors = TRUE` draws from the prior predictive.

## Changed defaults

| | 3.11 | 4.0 | for the 3.11 behaviour |
|---|---|---|---|
| backend | Stan | Julia | `backend = 'stan'` |
| random effects when optimising, `intoverpop = 'auto'` or `TRUE` | augmented filter | Laplace (Julia) | `intoverpop = 'augmented'` |
| covariance parameterisation | `'rawcorr'` | matrix logarithm (`'z'`) | `model$covmattransform <- 'rawcorr'` |
| `cores` for `ctFitUncertainty()`, `ctFitCheckCov()`, `ctKalmanArray()` | 1 | `getOption('mc.cores', 2)` | `cores = 1` |
| `cores` when sampling | one per chain | `getOption('mc.cores', 2)` | `cores = 'maxneeded'` |
| ACF plots | quantile regression | precision weighted spline | `method = 'quantile'` |
| subjects in `ctDiscretePars()`, `ctPredictTIP()`, `ctExtract()`, `ctKalmanArray()` | position | id in the data | `realid = FALSE` |
| importance sampling target | 100 effective draws | 200 | `optimcontrol$isESS = 100` |

A fit saved by 3.x and passed to `ctFitUpdate()` stays on Stan.

## Renamed and removed

| old | now |
|---|---|
| `ctOptimUncertainty()` | `ctFitUncertainty()` (old name kept) |
| `ctFitCovCheck()` | `ctFitCheckCov()` (old name kept) |
| `ctGenerateFromPriors()` | `ctGenerate(fromPriors = TRUE)` (deprecated) |
| `ctModelCoverage_check()` | `ctCoverageCheck()` |
| `ctCheckFit()` | `ctPostPredPlots()` and `ctFitCheckCov()` |
| `ctFitUpdateModel()`, `ctStanUpdModel()` | `ctFitUpdate()` |
| `ctGenerate(n.subjects = )` | `ctGenerate(n = )` |
| `ctDiscretePars(observational = )` | `impulseType = 'observed'` or `'unit'` |
| `ctFit(iter = , chains = , control = )` | entries of `sampleControl` |
| `uncertainty = 'surrogate'` | `'checked'` |
| `ctFit(vb = TRUE)` | removed |

## Later 4.x releases

Nothing yet.
