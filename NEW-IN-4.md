# New in ctsem 4

Version 4 adds an estimation engine written in Julia, now the default,
and most new features build on it. Most code written for version 3 still
runs. The tables at the end list the changed defaults and the renamed
and removed functions.
[NEWS.md](https://cdriveraus.github.io/ctsem/NEWS.md) is the change log.

## Julia engine

``` r

ctJuliaInstall()   # once
fit <- ctFit(data, model)
```

[`ctJuliaInstall()`](https://cdriveraus.github.io/ctsem/reference/ctJuliaInstall.md)
sets up Julia and the engine, asking before each download. The first fit
after installing or updating ctsem takes a few minutes longer while the
engine compiles.

- Fits are two to three times faster than with Stan, and use multiple
  cores (`cores`).
- New models need no compilation, so a nonlinear model starts fitting at
  once. With Stan, each new nonlinear model took minutes to compile.
- A saved fit can be summarised without Julia:
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`ctSummaryMatrices()`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md),
  [`ctDiscretePars()`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md),
  [`ctExtract()`](https://cdriveraus.github.io/ctsem/reference/ctExtract.md),
  [`ctModelLatex()`](https://cdriveraus.github.io/ctsem/reference/ctModelLatex.md)
  and
  [`ctNetwork()`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md)
  work on it. `ctsemtestfit` is an example fit.
- `backend = 'stan'` fits as ctsem 3 did.

## Random effects

- Random effects are integrated with a Laplace approximation. Version 3
  carried them as extra states of the Kalman filter, which
  underestimates individual differences in the dynamics. In linear
  models, random intercepts are integrated exactly.
  `intoverpop = 'augmented'` gives the old method, and
  [`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md)
  warns when that method meets nonlinear dynamics.
- Each fit checks the approximation and says if it affects the
  estimates, with a warning when they would move by a standard error or
  more. `intoverpop = 'quadrature'` then refines the fit, and
  `ctFitUncertainty(fit, 'sample')` samples the full posterior instead.
- `intoverpop = 'variational'` (experimental) integrates the random
  effects with an approximation to each subject’s posterior that is
  fitted during the optimisation. On some models this gives a higher log
  likelihood than the Laplace approximation, and fits take about twice
  as long.
- Random effects can be at multiple nested levels, such as people within
  couples: `id = c('person', 'couple')`, with each parameter’s
  `indvarying` naming the levels it varies at. See the [multilevel
  vignette](https://cdriveraus.github.io/ctsem/articles/laplace-multilevel.html).
- `poprank`, `popunique` and `popblocks` give the random effects
  covariance fewer parameters when many parameters vary.

## Measurement models

Observed variables can be continuous, binary, ordinal, counts, censored
or a choice with its response time, in any mix, set for each variable
with `manifesttype`.

| `manifesttype` | type | also needed |
|----|----|----|
| 0 | continuous |  |
| 1 | binary | optional `ASYMPTOTES` for guessing and ceiling |
| 2 | ordinal | `ncategories` |
| 3 | count |  |
| 4 | censored | `censormin`, `censormax` |
| 5 | choice and response time (drift diffusion model) | optional `WIENER` for boundary, nondecision time and bias |

Nonlinear measurement and dynamics are filtered more accurately than in
version 3. See the vignettes on [measurement
models](https://cdriveraus.github.io/ctsem/articles/measurement-models.html),
[item response
models](https://cdriveraus.github.io/ctsem/articles/irt-items.html) and
[ordinal or continuous
data](https://cdriveraus.github.io/ctsem/articles/ordinal-or-continuous.html).

[`ctModelDDM()`](https://cdriveraus.github.io/ctsem/reference/ctModelDDM.md)
builds a drift diffusion model of choices and response times with the
drift as a latent process, constant or changing over the session.

## Covariates and weights

- Smooth effects of time independent covariates:
  `TIpredNames = ~ s(age) + sex`.
- Subject weights, for population estimates from a sample that does not
  represent the population: `ctFit(..., weights = w)`. See the [sampling
  weights
  vignette](https://cdriveraus.github.io/ctsem/articles/sampling-weights.html).

## Uncertainty and sampling

- Standard errors come from the Hessian, as before, and each fit checks
  whether that approximation holds.
- [`ctFitUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md),
  which replaces
  [`ctOptimUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md),
  computes the uncertainty by other methods, such as a sandwich
  estimator, the bootstrap or importance sampling.
- `ctFit(optimize = FALSE)`, or `ctFitUncertainty(fit, 'sample')` after
  a fit, samples the posterior until the effective sample size reaches
  its target. See the [uncertainty
  vignette](https://cdriveraus.github.io/ctsem/articles/uncertainty.html).
- When sampling, ctsem checks whether the first draws are centred away
  from the estimate, which can mean the optimiser stopped at a local
  maximum. If optimising from there reaches a higher log posterior,
  `ctFit(optimize = FALSE)` fits again from that point and
  `ctFitUncertainty(fit, 'sample')` gives a warning.
- If the approximation for the random effects would move an estimate by
  a standard error or more, sampling refines the fit first and starts
  there.

## Checking and improving models

| function | use |
|----|----|
| `ctReport(fit)` | run the main summaries, plots and checks, and save them to a folder |
| `ctIdentify(data, model)` | see which parameters the data can inform, before fitting |
| `ctScoreTest(fit)` | test whether the model is missing random effects, covariate effects or nonlinearity, without refitting ([vignette](https://cdriveraus.github.io/ctsem/articles/score-tests.html)) |
| `ctFitAdequacy(fit)` | test the model’s one-step-ahead predictions |
| `ctModelImprove(fit)` | search for a better model step by step, guided by score tests |
| `ctModelImprove(fit, reduceOnly = TRUE)` | simplify the parts of a model the data do not identify: parameters the data do not determine, or a random effects covariance of too high a rank |
| `ctLOO(fit, method = 'psis')` | cross validation without refitting |
| `ctFitProfile(fit)` | profile likelihood |
| `ctFitFlatRetry(fit)` | fit again from the other end of the range of a parameter the fit left where the data say little about it, and keep the better fit |
| [`anova()`](https://rdrr.io/r/stats/anova.html), [`coef()`](https://rdrr.io/r/stats/coef.html), [`vcov()`](https://rdrr.io/r/stats/vcov.html), [`logLik()`](https://rdrr.io/r/stats/logLik.html) | the standard R generics |

[`ctModelHigherOrder()`](https://cdriveraus.github.io/ctsem/reference/ctModelHigherOrder.md),
[`ctModelSlowComponent()`](https://cdriveraus.github.io/ctsem/reference/ctModelSlowComponent.md),
[`ctModelPeriodicComponent()`](https://cdriveraus.github.io/ctsem/reference/ctModelPeriodicComponent.md)
and
[`ctModelUniqueProcess()`](https://cdriveraus.github.io/ctsem/reference/ctModelUniqueProcess.md)
extend a model with second order dynamics, a slow process, a cycle, or a
process specific to one indicator.

## Plots

[`ctNetwork()`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md)
and
[`ctNetworkPlot()`](https://cdriveraus.github.io/ctsem/reference/ctNetworkPlot.md)
draw temporal and contemporaneous networks,
[`ctPhasePortrait()`](https://cdriveraus.github.io/ctsem/reference/ctPhasePortrait.md)
a phase portrait, and
[`ctStateDependencePlot()`](https://cdriveraus.github.io/ctsem/reference/ctStateDependencePlot.md)
how a nonlinear model’s parameters change with the state.
[`ctVarianceDecomposition()`](https://cdriveraus.github.io/ctsem/reference/ctVarianceDecomposition.md)
splits variance into between subject, process and measurement parts, and
[`ctTracePlot()`](https://cdriveraus.github.io/ctsem/reference/ctTracePlot.md)
shows the path of the optimiser or sampler.

## Model specification

- A model can be written in lavaan, Mplus or brms style syntax with
  [`ctModelLavaan()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md),
  [`ctModelMplus()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md)
  and
  [`ctModelFormula()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md),
  and
  [`ctModelSyntax()`](https://cdriveraus.github.io/ctsem/reference/ctModelSyntax.md)
  translates between them
  ([vignette](https://cdriveraus.github.io/ctsem/articles/model-syntax.html)).
- Parameter options can be named, as in
  `'cint, indvarying=TRUE, tipreds=c(age)'`.
  [`ctPar()`](https://cdriveraus.github.io/ctsem/reference/ctPar.md)
  writes these, and
  [`ctCov()`](https://cdriveraus.github.io/ctsem/reference/ctCov.md)
  writes covariance matrices.
- Covariance matrices use a new parameterisation (a matrix logarithm),
  so their off-diagonal parameters are not correlations.
  [`summary()`](https://rdrr.io/r/base/summary.html) reports the
  correlation matrices (`DIFFUSIONcor`, `T0cor`, `MANIFESTcor`). An
  off-diagonal parameter fixed to zero gives a covariance of zero, but
  an estimate near zero does not show a covariance near zero with three
  or more variables. Read covariances and correlations from the reported
  matrices, not from the parameters.
- [`ctGenerate()`](https://cdriveraus.github.io/ctsem/reference/ctGenerate.md)
  simulates data from any model, including random effects and
  non-Gaussian indicators.

## Changed defaults

|  | ctsem 3.11 | ctsem 4 | for the old behaviour |
|----|----|----|----|
| backend | Stan | Julia | `backend = 'stan'` |
| random effects when optimising | Kalman filter states | Laplace approximation | `intoverpop = 'augmented'` |
| covariance parameterisation | `'rawcorr'` | matrix logarithm | `model$covmattransform <- 'rawcorr'` |
| `cores` in [`ctFitUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md), [`ctFitCheckCov()`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md), [`ctKalmanArray()`](https://cdriveraus.github.io/ctsem/reference/ctKalmanArray.md) | 1 | `getOption('mc.cores', 2)` | `cores = 1` |
| `cores` when sampling | one per chain | `getOption('mc.cores', 2)` | `cores = 'maxneeded'` |
| ACF plots | quantile regression | weighted spline | `method = 'quantile'` |
| subjects in [`ctDiscretePars()`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md), [`ctPredictTIP()`](https://cdriveraus.github.io/ctsem/reference/ctPredictTIP.md), [`ctExtract()`](https://cdriveraus.github.io/ctsem/reference/ctExtract.md), [`ctKalmanArray()`](https://cdriveraus.github.io/ctsem/reference/ctKalmanArray.md) | by position | by id | `realid = FALSE` |

A fit saved with ctsem 3 and passed to
[`ctFitUpdate()`](https://cdriveraus.github.io/ctsem/reference/ctFitUpdate.md)
is refitted with Stan.

## Renamed and removed

| old | new |
|----|----|
| [`ctOptimUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md) | [`ctFitUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md) (old name still works) |
| [`ctFitCovCheck()`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md) | [`ctFitCheckCov()`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md) (old name still works) |
| `ctstantestdat` | `ctsemtestdat` (old name still works) |
| [`ctGenerateFromPriors()`](https://cdriveraus.github.io/ctsem/reference/ctGenerateFromPriors.md) | `ctGenerate(fromPriors = TRUE)` |
| [`ctModelCoverage_check()`](https://cdriveraus.github.io/ctsem/reference/ctsem-defunct.md) | [`ctCoverageCheck()`](https://cdriveraus.github.io/ctsem/reference/ctCoverageCheck.md) |
| [`ctCheckFit()`](https://cdriveraus.github.io/ctsem/reference/ctsem-defunct.md) | [`ctPostPredPlots()`](https://cdriveraus.github.io/ctsem/reference/ctPostPredPlots.md) and [`ctFitCheckCov()`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md) |
| [`ctFitUpdateModel()`](https://cdriveraus.github.io/ctsem/reference/ctsem-defunct.md), [`ctStanUpdModel()`](https://cdriveraus.github.io/ctsem/reference/ctsem-defunct.md) | [`ctFitUpdate()`](https://cdriveraus.github.io/ctsem/reference/ctFitUpdate.md) |
| `ctGenerate(n.subjects = )` | `ctGenerate(n = )` |
| `ctDiscretePars(observational = )` | `impulseType = 'observed'` or `'unit'` |
| `ctFit(iter = , chains = , control = )` | `sampleControl` |
| `uncertainty = 'surrogate'` | `'checked'` |
| `ctFit(vb = TRUE)` | removed |
| `ctEmpiricalBayesFit()` | [`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) with random effects, and [`ctSubjectPars()`](https://cdriveraus.github.io/ctsem/reference/ctSubjectPars.md) for each subject’s estimates |

## Later 4.x releases

Nothing yet.
