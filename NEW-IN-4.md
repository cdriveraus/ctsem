# New in ctsem 4

Version 4 adds an estimation engine written in Julia, now the default, and most
new features build on it. Most code written for version 3 still runs. The
tables at the end list the changed defaults and the renamed and removed
functions. [NEWS.md](NEWS.md) is the change log.

## Julia engine

```r
ctJuliaInstall()   # once
fit <- ctFit(data, model)
```

`ctJuliaInstall()` sets up Julia and the engine, asking before each download.
The first fit after installing or updating ctsem takes a few minutes longer
while the engine compiles.

- Fits are two to three times faster than with Stan, and use multiple cores
  (`cores`).
- New models need no compilation, so a nonlinear model starts fitting at once.
  With Stan, each new nonlinear model took minutes to compile.
- `backend = 'stan'` fits as ctsem 3 did.

## Random effects

- Random effects are integrated with a Laplace approximation. Version 3 carried
  them as extra states of the Kalman filter, which underestimates individual
  differences in the dynamics. In linear models, random intercepts are
  integrated exactly. `intoverpop = 'augmented'` gives the old method, and
  `ctFit()` warns when that method meets nonlinear dynamics.
- Each fit checks the approximation and warns if it affects the estimates.
  `intoverpop = 'quadrature'` then refines the fit, and
  `ctFitUncertainty(fit, 'sample')` samples the full posterior instead.
- Random effects can be at multiple nested levels, such as people within
  couples: `id = c('person', 'couple')`, with each parameter's `indvarying`
  naming the levels it varies at. See the
  [multilevel vignette](https://cdriveraus.github.io/ctsem/articles/laplace-multilevel.html).
- `poprank`, `popunique` and `popblocks` give the random effects covariance
  fewer parameters when many parameters vary.

## Measurement models

Observed variables can be continuous, binary, ordinal, counts or censored, in
any mix, set for each variable with `manifesttype`.

| `manifesttype` | type | also needed |
|---|---|---|
| 0 | continuous | |
| 1 | binary | optional `ASYMPTOTES` for guessing and ceiling |
| 2 | ordinal | `ncategories` |
| 3 | count | |
| 4 | censored | `censormin`, `censormax` |

Nonlinear measurement and dynamics are filtered more accurately than in version
3. See the vignettes on
[measurement models](https://cdriveraus.github.io/ctsem/articles/measurement-models.html),
[item response models](https://cdriveraus.github.io/ctsem/articles/irt-items.html)
and [ordinal or continuous data](https://cdriveraus.github.io/ctsem/articles/ordinal-or-continuous.html).

## Covariates and weights

- Smooth effects of time independent covariates:
  `TIpredNames = ~ s(age) + sex`.
- Subject weights, for population estimates from a sample that does not
  represent the population: `ctFit(..., weights = w)`. See the
  [sampling weights vignette](https://cdriveraus.github.io/ctsem/articles/sampling-weights.html).

## Uncertainty and sampling

- Standard errors come from the Hessian, as before, and each fit checks whether
  that approximation holds.
- `ctFitUncertainty()`, which replaces `ctOptimUncertainty()`, computes the
  uncertainty by other methods, such as a sandwich estimator, the bootstrap or
  importance sampling.
- `ctFit(optimize = FALSE)`, or `ctFitUncertainty(fit, 'sample')` after a fit,
  samples the posterior until the effective sample size reaches its target. See
  the [uncertainty vignette](https://cdriveraus.github.io/ctsem/articles/uncertainty.html).
- When sampling, ctsem checks whether the first draws are centred away from
  the estimate, which can mean the optimiser stopped at a local maximum. If
  optimising from there reaches a higher log posterior,
  `ctFit(optimize = FALSE)` fits again from that point and
  `ctFitUncertainty(fit, 'sample')` gives a warning.

## Checking and improving models

| function | use |
|---|---|
| `ctReport(fit)` | run the main summaries, plots and checks, and save them to a folder |
| `ctIdentify(data, model)` | see which parameters the data can inform, before fitting |
| `ctScoreTest(fit)` | test whether the model is missing random effects, covariate effects or nonlinearity, without refitting ([vignette](https://cdriveraus.github.io/ctsem/articles/score-tests.html)) |
| `ctFitAdequacy(fit)` | test the model's one-step-ahead predictions |
| `ctModelImprove(fit)` | search for a better model step by step, guided by score tests |
| `ctModelImprove(fit, reduceOnly = TRUE)` | simplify the parts of a model the data do not identify: parameters the data do not determine, or a random effects covariance of too high a rank |
| `ctLOO(fit, method = 'psis')` | cross validation without refitting |
| `ctFitProfile(fit)` | profile likelihood |
| `anova()`, `coef()`, `vcov()`, `logLik()` | the standard R generics |

`ctModelHigherOrder()`, `ctModelSlowComponent()`, `ctModelPeriodicComponent()`
and `ctModelUniqueProcess()` extend a model with second order dynamics, a slow
process, a cycle, or a process specific to one indicator.

## Plots

`ctNetwork()` and `ctNetworkPlot()` draw temporal and contemporaneous networks,
`ctPhasePortrait()` a phase portrait, and `ctStateDependencePlot()` how a
nonlinear model's parameters change with the state. `ctVarianceDecomposition()`
splits variance into between subject, process and measurement parts, and
`ctTracePlot()` shows the path of the optimiser or sampler.

## Model specification

- A model can be written in lavaan, Mplus or brms style syntax with
  `ctModelLavaan()`, `ctModelMplus()` and `ctModelFormula()`, and
  `ctModelSyntax()` translates between them
  ([vignette](https://cdriveraus.github.io/ctsem/articles/model-syntax.html)).
- Parameter options can be named, as in
  `'cint, indvarying=TRUE, tipreds=c(age)'`. `ctPar()` writes these, and
  `ctCov()` writes covariance matrices.
- Covariance matrices use a new parameterisation (a matrix logarithm), so their
  off-diagonal parameters are not correlations. `summary()` reports the
  correlation matrices (`DIFFUSIONcor`, `T0cor`, `MANIFESTcor`).
- `ctGenerate()` simulates data from any model, including random effects and
  non-Gaussian indicators.

## Changed defaults

| | ctsem 3.11 | ctsem 4 | for the old behaviour |
|---|---|---|---|
| backend | Stan | Julia | `backend = 'stan'` |
| random effects when optimising | Kalman filter states | Laplace approximation | `intoverpop = 'augmented'` |
| covariance parameterisation | `'rawcorr'` | matrix logarithm | `model$covmattransform <- 'rawcorr'` |
| `cores` in `ctFitUncertainty()`, `ctFitCheckCov()`, `ctKalmanArray()` | 1 | `getOption('mc.cores', 2)` | `cores = 1` |
| `cores` when sampling | one per chain | `getOption('mc.cores', 2)` | `cores = 'maxneeded'` |
| ACF plots | quantile regression | weighted spline | `method = 'quantile'` |
| subjects in `ctDiscretePars()`, `ctPredictTIP()`, `ctExtract()`, `ctKalmanArray()` | by position | by id | `realid = FALSE` |

A fit saved with ctsem 3 and passed to `ctFitUpdate()` is refitted with Stan.

## Renamed and removed

| old | new |
|---|---|
| `ctOptimUncertainty()` | `ctFitUncertainty()` (old name still works) |
| `ctFitCovCheck()` | `ctFitCheckCov()` (old name still works) |
| `ctGenerateFromPriors()` | `ctGenerate(fromPriors = TRUE)` |
| `ctModelCoverage_check()` | `ctCoverageCheck()` |
| `ctCheckFit()` | `ctPostPredPlots()` and `ctFitCheckCov()` |
| `ctFitUpdateModel()`, `ctStanUpdModel()` | `ctFitUpdate()` |
| `ctGenerate(n.subjects = )` | `ctGenerate(n = )` |
| `ctDiscretePars(observational = )` | `impulseType = 'observed'` or `'unit'` |
| `ctFit(iter = , chains = , control = )` | `sampleControl` |
| `uncertainty = 'surrogate'` | `'checked'` |
| `ctFit(vb = TRUE)` | removed |

## Later 4.x releases

Nothing yet.
