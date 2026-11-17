# New in ctsem 4

ctsem 4 adds an estimation engine written in Julia, now the default, and most
new features build on it. Models and code written for ctsem 3 still work. This
page summarises what changed; [NEWS.md](NEWS.md) lists every change.

## Julia engine

```r
ctJuliaInstall()   # once
fit <- ctFit(data, model)
```

`ctJuliaInstall()` sets up Julia and the engine, asking before each download.
The first fit after installing or updating ctsem takes a few extra minutes
while the engine is compiled.

- Fits are two to three times faster than with Stan on one core, and use
  several cores at once (`cores`).
- Models are not compiled one by one, so a nonlinear model starts fitting
  straight away. With Stan, each new nonlinear model took several minutes to
  compile.
- `backend = 'stan'` still fits as ctsem 3 did.

## Random effects

- Random effects are integrated with a Laplace approximation. Version 3 carried
  them as extra states of the Kalman filter, which underestimates individual
  differences in the dynamics. Random intercepts in linear models are still
  integrated exactly, also beside random dynamics, where Laplace then
  approximates only the integral over the rest
  (`optimcontrol$laplace_onestep = TRUE` approximates them all together).
  `intoverpop = 'augmented'` gives the old method, and `ctFit()` warns when
  that method meets nonlinear dynamics.
- Each fit checks the approximation and warns if it affects the estimates.
  `intoverpop = 'quadrature'` then refines the fit, and
  `ctFitUncertainty(fit, 'sample')` samples the full posterior instead.
- Random effects can sit at several nested levels, such as people within
  couples: `id = c('person', 'couple')`, with each parameter's `indvarying`
  naming the levels it varies at. See the
  [multilevel vignette](https://cdriveraus.github.io/ctsem/articles/laplace-multilevel.html).
- `poprank`, `popunique` and `popblocks` simplify the random effects
  covariance when many parameters vary.

## Measurement models

Observed variables can be continuous, binary, ordinal, counts or censored, set
for each variable with `manifesttype`, and types can be mixed in one model.

| `manifesttype` | type | also needed |
|---|---|---|
| 0 | continuous | |
| 1 | binary | optional `ASYMPTOTES` for guessing and ceiling |
| 2 | ordinal | `ncategories` |
| 3 | count | |
| 4 | censored | `censormin`, `censormax` |

Measurement and dynamics that depend on the latent state are also handled more
accurately. See the vignettes on
[measurement models](https://cdriveraus.github.io/ctsem/articles/measurement-models.html),
[item response models](https://cdriveraus.github.io/ctsem/articles/irt-items.html)
and [ordinal or continuous data](https://cdriveraus.github.io/ctsem/articles/ordinal-or-continuous.html).

## Covariates and weights

- Smooth effects of time independent covariates:
  `TIpredNames = ~ s(age) + sex`.
- Missing values of time independent covariates are handled automatically.
- Subject weights, to estimate population values from an unrepresentative
  sample: `ctFit(..., weights = w)`. See the
  [sampling weights vignette](https://cdriveraus.github.io/ctsem/articles/sampling-weights.html).

## Uncertainty and sampling

- Standard errors come from the Hessian, and ctsem checks whether that
  approximation holds.
- `ctFitUncertainty()` recomputes the uncertainty by other methods, such as
  the bootstrap or importance sampling. It replaces `ctOptimUncertainty()`.
- `ctFit(optimize = FALSE)`, or `ctFitUncertainty(fit, 'sample')` after a fit,
  samples the posterior until the effective sample size is large enough. See
  the [uncertainty vignette](https://cdriveraus.github.io/ctsem/articles/uncertainty.html).

## Checking and improving models

| function | use |
|---|---|
| `ctReport(fit)` | run the main summaries, plots and checks, and save them to a folder |
| `ctIdentify(data, model)` | see which parameters the data can inform, before fitting |
| `ctScoreTest(fit)` | test whether the model is missing random effects, covariate effects or nonlinearity, without refitting ([vignette](https://cdriveraus.github.io/ctsem/articles/score-tests.html)) |
| `ctFitAdequacy(fit)` | test the model's one-step-ahead predictions |
| `ctModelImprove(fit)` | search for a better model, step by step, guided by those tests |
| `ctLOO(fit, method = 'psis')` | cross validate without refitting |
| `ctFitProfile(fit)` | profile likelihood |
| `anova()`, `coef()`, `vcov()`, `logLik()` | the usual R generics |

`ctModelHigherOrder()`, `ctModelSlowComponent()`, `ctModelPeriodicComponent()`
and `ctModelUniqueProcess()` extend a model with a second order, a slow
process, a cycle, or a process specific to one indicator.

## Plots and summaries

`ctNetwork()` and `ctNetworkPlot()` draw temporal and contemporaneous networks,
`ctPhasePortrait()` a phase portrait, `ctStateDependencePlot()` how a nonlinear
model's parameters change with the state, and `ctVarianceDecomposition()`
splits variance into between person, process and measurement parts.
`ctTracePlot()` shows the optimisation or sampling. `summary()` also gives the
correlation matrices (`DIFFUSIONcor`, `T0cor`, `MANIFESTcor`).

## Model specification

- Parameter options can be named, as in
  `'cint, indvarying=TRUE, tipreds=c(age)'`. `ctPar()` writes these, and
  `ctCov()` writes covariance matrices.
- Covariance matrices use a new parameterisation (a matrix logarithm), so their
  off-diagonal parameters are not correlations. The correlation matrices in
  `summary()` give the correlations.
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
