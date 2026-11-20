# Empirical Bayes subject-wise ctsem fits

Fits one ctsem model per subject using the model prior, estimates the
empirical marginal distribution of the raw parameters, then fits each
subject again using the resulting empirical Bayes prior.

## Usage

``` r
ctEmpiricalBayesFit(
  datalong,
  model,
  subjects = "all",
  priors = TRUE,
  optimize = TRUE,
  cores = getOption("mc.cores", 2L),
  subjectFitArgs = list(),
  Npasses = 2,
  ebUse = c("rawest", "rawposterior"),
  ebRobust = TRUE,
  ebOutlierMAD = 6,
  ebOutlierQuantiles = c(0.025, 0.975),
  ebWinsorize = TRUE,
  minsd = 1e-06,
  verbose = 0,
  progress = TRUE,
  ...
)
```

## Arguments

- datalong:

  Long format data containing multiple subjects.

- model:

  Model object from
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md).
  Time independent predictors are not supported. A random-effect-free
  copy of this model is used for the per-subject fits.

- subjects:

  Vector of subject identifiers to fit, or `'all'`.

- priors:

  Logical. Passed to
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md);
  defaults to `TRUE`.

- optimize:

  Logical. Passed to
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md);
  defaults to `TRUE`.

- cores:

  Number of subjects to fit in parallel. Each individual subject-level
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) call
  uses one core.

- subjectFitArgs:

  Named list of additional arguments passed to each
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) call.
  For optimized fits, `optimcontrol$stochastic` defaults to `FALSE` for
  all EB passes unless supplied here or in `...`. First-pass fits force
  `optimcontrol$estonly=TRUE`.

- Npasses:

  Total number of subject-wise fitting passes. The default `2` fits once
  with the model prior, builds a marginal empirical Bayes prior, then
  fits once with that prior. Values above `2` repeatedly map the
  previous pass estimates back to the original raw scale, rebuild the
  marginal EB prior, and refit.

- ebUse:

  `'rawest'` to build the empirical Bayes prior from first pass point
  estimates, or `'rawposterior'` to pool raw posterior samples.

- ebRobust:

  Logical. If TRUE, the empirical Bayes prior is built from robust raw
  summaries after outlier handling.

- ebOutlierMAD:

  Positive numeric. Raw values farther than this many MAD-scaled
  deviations from the median are treated as outliers. Use `Inf` or
  `NULL` to disable this rule.

- ebOutlierQuantiles:

  Length two numeric vector of lower and upper quantiles used to bound
  first-pass raw values, or `NULL` to disable.

- ebWinsorize:

  Logical. If TRUE, outlying first-pass raw values are clamped to the
  outlier bounds before computing the EB prior. If FALSE, they are set
  to missing for EB prior construction.

- minsd:

  Lower bound used for empirical raw SDs before model adjustment.

- verbose:

  Integer from 0 to 2. Passed to
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).

- progress:

  Logical. If TRUE, report the current EB fitting stage and overwrite a
  single console line with the subject fitting percentage.

- ...:

  Additional arguments passed to each
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) call.

## Value

Object of class `ctEmpiricalBayesFit`, containing the subject fit lists
and metadata. `$initialfits` contains the first-pass model prior fits,
`$fits` contains the final empirical Bayes prior fits, and `$passfits`
contains every pass. Use
[`summary()`](https://rdrr.io/r/base/summary.html) to compute final
transformed-parameter means, SDs, covariances, correlations, and outlier
diagnostics.

## Details

Works with either `backend='stan'` or `backend='julia'`, passed via
`subjectFitArgs` or `...`. Every per-subject
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) call
uses the same backend. `ebUse='rawposterior'` needs each subject fit to
carry posterior draws, which optimised fits on both backends do unless
`optimcontrol$estonly=TRUE`, or a sampled fit (`optimize=FALSE`).

## Examples

``` r
# \donttest{
if (isTRUE(ctJuliaStatus()$available)) {
  model <- ctModel(type='ct', manifestNames='Y1', LAMBDA=matrix(1))
  eb <- ctEmpiricalBayesFit(ctstantestdat, model, cores=2,
    subjectFitArgs=list(optimcontrol=list(finishsamples=20)))
  summary(eb)
}
# }
```
