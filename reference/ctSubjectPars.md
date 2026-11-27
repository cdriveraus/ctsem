# Extract an array of subject specific parameters from a fit

Extract an array of subject specific parameters from a fit from
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), on
either backend. `ctStanSubjectPars` is maintained as a
backward-compatible alias.

## Usage

``` r
ctSubjectPars(
  fit,
  pointest = TRUE,
  cores = getOption("mc.cores", 2L),
  nsamples = "all",
  level = NULL,
  uncertainty = NULL
)
```

## Arguments

- fit:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).

- pointest:

  if TRUE, returns only the set of individual difference parameters
  based on the max a posteriori estimate (or the median if sampling
  approaches were used).

- cores:

  Number of cores to use. `backend='stan'` only, and an error for a
  julia fit, which already holds the subject matrices.

- nsamples:

  Number of samples to calculate parameters for. Not used if
  pointest=TRUE.

- level:

  For a model with nested ids (`id = c('person', 'couple')` in
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)),
  the level to report. The default is the subject level (the first id).
  A grouping level above it gives each group's own parameter values,
  which are the population values shifted by that group's random effects
  and those of the levels above it, without subject level effects or TI
  predictor effects. Grouping levels need a `backend='julia'` fit, and
  their draws follow `uncertainty` as for subjects.

- uncertainty:

  For `backend='julia'` fits with `pointest=FALSE`, what the draws
  include. `'total'` (the default) pairs each draw of the population
  parameters with a draw of each subject's random effects from their
  posterior given that subject's data. The spread then includes both the
  uncertainty of the population parameters and how well the subject's
  own data determine its parameters. `'conditional'` keeps the
  population parameters at the estimate and draws only the subject's
  effects. `'population'` varies only the population parameters, with
  each subject at its estimate, which understates a subject's
  uncertainty when its own data are weak. A fit that sampled its random
  effects (`optimize=FALSE`) only supports `'total'`.

## Value

an nsamples by nsubjects by nparams array, with the subjects (or groups,
for a grouping `level`) labelled by their ids in the data.

## Details

This function returns the estimates of individual parameters, taking
into account any covariates and random effects, at every level of a
nested model.

A julia fit also stores, for every level with random effects, each
group's posterior for its effects on the raw parameter scale, in
`fit$randomeffects`: one entry per level, with each group's `mean`
deviation from the population value and its covariance packed column by
column (`cov`). The draws with `uncertainty='total'` come from this
posterior.

## Examples

``` r
# Each subject's values of the parameters that vary over subjects, at the
# fit's estimate. The array has one row for the estimate, a column for each
# subject (named by its id) and a slice for each parameter.
indpars <- ctSubjectPars(ctstantestfit)
dimnames(indpars)$param
indpars[1, 1:5, ]

# The continuous time intercepts of the two processes, one point per
# subject.
plot(indpars[1, , 'cint1'], indpars[1, , 'cint2'],
  xlab = 'cint1', ylab = 'cint2')

# \donttest{
if (isTRUE(ctJuliaStatus()$available)) {
model <- ctModel(LAMBDA = diag(2),
  manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
  TDpredNames = 'TD1')
fit <- ctFit(ctsemtestdat, model, priors = TRUE)

# Draws show how well each subject's values are known, with the draw in
# the first dimension. The standard deviation of the draws of mm_Y1 for
# the first five subjects:
draws <- ctSubjectPars(fit, pointest = FALSE, nsamples = 100)
print(apply(draws[, 1:5, 'mm_Y1'], 2, sd))

# The default uncertainty = 'total' includes the uncertainty of the
# population parameters. 'conditional' holds them at the estimate, so the
# spread is smaller.
conditional <- ctSubjectPars(fit, pointest = FALSE, nsamples = 100,
  uncertainty = 'conditional')
print(apply(conditional[, 1:5, 'mm_Y1'], 2, sd))
}
# }
```
