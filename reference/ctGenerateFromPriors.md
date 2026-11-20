# Generate data from a model's priors (deprecated)

Deprecated. Use
[`ctGenerate`](https://cdriveraus.github.io/ctsem/reference/ctGenerate.md)`(fromPriors = TRUE)`,
which does the same. `ctStanGenerate` is an alias and is also
deprecated.

## Usage

``` r
ctGenerateFromPriors(
  cts,
  datastruct = NA,
  is = FALSE,
  fullposterior = TRUE,
  nsamples = 200,
  parsonly = FALSE,
  cores = getOption("mc.cores", 2L),
  backend = c("auto", "julia", "stan"),
  n.subjects = 20,
  Tpoints = NULL,
  dtmean = 1
)
```

## Arguments

- cts:

  A model from
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
  or
  [`ctModelConvertOMX`](https://cdriveraus.github.io/ctsem/reference/ctModelConvertOMX.md),
  or a fit from either backend, which supplies its model and design.

- datastruct:

  Long format data structure as used by ctsem, giving the subjects,
  times and missing observations to generate. Ignored when `cts` is a
  fit, which uses its own data. If NA and `cts` is a model, a balanced
  design of `n.subjects` subjects observed at `Tpoints` occasions
  `dtmean` apart is used.

- is:

  Deprecated and ignored, with a warning if set.

- fullposterior:

  Generate from the full prior, or from its mode (the raw origin)?

- nsamples:

  How many samples to generate?

- parsonly:

  If TRUE, only return samples of raw parameters, don't generate data.

- cores:

  Number of cpu cores to use.

- backend:

  Which engine generates data from each parameter draw: `'auto'` (the
  default) uses `'julia'` when a julia session is available and `'stan'`
  otherwise, or `'julia'` or `'stan'`. The parameter draws come from the
  model's prior and are the same either way. Not used with
  `parsonly=TRUE`.

- n.subjects, Tpoints, dtmean:

  The design used when `cts` is a model and no `datastruct` is given.
  `Tpoints` defaults to the model's own when it has one. The default
  `n.subjects` is smaller than
  [`ctGenerate`](https://cdriveraus.github.io/ctsem/reference/ctGenerate.md)'s
  because `nsamples` datasets are returned.

## Value

List containing Y, an array of nsamples by data rows by manifest
variables, and llrow, an array of nsamples by data rows log likelihoods.
With `parsonly=TRUE`, the prepared model with the prior draws in
`$stanfit$rawposterior` and `$stanfit$transformedpars`.

## Examples

``` r
# \donttest{
#generate and plot samples from prior predictive
priorpred <- ctGenerateFromPriors(cts = ctstantestfit,cores=2,nsamples = 50)
# backend = 'auto' generates with julia where it is available. Julia keeps
# running between calls, holding what it has compiled; this ends it.
if (isTRUE(ctJuliaStatus()$available)) ctJuliaKill("session")
# }
```
