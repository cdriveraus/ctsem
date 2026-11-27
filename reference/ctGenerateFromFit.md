# Generate data from a fitted model

Adds a `$generated` object to a fit, holding random data generated from
the fit's posterior (or its estimate), on either backend.
`ctStanGenerateFromFit` is maintained as a backward-compatible alias.

## Usage

``` r
ctGenerateFromFit(
  fit,
  nsamples = 200,
  fullposterior = FALSE,
  verboseErrors = FALSE,
  cores = getOption("mc.cores", 2L),
  intoverstates = "fit",
  effects = c("population", "fitted")
)
```

## Arguments

- fit:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), on
  either backend.

- nsamples:

  Positive integer specifying number of datasets to generate.

- fullposterior:

  Logical indicating whether to sample from the full posterior (original
  nsamples) or the posterior mean.

- verboseErrors:

  if TRUE, print verbose output when errors in generation encountered.

- cores:

  Number of cpu cores to use: R worker processes for a stan fit, and for
  a julia fit the engine threads its draws are generated on.

- intoverstates:

  For `backend='julia'` fits, how the latent states are generated.
  `'fit'` (the default) follows the fit. A fit with explicit states
  (`intoverstates=FALSE`) then generates states explicitly, and every
  other fit generates each observation from the Kalman filter's
  one-step-ahead prediction, the density the fit used. A generated
  dataset's `llrow` is then the likelihood under which it was generated.

  `FALSE` simulates the model directly, whatever the fit used: each
  subject's trajectory from the process at its own parameters, and then
  each observation from its measurement model given the state, as
  [`ctGenerate`](https://cdriveraus.github.io/ctsem/reference/ctGenerate.md)
  does. The two differ for non-Gaussian indicators and for state
  dependent dynamics, where the filter is an approximation. Comparing
  data generated each way, with the same `effects`, shows the effect of
  that approximation, separately from any misfit of the model to the
  observed data. `TRUE` always generates through the filter.

  In all cases the states are generated anew.

- effects:

  For `backend='julia'` fits with random effects. `'population'` (the
  default) draws new subjects, and groups at every level, from the
  fitted population distribution at each parameter draw, so the
  generated data are the model's prediction for a new sample with the
  same design. `'fitted'` keeps the fitted subjects, each at its sampled
  effects when the fit sampled and saved them, and otherwise at its
  conditional mode. A fit with `intoverpop='augmented'` has no
  per-subject effects, so only `'population'` applies. Stan fits
  generate as the fit was made, and give an error if the argument is
  set.

## Value

The fit, with a `$generated` element: `Y`, the generated data as an
array of sample by data row by manifest variable, on the original time
and missingness structure, and `llrow`, each generated row's log
likelihood.

## Examples

``` r
gen <- ctGenerateFromFit(ctstantestfit, nsamples=3,fullposterior=TRUE,cores=1)
plot(gen$generated$Y[3,,2],type='l') #Third random data sample, 2nd manifest var, all time points. 
```
