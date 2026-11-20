# Profile likelihood for a ctsem julia fit

Fixes one raw parameter at a sequence of values, re-optimising all the
others at each, and reports how the likelihood changes. This shows
whether the data determine the parameter, and gives likelihood based
intervals. The curvature at the estimate can mislead about this along a
flat direction (see Directions with no curvature in
[`ctFitUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)).

## Usage

``` r
ctFitProfile(
  fit,
  parameters = NULL,
  points = 8L,
  step = NULL,
  growth = 1.6,
  level = 0.95,
  maxiter = 500L,
  verbose = 0,
  cores = getOption("mc.cores", 2L)
)
```

## Arguments

- fit:

  a `ctJuliaFit`.

- parameters:

  which parameters to profile: raw parameter names, indices, `NULL` for
  all of them, or `"flagged"` for only those the fit's own diagnostics
  question (parameters at the edge of their transform, unidentified
  directions, and intervals wider than the curvature supports).
  Profiling every parameter of a large model takes many optimisations.

- points:

  maximum points per side.

- step:

  first displacement, in raw units. `NULL` takes half of each
  parameter's standard error where the fit has one, so an identified
  parameter gets several points inside its interval, and 0.5 otherwise.

- growth:

  factor the step is multiplied by after a point that changed the
  likelihood by less than the target, so a flat direction is walked out
  quickly.

- level:

  confidence level setting the likelihood-ratio bar.

- maxiter:

  iteration cap for each constrained optimisation.

- verbose:

  0 for silence, 1 to report each point as it is computed.

- cores:

  The most CPU cores the profile uses at once. By default the two sides
  of each parameter's profile run one after another in this session, on
  that many engine threads. When the fit itself took more than five
  minutes to optimise, they run up to `cores` at a time instead, each in
  a background R session with its own Julia session. This needs the
  future package. The sessions stay open for the next call, and
  [`ctJuliaKill`](https://cdriveraus.github.io/ctsem/reference/ctJuliaKill.md)`('workers')`
  ends them. The profile is the same either way, except that with
  several sessions every side finishes before `$better` is reported.

## Value

An object of class `ctFitProfile`, with `print` and `plot` methods.
`$profile` has one row per point (parameter, raw value, transformed
value where there is one, log likelihood, and drop from the estimate).
`$summary` has one row per parameter (verdict, limits, and how far each
side was walked). `$better` is `NULL` unless a constrained fit beat the
estimate, in which case it holds that point and the run stopped there.

## Details

Each point is a constrained optimisation, so a profile can find that the
fit was not at a maximum. When that happens the run stops and the better
point is returned in `$better`. Refit from there and profile again.

## Examples

``` r
# \donttest{
# Needs a working julia backend, so the example is inert where Julia is
# absent -- including on CRAN, whose check machines have none.
if (isTRUE(ctJuliaStatus()$available)) {
  model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
    LAMBDA = matrix(1))
  model$pars$indvarying <- FALSE
  fit <- ctFit(ctstantestdat, model, backend = 'julia', cores = 1)

  # Two named parameters; parameters = "flagged" profiles only those the
  # fit's own diagnostics doubt, and NULL profiles every one.
  result <- ctFitProfile(fit, parameters = c('drift_eta1', 'mm_Y1'))
  print(result)
  plot(result)
  # Julia keeps running between calls, holding what it has compiled;
  # this ends it.
  ctJuliaKill("session")
}
# }
```
