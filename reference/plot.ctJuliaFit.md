# Plots for ctJuliaFit objects

Plots for a fit from `backend='julia'`, as
[`plot.ctStanFit`](https://cdriveraus.github.io/ctsem/reference/plot.ctStanFit.md)
gives for stan fits. `'regression'` plots the model implied regression
coefficients against the time interval
([`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md)),
`'kalman'` plots expectations
([`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md)),
`'priorcheck'` plots each population parameter's prior and posterior
([`ctPlotPosterior`](https://cdriveraus.github.io/ctsem/reference/ctPlotPosterior.md)),
and `'trace'` plots the optimiser's path to the estimate, or the chains
of a sampled fit
([`ctTracePlot`](https://cdriveraus.github.io/ctsem/reference/ctTracePlot.md)).
The stan plot types `'density'` and `'intervals'` are not available, and
`'priorcheck'` shows the posterior densities instead.

## Usage

``` r
# S3 method for class 'ctJuliaFit'
plot(x, types = "all", wait = TRUE, ...)
```

## Arguments

- x:

  Fit object from `ctFit(..., backend='julia')`.

- types:

  Vector of plot types: 'all', 'regression', 'kalman', 'priorcheck',
  'trace'. 'all' leaves out 'regression' for a discrete time model,
  'priorcheck' for a fit with no posterior draws, and 'trace' for a fit
  with nothing to trace.

- wait:

  Logical. Pause between plots?

- ...:

  Passed through to
  [`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md)
  and
  [`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md).
  Beware clashes when `types='all'`.

## Value

Nothing. Generates plots.

## Examples

``` r
# \donttest{
# Needs a working julia backend, so the example is inert where Julia is
# absent -- including on CRAN, whose check machines have none.
if (isTRUE(ctJuliaStatus()$available)) {
  data <- ctsemtestdat
  model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
    LAMBDA = matrix(1))
  model$pars$indvarying <- FALSE
  fit <- ctFit(data, model, backend = 'julia', cores = 1)
  plot(fit, wait = FALSE)
  # Julia keeps running between calls, holding what it has compiled;
  # this ends it.
  ctJuliaKill("session")
}
# }
```
