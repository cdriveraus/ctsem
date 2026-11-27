# Plot how a fit reached its estimate

For an optimized fit, draws the log posterior and the gradient norm
against iteration, from the record the engine kept while fitting. For a
sampled fit, draws each chain's draws of every population mean and
standard deviation, on the scale
[`summary()`](https://rdrr.io/r/base/summary.html) reports them on, so
chains that have not settled on one distribution can be seen.

## Usage

``` r
ctTracePlot(fit, which = NULL, plot = TRUE)
```

## Arguments

- fit:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) with
  the Julia backend.

- which:

  Which panels to draw: trace columns (`'objective'`, `'gradient_norm'`,
  `'inner_converged'`) for an optimized fit, parameter names for a
  sampled one. Defaults to all of them.

- plot:

  If `FALSE`, return the data instead of plotting it. For an optimized
  fit this is the trace, and for a sampled fit one row per draw and
  parameter.

## Value

A ggplot, or with `plot = FALSE` a data frame.

## Details

The Julia backend records the optimisation trace during the fit and
returns it as `fit$optim$trace`. For output while a fit is still
running, see `optimcontrol$callback` in
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md). The
chains of a sampled stan fit are drawn by `plot(fit, types = 'trace')`.

## Examples

``` r
# \donttest{
# Needs a working julia backend, so the example is inert where Julia is
# absent -- including on CRAN, whose check machines have none.
if (isTRUE(ctJuliaStatus()$available)) {
  model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
    LAMBDA = matrix(1))
  model$pars$indvarying <- FALSE
  fit <- ctFit(ctsemtestdat, model, backend = 'julia', cores = 1)
  print(ctTracePlot(fit))
  # Julia keeps running between calls, holding what it has compiled;
  # this ends it.
  ctJuliaKill("session")
}
# }
```
