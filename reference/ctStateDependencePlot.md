# How model matrix cells vary with a latent state

Plots every context-dependent matrix cell against one latent process,
over the range that process covers, with the posterior spread at each
point.

[`summary()`](https://rdrr.io/r/base/summary.html) reports each cell as
a single number at a single state. This plot shows the same cell as a
function of the state, with the summary's number as one point on the
curve.

## Usage

``` r
ctStateDependencePlot(
  fit,
  along = 1,
  extent = "data",
  gridsize = 25,
  nsamples = 20,
  quantiles = c(0.025, 0.5, 0.975),
  state = "asymptotic",
  tipreds = NULL,
  tolerance = 1e-06,
  plot = TRUE
)
```

## Arguments

- fit:

  A ctJuliaFit.

- along:

  The latent process to vary, by name or index.

- extent:

  `'data'` (the default) spans the smoothed latent states, `'sd'` spans
  two stationary standard deviations either side of T0MEANS, or supply a
  numeric range of length two.

- gridsize:

  Points along the range. One engine call each.

- nsamples:

  Posterior draws for the interval, or 1 for the point estimate alone.

- quantiles:

  Three quantiles: lower ribbon, line, upper ribbon.

- state:

  Where the other latent processes are held: `'asymptotic'`, `'mean'`,
  `'T0MEANS'`, or a numeric state.

- tipreds:

  Time independent predictor values, or NULL for the population.

- tolerance:

  A cell is treated as varying when its range across the grid exceeds
  this fraction of its own magnitude.

- plot:

  If FALSE, return the values as a data frame rather than a plot.

## Value

A ggplot, or a data frame when `plot=FALSE`.

## Details

Every model matrix is evaluated over the grid, and a cell is shown when
its value changes. This includes cells that depend on the state
indirectly, through a PARS entry or through a covariance built from a
state dependent cell.

For a model with no state dependent cells, the function says so and
returns no plot.

## See also

[`ctContextDependence`](https://cdriveraus.github.io/ctsem/reference/ctContextDependence.md),
[`ctPhasePortrait`](https://cdriveraus.github.io/ctsem/reference/ctPhasePortrait.md)

## Examples

``` r
# \donttest{
# Needs a working julia backend, so the example is inert where Julia is
# absent -- including on CRAN, whose check machines have none.
if (isTRUE(ctJuliaStatus()$available)) {
  set.seed(1)
  generating <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2), DRIFT = matrix(c(-.4, .1, 0, -.3), 2, 2),
    CINT = matrix(c(.2, .1), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c(.5, 0, 0, .4), 2, 2)))
  datalong <- as.data.frame(suppressMessages(ctGenerate(generating, n = 25,
    burnin = 5, dtmean = 1, logdtsd = .1, wide = FALSE, Tpoints = 12)))

  # eta1's own decay depends on where eta2 is -- a state-dependent DRIFT cell.
  model <- suppressMessages(ctModel(type = 'ct', n.latent = 2, n.manifest = 2,
    manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
    LAMBDA = diag(2), PARS = c('dr11|-log1p_exp(param)'),
    DRIFT = matrix(c('dr11 * (1 + 0.2 * eta2)', 'd21', 0, 'd22'), 2, 2),
    CINT = matrix(c('c1', 'c2'), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(.2, 2), DIFFUSION = matrix(c('df1', 0, 0, 'df2'), 2, 2)))
  model$pars$indvarying <- FALSE

  fit <- ctFit(datalong, model, backend = 'julia', cores = 1, verbose = 0)
  ctStateDependencePlot(fit, along = 'eta2', gridsize = 11, nsamples = 5)
  # Julia keeps running between calls, holding what it has compiled;
  # this ends it.
  ctJuliaKill("session")
}
# }
```
