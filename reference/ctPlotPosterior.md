# Prior and posterior distributions of a fit's parameters

Plots, for each free population parameter of a fit, its prior and its
posterior on the natural scale. This covers the population mean, and for
an individually varying parameter the population standard deviation at
each level it varies at. Works for fits from either backend. The prior
is ctsem's `normal(0,1)` raw prior passed through each parameter's
transform, which is the prior used with `priors = TRUE`. For a fit made
without priors, the plot shows how far the data moved each parameter
from where that prior would put it. The posterior comes from the fit's
draws, which are posterior draws for a sampled fit, and approximate ones
from
[`ctFitUncertainty`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
for an optimized fit. `ctStanPlotPost` is an alias.

## Usage

``` r
ctPlotPosterior(
  obj,
  rows = "all",
  npp = 4,
  priorwidth = TRUE,
  smoothness = 1,
  priorsamples = 10000,
  plot = TRUE,
  wait = FALSE,
  ...
)
```

## Arguments

- obj:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), on
  either backend. A model on its own has no posterior; `plot(model)`
  plots its priors.

- rows:

  Character string 'all', or which parameters to plot, by position in
  the fit's free population parameters, in parameter order.

- npp:

  Integer number of parameters to show per page.

- priorwidth:

  if TRUE, plots will be scaled to show bulk of both the prior and
  posterior distributions. If FALSE, scale is based only on the
  posterior.

- smoothness:

  Positive numeric – multiplier to modify smoothness of density plots,
  higher is smoother but can cause plots to exceed natural boundaries,
  such as standard deviations below zero.

- priorsamples:

  number of samples from prior to use. More is slower.

- plot:

  Logical, if FALSE, ggplot objects are returned in a list instead of
  plotting.

- wait:

  If true, user is prompted to continue before plotting next graph. If
  false, graphs are plotted one after another without waiting.

- ...:

  For a stan fit, passed to the preparation of the prior draws:
  `cores = x` will speed things up, where x is the number of cpu cores
  to use. A julia fit takes nothing here.

## Value

With `plot = TRUE`, `NULL` invisibly, having drawn the plots; with
`plot = FALSE`, a list of ggplot objects.

## Examples

``` r
# \donttest{
ctPlotPosterior(ctstantestfit, rows=3:4)
# }
```
