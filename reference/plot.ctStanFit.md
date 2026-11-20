# plot.ctStanFit

Plots for ctStanFit objects

## Usage

``` r
# S3 method for class 'ctStanFit'
plot(x, types = "all", wait = TRUE, ...)
```

## Arguments

- x:

  Fit object from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).

- types:

  Vector of character strings defining which plots to create. 'all'
  plots all possible types, including: 'regression', 'kalman',
  'priorcheck', 'trace', 'density','intervals'.

- wait:

  Logical. Pause between plots?

- ...:

  Arguments to pass through to the specific plot functions. Bewar of
  clashes may occur if types='all'. For details see the specific
  functions generating each type of plot.

## Value

Nothing. Generates plots.

## Details

This function is just a wrapper calling the necessary functions for
plotting - it may be simpler in many cases to access those directly.
They are:
[`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md),[`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md),
[`ctPlotPosterior`](https://cdriveraus.github.io/ctsem/reference/ctPlotPosterior.md),`stan_trace`,
`stan_dens`,`stan_plot` rstan offers many plotting possibilities not
available here, to use that functionality one must simply call the
relevant rstan plotting function. Use `x$stanfit` as the stan fit object
(where x is the name of your ctStanFit object). Because a ctStanFit
object has many parameters, the additional argument
`pars=ctRawParnames(x,'pop_')` is recommended. This denotes population
means, but see
[`ctRawParnames`](https://cdriveraus.github.io/ctsem/reference/ctRawParnames.md)
for other options.

## Examples

``` r
# \donttest{
plot(ctstantestfit,types=c('regression','kalman','priorcheck'), wait=FALSE)
# }
```
