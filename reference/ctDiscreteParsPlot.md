# ctDiscreteParsPlot

Plots the output from
[`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md),
for model implied regression strengths at specified times for continuous
time models fit with
[`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md).
`ctStanDiscreteParsPlot` is maintained as a backward-compatible alias.

## Usage

``` r
ctDiscreteParsPlot(
  x,
  indices = "all",
  quantiles = c(0.025, 0.5, 0.975),
  latentNames = "auto",
  ylab = "Coefficient",
  xlab = "Time interval",
  ylim = NA,
  facets = NA,
  splitSubjects = TRUE,
  colour = "Effect",
  title = "auto",
  polygonalpha = 0.1,
  ggcode = FALSE,
  line = c("estimate", "median")
)
```

## Arguments

- x:

  list object returned from
  [`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md).

- indices:

  Either a string specifying type of plot to create, or an n by 2 matrix
  specifying which indices of the output matrix to plot. 'AR' specifies
  all diagonals, for discrete time autoregression parameters. 'CR'
  specifies all off-diagonals,for discrete time cross regression
  parameters. 'all' plots all AR and CR effects at once.

- quantiles:

  numeric vector of length 3, with values between 0 and 1, specifying
  which quantiles to plot. The default plots 95% credible intervals as
  the band, and the posterior median at 50% as the line where `line`
  does not draw the estimate.

- latentNames:

  Vector of character strings denoting names for the latent variables.
  'auto' just uses eta1 eta2 etc.

- ylab:

  y label.

- xlab:

  x label.

- ylim:

  Two values, the lower and upper limits of the y axis. Intervals that
  extend past them are cut at the edge.

- facets:

  May be 'Subject' or 'Effect'.

- splitSubjects:

  if TRUE, subjects are plotted separately, if FALSE they are combined.

- colour:

  Character string denoting how colour varies. 'Effect' or 'Subject'.

- title:

  Character string. 'auto' generates automatically, NULL can be used to
  disable title.

- polygonalpha:

  Numeric between 0 and 1 to multiply the alpha of the fill.

- ggcode:

  if TRUE, returns a list containing the data.table to plot (for each
  curve and time interval, the outer `quantiles` of the draws as `ymin`
  and `ymax`, and the line as `value`), and a character string that can
  be evaluated (with the necessary arguments such as ylab etc filled
  in). For modifying plots.

- line:

  What the line shows. `'estimate'` (the default) draws the effects at
  the fit's estimate when `x` carries them (`attr(x, 'estimate')`, for
  an optimised fit), and otherwise the middle of `quantiles`. `'median'`
  draws the middle of `quantiles`. The band shows the quantiles of the
  draws.

## Value

A ggplot2 object. This can be modified by the various ggplot2 functions,
or displayed using print(x).

## Details

For a model whose dynamics depend on the state, the curves are computed
at one point in the state space, and the subtitle names the point. The
default, the population initial state with time dependent predictors at
zero, can be far from the data. Choose the point, and whether to
linearise or simulate, with `ctDiscretePars(state = , method = )`.

## Examples

``` r
data.table::setDTthreads(1) #ignore this line
x <- ctDiscretePars(ctsemtestfit)
ctDiscreteParsPlot(x, indices='CR')

#to modify plot:
g <- ctDiscreteParsPlot(x, indices='CR') +
 ggplot2::labs(title='My ggplot modification')
print(g)
```
