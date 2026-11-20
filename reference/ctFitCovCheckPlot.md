# ctFitCovCheckPlot

Plot the results of ctFitCheckCov.

## Usage

``` r
ctFitCovCheckPlot(x, maxlag = 10, vars = NA, splitvar = NA, cor = FALSE, ...)
```

## Arguments

- x:

  Data output from \`ctFitCheckCov(..., plot = FALSE)\`. If a ggplot
  object or list of ggplot objects is supplied, it is returned
  unchanged. Data with an \`x\` column (from \`ctFitCheckCov(x = )\`
  other than \`"lag"\`) are drawn against it, a panel per lag and
  variable.

- maxlag:

  Maximum lag to plot.

- vars:

  Optional character vector of variable names to plot. If \`NA\`, all
  variables are plotted.

- splitvar:

  Optional character string specifying a variable to split the plot by.

- cor:

  Logical; if \`TRUE\`, plot correlations instead of covariances. (label
  change)

- ...:

  not used.

## Value

ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
ctFitCovCheckPlot(ctFitCheckCov(ctstantestfit,cor=TRUE,plot=FALSE),maxlag=3)
} # }
```
