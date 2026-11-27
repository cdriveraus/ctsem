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
# \donttest{
# The lagged correlations in the data, and the same correlations in data
# generated from the fit, with one row per pair of variables and lag.
check <- ctFitCheckCov(ctstantestfit, cor = TRUE, plot = FALSE)
head(check[, c('rowvar', 'colvar', 'lag', 'empirical', 'q025', 'q975')])

# One plot per variable, with its correlation with each variable at lags
# 0 to 3. The data's correlation is a point, and the 95 per cent interval
# of the generated data's is a bar. Points outside their bar show where
# the model does not reproduce the data.
ctFitCovCheckPlot(check, maxlag = 3)

# Only the plot for Y1.
ctFitCovCheckPlot(check, maxlag = 3, vars = 'Y1')
# }
```
