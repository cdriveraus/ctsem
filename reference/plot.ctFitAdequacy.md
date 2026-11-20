# Plot the adequacy checks

Without `family`, plots each family's Holm adjusted p-value on a log
scale, against `alpha`. With a family, plots its statistics. For the
families over time gaps (`residuals`, `latent`, `volatility`, `scores`)
this is each pair's estimate by gap with a 95 per cent interval, for
`periodic` each variable's z by period, and for the others the 40
largest z. The intervals and z values are for each statistic alone. The
family's p-value, which allows for the number of statistics, is the
test.

## Usage

``` r
# S3 method for class 'ctFitAdequacy'
plot(x, family = NULL, ...)
```

## Arguments

- x:

  A `ctFitAdequacy`.

- family:

  One family, or NULL for all of them.

- ...:

  Unused.

## Value

A ggplot.
