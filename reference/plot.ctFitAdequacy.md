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

## Examples

``` r
# \donttest{
if (isTRUE(ctJuliaStatus()$available)) {
model <- ctModel(LAMBDA = diag(2),
  manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
  TDpredNames = 'TD1')
fit <- ctFit(ctsemtestdat, model, priors = TRUE)
adequacy <- ctFitAdequacy(fit)

# One point per family of checks at its adjusted p-value. A family left of
# the dashed line at alpha, drawn in red, detects misfit.
print(plot(adequacy))

# The correlations of the one-step-ahead errors of each pair of variables,
# by the time gap between rows, each with its 95 per cent interval.
# Correlations near zero mean the model accounts for how the variables
# follow each other over time. A red point's interval excludes zero, but
# the family's adjusted p-value, which allows for the number of
# correlations, is the test.
print(plot(adequacy, family = 'residuals'))
}
# }
```
