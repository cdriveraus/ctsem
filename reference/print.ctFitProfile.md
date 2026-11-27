# Print a profile likelihood

Print a profile likelihood

## Usage

``` r
# S3 method for class 'ctFitProfile'
print(x, digits = 3, ...)
```

## Arguments

- x:

  a `ctFitProfile`.

- digits:

  significant digits for the reported values.

- ...:

  ignored.

## Value

`x`, invisibly.

## Examples

``` r
# \donttest{
if (isTRUE(ctJuliaStatus()$available)) {
model <- ctModel(type = 'ct', manifestNames = 'Y1', latentNames = 'eta1',
  LAMBDA = matrix(1))
model$pars$indvarying <- FALSE
fit <- ctFit(ctsemtestdat, model)
result <- ctFitProfile(fit, parameters = c('drift_eta1', 'mm_Y1'))

# Each parameter's estimate and the limits of its 95 per cent profile
# likelihood interval, on the raw scale, and whether the data identify it.
print(result)
}
# }
```
