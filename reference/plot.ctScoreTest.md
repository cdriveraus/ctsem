# Plot score tests

Each family's strongest relaxations: their p-values adjusted within the
family (Benjamini-Hochberg), on a log scale, against .05.

## Usage

``` r
# S3 method for class 'ctScoreTest'
plot(x, family = NULL, top = 10, ...)
```

## Arguments

- x:

  A `ctScoreTest`.

- family:

  Families to show; by default every family tested.

- top:

  Relaxations shown per family, the strongest.

- ...:

  Unused.

## Value

A ggplot.

## Examples

``` r
# \donttest{
if (isTRUE(ctJuliaStatus()$available)) {
# A model in which the innovations of the two processes are uncorrelated.
model <- ctModel(LAMBDA = diag(2),
  manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
  TDpredNames = 'TD1',
  DIFFUSION = matrix(c('diff_eta1', 0, 0, 'diff_eta2'), 2, 2))
fit <- ctFit(ctsemtestdat, model, priors = TRUE)
st <- ctScoreTest(fit)

# A panel per family, with a point per relaxation at its adjusted
# p-value. Points left of the dashed line at .05, drawn in red, are the
# relaxations the data support, here freeing DIFFUSION[eta2,eta1].
print(plot(st))

# Two of the families, with the three strongest relaxations of each.
print(plot(st, family = c('fixed', 'heterogeneity'), top = 3))
}
# }
```
