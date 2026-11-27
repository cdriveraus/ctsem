# Write a covariance matrix as the cells a ctsem model reads

ctsem's covariance matrices (`T0VAR`, `DIFFUSION`, `MANIFESTVAR`, and
the population spreads `RAWPOPVAR` and `RAWPOPVAR_<idname>`) are
specified as a lower triangle, with standard deviations on the diagonal
and, below it, values that the model's `covmattransform` turns into
correlations. Except for two variables under `'z'`, where each value is
Fisher's z of the correlation, the correlations depend on all the values
together, so the cells for a given covariance are hard to work out by
hand. `ctCov()` computes them.

## Usage

``` r
ctCov(cov, covmattransform = "z")
```

## Arguments

- cov:

  A symmetric positive definite covariance matrix.

- covmattransform:

  The construction to write the cells for: `'z'` (the default, and what
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
  uses), `'rawcorr'` or `'cholesky'`. A `ctCov()` matrix assigned into a
  model is rewritten for that model's own construction, so this matters
  only when the cells are read directly.

## Value

The cells, a lower triangular numeric matrix of class `ctCov`, with the
covariance as attribute `covariance`. `'rawcorr'` cannot reach every
correlation matrix, and a covariance it cannot represent is refused.

## Details

Under `'z'`, a cell of `T0VAR`, `DIFFUSION` or `MANIFESTVAR` that is
fixed to 0 gives a covariance of zero, and `ctCov()` writes a zero
covariance as a cell of 0. Under `'rawcorr'` and `'cholesky'` a cell of
0 gives a zero covariance only in some positions. The cells of
`RAWPOPVAR` are not treated this way, and `backend = 'stan'` stops for a
model with a cell fixed to 0 between variables that other covariances
connect.

An estimated cell that is 0, or near 0, does not show that the
covariance is zero when there are three or more variables. Read
covariances and correlations from the matrices that
[`summary.ctStanFit`](https://cdriveraus.github.io/ctsem/reference/summary.ctStanFit.md)
and
[`ctSummaryMatrices`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md)
report (`DIFFUSIONcov`, `DIFFUSIONcor`, `T0cov`, `T0cor`, `MANIFESTcov`
and `MANIFESTcor`), and not from the cells.

Use the result wherever a fixed matrix can go, such as
`ctModel(DIFFUSION = ctCov(S))`, `m$matrices$T0VAR <- ctCov(S)`, or
`m$RAWPOPVAR <- ctCov(S)`, always for the whole matrix. The model also
stores the covariance, so the matrix keeps its meaning if
`covmattransform` is changed later, and the R generator of
[`ctGenerate`](https://cdriveraus.github.io/ctsem/reference/ctGenerate.md),
which reads covariance cells as a Cholesky factor, generates from the
same covariance as the julia generator. A cell edited by hand afterwards
replaces what `ctCov()` wrote.

For `RAWPOPVAR`, the covariance is that of the raw, untransformed
parameters across subjects (or groups, for `RAWPOPVAR_<idname>`), which
is the scale a fit reports it on.

## Examples

``` r
S <- matrix(c(1, .6, .3,
              .6, 2, .5,
              .3, .5, 1.5), 3)
ctCov(S)

m <- ctModel(LAMBDA = diag(3), DIFFUSION = ctCov(S), Tpoints = 5)
```
