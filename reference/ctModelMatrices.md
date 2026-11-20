# Matrix view for ctStanModel objects

Access or replace the matrix representation of a `ctStanModel` while
keeping the `pars` data frame as the canonical model specification.

## Usage

``` r
ctModelMatrices(x)

ctModelMatrices(x) <- value
```

## Arguments

- x:

  A `ctStanModel` object.

- value:

  A named list of matrices, typically copied from `ctModelMatrices(x)`
  or `x$matrices`.

## Value

`ctModelMatrices()` returns a named list of matrices. The replacement
form returns the updated `ctStanModel`.

## Details

`ctModelMatrices(x)` reconstructs matrices from `x$pars`. The
replacement form updates matching `matrix`, `row`, and `col` entries in
`x$pars`. Numeric matrix entries become fixed values; non-numeric
entries become parameter labels. Metadata such as transforms and
individual variation settings is preserved when possible and reset to
fit-ready defaults when fixed/free status changes.

The same view is available as `x$matrices`. Direct replacements such as
`x$matrices$DRIFT[1, 2] <- "cross"` update `x$pars`; detached copies
must be assigned back with `x$matrices <- mats`. A placeholder
`matrices` element is stored in the object so the view is visible in
`names(x)` and printed objects, but the matrix values are always derived
from `x$pars`.

The view also holds the population spread of the individually varying
parameters: `RAWPOPVAR` for the subject level and `RAWPOPVAR_<idname>`
for each grouping level, one row and column per parameter varying at
that level, also reachable as `x$RAWPOPVAR`. The diagonal holds standard
deviations, and the lower triangle values that determine the
correlations, both for the raw, untransformed parameters, the scale a
fit reports them on. A number fixes a cell, and a label leaves it free.
[`ctCov`](https://cdriveraus.github.io/ctsem/reference/ctCov.md) writes
a whole matrix from a covariance matrix, and can also be used to give
`T0VAR`, `DIFFUSION` or `MANIFESTVAR` as a covariance matrix.
