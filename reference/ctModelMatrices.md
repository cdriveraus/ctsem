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

## Examples

``` r
model <- ctModel(LAMBDA = diag(2), manifestNames = c('Y1', 'Y2'),
  latentNames = c('eta1', 'eta2'))

# The model as matrices. A label is a free parameter and a number a fixed
# value. DRIFT[1, 2] is the effect of eta2 on eta1.
mats <- ctModelMatrices(model)
mats$DRIFT

# Fix the effect of eta2 on eta1 to zero, and give both auto effects the
# same label so that they are estimated as one parameter.
mats$DRIFT[1, 2] <- 0
mats$DRIFT[2, 2] <- 'drift_eta1'
ctModelMatrices(model) <- mats
model$pars[model$pars$matrix == 'DRIFT', c('row', 'col', 'param', 'value')]

# A single cell can be changed through model$matrices.
model$matrices$DRIFT[1, 2] <- 'drift_eta1_eta2'

# The spread over subjects of the individually varying parameters, here
# the initial means and the manifest means. Standard deviations are on
# the diagonal and correlation parameters below it.
model$RAWPOPVAR
```
