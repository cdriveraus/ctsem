# Write a parameter specification

Builds the `|`-separated cell specification that ctsem's model matrices
use, from named arguments.

## Usage

``` r
ctPar(
  name,
  transform = NA,
  indvarying = NA,
  sdscale = NA,
  tipreds = NULL,
  mean = NA
)
```

## Arguments

- name:

  Parameter name. Two cells given the same name are the same parameter,
  which is how equality constraints are written.

- transform:

  Transform expression in terms of `param`, e.g. `'exp(param)'`. If
  omitted, the default for that matrix and cell is used, which differs
  between, for example, variances, drift diagonals and means.

- indvarying:

  Where the parameter varies: `TRUE` at every level of the model's `id`
  hierarchy, `FALSE` at none, or the names of the levels it varies at.
  With `id = c('person', 'couple')` in
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md),
  `indvarying = 'couple'` gives a parameter that differs between couples
  but not between the persons within a couple, and `TRUE` is the same as
  `c('person', 'couple')`. With a single id the two forms mean the same
  thing. The model's `pars$indvarying` column holds the same values and
  may be edited directly.

- sdscale:

  Multiplier on the population standard deviation's prior, at every
  level where the parameter varies. Only meaningful when it varies.

- tipreds:

  Time independent predictors acting on this parameter, as a character
  vector of names for effects to estimate, e.g. `'age'` or
  `c('age', 'sex')`. Naming any predictor turns the others off for this
  cell, so the list is the full set of effects on it. Each name must be
  one of the model's `TIpredNames`, or building the model gives an
  error. `character(0)` (or `tipreds=` with nothing after it, in a cell)
  means no predictor acts on this cell. Leaving the argument out keeps
  the model's default.

- mean:

  A fixed population mean for a parameter that varies (`indvarying`), on
  the scale [`summary()`](https://rdrr.io/r/base/summary.html) reports
  population means on. Only the random effects around it are then
  estimated, as for a person-level intercept with mean zero. Written
  into the cell as `'name=value|...'`, and stored on the model in
  `fixedpopmeans`, a named vector that may also be set directly. Needs
  `backend = 'julia'`.

## Value

A single character string, suitable for a model matrix cell.

## Details

The output is the same string the `|` syntax produces. Fields left at
their default are left empty, which the parser reads as unspecified.

The same named fields can also be written directly into a model matrix
cell. A cell containing `|` is read as a `|`-separated specification.
Otherwise a cell with a field of the form `name=value` is read as a
named specification, and anything else as a plain parameter name. These
three specify the same parameter:


      MANIFESTMEANS = matrix('mm||TRUE|0.5')
      MANIFESTMEANS = matrix(ctPar('mm', indvarying = TRUE, sdscale = 0.5))
      MANIFESTMEANS = matrix('mm, indvarying=TRUE, sdscale=0.5')
      

In the named form, fields are separated by commas, the parameter name
comes first, and several predictors or levels are given as
`tipreds=c(age, sex)` or `indvarying=c(person, couple)`. Commas inside
brackets belong to the value, so `transform=pnorm(param, 0, 1)` is one
field. In the `|` form, several values are separated by commas, as in
`'cross||person,couple'`. An unknown field name is an error that names
the nearest valid field, so a misspelling such as `topreds=age` is
caught. Comparisons in a state dependent cell (`==`, `>=`, `<=`, `!=`)
and named arguments inside a function call are read as part of the
value.

## See also

[`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md) for
the matrices these cells go in.

## Examples

``` r
ctPar('mm', indvarying = TRUE)
ctPar('drift11', transform = '-exp(param)')
ctPar('mm', indvarying = TRUE, tipreds = c('age', 'sex'))
# varying between couples only, in a model with id = c('person', 'couple')
ctPar('cross', indvarying = 'couple')
# a person-level intercept whose population mean is zero
ctPar('lev', indvarying = TRUE, mean = 0)

if (FALSE) { # \dontrun{
# These two models are identical.
model <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
  LAMBDA = matrix(1),
  MANIFESTMEANS = matrix(ctPar('mm', indvarying = TRUE, sdscale = 0.5)))

model <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
  LAMBDA = matrix(1),
  MANIFESTMEANS = matrix('mm, indvarying=TRUE, sdscale=0.5'))
} # }
```
