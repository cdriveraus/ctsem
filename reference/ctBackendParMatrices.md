# Model-implied parameter matrices from a julia backend fit

Materialise every model matrix – and the covariance and asymptotic
matrices derived from them – from a raw (unconstrained) parameter
vector, using the same engine code the likelihood uses.

## Usage

``` r
ctBackendParMatrices(
  fit,
  raw = NULL,
  tipreds = NULL,
  state = NULL,
  filterstate = NULL,
  time = 0,
  dt = 0,
  trim = TRUE
)
```

## Arguments

- fit:

  A `ctJuliaFit`, or a prepared model from `ctFit(..., fit=FALSE)`.

- raw:

  Raw parameter vector. Defaults to the fitted estimate.

- tipreds:

  Time-independent predictor values for the subject to materialise.
  Defaults to all zero, i.e. the population values.

- state:

  Latent state at which to evaluate state-dependent cells: one entry per
  latent process. Defaults to the population `T0MEANS`. The carrier
  states an `intoverpop` model appends hold the individually varying
  parameters and are not part of this – they always come from `raw`, so
  a cell reading one still varies across posterior draws when a state is
  supplied.

- filterstate:

  A whole state vector as the filter carries it, carriers included, for
  reporting one person's matrices at a state that person's own random
  effects produced. Mutually exclusive with `state`, which is what a
  caller choosing an evaluation point wants.

- time, dt:

  Time and time interval passed to state-dependent expressions that use
  them.

- trim:

  Report matrices over the real latent processes (the default), as
  ctsem's summaries do, rather than over the augmented state an
  `intoverpop` model filters in. `trim=FALSE` returns the matrices the
  filter actually uses, carrier states included.

## Value

A named list of matrices, with an attribute `stateDependent` giving the
cells whose values are conditional on `state`.

## Details

Internal. `ctBackend*` means internal throughout the package, so this is
not exported; the public, backend-neutral accessor for the population
matrices is
[`ctSummaryMatrices`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md).

For a linear model no cell depends on the state and the result is exact.
For a nonlinear one, the state-dependent cells are evaluated at `state`
and are conditional on it; `attr(x, 'stateDependent')` names them.

## Examples

``` r
# \donttest{
# ctBackendParMatrices(fit)$DRIFT
# }
```
