# Does a model or fit have context-dependent (nonlinear) matrix cells?

A ctsem matrix cell may be an expression of a latent process or a time
dependent predictor, such as `DRIFT[1,1] = '-log1p(exp(param)) * eta2'`.
Such a cell has a value only at a given state and set of predictor
values, so summaries and plots report the point at which they evaluated
it.

## Usage

``` r
ctModelIsNonlinear(x)
```

## Arguments

- x:

  A `ctStanModel`, `ctStanFit` or `ctJuliaFit`.

## Value

Logical.

## Details

References to the state that carries an individually varying parameter
(with `intoverpop='augmented'`) do not count, since these are
parameters, reported exactly.

## Examples

``` r
# FALSE: no cell of the bundled fit's model reads a latent state or a time
# dependent predictor.
ctModelIsNonlinear(ctsemtestfit)

# TRUE for a model in which the auto effect of eta1 changes with eta2.
model <- ctModel(LAMBDA = diag(2), manifestNames = c('Y1', 'Y2'),
  latentNames = c('eta1', 'eta2'), PARS = c('a|-log1p_exp(param)', 'b'),
  DRIFT = matrix(c('a + b * eta2', 'drift21', 0, 'drift22'), 2, 2))
ctModelIsNonlinear(model)
```
