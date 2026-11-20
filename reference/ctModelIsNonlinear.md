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
ctModelIsNonlinear(ctstantestfit)
```
