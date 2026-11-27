# log1p_exp

Maps the stan function so the same code works in R.

## Usage

``` r
log1p_exp(x)
```

## Arguments

- x:

  value to use.

## Value

`log(1 + exp(x))`, elementwise, finite where `exp(x)` overflows.

## Examples

``` r
log1p_exp(-3)
```
