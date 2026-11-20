# Inverse logit

Maps the stan function so the same code works in R.

## Usage

``` r
inv_logit(x)
```

## Arguments

- x:

  value to calculate the inverse logit for.

## Value

`1 / (1 + exp(-x))`, elementwise.

## Examples

``` r
inv_logit(-3)
```
