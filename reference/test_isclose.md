# Tests if 2 values are close to each other

Tests if 2 values are close to each other

## Usage

``` r
test_isclose(..., tol = 1e-08)
```

## Arguments

- ...:

  values to compare. At least two, and none of them empty. Comparing a
  zero-length or NULL value is an error.

- tol:

  tolerance

## Value

Logical or testthat output.

## Examples

``` r
test_isclose(1,1.0000001, tol=1e-4)
```
