# Print a ctsem model

Prints a summary of the model: its type, the parameters it will
estimate, and its individually varying parameters.

## Usage

``` r
# S3 method for class 'ctStanModel'
print(x, matrices = NULL, ...)
```

## Arguments

- x:

  A model from
  [`ctModel`](https://cdriveraus.github.io/ctsem/reference/ctModel.md).

- matrices:

  Whether to print the model matrices. Defaults to TRUE for models small
  enough to read.

- ...:

  Ignored.

## Value

`x`, invisibly.

## Details

The full specification is in `x$pars`, which can be printed directly.
[`ctModelMatrices`](https://cdriveraus.github.io/ctsem/reference/ctModelMatrices.md)
gives it as matrices.
