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

## Examples

``` r
model <- ctModel(LAMBDA = diag(2), manifestNames = c('Y1', 'Y2'),
  latentNames = c('eta1', 'eta2'))

# The number of free parameters, which of them vary over subjects, and each
# matrix with its parameter labels and fixed values.
model

# The same without the matrices.
print(model, matrices = FALSE)

# Every cell with its transform and other settings.
head(model$pars)
```
