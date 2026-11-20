# Print ctStanFit summaries

Print ctStanFit summaries

## Usage

``` r
# S3 method for class 'summary.ctStanFit'
print(
  x,
  width = getOption("width"),
  sections = names(x),
  row.names = TRUE,
  ...
)
```

## Arguments

- x:

  Object returned by
  [`summary.ctStanFit`](https://cdriveraus.github.io/ctsem/reference/summary.ctStanFit.md).

- width:

  Console width to use for wrapping note text and formatting tables.

- sections:

  Optional character vector of summary sections to print.

- row.names:

  Logical. Print row names for data frame sections?

- ...:

  Additional arguments passed to
  [`print`](https://rdrr.io/r/base/print.html) for individual sections.

## Value

Invisibly returns `x`.
