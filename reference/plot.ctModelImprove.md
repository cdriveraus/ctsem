# Plot a model improvement search

Plots the search's path, one row per refit in the order tried, from the
log likelihood its round started at to the one the refit reached, both
as the gain over the starting fit. Kept changes are filled and rejected
ones open. Pruning decisions are in the subtitle.

## Usage

``` r
# S3 method for class 'ctModelImprove'
plot(x, ...)
```

## Arguments

- x:

  A `ctModelImprove`.

- ...:

  Unused.

## Value

A ggplot.
