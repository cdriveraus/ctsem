# Plot score tests

Each family's strongest relaxations: their p-values adjusted within the
family (Benjamini-Hochberg), on a log scale, against .05.

## Usage

``` r
# S3 method for class 'ctScoreTest'
plot(x, family = NULL, top = 10, ...)
```

## Arguments

- x:

  A `ctScoreTest`.

- family:

  Families to show; by default every family tested.

- top:

  Relaxations shown per family, the strongest.

- ...:

  Unused.

## Value

A ggplot.
