# Plot a profile likelihood

One panel per parameter, showing the drop in log likelihood against the
value the parameter was held at, with the likelihood-ratio bar drawn
across. A curve that crosses the bar on both sides gives a confidence
interval, and one that stays near zero shows a direction the data do not
determine.

## Usage

``` r
# S3 method for class 'ctFitProfile'
plot(x, parameters = NULL, ...)
```

## Arguments

- x:

  a `ctFitProfile`.

- parameters:

  which to draw; defaults to all that were profiled.

- ...:

  passed to `plot`.

## Value

`x`, invisibly.
