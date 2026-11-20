# Plots uncertainty bands with shading

Plots uncertainty bands with shading

## Usage

``` r
ctPoly(x, y, ylow, yhigh, steps = 20, ...)
```

## Arguments

- x:

  x values

- y:

  y values

- ylow:

  lower limits of y

- yhigh:

  upper limits of y

- steps:

  number of polygons to overlay - higher integers lead to smoother
  changes in transparency between y and yhigh / ylow.

- ...:

  arguments to pass to polygon()

## Value

Nothing. Adds a polygon to existing plot.

## Examples

``` r
plot(0:100,sqrt(0:100),type='l')
ctPoly(x=0:100, y=sqrt(0:100), 
yhigh=sqrt(0:100) - runif(101), 
ylow=sqrt(0:100) + runif(101),
col=adjustcolor('red',alpha.f=.1))
```
