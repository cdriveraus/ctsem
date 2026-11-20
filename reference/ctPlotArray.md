# Plots three dimensional y values for quantile plots

1st margin of \$Y sets line values, 2nd sets variables, 3rd quantiles.

## Usage

``` r
ctPlotArray(
  input,
  grid = FALSE,
  add = FALSE,
  colvec = "auto",
  lwdvec = "auto",
  ltyvec = "auto",
  typevec = "auto",
  plotcontrol = list(ylab = "Array values", xaxs = "i"),
  legend = TRUE,
  legendcontrol = list(),
  polygon = TRUE,
  polygonalpha = 0.1,
  polygoncontrol = list(steps = 25)
)
```

## Arguments

- input:

  list containing 3 dimensional array to use for Y values, `$y` and
  vector of corresponding x values `$x`.

- grid:

  Logical. Plot with a grid?

- add:

  Logical. If TRUE, plotting is overlayed on current plot, without
  creating new plot.

- colvec:

  color vector of same length as 2nd margin.

- lwdvec:

  lwd vector of same length as 2nd margin.

- ltyvec:

  lty vector of same length as 2nd margin.

- typevec:

  type vector of same length as 2nd margin.

- plotcontrol:

  list of arguments to pass to plot.

- legend:

  Logical. Draw a legend?

- legendcontrol:

  list of arguments to pass to
  [`legend`](https://rdrr.io/r/graphics/legend.html).

- polygon:

  Logical. Draw the uncertainty polygon?

- polygonalpha:

  Numeric, multiplier for alpha (transparency) of the uncertainty
  polygon.

- polygoncontrol:

  list of arguments to pass to
  [`ctPoly`](https://cdriveraus.github.io/ctsem/reference/ctPoly.md)

## Value

Nothing. Generates plots.

## Examples

``` r
#'
input<-ctTIpredEffects(ctstantestfit, plot=FALSE, whichpars='CINT',
 nsamples=10,nsubjects=10)
    
ctPlotArray(input=input)
```
