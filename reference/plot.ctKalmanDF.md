# Plots prediction output from ctPredict.

Plots prediction output from ctPredict.

## Usage

``` r
# S3 method for class 'ctKalmanDF'
plot(
  x,
  subjects = unique(x$Subject),
  kalmanvec = c("y", "yprior"),
  errorvec = "auto",
  errormultiply = 1.96,
  plot = TRUE,
  elementNames = NA,
  polygonsteps = 10,
  polygonalpha = 0.2,
  facets = "Variable",
  ...
)
```

## Arguments

- x:

  Output from
  [`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md).
  In general it is easier to call
  [`ctPredict`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md)
  directly with the `plot=TRUE` argument, which calls this function.

- subjects:

  which subjects to plot, as they appear in the Subject column of `x`:
  the ids in the data, unless `ctPredict` read them as ctsem's
  numbering.

- kalmanvec:

  string vector of the elements of the output to plot: 'y', the data, or
  'y', 'eta', 'err' or 'errstd' followed by 'prior', 'upd' or 'smooth'.
  'y' is the observed variables and 'eta' the latent states; 'prior' is
  conditional on all previous data, 'upd' on all data up to the current
  time point, and 'smooth' on all data. 'err' is the data minus the
  estimate, and 'errstd' its standardised form, present when the output
  was made with `standardisederrors=TRUE`. The default, 'y' and
  'yprior', plots the data against the one step ahead predictions.

- errorvec:

  vector of names indicating which kalmanvec elements to plot
  uncertainty bands for. 'auto' plots all possible.

- errormultiply:

  Numeric denoting the multiplication factor of the std deviation of
  errorvec objects. Defaults to 1.96, for 95% intervals.

- plot:

  if FALSE, plots are not generated and the ggplot object is simply
  returned invisibly.

- elementNames:

  if NA, all relevant object elements are included – e.g. if yprior is
  in the kalmanvec argument, all manifest variables are plotted, and
  likewise for latent states if etasmooth was specified. Alternatively,
  a character vector specifying the manifest and latent names to plot
  explicitly can be specified.

- polygonsteps:

  Number of steps to use for uncertainty band shading.

- polygonalpha:

  Numeric for the opacity of the uncertainty region.

- facets:

  when multiple subjects are included in multivariate plots, the default
  is to facet plots by variable type. This can be set to NA for no
  facets, or "Subject" for facetting by subject.

- ...:

  not used.

## Value

A ggplot2 object. Side effect – Generates plots.

## Examples

``` r
### Get output from ctPredict
x<-ctPredict(ctstantestfit,subjects=2,timestep=.01)

### Plot with plot.ctKalmanDF
plot(x, subjects=2)

###Single step procedure:
ctPredict(ctstantestfit,subjects=2,
  kalmanvec=c('y','yprior'),
  elementNames=c('Y1','Y2'), 
  plot=TRUE,timestep=.01)
```
