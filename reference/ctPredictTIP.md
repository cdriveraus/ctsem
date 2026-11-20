# ctPredictTIP

Outputs the estimated effect of time independent predictors (covariate
moderators) on the expected observations.

## Usage

``` r
ctPredictTIP(
  sf,
  tipreds = "all",
  subject = NULL,
  timestep = "auto",
  doDynamics = TRUE,
  plot = TRUE,
  quantiles = c(0.16, 0.5, 0.84),
  discreteTimeQuantiles = c(0.025, 0.5, 0.975),
  dynamicsControl = list(),
  showUncertainty = TRUE,
  TIPvalues = NA,
  realid = TRUE
)
```

## Arguments

- sf:

  A fitted object from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md), from
  any backend.

- tipreds:

  A character vector specifying which time independent predictors to
  use. Default is 'all', which uses all time independent predictors in
  the model.

- subject:

  The subject for which predictions are made, as an id in the data (see
  `realid`), or NULL (the default) for the first subject. This is
  relevant only when time dependent predictors are also included in the
  model.

- timestep:

  A numeric value specifying the time step for predictions. Default is
  'auto', which tries to automatically determine an appropriate time
  step.

- doDynamics:

  A logical value indicating whether to plot the effects of time
  independent predictors on the dynamics of the system. Default is TRUE.
  Can be problematic for systems with many dimensions.

- plot:

  A logical value indicating whether to ggplot the results instead of
  returning a data.frame of predictions. Default is TRUE.

- quantiles:

  A numeric vector specifying the quantiles of the time independent
  predictors to plot. Default is c(.16,.5,.84): the median, and one
  standard deviation either side for a normally distributed predictor.

- discreteTimeQuantiles:

  a numeric vector of length 3 specifying the quantiles of the discrete
  time points to plot, when showUncertainty is TRUE.

- dynamicsControl:

  A named list of additional arguments to pass to
  [`ctDiscretePars`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md)
  and
  [`ctDiscreteParsPlot`](https://cdriveraus.github.io/ctsem/reference/ctDiscreteParsPlot.md)
  when `doDynamics=TRUE`. Arguments matching `ctDiscreteParsPlot` are
  automatically routed to the plot call. Internally controlled arguments
  `fit`, `ctstanfitobj`, `plot`, `subjects`, `realid`, `impulseType`
  (and its deprecated spelling `observational`), `x`, `quantiles`, and
  `splitSubjects` are ignored if supplied.

- showUncertainty:

  A logical value indicating whether to plot the uncertainty of the
  predictions. Default is TRUE.

- TIPvalues:

  An nvalue \* nTIpred numeric matrix specifying the fixed values for
  each time independent predictor effect to plot. Default is NA, which
  instead relies on the quantiles specified in the quantiles argument.

- realid:

  If TRUE (the default), subject is read as an id in the data, and an
  integer that is not an id there as ctsem's numbering of the subjects,
  1 to N, with a message. If FALSE, it is read as that numbering.

## Value

If plot is TRUE, a list of ggplot objects showing the estimated effects
of covariate moderators. Otherwise, a data frame with the predictions.

## Details

This function estimates the effects of covariate moderators on the
expected process and observations for a specified subject in a dynamic
system. The covariate moderators are defined at the specified quantiles,
and their effects on the trajectory are plotted or returned as a data
frame.

## Examples

``` r
# Example usage:
ctPredictTIP(ctstantestfit, tipreds='all', doDynamics=FALSE, plot=TRUE)
```
