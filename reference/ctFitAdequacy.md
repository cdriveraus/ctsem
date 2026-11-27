# Adequacy of a fitted model, as tests

Tests a fitted model in families, each aimed at one way a continuous
time model can fall short, and gives each family a p-value. Most
families use the one-step-ahead prediction errors, which under the model
are independent and standard normal however the observations are spaced.
Each variable's error is taken given the others observed at the same
time, so a result can be tied to a variable.

## Usage

``` r
ctFitAdequacy(
  fit,
  covariates = NULL,
  families = c("residuals", "latent", "volatility", "pairs", "periodic", "distribution",
    "fit", "dependence", "effects", "tvcovariates", "subjects", "nonlinear",
    "covariates"),
  gaps = NULL,
  periods = NULL,
  tvcovariates = NULL,
  close = 0.05,
  alpha = 0.05,
  projection = c("auto", "derivative", "plugin", "opg", "none"),
  parameters = NULL,
  nsamples = 100,
  lags = 0:3,
  nboot = 1999L,
  cores = getOption("mc.cores", 2L)
)
```

## Arguments

- fit:

  A fit from
  [`ctFit`](https://cdriveraus.github.io/ctsem/reference/ctFit.md),
  either backend. The `projection` needs `backend = 'julia'`.

- covariates:

  Optional subject level covariates, as for
  [`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md),
  for the `covariates` family.

- families:

  Which families to run; see Details. `scores` and `lagged` run only
  when named here.

- gaps:

  Break points, in the data's time units, for the time gap bins of the
  residual correlations. By default the eighth, quarter and half
  quantiles of the time between every pair of a subject's rows.

- periods:

  Periods, in the data's time units, for the `periodic` family, such as
  24 for a daily cycle with time in hours. By default a range of periods
  is scanned.

- tvcovariates:

  Optional time varying covariates for the `tvcovariates` family: a data
  frame with the model's subject id and time columns and one or more
  numeric columns, matched to the fit's rows by id and time.

- close:

  The root mean square residual correlation regarded as a close fit.

- alpha:

  Level for `adequate` and `close_fit`.

- projection:

  How the tests based on the errors account for the estimation of the
  parameters. `'derivative'` differentiates the statistics with respect
  to every parameter (two passes of the filter per parameter) and uses
  the fit's Hessian. `'plugin'`, on a Laplace fit, does the same more
  cheaply, with each subject's random effects estimated once from the
  whole record, and is approximate. `'opg'` uses the outer product of
  the scores, which needs many more subjects than parameters. `'none'`
  ignores estimation, which is conservative and can hide a strong misfit
  in the `distribution` family. When random effects affect the dynamics
  or the variances, `'derivative'` and `'plugin'` correct only the
  `distribution`, `fit`, `effects` and `scores` families. `'auto'` (the
  default) takes the first that is affordable: `'derivative'` when the
  fit has a Hessian and the passes take under
  `getOption('ctsem.adequacy.budget', 120)` seconds, then `'plugin'` on
  the same terms, `'opg'` with at least ten subjects per parameter, and
  otherwise `'none'`.

- parameters:

  For the `scores` family, the raw parameters whose row scores are used.
  By default up to 16 of the dynamic and measurement parameters, those
  in DRIFT, DIFFUSION and CINT first.

- nsamples:

  Number of generated datasets, for the families that need them.
  Generated data already on `fit` are reused if there are enough.

- lags:

  Observation lags for the `lagged` family.

- nboot:

  Bootstrap draws for each family's calibration.

- cores:

  Cores for the engine and for generating data.

## Value

An object of class `ctFitAdequacy`, containing: `families`, one row per
family with the number of statistics, the worst one, its z, the family's
p-value and the Holm adjusted p-value; `statistics`, every statistic,
with what it is about in the columns `variable`, `with` (the other
variable of a pair, or what the variable is correlated with), `gap` (the
time gap bin or lag), `period` and `covariate`, each NA where it does
not apply, and its observed value, z and two-sided p; `size`, the root
mean square residual correlations with their intervals and close fit
tests; `errors`, the rows the error based families use, in long form,
with `Subject`, `Time`, `Element` (`'error'` for each variable's
standardised one-step-ahead error, `'update'` for each latent's
standardised update, `'etaprior'` for its predicted mean), `Row`, the
variable and `value`, or NULL when no such family ran; `adequate`;
`close_fit`; `informative`, whether the `tvcovariates` predict the
errors at `alpha` (NA without them); `projection`, the method used;
`degenerate`, the axes the `effects` family left out; and `fit`, with
any generated data, for reuse.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
families, or one family's statistics.

## Details

- residuals:

  Correlations between the errors of a subject's rows, for every pair of
  rows, by time gap (`gaps`): the residual autocorrelation and
  cross-correlation in continuous time. Slow structure missing from the
  model shows as correlation over the gaps it spans.

- latent:

  The same for each latent process, using the standardised update each
  row makes to the predicted state. There is one set per process, so
  misfit in the dynamics is not diluted by item specific structure.

- volatility:

  The same for the squared errors of each variable: variance that
  clusters in time.

- pairs:

  The correlation of each pair of variables' errors in the same row,
  against the correlation the model implies. This detects a wrong
  measurement error covariance or loading.

- periodic:

  A cycle missing from the model, such as a daily rhythm. Each
  variable's errors, and each latent's updates, are correlated with the
  same subject's at other times, weighted by the cosine of the time
  between them over the period. This is the periodogram of the errors at
  that period. It depends only on the time between rows, so it finds a
  cycle whatever its phase in each subject. It uses the given `periods`,
  or by default scans from twice the typical interval between rows to
  half the typical record length. The estimate is the errors'
  autocorrelation at the period, and a cycle with amplitude a, in error
  standard deviations, gives a^2 / 2.

- distribution:

  Skewness and kurtosis of each variable's errors, and each row's
  squared standardised distance against the mean and variance of its
  chi-square distribution, with Anderson-Darling tests against the
  normal and chi-square.

- fit:

  How unevenly subjects are fitted: the spread of each subject's summed
  distance against the chi-square's, and the worst subject.

- dependence:

  The distance correlation of each variable's error, and of its square,
  with the predicted latent states, the previous errors and the time
  interval, tested by permutation. This detects a conditional mean or
  variance the model gets wrong, in any form.

- effects:

  On a single level Laplace fit, whether a normal distribution suits the
  random effects. At cuts along each axis, each subject's smoothed
  posterior probability that its standardised effect lies below the cut
  is compared with the normal's. Under the model these differences have
  mean zero, and a skewed or bimodal distribution, or a subgroup, moves
  them (based on the gradient function of Verbeke & Molenberghs, 2013,
  and Drikvandi, Verbeke & Molenberghs, 2017). An effect the data leave
  at its prior in every subject is left out (`degenerate`). Fits with
  more than four random effects per subject are not tested, with a
  message.

- tvcovariates:

  With `tvcovariates`, whether time varying variables not in the model
  predict its errors. Each variable's error, and each latent's update,
  is regressed on each standardised covariate at the same row and at the
  subject's previous row. The previous row tests prediction, in
  Granger's sense. The same row tests association, which also shows when
  the outcome itself moves the covariate. A correct model of the
  outcomes fails this test when a covariate is another measure of the
  same latent states, since it sharpens the estimate of the state. So a
  failure says that a model using these variables would predict better,
  and not that the model of the outcomes is wrong. This family is not
  counted in `adequate`.

- scores:

  Run only when named in `families`. White's (1987) dynamic information
  matrix test: the score of each row's one-step log likelihood for the
  dynamic and measurement parameters (`parameters`), correlated with the
  scores of the same subject's other rows, by time gap. Under the model
  these are uncorrelated, so this checks for dynamic misspecification of
  any kind. Because it is so broad, it is the least powerful family, and
  [`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
  detects the departures it covers far more often. It can take hours on
  large data.

- subjects, nonlinear, covariates, lagged:

  Statistics of the data compared with datasets generated from the fit:
  how much subjects differ in mean, standard deviation and lag one
  correlations, whether change depends nonlinearly on level, with
  `covariates` whether subjects' statistics depend on them, and the
  data's own lagged correlations (`lags`). These use the data twice and
  are conservative. The lagged correlations are usually reproduced by
  the fit, so `lagged` runs only when named in `families`.

The tests based on the errors account for the estimation of the
parameters (`projection`), and take their variance from the errors
themselves, so they hold with heavy tails or changing variance. The
`periodic` statistics are skewed, and are referred to a chi-square
matched to their skewness. The `distribution` and `fit` families use a
studentised bootstrap over subjects. Within a family, the largest
statistic and the sum of the squared statistics are each tested, taking
the statistics' correlation into account, and combined with the Cauchy
combination test. Families are combined by Holm's method, and `adequate`
is `TRUE` when no family's adjusted p-value is below `alpha`. The
`tvcovariates` family is left out of this, and its p-value is reported
in `informative`.

On a Laplace fit (`intoverpop = 'laplace'`) whose random effects affect
the dynamics or the variances, each row's error comes from the model's
predictive distribution given the subject's earlier rows, with the
random effects integrated by importance sampling
(`getOption('ctsem.sequential.particles', 512)` quasi-random draws per
subject). The draws are fixed by `getOption('ctsem.sequential.seed')`,
so a call can be repeated, and another seed or more draws shows how much
of a result is due to the sampling. With random effects at more than one
level this is not available, so the row distances are left out and the
latent family's p-values are approximate, with a message.

For a binary, ordinal, count or censored item, the error is replaced by
its randomised probability integral transform, mapped through the normal
quantile function. Under the model these are independent standard
normal, so every family uses them as it uses Gaussian errors. The latent
family's updates on such rows rest on the filter's Gaussian
approximation for the item, so its p-values are approximate there, with
a message. The `pairs` family does not pair a non-Gaussian item with a
Gaussian one.

With many observations every misfit becomes significant, so the size of
the residual, latent and within row correlations the model leaves is
reported too, as a bias corrected root mean square (as for the SRMR,
Maydeu-Olivares 2017), with a 95 per cent interval and `p_close`, the
p-value for the test that it is at least `close`. `close_fit` is `TRUE`
when every failing family is one of these three, and each of their sizes
is shown to be below `close` at `alpha`. A size whose interval reaches
`close` is not a close fit. A close fit refers to these averages over
many correlations, and one large discrepancy, or a biased estimate of a
quantity of interest, can still hide under a small average.

A failure in one variable's equation also shows in the other variables'
errors, so a failing statistic does not always point to where the
problem is.
[`ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
is more sensitive to the departures it covers, and
[`ctModelImprove`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
uses both.

## References

Newey, W. K. (1985). Maximum likelihood specification testing and
conditional moment tests. *Econometrica, 53*, 1047-1070.

White, H. (1987). Specification testing in dynamic models. In T. Bewley
(Ed.), *Advances in Econometrics: Fifth World Congress*. Cambridge
University Press.

Verbeke, G., & Molenberghs, G. (2013). The gradient function as an
exploratory goodness-of-fit assessment of the random-effects
distribution in mixed models. *Biostatistics, 14*, 477-490.

Drikvandi, R., Verbeke, G., & Molenberghs, G. (2017). Diagnosing
misspecification of the random-effects distribution in mixed models.
*Biometrics, 73*, 63-71.

Chang, J., Yao, Q., & Zhou, W. (2017). Testing for high-dimensional
white noise using maximum cross-correlations. *Biometrika, 104*,
111-127.

Chen, S. X., & Qin, Y.-L. (2010). A two-sample test for high-dimensional
data with applications to gene-set testing. *The Annals of Statistics,
38*, 808-835.

Liu, Y., & Xie, J. (2020). Cauchy combination test. *Journal of the
American Statistical Association, 115*, 393-402.

Maydeu-Olivares, A. (2017). Assessing the size of model misfit in
structural equation models. *Psychometrika, 82*, 533-558.

## See also

[`ctPostPredPlots`](https://cdriveraus.github.io/ctsem/reference/ctPostPredPlots.md)
and
[`ctFitCheckCov`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md)
to see the data against the fit.

## Examples

``` r
# \donttest{
if (isTRUE(ctJuliaStatus()$available)) {
model <- ctModel(LAMBDA = diag(2),
  manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
  TDpredNames = 'TD1')
fit <- ctFit(ctsemtestdat, model, priors = TRUE)
adequacy <- ctFitAdequacy(fit)
print(adequacy)     # one row per family of tests
print(plot(adequacy))
}
# }
```
