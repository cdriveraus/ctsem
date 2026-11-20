# What is the model missing? Score tests after a fit

A fitted model can be wrong in ways its estimates do not show. Subjects
may differ in a parameter the model holds fixed, a parameter may depend
on a covariate that is not in the model, or the dynamics may depend on
the level of the process.
[`ctScoreTest()`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
tests a fitted model, without refitting it, for whether such extensions
would fit better:

- **heterogeneity**: does a parameter differ between subjects, so that
  it should be individually varying?
- **covariates**: does a parameter depend on a subject level covariate
  that is not in the model? This includes the size of a random effect as
  well as its mean.
- **state dependence**: does a drift coefficient change with the level
  of a latent state, so that the dynamics are nonlinear?
- **volatility**: does the size of a process’s random shocks change with
  the level of a state, for example negative affect being more variable
  when it is high?

Each is a score (Lagrange multiplier) test, based on the slope of the
log likelihood at the fitted model in the direction of the extension. It
needs a few evaluations of the fitted model’s gradient and no new fit
for each candidate, so many candidates can be screened after one fit.
Each test also gives a rough estimate of the size of the effect. The
tests are screens that point to what to add, and refitting shows how
much it matters.

The three examples below each generate data with something the fitted
model lacks, fit the model without it, run
[`ctScoreTest()`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md),
and refit with what it found. Later sections cover a fourth family of
tests (cells the model fixes at zero), posterior predictive checks as
tests
([`ctFitAdequacy()`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)),
and a search that combines both and does the refitting
([`ctModelImprove()`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)).

``` r

library(ctsem)
set.seed(3)
```

## The process

Two latent processes, each measured by one indicator, observed 30 times
at irregular intervals for 150 subjects. `eta1` returns towards its
baseline at rate 0.6 and `eta2` at 0.4, and `eta1` pushes `eta2` up,
with a cross effect of 0.25. Each subject has their own measurement
intercepts. The generating model is written once, and each example
changes one cell of it.

``` r

generator <- function(DRIFT, MANIFESTMEANS = c("mm1|param|TRUE|1", "mm2|param|TRUE|1"),
  TIpredNames = NULL) {
  ctModel(type = "ct", silent = TRUE,
    manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), DRIFT = DRIFT, CINT = matrix(0, 2, 1),
    DIFFUSION = matrix(c(0.8, 0, 0, 0.6), 2, 2), MANIFESTVAR = diag(0.5, 2),
    MANIFESTMEANS = matrix(MANIFESTMEANS, 2, 1),
    T0MEANS = matrix(0, 2, 1), T0VAR = diag(1, 2), TIpredNames = TIpredNames)
}
```

The model fitted in every example is ctsem’s default for two latents
with one indicator each. Every drift, diffusion and measurement error
parameter is free and shared by all subjects, and the measurement
intercepts and initial states vary over subjects.

``` r

model <- ctModel(type = "ct", silent = TRUE,
  manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
  LAMBDA = diag(2))
```

## 1. Subjects who differ in a parameter

Here the subjects differ in how quickly `eta1` returns to baseline. The
drift `DRIFT[1,1]` varies between subjects around -0.6, with a standard
deviation of about 0.25.

``` r

gen1 <- generator(DRIFT = matrix(c("drift11|-log1p_exp(-param)|TRUE|0.6", 0.25,
  0, -0.4), 2, 2))
d1 <- ctGenerate(gen1, n = 150, Tpoints = 30, dtmean = 1, logdtsd = 0.5,
  backend = "julia", values = c(drift11 = -0.6, mm1 = 0, mm2 = 0))
truedrift <- attr(d1, "subjectpars")$drift11
c(mean = mean(truedrift), sd = sd(truedrift))
#>       mean         sd 
#> -0.6345464  0.2660944
d1 <- data.frame(d1)
```

The default model shares one drift among all subjects.

``` r

fit1 <- ctFit(d1, model, backend = "julia", cores = 2)
st1 <- ctScoreTest(fit1, tests = "heterogeneity")
st1
#> Score tests at the fitted model (laplace route, 150 subjects)
#> 
#> Unmodelled heterogeneity (LM: the variance with its covariances, against an
#> equal mixture of chi-square on df - 1 and df; z: the variance alone, one sided;
#> p: the stronger of the two; sd: the one-step estimate, raw scale):
#>        parameter    LM df      z      p  p_adj      sd sd_natural
#>       drift_eta1 19.40  5  2.580 0.0055 0.0292 0.19100    0.18400
#>  drift_eta1_eta2 16.30  5  1.160 0.0065 0.0292 0.16300    0.16300
#>        diff_eta2 14.20  5 -0.104 0.0175 0.0525 0.00176    0.00204
#>        diff_eta1  9.59  5  1.790 0.0420 0.0918 0.10400    0.15600
#>           mvarY2 12.00  5  1.600 0.0510 0.0918 0.05810    0.05690
#>   diff_eta2_eta1  5.00  5  0.776 0.3330 0.5000 0.09630    0.09630
#>           mvarY1  2.61  5  0.656 0.4120 0.5300 0.05140    0.04990
#>  drift_eta2_eta1  2.81  5 -0.157 0.7920 0.8900      NA         NA
#>       drift_eta2  2.48  5 -1.410 0.9020 0.9020      NA         NA
#> Covariances with: T0m_eta1, T0m_eta2, mm_Y1, mm_Y2
#> 
#> A new dimension of variation: largest eigenvalue 6.49, p = 0.0565; direction drift_eta1_eta2 1, drift_eta1 1, drift_eta2_eta1 0.54, diff_eta2_eta1 0.39, diff_eta1 0.38
```

Each row tests whether that parameter varies over subjects, in the way
[`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) would
fit it if it were made individually varying. That means its variance
together with its covariances with the effects already in the model,
here each subject’s initial state and measurement intercepts (named in
the line under the table). A variance cannot be negative, so that part
of the test is one sided, and `LM` is compared with an equal mixture of
chi-square distributions on `df - 1` and `df` degrees of freedom. The
covariances matter most for an intercept. An intercept that varies
independently of where each subject starts makes subjects drift apart
from a common start, while the same intercept correlated with the
initial state is a random intercept. `z` tests the variance alone, and
`p` is the smaller of the two tests’ p-values.

The drift of `eta1` is flagged on its variance alone (`z`). Other rows
are flagged too, mostly through their covariances (`LM` large where `z`
is small): the cross effect of `eta2` on `eta1`, and the diffusion of
`eta2`. This is expected. When a drift varies between subjects but is
fixed in the model, how each subject’s data relate to their starting
point and intercepts changes, and the scores of other parameters pick up
some of that. The test shows that the subjects differ and roughly where,
but cannot always say which of several related parameters is the source.
To find out, free the strongest candidate and test again.

The line under the table tests all directions at once, for subjects
varying along some combination of these parameters, with the covariances
with existing effects held at zero. This can find variation spread over
several parameters, which no single row isolates. Here one parameter
carries the variation, and its own row is the sharper test. When the
population covariance has reduced rank (`poprank`), the rows test
variation outside the shared dimensions and the line is not given,
because in simulations the rows found a missing dimension far more
often.

``` r

model1b <- model
model1b$pars$indvarying[model1b$pars$param %in% "drift_eta1"] <- TRUE
fit1b <- ctFit(d1, model1b, backend = "julia", cores = 2)
c(loglik_without = fit1$estimate$loglik, loglik_with = fit1b$estimate$loglik)
#> loglik_without    loglik_with 
#>      -11177.43      -11165.35
ctScoreTest(fit1b, tests = "heterogeneity")
#> Score tests at the fitted model (laplace route, 150 subjects)
#> 
#> Unmodelled heterogeneity (LM: the variance with its covariances, against an
#> equal mixture of chi-square on df - 1 and df; z: the variance alone, one sided;
#> p: the stronger of the two; sd: the one-step estimate, raw scale):
#>        parameter    LM df      z      p p_adj     sd sd_natural
#>        diff_eta2 14.30  6 -0.133 0.0345 0.276     NA         NA
#>  drift_eta1_eta2 12.00  6  0.420 0.0915 0.277 0.1130     0.1130
#>           mvarY2 11.70  6  1.610 0.1040 0.277 0.0580     0.0568
#>   diff_eta2_eta1  6.93  6  0.902 0.3120 0.623 0.1250     0.1250
#>        diff_eta1  2.59  6  0.660 0.4100 0.657 0.0356     0.0548
#>           mvarY1  4.93  6  0.309 0.6220 0.829 0.0200     0.0191
#>  drift_eta2_eta1  1.28  6 -0.561 0.9080 0.960     NA         NA
#>       drift_eta2  2.40  6 -1.410 0.9600 0.960     NA         NA
#> Covariances with: T0m_eta1, T0m_eta2, drift_eta1, mm_Y1, mm_Y2
#> 
#> A new dimension of variation: largest eigenvalue 2.57, p = 0.828; direction diff_eta2_eta1 1, mvarY2 -0.46, drift_eta1_eta2 0.28, mvarY1 0.21, drift_eta2_eta1 0.21
```

With the drift free to vary, the log likelihood rises by about 12. Once
the tests are adjusted for their number (`p_adj`) nothing else stands
out, so the other rows’ signals came from the drift. (A drift that
varies between subjects enters the filter nonlinearly, so
[`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) used
the Laplace approximation for the refit. The tests work with either.)

## 2. A covariate the model does not contain

Now each subject has a standardised `age`, which matters in two ways.
Older subjects have a higher intercept on `Y1`, and in older subjects
`eta1` pushes `eta2` harder (the cross effect `DRIFT[2,1]` rises by 0.3
per standard deviation of age). A second covariate, `w`, is unrelated to
anything.

``` r

gen2 <- generator(DRIFT = matrix(c(-0.6, "drift21||||age=0.3", 0, -0.4), 2, 2),
  MANIFESTMEANS = c("mm1|param|TRUE|1|age=0.6", "mm2|param|TRUE|1"),
  TIpredNames = c("age", "w"))
d2 <- data.frame(ctGenerate(gen2, n = 150, Tpoints = 30, dtmean = 1,
  logdtsd = 0.5, backend = "julia", values = c(drift21 = 0.25, mm1 = 0, mm2 = 0)))
covariates <- unique(d2[, c("id", "age", "w")])
head(covariates)
#>     id        age           w
#> 1    1  0.9272212  0.03814355
#> 31   2 -0.4992809 -1.77576431
#> 61   3 -1.6353253  0.67181018
#> 91   4 -0.9785608 -0.03100340
#> 121  5 -0.0775675  2.76565752
#> 151  6 -1.1776122  0.72924918
```

The model is fitted without either covariate, and the covariates are
given to
[`ctScoreTest()`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
afterwards, one row per subject.

``` r

fit2 <- ctFit(d2[, c("id", "time", "Y1", "Y2")], model, backend = "julia", cores = 2)
st2 <- ctScoreTest(fit2, covariates = covariates, tests = "covariates")
st2
#> Score tests at the fitted model (laplace route, 150 subjects)
#> 
#> Covariate age (linear): joint LM = 65.2, df = 17, p = 0.009
#>        parameter estimate      se estimate_natural      p   p_adj class
#>  drift_eta2_eta1  0.36900 0.06440           0.3690 <5e-04 0.00425  <NA>
#>            mm_Y1  0.05240 0.00958           0.5240 <5e-04 0.00425  <NA>
#>        diff_eta2 -0.06790 0.03770          -0.0888 0.0785 0.44500     B
#>           mvarY2  0.04780 0.03130           0.0426 0.1240 0.44500     B
#>   popsd_T0m_eta2 -0.05730 0.03700               NA 0.1310 0.44500  <NA>
#>           mvarY1 -0.05450 0.04070          -0.0483 0.2040 0.53300     A
#>       drift_eta1 -0.06970 0.05610          -0.0734 0.2310 0.53300  <NA>
#>         T0m_eta2  0.00982 0.00882           0.0982 0.2720 0.53300  <NA>
#>        diff_eta1  0.04200 0.03730           0.0666 0.2820 0.53300     A
#>       drift_eta2  0.05050 0.05110           0.0394 0.3310 0.55600  <NA>
#> 
#> Covariate w (linear): joint LM = 27.6, df = 17, p = 0.246
#>        parameter estimate      se estimate_natural      p p_adj
#>            mm_Y2   0.0149 0.00768           0.1490 0.0585 0.544
#>           mvarY1   0.0437 0.02310           0.0387 0.0640 0.544
#>        diff_eta1  -0.0338 0.02350          -0.0537 0.1620 0.680
#>            mm_Y1   0.0127 0.00971           0.1270 0.1860 0.680
#>  drift_eta2_eta1   0.0540 0.04500           0.0540 0.2340 0.680
#>   popsd_T0m_eta2   0.0406 0.04030               NA 0.3060 0.680
#>   diff_eta2_eta1  -0.0396 0.03970          -0.0396 0.3120 0.680
#>  drift_eta1_eta2   0.0290 0.02890           0.0290 0.3200 0.680
#>       drift_eta2   0.0429 0.04940           0.0335 0.3950 0.688
#>   popsd_T0m_eta1  -0.0361 0.04200               NA 0.4050 0.688
```

For each covariate, the joint test asks whether it matters for any
parameter. Each parameter’s row then gives the effect a refit with that
covariate would be expected to estimate, per unit of the covariate.
`age` is found, on the intercept of `Y1` and on the cross effect, at
about the sizes it was generated with, and `w` is not. The effects are
estimated jointly, so a covariate that moves one parameter is not
attributed to another parameter correlated with it. They are also
estimated with the other parameters held at their fitted values. The
model without `age` has already absorbed some of its effect, as the
intercepts’ standard deviation has grown to take in the age differences,
and letting those parameters move as well would extrapolate the model’s
curvature over that whole change.

The tests also cover the standard deviations of the random effects (the
`popsd_` rows), that is, whether subjects differ more at some ages than
others. A covariate can also be tested as a step at an unknown threshold
with `form = "threshold"`, for an effect that switches on above some
age, for example. This scans the possible cuts and allows for having
chosen the best one.

The refit includes `age` as a time independent predictor, with effects
on every parameter:

``` r

model2b <- ctModel(type = "ct", silent = TRUE,
  manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
  LAMBDA = diag(2), TIpredNames = "age")
fit2b <- ctFit(d2, model2b, backend = "julia", cores = 2)
effects <- summary(fit2b)$tipreds
effects[grepl("mm_Y1|drift_eta2_eta1", rownames(effects)), ]
#>                          mean    sd  2.5%   50% 97.5%  mode     z
#> tip_age_drift_eta2_eta1 0.340 0.038 0.267 0.341 0.414 0.342 9.014
#> tip_age_mm_Y1           0.503 0.086 0.345 0.500 0.676 0.507 5.831
```

## 3. Dynamics that depend on the state

In the last example `eta2` moderates how `eta1` returns to baseline.
When `eta2` is high, `eta1` returns more slowly. The drift of `eta1` is
`-0.6 + 0.3 * tanh(eta2)`, so it ranges between -0.9 and -0.3, and near
`eta2 = 0` it rises by about 0.3 per unit of `eta2`.

``` r

gen3 <- generator(DRIFT = matrix(c("-0.6 + 0.3 * tanh(eta2)", 0.25, 0, -0.4), 2, 2))
d3 <- data.frame(ctGenerate(gen3, n = 150, Tpoints = 30, dtmean = 1,
  logdtsd = 0.5, backend = "julia", values = c(mm1 = 0, mm2 = 0)))
```

The fitted model is linear, with each drift coefficient a constant.

``` r

fit3 <- ctFit(d3, model, backend = "julia", cores = 2)
st3 <- ctScoreTest(fit3, tests = "statedependence")
st3
#> Score tests at the fitted model (laplace route, 150 subjects)
#> 
#> State dependence of the drift:
#>  equation    LM df      p
#>      eta1 31.60  3 <5e-04
#>      eta2  4.39  3 0.2360
#>       all 39.50  6 <5e-04
#> 
#>       term equation     LM        p    p_adj estimate     se
#>  eta1*eta2     eta1 25.700   <5e-04 2.38e-06   0.2550 0.0504
#>  eta1*eta1     eta1 18.100   <5e-04 6.31e-05   0.1020 0.0240
#>  eta2*eta2     eta1  5.590 1.81e-02 3.62e-02   0.0874 0.0370
#>  eta2*eta2     eta2  3.180 7.46e-02 1.12e-01  -0.0346 0.0194
#>  eta1*eta1     eta2  0.302 5.83e-01 6.17e-01   0.0123 0.0225
#>  eta1*eta2     eta2  0.251 6.17e-01 6.17e-01       NA     NA
#>                                                             reads
#>  DRIFT[eta1,eta1] varies with eta2, or DRIFT[eta1,eta2] with eta1
#>                                 DRIFT[eta1,eta1] varies with eta1
#>                                 DRIFT[eta1,eta2] varies with eta2
#>                                 DRIFT[eta2,eta2] varies with eta2
#>                                 DRIFT[eta2,eta1] varies with eta1
#>  DRIFT[eta2,eta1] varies with eta2, or DRIFT[eta2,eta2] with eta1
```

Each term is one drift coefficient changing linearly with one state. The
term `eta1*eta2` in the equation for `eta1` stands out. Either the drift
of `eta1` changes with `eta2`, which is how the data were generated, or
the effect of `eta2` on `eta1` changes with `eta1`. To first order these
are the same term, so the table names both, and which to fit is a
substantive choice. The estimate is the change in the drift coefficient
per unit of the state. With the dependence written into the model:

``` r

model3b <- ctModel(type = "ct", silent = TRUE,
  manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
  LAMBDA = diag(2),
  DRIFT = matrix(c("dr11 + dstate * eta2", "dr21", "dr12", "dr22"), 2, 2),
  PARS = c("dr11|-log1p_exp(-param)", "dstate"))
fit3b <- ctFit(d3, model3b, backend = "julia", cores = 2)
c(loglik_without = fit3$estimate$loglik, loglik_with = fit3b$estimate$loglik)
#> loglik_without    loglik_with 
#>      -11224.14      -11202.13
summary(fit3b)$popmeans[c("dr11", "dstate"), ]
#>          mean    sd  2.5%    50%  97.5%   mode
#> dr11   -0.585 0.055 -0.70 -0.583 -0.485 -0.584
#> dstate  0.220 0.035  0.15  0.220  0.287  0.220
```

The refit’s `dstate` is a little smaller than the screen’s estimate,
because the generating `tanh` flattens away from zero and a linear term
averages over the whole range of `eta2`. The first screen also flagged
two weaker terms in the equation for `eta1`, `eta1*eta1` and
`eta2*eta2`. These are probably echoes of the same moderation, since
`eta1` and `eta2` are correlated through the cross effect. Testing the
refitted model shows nothing left once the tests are adjusted for their
number:

``` r

ctScoreTest(fit3b, tests = "statedependence")
#> Score tests at the fitted model (laplace route, 150 subjects)
#> 
#> State dependence of the drift:
#>  equation   LM df      p
#>      eta1 2.34  3 0.5180
#>      eta2 8.70  3 0.0395
#>       all 9.74  6 0.1930
#> 
#>       term equation       LM      p  p_adj estimate     se
#>  eta2*eta2     eta2 6.45e+00 0.0111 0.0665 -0.04690 0.0185
#>  eta1*eta1     eta2 3.77e+00 0.0520 0.1560 -0.04390 0.0226
#>  eta1*eta2     eta2 3.01e+00 0.0829 0.1660 -0.05210 0.0300
#>  eta1*eta1     eta1 2.16e+00 0.1420 0.2130  0.03480 0.0237
#>  eta2*eta2     eta1 1.90e-02 0.8900 1.0000 -0.00504 0.0365
#>  eta1*eta2     eta1 4.02e-08 1.0000 1.0000       NA     NA
#>                                                             reads
#>                                 DRIFT[eta2,eta2] varies with eta2
#>                                 DRIFT[eta2,eta1] varies with eta1
#>  DRIFT[eta2,eta1] varies with eta2, or DRIFT[eta2,eta2] with eta1
#>                                 DRIFT[eta1,eta1] varies with eta1
#>                                 DRIFT[eta1,eta2] varies with eta2
#>  DRIFT[eta1,eta1] varies with eta2, or DRIFT[eta1,eta2] with eta1
```

The term now in the model has no score left, and no other term stands
out.

The test evaluates each candidate term at the state predicted from
earlier observations. Because that prediction has not seen the
observation being scored, the observation’s contribution to the score
has mean zero when the model is right. Smoothed or updated state
estimates do not have this property, and in simulations, tests using
them as covariates rejected correct models nearly every time.

The polynomial terms assume the dependence is smooth enough to show in a
linear term. `tests = "instability"` makes no assumption about its form.
It takes each observation’s score, which says how much that observation
pulls on each parameter, and orders the observations by time and by each
latent’s state as predicted from earlier observations. Under the model,
the running sum of the scores behaves like a Brownian bridge in any
order. A parameter that changes over the study, or differs between low
and high levels of a state, bends it.

``` r

ctScoreTest(fit3, tests = "instability")
#> Score tests at the fitted model (laplace route, 150 subjects)
#> 
#> Stability across the data (Nyblom), by ordering:
#>    ordering    L df     p
#>        time 1.56 11 0.368
#>  state eta1 2.29 11 0.104
#>  state eta2 3.25 11 0.006
#> 
#>    ordering       parameter     L      p   p_adj
#>  state eta2      drift_eta1 1.880 <5e-04 0.00275
#>  state eta2       diff_eta1 1.800 <5e-04 0.00275
#>  state eta1           mm_Y1 0.156 0.0135 0.12600
#>  state eta1      drift_eta1 0.645 0.0275 0.12600
#>  state eta1       diff_eta1 0.508 0.0345 0.12600
#>  state eta1 drift_eta1_eta2 0.447 0.0610 0.16200
#>  state eta2           mm_Y1 0.101 0.0685 0.21100
#>  state eta1           mm_Y2 0.116 0.0735 0.16200
#>  state eta2          mvarY1 0.374 0.0905 0.21100
#>  state eta2       diff_eta2 0.356 0.0960 0.21100
```

Ordered by the predicted state of `eta2`, the drift and the diffusion of
`eta1` stand out. The moderation that generated the data is found
without specifying its form. Ordered by time, nothing stands out. The
test needs each observation’s score, which takes a pass of the filter
for every eight parameters, so it only runs when asked for.

## Everything at once

By default
[`ctScoreTest()`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
runs every family: the three above (covariates only when given) and a
fourth, which tests each cell the model fixes at zero as if it were
free, here the measurement error covariance and the cross loadings. On
the last data set:

``` r

ctScoreTest(fit3)
#> Score tests at the fitted model (laplace route, 150 subjects)
#> 
#> Unmodelled heterogeneity (LM: the variance with its covariances, against an
#> equal mixture of chi-square on df - 1 and df; z: the variance alone, one sided;
#> p: the stronger of the two; sd: the one-step estimate, raw scale):
#>        parameter    LM df      z     p p_adj     sd sd_natural
#>  drift_eta1_eta2 11.00  5  2.080 0.032 0.212 0.1990     0.1990
#>        diff_eta2  4.42  5  1.870 0.047 0.212 0.0859     0.0957
#>   diff_eta2_eta1 10.70  5  0.604 0.073 0.219 0.0913     0.0913
#>       drift_eta1  3.95  5  1.610 0.104 0.229 0.1240     0.1120
#>       drift_eta2  4.04  5  1.390 0.142 0.229 0.1230     0.0819
#>        diff_eta1  3.73  5  1.340 0.153 0.229 0.0515     0.0782
#>  drift_eta2_eta1  6.62  5 -0.144 0.333 0.426     NA         NA
#>           mvarY2  4.18  5  0.738 0.379 0.426 0.0390     0.0377
#>           mvarY1  2.91  5  0.273 0.602 0.602 0.0274     0.0269
#> Covariances with: T0m_eta1, T0m_eta2, mm_Y1, mm_Y2
#> 
#> A new dimension of variation: largest eigenvalue 6.6, p = 0.083; direction drift_eta2 1, drift_eta1_eta2 0.94, diff_eta2_eta1 0.72, diff_eta2 0.51, drift_eta2_eta1 0.48
#> 
#> State dependence of the drift:
#>  equation    LM df      p
#>      eta1 31.60  3 <5e-04
#>      eta2  4.39  3 0.2340
#>       all 39.50  6 <5e-04
#> 
#>       term equation     LM        p    p_adj estimate     se
#>  eta1*eta2     eta1 25.700   <5e-04 2.38e-06   0.2550 0.0504
#>  eta1*eta1     eta1 18.100   <5e-04 6.31e-05   0.1020 0.0240
#>  eta2*eta2     eta1  5.590 1.81e-02 3.62e-02   0.0874 0.0370
#>  eta2*eta2     eta2  3.180 7.46e-02 1.12e-01  -0.0346 0.0194
#>  eta1*eta1     eta2  0.302 5.83e-01 6.17e-01   0.0123 0.0225
#>  eta1*eta2     eta2  0.251 6.17e-01 6.17e-01       NA     NA
#>                                                             reads
#>  DRIFT[eta1,eta1] varies with eta2, or DRIFT[eta1,eta2] with eta1
#>                                 DRIFT[eta1,eta1] varies with eta1
#>                                 DRIFT[eta1,eta2] varies with eta2
#>                                 DRIFT[eta2,eta2] varies with eta2
#>                                 DRIFT[eta2,eta1] varies with eta1
#>  DRIFT[eta2,eta1] varies with eta2, or DRIFT[eta2,eta2] with eta1
#> 
#> State dependence of the innovation sds:
#>  equation    LM df      p
#>      eta1 30.80  2 <5e-04
#>      eta2  5.99  2 0.0580
#>       all 39.50  4 <5e-04
#> 
#>                      term     LM        p    p_adj estimate     se
#>  sd(eta1) x exp(v * eta2) 28.500   <5e-04 3.77e-07   0.1900 0.0356
#>  sd(eta1) x exp(v * eta1)  9.610 1.93e-03 3.87e-03   0.1120 0.0361
#>  sd(eta2) x exp(v * eta2)  5.900 1.52e-02 2.02e-02  -0.0973 0.0401
#>  sd(eta2) x exp(v * eta1)  0.722 3.95e-01 3.95e-01  -0.0315 0.0371
#> 
#> Fixed cells, freed: joint LM = 12.8, df = 1, p = 0.002
#>                cell   LM      p    p_adj estimate     se
#>  MANIFESTVAR[Y2,Y1] 12.8 <5e-04 0.000339    0.196 0.0548
#> 
#> Parameters depending on a subject's level (a random intercept):
#>  level    LM df     p
#>  mm_Y1  8.82  9 0.514
#>  mm_Y2  9.89  9 0.458
#>    all 24.00 18 0.371
#> 
#>                     term    LM      p p_adj estimate     se
#>  drift_eta1_eta2 ~ mm_Y2 2.730 0.0982 0.873   0.0480 0.0290
#>           mvarY2 ~ mm_Y2 2.020 0.1550 0.873   0.0153 0.0108
#>       drift_eta2 ~ mm_Y2 1.480 0.2230 0.873  -0.0181 0.0148
#>  drift_eta1_eta2 ~ mm_Y1 1.400 0.2370 0.873  -0.0343 0.0290
#>       drift_eta1 ~ mm_Y2 1.060 0.3030 0.873   0.0270 0.0262
#>           mvarY2 ~ mm_Y1 0.813 0.3670 0.873       NA     NA
#>   diff_eta2_eta1 ~ mm_Y2 0.744 0.3880 0.873       NA     NA
#>        diff_eta1 ~ mm_Y2 0.721 0.3960 0.873       NA     NA
#>  drift_eta2_eta1 ~ mm_Y1 0.527 0.4680 0.873       NA     NA
#>   diff_eta2_eta1 ~ mm_Y1 0.487 0.4850 0.873       NA     NA
```

The state dependence stands out, in the drift of `eta1` and in the size
of its innovations. So does the measurement error covariance between
`Y1` and `Y2`, although the data were generated without one, probably as
another effect of the missing moderation. The heterogeneity tests stay
quiet here, but need not. A drift that changes with the state can also
look like a drift that differs between subjects, and in simulations a
missing nonlinearity did affect them. A badly misspecified model can
show up in several tests at once, so free what stands out most, refit,
and test again.

## Do the model’s predictions look like the data?

The usual check of a dynamic model compares the data with data generated
from the fitted model: the autocorrelations and cross-correlations (ACF
and CCF), how much subjects differ from each other, and so on.
[`ctFitAdequacy()`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)
turns these comparisons into tests, grouped in families aimed at
different kinds of misfit.

Most families use the one-step-ahead prediction errors. Under the model
these are independent and standard normal however the observations are
spaced, so no generated data are needed. Each variable’s error is taken
given the other variables observed at the same time, so a result points
to a variable.

- `residuals`: correlations of the errors within subjects, over every
  pair of a subject’s observations, binned by the time between them.
  This is the residual ACF and CCF in continuous time. A slow process
  missing from the model shows as correlation over the time gaps it
  spans.
- `latent`: the same for each latent process, using the standardised
  change in the predicted state at each observation, correlated with the
  same subject’s changes at other times. With many indicators per
  process, each item’s errors can be nearly uncorrelated while the
  processes themselves are misspecified, for example with a cross effect
  missing. One statistic per process detects this more easily than one
  per item.
- `volatility`: the same for the squared errors, for variance that
  clusters in time.
- `pairs`: whether two variables’ errors at the same time correlate as
  the model says, for a wrong measurement error covariance or loading.
- `periodic`: a cycle the model lacks, such as a daily rhythm, from the
  periodogram of each variable’s errors and each process’s updates. It
  depends only on the time between observations, so it finds a rhythm
  whatever its phase in each subject. Give the periods that matter
  (`periods = 24` with time in hours), or by default it scans a range.
- `tvcovariates` (when `tvcovariates` are given): whether time varying
  variables the model does not use predict its errors, at the same and
  the previous observation.
- `distribution`: skewness, kurtosis and outlying occasions.
- `fit`: whether some subjects fit much worse than the model allows.
- `dependence`: whether an error or its size depends in any way on what
  came before: the predicted states, the previous errors, or the time
  since the last observation.
- `effects`: for random effects integrated with
  `intoverpop = "laplace"`, whether their normal distribution fits. At
  cuts along each axis, each subject’s posterior probability that its
  effect lies below the cut is averaged over subjects. Under the model
  this average matches the normal distribution, and a skewed or bimodal
  distribution of subjects moves it.
- `scores` (when asked for): whether each observation’s score is
  uncorrelated with the scores of the same subject’s other observations,
  which holds under the model. This checks for any kind of dynamic
  misspecification, but is weak because it is spread over many
  statistics. The targeted tests above detect what they cover far more
  often.

The error statistics allow for the parameters having been estimated from
the same data, by projecting each subject’s contribution on its scores,
as the score tests do. Three further families compare statistics of the
data with generated data: the spread over subjects of their own means,
standard deviations and lag one correlations, whether change depends
nonlinearly on level, and, when covariates are given, whether the
subjects’ statistics depend on them. With many observations every misfit
becomes significant, so the size of the remaining misfit is also
reported, as the root mean square residual correlation of the items and
of the latent processes, with intervals and a test of whether each is at
most `close` (0.05 by default).

On the fit from the second example, which lacks `age`, the checks fail
through the subjects’ spread and its relation to age:

``` r

ctFitAdequacy(fit2, covariates = covariates)
#> Adequacy: misfit detected by 2 of 12 families at alpha = 0.05
#> (estimation removed by: derivative; 100 generated datasets)
#> 
#>        family statistics                   worst worst_z        p    p_adj
#>     residuals         16         Y2~Y1 gap >10.2   1.620 8.49e-01 1.00e+00
#>        latent         16     eta2~eta1 gap >10.2   1.630 8.30e-01 1.00e+00
#>    volatility          8     Y2^2 gap 4.88-10.20   1.790 6.75e-01 1.00e+00
#>         pairs          1                   Y1:Y2  -0.027 9.63e-01 1.00e+00
#>      periodic         60 latent eta1 period 16.5   2.470 3.48e-01 1.00e+00
#>  distribution          8            kurtosis(Y1)  -1.160 9.17e-01 1.00e+00
#>           fit          2      spread of unit fit   2.410 3.20e-02 3.20e-01
#>       effects         20             P(u3 < 1.5)  -3.110 1.64e-01 1.00e+00
#>    dependence          4          dcor(Y2, past)   1.120 7.07e-01 1.00e+00
#>      subjects          8          sd(lag1_Y1_Y2)   5.510 2.57e-06 2.83e-05
#>     nonlinear         10           cor(dY1^2,Y2)   1.310 7.23e-01 1.00e+00
#>    covariates         16     cor(age,lag1_Y1_Y2)   7.850 9.94e-11 1.19e-09
#> 
#> Size of what remains (root mean square, bias-corrected; p_close tests that it is at least close):
#>                                     measure estimate lower   upper close
#>              residual correlation, all gaps        0     0 0.00355  0.05
#>         residual correlation, shortest gaps        0     0 0.00654  0.05
#>       latent residual correlation, all gaps        0     0 0.00277  0.05
#>  latent residual correlation, shortest gaps        0     0 0.00675  0.05
#>                      within-row correlation        0     0 0.00000  0.05
#>  p_close
#>        0
#>        0
#>        0
#>        0
#>        0
```

On the linear fit to the state dependent data:

``` r

ctFitAdequacy(fit3)
#> Adequacy: misfit detected by 1 of 11 families at alpha = 0.05, each shown within the close-fit bound (remaining correlations below 0.05)
#> (estimation removed by: derivative; 100 generated datasets)
#> 
#>        family statistics                    worst worst_z        p   p_adj
#>     residuals         16     Y2~Y2 gap 2.59- 4.93   -3.11 0.042300 0.42300
#>        latent         16 eta2~eta2 gap 2.59- 4.93   -2.70 0.158000 0.82000
#>    volatility          8          Y2^2 gap <=2.59    1.77 0.552000 1.00000
#>         pairs          1                    Y1:Y2    3.80 0.000415 0.00457
#>      periodic         60           Y2 period 2.52    2.41 0.204000 0.82000
#>  distribution          8                 skew(Y1)    1.89 0.779000 1.00000
#>           fit          2       spread of unit fit    1.31 0.228000 0.82000
#>       effects         20             P(u2 < -1.5)    3.09 0.079400 0.63600
#>    dependence          4         dcor(Y1^2, past)    2.09 0.137000 0.82000
#>      subjects          8           sd(lag1_Y2_Y1)    2.73 0.117000 0.81800
#>     nonlinear         10            cor(dY1^2,Y2)    2.77 0.059000 0.53100
#> 
#> Size of what remains (root mean square, bias-corrected; p_close tests that it is at least close):
#>                                     measure estimate lower   upper close
#>              residual correlation, all gaps  0.00666     0 0.01220  0.05
#>         residual correlation, shortest gaps  0.00332     0 0.01040  0.05
#>       latent residual correlation, all gaps  0.00347     0 0.00887  0.05
#>  latent residual correlation, shortest gaps  0.00000     0 0.00725  0.05
#>                      within-row correlation  0.01530     0 0.02230  0.05
#>  p_close
#>  0.0e+00
#>  0.0e+00
#>  0.0e+00
#>  0.0e+00
#>  1.2e-64
```

Only the `pairs` family fails, on the same measurement error covariance
the score tests flagged, and the remaining correlations are within the
close-fit bound. The state dependence that the score test found
(adjusted p of 2e-6) does not show. A moderated drift, averaged over the
data, still decays much like a single exponential, so the errors it
leaves are nearly uncorrelated and the root mean square residual
correlation is close to zero. The checks are useful for things no score
test covers: a missing slow process, variance that clusters, errors with
heavy tails, or a process of the wrong order. In simulations the
residual correlations over time gaps found a missing slow process (with
a time constant ten times the fast one’s) in 65 of 100 data sets, where
correlations at a few lags found it in 5. Passing the checks does not
show that the model is right, which is why the search in the next
section screens with the score tests and stops on the checks.

## Automating the search

[`ctModelImprove()`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
runs the tests in a loop. It screens the fit, refits with the most
strongly indicated changes one at a time (starting from the current
estimates), and keeps the best if the refit’s likelihood confirms it. It
repeats this until nothing more is indicated, and then checks the model
with
[`ctFitAdequacy()`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md).
If a check fails, it tries structural changes the score tests cannot
propose: a slow process driving a latent
([`ctModelSlowComponent()`](https://cdriveraus.github.io/ctsem/reference/ctModelSlowComponent.md)),
tried first for a latent whose updates stay correlated over long gaps; a
daily or other cycle at the period the `periodic` family points to
([`ctModelPeriodicComponent()`](https://cdriveraus.github.io/ctsem/reference/ctModelPeriodicComponent.md));
a separate process for an item whose errors stay correlated over time
([`ctModelUniqueProcess()`](https://cdriveraus.github.io/ctsem/reference/ctModelUniqueProcess.md));
a second order process
([`ctModelHigherOrder()`](https://cdriveraus.github.io/ctsem/reference/ctModelHigherOrder.md));
or models you supply. Finally it tries removing each change it kept, in
case a later change made it unnecessary.

Here the data have two of the problems above. `age` affects the
intercept of `Y1` and the cross effect, `eta2` moderates how `eta1`
returns to baseline, and `w` is unrelated.

``` r

gen4 <- generator(DRIFT = matrix(c("-0.6 + 0.3 * tanh(eta2)", "drift21||||age=0.3",
  0, -0.4), 2, 2),
  MANIFESTMEANS = c("mm1|param|TRUE|1|age=0.6", "mm2|param|TRUE|1"),
  TIpredNames = c("age", "w"))
d4 <- data.frame(ctGenerate(gen4, n = 150, Tpoints = 30, dtmean = 1,
  logdtsd = 0.5, backend = "julia", values = c(drift21 = 0.25, mm1 = 0, mm2 = 0)))
fit4 <- ctFit(d4[, c("id", "time", "Y1", "Y2")], model, backend = "julia", cores = 2)
```

``` r

imp <- ctModelImprove(fit4, covariates = unique(d4[, c("id", "age", "w")]))
imp
#> Model improvement: stopped because adequate
#> Log likelihood -11294.2 -> -11161.75 (23 -> 26 parameters)
#> 
#>  round                                                             change
#>      1                                        age as a predictor of mm_Y1
#>      2                                     eta1*eta2 in the drift of eta1
#>      3      drift_eta2_eta1 individually varying, with T0m_eta1, T0m_eta2
#>      4                         age also as a predictor of drift_eta2_eta1
#>  prune                                   drop age as a predictor of mm_Y1
#>  prune                                drop eta1*eta2 in the drift of eta1
#>  prune drop drift_eta2_eta1 individually varying, with T0m_eta1, T0m_eta2
#>  prune                    drop age also as a predictor of drift_eta2_eta1
#>  screen_p dloglik dpar  done                                         note
#>  5.50e-03   20.90    1  TRUE               LR p = 1.03e-10 against 0.0167
#>  2.17e-06   11.70    1  TRUE               LR p = 1.25e-06 against 0.0167
#>  5.00e-04   33.30    5  TRUE               LR p = 5.36e-13 against 0.0167
#>  5.00e-04   67.10    1  TRUE                 LR p = 4.93e-31 against 0.05
#>        NA  -20.30   -1 FALSE                 LR p = 1.78e-10 against 0.05
#>        NA      NA   -1 FALSE kept without a refit: Wald 21.4 against 7.68
#>        NA   -1.12   -5  TRUE                    LR p = 0.815 against 0.05
#>        NA  -99.30   -1 FALSE                 LR p = 4.42e-45 against 0.05
#> 
#> Adequacy: no family detects misfit at alpha = 0.05
#> (estimation removed by: plugin; 100 generated datasets)
#> 
#>        family statistics                worst worst_z      p p_adj
#>     residuals         16     Y2~Y2 gap <=2.54   2.960 0.0728 0.728
#>        latent         16 eta2~eta2 gap <=2.54   2.870 0.0635 0.699
#>    volatility          8       Y2^2 gap >10.1  -1.870 0.2340 1.000
#>         pairs          1                Y1:Y2  -0.991 0.2270 1.000
#>      periodic         60       Y2 period 2.29   2.390 0.8900 1.000
#>  distribution          9         kurtosis(Y2)  -0.988 0.8860 1.000
#>           fit          2   spread of unit fit  -2.500 0.0238 0.286
#>       effects         20        P(u4 < -0.75)   2.820 0.1240 1.000
#>    dependence          4     dcor(Y2^2, past)   1.910 0.1770 1.000
#>      subjects          8       sd(lag1_Y1_Y2)  -1.830 0.8080 1.000
#>     nonlinear         10        cor(dY2^2,Y1)   0.873 0.5650 1.000
#>    covariates          8    cor(w,lag1_Y2_Y1)   0.588 0.8690 1.000
#> 
#> Size of what remains (root mean square, bias-corrected; p_close tests that it is at least close):
#>                                     measure estimate lower  upper close
#>              residual correlation, all gaps  0.00524     0 0.0104  0.05
#>         residual correlation, shortest gaps  0.01030     0 0.0171  0.05
#>       latent residual correlation, all gaps  0.00691     0 0.0127  0.05
#>  latent residual correlation, shortest gaps  0.01060     0 0.0181  0.05
#>                      within-row correlation  0.00000     0 0.0144  0.05
#>    p_close
#>   0.00e+00
#>  1.31e-139
#>   0.00e+00
#>  1.87e-104
#>  1.26e-122
```

Each round shows what was tried, the score test that proposed it, and
the result of the refit. The first round added `age` as a predictor of
the intercept of `Y1`, and the second the moderation. The third let the
cross effect vary between subjects, together with its covariances with
the initial states. Once `age` is in the model, its effects on the other
parameters are tested too, and the fourth round added its effect on the
cross effect. When the search then tried removing each kept change, the
random effect on the cross effect was no longer needed, because age
explains how subjects differ in it, and it was dropped. The other
changes were kept, and every adequacy check passed. The final model
includes each effect that generated the data, with age effects of 0.58
on the intercept of `Y1` and 0.32 on the cross effect, against 0.6 and
0.3 in the generating model. A search that decides at `alpha = 0.05`
will sometimes keep effects that are not there, and `alpha = 0.01` is
stricter.

The final fit is an ordinary fit:

``` r

s <- summary(imp$fit)
s$tipreds
#>                          mean    sd  2.5%   50% 97.5%  mode      z
#> tip_age_drift_eta2_eta1 0.324 0.027 0.269 0.324 0.375 0.323 12.026
#> tip_age_mm_Y1           0.577 0.083 0.416 0.572 0.741 0.578  6.986
s$popmeans[grepl("^sdep", rownames(s$popmeans)), ]
#>                      mean    sd  2.5%   50% 97.5%  mode
#> sdep_eta1_eta1_eta2 0.131 0.028 0.077 0.131 0.183 0.131
```

The model was chosen using the data it is fitted to, so its standard
errors and tests are too optimistic about the structure the search
added. Treat that structure as a hypothesis to confirm on other data.
Within one data set, `holdout = 0.3` sets aside 30 per cent of the
subjects. The search runs on the rest, and each change it keeps is then
tested on the held-out subjects alone. With a large data set,
`close = 0.05` also stops the search once the remaining residual
correlations are not significantly larger than 0.05, since otherwise
every small misfit would be significant.

## How the tests work, briefly

- The scores are those of the objective the model was fitted with, so
  they integrate over random effects in the same way as the fit. A
  covariate’s effect adds to a parameter’s raw value, so its score is
  the covariate times each subject’s score for that parameter. Any
  number of covariates therefore costs a few matrix products after one
  set of per-subject scores and Hessians.
- At zero, a random effect’s variance has the score
  $`\frac{1}{2}\sum_i (s_{ij}^2 + H_{i,jj})`$, and its covariance with
  an existing effect has the score
  $`\sum_i (s_{ij} s_{ir} + H_{i,jr})`$. This is the information matrix
  test in those parameters, which is the score test for neglected
  heterogeneity. Every subject already has one such effect, their
  initial state, so the covariances are tested even in a model with no
  random effects.
- State dependence is tested on the fitted model extended with the drift
  terms at zero, which gives exactly the same fit.
- The fitted model’s curvature is used to project the fitted parameters
  out of each test and to size the one-step estimates. Each test’s
  variance comes from the spread of the subjects’ contributions, and its
  p-value from resampling subjects. This keeps a test valid when the
  model is wrong in some other way, which is when a screen is used. A
  variance based on the model’s curvature would make unrelated tests
  reject too.

The tests need subjects to be independent units, so a model with levels
above the subject can only be tested for state dependence. See
[`?ctScoreTest`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
for details.
