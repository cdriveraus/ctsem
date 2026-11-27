# Drift diffusion models

``` r

library(ctsem)
library(ggplot2)
```

The drift diffusion model describes a choice between two responses and
the time it takes. Evidence accumulates at a mean rate, the drift, until
it reaches one of two boundaries. The boundaries are `boundary` apart,
the accumulation starts at `bias` times that distance from the lower
boundary, and a nondecision time `ndt` is added for encoding and
responding.

In ctsem the drift is a latent process. It can be constant within a
person, change over a session, or depend on other processes in the
model, and each trial’s likelihood is exact. Boundary, nondecision time
and bias can differ between people and change with trial level
predictors.

A choice and response time variable holds the response time in seconds,
made negative for a response at the lower boundary. For accuracy coding,
where the upper boundary is the correct response:

``` r

dat$rt <- ifelse(dat$correct == 1, dat$rt, -dat$rt)
```

## A drift diffusion model

We simulate 20 people with 200 trials each. Half the trials are hard,
which lowers the drift, and in alternate blocks of 50 trials people are
told to respond quickly, which lowers the boundary. People differ in
their drift, boundary and nondecision time. The time column is the trial
number.

``` r

set.seed(1)
design <- do.call(rbind, lapply(1:20, function(i) data.frame(id = i,
  time = 1:200, hard = rbinom(200, 1, 0.5), speed = rep(c(0, 1), each = 50,
    times = 2))))

model <- ctModelDDM(name = 'rt', drift = ~ hard, boundary = ~ speed)
truth <- list(popmeans = c(drift_mean = 1.5, drift_hard = -0.8, boundary = 1.4,
  boundary_speed = -0.4, ndt = 0.3),
  popsd = c(drift_mean = 0.4, boundary = 0.2, ndt = 0.04))
dat <- data.frame(ctGenerate(model, datastruct = design, backend = 'julia',
  values = truth))
head(dat)
#>   id time hard speed         rt
#> 1  1    1    0     0  0.6209449
#> 2  1    2    0     0  0.6788036
#> 3  1    3    1     0  0.5821203
#> 4  1    4    1     0  0.4483293
#> 5  1    5    0     0  1.0161920
#> 6  1    6    1     0 -0.6448534
```

[`ctModelDDM()`](https://cdriveraus.github.io/ctsem/reference/ctModelDDM.md)
builds the model from the predictors of each parameter. The drift’s
predictors shift the drift on their own trial. The boundary’s and the
nondecision time’s act on the log scale, so `boundary_speed = -0.4`
multiplies the boundary by `exp(-0.4)` under speed instructions. By
default the drift, boundary and nondecision time vary over people. The
model is an ordinary ctsem model:

``` r

model
#> ctsem model (continuous time)
#>   1 latent (drift), 1 manifest (rt)
#>   2 time dependent predictors: hard, speed
#>   non-Gaussian indicators: rt (choice and response time, drift diffusion)
#> 
#> 7 free parameters
#> 
#> 3 individually varying: drift_mean, boundary, ndt
#>   This is a multilevel model. It estimates a 3 x 3 population covariance
#>   over those parameters -- 3 standard deviations and 3 correlations on top of the 7 above.
#>   Set indvarying = FALSE in a cell (see ?ctPar) for any that should not vary.
#> 
#> T0MEANS
#>                 
#> drift drift_mean
#> 
#> T0VAR
#>       drift
#> drift 0.001
#> 
#> DRIFT
#>       drift
#> drift     0
#> 
#> DIFFUSION
#>       drift
#> drift     0
#> 
#> CINT
#>        
#> drift 0
#> 
#> LAMBDA
#>    drift
#> rt     1
#> 
#> MANIFESTMEANS
#>                     
#> rt drift_hard * hard
#> 
#> MANIFESTVAR
#>    rt
#> rt  0
#> 
#> TDPREDEFFECT
#>       hard speed
#> drift    0     0
#> 
#> PARS
#>                [,1]
#> [1,]       boundary
#> [2,]     drift_hard
#> [3,] boundary_speed
```

``` r

fit <- ctFit(dat, model, cores = 4)
s <- summary(fit)
data.frame(estimate = s$popmeans[names(truth$popmeans), 'mode'],
  sd = s$popmeans[names(truth$popmeans), 'sd'], generated = truth$popmeans)
#>                estimate    sd generated
#> drift_mean        1.593 0.102       1.5
#> drift_hard       -0.817 0.061      -0.8
#> boundary          1.325 0.044       1.4
#> boundary_speed   -0.407 0.013      -0.4
#> ndt               0.298 0.010       0.3
data.frame(estimate = s$popsd[names(truth$popsd), 'mode'],
  generated = truth$popsd)
#>            estimate generated
#> drift_mean    0.409      0.40
#> boundary      0.192      0.20
#> ndt           0.048      0.04
```

The estimates are close to the values the data were generated from. Each
person’s own values come from `ctSubjectPars(fit)`.

To let the drift change over the session, use `dynamic = TRUE`. The
drift then follows a mean reverting process around each person’s level,
with an auto effect `drift_auto` and a diffusion `drift_diffusion` on
the scale of the time column. `sv = TRUE` adds trial to trial
variability of the drift, and `bias = ~ 1` estimates the start point.

## A drift diffusion model related to another process

Here a person’s alertness drifts over a ten minute session, and the
drift of their decisions rises and falls with it. Alertness is also
rated now and then. Both measures load on the same latent process. The
rating loads on it directly, and the drift through a loading
`drift_alert`. A model with more than one kind of indicator is written
with
[`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md),
which takes the choice and response time indicator as
`manifesttype = 5`, with its boundary, nondecision time and bias in
`WIENER`.

``` r

model2 <- ctModel(type = 'ct', manifestNames = c('rating', 'rt'),
  manifesttype = c(0, 5), latentNames = 'alert',
  LAMBDA = matrix(c(1, 'drift_alert'), 2, 1),
  MANIFESTMEANS = matrix(c('rating_mean||FALSE', 'drift_mean||TRUE'), 2, 1),
  MANIFESTVAR = matrix(c('rating_sd', 0, 0, 0), 2, 2),
  DRIFT = matrix('alert_auto'), DIFFUSION = matrix('alert_diffusion'),
  CINT = matrix(0), T0MEANS = matrix('stationary'), T0VAR = matrix('stationary'),
  WIENER = c('boundary||TRUE', 'ndt||TRUE', 0.5))
```

The drift of the rt indicator is its linear predictor,
`drift_mean + drift_alert * alert`. Alertness has mean zero, and its
scale is set by the rating’s loading of one. The time column is in
minutes, with a trial every two to four seconds and a rating after every
20th trial.

``` r

set.seed(2)
design2 <- do.call(rbind, lapply(1:20, function(i) {
  time <- cumsum(runif(200, 2, 4)) / 60
  data.frame(id = i, time = time, rating = ifelse(1:200 %% 20 == 0, 0, NA),
    rt = 0)
}))
truth2 <- list(popmeans = c(alert_auto = -0.4, alert_diffusion = 0.6,
  rating_mean = 0, rating_sd = 0.3, drift_mean = 1.2, drift_alert = 0.8,
  boundary = 1.3, ndt = 0.3),
  popsd = c(drift_mean = 0.3, boundary = 0.2, ndt = 0.04))
gen2 <- ctGenerate(model2, datastruct = design2, backend = 'julia',
  values = truth2)
dat2 <- data.frame(gen2)
```

``` r

fit2 <- ctFit(dat2, model2, cores = 4)
s2 <- summary(fit2)
data.frame(estimate = s2$popmeans[names(truth2$popmeans), 'mode'],
  sd = s2$popmeans[names(truth2$popmeans), 'sd'], generated = truth2$popmeans)
#>                 estimate    sd generated
#> alert_auto        -0.276 0.077      -0.4
#> alert_diffusion    0.495 0.064       0.6
#> rating_mean        0.111 0.104       0.0
#> rating_sd          0.334 0.042       0.3
#> drift_mean         1.235 0.113       1.2
#> drift_alert        0.822 0.086       0.8
#> boundary           1.336 0.051       1.3
#> ndt                0.311 0.012       0.3
```

The loading `drift_alert` says how much the drift changes per unit of
alertness. Because the drift is filtered together with the ratings,
every trial also informs the estimate of alertness, between ratings as
well as at them. The smoothed alertness of the first person, against the
alertness the data were generated from and the ratings:

``` r

k <- ctKalmanArray(fit2, timestep = 'asdata', subjects = 1)
states <- attr(gen2, 'states')
one <- data.frame(time = k$time[, 1], smooth = k$etasmooth[1, , 1],
  sd = sqrt(k$etasmoothcov[1, , 1, 1]),
  generated = states$alert[states$id == 1], rating = dat2$rating[dat2$id == 1])
ggplot(one, aes(time)) +
  geom_ribbon(aes(ymin = smooth - 1.96 * sd, ymax = smooth + 1.96 * sd),
    fill = 'steelblue', alpha = 0.25) +
  geom_line(aes(y = smooth), colour = 'steelblue') +
  geom_line(aes(y = generated), colour = 'grey30', linetype = 'dashed') +
  geom_point(aes(y = rating), colour = 'firebrick') +
  theme_bw() + labs(x = 'minutes', y = 'alertness')
```

![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAABUAAAAKgCAMAAABtOycFAAABQVBMVEUAAAAAACsAAFUAKysAK1UAK4AAVYAAVaorAAArKwArKysrK1UrK4ArVVUrVYArVaorgIArgNQzMzNGgrRNTU1NTWtNTYhNa4hNa6ZNiMRVAABVACtVKwBVKytVK4BVVStVVVVVVYBVgFVVgIBVgKpVgNRVqtRVqv9rTU1ra4hriKZriMRrpuGAKwCAKyuAVSuAVVWAVYCAqtSA1P+ITU2Ia02Ia2uIxOGIxP+ma02miGumxOGm4f+qVQCqVSuqVVWqgCuqgFWqqqqqqtSq1NSq1P+q//+yIiLC0d3EiE3EiGvExKbE4f/E///R4OzUgCvUgFXUqlXUqoDUqqrUqtTU1KrU1NTU1P/U///hpmvhxIjhxKbh///r6+v/qlX/xIj/1ID/1Kr/1NT/4ab/4cT/4eH//6r//8T//9T//+H///999sIfAAAACXBIWXMAAB2HAAAdhwGP5fFlAAAgAElEQVR4nO2dC7ssR3WeR4CiAzEQ0GBAxgHbBBSwDTEiAWwRsA37yMEjQxIjGR9dQJGO9v//AZn79KVqVa2qtaq+7vne5wGds/dMz+pvrXpPz0xPz+aREEJIEZveBRBCyFKhQAkhpBAKlBBCCqFACSGkEAqUEEIKoUAJIaQQCpQQQgqhQAkhpBAKlBBCCqFACSGkEAqUEEIKoUAJIaQQCpQQQgpxEegXCSFklbQQqP4uz57Zl1EDVjlo6TAeGax6GI+ILh4KNA+sctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSYTwyWPUwHhEK1AOsctDSmcez3fYo4wx+PF1hPCIUqAdY5aClE4inp0EXEE9PGI8IBeoBVjlo6VCgMlj1MB4RCtQDrHLQ0qFAZbDqYTwiFKgHWOWgpUOBymDVw3hEKFAPsMpBS4cClcGqh/GIUKAeYJWDlg4FKoNVD+MRoUA9wCoHLZ1RPCd1UqADsOphPCIUqAdY5aClM4yn6xmgJ5DjAYDxiFCgHmCVg5bOIB4AfyLHgwDjEaFAPcAqBy0dClQGqx7GI0KBeoBVDlo6FKgMVj2MR4QC9QCrHLR0KFAZrHoYjwgF6gFWOf3TmWiSAhXBqofxiFCgHmCV0z2d7cSTvesZ0z2eKVj1MB4RCtQDrHJ6pzP1J+ORwaqH8YhQoB5gldM5nZk/+UkkGax6GI8IBeoBVjl905n7MxRPR4OiDQ+nRwarHgrUA6xyegt09iMKVASrHsYjQoF6gFUOWjoUqAxWPYxHhAL1AKsctHQoUBmsehiPCAXqAVY5AOmM/UiBimDVw3hEKFAPsMpBSGckyGfzn1KgN7DqYTwiFKgHWOUgpEOB5oNVD+MRoUA9wCoHIZ2UQDuCEM8IrHoYjwgF6gFWOR3TuRoyV6APe7yLmoA2PJweGax6KFAPsMpZkEAfHtobFG14OD0yWPVQoB5glbMcgT48dDAo2vBwemSw6qFAPcAqZzECfXjoYVC04eH0yGDVQ4F6gFUOgkBHUKAiWPUwHhEK1AOscnAFOoECPYJVD+MRoUA9wCqHAhVBGx5OjwxWPRSoB1jlLEGgpxvyTaQDWPUwHhEK1AOsctDSCcYzNGjjapYQTz8YjwgF6gFWOWjpCALlifSPnJ4EWPVQoB5glYOWjiTQDiwinn4wHhEK1AOsctDSoUBlsOphPCIUqAdY5UCkMxQkBSqCVQ/jEaFAPcAqByOdgSGfzX5CgQ7AqofxiFCgHmCVg5EOBZoLVj2MR4QC9QCrHIx0EgLtB0Y8A7DqYTwiFKgHWOX0SmcsSAo0F6x6GI8IBeoBVjkUqAja8HB6ZLDqoUA9wCqHAhVBGx5OjwxWPRSoB1jlYAh08FcKVASrHsYjQoF6gFUOiEBvUKAiWPUwHhEK1AOsclAFCgLa8DAeGax6KFAPsMqhQEXQhofxyGDVQ4F6gFXOIgTKE+mvYNXDeEQoUA+wykFLJxxPN4MuI55uMB4RRIE+I3fIdtu7AkKc4RFoHljloKXDI1AZrHoYjwjiEaj+LmhNxioHLR0KVAarHsYjQoF6gFUOTjpnSVKgIlj1MB4RCtQDrHJw0hkJNPpBpcbgxHMGqx7GI0KBeoBVDk46FGgOWPUwHhEK1AOscnDS2Z4sGRRoN3DiOYNVD+MRoUA9wCqnTzpBP1KgGWDVw3hEKFAPsMqBEujhIJQCFcGqh/GIUKAeYJUDJNDTk3gKVASrHsYjQoF6gFUOkkCPUKAiWPUwHhEK1AOscihQEbTh4fTIYNVDgXqAVQ6oQFFAGx7GI4NVDwXqAVY5FKgI2vAwHhmseihQD7DKWYhAez2lRxseTo8MVj0UqAdY5aClQ4HKYNXDeEQoUA+wykFLhwKVwaqH8YhQoB5glYOWDgUqg1UP4xGhQD3AKgctHQpUBqsexiNCgXqAVQ5aOhSoDFY9jEeEAvUAqxy0dChQGax6GI8IBeoBVjlo6VCgMlj1MB4RCtQDrHLQ0mE8Mlj1MB4RCtQDrHJ6pCMeTjIeEax6GI8IBeoBVjkUqAja8DAeGax6KFAPsMqhQEXQhofxyGDVQ4F6gFUOBSqCNjyMRwarHgrUA6xyKFARtOFhPDJY9VCgHmCVQ4GKoA0P45HBqocC9QCrHApUBG14GI8MVj0UqAdY5SxIoD1OpkcbHk6PDFY9FKgHWOWgpUOBymDVw3hEKFAPsMpBS4cClcGqh/GIUKAeYJWDlk6GQFuKdEHx9IDxiFCgHmCVg5aOKNDt6T/NallUPD1gPCIUqAdY5aClI8VzFmirSg4sKZ4OMB4RCtQDrHLQ0kkKtO0roUuKpwOMR4QC9QCrHLR0UgJt/E7SkuLpAOMRoUA9wCoHLR1RoI+Nn8AvK54OMB4RCtQDrHLap5M4hrz3eBJg1cN4RChQD7DKWaJAGz6NRxueu5+eBFj1UKAeYJWzRIE2NCja8Nz99CTAqocC9QCrHApUBG147n56EmDVQ4F6gFUOBSqCNjx3Pz0JsOqhQD3AKocCFUEbnrufngRY9VCgHmCV01eguxkUqAxWPYxHhAL1IM8QrRTRTaBzd554usveRAPQhmfRhmgAVj0UqAeZAm3kiC7pxOR5EuguqVAKFATGI0KBepBXTiuFNk9HkudFoGmFtgJteBZtiAZg1UOBepBZzraNQk/ptHis4wOk9HkRKIpC0YZn0YZoAFY9FKgH2eU0OQQ9p9NA1/vtp/V5EyiGQdGGZ9GGaABWPRSoB1jlXNNxN+h2m+PPm0AhDIo2PKjTgwJWPRSoB1jlDNJBEyiCQdGGB3d6MMCqhwL1AKucZgLd6gUKIFG04cGdHgyw6qFAPUiV0+2Kwa4PvMv051SgnR2KNjyLNkQDsOqhQD1ICnSF31mx1+BJoOnj0LlAeyoUbXgWbYgGYNVDgXpwfwLd7SoFOjPoij+olQCrHsYjQoF6kChn9ma4ryv80zk78GzOMoFODUqBYsB4RChQD5ICTf3AlH06rg8wNWGhQCcGpUAxYDwiFKgHWoH62sJZoHMTFgp07FAKFAPGI0KBeiCXEzid3VmgnifQJ2SpE+jAoRQoBoxHhAL1ICXQvJ9Z4fkBpAJ9JgR6USgFigHjEaFAPdCX46iLiT9NH6nInymB7uzrlEAbnkUbogFY9VCgHiAJdK/PZ+O/m226TJ9pge6M65RBG55FG6IBWPVQoB4ACfTgTyeBlvozLdCmp9WjDc+iDdEArHooUA+Aytn70kegU+dlfowzT6AtDYo2PEjT88h4ElCgHmCVMxWohUHnysv9IHyeQBsaFG14oKcHAKx6KFAPsMqZpGMh0IDxjAXazqBow4M9Pf3BqocC9QCrnJRA1ZeqDwrPWqDNDIo2PNjT0x+seihQD7DKyRKoQqERf1KgRmDVw3hEKFAPxHLaXorpcS7QyX/Pf8utKqI7hT8zBdrKoGjDs2hDNACrHgrUA2iBXsvYjv+at7GY7RwE2sigaMOzaEM0AKseCtSDBQh0+py90p8OAt1vsX7fM0AbHnBD9P3X7hE9HhkKNI8CgXpo9bzJiEBLtpgvyXqB7hr9S4M2PMCGkPvVqqBGj5MHBepByRGogy7Omwyl09Of+QJtsibRhgfUEDkdi21D+p2+ILMtWUCBegAv0LKz6cuVOSZfoC0MijY8mIbIbG1oC1k3UhRUd3djKFAPMAS6jQq0rz8pUBmsek7xZPd2evfc2ykKKr+rAxSoBygCPVcTEGjBEGevoSQKgTYwKNrwIBpC0dzhffNvqSmofp8MoUA9KCnHXKDbgEBrplmxiBJoBOpvULThwTOEsr3Xu2puqyjIbt8MoEA9ABHo+Q/XdJLTLNWgXEUSFKgIVj3Pnqn7e7qj7taKgqz3sAoK1INTOTolugs0S1nRralXkUDuifQnjGOZgTY8aIbQdevSMr8Wg8XTVaAfvfbVwIPot4O2BkoEas5IoHmDfD0bfrax7HvmoFySzjmhDQ+YIZ4WCFSPoiCsePoK9I0tBepP/pOw0RVBhpvIuKditWiXpHc8vttXA1XPro1AFT2GiqerQD9+Y7tygWIYVLEEtuOLKl22kHNHxWKhQEWA6jk0q41A85sMFM9jV4H+7ofbdQsU4hBUtwS22+ll6c7bSN9RsVbUS9I1IbThgTHEuVmNBJrdZJR4TvQT6Nvb7Xf/be0C7W3Q4iWgEqL29vp6PDNCGx4QQ1yb1UqguV3GiOdCR4F++W8f3123QMOHoO2sWrEEKNB+INQzaFY7gea1GSGeG33fRBoJ9IsXnq2B7fH/tqHfhH7owNMatlvf22u32iYyUjc19fTefWco0GyOmgyqsoFA6wd5qzSit0DXvrIw8Oihkt4R+OIq0KtI9dtBexb2THqiHv2N6Te2D55SFz0J0z6H15BdT+ysKlPQhqfXc9RYsxo+hT+QrBOrXUBP4a8Pot8O2hp4Jr1/VPIbBec5HLinUKB+Bi0SqJtB0YanjyHizWos0GSfsdpFgXrwTHKhq0CvY1grUM9DUIVAGxgUbXh6GEJqVmuBpvqM1S4K1INn0glMjgIdTOEqBNrkEBRteJobItGs5gJNNBqrXRSoB88kFfoJdDCD23qB1i6DOIUCdTIo2vC0NkSqWe0FKjcaq10UqC2nQ8+BQLcnMu9c9djDEdzWvonkSalAfQyKMzxn2taTbFaP6ZEKxmoXBWrLSZY9BDqewLUIdPL3inyi4AzPmab1pJvVZXqEirHaRYHacj4Cnf00887lDzwZwJUIdEp5PnFwhudMy3oymtVneuIlY7WLArXl5MCZQL0fdjZ/bQWqfMOpoh6H7HCG50zDenKa1emfX4R4MuAV6W0JC9SX1CiuSaAOBsUZnjPt6slqVq/nL/3jyYECtaW1QHMmkQIVwRmeM83qyWtWtxeAeseTBQXqQaty8gZxVQK1Nyja8Nzd9ETpG08eFKgH+nKkiRHvkkH5EsgVY0uBmhsUbXhaGSK3WR3fguwZTyYUqAdSOZe32jUzEyJ/DClQEbThaWSI7Gb1PIejXzy5UKAeJARaMDUzFFO4bIHOt23Toytow9PGEPnN6noSXK94sqFAPZAFWjY3YzRDuGyBuhsUbXiaGELRrL5nEfeJJx8K1INQOZen7qJsbr+UHyBx7/Fj1AjU5RPxtQK1NSja8LQwhKZZnT+G0SMeBRSoB4JA5aO1zCsIx+4d2kzVEgAVqKlB0YangSFUzer9Obb28WigQD0IC3R7nNxsgcY1kbwzBappFtjw+BtC16zeAp02G6tdFGgZD3tmP7x+YjNYzlkECoFGRBG983UjqxeopUHvTqDKZnUX6KTZWO2iQIt4eAgZ9HoxkAqBpg0av+9tI2YCLb6nhIFADQ16bwLVNqu/QMfNxmoXBVrCw0PQoCmBZugo/ZZz/K6xjSAsgREWArUzKAUqNwtgelrGo4QCLeDhIWjQrSjQvOfDyRf8VPc8gbAERtjUY9XNOxOoulkI09MuHi0UaAFRgV7rmd9nl/d8OPV8Vbhn9DcQS2AIBSriWo++WRDTk47H6kvBlVCgBRQINHdQYhrUbWWM9xJQv9NkVI9RN+9KoPqUMQQ6aHYsnj4GpUALUAvUaIIK77lWgRoZ9J4EWhAyBSpCgZYQfhPJV6DlrFagNga9I4GWZAwi0FuvKdDUg+jvAnIa05VJOZ1HjwIVoUDlZoEI9NpsCjT1IPq7gJxIf2VcTuHAmJ2DWbcE0nrsJ1ATg96PQIsipkBFKFAPRuUUzovdhTxWLFALg6INDwUaIREPBXp5EP1d0NbAsJzCaTH8EKWzQPWm19YjPIBFs8CGx02gytQvzYIR6E6OhwK9PIj+LmhrYFBO6bBopJQwmLtAtZtU1+NqULThoUBjSPF0Og2UAnXhVo5+Sk6yyLfSdjv97PuUWoFGNl5+iGwp0HqDog2PpUCHWiltFo5Ad0I8FOjtQfR3QVsD13IKhuQs0Oxbp/RZvQSiAi02qKlAqw2KNjymR6A3rxQ3iwKVoEA9uJRTMiQqL6XtuaNAU80CGx4XgZa2Ckugu3g8FOjtQfR3QVgDww6eyymaEftLcPoItOIsAVuB1hoUYXhGmAu0tE9HoAS6i8ZDgd4eRH8XhDUwF2jZiNh/D1GtQCM/BnkNdFdrUIThGWEs0OI2naBARShQK6YCLZ4RNIGGqfG8eT1VjUMYnhG2Aq3NFkugOwo0/SD6uyCsgbFAK0ZkOPMmMnUSaPl9KVARw3oMBqiPQON1U6DJB9HfBWENDDtoduKlyfP59Qt0V+NQhOEZYVaPSbDQAh0uOp5If3kQ/V0Q1sCgf3WGGI0OBZpNaeMQhmeEUT1GsaIJdBcTaB+DUqBWDM+4s5u58CDp5AX2KpaXQEsVijA8I0zq2VldjQZOoE8He0mBhh5EfxeENXBt366BQFVr424EWmZQhOEZYVDPMQ2bdyMpUBEK1JhTjw0HKTRJlx/mLZH7EWiRQYGG50R1Pecw1irQQZe33V8EpUBtufTYd5JaC9T4zCrPJalvGc7wnKms5xrFHQi0/7tIFKgp1x67TtK2p0At3tWqLEFE3TOY4blQV88gsQUL9Eh4B57uRoeg4T83gwK1ZNBjsyGSBquDQCuX5ct7hHiOv06VkEDbNJThuVJTT1FiIogCPfd4b0wKNPAg+ruArIFRj10H6zxaXQRasoGXA0R/nyohia5rIMNzQ6hnm5DEPLFqjUIKdHdO43F7+UMvKFA7Jj12HSzN6rAQ6Pb2R+V959YM2XTw+7BBVY+rahvG8AwoFmgwsVqDYgv0cfyH9lCgZkx77DpYbQU6fNFAsR5jB5VPJ7+d/D5oUN1LB5q+QQzPkHg9W1Ggs8SG/ymmYHoO37hY96Anwj0/1XOK45qLZf4qKFAr5j32pJ9AM+8hPSFP1VN1CHq6naJxCMMzQhRo9HdyHuXop+f0nd91j3pCEOjukQKNPYj+Lv3XwLjVrQSaR3OBii9m7jIEGj4EzXjkwQ2zW9d/eCZIAo39JhlHMerpeXgwMygFWgIFWjpaESwEqnhcWZ7HgtKbCJSgFGi2QvsPz4RJPQM7REWRE0chwvQEN/3wkGnQjMIkge4o0NiD6O/Sew1MW706geaTtOcuJ57AJhTfsnf9Y1bzeg/PjHE9GXLIjaOI+PSEX6LMFajwbTS35zvBgi5/oEAjD6K/S+c1sPMVaPhc+vz7NxRojj6z4klvJMoomozuoQs0eXtFHAXEpicmQEuBhgu6/umWEQU6fBD9XfqugXnP/QWqoZlA8/SZJ9Byg47jSrcPW6BnN+jfPLIiMj1R/2U/hY/Pda5AEU5mokCrCfT8PgWaqc+8eKwEmlYotEAvUojJoTilbMLTI5xVlv0mUr1AxyczdTEoBVpLqOfWAh1fY1l7/zYCzdZnpkArnsRPSDQQWaADOQRvGt1nu6u/BKdHPCu3+k14CrSKZQk02HNjY00Eql0cdQLNWwy5z95PBeXcyE6gCYMuRKCh9+Lju2x3/ayIQKW71J7ElC/QXSyjZlCgdWT0uJ7xv/fqk9qrBJp1OKHSZ65AiwwajkRsIbJAB8wFKubgK1BfFALdDT9dQIGeHkR/l15rIK/H9fQT6O0FrfijKfWZG0+RQSNVSk1cikCnBi3JoQB4gfZ9F4kCLSe3x/WMrsfZUqCDt1Rjb7uq9ZkdT8m5oLFMhDYuRKDTQ9BEDM0Eqn8k1R3mm5/UQ4HOHkR/l+NeVH49uBJxRM0FOvyz+sJIZgKd/z79qaNgQXk3m284vcfRW8Q7uRSB3jiIojAGPdYCFc7/jNx8WtDlN6f/rFmgz//33/+DtqAagbY0qDgHns96bhfZ6S7QMn3mx6M/BBWaEu3kIgWazCH2FzXGAlX6M7D9p5NfZHzUwBEngb7/3956fPzwa5vN5lM/0JZUIdB2Cj0OQnQSGghUMYlWAh3/qtSeO4VA1Yeg0u/F4UEiQ6Cp5OwFGp22Wh+q7zAVaNMnnjN8BPrm5hO/fHx8fXPg8CcVVQJtk+apfV2PQNsINHpWdIU+FfGoD0GPkShfB6VARU7TIwyb8hm59q73KNB3jtp878nmxbfef2XzBWVJdQJtEmeq9y3euGwj0OBpTFX23GkEqj0EPf1aZ9DFCTQ3vVsoNQY9Ts9wE9ONabY+uW2eQCc3ejr7eZOuRHAR6Ot7cx40+sIPDv9/+LPqQXQ3PzAUaINAU3PRSKC5t7Q9kT5xrc+sgrJvqT0ElQUaHgxggQZf3csOb5JKIYfpGUlsujWVQAvuGhXo9SfNOhPAQ6DP/+pgzrNGP3xF+xy+XqDOmWb32JVWAh1Tb8+dSqAFb8SvSaABg2Znp8ssyn56Ih/kSG49+bhGAu1pUA+Bnpz54Sublx67CdQx1Pweu9JBoCb23KnisTybfhceimUJtCCOuufwR4FmbGz63YCn2xYUNnuJYFrQfNstuzPBT6DvPdl8+7GjQL1SDQ9CsMeWzEatuUCt9KmLx1agoaFYlEAL0pADSfL0aeKNn/DXq56GpeSsidR9AgLtaFC/p/BvHl8C7fIaqGeqoYbHXue2pGINmAjUzJ47pUALHlR3MhOuQAPfIacP4xhI2d2O7AU63tbor9EvVj2zldtXIdARl7zaNemC05tILx3efj+Ys/278PNcLQkOQfRMC0P6CtRSn8p47lqg01/ps6hmMj2jYU89KdkLdCvewkqg/QzqdRrTgW8/Pv/e5nQcqnoQ5e0f4wI1V6jcWrnHdYgn7svUC9RUn1qBWj7yLjARCxJo5i7afY5zFxFo8MAz+OG4k0FjAxQ4OSo55ncg0P3T9z0vHd9IeuHb2pJMBWpr0OGsqHtcSfk7AbUCNdanNh5jgc4mYjkCve2DPAuGn4QPTU/kabvw6eLTk3mxkdd7Jot/GrxBMLAGeH2U8/t/9jf7/3z4J5/7hbokW4EaKjTU7kiPEzNQQvmqqBOo7bP3I0qBZn/HcV5CweFBIibQdAA5v9QynR7hcFK6PMM298l8hkBDtwgF1oJVXY0pglFdwXZHepyYgbbUCNRBn+p4ggItefXsTGh4kIjUk7+rNgI9ZzyaHvFFT0mgh21Jd1ZcpDEs0Mm3ezTjHgSqU2ikCdN+yz1OzEBbygXqok+9QHMNmiuOwPAgUStQo4uBnl91H0xP1IDXD4CNPvIZOPteOno9/3/yX8bI7h0iWptAm1/OLopmY8EuSB2ds3iBxl7mskFbT2TRBt6wyNteYHiQCNeTv6uR3xVo9aDQy/TIx4+X/48I9Pbn4FYu/xzGShQ/CXXhcU0C7XY5uyj5G8v5+Me6Bepqz11BPMJBy+gHuYqYDw8SwXpSO5/+Xclx6f4+p+mR5yEl0Pm5T9P72wi0A6u7nF2U7I3lnDyyYoG6qvOEXqDx543SD6LMhweJUD3JnU//rvCJ/XF6EiORFujo5sJLARmlR/ejeaMOrO9ydlGytrSVTh7JaOKJ5Qq0gT5L4gmuttkPtvM/RZgNDxKBeuY7XyLQIoPupyc5EoO3gJJlHImNWU7p8e22b9VKL2enynfsy/3fhJNHclmgQH1f9RxT8ppsxo3ihzwzZsODxLwedV6pfBQ8fZoxFAXnRMReDE1uX5ieHs1a6eXs8vOVBZreaIilCbSlPXdF8eRUphDocBRQBXqbQ31ckXxKDJo1FUUnlYUGrk6gHQy63svZ5eYrCjRjmyHQBfrynKYFFdxHV2BaFZPhQcJLoEWHoHmzUXhW7nzyKNBHmMvZZQY8uuTN1uDwM9Hj9kwE2tmeu0KBqopclUAL0ipNZcJxOHK6VXzi6XT8Ym/DD/4s1tO+WSu+nF1mwmOB3v6ct70gsALtI8wpRfEYC/Q2COgCVex0Yrc1mrt+csjy+wxCqAdSrqd5s1Z9ObusiIdHnUb/6oMKFMKeu1KBagrPOdYaDQ8SpQK9fHbIhOtH170FqlaoVM9+91s3C/Byds/KeJrN/sbb7e2Ogz+f/5i/qaVwebLeu44KVMUfPkGToHDMmlEyi9vTZ4fqOV7I02JDeY9mNZn7ne/cNZkFXs4uyPGw8/p20eQIVLGdGGhHoJ1e64xRGE/+eYN5DIcHidERqHanag9BL6Ny3U7uEWj+1UAiD5v1ILvkEWjjQ9A7uJxdkKNBz+LUnveZnhAogUK580ipQANfcrx2gcarN71m3YHx+4qNBZpp0OPVTVJFtG3WXVyNKUDsygOZXZQBEiiePsvjkb7kuGTVDoYHiUE9UvWlBt1uQ6+VTp+ktBZolkG3dyXQ57/Vl3N8EP1dCgQayTnnjssRaJv3UdWUX15v+gPFpzcDDIYHiVs9YvWqPX7Yc73fTKCBf2WbCzRp0O35cqNJgTY1qJdAf/Onh8/Df/gn31Sew/TYQKCXBgc2NOhCqksSEMa6HlKsSKBxgxYt29vwIHGtRy5es8cPDwODTgk+SdEKdHhF5G3ZE4LEU6WL9tMCbWlQH4E+/+npOkwfvrL5lPZiTA0Eeu3uZDPzGyTuH8PBWMICmDM+Q349Ag28j3RuRtn5O7fhQeJST6J4xawA61cAACAASURBVB4/PMQNmnqJRy3Q6VWZ8skw6O4+BPr6ZvOpP3/yiV8+/+uN+jz6lgIdRx28QeL+YeyNJR5CjJm95b4mgcYMWvi88To8SJzqSdY+3eV4Ag8PcYMmXyHPfgo/KszxhdCceho2y+k80G+dP8P569MHOlUPorz9Y41Ab1lHb5C4fwhzY0mHEGMCz8dWJNDY6fThA9C8jyNhCjQdxUyg0d2NCzTjDUb99FQINMegWfW0a5bPJ5G+cP0Q/JvHS4qoHkR5+0elQNPP9uAEKh1CDAm/5b4mgUY+0RluSd4n4iEFWpCMWqB5J2gUCLT42qO7HIOuX6Cnz8KfBfreE6SLiRyIXDlGczG0FH0EGlsQ6xJo3td7CD8d8Ygp0JJklALNPb+taHoqVlCyqLx6mjXL8XJ2R3NiXY3pQFigede8zqSHQOPrYVUCzfx6D/HHAx4hBVqSi3TQF/Zn3pi3FmjyZYXMelo1iwI9/WybuIWKDgIVxm5dAtVclmncyG3oe8sRBTqrXL+vE6bDc5wWUIGmDJpbT6NmuTyFP7xxdDYn2uXsii/+qqHJm0gvT4hXszKB5hs0V6CTn/TlGI9+HuURDvgTVqCJp/HZ9TTplsubSG+eLiRyEOhepihvIpl/eDiO+2lMU3nKE7cugSoMmiHQcDyXwRr/rQnHeApmNf8ul2mJP+kfXQC+x/SIxwO5G2nSLheBvvdk8/m3jgJ9/2sFl7NT3v4x9wi0nUGdT6RXfr59bQLNN2ixQOfoh7KMQ0FFk5p9n9vs4ApUMqhYz2iXWvTL63vhN5vPPHnhs5/e/1d7PWUngUZWjw/lV8vIMaPKnrsVClT5BUkXjiMwnwFlPPrp1HAsSPEEXj3Sw+GJHoL2F6gw5PkCbWFQp8/C/8uTzRm1P+9ToPJz8vyn68Fq1ifQIoNGTvAujEc/pUnOQnuqGFTtTI+nB1igcYOG6zk/Q5vskUOTJnhdTOT3P//M3p6fLLgcKLhAszainLmJGAMyrfHnCgVaZlBbgR7QT6rI5cPemjlVzfR0eJAFGn0aH6zn8h7BdI+MOzTnbq4HCifQpBPLjTmtZn0CLXoSby/QA4aTfzkC1e5ULrNRwn0N9Eh48kP1XM9SoUDvQ6BVR5RaVinQwtCK30RKUj/4l9nS1aN6971gkz2nJ1hxoJ7bedKzNOrbIuMn0H+9oL2usqtA86ZNGEoLgTbR5o01CrTYoHNM49HP7pVdmUCzn8TnjxyMQIM16wTqbVAfgb7/vc0NjE8inZPN/Pc6frNqgTa2526lArUzqHU8+vE9USzQvJEuHbq+0xOoWilQZ4O6CPTDVzYUaJj2+uy9BALY1IMq0MIlOyio4CGTY1k8dZ2nZ163KNDAFoq6kY3TJ5E2n/qLv7/wD0gf5ews0Pb23HVfAnOMBJr7NY4JPOLRj/CwoIIHTO1n+dz1np5Z5eKbSCH0zVDg9Fl47cc3Rw+iv0u2QHOvU5j1GY04kZnros/+S2CGUT05YXYS6E67bEcFFTxaYj8rBq/79ExrF09jCqIXSj4+V2NSf3xz9CD6uygEmte08gvCHgnPXCd/9l8CU6zqgRaoat2OCyp4LHk/awav//RMqhdPpA+iF0o+fpezKwZBoJWnPIV63EufCEtggplAc75JPHkTx3iyB3hcUMEjifsZn7y+8eQyrr+gHr1RsvG7In0x6xRoP38iLIExZvXkCHR7+cNu8odrNZ7x5M3vpKCCx5HmVZq89JwjTM9oD0rq0StFpZ5sst9E0n8CfvAg+rs4CLSgSzdmPe6oT4wlMMJOoPmHoGUCrf/sRdb8TgoqeBihUHH0liHQ0T7cgUD3h6DfKi2owbdy+jPtcVd/YiyBIXb1pGPtLdCctTstqOBR4oXKoxe43+RHGNNTKVA/g7o8hf/+1zebFz77jTN/hnQaUxsmPe6qT5QlMMCwHgSBJm6THt/JY5m2KzF7SxHowKByPbGO6aWiUE82uW8ibcBOpG94Mfoj4x539ifKErhhKdBUtHNvTochV6BZV88MkZre6WYs25WavcUI9NZpCjT1ILqbH4AWaG9/wiyBK2I9ymcS+gtk6AS6Df45sukI8vDONmP670sinuUIdJcr0G24UXqr5Ksnm6VejclKoLmbGfa4uz9xlsAFqR71SzGZ8ZYJdHS3yBRlnDMsDu/soSxfIk4eoC9IoOd9SdRDgSILNHM7gx739yfOErgg1PMoNzJEO4GmbhJHmN35dgxPUii9bOoAnOk5701pPWU2ylJPNrlvIg3eN3rv63/U/02kbgLt+/b7GZwlcCZez7Gbyq25Xm867t3ATQTiszvfjt3nDIq/+WQA0PSc9mf9Ah19Ekn/sSQYgc7vpRUohD+RlsCJaD3ndio3V5lxLJ5js7fjv0ZvlyQ2uoHtmF0qoPi794YgTU+VQH0M6i7Q9560EmjtVZAD99omfhDj1GMMfWItgSOxei7t1G6vLuaoQDPvXyXQ0HasLlZlM3xI03PcpfJ69GLJVE82SYFO3oA/8mKjp/CC3gqfws8FmnnHp30u/BkBaQkcidRz66d2g1VJR+JRzUzGjcOTGyxI8cBRzIYPanoOO1VRj94seerJJn0E+s5coN9WllT8FL7yEkpzphvM3n6b7zrKBWoJHAjXM+indoNS2MmmxQSqKaD0jaRwQZpHDmM4fVjTs9+tdQv0+f/8xje+/uT2MaRv/MU/aksqfw3U/HTPIoFi2XOHtgR2EUOMGqrepJB4qmvheHSzVPpGUrggzSMHsZw/sOnJFmijc5nWdTk7+0PQdE8moNlzB7cEwoYYN1S/zXjsRQJVTlLh66CRglQPPcd2AMGm5+WXawRqb1D305j0VAnU1qD67Z2GF2vmwJZAKJ5ZRwu2GvNGmUCDNxXeiS86BI0UlLEpAaU/yw7Q+1H03cxXCnyUoZ5sMk+kf/1zvyit57HuNCZ7g6pufpldrJlDWwLzeOYdLdpuWB2zHlacp2N7LlPkVnXtUh9+JkqHmx71R89GlBgpqZ5ssj8Lr33jaPQg+rvczgO1fg6vuvV1AWPNHNwSmNYTamnRhoMGnf+jOv5BMp7R5+HTt5HI2seqdumfvi9NoE/z9jB2LFVgpKR6skF/DbQrg9WLNXN4S2D0t3BLCzcdUGjq8xAZAo1eCy/9OaUJObtY066Clz/HV5yaV4M3PbkGDf64wEhJ9WSz+K/0cGP83hHWzAEugQGxnhZuO+cVQL1Aw/cc/D37qUp0B2+brhOo+i7j66XMdgRwevL+lWhzYTunr/R4seJFUHOB+l/LbvrWO9bMAS6BA4meFm9dewmiCoEWvOQe3UETgZa8/X554MhljBCnJ8ug0d7o/ZJSTzaZAv39zzabT+Jckd7ijSVhE4ETl7BmDnEJZDS1ePPJg1Cz10BLJiu2exYCLTp/6fzAscv1IU5P1o5Ge6P3S0o92SzzgsoGAg0M18tDJr/Dmjm0JfCYN3PlD5Ay6KmZl5ZWvAtfNFmR3TMQaNXl62L7gjY9p3iqTnXVCyahnmzuVaChf5yj9txRoHFO7cpqas3DyAo9XWfp0tRxPJphKTxpLrx39QL1+QAH0vQcMbhSj94wsnqyWeYFlasFGvancAesmQNZAoN25XW15sHSnwcLC9T6EiJBgjtXK1CvD8CBTM+Ncz01+6s3jKyebJYhUIvnWaP7h7YgNg9r5gCWwKRdeV2te8gchR7/M4zH4RPwIUL7VilQtw8QA0zPGIurlesVI6onG5VAn/9WX87xQfR3kQVaGHJse0fk1mHNXPclMGtXZlsrHzal0JBAVQ9g+rG3OoH6XX+h+/RMudWDYFAvgf7mTw8vfn74J9/UfyreWKDml2g6kGgc1sz1XQKhdmW2tfqxsy7sMoinalYCH3gq25y+XZ7XrwEWaMW/G3rHSOrJJlegz396evfow1c2n1J/KMlAoMPBdTkLNNE2rJnrtgSi7crta30JGQq9xVP3b+38HNFCg6rb5Xr9L2iB9jeoj0Bf32w+9edPPvHL53+tvyC9wZtILgedQ5Jr0vnxdXRZAmK7cvtqUUjSLgOBVj1Q6CT7nI/RzwvSPa7z5ROxBSouRTFyvWXCs+wh0Hc2m2+dPxH/6yftrkg/CM7ZoKmBxZq55ksg2a7cvtqUkxCMMp7SS9p5CdT78rPgApUMuliBvr75wvWSIm9uXlKWZHAak7NAkxOLNXNtl0BOu7Iba1OS7BgrgUq/EX87K0hRjf/Vu/EFWvhtLnrNBGfZ7WIiZ4G2+1bO7OSqSb8v4frwWhougcx2ZTfWqCxRM4YClS525yLQBt9+gC5QwaCJyPWeCc2y2+XszgLt8r3wnQ9A71ague3K76xVZYJozAQqXy00/3Wl/HpafHkMvEDj6zERud4zoVlepUB9SQ8t1sy1WgL57crvrFlxcdeoBXq5elHgV9Jd7AXa5Mu3SqbHdTXO6+l5COrzFP7wxtHZnO80+154g95kkTG19yhQTbsUtzWrL3yx+p0+ntvl32a/Ee8hrubhL3PraeLPgukxbtyUQD2xIBYq0NMbRyeB7mXa4U0kVzKm9r4Eqm+X4rZ2ZcYMeoin5BPw8ePN8E+sBdrGn+rpcejcmFA9hV/Iqp7cwCx7CPS9J5vPv3UU6Ptf26ivTm8oUI83k3LG9o4Eqm+WTqC+Br0KVHUNpvN/YgKd32Eb/Hlgo0fy2tXIn8rp8encuKDQDyMCTW2rZHons+xyIv2bm83mM09e+Oyn9//9grYkS4HqWpNFztiuS6D6fqTapbmxTQZHwt+XtI9HdeJwwQFoeutqgbbyp2p6/Fo3LCj0w8I0DGbZ57Pw//LkcjVQtT8NBdrrAJQCTbRLc2ObDE7MrHMVqH5bkfsUzpxWoM38qZge397dCgr+tJdBvS4m8vuff2Zvz0+WfEG8lUC3HgegeX1alUD17Ui2S3VroxSOTL3jINCSstQCbefP7Onxb96loPCPOxl0jdcDPePxgc7MwV2TQM3721Og0waeBFoyKZH7lM6cTqAN/Zk7PW3adywo8nMK9PIg+ruEn8J38+eaBGre3ketQD0NehRo0b+0tgegOoG29Gfe9DRr304QaBeD2gr0+fe/Maf7t3Lakduj9QhU34ucdulubpfEgZl9nrpfvEtJql1N/ZkzPW0bGK2ny4XtbAU6+To5kC+VMyN7cFcjUH0rstqlvL1hFru5f/A/qzimrT/T8bRuYFygPQxKgeaTP7lYS7LcEPpO5LVLeXvLMHazhdZQoHnHunI9jf2ZjKd5C+P1BILJSlw9wAPW/BqoMYrJXYlAzVt7aZf2DqZxTA26LIG29mcing4tFARaaFDtOA6hQHPRTO46BGre2Wu7tHewzWPSykUJtLk/xXi69FCop8MhKAWaiWpyVyFQ88be2qW+h3Eio2YuSaDt/SnE06mJYjyzH+VFrh7IK24Cff6vZ37zn9fwGqhuctcgUH0X8tulvod1JMN2NhVoznKWnqM292c0noK+GxUk/C58tYMMCvbmhI9A3//e2t5E0k3uCgSqb4KiXep7mGcyUFHLd+HrBNrDn+F4SjtvUpD0y1KBFu+Ri0DHb8a/uAKBKid38QLVt0DVLv1dzEO5yQhAoOMfm5+nU0UgnorWWxQk/TJ4uZg8ynbI6Xqgmxc++/Unh/9tXviWtiRAgWr/6V+6QPUd0LWr4D7mqVx7uhiBdvHnLJ7K3hsUJP62XKBle+Z0RfoX3zpfl/7NHt8Lb436qdPCBapvgLJdBfexj6WLQCM/zhFolyfw03gMul9dkPjbqpAK9sbnO5GOF1F+83gpu9ej3wv/8U9e3W6/86v5g2gqOuErUP3kLlug+vy17Sq5k30u57YCfBIpR6Cd/DmKx6b9tQXJv65KSb8zfl8q9/jO8cs83ol9pcdHr20P/PE/zx5EU9EJV4EWTG7/JTlEaQh9/Op2Fd3LPphTY5ch0F7+HMRj1v/KguRfV+ak3RdXgR6evX/4SuQ5/Bvbr/zq8YMfb7/y79MH0VR0wlOgJZPbf0kOKf9SBidQBLoDEeg2LdBu/rzGYzoBVQUlfl8ZlHJPfF4DPT6Ff+/J4NuNZ/zh1eOx50evfenvpg+iqeiEo0CLJrf7khyhMoQ+/IJ2ld3NPpljbxEEOvqr4lvTGnCOx3YCanqZvlhVxcZ32j11eRf+9eOrn6eXQk8anfP29qvn//7l9EE0FZ3wE2jZv/zdl+SIim+18QFGoLt+Ah1KMy3Qfv48xWPb/7pupi+XmrWZ+Nvzqt3w+lbOz/3i8Db8Fw4yDT+Ff2P7o+N/3z2LdPAgmopOuAm08JnTUgWqT76IQoE6HYJ2Eeho+SYF2u8J/DEe0+ZXdzPjetM5m9nGPxam2QufTyK9fvz80TubzQtPIl8r9/GPz0/d//Dq5UXQL154VsZTe/aD67BVVApzb4fDPndp8Hb7dKu5fecxXFA7TyTSOme/F2jsFl57/EzxrZx7gT5//fhBpOAB6CIEek/+NJ0TLxz2u0ePt8LaDbBSf/r088A4r1nW1x9IXfDZYcXFRP7P3pvPf/2Zz3wz/B78QKDTE5lgnsKXP3Fa3lN4feYVlD6Fd3kVtMfTY/njMZN29Xz+vh8M3XPUBg3N+dLS4V9mYV9/ILchbw96Xc4ucAR6fRD91jwEWrOylibQ4jaWUb4kPeLpotAjkQU8blc/f56b5SnQkoZmfevz4C+lAs1bFRRohKpltSyBFjexFCyBHp4i95FUpkBblDLn2ixXgRa0NGdxSQLd5go0a2F0u6Ay+LvwdUtqIQIt714dFUvSJ56X+zg0R6Cd/DlolrNA1T3NEuggtem77baXF+km0Mv5n5jngVYupyUItLx11cAJdNdHoRkCbV7UvFnuAlU2NWtxCYegmm+xTtfeTaDIn0SqXkvoAi1vmwk1S9IvnvYKPR4bzZfzoF0tK4o1y1+guq7mCXR0CDr6lUag6ZXSTaAf/3j7ZdDPwtcvJFiBlvfLEkyBdlDoUaCz9XxrV5tyEs1qIVBNW/MW1yC5QoGebpcqvN+Xyn2AejUmg6kFFWhFt0ypWpJ+8eyaK1QWaINacprVRKCKtmYuroFAlaGc73FujGk8pt/K+cFP9v78zvT4s7tALaYWUqA1vbIFV6CNT2oSBepcSXaz2gg0v6+5Ai0P7/S206UxlvHcwdcam0yth0DL9/BgCOsmVVG3JE1T3c1fIm6o0G3wizqfXsrweVBtsxoJNLuvuYurPL/t8Ag0ERgFOk3dYmztBXraz7L7Pm21BHLBFmj7g9AJl7MC7B+sqFntpidvJ7IXV7VAcwxKgY4ztxlba4He9rTkrg2XQB6V9RiHGzjLq+VBaEigHi/FFjer4fRk7Ui+QJMRht9PupiTAlVjNbbGAh3uasG91iZQY4OGTpN9+YLtQwUICdTlgYub1XJ6cnYkf3GNUzy+Zz37VFLgbpcf3n4pVEyBjgI3mltTgYZ2WHUXClQk/DmDlxs5dL6EcZ69n5rVdHoy9kSxuMan05cLVEhvrQItmUK7wbUUaHSns2+/OoHaGlS6VECr49DJIzpst7xZbacnvScagQ5Ppw/dQBZoTnwrEujLMbLTNhtcO4Gm9jzn1usTqKlBE9daaavQw2M5ncNR2qzW05PaFU08p8bFT53P/1RSrNq1CvSUnkailqvEagnk7X3qxhSoSPpqf80cenoYB4HWNKv59CT2RRXPqWtRT9ZfWGRFAo2HmLMATFeIyRIwC3SFArU0aM71pvXPZ0o4b/7eBZpork6g8iEoBXoi/SZSavptV0bdErAOdI0CNTRo5pfKlb4ylM9lm15nERc2q8f0SDuji8dMoCavjy1XoDvpKMJ8PVQsAfM4KdBUs5Qr0kmhtw2aC7SuWV2mR9gbbbt2ygswRQkWekcCPcVZ+U5THqVLwDzLWzpImNRj16yybhlPzWBrFOij1F5lPIdYKdAE2vNAvZ+NFS0B8yBH6SBhU49Zs0qFZTk7ww35fY6tqFmdpie6O9p4zAQaDPI+BeqNfgmYpzhNBwkcgT7sKRZo6AlN+YZuf3H8HFtJs3pNT2x/1AKt/Ez8gECVFKgHyh6bRxhIBwmjeur79HCkTlgWLwiN70GBnojsjzoeTTfGzpwaNFAkBeqBosfm8UXSQcKqnto2PZyp3c6AIolObmwr0NqQO05PeIf0AlW0InHx+nmNFKgHuT02zy6eDhIgAn14cDDoCZVEp7dzv5KCiq7TE9ojfTz5Bk1++8esQgrUg5wem+cmp4OEWT11XXIU6IFchc5u1ORKCtn0nZ7ALhXEky/QyV8p0D4kemweWUY6SNjVU9UlZ4HuMi9pM7sJBTpgvkslAo12QX7NM/Du/bQ+CtQDucfmieWkg8TdCDRHofPfGwrUIOLe0zPbp5J4cgUq/jYYKQXqgdRj87zy0kHCsJ6qNvn7M32+feC3FOgYg3iiLRgpcuZLCrQT8R6bp5WbDhIoAjU5jSmJqNDQr9pdDDEHgOkxiCfrEDTvdPtxbRSoB5Eem0elSAcJy3rqGlV5In0mcYUGf25Wj0nACNNTH0/eIWhBqhSoB+EemyelSQcJ03rqm+Uu0OhLoeGfUqAz6uOxE+g4VgrUg1CPzXPSpYOEbT3VzWog0LBCI8elbS/HnQJjeqrjiR2ClnxMflgYBerBvMfmKWnTQcK4ntpmNRFoyKCRRW1Uj1G6INNTHY+hQIfRUqAeTHtsnpE+HSSs66lsVhuBzg84YwdFFGiI2nhMr7l2K4sC9WDSY/OICtJBwryeuma1EujkIDT63rxNPVbZwkxPbTwUaBxogZoHVJQOEvcq0JFB4+c28Ru1wlTGIx6CKp/J34qiQD0Y9Ng8Hj04S+CMQz01zWon0NsZTdK5oRb12CULND2V8UgGjQk08vNbTRRoKXnfOmCeTgFAS+CERz3ljWwq0MH17qI3MajHMFik6amMp0CgkV/cSqJAS3kUHvLp4Db9QVoCR1zqKW5kW4Hu0p/urK/HMleo6amLR0idAm0t0PMjh3+JNHNgS+CATz2lnWwu0BTV9ZjGijU9dfHkXpUp+YtbQRRoIZeHDsYKNXNgS+DRL56yTq5OoLahgk1PXTz6T3RSoD4k6inacS/AloBnPCWtXJtAjSNFm566eGYCPQuSAqVAo6AtAc94ClpJgYrATU+dQCMGpUAp0ChwS8C1HnUrVyZQ6zzxpqcqnplBE2eAUqAupOop2W838JaA69a1vVyXQM3jBJyeqnSVAk2lTIEWkaqnZL/dAFwCvuh6uSqB2oeJOT3lCU0PQSnQx/YCTdZTsNt+YC4BT1TNpEBFUKenOCIKdA4FKoG6BBzRNHNNAnWIEnh6yjKaHIJSoI8UqAzwEnBD0cwVCdQjSejpKUppbFAK9LG5QNP16HfBEegl4EduN9cjUJcYwaenJKf6C9sNqqFAa+KL1aPfBUfAl4AjWd1cjUB9MoSfHn1Q9ddWHlRDgdbEF6tHvwuOwC8BRzK6uRaBOiW4gOlRR1V9KtOgGgpUTUY9+l1wZAFLwJFkO1ciUK/8FjE9yqw0pzIFfzeohgJVk1GPfhccWcQScCTRznUI1C29ZUyPMi3NO/GhXw6qoUDVZNSj3wVHlrEEPBHbuQqB+mW3lOnR5UWBjmgq0Jx69LvgyFKWgCdCP9cgUMfkFjM9qsAUpzJRoMbk1KPfBUcWswQ8ifdzBQL1DG4506OKjAIdQoFKLGcJuBLr5/IF6hrbcqZHldntEHQvSAq0oUCz6tHvgiPLWQLOhBu6eIH6hrag6VGlNjQoBUqBCixoCTgTbOjSBeqc2YKmR5fb8BBUggI1Ja8e/S44sqAl4E6gowsXqHdiC5oeXczXQ9CCT8MPqqFAS6OT6tHvgiMLWgINmHV00QL1z2tJ06PL+WJQCrShQDPr0e+CI0taAg2YtnTJAm0Q15KmRxk0BXqBApVY0hJowrilCxZoi7AWNT26oM+HoBRoO4Hm1qPfBUcWtQSaMOrpcgXaJKtFTY8y6ZNBKVAKVGRRS6ANw54uVqBtolrU9GijLr2u3aAaCrQot0Q9+l1wZFFLoBW3plKgIsuaHmXUpVcGHVRDgRbllqhHvwuOLGsJNOPS1KUKtFFMy5oebdaFBh1UQ4GWxJaqR78LjixrCbTj3NWFCrRVSgubHm3YGQINvEY6qIYCzUZRj34XHFnYEmjHqa3LFGizkBY2Pdqwcw5B5x/1HFRDgWajqEe/C44sbAk05NjWRQq0XUYLmx512hkG3W6nH5cfVEOB5qKpR78LjixsCbTk0NclCrRhREubHnXcOa+CUqD1qOrR74IjS1sCLdktUqAtE1ra9KjjznsfiQKtRVWPfhccWdoSaMoSBdo0oKVNjz7vgnfiB9VQoNrIcurR74IjS1sCbVmeQNvms7jp0QeuN+igGgpUmVhWPfpdcGRxS6AtSxNo43gWNz36wF9WG3RQDQWqTCyrHv0uOLK4JdCWZ8sSaPN4sLqVMT36xNUGHVRDgeoCy6tHvwuOLG8JNOXZM8/JKUAUaId4mj+kjINA1QYdVEOBqvLKrEe/C44sbwk0ZR+P4+gUIAm0SzxYeAhUa9BBNRSoKq/MevS74MjylkBTDvE4zo4eQaCd4oEiXU9J6C+rFDqohgLVxJVbj34XHFngEmhJ0+/EziEu0G7xIOEjUJ1BB9VQoIq0suvR38WRBS6BlpzicZsePVGBdowHiIx6ynJPKHR4Kv2gGgo0SUE9+rs4ssQl0JBzPG7joyYm0K7x4OAmUIVBB9VQoCkKdgFr5pa4BBpyicdrftREBNo5Hhj8BJpQKAVaRsEuYM3cEpdAQyjQrHhgyKmnOHvJoIsR6LMynvpQWA1ZHk4TZEPvcBZFRc4HhYZ/s91WN+MOj0CLdgHrH+1FHkO0YxCPzwhpCR6BIsSDQVY9FfFHD0IXcwSqvwsFKrHMJdCMYTw+M6Qk2wA6GgAAEe5JREFUJFCMeCDwFmhUoRRoAWW7gDVzy1wCzRjF4zJESihQEX+BHg06dygFqqdwF7BmbplLoBnjeDymSElAoL2yeVzs9FT2IGjQUD8oUJHCXcCauYUugVZM4vEYIx1zgXZK5shSp6e6DZJBB9VQoAKlu4A1c0tdAo2YxuMwRzpmAu2Ty5mlTk99H4SD0EE1FGic4l3AmrmlLoFGzOKxH6QJD3uEX08F2iWVK4udnvo+xQ06qIYCjVO8C1gzt9gl0IZ5PPaTNOLhQTboRKA9Mhmw2OmxaFXMoINqKNAo5buANXOLXQJtCMRjPkpDHh4SBh0LtEMiI5Y7PRbNihh0UA0FGqNiF7BmbrlLoAmheKxnacDDQ8qgI4G2z2PCcqfHpF3hp/GDaijQCDW7gDVzy10CTQjGYzxMA3QCbZ7GjAVPj0m/ggYdVEOBhqnaBayZW/ASaEE4HttpGqASaOssAix5emw6NjTo+WT6QTUUaJC6XcCauSUvgQZE4jEdpwEagTZOIsiSp8eoZQODUqB5VO4C1swteQk0IB6P5UTdyH8TqWUKURY9PUYtuxl0ezLooBoKNEDtLmDN3KKXgD9SPIYzdSP3NKZ2EUgsenqsWjY06Lg1FGiA6l3AmrlFLwF/5HjspupG1on0rXY/xbKnx6plV4NSoGnqdwFr5pa9BNxJxWM2V7k8xbHn49Knx6wpF4OensMPqqFApxjsAtbMLXsJuJOOx2qyMudvafE0RlePWVtuBt1RoCIGu4A1cwtfAt7kxGM1W1nTt7x4mqKsx6wzF4Geu3SphgKdYLELWDO39CXgTF48RtOVM3tLjKchvQS6G5wOOqiGAh1jsgtYM7f0JeBMbjw285UxesuMpxnaesz6MzgddFANBTrCZhewZm7xS8CX7HhMBixj8hYaTyu6CXRg0EE1FOgIm13AmrnFLwFflPGYjJk4douOxx91PXZtuhp0UA0FOsRoF7BmbvlLwJWCeCxGLT52i4/Hl44C3VGgMla7gDVzy18CrhTGYzBukalbRTx+dP337mzQQTUU6A2zXcCauRUsAU8q4qkfudDQrSYeHzoL9OVx0yjQAWa7gDVzK1gCnlTGYz9yq4rHnoJ66tVw5WTQQTUU6BW7XcCauTUsAUcc4ymbuLuJp4y+Aj0ZdFANBXrBcBewZm4NS8CRNvEo5u0e41FQUk+tG4YcDDqohgI9Y7kLWDO3iiXgB+MRWUU8Jua8sDfooBoK9ITpLmDN3CqWgB+MR2Qd8dio88zQoBToCdtdwJq5dSwBNxiPyDriMVLnGR6BTjHeBayZW8cScIPxiKwkHit3Tn1BgY7zsAFr5layBLxgPCIricfMnRNhUKA7e39yCchg1cN4RNYSj509dxRoLA0rsGZuLUvACcYjspZ4DPVJgfrqk0sgAVY9jEdkNfFQoHGKBWpe/Lkep+2WsZol4APjEVlNPBRonDKBmld+q8dv0wWsZgn4wHhEVhMPBRqnQKDmZY/qcd26ltUsAR8Yj8h64qFAo2gFal7ztB7vB1CxniXgAuMRWVE8FGgMpUDNS57V4/4IGla0BDxgPCIriocCjaESqHnBgXoaPEY+K1oCHjAekTXFQ4FGUAjUvNxgPU0eJZc1LQEHGI/IquKhQMNkCdS8UKGeho+VZlVLwB7GI7KqeCjQMBkCNS9TrKfpo6VY1RKwh/GIrCseCjRIUqDmRSbqafx4MutaAuYwHpGVxUOBhpAFal5hup72DymwsiVgDeMRWVk8FGgISaDm9eXU0+NBo6xsCVjDeETWFg8FGiAqUPPiMuvp9Lhh1rYEjGE8ImuLx0CfD3uu1axYoOaV5dfT76EDrG0JGMN4RFYXj4k/rwZdr0DN69LU0/PBZ6xuCdjCeETWF4+JPy8GXatAzavS1dP34SesbwmYwnhE1hePjT/PBl2pQM2LUtbT+fHHrG8JmMJ4RFYYDwU6ZSpQ85q09fQuYMQKl4AljEdkhfFQoFMmAjUvSV1P7wJGrHAJWMJ4RNYYDwU6YSxQ84rUYM3cGpeAIYxHZI3x1Ah09W8imddTANbMrXEJGMJ4RFYZj4FBL9WsTKDm1RSBNXOrXAJ2MB6RdcZTbdBrNasSqHkthWDN3DqXgBmMR2Sd8VQJdKUf5QQCqxy0dBiPDFY9K42HAh2C1mSsctDSYTwyWPWsNR4KdABak7HKQUuH8chg1bPWeCjQAWhNxioHLR3GI4NVz2rjoUBvoDUZqxy0dBiPDFY9642HAr2C1mSsctDSYTwyWPWsNx4K9Apak7HKQUuH8chg1bPieCjQC2hNxioHLR3GI4NVz4rjoUAvoDUZqxy0dBiPDFY9a46HAj2D1mSsctDSYTwyWPWsOh4K9ARak7HKQUuH8chg1bPqeCjQE2hNxioHLR3GI4NVz7rjoUCPoDUZqxy0dBiPDFY9646HAj2C1mSsctDSYTwyWPWsPB4K9ABak7HKQUuH8chg1bPyeCjQA2hNxioHLR3GI4NVz9rjoUAf8ZqMVQ5aOoxHBque1cdDgeI1GasctHQYjwxWPauPhwLFazJWOWjpMB4ZrHrWHw8FCtdkrHLQ0mE8Mlj13EE8FChak7HKQUuH8chg1XMH8SxXoB+99tXAg+i3g9ZkrHLQ0mE8Mlj13EM8ixXoG1sK1B+0dBiPDFY9dxHPMgX68RtbCrQBaOkwHhmseu4jniUK9Hc/3FKgLUBLh/HIYNVzJ/EsT6Bvb7ff/TcKtAFo6TAeGax67iaexQn0y3/7+C4F2gC0dBiPDFY99xPPwgR6YCTQL154RgghzXmaS+kDUKCEkNWyaIFeRarfDtrTDKxy0NJhPDJY9dxZPAt+Cn99EP120JqMVQ5aOoxHBquee4sHXKDvbo/86PI3CtQftHQYjwxWPfcXDwXaHaxy0NJhPDJY9dxhPMgCHUOBtgAtHcYjg1XPPcZDgXYGqxy0dBiPDFY9dxoPBdoTrHLQ0mE8Mlj13G08FGg/sMpBS4fxyGDVc8/xUKCdwCoHLR3GI4NVz33HAy7Q8IPo74LWZKxy0NJhPDJY9dx7PBRoB7DKQUuH8chg1cN4KNDmYJWDlg7jkcGqh/EEDDqohgJ1AKsctHQYjwxWPYznca7QQTUUqANY5aClw3hksOphPAco0KZglYOWDuORwaqH8RyhQFuCVQ5aOoxHBqsexnOCAm0IVjlo6TAeGax6GM8FCrQZWOWgpcN4ZLDqYTxXKNBWYJWDlg7jkcGqh/HcoEAbgVUOWjqMRwarHsYzgAJtA1Y5aOkwHhmsehjPEAq0CVjloKXDeGSw6mE8IyjQFmCVg5YO45HBqofxjKFAG4BVDlo6jEcGqx7GM4EC9QerHLR0GI8MVj2MZwYF6g1WOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwgF6gFWOWjpMB4ZrHoYjwiiQAkhZJVQoIQQUkgDgRYwrYsMYToijEeE8YjUxUOBLgGmI8J4RBiPCAW6fpiOCOMRYTwiFOj6YToijEeE8YhQoOuH6YgwHhHGI0KBrh+mI8J4RBiPCAW6fpiOCOMRYTwi6xAoIYQsDgqUEEIKoUAJIaQQCpQQQgqhQAkhpBAKlBBCCqFACSGkEAqUEEIK6SnQj3/y6nb7nV8lfnSv/O6/b7dfGkXx0WvbI3/8z92KgiGUBYfnzMc/3l645cPpOfPRa189/8lCQB0Feu7osKGBH90r/3Qa9y/93e1Hf3iVS+BCIAsOz4WgQDk9Z97YngVqIqCOAn1j+5VfPX7w4+1X/l360Z3y7vZL/+PxEMWgl+9eOk9CWXB4pvzh1cG/v5yeIx+/sb0EYSKgfgL9w6tHN3z02q3JgR/dKftjiB8d/rv/F/FH1x++MfjzvTPPgsMzZT9Ef3n7G6fnwO9+uL0I1EZA/QT69nk/3r51OfCjO+Wj185Hnm/covj4x1TDhUAWHJ4pbw8PpDg9B97ebr/7b9dBsRBQP4Fe/kUcPLUI/OjeGQj0o9e+/L/2/37+V75LEsyCwzNh9OyF03Pk7S//7XVAbATUTaDXfxH/8Orl38nAj+6d8fOL8xsDfCoWyILDM+XtkQU4PRfOdjQSEAUKzXARvLt/+vHvj//vJ1s+FwtkweGZMHkhj9NzYX0CvbzRHPjRnfPucNzfvj714It88yw4PBPeHb+VzOm5MBdojYAQBMoj0AjvvvqlwBOud3mizpVbFhyeMZfzOKZwetZ3BMo1EObt8NMtHmLdCB49cHge41PC6VmJQPkufIp/irxcRUPcGGTB4RkRi4HTs5J34a/nWo1Ow5r96G75+I3tl0eHCtfnZDREMAsOz4jJa52cnivvTk72rBMQP4kEyhuz16reuDz14JkooSw4PENm581zei5cBLr0TyLtW/nlyedOAz+6V96eZ/CHVw8nonzwQ6YTzILDM+T6UbYLnJ4LF4HaCKjjxUQ+GFz55HzZgw/UF0NZKZdrj50+uHtO5+3z5XTu+8MkJ4ZZcHjmDF7q5PSMub6IYSKgntcD/eAn+1q/c+zz5boxgx/dNe9e/TkQ6OMHh0uEfpfpHBhkweGZMzhbidMz5vYqsIWAeEV6QggphAIlhJBCKFBCCCmEAiWEkEIoUEIIKYQCJYSQQihQQggphAIlhJBCKFBCCCmEAiWEkEIoUEIIKYQCJYSQQihQQggphAIly+D1zYtv9a6BkAkUKFkGuQL99R/Rs6QZFChZBpkC5YEqaQkFSlYFBUpaQoGSVUGBkpZQoGRVUKCkJRQoQeL5X21eevzN1zabFz6/9+Bv/nSz2XzuF8ffnMz44SubLzz/+ac3m09+8+jJw99P93xz84lf7v/vyPFH73/vyX4z53uf/7r5j3/Tfp/IiqFACRJ7gf6H108W/MQvf3r6wws/OPzmKtD/9LXTjz/1y0dJoOc/bjafP/7615e/Hu9GiA0UKEFiL9C98t56fL6X5yc3n/vt4/OfbU5Pyq8C3bzwX956fP+vTpqcCPT2FP7N06Hr7396Muh7TzYv7v962NpLvfaNrBAKlCDx/CzGoylPqnt9aMbDj799vuHlKf3pnmOB7oX50uXnhyPY82+vWyPEBAqUILH34tV0p6fuZwPeBHp+k+hkwphAr8I8vap62QghtlCgBInzgeXj4FDxnYlAr0eWgkBP1nw8b2j/k3cO7yf9Y8s9IXcBBUqQuJkvKtCRMCMCPTzTv3H4xekNqU9+87fNd4msGQqUIOEm0NMpUYd34XkiE7GDAiVI2Ak08Gb7//3+pw8K/bb7XpC7gQIlSJQL9PXJa6DhDyRdz4oixAIKlCBRINCXLnccvgt/fQ//5NLBm0o8j4kYQoESJLQCvR5qvrOZngd6PtB85/ic/XoaU/TYlJACKFCChFagh08cvfiLx+c/f7K5/P2FHzz/7fHnh7eLLk/Z9/d74Vv7/x4+wcTXQIkZFChBQi3Q42c/97z4s9Pf3zn85aXHwWfhTwec7z25/P0LzXeKrBcKlCChFujj6dpM37r+/dd7U750cOb73/v08PJLz3/+mcOJoJ/niaDEEAqUEEIKoUAJIaQQCpQQQgqhQAkhpBAKlBBCCqFACSGkEAqUEEIKoUAJIaQQCpQQQgqhQAkhpBAKlBBCCqFACSGkEAqUEEIKoUAJIaQQCpQQQgqhQAkhpBAKlBBCCqFACSGkEAqUEEIKoUAJIaSQ/w9zS7O53CwR/wAAAABJRU5ErkJggg==)

The band is the 95% interval of the smoothed alertness, the dashed line
the generating alertness and the points the ratings.

The same structure links decisions to any process ctsem can model, such
as a mood or stress rating, a physiological signal, a process with its
own predictors, or several drift diffusion tasks loading on shared
latent processes.
