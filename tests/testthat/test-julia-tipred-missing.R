# Sampling missing time-invariant (TI) predictor values, julia backend.
#
# SPEC-tipred-sampling.md (review/, not shipped with the package) is the
# design document and carries Charles's 2026-09-03 decisions, which this file
# verifies directly. The engine-side correctness argument -- the closed-form
# check at two draw counts, the gradient against an analytic derivative and
# against FiniteDiff, the adjoint guard -- lives in
# inst/julia/ContinuousTimeSEM/test/test_ti_missing_predictor.jl; this file is
# the R side: that `.ctJuliaPrepare()` builds the right spec, that the
# fallback rule and its warning fire correctly, that `optimcontrol$tipredMissing`
# chooses between imputing and sampling on every route, that `ctFit()`
# actually samples end to end, and one targeted comparison against Stan's
# independent implementation on the smallest model that exercises the path.

.tipred_missing_model <- function() {
  model <- suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix("t0m", 1, 1), n.TIpred = 1, TIpredNames = "group",
    tipredDefault = FALSE
  ))
  model$pars$group_effect[model$pars$param == "t0m"] <- TRUE
  model
}

test_that("a sampled fit places with a missing TI predictor imputed, then samples it", {
  model <- .tipred_missing_model()
  dat <- data.frame(id = rep(1:3, each = 3), time = rep(0:2, 3), Y1 = 0,
    group = rep(c(-1, 2, NA), each = 3))

  prepped <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    optimize = FALSE, fit = FALSE, intoverpop = "augmented")))
  expect_s3_class(prepped, "ctJuliaModel")
  # The placement optimises with the value at its imputed mean. Single
  # predictor, includeOutcome defaults TRUE but there is nothing to condition
  # on here (constant Y1) -- falls through to marginal: mean and SD of the two
  # complete cases.
  expect_null(prepped$ti_missing)
  expect_equal(nrow(prepped$ti_imputed), 1L)
  expect_equal(prepped$ti_imputed$subject, 3L)
  expect_equal(prepped$ti_imputed$predictor, 1L)
  expect_equal(prepped$ti_imputed$mu, mean(c(-1, 2)), tolerance = 1e-8)
  expect_equal(prepped$ti_imputed$sigma, sd(c(-1, 2)), tolerance = 1e-8)
  expect_equal(prepped$tipred_data[3, 1], mean(c(-1, 2)))

  # The sample makes it a parameter, its raw index past every other block:
  # the model parameters (t0m is indvarying, so the augmented route also
  # carries its population SD), then the one TI-effect coefficient, then this.
  sampled <- ctsem:::.ctJuliaTIMissingSampled(prepped)
  expect_null(sampled$ti_imputed)
  expect_equal(sampled$ti_missing$parameter, max(prepped$ti_effects$coefficient) + 1L)
  expect_equal(sampled$ti_missing$mu, prepped$ti_imputed$mu)
  # Never NA in what reaches the engine -- the missing cell holds the
  # placeholder, always overwritten before it is read.
  expect_false(anyNA(sampled$tipred_data))
  expect_equal(sampled$tipred_data[3, 1], 99999)

  # Complete data still prepares on the same path with nothing new, so the
  # new machinery is about the missing cell specifically.
  dat$group[dat$id == 3] <- .5
  prepared <- suppressMessages(ctFit(dat, model, backend = "julia",
    optimize = FALSE, fit = FALSE, intoverpop = "augmented"))
  expect_null(prepared$ti_missing)
  expect_null(prepared$ti_imputed)
  expect_equal(as.numeric(prepared$tipred_data), c(-1, 2, .5))
})

test_that("fitting a missing TI predictor with the default (adjoint) gradient now works, same as explicit forward", {
  skip_without_julia()
  # The reverse pass initially had no cotangent for a sampled TI predictor
  # value, so a fit here briefly refused unless gradient='forward' was
  # requested by name (see the git history of .ctFitJuliaBackendImpl). The
  # Julia engine's adjoint now covers this (test_ti_missing_predictor.jl,
  # cross-checked against ForwardDiff and FiniteDiff), so ctFit() no longer
  # singles this model out: the default 'adjoint' and an explicit 'forward'
  # both fit it without error, and are two gradient methods for the same
  # posterior rather than one being the only offer.
  model <- .tipred_missing_model()
  dat <- data.frame(id = rep(1:3, each = 3), time = rep(0:2, 3),
    Y1 = rnorm(9), group = rep(c(-1, 2, NA), each = 3))

  fit_default <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    optimize = FALSE, intoverpop = "augmented", chains = 1, iter = 10,
    cores = 1, sampleControl = list(warmup = 5))))
  expect_s3_class(fit_default, "ctJuliaFit")

  fit_forward <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    optimize = FALSE, intoverpop = "augmented", chains = 1, iter = 10,
    cores = 1, sampleControl = list(warmup = 5),
    optimcontrol = list(gradient = "forward"))))
  expect_s3_class(fit_forward, "ctJuliaFit")
})

# Missing TI predictor values are imputed when optimising and sampled when
# sampling, on every route, unless `optimcontrol$tipredMissing` says otherwise.
test_that("optimcontrol$tipredMissing chooses imputing or sampling, on every route", {
  model <- .tipred_missing_model()
  dat <- data.frame(id = rep(1:3, each = 3), time = rep(0:2, 3), Y1 = 0,
    group = rep(c(-1, 2, NA), each = 3))
  prep <- function(...) suppressWarnings(suppressMessages(ctFit(dat, model,
    backend = "julia", fit = FALSE, ...)))

  # Sampling needs no particular route: 'auto' samples the effects.
  defaulted <- prep(optimize = FALSE)
  expect_identical(defaulted$args$resolved$intoverpop, "none")
  expect_identical(defaulted$tipredMissing, "auto")
  expect_true(ctsem:::.ctJuliaTIMissingSamples(defaulted))
  expect_false(ctsem:::.ctJuliaTIMissingSamples(defaulted, "impute"))
  expect_false(ctsem:::.ctJuliaTIMissingSamples(prep(optimize = FALSE,
    optimcontrol = list(tipredMissing = "impute"))))

  # Optimising imputes unless asked to estimate the values as parameters, and
  # a fit that estimated them cannot then be sampled with them imputed.
  expect_null(prep(intoverpop = "laplace")$ti_missing)
  estimated <- prep(intoverpop = "laplace",
    optimcontrol = list(tipredMissing = "sample"))
  expect_equal(nrow(estimated$ti_missing), 1L)
  expect_null(estimated$ti_imputed)
  expect_error(ctsem:::.ctJuliaTIMissingSamples(estimated, "impute"),
    "estimated its missing TI")

  # Refused by name: a value that is not one, the setting on stan, and the
  # setting in the sampler's list, where ctFit() does not read it.
  expect_error(prep(optimcontrol = list(tipredMissing = "drop")),
    "must be 'auto', 'impute' or 'sample'")
  expect_error(suppressWarnings(suppressMessages(ctFit(dat, model,
    backend = "stan", fit = FALSE, optimcontrol = list(tipredMissing = "sample")))),
    "tipredMissing")
  expect_error(ctsem:::.ctBackendSampleCheckControl(list(tipredMissing = "impute"),
    "sampleControl"), "optimcontrol$tipredMissing", fixed = TRUE)
})

test_that("the Laplace objective's gradient carries a sampled TI predictor value", {
  skip_without_julia()
  # The value moves the subject's raw parameters through its TI effect, and
  # with them the subject's mode and curvature: central differences against
  # the engine's gradient, on every coordinate, at a point off the imputed
  # mean so the imputation density's own term is not zero.
  model <- .tipred_missing_model()
  model$pars$group_effect[model$pars$param == "drift"] <- TRUE
  set.seed(5)
  dat <- data.frame(id = rep(1:5, each = 4), time = rep(c(0, .7, 1.9, 3.2), 5),
    Y1 = rnorm(20), group = rep(c(-1, .4, NA, 1.2, -.3), each = 4))
  spec <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE)))
  spec <- ctsem:::.ctJuliaTIMissingSampled(spec)
  npar <- ctsem:::.ctBackendNpar(spec)
  x <- seq(-.2, .2, length.out = npar)
  cell <- spec$ti_missing$parameter
  x[cell] <- spec$ti_missing$mu + .8
  value <- function(x) ctJuliaEvaluate(spec, x, gradient = FALSE, cores = 1L)$value
  g <- as.numeric(ctJuliaEvaluate(spec, x, gradient = TRUE, cores = 1L)$gradient)
  h <- 1e-5
  fd <- vapply(seq_len(npar), function(j) {
    e <- replace(numeric(npar), j, h)
    (value(x + e) - value(x - e)) / (2 * h)
  }, numeric(1))
  expect_equal(g, fd, tolerance = 1e-5)
  # And not the imputation term alone: the effect reaches the likelihood.
  expect_gt(abs(g[cell] + .8 / spec$ti_missing$sigma^2), 1e-3)
})

# The optimising path imputes a missing TI predictor from the other variables,
# as the Stan path does, and the new missingness detection must not reach it:
# SPEC-tipred-sampling.md is explicit that this path is not to change. What
# the warning has to carry is that the value was manufactured and that it is
# the predictor's own effect estimate that pays for it.
test_that("the julia optimising path imputes a missing TI predictor, with a warning", {
  model <- .tipred_missing_model()
  set.seed(20260902)
  dat <- data.frame(id = rep(1:4, each = 3), time = rep(0:2, 4),
    Y1 = stats::rnorm(12), group = rep(c(-1, 2, NA, .5), each = 3))

  warned <- character()
  prepped <- withCallingHandlers(suppressMessages(ctFit(dat, model,
    backend = "julia", optimize = TRUE, fit = FALSE)),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  told <- grep("TIpreds", warned, value = TRUE)
  expect_length(told, 1L)
  expect_match(told, "1 in group", fixed = TRUE)
  expect_match(told, "imputed", fixed = TRUE)
  expect_match(told, "overly confident", fixed = TRUE)

  # Imputed, not sampled: no sampled cell, and a value that is neither the
  # sentinel nor missing.
  expect_null(prepped$ti_missing)
  values <- as.numeric(prepped$tipred_data)
  expect_equal(values[c(1, 2, 4)], c(-1, 2, .5))
  expect_true(is.finite(values[3]))
  expect_false(values[3] == 99999)
})

test_that("julia optimises with the imputed TI predictor values", {
  skip_without_julia()

  model <- suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix(-.4, 1, 1),
    DIFFUSION = matrix(.5, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix("t0m", 1, 1), n.TIpred = 1, TIpredNames = "group",
    tipredDefault = FALSE
  ))
  model$pars$group_effect[model$pars$param == "t0m"] <- TRUE
  set.seed(20260902)
  group <- c(-1, -.5, 0, .5, 1, NA)
  dat <- do.call(rbind, lapply(seq_along(group), function(i) {
    data.frame(id = i, time = 0:3,
      Y1 = stats::rnorm(4, ifelse(is.na(group[i]), .5, group[i]), .5),
      group = group[i])
  }))

  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    optimize = TRUE, cores = 1, savescores = FALSE)))
  expect_s3_class(fit, "ctJuliaFit")
  expect_true(is.finite(fit$estimate$loglik))
  # The imputed value is what the fit conditioned on -- not a sentinel, and not
  # dropped.
  expect_false(any(fit$model_spec$tipred_data == 99999))
  expect_false(anyNA(fit$model_spec$tipred_data))
  expect_equal(as.numeric(fit$model_spec$tipred_data)[1:5], group[1:5])
})

test_that("the imputation fallback rule warns naming the fallback taken, and only when it fires", {
  model <- suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix("t0m", 1, 1), n.TIpred = 2,
    TIpredNames = c("group", "age"), tipredDefault = FALSE
  ))
  model$pars$group_effect[model$pars$param == "t0m"] <- TRUE
  set.seed(7)
  n <- 13
  dat <- data.frame(id = rep(1:n, each = 3), time = rep(0:2, n),
    Y1 = rnorm(3 * n),
    group = rep(c(rnorm(n - 1), NA), each = 3),
    age = rep(rnorm(n), each = 3))

  # `n - 1` = 12 complete cases. The "full" tier needs one coefficient per
  # other predictor plus per-subject manifest mean and SD (3 here), so
  # ~15 complete cases by the ~5-per-coefficient rule -- 12 is short, so it
  # falls back (to "predictors": 2 coefficients, ~10 needed, 12 available,
  # which succeeds) and warns naming that.
  warned <- character()
  withCallingHandlers({
    prepped <- suppressMessages(ctFit(dat, model, backend = "julia",
      optimize = FALSE, fit = FALSE, intoverpop = "augmented"))
  }, warning = function(w) { warned[[length(warned) + 1L]] <<- conditionMessage(w); invokeRestart("muffleWarning") })
  fallback_warning <- grep("not enough complete cases", warned, value = TRUE)
  expect_length(fallback_warning, 1L)
  expect_match(fallback_warning, "TI predictor 'group'", fixed = TRUE)

  # With includeOutcome=FALSE the "full" tier is never attempted, so the
  # "predictors" tier (12 complete cases for 2 coefficients, comfortably
  # over the ~5-per-coefficient rule) succeeds outright -- no warning.
  warned2 <- character()
  withCallingHandlers({
    prepped2 <- suppressMessages(ctFit(dat, model, backend = "julia",
      optimize = FALSE, fit = FALSE, intoverpop = "augmented",
      optimcontrol = list(tipredMissingIncludeOutcome = FALSE)))
  }, warning = function(w) { warned2[[length(warned2) + 1L]] <<- conditionMessage(w); invokeRestart("muffleWarning") })
  expect_length(grep("not enough complete cases", warned2, value = TRUE), 0L)
  expect_equal(prepped2$ti_imputed[, c("subject", "predictor")],
    prepped$ti_imputed[, c("subject", "predictor")])
})

# The small sampled fit the three tests below ask things of: that it samples
# the missing cell at all, and that two consumers handle it. It was
# fitted three times, identically -- same seed, data and settings -- and is
# fitted once now (`fit_cached()`, helper-julia.R).
.tipred_missing_sampled <- function() fit_cached("tipred_missing_sampled", {
  model <- .tipred_missing_model()
  set.seed(1)
  dat <- data.frame(id = rep(1:6, each = 4), time = rep(0:3, 6),
    Y1 = rnorm(24, 0, 1),
    group = rep(c(-1, 2, 0.3, -0.7, 1.5, NA), each = 4))
  suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    optimize = FALSE, intoverpop = "augmented", chains = 1, iter = 60,
    cores = 1, sampleControl = list(warmup = 30),
    optimcontrol = list(gradient = "forward"))))
})

test_that("a small julia fit actually samples a missing TI predictor value end to end", {
  skip_without_julia()
  fit <- .tipred_missing_sampled()
  expect_s3_class(fit, "ctJuliaFit")
  idx <- fit$model_spec$ti_missing$parameter
  expect_length(idx, 1L)
  draws <- fit$estimate$rawposterior[, idx]
  expect_true(all(is.finite(draws)))
  # Not asserting a tight value -- the closed-form check for the sampler's
  # correctness lives in the julia engine suite. This just confirms the R
  # side actually reaches a real posterior over the sampled cell, with
  # variation across draws (a mean-imputation bug would collapse this to a
  # single repeated number).
  expect_gt(sd(draws), 0)
})

test_that("closed form via ctFit(): posterior of an isolated missing predictor recovers its prior", {
  skip_without_julia()
  # No TI effect at all (the `group_effect` column is never set), so the
  # sampled value has no process-likelihood contribution and its marginal
  # posterior is exactly its Normal(mu, sigma) imputation conditional -- the
  # same construction as the julia engine suite's closed-form test, run here
  # through the real R -> Julia bridge. `t0m` is left free -- a model with
  # zero free parameters at all hits an unrelated, pre-existing
  # `.ctModelMatSetup()` error (`undefined columns selected`), reproduced
  # independently of backend and of this feature; one free parameter with no
  # TI effect on it keeps the raw vector at "one ordinary parameter plus one
  # sampled predictor value" without touching that.
  model <- suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix(-1.5, 1, 1),
    DIFFUSION = matrix(.3, 1, 1), MANIFESTVAR = matrix(.4, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(.5, 1, 1),
    T0MEANS = matrix("t0m", 1, 1), n.TIpred = 1, TIpredNames = "group",
    tipredDefault = FALSE
  ))
  set.seed(3)
  dat <- data.frame(id = rep(1:6, each = 4), time = rep(0:3, 6),
    Y1 = rnorm(24, 0, 1),
    group = rep(c(-2, -1, 1, 2, 0.5, NA), each = 4))
  mu <- mean(c(-2, -1, 1, 2, 0.5))
  sigma <- sd(c(-2, -1, 1, 2, 0.5))

  # On the default route, which samples the random effects jointly.
  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    optimize = FALSE, chains = 1, cores = 1,
    sampleControl = list(iter = 300L, warmup = 150L))))
  expect_identical(fit$args$resolved$intoverpop, "none")
  draws <- fit$estimate$rawposterior[, fit$model_spec$ti_missing$parameter]

  # One draw count, not two. The error shrinking with the draws -- the error at
  # 4000 draws below the error at 200 -- is the sampler's property, and the
  # julia engine suite (test_ti_missing_predictor.jl) shows it, with several
  # counts and a fixed seed. What this test adds is that the R side reaches
  # that posterior through `ctFit()`, which one count shows. It ran a second
  # chain of 1500 draws against the same two bounds, 90% of this file's time,
  # and no comparison between the two was ever asserted: one realisation per
  # count cannot order them. The bounds catch what the R side can get wrong --
  # the imputation's spread handed to the engine doubled, or its mean moved by
  # a unit, fail them at 150 draws.
  expect_lt(abs(mean(draws) - mu), 0.6)
  expect_lt(abs(var(draws) - sigma^2), sigma^2)
})

test_that("targeted stan comparison: same model, same gap, agreeing posteriors", {
  skip_without_julia()
  # The smallest model that exercises the path: one latent, one manifest, one
  # TI predictor with one missing cell, a TI effect on T0MEANS. Symmetric
  # complete-case predictor values (mean exactly 0) and tipredsimputedscale
  # set to their SD make Stan's own imputation prior -- fixed at
  # Normal(0, tipredsimputedscale) -- coincide with this feature's marginal
  # fallback conditional (single predictor, includeOutcome=FALSE), so the two
  # backends are targeting the same posterior rather than two different
  # priors by construction.
  model <- .tipred_missing_model()
  groupvals <- c(-2, -1, 1, 2, NA)
  sigma0 <- sd(groupvals[!is.na(groupvals)])
  model$tipredsimputedscale <- sigma0
  set.seed(42)
  dat <- data.frame(id = rep(1:5, each = 4), time = rep(0:3, 5),
    Y1 = rnorm(20, 0, 1), group = rep(groupvals, each = 4))

  jfit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    optimize = FALSE, intoverpop = "augmented",
    cores = 1, sampleControl = list(chains = 1, iter = 800L, warmup = 300L),
    optimcontrol = list(gradient = "forward",
      tipredMissingIncludeOutcome = FALSE))))
  jidx <- jfit$model_spec$ti_missing$parameter
  jdraws <- jfit$estimate$rawposterior[, jidx]

  sfit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "stan",
    optimize = FALSE, iter = 800L, chains = 1, cores = 1, verbose = 0)))
  sdraws <- as.numeric(rstan::extract(sfit$stanfit$stanfit)$tipredsimputed)

  # Observed on this machine: julia mean 0.150 sd 1.738 (500 draws), stan
  # mean 0.556 sd 1.452 (400 draws) -- both centred near the shared
  # Normal(0, 1.826) prior with a similar, modest pull from the process
  # likelihood, and within single-chain Monte Carlo noise of one another.
  # The tolerance below is wide enough to absorb that noise but would still
  # catch a genuinely wrong conditional: a wrong sign, a mean many SDs off,
  # or a scale wrong by a factor of two.
  expect_lt(abs(mean(jdraws) - mean(sdraws)), 1.2)
  expect_gt(sd(jdraws) / sd(sdraws), 0.5)
  expect_lt(sd(jdraws) / sd(sdraws), 2.0)
})

# Auxiliary/uncertainty functions on this same sampled-missing-TI-predictor
# fit -- the 2026-09 sampled-fit review (review/J7-sampled-fit-support.md).
test_that("per-subject scores carry a sampled TI predictor value on its own subject's row", {
  skip_without_julia()
  fit <- .tipred_missing_sampled()
  cells <- fit$model_spec$ti_missing
  expect_equal(nrow(cells), 1L)
  est <- as.numeric(fit$estimate$raw)
  scores <- ctsem:::.ctBackendScoreMatrix(fit, est)
  gradient <- as.numeric(ctJuliaEvaluate(fit, est, gradient = TRUE, cores = 1L)$gradient)
  expect_equal(colSums(scores), gradient, tolerance = 1e-8)
  expect_equal(which(scores[, cells$parameter] != 0), cells$subject)
})

test_that("ctGenerateFromFit() generates from a fit that sampled a TI predictor value", {
  skip_without_julia()
  fit <- .tipred_missing_sampled()
  expect_no_error(suppressMessages(ctGenerateFromFit(fit, nsamples = 5, cores = 1)))
})

test_that("a rebuild over some subjects takes the fit's values for its missing TI predictors", {
  model <- .tipred_missing_model()
  dat <- data.frame(id = rep(1:4, each = 3), time = rep(0:2, 4), Y1 = 0,
    group = rep(c(-1, 2, NA, .5), each = 3))
  fitted <- suppressWarnings(suppressMessages(ctFit(dat, model,
    backend = "julia", fit = FALSE)))
  # Subjects 3 and 4: one observed value, too few to impute from on its own.
  rebuild <- function() suppressMessages(ctsem:::.ctJuliaPrepare(
    dat[dat$id %in% 3:4, ], fitted$model, intoverpop = fitted$intoverpop))
  expect_true(is.na(rebuild()$ti_imputed$mu))

  carried <- ctsem:::.ctJuliaCarryTIMissing(rebuild(), fitted)
  expect_equal(carried$tipred_data[, 1], c(fitted$ti_imputed$mu, .5))
  expect_null(carried$ti_missing)

  # From a fit that sampled it: the same raw parameter, so the fit's raw
  # vectors apply to the rebuild unchanged.
  sampled <- ctsem:::.ctJuliaTIMissingSampled(fitted)
  carried <- ctsem:::.ctJuliaCarryTIMissing(rebuild(), sampled)
  expect_null(carried$ti_imputed)
  expect_equal(carried$ti_missing$parameter, sampled$ti_missing$parameter)
  expect_equal(carried$ti_missing$subject, 1L)
  expect_equal(carried$tipred_data[1, 1], 99999)

  # A subject the fit never had, with nothing to impute from, is refused.
  stranger <- dat[dat$id %in% 3:4, ]
  stranger$id <- stranger$id + 10
  expect_error(ctsem:::.ctJuliaCarryTIMissing(suppressMessages(
    ctsem:::.ctJuliaPrepare(stranger, fitted$model, intoverpop = fitted$intoverpop)),
    fitted), "fewer than two subjects")
})
