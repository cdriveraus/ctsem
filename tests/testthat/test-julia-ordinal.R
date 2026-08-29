# Ordinal observations on the julia backend.
#
# The measurement update integrates the cumulative-logit likelihood against the
# predicted state, exactly as the binary path integrates the Bernoulli one --
# the only difference is the likelihood inside the node loop. So the tests that
# matter are the same ones: a likelihood that can be written down, a gradient
# finite differences can confirm, and the fact that the binary path still gives
# its own answer through the shared code.

.jord_thresholds <- c(-1.0, 0.4, 1.9)   # four categories

.jord_draw <- function(eta, tau) {
  invlog <- function(x) 1 / (1 + exp(-x))
  cumulative <- sapply(seq_along(tau), function(k) invlog(tau[k] - eta))
  p <- cbind(cumulative, 1)
  p <- cbind(p[, 1, drop = FALSE], t(apply(p, 1, diff)))
  apply(p, 1, function(pr) sample.int(length(pr), 1, prob = pmax(pr, 0)))
}

.jord_data <- function(nsubjects = 40, nobs = 12, nindicators = 2, seed = 7) {
  set.seed(seed)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(0.001),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = nobs))
  d <- data.frame(ctGenerate(gen, n.subjects = nsubjects, Tpoints = nobs,
    backend = "r"))
  for (i in seq_len(nindicators)) {
    d[[paste0("o", i)]] <- .jord_draw(d$eta, .jord_thresholds)
  }
  d$eta <- NULL
  d
}

.jord_model <- function(nindicators = 2, ncategories = 4) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = nindicators, manifestNames = paste0("o", 1:nindicators),
    latentNames = "eta1", manifesttype = rep(2L, nindicators),
    ncategories = rep(ncategories, nindicators),
    LAMBDA = matrix(1, nindicators, 1),
    MANIFESTMEANS = matrix(0, nindicators, 1), CINT = matrix(0),
    T0MEANS = matrix(0), MANIFESTVAR = diag(0, nindicators))))
  m$pars$indvarying <- FALSE
  m
}

.jord_spec <- function(d, m) {
  spec <- ctsem:::.ctJuliaPrepare(d, m, priors = FALSE,
    intoverpop = "augmented")
  structure(spec, class = c("ctJuliaModel", "ctFitModel"))
}

test_that("ctModel refuses an ordinal specification it cannot act on", {
  base <- list(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "o1", latentNames = "eta1", LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(0), CINT = matrix(0), T0MEANS = matrix(0),
    MANIFESTVAR = matrix(0))
  build <- function(...) suppressWarnings(suppressMessages(
    do.call(ctModel, c(base, list(...)))))
  # The category count cannot be inferred: ctModel never sees the data, and
  # ctGenerate needs it as much as ctFit does.
  expect_error(build(manifesttype = 2L), "ncategories")
  expect_error(build(manifesttype = 2L, ncategories = 2L), "at least 3")
  expect_error(build(manifesttype = 3L), "must be 0")
  m <- build(manifesttype = 2L, ncategories = 4L)
  # ctModel returns the converted model, so the matrix is in `pars` rather
  # than sitting on the object under its own name.
  thr <- m$pars[m$pars$matrix %in% "THRESHOLDS", ]
  expect_equal(nrow(thr), 3L)
  expect_equal(max(thr$col), 3L)
  expect_equal(m$ncategories, 4L)
})

test_that("THRESHOLDS holds a free first threshold and positive gaps", {
  m <- .jord_model(nindicators = 1, ncategories = 4)
  thr <- m$pars[m$pars$matrix %in% "THRESHOLDS", ]
  expect_equal(nrow(thr), 3L)
  expect_true(all(is.na(thr$value)))
  # Column 1 is an unconstrained threshold; the rest are gaps, and a gap that
  # could go negative would let the thresholds cross.
  expect_false(grepl("log1p_exp", thr$transform[thr$col == 1]))
  expect_true(all(grepl("log1p_exp", thr$transform[thr$col > 1])))
})

test_that("a ragged model leaves the unused threshold cells fixed", {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("o1", "o2"), latentNames = "eta1",
    manifesttype = c(2L, 2L), ncategories = c(4L, 3L),
    LAMBDA = matrix(1, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix(0), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 2))))
  thr <- m$pars[m$pars$matrix %in% "THRESHOLDS", ]
  # o1 needs three thresholds and o2 two, so the matrix is 2x3 with one cell
  # unused -- fixed, and never read by the engine.
  expect_equal(sum(is.na(thr$value)), 5L)
  spare <- thr[thr$row == 2 & thr$col == 3, ]
  expect_equal(as.numeric(spare$value), 0)
})

test_that("one ordinal observation has the likelihood it can be shown to have", {
  skip_on_cran()
  skip_without_julia()
  # T0MEANS 0, T0VAR 1, one observation in the bottom category of a
  # three-category variable with both thresholds fixed at zero. The bottom
  # category is then P(y=1|eta) = inv_logit(0 - eta), and its marginal is
  # int inv_logit(-eta) N(eta;0,1) deta = 0.5 by symmetry.
  d <- data.frame(id = 1L, time = 0, o1 = 1L, x = NA_real_)
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("o1", "x"), latentNames = "eta1",
    manifesttype = c(2L, 0L), ncategories = c(3L, 0L),
    LAMBDA = matrix(c(1, 0), 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix(0), T0MEANS = matrix(0), DRIFT = matrix(-0.3),
    DIFFUSION = matrix("dif"), T0VAR = matrix(1),
    MANIFESTVAR = diag(c(0, 1), 2))))
  m$pars$indvarying <- FALSE
  m$pars$value[m$pars$matrix %in% "THRESHOLDS"] <- 0
  m$pars$param[m$pars$matrix %in% "THRESHOLDS"] <- NA
  m$pars$transform[m$pars$matrix %in% "THRESHOLDS"] <- NA
  handle <- .jord_spec(d, m)
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  value <- ctJuliaEvaluate(handle, rep(0, npar), gradient = FALSE)$value
  expect_equal(as.numeric(value), log(0.5), tolerance = 1e-7)
})

test_that("the gradient matches finite differences, thresholds included", {
  skip_on_cran()
  skip_without_julia()
  # The reverse pass is hand written and the threshold cotangent has to travel
  # back through the cumulative sum as well as the quadrature; a reverse pass
  # that misses a path returns a confident wrong number rather than failing.
  handle <- .jord_spec(.jord_data(nsubjects = 20), .jord_model())
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  set.seed(4)
  at <- stats::rnorm(npar, 0, 0.4)
  adjoint <- ctJuliaEvaluate(handle, at, gradient = TRUE)$gradient
  eps <- 1e-6
  fd <- vapply(seq_len(npar), function(i) {
    up <- at; up[i] <- up[i] + eps
    dn <- at; dn[i] <- dn[i] - eps
    (ctJuliaEvaluate(handle, up)$value - ctJuliaEvaluate(handle, dn)$value) /
      (2 * eps)
  }, numeric(1))
  expect_equal(as.numeric(adjoint), fd, tolerance = 1e-5)

  # Every threshold must actually be in the gradient: a cotangent silently
  # dropped would leave a zero here and an unmoved parameter in a fit.
  tab <- handle$parameter_table
  thresholds <- unique(tab$parnumber[tab$matrix %in% "THRESHOLDS" &
      !is.na(tab$parnumber)])
  expect_gt(length(thresholds), 0L)
  expect_true(all(abs(adjoint[thresholds]) > 1e-8))
})

test_that("the ordinal adjoint is exact with more than one latent state", {
  skip_on_cran()
  skip_without_julia()
  # The companion of the binary case in test-julia-binary.R, and it found the
  # same bug harder: with a random effect expanding the state, the adjoint was
  # out by 165% before the covariance cotangent stopped being assumed
  # symmetric. Forward mode rather than a finite difference, because it is
  # exact and so leaves nowhere for the disagreement to hide.
  d <- .jord_data(nsubjects = 15, nobs = 12, nindicators = 2)
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("o1", "o2"), latentNames = "eta1",
    manifesttype = c(2L, 2L), ncategories = c(4L, 4L),
    LAMBDA = matrix(1, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 2))))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% "cint"] <- TRUE
  handle <- .jord_spec(d, m)
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  set.seed(3)
  at <- stats::rnorm(npar, 0, 0.3)
  adjoint <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
    gradient_method = "adjoint")$gradient)
  forward <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
    gradient_method = "forward")$gradient)
  expect_equal(adjoint, forward, tolerance = 1e-9)
})

test_that("a fit recovers the thresholds it generated from", {
  skip_on_cran()
  skip_without_julia()
  d <- .jord_data(nsubjects = 120, nobs = 15, nindicators = 3, seed = 5)
  m <- .jord_model(nindicators = 3)
  fit <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    cores = 1, optimcontrol = list(estonly = TRUE))))
  est <- summary(fit)$popmeans
  # THRESHOLDS holds the first threshold then gaps, so that is what the truth
  # has to be expressed as too.
  gaps <- c(.jord_thresholds[1], diff(.jord_thresholds))
  for (i in 1:3) {
    for (k in 1:3) {
      name <- paste0("threshold_o", i, "_", k)
      expect_true(name %in% rownames(est))
      expect_equal(unname(est[name, "mean"]), gaps[k], tolerance = 0.25)
    }
  }
  expect_equal(unname(est["drift_eta1", "mean"]), -0.3, tolerance = 0.15)
  expect_equal(unname(est["diff_eta1", "mean"]), 0.8, tolerance = 0.2)
})

test_that("generation draws categories in the proportions the model implies", {
  skip_on_cran()
  skip_without_julia()
  # Every parameter fixed. With free ones the values come from a prior draw,
  # which can put the latent variance so high that only the end categories are
  # ever reached -- true of the model that was drawn, but no test of the code.
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "o1", latentNames = "eta1",
    manifesttype = 2L, ncategories = 4L, LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(0), CINT = matrix(0), T0MEANS = matrix(0),
    MANIFESTVAR = matrix(0), DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8),
    T0VAR = matrix(1))))
  m$pars$indvarying <- FALSE
  gaps <- c(.jord_thresholds[1], diff(.jord_thresholds))
  sel <- m$pars$matrix %in% "THRESHOLDS"
  m$pars$value[sel] <- gaps[m$pars$col[sel]]
  m$pars$param[sel] <- NA
  m$pars$transform[sel] <- NA

  set.seed(1)
  g <- suppressWarnings(suppressMessages(
    ctGenerate(m, n.subjects = 300, Tpoints = 12, backend = "julia")))
  values <- g[, "o1"]
  expect_true(all(values == round(values)))
  expect_true(all(values >= 1 & values <= 4))
  expect_equal(sort(unique(values)), c(1, 2, 3, 4))

  # The same latent process, categorised in R, says what the proportions
  # should be -- a generator that ignored the thresholds would still pass the
  # checks above but not this one.
  set.seed(2)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(1e-8),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = 12))
  eta <- data.frame(ctGenerate(gen, n.subjects = 300, Tpoints = 12,
    backend = "r"))$eta
  invlog <- function(x) 1 / (1 + exp(-x))
  cum <- sapply(.jord_thresholds, function(t) invlog(t - eta))
  p <- cbind(cum, 1)
  p <- cbind(p[, 1, drop = FALSE], t(apply(p, 1, diff)))
  expect_equal(as.numeric(prop.table(table(values))), colMeans(p),
    tolerance = 0.05)
})

test_that("stan refuses an ordinal model rather than ignoring the thresholds", {
  skip_on_cran()
  d <- .jord_data(nsubjects = 10, nobs = 5, nindicators = 1)
  m <- .jord_model(nindicators = 1)
  expect_error(suppressWarnings(suppressMessages(
    ctFit(d, m, backend = "stan", cores = 1,
      optimcontrol = list(estonly = TRUE)))), "julia")
})

test_that("ordinal data are checked against the model's category count", {
  skip_on_cran()
  d <- .jord_data(nsubjects = 10, nobs = 5, nindicators = 1)
  m <- .jord_model(nindicators = 1, ncategories = 4)
  bad <- d
  bad$o1[1] <- 9
  expect_error(suppressWarnings(suppressMessages(ctFit(bad, m,
    backend = "julia", cores = 1, optimcontrol = list(estonly = TRUE)))),
    "ncategories")
  zerobased <- d
  zerobased$o1 <- zerobased$o1 - 1L
  expect_error(suppressWarnings(suppressMessages(ctFit(zerobased, m,
    backend = "julia", cores = 1, optimcontrol = list(estonly = TRUE)))),
    "consecutive integers")
})

test_that("mixed ordinal, binary and Gaussian indicators fit together", {
  skip_on_cran()
  skip_without_julia()
  set.seed(6)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(0.001),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = 12))
  d <- data.frame(ctGenerate(gen, n.subjects = 50, Tpoints = 12,
    backend = "r"))
  d$o1 <- .jord_draw(d$eta, .jord_thresholds)
  d$b1 <- stats::rbinom(nrow(d), 1, 1 / (1 + exp(-d$eta)))
  d$y1 <- d$eta + stats::rnorm(nrow(d), 0, 0.5)
  d$eta <- NULL
  mvar <- diag(0, 3); mvar[3, 3] <- "mvar"
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 3, manifestNames = c("o1", "b1", "y1"), latentNames = "eta1",
    manifesttype = c(2L, 1L, 0L), ncategories = c(4L, 0L, 0L),
    LAMBDA = matrix(1, 3, 1), MANIFESTMEANS = matrix(0, 3, 1),
    CINT = matrix(0), T0MEANS = matrix(0), MANIFESTVAR = mvar)))
  m$pars$indvarying <- FALSE
  fit <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    cores = 1, optimcontrol = list(estonly = TRUE))))
  expect_true(is.finite(fit$estimate$loglik))
  est <- summary(fit)$popmeans
  expect_true("threshold_o1_1" %in% rownames(est))
  expect_true("mvar" %in% rownames(est))
  expect_equal(unname(est["drift_eta1", "mean"]), -0.3, tolerance = 0.2)
})

test_that("post-fit functions work on an ordinal fit", {
  skip_on_cran()
  skip_without_julia()
  d <- .jord_data(nsubjects = 30, nobs = 8, nindicators = 2)
  fit <- suppressWarnings(suppressMessages(ctFit(d, .jord_model(),
    backend = "julia", cores = 1, optimcontrol = list(estonly = TRUE))))
  expect_s3_class(ctKalman(fit, subjects = 1), "data.frame")
  expect_s3_class(ctPredict(fit, subjects = 1), "data.frame")
  expect_error(ctACFresiduals(fit), NA)
  expect_error(summary(fit), NA)
})

test_that("ctModelLatex renders the ordinal link", {
  m <- .jord_model(nindicators = 2)
  tex <- paste(as.character(ctModelLatex(m, compile = FALSE)), collapse = "\n")
  expect_true(grepl("operatorname{logit}", tex, fixed = TRUE))
  expect_true(grepl("tau_{1,k}", tex, fixed = TRUE))
  # An ordinal observation is integrated, not given a Gaussian error, and the
  # rendering has to say so rather than showing a covariance that is not used.
  expect_true(grepl("carry no measurement error", tex, fixed = TRUE))
})

test_that("the print method names the ordinal variables and their categories", {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("o1", "b1"), latentNames = "eta1",
    manifesttype = c(2L, 1L), ncategories = c(5L, 0L),
    LAMBDA = matrix(1, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix(0), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 2))))
  out <- paste(utils::capture.output(print(m)), collapse = "\n")
  expect_true(grepl("o1 (ordinal, 5 categories)", out, fixed = TRUE))
  expect_true(grepl("b1 (binary)", out, fixed = TRUE))
})

test_that("a degenerate predicted variance costs its likelihood, not nothing", {
  skip_on_cran()
  skip_without_julia()
  # T0VAR fixed at zero makes the first occasion's linear predictor known
  # exactly, so the observation there contributes `log P(y | eta)` and moves
  # nothing. Returning zero instead -- as a likelihood of one -- is not a
  # harmless edge case: it makes a collapsing variance *pay*, so an optimiser
  # with a free T0VAR is rewarded for driving it to zero and the objective is
  # discontinuous at the point it is driving towards.
  #
  # One subject, one occasion, everything fixed, so the whole log likelihood is
  # that single term and can be written down.
  d <- data.frame(id = 1L, time = 0, o1 = 2L)
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "o1", latentNames = "eta1",
    manifesttype = 2L, ncategories = 4L, LAMBDA = matrix(1),
    MANIFESTMEANS = matrix(0), CINT = matrix(0), T0MEANS = matrix(0),
    MANIFESTVAR = matrix(0), DRIFT = matrix(-0.3), DIFFUSION = matrix(0.8),
    T0VAR = matrix(0))))
  m$pars$indvarying <- FALSE
  gaps <- c(.jord_thresholds[1], diff(.jord_thresholds))
  sel <- m$pars$matrix %in% "THRESHOLDS"
  m$pars$value[sel] <- gaps[m$pars$col[sel]]
  m$pars$param[sel] <- NA
  m$pars$transform[sel] <- NA

  handle <- .jord_spec(d, m)
  value <- as.numeric(ctJuliaEvaluate(handle, numeric(0))$value)
  # eta is exactly zero, so category 2 has probability
  # plogis(tau2) - plogis(tau1).
  expected <- log(stats::plogis(.jord_thresholds[2]) -
      stats::plogis(.jord_thresholds[1]))
  expect_equal(value, expected, tolerance = 1e-10)
  expect_lt(value, 0)
})

test_that("the gradient stays finite as a threshold gap closes", {
  skip_on_cran()
  skip_without_julia()
  # The interval probability is computed as a product of tails rather than a
  # difference of CDFs, so a gap far below the point where the difference would
  # lose its digits still gives a usable score. A gap is a free parameter under
  # a positive transform, so an optimiser reaches these values by itself.
  d <- .jord_data(nsubjects = 15, nobs = 8, nindicators = 1)
  m <- .jord_model(nindicators = 1)
  handle <- .jord_spec(d, m)
  tab <- handle$parameter_table
  npar <- max(tab$parnumber, na.rm = TRUE)
  gappar <- unique(tab$parnumber[tab$matrix %in% "THRESHOLDS" & tab$col > 1 &
      !is.na(tab$parnumber)])
  expect_gt(length(gappar), 0L)
  # Down to -15, which is a gap of about 1e-13. Below that the *cumulated*
  # thresholds stop resolving the gap at all -- `tau + 1e-35` is `tau` -- and
  # the category between two numerically identical thresholds genuinely has
  # probability zero. `_CTSEM_SATURATION` already reports a raw magnitude of 20
  # as not converged, so the optimiser is told about that region rather than
  # being expected to work in it.
  for (raw in c(-5, -10, -15)) {
    at <- rep(0.1, npar)
    at[gappar] <- raw          # gaps of roughly 9e-5 down to 2e-13
    got <- ctJuliaEvaluate(handle, at, gradient = TRUE)
    expect_true(is.finite(as.numeric(got$value)), info = paste("raw", raw))
    expect_true(all(is.finite(as.numeric(got$gradient))),
      info = paste("raw", raw))
  }
})

test_that("an unlikely observation is a large penalty, not an impossible row", {
  skip_on_cran()
  skip_without_julia()
  # Push the first threshold to about 300, so every observed category above the
  # first has a probability around exp(-300). Computed as a probability that
  # underflows to zero, the row becomes impossible, the subject's likelihood is
  # -Inf and the whole trial point is invalid -- which is what a Laplace inner
  # mode solve then has nothing to work with. Computed as a log likelihood it is
  # merely a large negative number, and the optimiser can walk back out.
  d <- .jord_data(nsubjects = 10, nobs = 8, nindicators = 1)
  m <- .jord_model(nindicators = 1)
  handle <- .jord_spec(d, m)
  tab <- handle$parameter_table
  npar <- max(tab$parnumber, na.rm = TRUE)
  first <- unique(tab$parnumber[tab$matrix %in% "THRESHOLDS" & tab$col == 1 &
      !is.na(tab$parnumber)])
  expect_equal(length(first), 1L)
  for (raw in c(5, 30, 70)) {
    at <- rep(0.1, npar)
    at[first] <- raw
    got <- ctJuliaEvaluate(handle, at, gradient = TRUE)
    expect_true(is.finite(as.numeric(got$value)), info = paste("raw", raw))
    expect_true(all(is.finite(as.numeric(got$gradient))),
      info = paste("raw", raw))
  }
  # And it really is a penalty: a threshold that far out must be much worse
  # than a sane one, or the test is passing on an answer that ignores the data.
  sane <- rep(0.1, npar); sane[first] <- 0
  far <- rep(0.1, npar); far[first] <- 30
  expect_lt(as.numeric(ctJuliaEvaluate(handle, far)$value),
    as.numeric(ctJuliaEvaluate(handle, sane)$value) - 100)
})

test_that("the adjoint is exact where the predicted variance is degenerate", {
  skip_on_cran()
  skip_without_julia()
  # T0VAR at zero makes the first occasion's linear predictor known exactly.
  # The forward update is then skipped and the observation contributes
  # `log P(y | eta)`, which still depends on the state, LAMBDA, MANIFESTMEANS
  # and the thresholds -- so the reverse pass has to contribute those terms
  # rather than skip the observation. It skipped it, which is a silently
  # missing piece of gradient exactly where a variance is collapsing.
  d <- .jord_data(nsubjects = 20, nobs = 8, nindicators = 2)
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("o1", "o2"), latentNames = "eta1",
    manifesttype = c(2L, 2L), ncategories = c(4L, 4L),
    LAMBDA = matrix(1, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix(0), T0MEANS = matrix(0), MANIFESTVAR = diag(0, 2),
    T0VAR = matrix(0))))
  m$pars$indvarying <- FALSE
  handle <- .jord_spec(d, m)
  npar <- max(handle$parameter_table$parnumber, na.rm = TRUE)
  set.seed(11)
  at <- stats::rnorm(npar, 0, 0.3)
  adjoint <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
    gradient_method = "adjoint")$gradient)
  forward <- as.numeric(ctJuliaEvaluate(handle, at, gradient = TRUE,
    gradient_method = "forward")$gradient)
  expect_true(all(is.finite(adjoint)))
  expect_equal(adjoint, forward, tolerance = 1e-9)
  # The thresholds must carry gradient here too, or the degenerate branch has
  # dropped exactly the term this test exists for.
  tab <- handle$parameter_table
  thresholds <- unique(tab$parnumber[tab$matrix %in% "THRESHOLDS" &
      !is.na(tab$parnumber)])
  expect_true(all(abs(adjoint[thresholds]) > 1e-8))
})
