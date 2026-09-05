# Starting values read off the data (R/ctDataStart.R).
#
# The load-bearing property is not that any particular number comes out, but
# that the number tracks the data's scale -- the old fixed start did not, and
# that is what let a badly scaled fit collapse the latent to zero variance.
# Most of this is therefore pure R against a known scale, with one Julia case
# for the failure that motivated it.

.datastart_gen <- function(scale = 1, nsubjects = 40, nobs = 10, seed = 3) {
  set.seed(seed)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 3,
    manifestNames = c("y1", "y2", "y3"), latentNames = "eta1",
    LAMBDA = matrix(c(1, 0.8, 1.2), 3, 1), DRIFT = matrix(-0.4),
    DIFFUSION = matrix(1.5), T0VAR = matrix(1.7), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(c(0, 0.5, -0.5), 3, 1),
    MANIFESTVAR = diag(0.8, 3), Tpoints = nobs))
  d <- data.frame(ctGenerate(gen, n.subjects = nsubjects, Tpoints = nobs,
    backend = "r"))
  d[, c("y1", "y2", "y3")] <- d[, c("y1", "y2", "y3")] * scale
  d
}

.datastart_model <- function() {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 3,
    manifestNames = c("y1", "y2", "y3"), latentNames = "eta1",
    LAMBDA = matrix(c(1, "lambda_y2", "lambda_y3"), 3, 1),
    CINT = matrix(0), T0MEANS = matrix(0)))
  m$pars$indvarying <- FALSE
  m
}

# The raw vector back on the parameter scale, so a test can say what the fit
# would actually start at rather than what its unconstrained coordinate is.
.datastart_values <- function(d, m) {
  spec <- ctsem:::.ctJuliaPrepare(d, m, priors = FALSE, intoverpop = "augmented")
  pt <- spec$parameter_table
  npar <- max(pt$parnumber, na.rm = TRUE)
  raw <- ctsem:::.ctDataStart(d, m, spec, npar)
  if (is.null(raw)) return(NULL)
  rows <- pt[!is.na(pt$parnumber) & !is.na(pt$transform) &
      pt$matrix %in% c("DRIFT", "DIFFUSION", "T0VAR", "MANIFESTVAR"), ,
    drop = FALSE]
  rows <- rows[!duplicated(rows$parnumber), , drop = FALSE]
  out <- stats::setNames(rep(NA_real_, nrow(rows)),
    paste0(rows$matrix, "[", rows$row, ",", rows$col, "]"))
  log1p_exp <- function(x) ifelse(x > 30, x, log1p(exp(x)))
  for (i in seq_len(nrow(rows))) {
    p <- raw[rows$parnumber[i]]
    if (!is.finite(p)) next
    text <- gsub("param[[][0-9]+[]]", "param", rows$transform[i])
    out[i] <- eval(parse(text = text), list(param = p, log1p_exp = log1p_exp))
  }
  out
}

test_that("the derived start follows the data's scale", {
  # The generating DIFFUSION is 1.5 * scale and the measurement sd 0.8 * scale.
  # The claim is only that the start lands in the right decade -- it is a
  # starting value, not an estimate.
  for (scale in c(0.01, 1, 100)) {
    v <- .datastart_values(.datastart_gen(scale), .datastart_model())
    expect_false(is.null(v), info = scale)
    expect_true(v[["DIFFUSION[1,1]"]] > 0.1 * 1.5 * scale, info = scale)
    expect_true(v[["DIFFUSION[1,1]"]] < 100 * 1.5 * scale, info = scale)
    expect_true(v[["MANIFESTVAR[1,1]"]] > 0.1 * 0.8 * scale, info = scale)
    expect_true(v[["MANIFESTVAR[1,1]"]] < 100 * 0.8 * scale, info = scale)
    # Mean reverting rather than a random walk or white noise.
    expect_true(v[["DRIFT[1,1]"]] < 0, info = scale)
    expect_true(v[["DRIFT[1,1]"]] > -3.01, info = scale)
  }
})

test_that("the start moves with the data where the fixed default cannot", {
  small <- .datastart_values(.datastart_gen(0.01), .datastart_model())
  large <- .datastart_values(.datastart_gen(100), .datastart_model())
  # Four orders of magnitude of data has to show up as orders of magnitude of
  # starting value; the old default was one number for both.
  expect_gt(large[["DIFFUSION[1,1]"]] / small[["DIFFUSION[1,1]"]], 100)
  expect_gt(large[["MANIFESTVAR[1,1]"]] / small[["MANIFESTVAR[1,1]"]], 100)
})

test_that("nothing derived escapes the clamps", {
  # A deliberately absurd scale: the answer must be bounded, not enormous.
  v <- .datastart_values(.datastart_gen(1e6), .datastart_model())
  expect_true(all(is.finite(v[!is.na(v)])))
  expect_lt(max(v[!is.na(v)]), 1e4)
  expect_true(v[["DRIFT[1,1]"]] >= -3.01 && v[["DRIFT[1,1]"]] < 0)
})

test_that("a count is read on the log rate, not on the rate", {
  # The distinction that matters: a Poisson indicator's observations are on the
  # rate scale while its latent is on the log rate, so sd(y) overstates the
  # latent scale badly -- measured, by a factor of 6.9 on a process of sd 1.13,
  # in the one direction that makes fits fail.
  y <- stats::rpois(500, exp(1.4 + 0.5 * stats::rnorm(500)))
  expect_equal(ctsem:::.ctDataStartManifestScale(y, 3L),
    stats::sd(log1p(y)))
  expect_false(isTRUE(all.equal(ctsem:::.ctDataStartManifestScale(y, 3L),
    stats::sd(y))))
  # Gaussian and censored are read directly; binary and ordinal say nothing,
  # because their link fixes the latent scale rather than the data.
  z <- stats::rnorm(500, 0, 2)
  expect_equal(ctsem:::.ctDataStartManifestScale(z, 0L), stats::sd(z))
  expect_equal(ctsem:::.ctDataStartManifestScale(z, 4L), stats::sd(z))
  expect_true(is.na(ctsem:::.ctDataStartManifestScale(z, 1L)))
  expect_true(is.na(ctsem:::.ctDataStartManifestScale(z, 2L)))
})

test_that("a process measured only through a logit gets the link's scale", {
  set.seed(4)
  d <- .datastart_gen(1)
  d$b1 <- stats::rbinom(nrow(d), 1, 0.5)
  d$b2 <- stats::rbinom(nrow(d), 1, 0.5)
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 2, manifestNames = c("b1", "b2"), latentNames = "eta1",
    manifesttype = c(1L, 1L), LAMBDA = matrix(1, 2, 1),
    MANIFESTMEANS = matrix(c("mm_b1", "mm_b2"), 2, 1),
    MANIFESTVAR = diag(0, 2), CINT = matrix(0), T0MEANS = matrix(0))))
  m$pars$indvarying <- FALSE
  v <- .datastart_values(d, m)
  # sd(y) for a 0/1 variable is at most 0.5 and means nothing about the latent;
  # the constant used instead is on the logit's own scale.
  expect_true(v[["T0VAR[1,1]"]] > 0.5)
  expect_true(v[["T0VAR[1,1]"]] < 5)
})

test_that("the discrete time diffusion is scaled by the autoregression", {
  d <- .datastart_gen(1)
  mdt <- suppressWarnings(suppressMessages(ctModel(type = "dt", n.latent = 1,
    n.manifest = 3, manifestNames = c("y1", "y2", "y3"), latentNames = "eta1",
    LAMBDA = matrix(c(1, "lambda_y2", "lambda_y3"), 3, 1),
    CINT = matrix(0), T0MEANS = matrix(0))))
  mdt$pars$indvarying <- FALSE
  v <- .datastart_values(d, mdt)
  expect_false(is.null(v))
  # A discrete time DRIFT diagonal is an autoregression in (0,1), not a
  # continuous time rate, and the innovation is scaled by sqrt(1-phi^2)
  # rather than by sqrt(2|a|) -- so it must come out below the process spread.
  expect_true(v[["DRIFT[1,1]"]] > 0 && v[["DRIFT[1,1]"]] < 1)
  expect_true(v[["DIFFUSION[1,1]"]] < v[["T0VAR[1,1]"]])
})

test_that("an unsolvable transform leaves the default rather than guessing", {
  # A target no transform can reach, and text that is not a transform at all.
  expect_true(is.na(ctsem:::.ctDataStartInvert(NA_character_, 1)))
  expect_true(is.na(ctsem:::.ctDataStartInvert("not a transform (", 1)))
  expect_true(is.na(ctsem:::.ctDataStartInvert("10 * log1p_exp(2 * param[1])",
    NA_real_)))
  # Reachable targets are solved to the value asked for.
  raw <- ctsem:::.ctDataStartInvert("10 * log1p_exp(2 * param[1])", 1.5)
  expect_equal(10 * log1p(exp(2 * raw)), 1.5, tolerance = 1e-6)
  # Out of reach within the bound clamps rather than failing.
  expect_equal(abs(ctsem:::.ctDataStartInvert("10 * log1p_exp(2 * param[1])",
    1e6)), 3)
})

test_that("the derived start rescues a badly scaled fit", {
  skip_without_julia()
  d <- .datastart_gen(0.01, nsubjects = 50)
  fit <- function(datastart) suppressWarnings(suppressMessages(
    ctFit(d, .datastart_model(), backend = "julia", cores = 2,
      optimcontrol = list(estonly = TRUE, datastart = datastart))))
  on_ <- fit(TRUE)
  off <- fit(FALSE)
  # Measured: 4565.60 against 4237.21, with the loadings running to -844 and
  # -1242 in the failing fit against generating values of 0.8 and 1.2.
  expect_true(isTRUE(on_$estimate$converged))
  expect_gt(as.numeric(on_$estimate$loglik), as.numeric(off$estimate$loglik) - 1)
  expect_lt(max(abs(on_$estimate$raw)), 10)
})
