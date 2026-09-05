# `ctIdentify()` answers a question nothing else in the package does: whether
# the data can inform the free parameters, before a fit is spent finding out.
# The two cases below are the ones that matter -- a model that is identified and
# one that is not, differing by a single freed loading.

.identify_data <- function(nsubjects = 40, nobs = 8) {
  set.seed(5)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6), MANIFESTVAR = matrix(0.3),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = nobs))
  ctGenerate(gen, n.subjects = nsubjects, Tpoints = nobs, backend = "r")
}

.identify_model <- function(lambda = 1) {
  suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = if (is.character(lambda)) matrix(lambda) else matrix(lambda),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(0)))
}

test_that("an identified model reports no uninformed directions", {
  skip_without_julia()
  result <- suppressMessages(ctIdentify(.identify_data(), .identify_model(),
    cores = 2))
  expect_s3_class(result, "ctIdentify")
  expect_length(result$parameters, 0L)
  expect_equal(result$nsubjects, 40L)
  expect_false(result$rankLimited)
  # The margin matters as much as the verdict: a threshold that only just
  # cleared would be a threshold tuned to this example.
  expect_gt(result$smallest, 1e-13)
})

test_that("a freed loading is reported with the parameter it trades against", {
  skip_without_julia()
  # LAMBDA free as well as DIFFUSION leaves only their product determined, so
  # the pair is unidentified however much data there is.
  result <- suppressMessages(ctIdentify(.identify_data(),
    .identify_model("lambda"), cores = 2))
  expect_gte(result$nweak, 1L)
  # Naming the parameters is the point -- a verdict alone would not tell anyone
  # which cell of which matrix to change.
  expect_true(all(c("lambda", "diff_eta1") %in% result$parameters))
  # Structurally flat means flat to machine precision, not merely small.
  expect_lt(result$smallest, 1e-13)
})

test_that("the information used is positive semi-definite away from any mode", {
  skip_without_julia()
  # The reason this uses scores rather than the Hessian: at an arbitrary point
  # the Hessian is indefinite and its flat directions describe the point, not
  # the data. The score information cannot be.
  spec <- .ctFitJuliaBackend(.identify_data(), .identify_model(), fit = FALSE,
    priors = FALSE, intoverpop = "augmented", cores = 1, verbose = 0)
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  set.seed(2)
  info <- ctsem:::.ctIdentifyInformation(spec, stats::rnorm(npar, 0, 1))
  values <- eigen(info$information, symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(values), -sqrt(.Machine$double.eps) * max(values))
  # Scaled to unit diagonal, so a threshold reads information rather than the
  # units each raw parameter happens to be on.
  expect_equal(diag(info$information), rep(1, npar), tolerance = 1e-8)
})

test_that("a model with no free parameters is refused rather than answered", {
  skip_without_julia()
  fixed <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6), MANIFESTVAR = matrix(0.3),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0)))
  expect_error(suppressMessages(ctIdentify(.identify_data(), fixed)),
    "no free parameters")
})
