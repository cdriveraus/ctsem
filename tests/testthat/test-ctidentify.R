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

# Random effects on a variance cell (DIFFUSION) and on DRIFT. Under
# `intoverpop='augmented'` the variance cell's effect is a state the
# observation mean function does not see, so the Kalman update never moves it
# and it is learned about only through its correlation with the drift effect:
# partially identified, not unidentified.
.identify_varying_model <- function() {
  model <- .identify_model()
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$matrix %in% c("DIFFUSION", "DRIFT")] <- TRUE
  model
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

test_that("two identical calls give identical output", {
  skip_without_julia()
  # The load-bearing one. Evaluation points 2 onwards were drawn from the
  # session stream, so the same call twice reported `always = (empty)` and
  # `always = popsd_df11` -- the reproducibility claim covered only the first
  # point. Different seeds set before each call, deliberately: if either
  # result depended on the session stream this would catch it.
  data <- .identify_data()
  model <- .identify_varying_model()
  set.seed(11)
  first <- suppressMessages(ctIdentify(data, model, cores = 1))
  set.seed(77)
  second <- suppressMessages(ctIdentify(data, model, cores = 1))
  expect_identical(first$parameters, second$parameters)
  expect_identical(first$partial, second$partial)
  expect_identical(first$sometimes, second$sometimes)
  # Every evaluation point, not just the reported summary.
  expect_identical(lapply(first$starts, `[[`, "at"),
    lapply(second$starts, `[[`, "at"))
  expect_identical(first, second)

  # And it must not move the caller's stream either, or a script that fits
  # after checking gets different starting values than one that does not.
  set.seed(5)
  before <- stats::rnorm(3)
  set.seed(5)
  invisible(suppressMessages(ctIdentify(data, model, cores = 1)))
  expect_equal(stats::rnorm(3), before)
})

test_that("a random effect on a variance cell reports the covariance as determined", {
  skip_without_julia()
  result <- suppressMessages(ctIdentify(.identify_data(),
    .identify_varying_model(), cores = 1))
  expect_gte(result$nweak, 1L)
  # The population sd of the diffusion effect and its correlation with the
  # drift effect, together, are the ridge -- and it is the *sd* the message is
  # about, named as the model names it rather than as the raw vector does.
  expect_identical(result$partial, "diff_eta1")
  expect_identical(result$partners, "drift_eta1")
  expect_length(result$structural, 0L)
  expect_true(all(c("popsd_diff_eta1", "rawcor_diff_eta1__drift_eta1") %in%
    result$ridge))
  # The advice, which is the whole point: the covariance is estimable, so
  # fixing one of the set throws it away.
  # `strwrap` folds the advice, so the comparison is on whitespace-normalised
  # text rather than on where the wrap happened to fall.
  printed <- gsub("[[:space:]]+", " ",
    paste(utils::capture.output(print(result)), collapse = " "))
  expect_match(printed, "only the covariances they generate are")
  expect_match(printed, "intoverpop='laplace'", fixed = TRUE)
  expect_false(grepl("not estimable from this data", printed))
  # Aggregating by name put the correlation here under text saying a fit may
  # still settle. It will not.
  expect_length(result$sometimes, 0L)
})

test_that("a coordinate in the flat subspace gets no reported spread, and is said to", {
  skip_without_julia()
  # Fit free, and deliberately so: the information matrix is the whole input to
  # the covariance a fit would report, so the fault is reachable in seconds
  # through the route `ctIdentify()` already uses rather than by spending a fit
  # to reproduce it. Same model as the partial-identification test above.
  spec <- suppressMessages(.ctFitJuliaBackend(.identify_data(),
    .identify_varying_model(), fit = FALSE, priors = FALSE,
    intoverpop = "augmented", cores = 1, verbose = 0))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  info <- ctsem:::.ctIdentifyInformation(spec, numeric(npar))
  parnames <- ctsem:::.ctBackendRawParameterNames(list(model_spec = spec), npar)

  # Exactly what a fit does with it: project the flat directions out, invert,
  # report the square roots of the diagonal.
  cov <- suppressWarnings(suppressMessages(
    ctsem:::ctOptimCovFromHessian(-info$information, warn = FALSE)))
  se <- sqrt(diag(cov))
  check <- ctsem:::.ctBackendIntervalCheck(-info$information, se, parnames)

  # The population sd of the diffusion effect is the coordinate the augmented
  # filter cannot see -- its carrier state enters only the predicted covariance,
  # so the Kalman update never moves it.
  expect_true("popsd_diff_eta1" %in% check$unidentified)
  expect_gte(check$nunidentified, 1L)
  mass <- check$table$nullmass[match(parnames, check$table$param)]
  names(mass) <- parnames
  # It is the whole of the flat direction here, and the identified coordinates
  # carry none of it -- not merely little, which is what makes the threshold a
  # threshold rather than a tuning.
  expect_equal(unname(mass["popsd_diff_eta1"]), 1, tolerance = 1e-8)
  expect_lt(max(mass[setdiff(parnames, "popsd_diff_eta1")]), 1e-10)

  # And what the user would otherwise read: a standard error of zero, which
  # prints as a point estimate with a zero-width interval rather than as a
  # parameter the data says nothing about.
  expect_equal(unname(se[match("popsd_diff_eta1", parnames)]), 0,
    tolerance = 1e-12)
  # The width check cannot reach it: a coordinate with no curvature of its own
  # has no ratio to be large.
  expect_false("popsd_diff_eta1" %in% check$parameters)

  # Said, rather than only computed -- in the brief register the warning uses,
  # since R truncates a warning at 1000 bytes.
  expect_warning(ctsem:::.ctBackendIdentifyWarn(NULL, data.frame(), check),
    "No width at all for")
})

test_that("a genuinely redundant parameter still gets the fix-or-remove advice", {
  skip_without_julia()
  # Not a ridge: LAMBDA free against the latent scale leaves nothing about the
  # pair determined, and the message must not soften.
  result <- suppressMessages(ctIdentify(.identify_data(),
    .identify_model("lambda"), cores = 1))
  expect_gte(result$nweak, 1L)
  expect_length(result$partial, 0L)
  expect_true(all(c("lambda", "diff_eta1") %in% result$structural))
  printed <- gsub("[[:space:]]+", " ",
    paste(utils::capture.output(print(result)), collapse = " "))
  expect_match(printed, "not estimable from this data as the model stands")
  expect_false(grepl("only the covariances they generate", printed))
})

test_that("the flat subspace is aggregated by subspace, not by name", {
  skip_without_julia()
  result <- suppressMessages(ctIdentify(.identify_data(),
    .identify_varying_model(), cores = 1))
  # The loading rotates along the ridge, so the names implicated at one point
  # are not the names implicated at the next. Intersecting them emptied
  # `$parameters` and moved everything into `$sometimes`; the union of what
  # lies in the flat subspace is the invariant statement.
  perpoint <- lapply(result$starts, function(start)
    sort(unique(unlist(lapply(start$directions, `[[`, "parameters")))))
  expect_true(result$rotating)
  expect_false(all(vapply(perpoint,
    function(names) setequal(names, result$parameters), logical(1))))
  expect_true(all(unlist(perpoint) %in% result$parameters))
})

test_that("the same advice is given after a fit as before one", {
  skip_without_julia()
  # `.ctBackendIdentifyWarn()` is the post-fit site and takes no model, so the
  # classification travels on the identifiability object. One helper writes the
  # text at both sites; this checks the post-fit one uses it -- in its `brief`
  # register, because a warning is truncated at 1000 bytes and the paragraph
  # `print.ctIdentify()` prints ran past that (the test above asserts the long
  # form on the print path).
  partial <- list(nweak = 1L, parameters = "popsd_diff_eta1", negative = 0L,
    directions = list(list(parameters = "popsd_diff_eta1",
      partial = list(parameters = "diff_eta1", partners = "drift_eta1"))))
  expect_warning(ctsem:::.ctBackendIdentifyWarn(partial, data.frame()),
    "only the covariance is")
  complete <- list(nweak = 1L, parameters = "lambda", negative = 0L,
    directions = list(list(parameters = c("lambda", "diff_eta1"))))
  expect_warning(ctsem:::.ctBackendIdentifyWarn(complete, data.frame()),
    "Not estimable as the model stands")
})

test_that("the post-fit warning fits inside R's warning length", {
  # R truncates a warning at getOption("warning.length"), 1000 bytes by
  # default, and it truncates the *end* -- which is where the advice and the
  # object to look at were. Both cases are measured here rather than trusted,
  # since the text is assembled from several pieces and any of them can grow.
  intervals <- list(nunidentified = 2L, unidentified = c("popsd_diff_eta1",
    "rawcor_diff_eta1__drift_eta1"))
  partial <- list(nweak = 1L, parameters = "popsd_diff_eta1", negative = 0L,
    directions = list(list(parameters = "popsd_diff_eta1",
      partial = list(parameters = "diff_eta1", partners = "drift_eta1"))))
  text <- tryCatch(
    ctsem:::.ctBackendIdentifyWarn(partial, data.frame(), intervals),
    warning = function(w) conditionMessage(w))
  expect_lt(nchar(text, type = "bytes"), 1000L)
  # The tail is the part that was being lost, so it is the part asserted.
  expect_match(text, "fit\\$identifiability")
  expect_match(text, "summary\\(\\) reports")

  complete <- list(nweak = 1L, parameters = "lambda", negative = 0L,
    directions = list(list(parameters = c("lambda", "diff_eta1"))))
  text <- tryCatch(
    ctsem:::.ctBackendIdentifyWarn(complete, data.frame(), intervals),
    warning = function(w) conditionMessage(w))
  expect_lt(nchar(text, type = "bytes"), 1000L)
  expect_match(text, "fit\\$identifiability")
})

test_that("intoverpop='laplace' works from ctIdentify", {
  skip_without_julia()
  # It is the route the partial-identification message recommends, and it used
  # to stop with "no parameters are marked indvarying" on a model that plainly
  # has them: `.ctJuliaLaplaceSpec()` reads the varying set from the prepared
  # data or from `modelmats`, and a model straight from `ctModel()` has
  # neither.
  result <- suppressMessages(ctIdentify(.identify_data(),
    .identify_varying_model(), cores = 1, intoverpop = "laplace"))
  expect_s3_class(result, "ctIdentify")
  expect_true(all(c("popsd_drift_eta1", "popsd_diff_eta1",
    "rawcor_diff_eta1__drift_eta1") %in% result$parnames))
  # And it does what the message says: the sd it could not separate under
  # 'augmented' is informed here.
  expect_length(result$partial, 0L)
  expect_false("popsd_diff_eta1" %in% result$parameters)
})

test_that("the laplace specification agrees with the one ctFit builds", {
  skip_without_julia()
  # The fallback added for the unprepared model must find the same random
  # effects `.ctPrepareData()` does, or `ctIdentify(intoverpop='laplace')`
  # would answer a question about a different model than the fit does.
  data <- .identify_data()
  model <- .identify_varying_model()
  bare <- suppressMessages(.ctFitJuliaBackend(data, model, fit = FALSE,
    intoverpop = "laplace", cores = 1, verbose = 0))
  prepared <- suppressMessages(ctFit(data, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE, cores = 1, verbose = 0))
  expect_equal(bare$laplace$nrandom, prepared$laplace$nrandom)
  expect_equal(sort(bare$laplace$param), sort(prepared$laplace$param))
  expect_equal(bare$laplace$sd_scale[order(bare$laplace$param)],
    prepared$laplace$sd_scale[order(prepared$laplace$param)])
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
