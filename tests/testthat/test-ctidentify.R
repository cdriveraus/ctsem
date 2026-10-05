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
  ctGenerate(gen, n = nsubjects, Tpoints = nobs, backend = "r")
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
  # which cell of which matrix to change. T0VAR is in the trade too, since the
  # latent's scale sets the initial variance as well as the diffusion, though
  # it carries only 0.002 to 0.065 of the flat direction across the three
  # evaluation points: a norm of 0.25 named it at one point, by 0.004.
  expect_true(all(c("lambda", "diff_eta1", "T0var_eta1") %in% result$parameters))
  # Named at every point, so the direction is not reported as turning between
  # them: it moves, but it keeps the same three coordinates.
  expect_false(result$rotating)
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

# A fit's report names what the likelihood measured flat, not only what the
# curvature has decayed far enough to call flat at the point the optimiser
# stopped. The stan-julia parity fixture's ridge of ten correlations sat at
# 1.7e-08 of the sharpest curvature at its earliest stopping point, above the
# eigenvalue rule's 1e-8, and was named at later points of the same ridge
# 4.4e-04 nats higher; its tenth correlation's loading drifted from 0.19 to
# 0.28 across them, either side of the absolute 0.25 that used to decide it.
test_that("the report names what the likelihood measured flat, at any curvature", {
  parnames <- c("a", "b", "c")
  flatline <- function(p) -0.5 * p[1]^2 - 0.25 * p[3]^2
  info <- diag(c(1, 3e-8, 0.5))
  screen <- ctsem:::.ctOptimFlatDirectionScreen(info, flatline, numeric(3))
  expect_equal(sum(screen$flat), 1L)

  # The curvature alone names nothing: 3e-8 is above its 1e-8.
  plain <- ctsem:::.ctBackendIdentifiability(-info, parnames)
  expect_equal(plain$nweak, 0L)
  # The likelihood's verdict does, raw or as the uncertainty stage stores it.
  measured <- ctsem:::.ctBackendIdentifiability(-info, parnames,
    screen = screen)
  expect_equal(measured$nweak, 1L)
  expect_identical(measured$parameters, "b")
  expect_identical(measured$directions[[1]]$evidence, "likelihood")
  expect_equal(measured$directions[[1]]$change, 0)
  expect_equal(measured$screen$candidates, 1L)
  stored <- list(n = 1L, bar = screen$bar, lengths = screen$lengths,
    rtol = screen$rtol, candidates = length(screen$candidates),
    evaluations = screen$evaluations, change = screen$change[screen$flat],
    eigenvalue = screen$eig$values[screen$flat] / max(screen$eig$values),
    vectors = screen$eig$vectors[, screen$flat, drop = FALSE])
  fromfit <- ctsem:::.ctBackendIdentifiability(-info, parnames, screen = stored)
  expect_identical(fromfit$parameters, "b")
  expect_equal(fromfit$directions[[1]]$relative, 3e-8)
  expect_equal(fromfit$screen$candidates, 1L)

  # Below the eigenvalue rule and not confirmed -- the walk rose past the bar --
  # it is still named, because a rise proves nothing; it says so.
  steep <- function(p) -0.5 * sum(p^2)
  curved <- diag(c(1, 1e-12, 0.5))
  unconfirmed <- ctsem:::.ctOptimFlatDirectionScreen(curved, steep, numeric(3))
  expect_false(any(unconfirmed$flat))
  kept <- ctsem:::.ctBackendIdentifiability(-curved, parnames,
    screen = unconfirmed)
  expect_identical(kept$parameters, "b")
  expect_identical(kept$directions[[1]]$evidence, "curvature")
  expect_true(is.na(kept$directions[[1]]$change))

  # Saturation is the reason when there is one, not the condition for naming.
  saturatedfit <- list(optim = list(saturated_parameters = "b"))
  reasoned <- ctsem:::.ctBackendIdentifiability(-info, parnames,
    fit = saturatedfit, screen = screen)
  expect_identical(reasoned$directions[[1]]$saturated, "b")
  expect_identical(measured$directions[[1]]$saturated, character())
})

test_that("a direction names every coordinate carrying a share of it", {
  # Ten coordinates along one flat direction, the smallest two at half the
  # largest. An absolute loading of 0.25 dropped them; a share of the direction
  # at the rounding floor does not, however many coordinates it is spread over.
  v <- c(0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2)
  v <- v / sqrt(sum(v^2))
  expect_lt(min(abs(v)), 0.25)
  information <- diag(10) - tcrossprod(v)
  report <- ctsem:::.ctBackendIdentifiability(-information, paste0("p", 1:10))
  expect_equal(report$nweak, 1L)
  expect_setequal(report$parameters, paste0("p", 1:10))
})

test_that("a fit names what lies in the flat subspace, in any basis, as its intervals do", {
  # The parity fixture's flat direction turns as the optimiser walks its ridge,
  # and the weaker of the ten correlations on it carry 0.005 to 0.07 of it
  # depending on where the fit stopped (test-stan-julia-parity.R). Named by a
  # third of the largest loading they came and went -- nine names at one
  # stopping point, seven at another; named by their share against
  # `.ctNullMassBar()`, which sits at rounding, they do not.
  parnames <- paste0("p", 1:6)
  v <- c(0.8, 0.6, 0.1, 1e-5, 0, 0)
  v <- v / sqrt(sum(v^2))
  information <- diag(6) - tcrossprod(v)
  report <- ctsem:::.ctBackendIdentifiability(-information, parnames)
  expect_equal(report$nweak, 1L)
  # p3 at an eighth of the largest loading is on the ridge and named, largest
  # share first; a leak of 1e-5 is rounding's size and is not.
  expect_identical(report$parameters, c("p1", "p2", "p3"))
  expect_identical(report$directions[[1]]$parameters, report$parameters)
  # The same coordinates as the intervals with no width. This direction's
  # curvature is zero, so the covariance drops it and the interval check
  # measures the same subspace against the same bar.
  check <- ctsem:::.ctBackendIntervalCheck(-information, rep(1, 6), parnames)
  expect_setequal(check$unidentified, report$parameters)

  # Two flat directions: any rotation of the pair is as good a basis as the
  # one a decomposition returns, and a coordinate's share of the subspace is
  # the same from all of them.
  w <- c(0, 0.1, 0.2, 0.9, 0, 0.3)
  w <- w - sum(w * v) * v
  w <- w / sqrt(sum(w^2))
  share <- ctsem:::.ctIdentifySubspaceShare(list(list(vector = v),
    list(vector = w)), 6L)
  expect_equal(sum(share), 2)
  turn <- 0.6
  rotated <- list(list(vector = cos(turn) * v + sin(turn) * w),
    list(vector = -sin(turn) * v + cos(turn) * w))
  expect_equal(ctsem:::.ctIdentifySubspaceShare(rotated, 6L), share)
  both <- diag(6) - tcrossprod(v) - tcrossprod(w)
  pair <- ctsem:::.ctBackendIdentifiability(-both, parnames)
  expect_equal(pair$nweak, 2L)
  expect_setequal(pair$parameters, parnames[share >= ctsem:::.ctNullMassBar()])

  # The same direction measured twice -- a screen's vector from a slightly
  # different matrix beside the eigenvector it did not quite match -- names
  # nothing that neither carries a share of. Orthogonalising the pair would
  # have made their small difference, spread over p5 and p6, a whole direction.
  again <- v + c(0, 0, 0, 0, 0.02, 0.02)
  remeasured <- ctsem:::.ctIdentifySubspaceShare(list(list(vector = v),
    list(vector = again)), 6L)
  expect_true(all(remeasured[5:6] < ctsem:::.ctNullMassBar()))
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

test_that("partial identification is a finding about the augmented route only", {
  skip_without_julia()
  # The mechanism is the augmented filter's -- a variance cell's carrier state
  # the update never moves. On the Laplace route the likelihood sees a level's
  # scales and correlations only through the covariance they build, so a flat
  # direction there that moves a variance is that variance undetermined. A
  # laplace fit whose flat sds had correlations near zero (AnomAuth S1) was
  # classified partial all the same, since a cross-covariance through a zero
  # correlation does not move with the scale, and was told it was under
  # intoverpop='augmented' and that 'laplace' would identify it. So the
  # classification has blocks on the augmented route and none on Laplace.
  data <- .identify_data()
  model <- .identify_varying_model()
  augmented <- suppressMessages(ctFit(data, model, backend = "julia",
    intoverpop = "augmented", fit = FALSE, cores = 1, verbose = 0))
  laplace <- suppressMessages(ctFit(data, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE, cores = 1, verbose = 0))
  blocks <- ctsem:::.ctIdentifyBlocks(augmented)
  expect_length(blocks, 1L)
  expect_identical(blocks[[1L]]$route, "augmented")
  expect_length(ctsem:::.ctIdentifyBlocks(laplace), 0L)
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
