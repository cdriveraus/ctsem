# The optimisation trace and the live callback.
#
# Two consumers of the same information with different costs: the trace records
# every iteration in Julia and crosses the bridge once, the callback crosses it
# while the fit runs and so is rate-limited. The tests that matter are that the
# trace is always there, that the callback is not tied to `verbose`, and that a
# callback which fails cannot take the fit with it.

.trace_data <- function(nsubjects = 60, nobs = 8) {
  set.seed(5)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.4), DIFFUSION = matrix(0.6), MANIFESTVAR = matrix(0.3),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix(0), Tpoints = nobs))
  ctGenerate(gen, n.subjects = nsubjects, Tpoints = nobs, backend = "r")
}

.trace_model <- function() {
  suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(0)))
}

.trace_fit <- function(...) {
  suppressWarnings(suppressMessages(ctFit(.trace_data(), .trace_model(),
    backend = "julia", cores = 2, verbose = 0,
    optimcontrol = c(list(estonly = TRUE), list(...)))))
}

test_that("the trace is recorded even with reporting off", {
  skip_on_cran()
  skip_without_julia()
  fit <- .trace_fit()
  trace <- fit$trace
  expect_s3_class(trace, "data.frame")
  expect_gt(nrow(trace), 1L)
  expect_named(trace, c("iteration", "objective", "gradient_norm"))
  # `verbose = 0` above: the fit whose trace turns out to be worth reading is
  # the one nobody thought to turn reporting on for.
  expect_true(all(diff(trace$iteration) > 0))
  # Ascent, and ending where the fit says it ended.
  expect_true(all(diff(trace$objective) >= -1e-8))
  expect_equal(trace$objective[nrow(trace)], fit$estimate$logposterior,
    tolerance = 1e-6)
  expect_lt(trace$gradient_norm[nrow(trace)], trace$gradient_norm[1])
})

test_that("the callback fires while the fit runs, rate limited, without verbose", {
  skip_on_cran()
  skip_without_julia()
  seen <- new.env()
  seen$iterations <- integer()
  fit <- .trace_fit(callback = function(iteration, total, objective, gradnorm) {
    seen$iterations <- c(seen$iterations, as.integer(iteration))
    NULL
  })
  # Not tied to `verbose`: a front end that draws rather than prints passes a
  # callback and leaves verbose at zero, which reported nothing at all when the
  # callback shared the printed line's enable flag.
  expect_gt(length(seen$iterations), 0L)
  # Rate limited on time, so fewer calls than iterations on any fit that is not
  # pathologically slow.
  expect_lte(length(seen$iterations), nrow(fit$trace))
  # The last iteration is always reported, forced past the cadence, because a
  # fit finishing inside one interval would otherwise report nothing.
  expect_equal(max(seen$iterations), fit$estimate$iterations)
})

test_that("a callback that errors is disabled and the fit survives", {
  skip_on_cran()
  skip_without_julia()
  # An error thrown out of an R callback never reaches the engine: it aborts
  # before a reply is sent and leaves the JuliaConnectoR bridge desynchronised,
  # so the fit is lost to a fault in a reporting function. It has to be caught
  # on the R side.
  reference <- .trace_fit()
  warned <- NULL
  fit <- withCallingHandlers(
    suppressMessages(ctFit(.trace_data(), .trace_model(), backend = "julia",
      cores = 2, verbose = 0, optimcontrol = list(estonly = TRUE,
        callback = function(...) stop("deliberate")))),
    warning = function(w) {
      if (grepl("progress callback failed", conditionMessage(w))) {
        warned <<- conditionMessage(w)
      }
      invokeRestart("muffleWarning")
    })
  expect_true(isTRUE(fit$estimate$converged))
  expect_equal(fit$estimate$loglik, reference$estimate$loglik, tolerance = 1e-6)
  expect_false(is.null(warned))
  expect_match(warned, "deliberate")
})

test_that("a callback does not change the answer", {
  skip_on_cran()
  skip_without_julia()
  plain <- .trace_fit()
  watched <- .trace_fit(callback = function(...) NULL)
  expect_equal(watched$estimate$raw, plain$estimate$raw, tolerance = 1e-6)
})

test_that("the laplace route traces the inner solve as well", {
  skip_on_cran()
  skip_without_julia()
  model <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix("mm")))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mm"] <- TRUE
  fit <- suppressWarnings(suppressMessages(ctFit(.trace_data(), model,
    backend = "julia", cores = 2, intoverpop = "laplace", verbose = 0,
    optimcontrol = list(estonly = TRUE))))
  # A Laplace fit that stalls usually stalls in the inner solve, which the
  # outer objective alone does not show.
  expect_true("inner_converged" %in% names(fit$trace))
  expect_gt(nrow(fit$trace), 1L)
})

test_that("ctTracePlot draws, and refuses a fit with no trace", {
  skip_on_cran()
  skip_without_julia()
  fit <- .trace_fit()
  file <- file.path(tempdir(), "ctsem-trace-test.png")
  on.exit(unlink(file), add = TRUE)
  grDevices::png(file, width = 600, height = 500)
  result <- ctTracePlot(fit)
  grDevices::dev.off()
  expect_true(file.exists(file))
  expect_gt(file.info(file)$size, 0)
  expect_s3_class(result, "data.frame")

  bare <- fit
  bare$trace <- NULL
  expect_error(ctTracePlot(bare), "no optimisation trace")
})

test_that("a callback that is not a function is refused before fitting", {
  skip_on_cran()
  skip_without_julia()
  expect_error(
    suppressMessages(ctFit(.trace_data(), .trace_model(), backend = "julia",
      cores = 1, optimcontrol = list(estonly = TRUE, callback = "nope"))),
    "must be a function")
})
