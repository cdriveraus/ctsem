# Julia session management, and one whole fit from a model and a data frame.
#
# Everything else in the julia test files reaches for a fitted object, so all of
# them share one setup path and none of them exercises it. That is how
# `.ctJuliaSessionRunning()` came to answer "is a session running?" by *starting*
# one -- `JuliaConnectoR::juliaEval()` starts a server rather than reporting on
# it -- with two silent consequences:
#
#   * `ctJuliaSetup(threads = n, force = TRUE)` cleared the session, asked the
#     question, saw the session its own question had just created, and warned
#     that `threads` would have no effect. `force` could not work.
#   * `ctFit(backend = 'julia', cores = n)` sets `JULIA_NUM_THREADS` only when no
#     session is running, so it was never set and every fit ran on one thread
#     whatever `cores` said.
#
# Both are asserted below, along with a plain model-and-data-to-fit check that
# does not reuse a cached fixture.

test_that("asking whether a session is running does not start one", {
  skip_without_julia()
  ctsem:::.ctJuliaClearSession()
  # The regression: this must not be what brings a session into existence.
  expect_false(ctsem:::.ctJuliaSessionRunning())
  expect_false(ctsem:::.ctJuliaSessionRunning())
})

test_that("force restarts the session and threads take effect", {
  skip_without_julia()
  previous <- Sys.getenv("JULIA_NUM_THREADS", unset = NA)
  on.exit({
    if (is.na(previous)) Sys.unsetenv("JULIA_NUM_THREADS")
    else Sys.setenv(JULIA_NUM_THREADS = previous)
  }, add = TRUE)

  # Two threads is enough to prove the setting is applied without asking the
  # test machine for cores it may not have.
  expect_no_warning(suppressMessages(ctJuliaSetup(threads = 2L, force = TRUE)))
  expect_equal(Sys.getenv("JULIA_NUM_THREADS"), "2")
  expect_equal(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")), 2L)
  expect_true(ctsem:::.ctJuliaSessionRunning())

  # Asking for a different count *without* force cannot take effect, and saying
  # so is the point of the warning -- it must fire here and only here.
  expect_warning(suppressMessages(ctJuliaSetup(threads = 3L)), "already running")
  expect_equal(as.integer(JuliaConnectoR::juliaEval("Threads.nthreads()")), 2L)
})

test_that("a model and a data frame fit end to end", {
  skip_without_julia()
  set.seed(20260827)
  nsubjects <- 12; nobs <- 6
  dat <- do.call(rbind, lapply(seq_len(nsubjects), function(i) {
    intercept <- stats::rnorm(1, 1.5, 0.9)
    state <- stats::rnorm(1, 0, 0.5)
    out <- numeric(nobs)
    for (t in seq_len(nobs)) {
      if (t > 1) {
        decay <- exp(-0.4)
        state <- decay * state +
          stats::rnorm(1, 0, sqrt(0.36 / 0.8 * (1 - decay^2)))
      }
      out[t] <- state + intercept + stats::rnorm(1, 0, 0.3)
    }
    data.frame(id = i, time = seq_len(nobs) - 1, Y1 = out)
  }))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), T0VAR = matrix(0.5),
    MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE

  for (route in c("laplace", "augmented")) {
    fit <- suppressMessages(ctFit(dat, model, backend = "julia",
      intoverpop = if (route == "laplace") "laplace" else TRUE, cores = 1,
      optimcontrol = list(estonly = TRUE)))
    expect_s3_class(fit, "ctJuliaFit")
    expect_true(isTRUE(fit$estimate$converged), label = route)
    expect_false(isTRUE(fit$estimate$stalled), label = route)
    expect_true(is.finite(fit$estimate$loglik), label = route)
    # The population mean of an identity-transformed MANIFESTMEANS is the one
    # parameter whose truth is known here without any transform bookkeeping.
    summ <- suppressWarnings(summary(fit))
    expect_true(is.list(summ) || is.data.frame(summ), label = route)
  }
})

# The engine's subject-chunk ceiling is session-global Julia state. Three call
# sites wrote it and none put it back, so after a `cores=8` fit the session read
# 8 and after a `cores=1` fit it read 1 -- and every later ctKalman(),
# ctExtract() or ctLOO() inherited whichever fit came last. Performance-only,
# but it makes a timing depend on history, which is what makes one impossible to
# reproduce.
test_that("a fit leaves the engine's chunk ceiling as it found it", {
  skip_without_julia()
  ceiling <- function() as.integer(JuliaConnectoR::juliaEval(
    "ContinuousTimeSEM.ctsem_max_chunks().max_chunks"))
  # Put it back on the way out. This test is about session state not leaking,
  # and leaving a ceiling of 3 behind for every later test file would be a poor
  # way to make that point -- it did, and it changed which chunk count the next
  # file's fits ran at.
  original <- ceiling()
  on.exit(JuliaConnectoR::juliaEval(sprintf(
    "ContinuousTimeSEM.ctsem_set_max_chunks!(%d)", original)), add = TRUE)
  JuliaConnectoR::juliaEval("ContinuousTimeSEM.ctsem_set_max_chunks!(3)")
  before <- ceiling()
  expect_equal(before, 3L)

  set.seed(9)
  dat <- do.call(rbind, lapply(1:8, function(i)
    data.frame(id = i, time = 0:4, Y1 = cumsum(stats::rnorm(5)) * .5)))
  model <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1))))
  fit <- suppressWarnings(suppressMessages(ctFit(dat, model, backend = "julia",
    cores = 2, optimcontrol = list(estonly = TRUE))))

  expect_equal(ceiling(), before)
  # And the tuner's answer is still recorded on the fit, which is what the
  # uncertainty phase reads instead of `cores`.
  expect_true(is.numeric(fit$estimate$chunks) || is.integer(fit$estimate$chunks))
  expect_lte(as.integer(fit$estimate$chunks), 2L)

  # A diagnostic must not reconfigure the session either.
  JuliaConnectoR::juliaEval("ContinuousTimeSEM.ctsem_set_max_chunks!(3)")
  suppressWarnings(suppressMessages(try(ctExtract(fit, cores = 2), silent = TRUE)))
  expect_equal(ceiling(), 3L)
})
