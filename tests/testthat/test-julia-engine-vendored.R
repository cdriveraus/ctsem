# The Julia engine ships inside ctsem, so `ctJuliaSetup()` needs no network
# access and no repository credentials.
#
# This exists because the previous arrangement -- `Pkg.add` from a URL pinned in
# inst/julia/engine.json -- pointed at a private GitLab repository, which meant
# backend='julia' could not be installed by anyone outside that project. That
# was invisible from inside it, and no test or CI job would have caught it,
# because every check ran against a checkout that was already there. These
# assertions are cheap and need no Julia at all; the one that does is skipped
# unless Julia is present.

test_that("the Julia engine is vendored inside the installed package", {
  engine <- system.file("julia", "ContinuousTimeSEM", package = "ctsem")
  expect_true(nzchar(engine))
  expect_true(file.exists(file.path(engine, "Project.toml")))
  expect_true(file.exists(file.path(engine, "src", "ContinuousTimeSEM.jl")))
  expect_true(file.exists(file.path(engine, "src", "r_interface.jl")))
})

test_that("the engine lock records provenance rather than an install source", {
  lock <- ctsem:::.ctJuliaEngineLock()
  # A commit, not a branch name: a branch would move under a released ctsem.
  expect_match(lock$revision, "^[0-9a-f]{40}$")
  expect_true(nzchar(lock$branch))
  expect_true(nzchar(lock$url))
})

test_that("the vendored engine declares no heavyweight dependencies", {
  # DataFrames, Distributions, StaticArrays and Reexport were loaded for one
  # row-container use and three zero-use `@reexport`s, and cost 44 packages,
  # 144 MB and ~5 s of load time in every R session. If one comes back, it
  # should be a deliberate decision rather than a silent regression.
  project <- readLines(system.file("julia", "ContinuousTimeSEM", "Project.toml",
    package = "ctsem"), warn = FALSE)
  deps <- project[seq(which(project == "[deps]") + 1L, length(project))]
  deps <- deps[seq_len(which(!nzchar(deps))[1L] - 1L)]
  names <- trimws(sub("=.*", "", deps))
  expect_false(any(c("DataFrames", "Distributions", "StaticArrays", "Reexport") %in% names))
})

test_that("ctJuliaSetup works from the vendored copy, with no project argument", {
  skip_without_julia()
  skip_on_cran()

  status <- ctJuliaSetup()
  expect_true(status$available)
  expect_match(status$revision, "^[0-9a-f]{40}$")

  # A model with no TI predictors, so the optional table columns are omitted
  # rather than sent as empty vectors -- JuliaConnectoR hangs on those.
  model <- ctModel(type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1))
  data <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, .5, 1.5), 2),
    Y1 = c(0, .1, .2, .1, 0, -.1))

  julia_spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE,
    priors = FALSE))

  # A one-element parameter vector also exercises the scalar-marshalling branch
  # of .ctJuliaVector, which is the one length that must go across as a list.
  julia_value <- ctJuliaEvaluate(julia_spec, -.7, gradient = TRUE)
  expect_true(is.finite(as.numeric(julia_value$value)))
  expect_length(as.numeric(julia_value$gradient), 1L)

  # And the vendored engine computes the same thing Stan does, which is what
  # makes it a usable copy rather than merely an importable one.
  skip_if_not_installed("rstan")
  stan_spec <- suppressMessages(ctFit(data, model, backend = "stan", fit = FALSE,
    priors = FALSE))
  stan_value <- rstan::log_prob(
    ctsem:::stan_reinitsf(ctsem:::stanmodels$ctsm, stan_spec$standata),
    upars = -.7, adjust_transform = FALSE, gradient = TRUE)
  expect_equal(as.numeric(julia_value$value), as.numeric(stan_value), tolerance = 1e-8)
  expect_equal(as.numeric(julia_value$gradient),
    as.numeric(attributes(stan_value)$gradient), tolerance = 1e-8)
})
