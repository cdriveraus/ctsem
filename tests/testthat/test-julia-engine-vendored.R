# The Julia engine is part of ctsem, so `ctJuliaSetup()` needs no network access
# and no repository credentials.
#
# This exists because the original arrangement -- `Pkg.add` from a URL pinned in
# a lock file -- pointed at a private GitLab repository, which meant
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

test_that("the engine version is a hash of the engine source", {
  # Identity by content, not by a maintained version string. The cached project
  # directory is keyed on this and is only populated when empty, so an
  # identifier that had to be updated by hand would leave anyone who had already
  # run the backend on a stale engine after an edit -- silently. Editing the
  # engine has to change this, which is what the second half asserts.
  version <- ctsem:::.ctJuliaEngineVersion()
  expect_match(version, "^[0-9a-f]{12}$")
  expect_true(grepl(version, ctsem:::.ctJuliaEnvDir(), fixed = TRUE))

  engine <- system.file("julia", "ContinuousTimeSEM", package = "ctsem")
  altered <- file.path(tempdir(), "ContinuousTimeSEM-altered")
  unlink(altered, recursive = TRUE)
  dir.create(altered, recursive = TRUE)
  file.copy(list.files(engine, full.names = TRUE), altered, recursive = TRUE)
  expect_identical(ctsem:::.ctJuliaEngineVersion(altered), version)
  cat("
# a change", file = file.path(altered, "src", "r_interface.jl"), append = TRUE)
  expect_false(identical(ctsem:::.ctJuliaEngineVersion(altered), version))
  unlink(altered, recursive = TRUE)
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
  expect_match(status$engine, "^[0-9a-f]{12}$")

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
