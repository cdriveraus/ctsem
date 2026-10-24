# intoverpop = 'auto': which route an optimised model with random effects takes.
#
# On the julia backend, when optimising, 'auto' and TRUE take the Laplace
# route, which integrates an effect that shifts a mean in a linear Gaussian
# model inside each subject's filter, exactly, and every other effect by the
# Laplace approximation; `.ctIntOverPopAuto()` records why. Stan has one route,
# the state-explicit target composes with the augmented route only, sampling
# integrates nothing, and a named route is taken as named.
#
# Everything here is resolved with fit = FALSE, which is pure R on both
# backends: no Julia session and no Stan compile.

.auto_data <- function(binary = FALSE) {
  set.seed(5)
  do.call(rbind, lapply(1:6, function(i) {
    y <- if (binary) stats::rbinom(5, 1, 0.5) else stats::rnorm(5)
    data.frame(id = i, time = 0:4, Y1 = y)
  }))
}

.auto_model <- function(varying, ...) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", ...)))
  m$pars$indvarying <- m$pars$param %in% varying
  m
}

# The route and reason ctFit() resolved, and any message it printed.
.auto_resolve <- function(model, backend = "julia", data = .auto_data(), ...) {
  said <- character(0)
  prepared <- withCallingHandlers(suppressWarnings(ctFit(data, model,
    backend = backend, fit = FALSE, ...)),
    message = function(m) {
      said <<- c(said, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(route = prepared$args$resolved$intoverpop,
    reason = prepared$args$resolved$intoverpopreason,
    said = paste(said, collapse = ""), prepared = prepared)
}

test_that("every model with random effects takes laplace, and the filter takes a mean's", {
  skip_on_cran()
  # Effects that shift a mean, with the affine transform: the filter
  # integrates them, and nothing is announced.
  for (m in list(
    .auto_model("mm", LAMBDA = matrix(1), MANIFESTMEANS = matrix("mm")),
    .auto_model("cint", LAMBDA = matrix(1), CINT = matrix("cint")),
    .auto_model("t0m", LAMBDA = matrix(1), T0MEANS = matrix("t0m")))) {
    r <- .auto_resolve(m)
    expect_identical(r$route, "laplace")
    expect_match(r$reason, "the filter integrates")
    expect_false(is.null(r$prepared$laplace$hybrid))
    expect_false(grepl("intoverpop", r$said, fixed = TRUE))
  }
  # Everything else is a Laplace coordinate: an effect in DRIFT, DIFFUSION,
  # MANIFESTVAR or LAMBDA, a mean cell with a nonlinear transform, a varying
  # parameter inside another cell's expression (a state-dependent one, and
  # through DRIFT from a cell that is otherwise a linear intercept, which a
  # name in an expression declared in PARS makes one parameter), and any
  # effect with non-Gaussian indicators.
  cases <- list(
    DRIFT = .auto_model("dr", LAMBDA = matrix(1), DRIFT = matrix("dr")),
    DIFFUSION = .auto_model("df", LAMBDA = matrix(1), DIFFUSION = matrix("df")),
    MANIFESTVAR = .auto_model("mv", LAMBDA = matrix(1), MANIFESTVAR = matrix("mv")),
    LAMBDA = .auto_model("lam", LAMBDA = matrix("lam")),
    transform = .auto_model("cint", LAMBDA = matrix(1),
      CINT = matrix("cint|log1p_exp(param)|TRUE")),
    state = .auto_model("b", LAMBDA = matrix(1), PARS = matrix("b"),
      CINT = matrix("b * eta1")),
    reach = .auto_model("cint", LAMBDA = matrix(1), PARS = matrix("cint"),
      CINT = matrix("cint"), DRIFT = matrix("-exp(cint)")))
  for (cell in names(cases)) {
    r <- .auto_resolve(cases[[cell]])
    expect_identical(r$route, "laplace", info = cell)
    expect_false(is.null(r$prepared$laplace), info = cell)
    expect_null(r$prepared$laplace$hybrid, info = cell)
  }
  m <- .auto_model("mm", LAMBDA = matrix(1), MANIFESTMEANS = matrix("mm"),
    manifesttype = 1L)
  r <- .auto_resolve(m, data = .auto_data(binary = TRUE))
  expect_identical(r$route, "laplace")
  expect_null(r$prepared$laplace$hybrid)
})

test_that("named routes, stan, sampling and the state-explicit target", {
  skip_on_cran()
  drift <- .auto_model("dr", LAMBDA = matrix(1), DRIFT = matrix("dr"))
  # Named routes are taken as named.
  expect_identical(.auto_resolve(drift, intoverpop = "augmented")$route, "augmented")
  expect_true(is.na(.auto_resolve(drift, intoverpop = "augmented")$reason))
  # Stan has one route.
  r <- .auto_resolve(drift, backend = "stan")
  expect_identical(r$route, "augmented")
  expect_match(r$reason, "backend='stan'")
  # Sampling integrates nothing, and the state-explicit target composes with
  # the augmented route only.
  expect_identical(.auto_resolve(drift, optimize = FALSE)$route, "none")
  # Asked of the helper, because ctFit() refuses this model's free Gaussian
  # MANIFESTVAR with intoverstates = FALSE before it gets that far.
  expect_identical(ctsem:::.ctIntOverPopAuto(drift, backend = "julia",
    intoverstates = FALSE)$route, "augmented")
  # The automatic substep mesh does not change the route: the laplace route
  # chooses its mesh at the random-effect modes.
  expect_identical(.auto_resolve(drift,
    nlcontrol = list(nsubsteps = "auto"))$route, "laplace")
  # T0VAR cannot vary on any route -- the flag is cleared with a message -- so
  # a flag there alone does not send the model to laplace, which would refuse
  # it for having nothing to integrate.
  t0 <- .auto_model(character(0), LAMBDA = matrix(1))
  t0$pars$indvarying <- t0$pars$matrix %in% "T0VAR" & is.na(t0$pars$value) &
    t0$pars$row == t0$pars$col
  expect_identical(.auto_resolve(t0)$route, "augmented")
})

test_that("an outer level takes laplace", {
  skip_on_cran()
  m <- .auto_model("mm", LAMBDA = matrix(1), MANIFESTMEANS = matrix("mm"))
  # mm's indvarying is TRUE, which with a study level is both levels.
  m$groupIDnames <- "study"
  auto <- ctsem:::.ctIntOverPopAuto(m, backend = "julia")
  expect_identical(auto$route, "laplace")
  # Whatever the backend: stan then refuses it by name.
  expect_identical(ctsem:::.ctIntOverPopAuto(m, backend = "stan")$route, "laplace")
  # Sampling integrates nothing at any level: the effects of every level are
  # sampled with the rest.
  nested <- .auto_data()
  nested$study <- ifelse(nested$id <= 3, 1L, 2L)
  expect_identical(.auto_resolve(m, data = nested, optimize = FALSE)$route, "none")
})

test_that("the variance-cell warning still fires for an explicit augmented route", {
  skip_on_cran()
  m <- .auto_model("df", LAMBDA = matrix(1), DIFFUSION = matrix("df"))
  expect_warning(suppressMessages(ctFit(.auto_data(), m, backend = "julia",
    intoverpop = "augmented", fit = FALSE)), "only partially identified")
})
