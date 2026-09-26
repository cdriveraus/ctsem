# intoverpop = 'auto': which route an optimised model with random effects takes.
#
# The augmented route -- each random effect a static latent state, integrated
# by the filter -- is exact and the cheapest there is for an effect that shifts
# a mean affinely with Gaussian indicators, and measurably the wrong estimator
# elsewhere; `.ctIntOverPopAuto()` records the measurements. So with
# backend = 'julia' and optimize = TRUE, 'auto' takes laplace wherever one of
# the triggers below holds, and says so in one line.
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

test_that("a linear-Gaussian random intercept stays on the augmented route", {
  skip_on_cran()
  for (m in list(
    .auto_model("mm", LAMBDA = matrix(1), MANIFESTMEANS = matrix("mm")),
    .auto_model("cint", LAMBDA = matrix(1), CINT = matrix("cint")),
    .auto_model("t0m", LAMBDA = matrix(1), T0MEANS = matrix("t0m")),
    # An affine transform of the user's own is still affine.
    .auto_model("cint", LAMBDA = matrix(1), CINT = matrix("cint|param*2+1|TRUE")))) {
    r <- .auto_resolve(m)
    expect_identical(r$route, "augmented")
    expect_match(r$reason, "affinely")
    expect_null(r$prepared$laplace)
    expect_false(grepl("intoverpop='auto'", r$said, fixed = TRUE))
  }
})

test_that("a random effect in DRIFT, DIFFUSION, MANIFESTVAR or LAMBDA takes laplace", {
  skip_on_cran()
  cases <- list(
    DRIFT = .auto_model("dr", LAMBDA = matrix(1), DRIFT = matrix("dr")),
    DIFFUSION = .auto_model("df", LAMBDA = matrix(1), DIFFUSION = matrix("df")),
    MANIFESTVAR = .auto_model("mv", LAMBDA = matrix(1), MANIFESTVAR = matrix("mv")),
    LAMBDA = .auto_model("lam", LAMBDA = matrix("lam")))
  for (cell in names(cases)) {
    r <- .auto_resolve(cases[[cell]])
    expect_identical(r$route, "laplace", info = cell)
    expect_match(r$reason, paste0("individual variation on ", cell), info = cell)
    expect_false(is.null(r$prepared$laplace), info = cell)
    # One line, naming the route and the reason.
    expect_match(r$said, paste0("intoverpop='auto' chose 'laplace': individual ",
      "variation on ", cell), fixed = TRUE, info = cell)
  }
})

test_that("a mean cell with a nonlinear transform takes laplace", {
  skip_on_cran()
  r <- .auto_resolve(.auto_model("cint", LAMBDA = matrix(1),
    CINT = matrix("cint|log1p_exp(param)|TRUE")))
  expect_identical(r$route, "laplace")
  expect_match(r$reason, "the transform of cint in CINT is not affine", fixed = TRUE)
  # The transform codes and text forms, read before they are turned into codes.
  affine <- ctsem:::.ctTransformAffine
  expect_true(affine(0)); expect_true(affine("0"))
  for (code in 1:5) expect_false(affine(code), info = code)
  expect_true(affine("param * 10")); expect_true(affine("param"))
  expect_true(affine(NA))
  expect_false(affine("exp(param)")); expect_false(affine("log1p_exp(param * -2) * -2"))
  expect_false(affine("1/(1+exp(-param))")); expect_false(affine("param^3"))
  # Text that cannot be evaluated is not called affine.
  expect_false(affine("notafunction(param)"))
})

test_that("a varying parameter inside another cell's expression takes laplace", {
  skip_on_cran()
  # Through a state-dependent expression in a mean cell.
  m <- .auto_model("b", LAMBDA = matrix(1), PARS = matrix("b"),
    CINT = matrix("b * eta1"))
  r <- .auto_resolve(m)
  expect_identical(r$route, "laplace")
  expect_match(r$reason, "enters a state-dependent expression in CINT", fixed = TRUE)
  # Through DRIFT, from a cell that is otherwise a linear intercept. A name in
  # an expression has to be declared in PARS, which makes the PARS cell and the
  # CINT cell one parameter.
  m <- .auto_model("cint", LAMBDA = matrix(1), PARS = matrix("cint"),
    CINT = matrix("cint"), DRIFT = matrix("-exp(cint)"))
  r <- .auto_resolve(m)
  expect_identical(r$route, "laplace")
  expect_match(r$reason, "cint reaches DRIFT", fixed = TRUE)
})

test_that("non-Gaussian indicators with any random effect take laplace", {
  skip_on_cran()
  m <- .auto_model("mm", LAMBDA = matrix(1), MANIFESTMEANS = matrix("mm"),
    manifesttype = 1L)
  r <- .auto_resolve(m, data = .auto_data(binary = TRUE))
  expect_identical(r$route, "laplace")
  expect_match(r$reason, "non-Gaussian indicators")
})

test_that("the routes 'auto' does not change", {
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
  # a flag there alone does not send the model to laplace.
  t0 <- .auto_model(character(0), LAMBDA = matrix(1))
  t0$pars$indvarying <- t0$pars$matrix %in% "T0VAR" & is.na(t0$pars$value) &
    t0$pars$row == t0$pars$col
  expect_identical(.auto_resolve(t0)$route, "augmented")
})

test_that("an outer level takes laplace without announcing it", {
  skip_on_cran()
  m <- .auto_model("mm", LAMBDA = matrix(1), MANIFESTMEANS = matrix("mm"))
  m$groupIDnames <- "study"
  m$pars$indvarying_study <- m$pars$param %in% "mm"
  auto <- ctsem:::.ctIntOverPopAuto(m, backend = "julia")
  expect_identical(auto$route, "laplace")
  expect_false(auto$announce)
  # Whatever the backend: stan then refuses it by name.
  expect_identical(ctsem:::.ctIntOverPopAuto(m, backend = "stan")$route, "laplace")
})

test_that("the variance-cell warning still fires for an explicit augmented route", {
  skip_on_cran()
  m <- .auto_model("df", LAMBDA = matrix(1), DIFFUSION = matrix("df"))
  expect_warning(suppressMessages(ctFit(.auto_data(), m, backend = "julia",
    intoverpop = "augmented", fit = FALSE)), "only partially identified")
})
