test_that("Julia backend preparation is serializable and does not start Julia", {
  model <- suppressWarnings(ctModel(
    type = "ct",
    LAMBDA = diag(1),
    DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1),
    MANIFESTVAR = matrix("residual", 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1),
    T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(0, 1, 1)
  ))
  dat <- data.frame(id = rep(1:2, each = 3), time = rep(0:2, 2), Y1 = 0)

  prepared <- suppressMessages(ctFit(
    dat, model, backend = "julia", fit = FALSE,
    backendcontrol = list(gradient = "forward"), nlcontrol = list(maxtimestep = .5)
  ))
  expect_s3_class(prepared, "ctJuliaModel")
  expect_true(is.data.frame(prepared$parameter_table))
  expect_equal(prepared$subject_starts, c(1L, 4L))
  expect_equal(prepared$max_timestep, .5)
  expect_true(all(c("PARS", "JAx", "Jy") %in% prepared$parameter_table$matrix))
  free_transforms <- prepared$parameter_table$transform[
    !is.na(prepared$parameter_table$parnumber)
  ]
  expect_true(all(grepl("param[", free_transforms, fixed = TRUE)))
  expect_silent(unserialize(serialize(prepared, NULL)))
})

test_that("Julia setup emits ASCII string literals", {
  expect_equal(ctsem:::.ctJuliaString("https://example.org/engine.git"),
    "\"https://example.org/engine.git\"")
  expect_equal(ctsem:::.ctJuliaString("C:/Julia/project"), "\"C:/Julia/project\"")
  expect_error(ctsem:::.ctJuliaString("bad\nvalue"), "one non-empty-line")
})

test_that("Julia default initial values match Stan's small random raw initialization", {
  set.seed(20260820)
  expected <- rnorm(4, 0, .01)
  set.seed(20260820)
  expect_equal(ctsem:::.ctJuliaInitialValues(4), expected)
  expect_equal(ctsem:::.ctJuliaInitialValues(2, c(.1, -.2)), c(.1, -.2))
})

test_that("Julia backend rejects unsupported capabilities before session startup", {
  model <- ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(0, 1, 1)
  )
  # Overriding one entry at a time, rather than passing through `...`, because
  # `...` would supply a second value for a formal the defaults already name.
  unsupported <- function(...) do.call(ctsem:::.ctJuliaUnsupported,
    utils::modifyList(list(model = model, optimize = TRUE, priors = FALSE,
      intoverpop = FALSE, vb = FALSE, gendata = FALSE, stanmodeltext = NA,
      compileArgs = list(), forcerecompile = FALSE), list(...)))
  binary <- model
  binary$manifesttype <- 1L
  expect_error(unsupported(vb = TRUE), "variational Bayes")
  expect_error(unsupported(gendata = TRUE), "generation")
  expect_error(unsupported(forcerecompile = TRUE), "Stan compilation controls")
  # Binary is supported now -- the filter integrates the observation rather
  # than linearising it -- so what is refused is anything beyond it.
  expect_no_error(unsupported(model = binary))
  # Ordinal is supported too, by the same quadrature with a cumulative logit
  # in place of the Bernoulli likelihood.
  ordinal <- model
  ordinal$manifesttype <- 2L
  expect_no_error(unsupported(model = ordinal))
  # Count (3) and censored (4) are supported now as well, so what is refused
  # is anything past them.
  count <- model
  count$manifesttype <- 3L
  expect_no_error(unsupported(model = count))
  beyond <- model
  beyond$manifesttype <- 5L
  expect_error(unsupported(model = beyond), "beyond censored")

  # `optimize = FALSE` used to be on that list and is not any more: the engine
  # has its own sampler, so the combination has to pass the capability check
  # rather than be rejected by it. Asserting the absence of the old error keeps
  # the lifted limitation from quietly coming back.
  expect_no_error(unsupported(optimize = FALSE))
})

test_that("Julia parameter preparation emits state and Jacobian expressions", {
  model <- ctModel(
    type = "ct", LAMBDA = matrix("1 + eta1", 1, 1),
    DRIFT = matrix(-.2, 1, 1), DIFFUSION = matrix(.2, 1, 1),
    MANIFESTVAR = matrix(.1, 1, 1), MANIFESTMEANS = matrix(0, 1, 1),
    T0VAR = matrix(1, 1, 1), T0MEANS = matrix(0, 1, 1)
  )
  table <- ctsem:::.ctJuliaParameterTable(model)
  expect_equal(table$updatetransform[table$matrix == "LAMBDA"], "1 + state[1]")
  expect_equal(table$value[table$matrix == "JAx"], -.2)
  expect_equal(table$updatetransform[table$matrix == "Jy"], "1 + 2 * state[1]")
})

test_that("Julia reuses the canonical raw transform for shared parameters", {
  model <- suppressWarnings(ctModel(
    type = "ct", n.latent = 2, LAMBDA = diag(2),
    PARS = matrix("cross||TRUE", 1, 1),
    DRIFT = matrix(c("d11", "PARS[1,1]", "d21", "d22"), 2, 2, byrow = TRUE),
    DIFFUSION = diag(c(.2, .15)), MANIFESTVAR = diag(c(.1, .1)),
    MANIFESTMEANS = matrix(0, 2, 1), T0VAR = diag(2), T0MEANS = matrix(0, 2, 1)
  ))
  data <- data.frame(id = rep(1:2, each = 2), time = rep(0:1, 2),
    Y1 = 0, Y2 = 0)

  prepared <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
  d22 <- subset(prepared$parameter_table,
    matrix == "DRIFT" & row == 2L & col == 2L)
  pars <- subset(prepared$parameter_table,
    matrix == "PARS" & row == 1L & col == 1L)

  expect_true(is.na(pars$parnumber))
  expect_equal(d22$parnumber, 3L)
  expect_equal(d22$transform, "0 + 1 * (param[3] * 1 + 0)")
})

test_that("Julia preparation retains missing manifest values", {
  model <- ctModel(type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1), MANIFESTVAR = matrix("residual", 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(0, 1, 1))
  dat <- data.frame(id = 1, time = 0, Y1 = NA_real_)
  prepared <- suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE))
  expect_true(is.na(prepared$manifest_data[1, 1]))
})

test_that("Julia preparation receives ctFit predictor arrays and TI effect mapping", {
  model <- suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(0, 1, 1), n.TDpred = 1, TDpredNames = "dose",
    TDPREDEFFECT = matrix("impulse", 1, 1), n.TIpred = 1,
    TIpredNames = "group", tipredDefault = FALSE
  ))
  model$pars$group_effect[model$pars$param == "drift"] <- TRUE
  dat <- data.frame(
    id = rep(1:2, each = 3), time = rep(c(0, 1, 3), 2), Y1 = 0,
    dose = c(1, NA, 0, 0, 1, 0), group = rep(c(-1, 2), each = 3)
  )

  prepared <- suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE))
  expect_equal(prepared$tdpred_data, matrix(c(1, 0, 0, 0, 1, 0), nrow = 1))
  expect_equal(unname(prepared$tipred_data), matrix(c(-1, 2), ncol = 1))
  expect_equal(prepared$ti_effects$parameter, 1L)
  expect_equal(prepared$ti_effects$predictor, 1L)
  expect_equal(prepared$ti_effects$coefficient,
    max(prepared$parameter_table$parnumber, na.rm = TRUE) + 1L)
})

test_that("Julia preparation expands individual differences into static states", {
  model <- ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(0, 1, 1), CINT = matrix("cint||TRUE", 1, 1)
  )
  dat <- data.frame(id = rep(1:2, each = 2), time = rep(0:1, 2), Y1 = 0)

  prepared <- suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE))
  expect_equal(prepared$nlatent, 1L)
  expect_equal(prepared$nlatent_augmented, 2L)
  expect_equal(prepared$dynamic_state_indices, 1L)
  expect_true(any(prepared$parameter_table$matrix == "DRIFT" &
    prepared$parameter_table$row == 2L & prepared$parameter_table$col == 2L))
  expect_true(any(prepared$random_effects$type == "sd"))
  expect_true(any(prepared$rewritten_cells$matrix == "CINT"))
})

test_that("Julia preparation rejects non-constant TI predictors", {
  model <- suppressWarnings(ctModel(
    type = "ct", LAMBDA = diag(1), DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix(.2, 1, 1), MANIFESTVAR = matrix(.1, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1), T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(0, 1, 1), n.TIpred = 1, TIpredNames = "group",
    tipredDefault = FALSE
  ))
  model$pars$group_effect[model$pars$param == "drift"] <- TRUE
  dat <- data.frame(id = c(1, 1), time = 0:1, Y1 = 0, group = c(0, 1))
  expect_error(
    suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE)),
    "constant within subject"
  )
})

test_that("Julia completes a full AnomAuth optimization", {
  skip_if_not_installed("JuliaConnectoR")
  skip_if(Sys.getenv("CTSEM_RUN_JULIA_E2E", unset = "") != "true",
    "Set CTSEM_RUN_JULIA_E2E=true to run the full Julia fit regression.")

  data(AnomAuth, package = "ctsem")
  model <- ctModel(LAMBDA = diag(2), n.latent = 2, n.manifest = 2,
    MANIFESTVAR = diag(0, 2), Tpoints = 5)
  model$pars$indvarying <- FALSE
  dat <- ctDeintervalise(ctWideToLong(AnomAuth, Tpoints = model$Tpoints,
    n.manifest = 2))

  set.seed(20260820)
  fit <- suppressMessages(ctFit(dat, model, backend = "julia", optimize = TRUE,
    savescores = FALSE, cores = 1))
  expect_s3_class(fit, "ctJuliaFit")
  expect_true(is.finite(fit$estimate$loglik))
  expect_true(fit$estimate$converged)
})

# Row 1 used to build the manifest covariance before the update-group transform
# that owns the cell had run, so a state-dependent MANIFESTVAR was read from a
# parameter slot nothing had written. It needs no Stan to pin: with one
# occasion per subject the state at the measurement is T0MEANS, which is fixed
# here, so the state-dependent cell must give exactly what the same model gives
# with that cell fixed at the transform's value there.
test_that("a state-dependent MANIFESTVAR is transformed before row 1 reads it", {
  skip_on_cran()
  skip_without_julia()

  t0 <- 1.5
  fixedsd <- .3 + .05 * t0   # what the expression yields at the initial state

  .m <- function(manifestvar) suppressWarnings(ctModel(
    type = "ct",
    LAMBDA = diag(1),
    DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1),
    MANIFESTVAR = matrix(manifestvar, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1),
    T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(t0, 1, 1)))

  statedep <- .m(paste0(".3 + .05 * eta1"))
  equivalent <- .m(as.character(fixedsd))

  # One row per subject, so every row is a row 1.
  set.seed(11)
  dat <- data.frame(id = 1:8, time = 0, Y1 = stats::rnorm(8, t0, 1))

  specA <- suppressMessages(ctFit(dat, statedep, backend = "julia", fit = FALSE))
  specB <- suppressMessages(ctFit(dat, equivalent, backend = "julia", fit = FALSE))

  npar <- max(specA$parameter_table$parnumber, na.rm = TRUE)
  expect_equal(npar, max(specB$parameter_table$parnumber, na.rm = TRUE))

  set.seed(12)
  raw <- stats::rnorm(npar, 0, .3)
  a <- ctJuliaEvaluate(specA, raw, gradient = FALSE)$value
  b <- ctJuliaEvaluate(specB, raw, gradient = FALSE)$value

  expect_true(is.finite(as.numeric(a)))
  expect_equal(as.numeric(a), as.numeric(b), tolerance = 1e-10)

  # And the transform must actually be doing something: a cell fixed at a
  # different value has to disagree, or the test above would pass on a model
  # where MANIFESTVAR never varied at all.
  specC <- suppressMessages(ctFit(dat, .m(as.character(fixedsd * 3)),
    backend = "julia", fit = FALSE))
  cc <- ctJuliaEvaluate(specC, raw, gradient = FALSE)$value
  expect_false(isTRUE(all.equal(as.numeric(a), as.numeric(cc), tolerance = 1e-6)))
})

# `covmattransform` never reached the engine: every call passed a literal 0 and
# the branch the other values select is commented out, so a 'cholesky' model
# was fitted with the default transform and measured 3.76 log units away from
# the same model on stan, silently. Refused rather than implemented.
test_that("a non-default covmattransform is refused on the julia backend", {
  skip_on_cran()

  .m <- function() suppressWarnings(ctModel(
    type = "ct",
    LAMBDA = diag(1, 2),
    DRIFT = matrix(c("drift11", 0, 0, "drift22"), 2, 2),
    DIFFUSION = matrix(c("diff11", "diff21", 0, "diff22"), 2, 2),
    MANIFESTMEANS = matrix(0, 2, 1),
    T0MEANS = matrix(0, 2, 1)))
  dat <- data.frame(id = rep(1:4, each = 2), time = rep(0:1, 4),
    Y1 = 0, Y2 = 0)

  for (tf in c("cholesky", "rawcorr_indep")) {
    model <- .m()
    model$covmattransform <- tf
    expect_error(
      suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE)),
      regexp = "covmattransform")
  }

  # The default still passes the guard. `fit = FALSE` stops before Julia is
  # needed, so this half does not depend on a Julia installation.
  model <- .m()
  expect_identical(model$covmattransform, "rawcorr")
  expect_error(suppressMessages(ctFit(dat, model, backend = "julia", fit = FALSE)),
    regexp = NA)
})

# The same row-1 ordering defect one hop further out. PARS lives in the
# *predict* transform group, and row 1 ran the td and the update group only, so
# an update-group cell written from a PARS cell -- LAMBDA, MANIFESTMEANS,
# MANIFESTVAR or Jy -- read a parameter slot the predict group had not filled
# and took the buffer's zero. Against stan that was 11.5 log units at one
# occasion per subject.
#
# Pinned like its sibling above and for the same reason: at one occasion the
# state at the measurement is T0MEANS. The PARS cell here is indvarying, so
# intoverpop rewrites it to the augmented population coordinate, whose T0MEANS
# is the raw parameter under the identity transform -- a known value. The
# consumer is MANIFESTVAR rather than LAMBDA because the manifest covariance
# does not enter Jy, so the two models agree exactly rather than to a
# linearisation.
test_that("a PARS cell is transformed before row 1's update group reads it", {
  skip_on_cran()
  skip_without_julia()

  t0 <- 1.5

  .m <- function(manifestvar) suppressWarnings(ctModel(
    type = "ct",
    LAMBDA = diag(1),
    PARS = matrix("mvp||TRUE", 1, 1),
    DRIFT = matrix("drift", 1, 1),
    DIFFUSION = matrix("diffusion", 1, 1),
    MANIFESTVAR = matrix(manifestvar, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1),
    T0VAR = matrix(1, 1, 1),
    T0MEANS = matrix(t0, 1, 1)))

  # One row per subject, so every row is a row 1.
  set.seed(11)
  dat <- data.frame(id = 1:8, time = 0, Y1 = stats::rnorm(8, t0, 1))

  specA <- suppressMessages(ctFit(dat, .m("PARS[1,1]"), backend = "julia", fit = FALSE))
  ptab <- specA$parameter_table
  npar <- max(ptab$parnumber, na.rm = TRUE)

  # The two halves of "the value the transform yields at row 1": PARS[1,1] is
  # the augmented state, and that state at row 1 is its own T0MEANS, which is
  # the raw parameter untransformed.
  augmented <- which(ptab$matrix == "T0MEANS" & ptab$row == 2)
  expect_identical(ptab$predicttransform[ptab$matrix == "PARS"][1], "state[2]")
  expect_identical(ptab$transform[augmented],
    paste0("param[", ptab$parnumber[augmented], "]"))

  raw <- rep(-.5, npar)
  raw[ptab$parnumber[augmented]] <- .4

  specB <- suppressMessages(ctFit(dat, .m("0.4"), backend = "julia", fit = FALSE))
  expect_equal(npar, max(specB$parameter_table$parnumber, na.rm = TRUE))

  a <- ctJuliaEvaluate(specA, raw, gradient = FALSE)$value
  b <- ctJuliaEvaluate(specB, raw, gradient = FALSE)$value

  expect_true(is.finite(as.numeric(a)))
  expect_equal(as.numeric(a), as.numeric(b), tolerance = 1e-10)

  # And the cell has to matter, or the equality above would also hold on a
  # model whose MANIFESTVAR never varied.
  specC <- suppressMessages(ctFit(dat, .m("1.2"), backend = "julia", fit = FALSE))
  cc <- ctJuliaEvaluate(specC, raw, gradient = FALSE)$value
  expect_false(isTRUE(all.equal(as.numeric(a), as.numeric(cc), tolerance = 1e-6)))
})
