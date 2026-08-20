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
  dat <- data.frame(id = 1, time = 0, Y1 = 0)
  expect_error(
    suppressMessages(ctFit(dat, model, backend = "julia", optimize = FALSE)),
    "optimize=FALSE"
  )
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

test_that("Julia completes the AnomAuth initial optimization step", {
  skip_if_not_installed("JuliaConnectoR")
  project <- Sys.getenv("CTSEM_JULIA_PROJECT", unset = "")
  skip_if(!nzchar(project) || !dir.exists(project),
    "Set CTSEM_JULIA_PROJECT to run Julia end-to-end fit tests.")

  data(AnomAuth, package = "ctsem")
  model <- ctModel(LAMBDA = diag(2), n.latent = 2, n.manifest = 2,
    MANIFESTVAR = diag(0, 2), Tpoints = 5)
  model$pars$indvarying <- FALSE
  dat <- ctDeintervalise(ctWideToLong(AnomAuth, Tpoints = model$Tpoints,
    n.manifest = 2))

  set.seed(20260820)
  fit <- suppressMessages(ctFit(dat, model, backend = "julia", optimize = TRUE,
    savescores = FALSE, cores = 1,
    backendcontrol = list(julia_project = project, maxiter = 2)))
  expect_s3_class(fit, "ctJuliaFit")
  expect_true(is.finite(fit$estimate$loglik))
  expect_gte(fit$estimate$iterations, 1L)
})
