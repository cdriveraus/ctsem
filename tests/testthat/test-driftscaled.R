# ctModel(driftScaled = ...) (R/ctModelDriftScaled.R): a process's diffusion
# sd and intercept parameterised by what they accrue over one time unit given
# its own drift diagonal. The same model in other coordinates, built from PARS
# cells and expressions.

.ds_model <- function(scaled = character(), diffusion = "auto", nlat = 1L,
  varying = "cint") {
  m <- if (nlat == 1L) {
    suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 2,
      manifestNames = c("y1", "y2"), latentNames = "eta1", LAMBDA = matrix(1, 2, 1),
      MANIFESTMEANS = matrix(0, 2, 1), CINT = matrix("cint"), T0MEANS = matrix(0),
      MANIFESTVAR = matrix(c("mvar", "0", "0", "mvar"), 2, 2), driftScaled = scaled)))
  } else {
    suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
      manifestNames = c("y1", "y2"), latentNames = c("a", "b"), LAMBDA = diag(2),
      MANIFESTMEANS = matrix(0, 2, 1), CINT = matrix(c("cint_a", "cint_b"), 2),
      DIFFUSION = diffusion, MANIFESTVAR = diag(0.3, 2), driftScaled = scaled)))
  }
  m$pars$indvarying <- m$pars$param %in% varying
  m
}
.ds_data <- function() {
  set.seed(1)
  do.call(rbind, lapply(1:12, function(i) {
    cint <- stats::rnorm(1, 0, .5); eta <- numeric(6); eta[1] <- stats::rnorm(1)
    for (t in 2:6) eta[t] <- .74 * eta[t - 1] + .86 * cint + stats::rnorm(1, 0, .6)
    data.frame(id = i, time = 0:5, y1 = eta + stats::rnorm(6, 0, .5),
      y2 = eta + stats::rnorm(6, 0, .5))
  }))
}

test_that("driftScaled moves each free parameter into PARS and writes its cell as an expression", {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 2,
    manifestNames = c("y1", "y2"), latentNames = c("a", "b"), LAMBDA = diag(2),
    CINT = matrix(c("cint_a", "cint_b"), 2), PARS = matrix(c("p1", "p2"), 1, 2),
    TIpredNames = "z", driftScaled = c("DIFFUSION", "CINT"))))
  p <- m$pars
  cell <- function(mat, r, c) p[p$matrix == mat & p$row == r & p$col == c, , drop = FALSE]
  expect_equal(m$driftScaled, c("DIFFUSION", "CINT"))
  # The user's own PARS cells are untouched; the moved ones follow them, the
  # rows padded to its width.
  expect_equal(cell("PARS", 1, 1)$param, "p1")
  expect_equal(cell("PARS", 1, 2)$param, "p2")
  moved <- p[p$matrix == "PARS" & p$row > 1 & p$col == 1, ]
  expect_setequal(moved$param, c("diff_a", "diff_b", "cint_a", "cint_b"))
  pad <- p[p$matrix == "PARS" & p$row > 1 & p$col == 2, ]
  expect_true(all(is.na(pad$param)) && all(pad$value == 0))
  # Name, transform, random effect and predictor effect move with the parameter.
  expect_equal(moved$transform[moved$param == "diff_a"], "10 * log1p_exp(2 * param) + 1e-10")
  expect_true(moved$indvarying[moved$param == "cint_a"])
  expect_true(as.logical(moved$z_effect[moved$param == "cint_a"]))
  # The cells read their PARS cell and their own drift diagonal.
  k <- moved$row[moved$param == "diff_b"]
  expect_equal(cell("DIFFUSION", 2, 2)$param,
    sprintf("PARS[%d,1] * sqrt(2 * DRIFT[2,2] / expm1(2 * DRIFT[2,2]))", k))
  k <- moved$row[moved$param == "cint_a"]
  expect_equal(cell("CINT", 1, 1)$param, sprintf("PARS[%d,1] * DRIFT[1,1] / expm1(DRIFT[1,1])", k))
  expect_false(cell("CINT", 1, 1)$indvarying)
  # Correlations stay where they are.
  expect_equal(cell("DIFFUSION", 2, 1)$param, "diff_b_a")
  # And the population sds keep their names.
  expect_true("popsd_cint_a" %in% diag(m$RAWPOPVAR))

  expect_error(ctModel(type = "dt", n.latent = 1, n.manifest = 1, LAMBDA = diag(1),
    driftScaled = "CINT"), "continuous-time")
  expect_error(ctModel(type = "ct", n.latent = 1, n.manifest = 1, LAMBDA = diag(1),
    driftScaled = "T0VAR"), "DIFFUSION")
})

test_that("a drift-scaled model is the default model in other coordinates", {
  skip_without_julia()
  d <- .ds_data()
  natural <- function(fit) {
    e <- ctBackendParMatrices(fit)
    c(loglik = fit$estimate$loglik, DRIFT = e$DRIFT[1, 1],
      DIFFUSIONcov = e$DIFFUSIONcov[1, 1], CINT = e$CINT[1, 1])
  }
  for (route in c("augmented", "laplace")) {
    fits <- lapply(list(character(), c("DIFFUSION", "CINT")), function(s) {
      set.seed(7)
      suppressWarnings(suppressMessages(ctFit(d, .ds_model(s), backend = "julia",
        intoverpop = route, cores = 1, verbose = 0)))
    })
    expect_equal(natural(fits[[2]]), natural(fits[[1]]), tolerance = 1e-4, info = route)
  }
})

test_that("a scaled diffusion starts where the default's does", {
  skip_without_julia()
  d <- .ds_data()
  spec <- suppressWarnings(suppressMessages(ctFit(d, .ds_model("DIFFUSION"), backend = "julia",
    intoverpop = "augmented", fit = FALSE)))
  start <- ctsem:::.ctDataStart(d, spec$model, spec, ctsem:::.ctBackendNpar(spec))
  pt <- spec$parameter_table
  k <- pt$parnumber[pt$matrix %in% "PARS" & pt$param %in% "diff_eta1"][1]
  expect_true(is.finite(start[k]))
})

test_that("a diffusion correlation fixed at zero stays zero in continuous time", {
  skip_without_julia()
  set.seed(2)
  d <- data.frame(id = rep(1:10, each = 8), time = rep(0:7, 10),
    y1 = stats::rnorm(80), y2 = stats::rnorm(80))
  m <- .ds_model(c("DIFFUSION", "CINT"), matrix(c("diff_a", 0, 0, "diff_b"), 2, 2),
    nlat = 2L, varying = character())
  set.seed(7)
  fit <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia", cores = 1,
    verbose = 0, optimcontrol = list(estonly = TRUE))))
  expect_equal(ctBackendParMatrices(fit)$DIFFUSIONcov[2, 1], 0)
})

test_that("a scaled diffusion needs a construction whose diagonal is an sd", {
  skip_without_julia()
  m <- .ds_model("DIFFUSION")
  m$covmattransform <- "cholesky"
  expect_error(suppressWarnings(suppressMessages(ctFit(.ds_data(), m, backend = "julia",
    fit = FALSE))), "covmattransform")
})

test_that("a drift-scaled random intercept keeps the Laplace units concave quadratics", {
  skip_without_julia()
  # Its effect sits on a PARS cell that only CINT reads, affine in it given the
  # drift, so each unit is still one Gaussian basin and the probes are skipped
  # (`.ctJuliaParsAffine()`); a random drift beside it makes the slope vary
  # with an effect, and they run again. Skipping changes nothing: values and
  # gradients as with the probes on.
  d <- .ds_data()
  spec <- suppressWarnings(suppressMessages(ctFit(d, .ds_model(c("DIFFUSION", "CINT")),
    backend = "julia", intoverpop = "laplace", fit = FALSE)))
  expect_true(isTRUE(spec$laplace$quadratic))
  both <- suppressWarnings(suppressMessages(ctFit(d, .ds_model(c("DIFFUSION", "CINT"),
    varying = c("cint", "drift_eta1")), backend = "julia", intoverpop = "laplace",
    fit = FALSE)))
  expect_false(isTRUE(both$laplace$quadratic))
  probed <- spec
  probed$laplace$quadratic <- FALSE
  handles <- lapply(list(spec, probed), function(s)
    structure(s, class = c("ctJuliaModel", "ctFitModel")))
  set.seed(9)
  for (i in 1:2) {
    x <- stats::rnorm(ctsem:::.ctBackendNpar(spec), 0, 0.5)
    a <- ctJuliaEvaluate(handles[[1]], x, gradient = TRUE, cores = 1L)
    b <- ctJuliaEvaluate(handles[[2]], x, gradient = TRUE, cores = 1L)
    expect_equal(a$value, b$value, tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-10)
  }
})
