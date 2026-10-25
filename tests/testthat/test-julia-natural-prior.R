# control$flat = 'natural' (R/ctNaturalPrior.R): the parameters without a
# prior flat on their natural scale in the density 'checked' and 'is'
# integrate over, by the change of variable log |f'| of each one's own
# transform.

.np_model <- function(drift = NULL) {
  args <- list(type = "ct", n.latent = 1, n.manifest = 2, manifestNames = c("y1", "y2"),
    latentNames = "eta1", LAMBDA = matrix(1, 2, 1), MANIFESTMEANS = matrix(0, 2, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0),
    MANIFESTVAR = matrix(c("mvar", "0", "0", "mvar"), 2, 2))
  if (!is.null(drift)) args$DRIFT <- matrix(drift)
  m <- suppressWarnings(suppressMessages(do.call(ctModel, args)))
  m$pars$indvarying <- m$pars$param %in% "cint"
  m
}
.np_data <- function() {
  set.seed(1)
  do.call(rbind, lapply(1:12, function(i) {
    cint <- stats::rnorm(1, 0, .5); eta <- numeric(6); eta[1] <- stats::rnorm(1)
    for (t in 2:6) eta[t] <- .74 * eta[t - 1] + .86 * cint + stats::rnorm(1, 0, .6)
    data.frame(id = i, time = 0:5, y1 = eta + stats::rnorm(6, 0, .5),
      y2 = eta + stats::rnorm(6, 0, .5))
  }))
}
.np_env <- function(x) list(param = x, inv_logit = stats::plogis,
  log1p_exp = function(z) ifelse(z > 30, z, log1p(exp(z))))

test_that("the engine's natural-scale prior is the terms R writes, on both routes", {
  skip_without_julia()
  handle <- function(s) structure(s, class = c("ctJuliaModel", "ctFitModel"))
  d <- .np_data()
  for (route in c("augmented", "laplace")) {
    spec <- suppressWarnings(suppressMessages(ctFit(d, .np_model(), backend = "julia",
      intoverpop = route, fit = FALSE, priors = FALSE)))
    nat <- ctsem:::.ctJuliaNaturalPrior(spec)
    names <- as.character(ctsem:::.ctBackendRawParameterNames(list(model_spec = spec),
      ctsem:::.ctBackendNpar(spec)))
    # Every sd and the drift, the population sd included (a RAWPOPVAR cell on
    # the augmented route, the level's own transform on Laplace); not the
    # intercept, whose transform is linear.
    expect_setequal(names[nat$index], c("drift_eta1", "diff_eta1", "mvar", "T0var_eta1",
      "popsd_cint"))
    spec2 <- spec
    spec2$naturalprior <- nat
    set.seed(3)
    x <- stats::rnorm(ctsem:::.ctBackendNpar(spec), 0, 0.7)
    e1 <- ctJuliaEvaluate(handle(spec), x, cores = 1L)
    e2 <- ctJuliaEvaluate(handle(spec2), x, cores = 1L)
    value <- sum(vapply(nat$value, function(t) eval(str2lang(t), .np_env(x)), numeric(1)))
    g <- numeric(length(x))
    g[nat$index] <- vapply(nat$slope, function(t) eval(str2lang(t), .np_env(x)), numeric(1))
    expect_equal(e2$value - e1$value, value, tolerance = 1e-10, info = route)
    expect_equal(as.numeric(e2$gradient) - as.numeric(e1$gradient), g, tolerance = 1e-10,
      info = route)
    module <- ctsem:::.ctJuliaModule(spec$project)
    h <- function(s) as.matrix(ctsem:::.ctBackendJuliaValue(module$ctsem_hessian(
      ctsem:::.ctJuliaObjective(handle(s), hybrid = TRUE), ctsem:::.ctJuliaNumericVector(x))))
    C <- matrix(0, length(x), length(x))
    C[cbind(nat$index, nat$index)] <- ctsem:::.ctNaturalCurvature(nat, x)
    expect_equal(h(spec2) - h(spec), C, tolerance = 1e-10, info = route)
  }
})

test_that("each parameter's Jacobian is its own transform's, whatever the user wrote", {
  skip_without_julia()
  d <- .np_data()
  term <- function(drift) {
    spec <- suppressWarnings(suppressMessages(ctFit(d, .np_model(drift), backend = "julia",
      intoverpop = "augmented", fit = FALSE, priors = FALSE)))
    nat <- ctsem:::.ctJuliaNaturalPrior(spec)
    k <- spec$parameter_table$parnumber[spec$parameter_table$matrix %in% "DRIFT"][1]
    nat[nat$index == k, , drop = FALSE]
  }
  # A drift free to be positive, written linearly: no Jacobian, so nothing
  # changes for it.
  expect_equal(nrow(term("drift|param")), 0L)
  # Negative by an exponential: log |f'| is the raw value itself.
  neg <- term("drift|-exp(param)")
  expect_equal(nrow(neg), 1L)
  expect_equal(eval(str2lang(neg$slope), .np_env(rep(0.3, 10))), 1)
})

test_that("'is' with flat = 'natural' samples the natural-flat posterior of a closed form", {
  skip_without_julia()
  # y ~ N(0, s^2), the one free parameter the measurement sd: under a flat
  # prior on s the posterior is s^-n exp(-S / (2 s^2)); under the raw prior
  # it is that over f'(raw). Both by quadrature over raw.
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "y", latentNames = "eta", LAMBDA = matrix(1), DRIFT = matrix(-1),
    DIFFUSION = matrix(0), CINT = matrix(0), T0MEANS = matrix(0), T0VAR = matrix(0),
    MANIFESTMEANS = matrix(0), MANIFESTVAR = matrix("mvar"))))
  set.seed(11)
  n <- 8
  d <- data.frame(id = rep(1:2, each = n / 2), time = rep(seq_len(n / 2) - 1, 2),
    y = stats::rnorm(n, 0, 0.7))
  S <- sum(d$y^2)
  set.seed(2)
  fit <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia", cores = 1,
    verbose = 0, priors = FALSE)))
  expect_equal(as.character(fit$model_spec$parameter_table$transform[
    fit$model_spec$parameter_table$param %in% "mvar"][1]), "1e-10 + 5 * log1p_exp(2 * param[1])")
  f <- function(raw) 1e-10 + 5 * ifelse(2 * raw > 30, 2 * raw, log1p(exp(2 * raw)))
  grid <- seq(-12, 6, length.out = 100001)
  postmean <- function(logdens) {
    w <- exp(logdens - max(logdens))
    sum(w * f(grid)) / sum(w)
  }
  loglik <- -n * log(f(grid)) - S / (2 * f(grid)^2)
  natural <- postmean(loglik + log(10 * stats::plogis(2 * grid)))
  raw <- postmean(loglik)
  expect_gt(natural - raw, 0.05)
  draw <- function(flat) {
    set.seed(1)
    suppressWarnings(suppressMessages(ctFitUncertainty(fit, "is", cores = 1,
      control = list(flat = flat, isESS = 2000))))
  }
  nat <- draw("natural")
  rw <- draw("raw")
  expect_lt(abs(mean(f(nat$estimate$rawposterior[, 1])) - natural), 0.025)
  expect_lt(abs(mean(f(rw$estimate$rawposterior[, 1])) - raw), 0.025)
  # The fit keeps its objective's curvature: the prior's is the sampled
  # density's, and the certification, the continuation and a later call read
  # the stored matrix.
  expect_equal(nat$uncertainty$hessian, rw$uncertainty$hessian)
})

test_that("R's log1p_exp is finite at any argument, as the transforms it evaluates need", {
  x <- c(-800, -40, -3, 0, 3, 40, 800)
  expect_equal(ctsem:::log1p_exp(x[3:5]), log1p(exp(x[3:5])), tolerance = 1e-15)
  expect_equal(ctsem:::log1p_exp(c(-800, 0, 40, 800)), c(0, log(2), 40, 800))
  expect_true(all(is.finite(ctsem:::log1p_exp(x))))
})

test_that("flat = 'natural' is refused where nothing integrates over a density", {
  expect_error(ctsem:::.ctBackendUncertainty(structure(list(estimate = list(raw = 1),
    args = list(resolved = list(intoverstates = TRUE))), class = "ctJuliaFit"),
    uncertainty = "hessian", draws = "normal", finishsamples = 10, cores = 1,
    control = list(flat = "natural"), verbose = 0), "applies to uncertainty")
  expect_error(ctsem:::.ctBackendUncertainty(structure(list(estimate = list(raw = 1),
    args = list(resolved = list(intoverstates = TRUE))), class = "ctJuliaFit"),
    uncertainty = "checked", draws = "imis", finishsamples = 10, cores = 1,
    control = list(flat = "log"), verbose = 0), "must be 'raw' or 'natural'")
})
