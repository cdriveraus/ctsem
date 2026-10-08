# The inner solve's choice of mode (`_laplace_probe_modes`). Where a unit's
# integrand has more than one maximum, Newton from the origin reaches whichever
# one the origin's basin holds; probes along the found mode's soft directions
# keep the candidate with the largest term, scored exactly as the evaluation
# scores it (`_laplace_unit_primal`).
#
# The fixture is a loading whose sign the data cannot fix: y = lambda * eta
# with eta centred, so a person's likelihood depends on lambda^2 and has a
# maximum either side of zero. At a population mean of zero the origin is the
# stationary point between them -- not a maximum -- which Newton accepts at
# once and the evaluation scores by the gated rule's soft direction.

.modes_sign_data <- function() {
  set.seed(3)
  n <- 12; nt <- 8
  d <- data.frame(id = rep(1:n, each = nt), time = rep(0:(nt - 1), n))
  lambda <- 1.5 * sample(c(-1, 1), n, TRUE)
  eta <- unlist(lapply(1:n, function(i) {
    e <- numeric(nt)
    e[1] <- rnorm(1)
    for (t in 2:nt) e[t] <- 0.6 * e[t - 1] + rnorm(1, 0, 0.8)
    e
  }))
  d$Y1 <- lambda[d$id] * eta + rnorm(n * nt, 0, 0.3)
  d
}

.modes_sign_model <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", id = "id",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix("lambda"),
    DRIFT = matrix(-0.5), DIFFUSION = matrix(0.8), MANIFESTVAR = matrix(0.3),
    MANIFESTMEANS = matrix(0), T0MEANS = matrix(0), T0VAR = matrix(1),
    CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% "lambda"
  m$pars$transform[m$pars$param %in% "lambda"] <- "param"
  m
}

test_that("the inner solve keeps the best mode it finds, scored as the evaluation scores it", {
  skip_without_julia()
  spec <- suppressMessages(ctFit(.modes_sign_data(), .modes_sign_model(),
    backend = "julia", intoverpop = "laplace", fit = FALSE))
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  objective <- ctsem:::.ctJuliaObjective(handle, hybrid = TRUE)
  probes <- function(on) invisible(ctsem:::.ctJuliaCall(
    "ContinuousTimeSEM.ctsem_set_mode_probes!", on))
  on.exit(probes(TRUE), add = TRUE)
  evaluate <- function(x) ctJuliaEvaluate(handle, x, gradient = TRUE,
    contributions = TRUE, cores = 1L)
  # Each unit's term at the modes the last evaluation left, by the function the
  # probes compare with: it must be the evaluation's own, the gated rule at a
  # curvature that does not factor included.
  ctsem:::.ctJuliaEval(paste(
    "function _ctsem_test_unit_terms(o, theta)",
    "  C = ContinuousTimeSEM",
    "  Ls = C._laplace_popchols(theta, o.spec)",
    "  aws = C._laplace_workspace!(o, Float64, length(theta))",
    "  [C._laplace_unit_term(o, U, theta, Ls, copy(o.modes[U]), aws)",
    "    for U in eachindex(o.units.members)]",
    "end", sep = "\n"))
  terms <- function(x) as.numeric(ctsem:::.ctJuliaCall("Main._ctsem_test_unit_terms",
    objective, ctsem:::.ctJuliaNumericVector(x)))

  # Population mean zero, raw sd 0.5: every origin a stationary point.
  x <- c(0, 0.5)
  probes(FALSE)
  off <- evaluate(x)
  expect_equal(terms(x), as.numeric(off$subject_loglik), tolerance = 1e-10)
  probes(TRUE)
  on <- evaluate(x)
  expect_equal(terms(x), as.numeric(on$subject_loglik), tolerance = 1e-10)
  gain <- as.numeric(on$subject_loglik) - as.numeric(off$subject_loglik)
  expect_true(all(gain >= -1e-10))
  expect_gt(max(gain), 0.1)

  # Away from zero the origin lies in a basin, and the probes' choice is a
  # smooth function there: the gradient is the value's.
  x <- c(0.3, 0.5)
  r <- evaluate(x)
  set.seed(1)
  v <- rnorm(2)
  v <- v / sqrt(sum(v^2))
  h <- 1e-5
  fd <- (evaluate(x + h * v)$value - evaluate(x - h * v)$value) / (2 * h)
  expect_equal(fd, sum(as.numeric(r$gradient) * v), tolerance = 1e-5)
})
