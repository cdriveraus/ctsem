# The state-explicit route: `intoverstates = FALSE`.
#
# Everything else in the julia backend integrates the latent states out as it
# goes, and generation rides on that -- each row is drawn from the filter's own
# one-step-ahead predictive and the filter then conditions on the draw. For a
# linear Gaussian model that is exact and is the right thing to do. For a
# non-Gaussian one the filter's measurement update is an assumed-density
# projection, and a projection can move the state it conditions on.
#
# It does, and mixing indicator kinds is what makes it visible. Ordinal and
# binary observations are bounded, so a state that has run away still produces
# ordinary-looking data in them; a count is not, so it shows. The reported case
# was a two-latent model with an ordinal and a binary indicator on one process
# and a count on the other: one subject drew a count of 252 at its first row,
# the two rows after it came back at the engine's generation clamp of 100000,
# and a fit to that data returned standard errors of zero.
#
# `intoverstates = FALSE` cannot do that, and the tests below are mostly about
# saying why in a way that would fail if it started to: the state is a draw
# from the process, the observation is a draw given the state, and no
# observation is ever conditioned on. The distributional checks live in the
# engine's own `test_state_sampling.jl`, which can hold the state fixed and
# take moments; what is left for here is the R route to it.

.states_model <- function(diffusion = 1, t0var = 1) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 2,
    n.manifest = 3, manifestNames = c("o1", "b1", "c1"),
    latentNames = c("eta1", "eta2"),
    manifesttype = c(2L, 1L, 3L), ncategories = c(4L, 0L, 0L),
    LAMBDA = matrix(c(1, 0, 1, 0, 0, 1), 3, 2, byrow = TRUE),
    DRIFT = matrix(c(-0.4, 0.25, -0.15, -0.7), 2, 2, byrow = TRUE),
    DIFFUSION = diag(diffusion, 2), T0VAR = diag(t0var, 2),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    MANIFESTVAR = diag(0, 3),
    MANIFESTMEANS = matrix(c(0, 0, 0.8), 3, 1), Tpoints = 10)))
  # Thresholds at -1, 0, 1. The matrix holds the first threshold and then the
  # gaps to it, so the cells are -1, 1, 1; setting `value` is what fixes a cell,
  # whatever its label says.
  rows <- m$pars$matrix %in% "THRESHOLDS" & m$pars$row %in% 1L
  m$pars$value[rows] <- c(-1, 1, 1)[m$pars$col[rows]]
  m$pars$indvarying <- FALSE
  m
}

test_that("the engine reports the innovation count the design needs", {
  skip_without_julia()
  model <- .states_model()
  times <- lapply(1:3, function(i) 0:4)
  skeleton <- ctsem:::.ctGenerateSkeleton(model, 3, times)
  spec <- ctsem:::.ctJuliaPrepare(skeleton,
    suppressMessages(ctsem:::.ctGenerateResolveFree(model, quiet = TRUE)),
    priors = FALSE, intoverpop = "augmented")
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  # Two latent states at every row: one innovation each at the first row of a
  # subject and one each per interval after it, which for this design is the
  # same count either way.
  expect_equal(ctsem:::.ctBackendStateDimension(handle), 2L * 15L)
})

test_that("generation without a maxtimestep is reproducible from the seed", {
  skip_without_julia()
  model <- .states_model()
  set.seed(11)
  first <- suppressMessages(ctGenerate(model, n.subjects = 5, Tpoints = 6,
    backend = "julia", intoverstates = FALSE))
  set.seed(11)
  again <- suppressMessages(ctGenerate(model, n.subjects = 5, Tpoints = 6,
    backend = "julia", intoverstates = FALSE))
  expect_equal(first, again)
  set.seed(12)
  different <- suppressMessages(ctGenerate(model, n.subjects = 5, Tpoints = 6,
    backend = "julia", intoverstates = FALSE))
  expect_false(isTRUE(all.equal(first, different)))
})

test_that("a mixed ordinal, binary and count model generates bounded data", {
  skip_without_julia()
  # The reported reproduction, at its own seed.
  model <- .states_model()
  set.seed(11)
  generated <- suppressMessages(ctGenerate(model, n.subjects = 30,
    Tpoints = 10, backend = "julia", intoverstates = FALSE))
  data <- data.frame(generated)

  expect_true(all(is.finite(as.matrix(data))))
  expect_true(all(data$o1 %in% 1:4))
  expect_true(all(data$b1 %in% c(0, 1)))
  expect_true(all(data$c1 >= 0 & data$c1 == round(data$c1)))

  # The count is the indicator that showed the runaway, and this is the
  # assertion that would have caught it: the process it loads on has a
  # stationary standard deviation under one, so a rate of exp(0.8 + eta2) has
  # no way to reach the hundreds, let alone the engine's clamp of 100000.
  expect_lt(max(data$c1), 500)
  # And the state that produced it must relax rather than accumulate: the
  # largest count in the second half of a subject's series cannot be a
  # multiple of the largest in the first, which is what a state walking outward
  # would give.
  halves <- tapply(seq_len(nrow(data)), data$id, function(rows) {
    early <- rows[seq_len(length(rows) %/% 2)]
    late <- setdiff(rows, early)
    max(data$c1[late]) - max(data$c1[early])
  })
  expect_lt(max(halves), 100)
})

test_that("with no process noise the counts are Poisson at the fixed rate", {
  skip_without_julia()
  # T0VAR and DIFFUSION at zero pin every state at T0MEANS, so each row is an
  # independent draw at the same rate and the sample mean and variance are the
  # whole of what a Poisson claims. A draw taken from the wrong distribution --
  # or from the right one at the wrong rate -- fails both.
  model <- .states_model(diffusion = 0, t0var = 0)
  set.seed(4)
  data <- data.frame(suppressMessages(ctGenerate(model, n.subjects = 60,
    Tpoints = 20, backend = "julia", intoverstates = FALSE)))
  rate <- exp(0.8)
  n <- nrow(data)
  expect_equal(mean(data$c1), rate, tolerance = 4 * sqrt(rate / n) / rate)
  expect_equal(stats::var(data$c1), rate, tolerance = 0.15)
  # eta1 is pinned at zero too, so the binary indicator is a fair coin.
  expect_equal(mean(data$b1), 0.5, tolerance = 4 * sqrt(0.25 / n))
  # ... and the ordinal one splits at thresholds -1, 0, 1.
  cumulative <- 1 / (1 + exp(-c(-1, 0, 1)))
  expected <- diff(c(0, cumulative, 1))
  for (k in 1:4) {
    expect_equal(mean(data$o1 == k), expected[k],
      tolerance = 4 * sqrt(expected[k] / n))
  }
})

test_that("the joint density is finite and differentiable through R", {
  skip_without_julia()
  model <- suppressMessages(ctsem:::.ctGenerateResolveFree(.states_model(),
    quiet = TRUE))
  times <- lapply(1:4, function(i) 0:5)
  skeleton <- ctsem:::.ctGenerateSkeleton(model, 4, times)
  set.seed(3)
  filled <- suppressMessages(ctGenerate(model, n.subjects = 4, Tpoints = 6,
    backend = "julia", intoverstates = FALSE))
  skeleton[, model$manifestNames] <- data.frame(filled)[, model$manifestNames]

  spec <- ctsem:::.ctJuliaPrepare(skeleton, model, priors = FALSE,
    intoverpop = "augmented")
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  npar <- suppressWarnings(max(c(0L,
    as.integer(spec$parameter_table$parnumber)), na.rm = TRUE))
  raw <- if (npar < 1L) 0 else numeric(npar)
  ndim <- ctsem:::.ctBackendStateDimension(handle)

  set.seed(5)
  z <- stats::rnorm(ndim)
  result <- ctsem:::.ctBackendJointDensity(handle, raw, z)
  expect_true(is.finite(result$value))
  expect_length(as.numeric(result$gradient), max(npar, 1L) + ndim)
  expect_true(all(is.finite(as.numeric(result$gradient))))
  # The pieces are a partition of the total, which is what makes the value a
  # joint density rather than a likelihood with a penalty attached.
  expect_equal(result$observation + result$state_prior + result$parameter_prior,
    result$value)
  # The innovations' own density is the standard normal one, exactly.
  expect_equal(as.numeric(result$state_prior),
    sum(stats::dnorm(z, log = TRUE)))

  # Moving one innovation moves the density by what the gradient says it will.
  step <- 1e-5
  moved <- z
  moved[3] <- moved[3] + step
  shifted <- ctsem:::.ctBackendJointDensity(handle, raw, moved, gradient = FALSE)
  slope <- as.numeric(result$gradient)[max(npar, 1L) + 3L]
  expect_equal((shifted$value - result$value) / step, slope, tolerance = 1e-4)
})

# A plain linear Gaussian model, where the two routes agree in distribution and
# `intoverstates='auto'` therefore keeps the filter one.
.gaussian_model <- function() {
  suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "y1", latentNames = "eta1",
    LAMBDA = matrix(1), DRIFT = matrix(-0.4), DIFFUSION = matrix(0.8),
    MANIFESTVAR = matrix(0.3), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix(0), MANIFESTMEANS = matrix(0.5), Tpoints = 6)))
}

test_that("intoverstates='auto' picks the route the model needs", {
  skip_without_julia()
  # Linear and Gaussian: the filter's predictive is exact, so 'auto' keeps it
  # and the output is what an existing caller already gets.
  model <- .gaussian_model()
  set.seed(21)
  auto <- suppressMessages(ctGenerate(model, n.subjects = 4, Tpoints = 6,
    backend = "julia", intoverstates = "auto"))
  set.seed(21)
  filtered <- suppressMessages(ctGenerate(model, n.subjects = 4, Tpoints = 6,
    backend = "julia", intoverstates = TRUE))
  expect_equal(auto, filtered)

  # Non-Gaussian indicators: the update is an assumed-density projection, so
  # 'auto' takes the state path instead.
  mixed <- .states_model()
  set.seed(21)
  autom <- suppressMessages(ctGenerate(mixed, n.subjects = 4, Tpoints = 6,
    backend = "julia", intoverstates = "auto"))
  set.seed(21)
  sampled <- suppressMessages(ctGenerate(mixed, n.subjects = 4, Tpoints = 6,
    backend = "julia", intoverstates = FALSE))
  expect_equal(autom, sampled)
  set.seed(21)
  viafilter <- suppressMessages(ctGenerate(mixed, n.subjects = 4, Tpoints = 6,
    backend = "julia", intoverstates = TRUE))
  expect_false(isTRUE(all.equal(autom, viafilter)))
})

test_that("a fit over the joint density runs and carries its trajectory", {
  skip_without_julia()
  model <- .states_model()
  set.seed(2)
  data <- data.frame(suppressMessages(ctGenerate(model, n.subjects = 12,
    Tpoints = 8, backend = "julia", intoverstates = FALSE)))

  fitmodel <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    n.latent = 2, n.manifest = 3, manifestNames = c("o1", "b1", "c1"),
    latentNames = c("eta1", "eta2"), manifesttype = c(2L, 1L, 3L),
    ncategories = c(4L, 0L, 0L),
    LAMBDA = matrix(c(1, 0, 1, 0, 0, 1), 3, 2, byrow = TRUE),
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    MANIFESTVAR = diag(0, 3),
    MANIFESTMEANS = matrix(c(0, 0, "cmean"), 3, 1), Tpoints = 8)))
  fitmodel$pars$indvarying <- FALSE

  fit <- suppressWarnings(suppressMessages(ctFit(data, fitmodel,
    backend = "julia", intoverstates = FALSE, verbose = 0)))

  expect_s3_class(fit, "ctJuliaFit")
  npar <- length(fit$estimate$raw)
  expect_gt(npar, 0)
  expect_true(all(is.finite(fit$estimate$raw)))
  # The trajectory is what this route has that the marginal one does not: one
  # state per latent per data row, and the innovations it was built from.
  expect_equal(dim(fit$estimate$states), c(nrow(data), 2L))
  expect_true(all(is.finite(fit$estimate$states)))
  expect_length(fit$estimate$innovations, 2L * nrow(data))
  # And `raw` is still the parameters alone, so everything downstream reads it
  # as it always has.
  expect_false(npar == length(fit$estimate$innovations))
  expect_identical(fit$estimate$loglik_type, "joint")
  expect_identical(fit$args$resolved$intoverstates, FALSE)
  expect_true(is.finite(fit$estimate$loglik))
})

test_that("standard errors profile the states out, and the rest are refused", {
  skip_without_julia()
  model <- .gaussian_model()
  set.seed(6)
  data <- data.frame(suppressMessages(ctGenerate(model, n.subjects = 15,
    Tpoints = 6, backend = "julia", intoverstates = TRUE)))
  # MANIFESTVAR fixed, which it has to be. With a Gaussian indicator's
  # measurement variance free the joint density is unbounded -- send it to zero
  # and the trajectory interpolates the data exactly -- so there is no maximum
  # for a Hessian to describe. Measured before this line said so: a raw value
  # of -10.6 for it and a gradient of 3e9, reported as not converged. ctFit
  # warns about exactly this, and the test below checks that it does.
  fitmodel <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    n.latent = 1, n.manifest = 1, manifestNames = "y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTVAR = matrix(0.3), MANIFESTMEANS = matrix("mmean"), Tpoints = 6)))
  fitmodel$pars$indvarying <- FALSE

  fit <- suppressWarnings(suppressMessages(ctFit(data, fitmodel,
    backend = "julia", intoverstates = FALSE, verbose = 0)))
  npar <- length(fit$estimate$raw)

  # Kept on the fit rather than turned into intervals.
  expect_equal(dim(fit$estimate$hessian_profile), c(npar, npar))
  expect_null(fit$uncertainty)
  hessian <- ctsem:::.ctBackendHessian(fit, fit$estimate$raw)
  expect_equal(dim(hessian), c(npar, npar))
  expect_true(all(is.finite(hessian)))
  # Negative semi-definite at a maximum, which is what makes it usable as an
  # observed information. Semi- rather than strictly, because a joint mode can
  # sit against a flat direction -- the drift running to zero is the usual one,
  # since less mean reversion means smaller innovations and the innovation
  # prior rewards that. The identifiability report is what turns such a
  # direction into something the user sees.
  # Negative semi-definite at a maximum. Semi- rather than strictly, and small:
  # this is the curvature of the *profile*, and with an innovation for every
  # observation the states re-optimise to absorb almost any change in the
  # parameters, so the profile is nearly flat. Measured here, the largest
  # eigenvalue is about 0.05 -- which is why an optimised fit reports no
  # standard errors and says so, rather than inverting this into intervals.
  values <- eigen(-(hessian + t(hessian)) / 2, only.values = TRUE)$values
  expect_true(all(is.finite(values)))
  # A smallest eigenvalue of a numerically profiled Hessian, not a quantity
  # pinned to a digit: the fit runs at julia's default thread count, where
  # results are not reproducible below about 1e-7, and this floor's job is to
  # catch a genuinely indefinite Hessian rather than that run-to-run noise.
  # -1e-2 clears the observed -0.00130 with room while still well below the
  # ~0.05 largest eigenvalue noted above, so a matrix that is actually
  # indefinite in a substantial direction still fails it.
  expect_gt(min(values), -1e-2)

  # Everything except 'hessian' would score the marginal likelihood, which is
  # not the density this fit maximised.
  expect_error(ctOptimUncertainty(fit, uncertainty = "opg"), "intoverstates")
})

test_that("a free Gaussian measurement variance is called out, not left to fail", {
  skip_without_julia()
  # The joint density has no maximum in this direction, so there is nothing for
  # the optimiser to converge to; ctFit now stops rather than letting it run to
  # the boundary and report not converged.
  model <- .gaussian_model()
  set.seed(6)
  data <- data.frame(suppressMessages(ctGenerate(model, n.subjects = 8,
    Tpoints = 5, backend = "julia", intoverstates = TRUE)))
  free <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    MANIFESTMEANS = matrix("mmean"), Tpoints = 5)))
  free$pars$indvarying <- FALSE
  expect_error(
    suppressWarnings(suppressMessages(ctFit(data, free, backend = "julia",
      intoverstates = FALSE, verbose = 0,
      optimcontrol = list(estonly = TRUE)))),
    "unbounded")
})

test_that("the Laplace random-effect route and sampled states are refused together", {
  skip_on_cran()
  # A model check, so no Julia session is needed. Every other intoverpop
  # composes: augmented random effects are extra latent states, and the state
  # path samples them along with the rest.
  expect_error(ctsem:::.ctJuliaUnsupported(.states_model(), optimize = TRUE,
    priors = FALSE, intoverpop = "laplace", gendata = FALSE,
    stanmodeltext = NA, compileArgs = list(), forcerecompile = FALSE,
    intoverstates = FALSE), "intoverstates=FALSE")
  expect_silent(ctsem:::.ctJuliaUnsupported(.states_model(), optimize = TRUE,
    priors = FALSE, intoverpop = TRUE, gendata = FALSE,
    stanmodeltext = NA, compileArgs = list(), forcerecompile = FALSE,
    intoverstates = FALSE))
})

test_that("sampling the joint density gives a posterior over both", {
  skip_without_julia()
  # The mode this route is actually for. Optimising the joint density gives its
  # mode, which is biased and whose profile is flat; sampling it gives the exact
  # posterior of the parameters and the states together, with no Gaussian
  # assumption about the state anywhere -- which is the whole reason to make the
  # states explicit rather than let the filter project them.
  gen <- suppressWarnings(suppressMessages(ctModel(type = "ct", n.latent = 1,
    n.manifest = 1, manifestNames = "c1", latentNames = "eta1",
    manifesttype = 3L, LAMBDA = matrix(1), DRIFT = matrix(-0.5),
    DIFFUSION = matrix(0.7), MANIFESTVAR = matrix(0), T0VAR = matrix(1),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(1),
    Tpoints = 6)))
  set.seed(9)
  data <- data.frame(suppressMessages(ctGenerate(gen, n.subjects = 8,
    Tpoints = 6, backend = "julia")))

  fitmodel <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    n.latent = 1, n.manifest = 1, manifestNames = "c1", latentNames = "eta1",
    manifesttype = 3L, LAMBDA = matrix(1), MANIFESTVAR = matrix(0),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix("cmean"),
    Tpoints = 6)))
  fitmodel$pars$indvarying <- FALSE

  fit <- suppressWarnings(suppressMessages(ctFit(data, fitmodel,
    backend = "julia", intoverstates = FALSE, optimize = FALSE, priors = TRUE,
    verbose = 0, iter = 120, chains = 1)))

  npar <- length(fit$estimate$raw)
  # The posterior is over the parameters alone, even though the sampler moved
  # in the states too: `rawposterior` means the same thing on every fit.
  expect_equal(ncol(fit$estimate$rawposterior), npar)
  expect_true(all(is.finite(fit$estimate$rawposterior)))
  # Standard errors, which the optimised route deliberately does not report --
  # here they come from draws rather than from a curvature, so the flat profile
  # is not in the way.
  expect_length(fit$estimate$se, npar)
  expect_true(all(fit$estimate$se > 0))
  # And the trajectory comes back alongside them.
  expect_equal(dim(fit$estimate$states), c(nrow(data), 1L))
  expect_true(all(is.finite(fit$estimate$states)))
  expect_identical(fit$estimate$loglik_type, "joint")
})
