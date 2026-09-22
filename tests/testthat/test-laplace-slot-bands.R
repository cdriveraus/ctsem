# The parallel width must not change the answer, and must not leak between
# routes.
#
# The Laplace path divides work on two nested axes -- units, and the members
# inside one unit -- and every task needs an adjoint workspace of its own.
# Those come from `laplace.workspaces`, indexed by a slot. A task holds a
# contiguous *band* of slots and lends disjoint sub-bands to anything it
# spawns, so two live tasks can never name the same one.
#
# Two things went wrong before that band was ambient rather than computed, and
# this file exists for both.
#
# The width leaked between routes. Slots were derived from the global chunk
# count, which the fit tuner *pins* for the rest of the session, so a fit
# followed by a `ctLaplaceCheck()` asked for slots the quadrature route had
# never budgeted. Nothing in the suite caught it, because the quadrature tests
# fit on one core and so never set a width above one.
#
# And the width changed the answer. A missing band at one call site left a
# whole phase serial; a shared local in a closure had two tasks writing one
# variable. Neither errored -- the first was merely slow, the second surfaced
# once as a `BoundsError` inside a spawned task, which the objective wrapper
# turns into a large finite penalty, so the optimiser sees a bad point rather
# than a bug. The assertion that catches both is that the gradient is the same
# serially and in parallel: the pool partitions members, so each worker's
# columns are summed in the same order whatever the width.

.band_model <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), T0MEANS = matrix(paste0("t0_", 1:2), 2, 1),
    T0VAR = diag(2), MANIFESTMEANS = matrix(paste0("mm_", 1:2), 2, 1),
    MANIFESTVAR = "diag", DRIFT = "auto", CINT = matrix(0, 2, 1),
    DIFFUSION = "auto", id = c("subject", "grp"), time = "time",
    tipredDefault = FALSE)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying_grp <- FALSE
  m$pars$indvarying[m$pars$matrix %in% c("T0MEANS", "MANIFESTMEANS")] <- TRUE
  m$pars$indvarying_grp[m$pars$matrix == "MANIFESTMEANS"] <- TRUE
  m
}

.band_data <- function(ngroup = 3L, nperson = 12L) {
  set.seed(4)
  d <- expand.grid(time = c(0, .2, .45, .7), person = seq_len(nperson),
    grp = seq_len(ngroup))
  d <- d[order(d$grp, d$person, d$time), ]
  d$subject <- paste0(d$grp, "_", d$person)
  d$Y1 <- rnorm(nrow(d))
  d$Y2 <- rnorm(nrow(d))
  d
}

test_that("the parallel width does not change the objective or its gradient", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(ctFit(.band_data(), .band_model(),
    backend = "julia", intoverpop = "laplace", cores = 4, verbose = 0,
    optimcontrol = list(estonly = TRUE, maxiter = 3))))

  at <- as.numeric(fit$estimate$raw)
  lpg <- ctsem:::.ctBackendLpgFunc(fit, gradient = TRUE)
  set_width <- function(n) JuliaConnectoR::juliaCall(
    "ContinuousTimeSEM.ctsem_set_max_chunks!", as.integer(n))
  original <- as.integer(JuliaConnectoR::juliaEval(
    "ContinuousTimeSEM.ctsem_max_chunks().max_chunks"))
  on.exit(set_width(original), add = TRUE)

  set_width(1L)
  serial <- lpg(at)
  set_width(4L)
  wide <- lpg(at)

  # Tight, and not `tolerance = 1e-6`. A difference at this size is a slot
  # collision or a shared local, not rounding -- see the header. The members a
  # worker owns are summed in their own order, so the last bit can move where
  # the partition does, and no more than that.
  expect_equal(as.numeric(wide), as.numeric(serial), tolerance = 1e-12)
  expect_equal(as.numeric(attr(wide, "gradient")),
    as.numeric(attr(serial, "gradient")), tolerance = 1e-10)
  expect_true(all(is.finite(attr(serial, "gradient"))))
})

test_that("a fit's pool width leaves the other routes working", {
  skip_without_julia()
  fit <- suppressWarnings(suppressMessages(ctFit(.band_data(), .band_model(),
    backend = "julia", intoverpop = "laplace", cores = 4, verbose = 0,
    optimcontrol = list(estonly = TRUE, maxiter = 3))))
  expect_true(is.finite(as.numeric(fit$estimate$loglik)))

  # The quadrature route claims its own band out of the pool the fit sized,
  # rather than inheriting a count it never budgeted for. Before the pool this
  # either indexed past the workspace vector or took a slot another chunk was
  # using.
  chk <- suppressWarnings(suppressMessages(
    ctLaplaceCheck(fit, nodes = 3L, correction = FALSE)))
  expect_s3_class(chk, "ctLaplaceCheck")
  expect_true(all(is.finite(unlist(chk[vapply(chk, is.numeric, logical(1))]))))
})
