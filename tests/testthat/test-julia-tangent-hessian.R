# The hand-coded Hessian (`tangent_hessian.jl`) against the dual route, its
# oracle (`ctsem_set_tangent_hessian!(false)`), on objectives only R builds:
# random intercepts the subject filter integrates -- a block of correlated
# effects, and own variances -- a random TD effect integrated beside them, TI
# predictors, which give each subject its own parameter layer, and a TD
# predictor. The engine's test_tangent_hessian.jl covers the sweep itself. The
# two routes compute the same derivative, so they agree to rounding.

.th_data <- function() {
  set.seed(5)
  n <- 30; t <- 8
  d <- data.frame(id = rep(seq_len(n), each = t))
  d$time <- ave(stats::rexp(n * t, 1), d$id, FUN = cumsum)
  d$ti1 <- rep(stats::rnorm(n), each = t)
  d$td1 <- stats::rbinom(n * t, 1, 0.3)
  u <- rep(stats::rnorm(n, 0, 0.6), each = t)
  for (v in c("y1", "y2", "y3")) {
    d[[v]] <- u + 0.3 * d$ti1 + stats::rnorm(n * t)
    d[stats::runif(n * t) < 0.15, v] <- NA
  }
  d
}

.th_model <- function() {
  m <- suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 3,
    manifestNames = c("y1", "y2", "y3"), latentNames = c("eta1", "eta2"),
    TIpredNames = "ti1", TDpredNames = "td1",
    LAMBDA = matrix(c(1, "lam21", 0, 0, 0, 1), 3, 2), silent = TRUE))
  m$pars$indvarying[m$pars$matrix == "TDPREDEFFECT" & m$pars$row == 1] <- TRUE
  m
}

# What the filter integrates, whether the hand route serves the objective and
# takes the point, and the largest difference from the dual's Hessian relative
# to its largest entry.
.th_compare <- function(fit) {
  f <- ctsem:::.ctJuliaEval(paste(
    "function (lap, x)",
    "  M = ContinuousTimeSEM",
    "  o = lap.objective",
    "  xe = M._laplace_extend(lap, collect(Float64, x))",
    "  h = o.params.hybrid",
    "  why = M._tan_refusal(o, M._HYBRID_FULL)",
    "  H = zeros(length(xe), length(xe))",
    "  taken = why === nothing && M._tan_hessian_columns!(H, o, xe,",
    "    collect(1:length(xe)), nothing, M._HYBRID_FULL)",
    "  M.ctsem_set_tangent_hessian!(false)",
    "  H0 = try M.ctsem_hessian(lap, x) finally M.ctsem_set_tangent_hessian!(true) end",
    "  H1 = M.ctsem_hessian(lap, x)",
    "  (block = M._hybrid_on(h) && h.block, items = M._hybrid_items(h),",
    "    kinds = sort(unique(h.kind)), served = why === nothing, taken = taken,",
    "    err = maximum(abs.(H1 .- H0)) / maximum(abs.(H0)))",
    "end", sep = "\n"))
  f(ctsem:::.ctJuliaObjective(fit, hybrid = TRUE),
    ctsem:::.ctJuliaNumericVector(as.numeric(fit$estimate$raw)))
}

test_that("the hand-coded Hessian is the dual's where the filter integrates deviations", {
  skip_without_julia()
  d <- .th_data()
  for (own in c(FALSE, TRUE)) {
    args <- list(d, .th_model(), backend = "julia", intoverpop = "laplace",
      cores = 1, optimcontrol = list(estonly = TRUE))
    if (own) args <- c(args, list(poprank = 0, popunique = TRUE))
    fit <- suppressWarnings(suppressMessages(do.call(ctFit, args)))
    r <- .th_compare(fit)
    info <- if (own) "own variances" else "a block"
    expect_equal(r$block, !own, info = info)
    # Initial means, manifest intercepts and the TD effect, all integrated.
    expect_equal(as.integer(r$kinds), c(1L, 2L, 4L), info = info)
    expect_true(r$served, info = info)
    expect_true(r$taken, info = info)
    expect_lt(r$err, 1e-10)
  }
})
