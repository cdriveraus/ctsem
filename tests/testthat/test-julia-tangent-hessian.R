# The hand-coded Hessian (`tangent_hessian.jl`) against the dual route, its
# oracle (the engine tests' dual_oracles.jl), on objectives only R builds:
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

# The engine tests' oracles, forward mode through the filter, into the session.
.th_oracles <- function() {
  path <- system.file("julia", "ContinuousTimeSEM", "test", "dual_oracles.jl",
    package = "ctsem")
  ctsem:::.ctJuliaEval(sprintf(
    'isdefined(Main, :_oracle_hessian) || include("%s")', path))
}

# What the filter integrates, whether the hand route serves the objective and
# takes the point, and the largest difference from the dual's Hessian relative
# to its largest entry.
.th_compare <- function(fit) {
  .th_oracles()
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
    "  H0 = _oracle_hessian(o, xe; hybrid_mode = M._HYBRID_FULL)",
    "  H1 = M.ctsem_hessian(lap, x)",
    "  (block = M._hybrid_on(h) && h.block, items = M._hybrid_items(h),",
    "    kinds = sort(unique(h.kind)), served = why === nothing, taken = taken,",
    "    err = maximum(abs.(H1 .- H0)) / maximum(abs.(H0)))",
    "end", sep = "\n"))
  f(ctsem:::.ctJuliaObjective(fit, hybrid = TRUE),
    ctsem:::.ctJuliaNumericVector(as.numeric(fit$estimate$raw)))
}

# The unit curvature by the hand-coded sweeps (`tangent_laplace.jl`) against the
# dual's, with the gated floor's cross curvature read off the same sweeps
# (`CTSEMHybridCross`), at each unit's mode: the largest differences in the
# curvature, and in the cross terms, mode and `M`, relative to their largest
# entries, the members with a cross term, and whether the hand route served.
.th_curvature <- function(fit) {
  .th_oracles()
  f <- ctsem:::.ctJuliaEval(paste(
    "function (lap, x)",
    "  M = ContinuousTimeSEM",
    "  xe = M._laplace_extend(lap, collect(Float64, x))",
    "  M.ctsem_laplace_evaluate(lap, xe; gradient = false)",
    "  Ls = M._laplace_popchols(xe, lap.spec)",
    "  errM = 0.0; errH = 0.0; errJ = 0.0; ncross = 0; units = 0",
    "  for U in eachindex(lap.units.members)",
    "    u = lap.modes[U]",
    "    isempty(u) && continue",
    "    units += 1",
    "    d, nm = length(u), length(lap.units.members[U])",
    "    c1 = M.CTSEMHybridCross{Float64}(d, nm)",
    "    A1 = M._laplace_block_dense(M._laplace_unit_curvature(lap, U, xe, Ls, u;",
    "      hook = c1), lap.units.blocks[U], d)",
    "    c0 = M.CTSEMHybridCross{Float64}(d, nm)",
    "    A0 = _oracle_unit_curvature(lap, U, xe, Ls, u; hook = c0)",
    "    errM = max(errM, maximum(abs.(A1 .- A0)) / maximum(abs.(A0)))",
    "    for m in 1:nm",
    "      c1.q[m] == c0.q[m] || return (errM = Inf, errH = Inf, errJ = Inf, ncross = -1,",
    "        units = units, served = false)",
    "      q = c0.q[m]",
    "      q == 0 && continue",
    "      ncross += 1",
    "      rel(a, b) = maximum(abs.(a .- b)) / max(1.0, maximum(abs.(b)))",
    "      errH = max(errH, rel(c1.H[m][1:q, :], c0.H[m][1:q, :]),",
    "        rel(c1.A[m][1:q, 1:q], c0.A[m][1:q, 1:q]), rel(c1.dhat[m][1:q], c0.dhat[m][1:q]))",
    "    end",
    "    # The joint curvature, the deviations fixed a little off their modes.",
    "    any(>(0), c0.q) || continue",
    "    layout = M._hybrid_layout(d, c0.q)",
    "    z = vcat(u, reduce(vcat, [c0.dhat[m][1:c0.q[m]] .+ 0.1 for m in 1:nm]))",
    "    J1 = M._hybrid_joint_curvature(lap, U, xe, Ls, z, layout).diag[1]",
    "    J0 = _oracle_joint_curvature(lap, U, xe, Ls, z, layout)",
    "    errJ = max(errJ, maximum(abs.(J1 .- J0)) / maximum(abs.(J0)))",
    "  end",
    "  served = M._tan_laplace_context(lap, xe, M._HYBRID_PROFILE) !== nothing &&",
    "    M._tan_laplace_context(lap, xe, M._HYBRID_FIXED) !== nothing",
    "  (errM = errM, errH = errH, errJ = errJ, ncross = ncross, units = units,",
    "    served = served)",
    "end", sep = "\n"))
  f(ctsem:::.ctJuliaObjective(fit, hybrid = TRUE),
    ctsem:::.ctJuliaNumericVector(as.numeric(fit$estimate$raw)))
}

test_that("the hand-coded unit curvature is the dual's, cross terms included", {
  skip_without_julia()
  d <- .th_data()
  m <- .th_model()
  # A random drift, a coordinate, beside the intercepts the filter integrates.
  m$pars$indvarying[m$pars$matrix == "DRIFT" & m$pars$row == 1 & m$pars$col == 1] <- TRUE
  for (own in c(FALSE, TRUE)) {
    args <- list(d, m, backend = "julia", intoverpop = "laplace", cores = 1,
      optimcontrol = list(estonly = TRUE))
    if (own) args <- c(args, list(poprank = 0, popunique = TRUE))
    fit <- suppressWarnings(suppressMessages(do.call(ctFit, args)))
    r <- .th_curvature(fit)
    info <- if (own) "own variances" else "a block"
    expect_true(r$served, info = info)
    expect_gt(r$units, 0)
    expect_gt(r$ncross, 0)
    expect_lt(r$errM, 1e-10)
    expect_lt(r$errH, 1e-10)
    expect_lt(r$errJ, 1e-10)
  }
})

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
