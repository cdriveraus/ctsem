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
  m$pars$indvarying <- ifelse(m$pars$matrix == "MANIFESTMEANS", TRUE,
    ifelse(m$pars$matrix == "T0MEANS", "subject", FALSE))
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
  # A session is only as wide as the call that started it, and this needs four
  # threads: at two, units split two ways and the members inside a unit never
  # run in parallel, so the nested axis would go untested and still pass. So
  # the fit restarts a narrower session rather than running at what it finds.
  fit <- withr::with_options(list(ctsem.julia.restart = TRUE),
    suppressWarnings(suppressMessages(ctFit(.band_data(), .band_model(),
      backend = "julia", intoverpop = "laplace", cores = 4, verbose = 0,
      optimcontrol = list(estonly = TRUE, maxiter = 3)))))
  expect_gte(as.integer(ctsem:::.ctJuliaEval("Threads.nthreads()")), 4L)

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

# A member under nested blocks takes its trace terms over its whole chain of
# blocks in one factor (`_laplace_member_chain_trace!`): here a person's four
# intercepts and their study's four drift effects, one jet a person, where the
# block-by-block assembly took four.
test_that("a nested unit's trace terms take each person's chain at once", {
  skip_without_julia()
  withr::local_options(ctsem.laplace.hybrid = FALSE)
  spec <- suppressMessages(ctFit(.nested_panel_data(2, 6, 8),
    .nested_panel_model("intercepts"), backend = "julia", intoverpop = "laplace",
    fit = FALSE, cores = 1))
  obj <- ctsem:::.ctJuliaObjective(spec, layer = TRUE)
  x <- rep(0, ctsem:::.ctBackendNpar(spec))
  for (lv in spec$laplace$levels) x[lv$sd_index[lv$sd_index > 0]] <- -1
  r <- ctsem:::.ctJuliaGet(ctsem:::.ctJuliaEval(paste(
    "(o, x) -> ContinuousTimeSEM._ctsem_with_max_chunks(1) do",
    "  C = ContinuousTimeSEM",
    "  C.ctsem_laplace_evaluate(o, x; gradient = false)",
    "  C.ctsem_reset_opcounts!()",
    "  r = C.ctsem_laplace_evaluate(o, x; gradient = true)",
    "  jets = C.ctsem_opcounts().sweep_jet",
    "  f = C._laplace_subject_gradients_forward(o, x)",
    "  theta = C._laplace_extend(o, collect(Float64, x))",
    "  Ls = C._laplace_popchols(theta, o.spec)",
    "  lam = map(eachindex(o.units.members)) do U",
    "    M = C._laplace_copy_blockmatrix(C._laplace_unit_curvature(o, U, theta, Ls, o.modes[U]))",
    "    blocks = o.units.blocks[U]",
    "    C._soft_smallest_eigenpair(M, blocks, C._soft_shifted_factor(M, blocks, 0.0)).lambda",
    "  end",
    "  (g = r.gradient, ref = vec(sum(f.scores; dims = 1)), jets = jets,",
    "   members = sum(length, o.units.members), units = length(o.units.members),",
    "   hard = all(>=(o.soft_blend_hi), lam))",
    "end", sep = "\n"))(obj, ctsem:::.ctJuliaNumericVector(x)))
  expect_equal(r$units, 2)
  # Every unit above the soft-direction quadrature's band, so the plain
  # assembly runs.
  expect_true(isTRUE(r$hard))
  expect_equal(as.numeric(r$g), as.numeric(r$ref), tolerance = 1e-8)
  expect_equal(r$jets, r$members)
})

# A nested unit's leaves integrate their own soft directions
# (`_laplace_block_soft_term`): each person's soft rule given the study at its
# mode, the whole unit's rule blended in as the study's own curvature, its
# Schur complement, turns soft. These Gaussian persons are never soft below
# the band's usual top (their curvature is the prior's and more), so the top is
# raised to 3 here, which puts every person's rule and the study's blend to
# work; the gradient is checked against forward mode over the term.
test_that("a nested unit's leaves integrate their own soft directions", {
  skip_without_julia()
  withr::local_options(ctsem.laplace.hybrid = FALSE)
  spec <- suppressMessages(ctFit(.nested_panel_data(3, 8, 4),
    .nested_panel_model("drift"), backend = "julia", intoverpop = "laplace",
    fit = FALSE, cores = 1))
  obj <- ctsem:::.ctJuliaObjective(spec, layer = TRUE)
  x <- rep(0, ctsem:::.ctBackendNpar(spec))
  for (lv in spec$laplace$levels) x[lv$sd_index[lv$sd_index > 0]] <- 1.5
  set.seed(1)
  x <- x + rnorm(length(x), 0, 0.5)
  r <- ctsem:::.ctJuliaGet(ctsem:::.ctJuliaEval(paste(
    "(o, x) -> begin",
    "  C = ContinuousTimeSEM",
    "  top = o.soft_blend_hi",
    "  o.soft_blend_hi = 3.0",
    "  try",
    "    theta = C._laplace_extend(o, collect(Float64, x))",
    "    g = C.ctsem_laplace_evaluate(o, x; gradient = true)",
    "    Ls = C._laplace_popchols(theta, o.spec)",
    "    soft = 0; w = Float64[]",
    "    for U in eachindex(o.units.members)",
    "      M = C._laplace_copy_blockmatrix(C._laplace_unit_curvature(o, U, theta, Ls, o.modes[U]))",
    "      blocks = o.units.blocks[U]; leaves = C._laplace_leaves(blocks)",
    "      soft += count(b -> C._laplace_leaf_soft_gate(o, M, b, o.soft_blend_hi), leaves)",
    "      push!(w, C._laplace_smoothstep_weight(C._laplace_rest_lowest(M, blocks, leaves).lambda,",
    "        o.soft_blend_lo, o.soft_blend_hi))",
    "    end",
    "    nested = C.ctsem_laplace_evaluate(o, x; gradient = true, nested_gradient = true)",
    "    was = C.ctsem_set_block_soft!(false)",
    "    whole = try C.ctsem_laplace_evaluate(o, x; gradient = false).value finally",
    "      C.ctsem_set_block_soft!(was) end",
    "    (g = g.gradient, value = g.value, nested = nested.gradient, soft = soft, w = w,",
    "     whole = whole)",
    "  finally",
    "    o.soft_blend_hi = top",
    "  end",
    "end", sep = "\n"))(obj, ctsem:::.ctJuliaNumericVector(x)))
  # Exercised: persons soft given the study, and the study's blend between
  # its ends.
  expect_gt(r$soft, 0)
  expect_true(any(r$w > 0 & r$w < 1))
  expect_equal(as.numeric(r$g), as.numeric(r$nested), tolerance = 1e-7)
  # The whole unit's rule is the oracle switch's, and another integral.
  expect_gt(abs(r$value - r$whole), 1e-3)
})
