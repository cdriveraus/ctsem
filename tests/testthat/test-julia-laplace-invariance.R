# Two classes of defect the Laplace route could have and none of the existing
# suite would catch, both borrowed from a bug found in the bigIRT package.
#
# 1. `rowsum(..., reorder=FALSE)` there grouped observations by the order they
#    were *encountered* rather than by their id's value, and it was invisible
#    whenever the data happened to arrive already sorted -- which is most of
#    the time, because most data is. The ctsem analogue would be building the
#    Laplace hierarchy (which subjects share a study, which raw index is
#    whose) from encounter order rather than from the id itself. `.ctJuliaPrepare`
#    sorts the data by (subject id, time) before anything else reads it
#    (R/ctJuliaBackend.R, `dat <- dat[order(dat[[model$subjectIDname]],
#    dat[[model$timeName]]), ]`), which should make every downstream
#    computation -- `subject_starts`, `.ctJuliaHierarchy`'s group vectors --
#    depend on id value rather than row position. The tests below assert that
#    invariant rather than assume it: shuffling which subject's block comes
#    first in the raw data must not move the objective, its gradient, or any
#    per-subject quantity once everything is remapped by id.
#
#    The nested case goes further than reshuffling row blocks (which the sort
#    above would neutralise even if the bug existed): it gives the *same*
#    eight subjects two different arrangements of id numbers -- one where
#    ascending id order groups each study contiguously, one where it
#    interleaves them -- so that the study-group index that
#    `.ctJuliaHierarchy` assigns by first encounter (`unique(labels)`,
#    R/ctJuliaBackend.R around line 1798) is a *different* integer for the
#    same real study in the two arrangements. If group membership were ever
#    read back inconsistently with how it was assigned, this is where it would
#    show: the model is otherwise identical, so the two must agree exactly.
#
# 2. A hand-assembled analytic gradient can have a term that is only wrong for
#    a correlation or off-diagonal loading away from zero, and the standard
#    diagnostic gap-at-zero hides it: several terms of a covariance gradient
#    are literally zero when the correlation is (see CLAUDE.md's account of
#    the same class of bug in bigIRT, where two of four gradient terms
#    vanished at rho=0). `.ctJuliaLaplaceSpec` lists the outer-gradient checks
#    it inherits from the engine's own suite (`test_laplace.jl` has a fixture
#    whose baseline correlation raw value is -0.15 and a test that perturbs it
#    further from zero -- see "the seeded gradient is exercised away from the
#    optimum"), but that suite has no fixture at all for a reduced-rank
#    `poprank` level, and the one R-side gradient test that does cover
#    reduced-rank levels (`test-poprank-levels.R`, "a reduced level's gradient
#    matches finite differences") evaluates at `fit$estimate$raw` from a
#    3-iteration fit started at `rnorm(npar, 0, .01)` -- so its loadings are
#    still close to their near-zero starting point when the gradient is
#    checked. The second half of this file adds a check at loadings and a
#    correlation coordinate set explicitly away from zero.

# --- helpers shared by both parts --------------------------------------------

# Each subject's content is generated from a seed keyed to its *role* (its
# position in a fixed underlying population), not to whatever numeric id or
# row position it is eventually given. That is what makes "the same subjects,
# relabelled" a meaningful thing to build: two arrangements can be given
# different id numbers and still be guaranteed to carry identical content.
.invar_role_series <- function(role, intercept, nobs = 4L, drift = -0.4,
    diffusion = 0.5) {
  set.seed(9000L + role)
  state <- stats::rnorm(1, 0, 0.4)
  out <- numeric(nobs)
  for (t in seq_len(nobs)) {
    if (t > 1L) {
      decay <- exp(drift)
      state <- decay * state + stats::rnorm(1, 0,
        sqrt(diffusion^2 / (-2 * drift) * (1 - decay^2)))
    }
    out[t] <- state + intercept + stats::rnorm(1, 0, 0.25)
  }
  data.frame(time = seq_len(nobs) - 1L, Y1 = out)
}

# The engine's own primitive behind `.ctBackendLaplaceSubjectPars`
# (R/ctBackendSummary.R): each subject's raw parameter vector after applying
# its own conditional-mode deviation, solved fresh from the given point. No
# fit object needed -- just the prepared spec and a raw vector -- which is
# what lets this run at a chosen point rather than at wherever an optimiser
# stopped.
.invar_subject_raw <- function(spec, raw) {
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(spec)
  out <- .ctBackendJuliaValue(module$ctsem_laplace_subject_values(objective,
    .ctJuliaNumericVector(raw)))
  matrix(as.numeric(out), ncol = length(raw))
}

# --- part 1: subject-order invariance ----------------------------------------

.invar_flat_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model
}

.invar_flat_data <- function(nsubjects = 10L, nobs = 4L) {
  do.call(rbind, lapply(seq_len(nsubjects), function(role) {
    set.seed(8000L + role)
    intercept <- stats::rnorm(1, 1.2, 0.8)
    d <- .invar_role_series(role, intercept, nobs)
    d$id <- role
    d
  }))
}

test_that("the laplace value, gradient and per-subject quantities do not depend on row order", {
  skip_without_julia()
  model <- .invar_flat_model()
  dat <- .invar_flat_data()

  # Shuffle which subject's block of rows comes first; ids and each subject's
  # own rows are untouched.
  set.seed(42)
  shuffled_ids <- sample(unique(dat$id))
  shuffled <- do.call(rbind, lapply(shuffled_ids, function(i) dat[dat$id == i, ]))
  rownames(shuffled) <- NULL
  # Confirm the shuffle actually moved something, or the test proves nothing.
  expect_false(identical(dat$id, shuffled$id))

  spec_a <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  spec_b <- suppressMessages(ctFit(shuffled, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  expect_equal(spec_b$laplace$npar, spec_a$laplace$npar)

  # A fixed, non-trivial point -- not the origin, where a term that only
  # depends on grouping could happen to be zero regardless of whether the
  # grouping is right.
  set.seed(11)
  raw <- stats::rnorm(spec_a$laplace$npar, 0, 0.3)

  ev_a <- ctJuliaEvaluate(spec_a, raw, gradient = TRUE, contributions = TRUE)
  ev_b <- ctJuliaEvaluate(spec_b, raw, gradient = TRUE, contributions = TRUE)

  expect_equal(as.numeric(ev_b$value), as.numeric(ev_a$value), tolerance = 1e-10)
  expect_equal(as.numeric(ev_b$gradient), as.numeric(ev_a$gradient),
    tolerance = 1e-8)

  # Per-subject contributions, mapped by id -- `labels` is the sorted-unique-id
  # vector `subject_loglik` is indexed by (R/ctJuliaBackend.R,
  # `.ctJuliaHierarchy`'s level-one `labels = subjects`).
  ids_a <- spec_a$laplace$levels[[1]]$labels
  ids_b <- spec_b$laplace$levels[[1]]$labels
  b_at_a <- match(ids_a, ids_b)
  expect_false(anyNA(b_at_a))
  expect_equal(as.numeric(ev_b$subject_loglik)[b_at_a],
    as.numeric(ev_a$subject_loglik), tolerance = 1e-10)

  # And the conditional-mode-driven per-subject parameter values: column
  # `re_index[1]` is "mmean"'s raw position, and each row is that subject's
  # own realised value of it.
  col <- spec_a$laplace$levels[[1]]$re_index[1]
  subj_a <- .invar_subject_raw(spec_a, raw)
  subj_b <- .invar_subject_raw(spec_b, raw)
  expect_equal(subj_b[b_at_a, col], subj_a[, col], tolerance = 1e-8)
})

# --- part 1, continued: a nested hierarchy under two different id labellings -

# Eight subjects, four per study. `study_of_role` fixes which study each role
# truly belongs to; `.invar_nested_build` then assigns the *numbers* used as
# subject id, which is the only thing that differs between the two
# arrangements below.
.invar_study_of_role <- c(rep("A", 4L), rep("B", 4L))

.invar_nested_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
    T0VAR = matrix(0.5), MANIFESTMEANS = matrix("mmean"),
    id = c("subject", "study"))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% "mmean"] <- TRUE
  model$pars$indvarying_study <- FALSE
  model$pars$indvarying_study[model$pars$param %in% "mmean"] <- TRUE
  model
}

.invar_nested_build <- function(id_for_role) {
  roles <- seq_along(id_for_role)
  rows <- lapply(roles, function(role) {
    study <- .invar_study_of_role[role]
    set.seed(7000L + role)
    intercept <- (if (study == "A") 0.6 else -0.6) + stats::rnorm(1, 0, 0.3)
    d <- .invar_role_series(role, intercept)
    d$subject <- id_for_role[role]
    d$study <- study
    d
  })
  out <- do.call(rbind, rows)
  out[order(out$subject, out$time), ]
}

test_that("a nested hierarchy agrees under two different id labellings of the same subjects", {
  skip_without_julia()
  model <- .invar_nested_model()

  # Contiguous: ascending id order visits study A's four subjects, then B's.
  contiguous_ids <- 1:8
  # Interleaved: ascending id order alternates B, A, B, A, ... To find each
  # role's id: role r's study is A for r in 1:4 and B for r in 5:8, and the id
  # sequence below places a B-role at every odd id and an A-role at every even
  # one, so the two arrangements assign the group index ("which integer means
  # study A") in the opposite order that `unique()` first encounters it.
  interleaved_ids <- c(2, 4, 6, 8, 1, 3, 5, 7)

  contiguous <- .invar_nested_build(contiguous_ids)
  interleaved <- .invar_nested_build(interleaved_ids)
  # Confirm the two arrangements really do put the studies in a different
  # order along ascending id -- otherwise this is the flat test again.
  expect_false(identical(
    unique(contiguous$study[order(contiguous$subject)]),
    unique(interleaved$study[order(interleaved$subject)])))

  spec_c <- suppressMessages(ctFit(contiguous, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  spec_i <- suppressMessages(ctFit(interleaved, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  expect_equal(spec_i$laplace$npar, spec_c$laplace$npar)
  expect_equal(vapply(spec_i$laplace$levels, `[[`, integer(1), "ngroups"),
    vapply(spec_c$laplace$levels, `[[`, integer(1), "ngroups"))

  set.seed(13)
  raw <- stats::rnorm(spec_c$laplace$npar, 0, 0.3)

  ev_c <- ctJuliaEvaluate(spec_c, raw, gradient = TRUE, contributions = TRUE)
  ev_i <- ctJuliaEvaluate(spec_i, raw, gradient = TRUE, contributions = TRUE)

  # The two describe the same eight subjects and studies under different
  # numbering, so the marginal likelihood and its gradient must agree exactly,
  # not merely be close.
  expect_equal(as.numeric(ev_i$value), as.numeric(ev_c$value), tolerance = 1e-10)
  expect_equal(as.numeric(ev_i$gradient), as.numeric(ev_c$gradient),
    tolerance = 1e-8)

  # Per-subject contributions, mapped by role (each arrangement's `id_for_role`
  # inverted through its own sorted-id `labels`, then matched role to role).
  role_of <- function(spec, id_for_role) match(spec$laplace$levels[[1]]$labels,
    id_for_role)
  role_c <- role_of(spec_c, contiguous_ids)
  role_i <- role_of(spec_i, interleaved_ids)
  expect_false(anyNA(role_c)); expect_false(anyNA(role_i))
  i_at_c <- match(role_c, role_i)
  expect_false(anyNA(i_at_c))
  expect_equal(as.numeric(ev_i$subject_loglik)[i_at_c],
    as.numeric(ev_c$subject_loglik), tolerance = 1e-10)
})

# --- part 2: gradients away from the zero-correlation special case ----------

.invar_corr_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    T0VAR = diag(0.5, 2),
    MANIFESTMEANS = matrix(c("mmean1", "mmean2"), 2, 1))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% c("mmean1", "mmean2")] <- TRUE
  model
}

.invar_corr_data <- function(nsubjects = 10L, nobs = 4L) {
  do.call(rbind, lapply(seq_len(nsubjects), function(role) {
    set.seed(6000L + role)
    i1 <- stats::rnorm(1, 1.0, 0.6); i2 <- stats::rnorm(1, -0.5, 0.6)
    d1 <- .invar_role_series(role, i1, nobs)
    d2 <- .invar_role_series(role + 1000L, i2, nobs)
    data.frame(id = role, time = d1$time, Y1 = d1$Y1, Y2 = d2$Y1)
  }))
}

# Central difference of the objective's *value*, matching the standard the
# suite already uses for the seeded Laplace gradient (test-poprank-levels.R,
# "a reduced level's gradient matches finite differences"): each perturbation
# re-solves the profiled inner modes, so a smaller step sharpens nothing and
# 1e-4 is about as good as this comparison gets.
.invar_central_difference <- function(spec, at, step = 1e-4) {
  vapply(seq_along(at), function(j) {
    h <- step * max(1, abs(at[[j]]))
    up <- at; up[j] <- up[j] + h
    lo <- at; lo[j] <- lo[j] - h
    up_v <- as.numeric(ctJuliaEvaluate(spec, up, gradient = FALSE)$value)
    lo_v <- as.numeric(ctJuliaEvaluate(spec, lo, gradient = FALSE)$value)
    (up_v - lo_v) / (2 * h)
  }, numeric(1))
}

test_that("the outer gradient matches finite differences with the correlation clearly away from zero", {
  skip_without_julia()
  model <- .invar_corr_model()
  dat <- .invar_corr_data()

  spec <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  laplace <- spec$laplace
  expect_equal(laplace$nrandom, 2L)
  expect_length(laplace$cor_index, 1L)

  set.seed(21)
  raw <- stats::rnorm(laplace$npar, 0, 0.2)
  # Population scales and the correlation coordinate set explicitly away from
  # zero and from each other -- the point the bigIRT-style defect hides at.
  raw[laplace$sd_index] <- c(0.4, -0.3)
  raw[laplace$cor_index] <- 0.55

  analytic <- as.numeric(ctJuliaEvaluate(spec, raw, gradient = TRUE)$gradient)
  expect_true(all(is.finite(analytic)))
  numeric_grad <- .invar_central_difference(spec, raw)

  expect_equal(analytic, numeric_grad, tolerance = 1e-3)
  # The correlation and scale coordinates specifically, not just the vector as
  # a whole: a term that only enters through `psi` or through `L` could cancel
  # in a coarser comparison.
  expect_equal(analytic[laplace$cor_index], numeric_grad[laplace$cor_index],
    tolerance = 1e-3)
  expect_equal(analytic[laplace$sd_index], numeric_grad[laplace$sd_index],
    tolerance = 1e-3)
})

# --- part 2, continued: a reduced-rank level's off-diagonal loadings --------

.invar_poprank_model <- function() {
  model <- suppressWarnings(suppressMessages(ctModel(
    type = "ct", manifestNames = c("Y1", "Y2", "Y3"),
    latentNames = c("eta1", "eta2", "eta3"), LAMBDA = diag(3),
    T0MEANS = matrix(0, 3, 1), CINT = matrix(0, 3, 1), T0VAR = diag(0.5, 3),
    MANIFESTMEANS = matrix(paste0("mmean", 1:3), 3, 1))))
  model$pars$indvarying <- FALSE
  model$pars$indvarying[model$pars$param %in% paste0("mmean", 1:3)] <- TRUE
  model
}

.invar_poprank_data <- function(nsubjects = 10L, nobs = 4L) {
  do.call(rbind, lapply(seq_len(nsubjects), function(role) {
    set.seed(5000L + role)
    ii <- stats::rnorm(3, c(1.0, -0.5, 0.3), 0.5)
    d1 <- .invar_role_series(role, ii[1], nobs)
    d2 <- .invar_role_series(role + 1000L, ii[2], nobs)
    d3 <- .invar_role_series(role + 2000L, ii[3], nobs)
    data.frame(id = role, time = d1$time, Y1 = d1$Y1, Y2 = d2$Y1, Y3 = d3$Y1)
  }))
}

test_that("a reduced level's gradient matches finite differences with its loadings away from zero", {
  skip_without_julia()
  model <- .invar_poprank_model()
  dat <- .invar_poprank_data()

  spec <- suppressMessages(ctFit(dat, model, backend = "julia",
    intoverpop = "laplace", poprank = 2, fit = FALSE))
  level <- spec$laplace$levels[[1]]
  expect_equal(level$rank, 2L)
  # k=3, rank=2: 3*2 - 2*1/2 = 5 loadings, some of them off the natural
  # diagonal -- the coordinates `test-poprank-levels.R`'s gradient check never
  # moves far from their near-zero starting values.
  expect_length(level$load_index, 5L)

  set.seed(23)
  raw <- stats::rnorm(spec$laplace$npar, 0, 0.2)
  # Loadings with mixed signs and magnitude clearly away from zero.
  raw[level$load_index] <- c(0.6, -0.5, 0.4, -0.45, 0.55)

  analytic <- as.numeric(ctJuliaEvaluate(spec, raw, gradient = TRUE)$gradient)
  expect_true(all(is.finite(analytic)))
  numeric_grad <- .invar_central_difference(spec, raw)

  expect_equal(analytic, numeric_grad, tolerance = 1e-3)
  expect_equal(analytic[level$load_index], numeric_grad[level$load_index],
    tolerance = 1e-3)
})
