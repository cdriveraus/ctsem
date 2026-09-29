# The categorical measurement update's moment Jacobian is written out by hand
# (`_quadrature_partials!`, binary_measurement.jl), and every Hessian and
# Laplace pass that runs through the reverse pass differentiates it again. A
# wrong partial there gives a plausible gradient and a plausible Hessian, so
# both orders are checked here against forward mode, which never reads it:
# the rule itself differentiated by ForwardDiff, and a two-state model's
# gradient and Hessian by forward-over-forward.
#
# The partials of the rule's placement -- where the mode sits and how wide
# the nodes are spread -- cancel to the rule's own error, since an exact
# integral does not depend on where it is centred. So an error there is
# invisible where the rule is accurate, and the diffuse cases below, where it
# is least accurate, are the ones that see it.

test_that("the rule's partials are forward mode's, to first and second order", {
  skip_without_julia()
  ctJuliaSetup()
  # Columns: the Jacobian and forward mode's, then its partials along a
  # direction -- as the Laplace seeds and the continuation Hessian take them --
  # and the moments' Hessians along the same direction.
  ctsem:::.ctJuliaEval('
    function _ctsem_test_catjac(a, b, y, tau, kind, direction)
      M = ContinuousTimeSEM
      FD = M.ForwardDiff
      a = Float64(a)
      b = Float64(b)
      nodes, weights = M._binary_rule()
      tau = collect(Float64, tau)
      f = x -> collect(M._binary_moments(x[1], sqrt(x[2]), y, nodes, weights,
        x[3:end], kind))
      x0 = [a; b; tau]
      _, _, _, J = M._binary_moment_jacobian(a, b, y, nodes, weights, tau, kind)
      D = FD.Dual{M._LaplaceSeedInner,Float64,1}
      d = direction[1:length(x0)]
      seed = (v, s) -> D(v, FD.Partials((s,)))
      _, _, _, JD = M._binary_moment_jacobian(seed(a, d[1]), seed(b, d[2]), y,
        nodes, weights, D[seed(tau[i], d[2 + i]) for i in eachindex(tau)], kind)
      second = [FD.partials(JD[r, j])[1] for r in 1:3, j in 1:length(x0)]
      reference = vcat([transpose(FD.hessian(x -> f(x)[r], x0) * d)
        for r in 1:3]...)
      return hcat(vec(J), vec(FD.jacobian(f, x0)), vec(second), vec(reference))
    end')
  kinds <- c(binary = 1L, ordinal = 2L, count = 3L)
  cases <- list(
    list("ordinal, interior", 0.7, 1.3, 3, c(-1, 0.4, 1.9), "ordinal"),
    list("ordinal, lowest", -0.4, 2.5, 1, c(-1, 0.4, 1.9), "ordinal"),
    list("ordinal, highest", 1.2, 0.4, 4, c(-1, 0.4, 1.9), "ordinal"),
    list("ordinal, narrow category", 0.2, 3, 2, c(-0.3, 0.2, 1.4), "ordinal"),
    list("binary, one", -1.1, 0.8, 1, numeric(0), "binary"),
    list("binary, zero", 0.9, 4, 0, numeric(0), "binary"),
    list("binary, diffuse", 0.4, 30, 1, numeric(0), "binary"),
    list("ordinal, diffuse", -0.5, 40, 2, c(-1, 0.4, 1.9), "ordinal"),
    list("count", 0.5, 0.6, 4, numeric(0), "count"))
  num <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
  for (cs in cases) {
    # A call written out in full rather than a juliaCall: evaluated as a new
    # top-level expression it sees the function just defined, and an empty
    # vector need not cross the bridge.
    got <- ctsem:::.ctJuliaEval(sprintf(
      "_ctsem_test_catjac(%s, %s, %d, Float64[%s], %d, [0.3, -0.2, 0.1, 0.25, -0.4])",
      num(cs[[2]]), num(cs[[3]]), as.integer(cs[[4]]), num(cs[[5]]),
      kinds[[cs[[6]]]]))
    expect_equal(got[, 1], got[, 2], tolerance = 1e-9, info = cs[[1]])
    expect_equal(got[, 3], got[, 4], tolerance = 1e-8, info = cs[[1]])
  }
})

# Two latent states with a cross-loaded ordinal indicator, so that `c = P λ`
# has two non-zero entries: with one state the covariance cotangent is a
# scalar, and a symmetry error in the reverse pass cannot show.
.jcatjac_data <- function(nsubjects = 6, nobs = 10, seed = 7) {
  set.seed(seed)
  A <- matrix(c(0.75, 0.1, -0.15, 0.8), 2, 2)
  rows <- lapply(seq_len(nsubjects), function(i) {
    eta <- matrix(0, nobs, 2)
    eta[1, ] <- stats::rnorm(2)
    for (t in 2:nobs) eta[t, ] <- A %*% eta[t - 1, ] + stats::rnorm(2, 0, 0.6)
    cut3 <- function(z) as.integer(cut(z + stats::rlogis(nobs),
      c(-Inf, -1, 0.3, 1.5, Inf)))
    data.frame(id = i, time = seq_len(nobs) - 1,
      o1 = cut3(eta[, 1] + 0.6 * eta[, 2]), o2 = cut3(eta[, 2]),
      b1 = stats::rbinom(nobs, 1, stats::plogis(eta[, 1])))
  })
  do.call(rbind, rows)
}

test_that("a two-state categorical model's gradient and Hessian are forward mode's", {
  skip_without_julia()
  m <- suppressMessages(ctModel(type = "ct", n.latent = 2, n.manifest = 3,
    manifestNames = c("o1", "o2", "b1"), latentNames = c("eta1", "eta2"),
    LAMBDA = matrix(c(1, "l12", 0, 1, 1, 0), 3, 2, byrow = TRUE),
    MANIFESTMEANS = matrix(0, 3, 1), CINT = matrix(0, 2, 1),
    T0MEANS = matrix(0, 2, 1), T0VAR = diag(1, 2),
    DIFFUSION = matrix(c(0.7, 0, 0.2, 0.6), 2, 2, byrow = TRUE),
    MANIFESTVAR = diag(0, 3), manifesttype = c(2L, 2L, 1L),
    ncategories = c(4L, 4L, 0L), silent = TRUE))
  m$pars$indvarying <- FALSE
  spec <- structure(ctsem:::.ctJuliaPrepare(.jcatjac_data(), m, priors = FALSE,
    intoverpop = "augmented"), class = c("ctJuliaModel", "ctFitModel"))
  npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
  set.seed(4)
  at <- stats::rnorm(npar, 0, 0.3)
  adjoint <- as.numeric(ctJuliaEvaluate(spec, at, gradient = TRUE,
    gradient_method = "adjoint")$gradient)
  forward <- as.numeric(ctJuliaEvaluate(spec, at, gradient = TRUE,
    gradient_method = "forward")$gradient)
  expect_equal(adjoint, forward, tolerance = 1e-9)
  # Forward over reverse differentiates the hand-written partials once more;
  # forward over forward never meets them.
  objective <- ctsem:::.ctJuliaObjective(spec)
  over_reverse <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_hessian",
    objective, ctsem:::.ctJuliaNumericVector(at))
  over_forward <- ctsem:::.ctJuliaCall("ContinuousTimeSEM.ctsem_hessian_forward",
    objective, ctsem:::.ctJuliaNumericVector(at))
  expect_equal(as.numeric(over_reverse), as.numeric(over_forward),
    tolerance = 1e-8)
})
