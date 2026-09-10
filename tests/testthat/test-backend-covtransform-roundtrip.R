# A reported covariance must be the construction of the reported raw matrix.
#
# This file exists because one defect shape has now appeared four times in this
# area, and every instance was invisible: a *reporting* path building a
# covariance by a route of its own, or with a construction hardcoded, while the
# likelihood used another. The estimates, the gradients and the fit were right
# each time; only what a user reads was wrong, and it looked entirely
# plausible.
#
#   - the stan `generated quantities` lookup that stopped matching, so every
#     population summary came out divided by its parameter's multiplier
#   - three hardcoded population-covariance sites that ignored
#     `covmattransform` while the filter honoured it
#   - `_ctsem_pack_matrices!` reporting a T0cov with no population block while
#     the filter's had one
#   - `_laplace_popchol` spelling out the code=0 construction by hand
#
# The check is a round trip rather than an independent reimplementation: for
# every covariance the summary reports, feed the *raw* matrix it reports
# alongside it back through the engine's own `sdcovsqrt2cov` at the model's own
# construction code, and require the two to agree to machine precision. That
# catches a reporting path taking a different route and a reporting path
# ignoring the code, which is all four incidents, and it does it for every
# construction the backend accepts rather than for whichever one the fixture
# happens to use. Whether `sdcovsqrt2cov` is itself correct is a separate
# question, answered against independent references in
# `test_constrain_cor_sqrt.jl` and `test_cov_expm.jl`.
#
# Three constructions, not four: the julia backend refuses 'rawcorr_indep',
# which selects the same construction as 'rawcorr' and differs only in the
# prior, so accepting it would accept a setting that does nothing. See the
# covmattransform gate in `.ctJuliaUnsupported`. Adding a code here is the
# natural place to notice that its reporting path was never checked.
#
# No random effects in the fixture, deliberately. The augmentation trims rows
# out of some reported arrays after they are built, so `construct(trim(raw))`
# and `trim(construct(raw))` are not the same object and a round trip over them
# would be comparing the wrong things. The population block's own round trip is
# the stan comparison in `test-backend-summary.R`.

.covrt_model <- function(transform) {
  model <- suppressWarnings(ctModel(type = "ct", n.latent = 3, n.manifest = 3,
    LAMBDA = diag(3), manifestNames = paste0("Y", 1:3),
    latentNames = paste0("eta", 1:3),
    # Full DIFFUSION and MANIFESTVAR: the off-diagonal coordinates are what the
    # construction differs on, so a diagonal fixture would pass under any of
    # the four codes and prove nothing.
    DRIFT = matrix(c("dr1", 0, 0, "cr21", "dr2", 0, 0, 0, "dr3"), 3, 3,
      byrow = TRUE),
    T0MEANS = matrix(0, 3, 1), MANIFESTMEANS = matrix(0, 3, 1),
    CINT = matrix(0, 3, 1)))
  model$covmattransform <- transform
  model
}

.covrt_data <- function() {
  set.seed(11)
  do.call(rbind, lapply(1:8, function(i) data.frame(id = i,
    time = c(0, .4, 1.1, 2.0),
    Y1 = stats::rnorm(4, 0, .5), Y2 = stats::rnorm(4, 0, .5),
    Y3 = stats::rnorm(4, 0, .5))))
}

.covrt_pointfit <- function(spec, model, raw) {
  structure(list(model_spec = spec, model = model, backend = "julia",
    estimate = list(raw = raw, loglik = NA_real_)),
    class = c("ctJuliaFit", "ctFit"))
}

test_that("every covmattransform reports covariances its own raw matrices rebuild", {
  skip_without_julia()
  data <- .covrt_data()
  # Which reported covariance comes from which reported raw matrix. Both pairs
  # are required below rather than merely attempted, so a rename that drops one
  # from the layout fails the test instead of quietly halving it.
  pairs <- list(c("DIFFUSION", "DIFFUSIONcov"), c("MANIFESTVAR", "MANIFESTcov"),
    c("T0VAR", "T0cov"))

  for (transform in c("rawcorr", "cholesky", "z")) {
    model <- .covrt_model(transform)
    spec <- suppressMessages(ctFit(data, model, backend = "julia", fit = FALSE))
    code <- if (is.null(spec$covmatcode)) 0L else as.integer(spec$covmatcode)
    npar <- max(spec$parameter_table$parnumber, na.rm = TRUE)
    set.seed(4)
    raw <- stats::rnorm(npar, 0, .4)
    fit <- .covrt_pointfit(spec, model, raw)
    pop <- ctsem:::.ctBackendPopArrays(fit)

    # The engine's own construction, reached directly. A closure rather than a
    # `juliaCall` per matrix: `sdcovsqrt2cov` returns a `Symmetric` wrapper the
    # bridge would hand back as a struct.
    construct <- JuliaConnectoR::juliaEval(
      "(m, c) -> Matrix(ContinuousTimeSEM.sdcovsqrt2cov(m, Int(c)))")

    for (pair in pairs) {
      rawname <- paste0("pop_", pair[1])
      covname <- paste0("pop_", pair[2])
      expect_true(rawname %in% names(pop), info = rawname)
      expect_true(covname %in% names(pop), info = covname)
      rawmat <- drop(pop[[rawname]])
      reported <- drop(pop[[covname]])
      rebuilt <- construct(rawmat, code)
      expect_equal(dim(reported), dim(rebuilt),
        info = paste(transform, covname))
      expect_equal(as.numeric(reported), as.numeric(rebuilt),
        tolerance = 1e-12, info = paste(transform, covname))
    }

    # And the two constructions actually differ, so the agreement above is not
    # one comparison twice over. 'z' is the code that would have exposed the
    # hardcoded sites, so it is the one asked to disagree with 0.
    if (identical(transform, "z")) {
      other <- construct(drop(pop$pop_DIFFUSION), 0L)
      expect_false(isTRUE(all.equal(as.numeric(drop(pop$pop_DIFFUSIONcov)),
        as.numeric(other), tolerance = 1e-6)))
    }
  }
})
