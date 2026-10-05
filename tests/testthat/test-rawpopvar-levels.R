# RAWPOPVAR at every level of the hierarchy, held on the Laplace route.
#
# `RAWPOPVAR` states the subject level's population spread and
# `RAWPOPVAR_<idname>` each grouping level's, as `indvarying_<idname>` and
# `sdscale_<idname>` already name a level's own columns. A number fixes that
# raw-scale sd or correlation coordinate. On the Laplace and 'none' routes a
# fixed entry takes no place in the raw vector: it reaches the engine as index
# 0 with its raw value, and `_laplace_popchol` builds the covariance from it.
# Before this, only the augmented route read a stated number, and a Laplace fit
# estimated the spread as though nothing had been said.

.rpv_model <- function() {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    id = c("id", "study"), manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), DRIFT = matrix(-0.5), DIFFUSION = matrix(.5),
    MANIFESTVAR = matrix(.3), MANIFESTMEANS = matrix("mm"),
    T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0))))
  m$pars$indvarying <- m$pars$param %in% "mm"
  m$pars$indvarying_study <- m$pars$param %in% "mm"
  m
}
.rpv_data <- function() data.frame(id = rep(1:8, each = 4),
  study = rep(1:2, each = 16), time = rep(0:3, 8), Y1 = rnorm(32))

test_that("each grouping level has its own RAWPOPVAR", {
  m <- .rpv_model()
  mats <- m$matrices
  expect_true(all(c("RAWPOPVAR", "RAWPOPVAR_study") %in% names(mats)))
  expect_equal(unname(mats$RAWPOPVAR_study["mm", "mm"]), "popsd_mm.study")
  mats$RAWPOPVAR_study["mm", "mm"] <- 0.3
  m$matrices <- mats
  expect_equal(unname(m$RAWPOPVAR_study["mm", "mm"]), "0.3")
  expect_equal(unname(m$RAWPOPVAR["mm", "mm"]), "popsd_mm")
  # A level nothing varies at has none.
  m$pars$indvarying_study <- FALSE
  expect_null(ctsem:::.ctModelRawPopVarSync(m)[["RAWPOPVAR_study"]])
})

test_that("a stated sd is held on the Laplace route, at its own level", {
  skip_without_julia()
  m <- .rpv_model()
  free <- suppressMessages(ctFit(.rpv_data(), m, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  m$RAWPOPVAR_study["mm", "mm"] <- 0.3
  spec <- suppressMessages(ctFit(.rpv_data(), m, backend = "julia",
    intoverpop = "laplace", fit = FALSE))
  levels <- spec$laplace$levels
  expect_gt(levels[[1]]$sd_index, 0L)
  expect_equal(levels[[2]]$sd_index, 0L)
  # One parameter fewer: the study sd is not estimated.
  expect_equal(ctsem:::.ctBackendNpar(spec), ctsem:::.ctBackendNpar(free) - 1L)
  # The engine builds exactly the stated sd, whatever the raw vector holds.
  handle <- structure(spec, class = c("ctJuliaModel", "ctFitModel"))
  module <- ctsem:::.ctJuliaModule(spec$project)
  popsd <- function(raw, level) sqrt(as.numeric(ctsem:::.ctBackendJuliaValue(
    module$ctsem_laplace_popcov(ctsem:::.ctJuliaObjective(handle),
      ctsem:::.ctJuliaNumericVector(raw), as.integer(level)))))
  npar <- ctsem:::.ctBackendNpar(spec)
  expect_equal(popsd(numeric(npar), 2L), 0.3, tolerance = 1e-8)
  expect_equal(popsd(rep(2, npar), 2L), 0.3, tolerance = 1e-8)
  # The subject level is still free, and moves with its own raw entry.
  expect_false(isTRUE(all.equal(popsd(numeric(npar), 1L), popsd(rep(2, npar), 1L))))
})

test_that("a level poprank reduces refuses a stated value", {
  skip_without_julia()
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = "Y1", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(-0.5), DIFFUSION = matrix(.5), MANIFESTVAR = matrix(.3),
    MANIFESTMEANS = matrix("mm"), T0VAR = matrix(1), T0MEANS = matrix(0),
    CINT = matrix("cint"))))
  m$pars$indvarying <- m$pars$param %in% c("mm", "cint")
  m <- ctsem:::.ctModelRawPopVarSync(m)
  m$RAWPOPVAR["mm", "mm"] <- 0.3
  expect_error(suppressMessages(ctFit(.rpv_data()[, c("id", "time", "Y1")], m,
    backend = "julia", intoverpop = "laplace", fit = FALSE, poprank = 1)),
    "RAWPOPVAR")
})
