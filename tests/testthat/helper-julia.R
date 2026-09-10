# Test guards. One spelling per axis, and nothing in tests/ reads an
# environment variable directly.
#
#   skip_without_julia()  this test drives the julia backend
#   skip_on_cran()        testthat's own: this test is too slow for CRAN
#   skip_on_32bit()       the stan models do not build in a 32-bit address space
#
# All three work in two places: inside `test_that()`, skipping that block, and
# at the top of a file before any setup code, skipping the rest of the file.
# Both are reported with their reason ("(code run outside of `test_that()`)"
# for the file-level form), which the `if (Sys.getenv("NOT_CRAN") == "true")`
# wrappers these replaced were not -- a wrapped-out file produced no results at
# all, so a whole file dropping out looked exactly like a file that passed.
#
# Why "needs julia" is one call and not two. The two axes used to be spelled
# three ways -- a file-level NOT_CRAN wrapper, `skip_on_cran()`, and
# `skip_without_julia()` -- with `skip_on_cran()` standing in for "needs julia"
# in some blocks and for "is slow" in others. Six test blocks and two examples
# that should have skipped were erroring instead, because the block they were
# in had only the CRAN guard. A julia test is always both things, so
# `skip_without_julia()` asserts both and is the only call such a test needs:
# add it and delete any `skip_on_cran()` beside it.
#
# An example cannot call any of these. It asks the same question through public
# API instead: wrap the body in `if (isTRUE(ctJuliaStatus()$available))`.

#' Skip unless a live Julia session is available (and this is not CRAN).
#'
#' The CRAN check comes first because the julia probe starts a Julia process.
skip_without_julia <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("JuliaConnectoR")
  # The discovery has to run ctsem's own bindir search first. `juliaSetupOk()`
  # looks for Julia on the PATH, and on Windows a juliaup install is behind a
  # shim that is not there, so a bare probe reports "Julia could not be found"
  # on a machine where `ctFit(backend = 'julia')` works perfectly well. Skipping
  # on that answer is how this suite came to be a no-op on Windows -- silently,
  # since a skip is not a failure.
  bin <- tryCatch(ctsem:::.ctJuliaBin(NULL), error = function(e) NULL)
  if (!is.null(bin) && nzchar(bin)) Sys.setenv(JULIA_BINDIR = bin)
  testthat::skip_if(
    !isTRUE(tryCatch(JuliaConnectoR::juliaSetupOk(), error = function(e) FALSE)),
    "Julia is not available.")
}

#' Skip on a 32-bit build.
#'
#' The stan models are large enough that a 32-bit toolchain runs out of address
#' space compiling them, so every test that fits one carries this beside
#' `skip_on_cran()`.
skip_on_32bit <- function() {
  testthat::skip_if(.Machine$sizeof.pointer == 4, "32-bit build.")
}


# ---------------------------------------------------------------------------
# Fitting one model on both backends
#
# The expensive files in this suite used to fit on stan, which is the default
# backend and the one being deprecated. They fit on julia now and cost a
# fraction of what they did. Setting CTSEM_TEST_STAN fits both and compares
# them, so the stan coverage those files used to give is a run away rather
# than deleted -- and the comparison is stronger than what was there before,
# which never checked one backend against the other at all.
#
#   Rscript -e 'Sys.setenv(NOT_CRAN="true", CTSEM_TEST_STAN="true"); ...'
#
# `test-stan-julia-parity.R` does this too, on small models built for it. The
# difference is scale: these are the suite's large designs -- 50 subjects x 50
# occasions, nonlinear LAMBDA, individually varying loadings, mixed binary and
# Gaussian indicators -- and a disagreement that only appears at that size
# would not show up there.
#
# Nothing below reads the environment except `.ctsem_test_stan()`, for the
# reason at the top of this file.

.ctsem_test_stan <- function() {
  v <- tolower(trimws(Sys.getenv("CTSEM_TEST_STAN", "")))
  v %in% c("true", "yes", "1", "t")
}

#' The backends a test that can run on either should fit.
#'
#' Always 'julia', first, so a test can write `fits$julia` unconditionally.
#' 'stan' as well when CTSEM_TEST_STAN is set.
test_backends <- function() {
  if (.ctsem_test_stan()) c("julia", "stan") else "julia"
}

#' Fit the same model on every backend `test_backends()` asks for.
#'
#' Arguments go to `ctFit()` unchanged, minus `backend`. `stanargs` are extra
#' arguments for the stan fit only, for the few settings that path needs and
#' julia has no equivalent of; passing one to julia is an error there rather
#' than a no-op, which is the house rule for an argument a backend cannot
#' honour.
#'
#' Returns a named list, always with `$julia`, `$stan` only when asked for.
fit_backends <- function(..., backends = test_backends(), stanargs = list()) {
  args <- list(...)
  if (!is.null(args$backend)) stop("fit_backends() chooses the backend")
  out <- list()
  for (be in backends) {
    a <- args
    if (identical(be, "stan")) a <- utils::modifyList(a, stanargs)
    a$backend <- be
    out[[be]] <- do.call(ctsem::ctFit, a)
  }
  out
}

# The summary sections both backends build, and that mean the same thing on
# each. `parmatrices` is deliberately not here: it is compared by the tests
# that care, row by matrix/row/col, because its row *order* differs between
# the two and a whole-table comparison would report that as a difference.
.CT_BACKEND_SECTIONS <- c("popmeans", "popsd", "rawpopcorr")

#' Compare the same model's summary across the backends that were fitted.
#'
#' A no-op with one backend, so a test can call it unconditionally.
#'
#' Compares the `mean` column of each section elementwise by row name, and the
#' log likelihood. Elementwise and by name, not as whole tables: the sections
#' can carry different row orders and different extra columns, and a norm over
#' the lot reports "not equal" without saying which parameter moved.
#'
#' `sections` names what to compare. `rawpopcorr` is the population
#' correlation, and a mismatch there is the first place the population
#' covariance shows up -- see the note in the files that pass
#' `sections = setdiff(.CT_BACKEND_SECTIONS, "rawpopcorr")`.
expect_backends_agree <- function(fits, tol = 1e-2, logliktol = 1e-1,
  sections = .CT_BACKEND_SECTIONS) {

  if (length(fits) < 2L) return(invisible(NULL))
  if (is.null(fits$julia) || is.null(fits$stan)) {
    stop("expect_backends_agree() wants $julia and $stan")
  }
  sj <- summary(fits$julia)
  ss <- summary(fits$stan)

  testthat::expect_equal(as.numeric(sj$loglik), as.numeric(ss$loglik),
    tolerance = logliktol)

  compared <- 0L
  for (sec in sections) {
    a <- sj[[sec]]
    b <- ss[[sec]]
    # A section absent on one side is not a failure -- a model with no random
    # effects has no popsd -- but absent on one and present on the other is.
    if (is.null(a) && is.null(b)) next
    testthat::expect_false(xor(is.null(a), is.null(b)),
      label = paste0(sec, " is present on one backend and not the other"))
    if (is.null(a) || is.null(b)) next

    nm <- intersect(rownames(a), rownames(b))
    # The rows have to line up by name or nothing is being compared, and a
    # comparison of nothing passes. This is the guard for that.
    testthat::expect_gt(length(nm), 0L)
    if (!length(nm)) next
    av <- stats::setNames(as.numeric(a[nm, "mean"]), nm)
    bv <- stats::setNames(as.numeric(b[nm, "mean"]), nm)
    testthat::expect_equal(av, bv, tolerance = tol)
    compared <- compared + length(nm)
  }
  # Same reason again, one level up: every section skipped is a green check
  # that checked nothing.
  testthat::expect_gt(compared, 0L)
  invisible(list(julia = sj, stan = ss))
}
