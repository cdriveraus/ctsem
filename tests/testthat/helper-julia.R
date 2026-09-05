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
