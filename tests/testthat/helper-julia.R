# Every backend test needs a live Julia session, and there is no point running
# any of them without one -- the engine *is* what they test. Kept in one place
# so the skip condition cannot drift between files.
#
# The discovery has to run ctsem's own bindir search first. `juliaSetupOk()`
# looks for Julia on the PATH, and on Windows a juliaup install is behind a
# shim that is not there, so a bare probe reports "Julia could not be found" on
# a machine where `ctFit(backend = 'julia')` works perfectly well. Skipping on
# that answer is how this suite came to be a no-op on Windows -- silently, since
# a skip is not a failure.
skip_without_julia <- function() {
  testthat::skip_if_not_installed("JuliaConnectoR")
  # Same search `ctJuliaSetup()` performs, and for the same reason: it knows
  # about juliaup and about the managed install, and `Sys.which` does not.
  bin <- tryCatch(ctsem:::.ctJuliaBin(NULL), error = function(e) NULL)
  if (!is.null(bin) && nzchar(bin)) Sys.setenv(JULIA_BINDIR = bin)
  testthat::skip_if(
    !isTRUE(tryCatch(JuliaConnectoR::juliaSetupOk(), error = function(e) FALSE)),
    "Julia is not available.")
}
