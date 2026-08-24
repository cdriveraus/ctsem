# Every backend test needs a live Julia session, and there is no point running
# any of them without one -- the engine *is* what they test. Kept in one place
# so the skip condition cannot drift between files.
skip_without_julia <- function() {
  testthat::skip_if_not_installed("JuliaConnectoR")
  testthat::skip_if(
    !isTRUE(tryCatch(JuliaConnectoR::juliaSetupOk(), error = function(e) FALSE)),
    "Julia is not available.")
}
