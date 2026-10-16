# Every help page parses to its end.
#
# A single unmatched quote inside \code{} -- `\code{s' cov s}` in ctLOO's
# roxygen -- makes the Rd parser read the rest of the page as one string, so
# its Examples and References vanish from the help, and nothing else that runs
# on a push notices: roxygen writes the page without complaint and the package
# still loads. The parser warns, and that warning is what this turns into a
# failure. Cheap and always runs: it reads the Rd sources, it fits nothing.

test_that("every Rd file in man/ parses without a warning", {
  mandir <- normalizePath(file.path(testthat::test_path(), "..", "..", "man"),
    mustWork = FALSE)
  skip_if(!dir.exists(mandir), "package sources not available")
  files <- list.files(mandir, pattern = "[.]Rd$", full.names = TRUE)
  expect_gt(length(files), 50L)   # found the help pages, not nothing
  problems <- character(0)
  for (f in files) {
    withCallingHandlers(tools::parse_Rd(f), warning = function(w) {
      problems <<- c(problems, paste0(basename(f), ": ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  }
  expect_identical(problems, character(0))
})
