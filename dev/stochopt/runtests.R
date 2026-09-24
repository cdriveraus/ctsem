# Run named test files against one tree and its installed library.
#   Rscript runtests.R <tree> <lib> <file> [<file> ...]
args <- commandArgs(TRUE)
tree <- path.expand(args[1]); lib <- path.expand(args[2]); files <- args[-(1:2)]
Sys.setenv(NOT_CRAN = "true")
.libPaths(c(lib, .libPaths()))
suppressMessages(devtools::load_all(tree, compile = FALSE, quiet = TRUE))
cat("LOADED\n"); flush(stdout())
for (f in files) {
  cat("==== FILE", f, "\n"); flush(stdout())
  r <- testthat::test_file(file.path(tree, "tests", "testthat", f),
    reporter = "summary", stop_on_failure = FALSE)
  d <- as.data.frame(r)
  cat("==== RESULT", f, "tests", nrow(d), "failed", sum(d$failed > 0),
    "errors", sum(d$error), "skipped", sum(d$skipped), "\n"); flush(stdout())
}
cat("DONE\n")
