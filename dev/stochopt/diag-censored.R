# Why does the censored recovery fit report converged = FALSE?
#   Rscript diag-censored.R <tree> <lib> [batch] [newton]
args <- commandArgs(TRUE)
tree <- path.expand(args[1]); lib <- path.expand(args[2])
Sys.setenv(NOT_CRAN = "true")
.libPaths(c(lib, .libPaths()))
suppressMessages(devtools::load_all(tree, compile = FALSE, quiet = TRUE))
# the test file's helpers: every top-level assignment to a dotted name
ex <- parse(file.path(tree, "tests/testthat/test-julia-censored.R"))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
  startsWith(as.character(e[[2]]), ".")) eval(e, globalenv())
for (oc in list(list(), list(batch = FALSE), list(newton = FALSE),
  list(batch = FALSE, newton = FALSE))) {
  fit <- suppressWarnings(suppressMessages(ctFit(.censored_data(nsubjects = 50,
    nobs = 10), .censored_model(), backend = "julia", intoverpop = "augmented",
    optimcontrol = c(list(estonly = TRUE), oc))))
  o <- fit$optim
  cat(sprintf("%-28s converged=%s it=%d newton=%d batch=%s gap_stop=%s pg=%.2e last=%.2e |g|=%.2e stalled=%s overshot=%s og=%.2e saturated=%s sat=%s escapes=%d ll=%.6f\n",
    paste(names(oc), unlist(oc), collapse = ","), o$converged, o$iterations,
    o$newton_steps, paste(o$batch_sizes, collapse = ","), o$stopped_by_gap,
    o$predicted_gain, o$last_gain, o$gradient_norm, o$stalled, o$overshot,
    o$overshoot_gain, o$saturated, paste(o$saturated_parameters, collapse = ","),
    o$stall_escapes, fit$estimate$loglik))
}
cat("DONE\n")
