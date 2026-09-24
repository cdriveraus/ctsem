# Which change moves the stan-julia parity fixture's julia fit?
#   Rscript dev/diag-parity.R <tree> <lib>
args <- commandArgs(TRUE)
tree <- path.expand(args[1]); lib <- path.expand(args[2])
Sys.setenv(NOT_CRAN = "true")
.libPaths(c(lib, .libPaths()))
suppressMessages(devtools::load_all(tree, compile = FALSE, quiet = TRUE))
model <- suppressWarnings(ctModel(
  type = "ct", n.latent = 2, LAMBDA = diag(1, 2),
  MANIFESTVAR = diag(c(.1, .1)), MANIFESTMEANS = matrix(0, 2, 1),
  T0VAR = diag(2),
  T0MEANS = c("t0a||TRUE", "t0b||TRUE"),
  CINT = c("B1||TRUE", "B2||TRUE"),
  DRIFT = matrix(c("auto1", "cross21||TRUE", "cross21||TRUE", "auto2"), 2, 2, byrow = TRUE),
  DIFFUSION = diag(c(.2, .15)),
  n.TDpred = 1, TDpredNames = "dose", TDPREDEFFECT = matrix(c("impulse", 0), 2, 1),
  n.TIpred = 1, TIpredNames = "group", tipredDefault = FALSE))
model$pars$group_effect[model$pars$param == "B1"] <- TRUE
set.seed(21)
data <- data.frame()
for (i in 1:6) data <- rbind(data, data.frame(id = i, time = c(0, .5, 1, 1.5),
  Y1 = rnorm(4, 0, .5), Y2 = rnorm(4, 0, .5), dose = c(0, 1, 0, 1),
  group = rep(rnorm(1), 4)))
for (oc in list(list(), list(restarts = 0L))) {
  jf <- suppressWarnings(suppressMessages(ctFit(data, model = model,
    backend = "julia", priors = FALSE, verbose = 0, optimcontrol = oc)))
  cat("==== optimcontrol:", paste(names(oc), unlist(oc)), "\n")
  cat("loglik", jf$estimate$loglik, "converged", jf$optim$converged,
    "status", jf$uncertainty$certification$status, "\n")
  cat("raw[14:17]", round(jf$estimate$raw[14:17], 3), "\n")
  cat("weak:", paste(jf$identifiability$parameters, collapse = ", "), "\n")
  h <- jf$optim$corrections
  if (length(h)) for (x in h) cat("  correction: predicted", signif(x$predicted, 3),
    "step_gain", signif(x$step_gain, 3), "escaped", isTRUE(x$escaped),
    "futile", isTRUE(x$futile), "resumed", isTRUE(x$resumed), "iterations",
    x$iterations, "\n")
  if (!is.null(jf$optim$restarts)) print(jf$optim$restarts)
}
cat("DONE\n")
