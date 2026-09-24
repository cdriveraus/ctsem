# ctsem_laplace_hessian(warm = TRUE) against warm = FALSE, through the engine.
#   Rscript engine-check.R <tree> <lib> <model> <seed>
args <- commandArgs(TRUE)
tree <- path.expand(args[1]); lib <- path.expand(args[2]); MODEL <- args[3]
SEED <- as.integer(args[4])
Sys.setenv(NOT_CRAN = "true")
.libPaths(c(lib, .libPaths()))
suppressMessages(devtools::load_all(tree, compile = FALSE, quiet = TRUE))
source(file.path(tree, "dev/stochopt/models.R"))
P <- make_problem(MODEL, SEED)
fit <- suppressWarnings(suppressMessages(ctFit(P$d, P$m, backend = "julia",
  intoverpop = "laplace", cores = 1, optimcontrol = list(estonly = TRUE))))
x <- fit$estimate$raw; npar <- length(x)
L <- .ctJuliaObjective(fit)
H <- function(warm) {
  t0 <- proc.time()[["elapsed"]]
  h <- JuliaConnectoR::juliaCall("ContinuousTimeSEM.ctsem_laplace_hessian", L,
    .ctJuliaVector(x), warm = warm)
  list(H = matrix(as.numeric(.ctBackendJuliaValue(h)), npar, npar),
    s = proc.time()[["elapsed"]] - t0)
}
invisible(H(TRUE)); invisible(H(FALSE))        # compile both
w <- H(TRUE); c0 <- H(FALSE)
cat(sprintf("ENGINE model=%s seed=%d cold %.1fs warm %.1fs (%.2fx) max|dH|/max|H| %.2e finite %s\n",
  MODEL, SEED, c0$s, w$s, c0$s / w$s,
  max(abs(w$H - c0$H)) / max(abs(c0$H)), all(is.finite(w$H))))
cat("DONE\n")
