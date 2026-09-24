# Does the engine's diagonal preconditioner slow its L-BFGS? Optim turns off
# the s'y/y'y initial scaling whenever a preconditioner is supplied.
#   Rscript precond.R <model> <seed>
args <- commandArgs(TRUE)
MODEL <- args[1]; SEED <- as.integer(args[2])
HERE <- Sys.getenv("STOCHOPT_DIR", "~/stoch/v2")
.libPaths(c(Sys.getenv("CTSEM_LIB", "~/dev/ctsemlib-stoch"), .libPaths()))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(ctsem))
source(file.path(HERE, "models.R"))
P <- make_problem(MODEL, SEED)
spec <- suppressWarnings(suppressMessages(ctFit(P$d, P$m, backend = "julia",
  fit = FALSE, intoverpop = P$route, cores = 1)))
npar <- ctsem:::.ctBackendNpar(spec)
set.seed(SEED)
start <- rnorm(npar, 0, .01)
derived <- try(ctsem:::.ctDataStart(P$d, P$m, spec, npar), silent = TRUE)
if (!inherits(derived, "try-error") && !is.null(derived)) {
  use <- is.finite(derived); start[use] <- derived[use]
}
invisible(ctsem:::.ctJuliaOptimise(spec, start, optimcontrol = list(maxiter = 3L)))
run <- function(oc, label) {
  t0 <- proc.time()[["elapsed"]]
  r <- ctsem:::.ctJuliaOptimise(spec, start, optimcontrol = oc)
  cat(sprintf("%-8s %-16s it=%4d f_calls=%4s %6.1fs f=%.6f\n", MODEL, label,
    as.integer(r$iterations), format(r$f_calls), proc.time()[["elapsed"]] - t0,
    as.numeric(r$maximum_loglik)[1])); flush(stdout())
}
for (rep in 1:2) {
  run(list(), "default")
  run(list(precondition = FALSE), "precond=FALSE")
  run(list(precondition = FALSE, lbfgs_memory = 10L), "noP,m=10")
}
cat("DONE\n")
