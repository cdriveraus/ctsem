# Where a warm 5000-subject panel fit spends its time, outside the optimiser.
#   Rscript overhead.R <lib> <seed>
args <- commandArgs(TRUE)
LIB <- args[1]; SEED <- as.integer(args[2])
HERE <- Sys.getenv("STOCHOPT_DIR", "~/stoch/v3")
.libPaths(c(LIB, .libPaths()))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(ctsem))
source(file.path(HERE, "models.R"))
P <- make_problem("panel5k", SEED)
el <- function(expr) { t0 <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t0 }
ids <- unique(P$d$id)
slice <- P$d[P$d$id %in% ids[1:500], ]
fitit <- function(d, oc) { set.seed(SEED)
  suppressWarnings(suppressMessages(ctFit(d, P$m, backend = "julia",
    intoverpop = P$route, cores = 1, optimcontrol = oc))) }
invisible(fitit(slice, list(estonly = FALSE)))            # warm
invisible(fitit(slice, list(estonly = FALSE, carefulfit = FALSE)))
t_spec <- el(spec <- suppressWarnings(suppressMessages(ctFit(P$d, P$m,
  backend = "julia", fit = FALSE, intoverpop = P$route, cores = 1))))
obj <- ctsem:::.ctJuliaObjective(spec)
npar <- ctsem:::.ctBackendNpar(spec)
x <- rnorm(npar, 0, .01)
t_grad <- min(sapply(1:3, function(i) el(ctJuliaEvaluate(spec, x, gradient = TRUE))))
t_hess <- el(ctsem:::.ctBackendHessianAt(spec, x))
t_full <- el(f1 <- fitit(P$d, list(estonly = FALSE)))
t_nocare <- el(f2 <- fitit(P$d, list(estonly = FALSE, carefulfit = FALSE)))
cat(sprintf("SPLIT seed=%d spec_build=%.1f gradient=%.2f hessian=%.1f full_fit=%.1f (it=%d) no_carefulfit=%.1f (it=%d) dll=%.2e\n",
  SEED, t_spec, t_grad, t_hess, t_full, f1$optim$iterations, t_nocare,
  f2$optim$iterations, f2$estimate$loglik - f1$estimate$loglik))
