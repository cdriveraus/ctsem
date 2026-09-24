# Cold against warm-started central-difference Laplace Hessians.
#   Rscript measure.R <tree> <lib> <model> <seed>
# model: ordinal (dev/stochopt/models.R, 300 subjects), panel (the same file's
# 1000-subject panel under intoverpop = 'laplace'), nonlinear (the
# test-julia-laplace.R fixture whose inner problem is multimodal away from the
# mode -- the case the warm start was switched off for).
args <- commandArgs(TRUE)
tree <- path.expand(args[1]); lib <- path.expand(args[2]); MODEL <- args[3]
SEED <- as.integer(args[4])
Sys.setenv(NOT_CRAN = "true")
.libPaths(c(lib, .libPaths()))
suppressMessages(devtools::load_all(tree, compile = FALSE, quiet = TRUE))
if (MODEL == "nonlinear") {
  ex <- parse(file.path(tree, "tests/testthat/test-julia-laplace.R"))
  for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
    startsWith(as.character(e[[2]]), ".laplace_")) eval(e, globalenv())
  set.seed(SEED)
  d <- .laplace_nonlinear_data(); m <- .laplace_nonlinear_model()
} else {
  source(file.path(tree, "dev/stochopt/models.R"))
  P <- make_problem(MODEL, SEED); d <- P$d; m <- P$m
}
fit <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
  intoverpop = "laplace", cores = 1, optimcontrol = list(estonly = TRUE))))
x <- fit$estimate$raw
npar <- length(x)
L <- .ctJuliaObjective(fit)
JuliaConnectoR::juliaEval(sprintf('include("%s")',
  normalizePath(file.path(tree, "dev/laplace-hessian/proto.jl"))))
jv <- function(v) .ctJuliaNumericVector(as.numeric(v))
module <- .ctJuliaModule(fit$model_spec$project)
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr
  list(v = v, s = proc.time()[["elapsed"]] - t0) }
# warm both once so compilation is out of the timings
invisible(JuliaConnectoR::juliaGet(JuliaConnectoR::juliaCall("LaplaceHess.point_costs", L, jv(x), 1L)))
costs <- JuliaConnectoR::juliaGet(JuliaConnectoR::juliaCall("LaplaceHess.point_costs", L, jv(x), 3L))
cold <- tm(matrix(as.numeric(.ctBackendJuliaValue(module$ctsem_laplace_hessian(L,
  .ctJuliaVector(x)))), npar, npar))
warmr <- tm(JuliaConnectoR::juliaGet(JuliaConnectoR::juliaCall(
  "LaplaceHess.warm_hessian", L, jv(x))))
warm <- matrix(as.numeric(warmr$v$H), npar, npar)
ok <- is.finite(cold$v) & is.finite(warm)
rel <- max(abs(warm[ok] - cold$v[ok])) / max(abs(cold$v[ok]))
se <- function(H) { e <- try(sqrt(diag(solve(-H))), silent = TRUE)
  if (inherits(e, "try-error")) rep(NA, nrow(H)) else e }
cat(sprintf("RESULT model=%s seed=%d npar=%d N=%d | gradient cold %.3fs warm %.3fs value %.3fs | inner its/unit cold %.1f warm %.1f | Hessian cold %.1fs warm %.1fs (%.2fx) = %.1f / %.1f gradients | max|dH|/max|H| %.2e | SE max rel diff %.2e | nonfinite cold %d warm %d | warm converged %s\n",
  MODEL, SEED, npar, length(unique(d$id)), costs$cold, costs$warm, costs$value,
  costs$cold_iterations, costs$warm_iterations, cold$s, warmr$s, cold$s / warmr$s,
  cold$s / costs$cold, warmr$s / costs$cold, rel,
  max(abs(se(warm) / se(cold$v) - 1), na.rm = TRUE), sum(!is.finite(cold$v)),
  sum(!is.finite(warm)), warmr$v$converged))
cat("DONE\n")
