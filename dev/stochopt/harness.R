# Harness for optimiser prototypes (dev/stochopt/stochopt.jl).
#
#   Rscript harness.R <model> <seed> <outdir> [phase]
#
# One model x one dataset per process, so a sweep is many processes and dev1's
# cores run them side by side (processes scale linearly there; threads do not).
# Everything is compared from the SAME start, with carefulfit off throughout:
# the prior warm-up is a basin-selection device and orthogonal to this.
#
# Cost is reported two ways: wall seconds inside this process, and counted
# engine calls converted with per-model primitive timings. The counted figure
# is the one that survives a busy machine.
args <- commandArgs(TRUE)
MODEL <- args[1]; SEED <- as.integer(args[2]); OUT <- args[3]
PHASE <- if (length(args) >= 4) args[4] else "endgame"
LIB <- Sys.getenv("CTSEM_LIB", "~/dev/ctsemlib-stoch")
HERE <- Sys.getenv("STOCHOPT_DIR", "~/stoch")
.libPaths(c(LIB, .libPaths()))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(ctsem))
cat("LOADED", MODEL, SEED, PHASE, "\n"); flush(stdout())

source(file.path(HERE, "models.R"))
P <- make_problem(MODEL, SEED)
cat("DATA", nrow(P$d), "rows", length(unique(P$d$id)), "subjects\n"); flush(stdout())

spec <- suppressWarnings(suppressMessages(ctFit(P$d, P$m, backend = "julia",
  fit = FALSE, intoverpop = P$route, cores = 1)))
npar <- ctsem:::.ctBackendNpar(spec)
start <- ctsem:::.ctJuliaInitialValues(npar, NULL, initsd = .01, spec = spec)
set.seed(SEED)
start <- rnorm(npar, 0, .01)
derived <- try(ctsem:::.ctDataStart(P$d, P$m, spec, npar), silent = TRUE)
if (!inherits(derived, "try-error") && !is.null(derived)) {
  use <- is.finite(derived); start[use] <- derived[use]
}
obj <- ctsem:::.ctJuliaObjective(spec)
JuliaConnectoR::juliaEval(sprintf('include("%s")',
  normalizePath(file.path(HERE, "stochopt.jl"))))
SO <- function(name) JuliaConnectoR::juliaFun(paste0("StochOpt.", name))
jv <- function(x) ctsem:::.ctJuliaNumericVector(as.numeric(x))
cat("PREPARED npar", npar, "\n"); flush(stdout())

# warm every code path once so no timing below includes compilation
invisible(ctsem:::.ctJuliaOptimise(spec, start, optimcontrol = list(maxiter = 3L)))
invisible(JuliaConnectoR::juliaGet(SO("primitive_times")(obj, jv(start), 1L)))
cat("WARM\n"); flush(stdout())

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- expr
  list(value = v, secs = proc.time()[["elapsed"]] - t0)
}

optimise <- function(innergaptol = NULL) {
  oc <- list()
  if (!is.null(innergaptol)) oc$innergaptol <- innergaptol
  r <- timed(ctsem:::.ctJuliaOptimise(spec, start, optimcontrol = oc))
  res <- r$value
  list(secs = r$secs, x = as.numeric(res$minimizer)[seq_len(npar)],
    f = as.numeric(res$maximum_loglik)[1L],
    iterations = as.integer(res$iterations),
    f_calls = as.integer(.ctsemOr(res$f_calls)),
    g_calls = as.integer(.ctsemOr(res$g_calls)),
    stopped_by_gap = isTRUE(res$stopped_by_gap),
    trace = ctsem:::.ctBackendTrace(res$trace))
}
.ctsemOr <- function(x) if (is.null(x)) NA else x

out <- list(model = MODEL, seed = SEED, npar = npar,
  nsubjects = length(unique(P$d$id)), nrows = nrow(P$d), route = P$route)

# ---- baseline: the engine optimiser as a fit runs it, then the certification
# Hessian every fit computes
base <- optimise()
h <- timed(ctsem:::.ctBackendHessianAt(spec, base$x))
Hb <- h$value
gb <- JuliaConnectoR::juliaGet(SO("evalgrad")(
  JuliaConnectoR::juliaCall("StochOpt.Ledger", 1L), obj, jv(base$x)))
gapb <- if (!is.null(Hb)) ctsem:::.ctBackendOptimGap(Hb, as.numeric(gb[[2]])) else NULL
base$hess_secs <- h$secs
base$final_gain <- if (!is.null(gapb)) gapb$gap else NA
base$trace <- base$trace  # kept: predicted gain per iteration
out$baseline <- base
cat(sprintf("BASE it=%d %.1fs + hess %.1fs f=%.6f gain=%.2g\n", base$iterations,
  base$secs, base$hess_secs, base$f, base$final_gain)); flush(stdout())

out$prim <- JuliaConnectoR::juliaGet(SO("primitive_times")(obj, jv(base$x), 3L))
cat("PRIM", paste(names(out$prim), signif(unlist(out$prim), 3)), "\n"); flush(stdout())

if (PHASE == "endgame") {
  runs <- list()
  N <- out$nsubjects
  variants <- list(
    list(curv = "exact"), list(curv = "chord"),
    list(curv = "subset", submax = max(20L, ceiling(N / 8))),
    list(curv = "bhhh"))
  for (tau in c(1, 1e-1, 1e-2)) {
    s1 <- optimise(innergaptol = tau)
    cat(sprintf("STAGE1 tau=%g it=%d %.1fs f=%.6f\n", tau, s1$iterations,
      s1$secs, s1$f)); flush(stdout())
    for (v in variants) {
      if (v$curv == "subset" && N < 40) next
      if (v$curv == "bhhh" && N < 2 * npar) next
      eg <- timed(JuliaConnectoR::juliaGet(SO("endgame")(obj, jv(s1$x),
        curvature = JuliaConnectoR::juliaEval(paste0(":", v$curv)),
        submax = as.integer(.ctsemOr(v$submax) %||% 0L))))
      e <- eg$value
      rec <- list(tau = tau, curvature = v$curv, submax = .ctsemOr(v$submax),
        stage1_secs = s1$secs, stage1_iterations = s1$iterations,
        stage1_f = s1$f, endgame_secs = eg$secs, status = e$status,
        f = e$f, iterations = e$iterations, refreshes = e$refreshes,
        final_gain = e$final_gain, extra_step = e$extra_step,
        grad = e$grad, value = e$value, hess = e$hess, scores = e$scores,
        engine_secs = e$seconds)
      runs[[length(runs) + 1L]] <- rec
      cat(sprintf("EG tau=%g %-6s %s it=%d refresh=%d df=%.2e gain=%.2g %.1f+%.1fs (base %.1f+%.1f)\n",
        tau, v$curv, e$status, e$iterations, e$refreshes, e$f - base$f,
        e$final_gain, s1$secs, eg$secs, base$secs, base$hess_secs)); flush(stdout())
    }
  }
  out$endgame <- runs
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, file.path(OUT, sprintf("%s-%s-%d.rds", PHASE, MODEL, SEED)))
cat("SAVED\n")
