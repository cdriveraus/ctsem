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

optimise <- function(innergaptol = NULL, precondition = NULL) {
  oc <- list()
  if (!is.null(innergaptol)) oc$innergaptol <- innergaptol
  if (!is.null(precondition)) oc$precondition <- precondition
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
# Cached per model and seed: the baseline does not depend on the prototype,
# and on the laplace model it is most of a run.
BASECACHE <- Sys.getenv("STOCHOPT_BASECACHE", file.path(HERE, "basecache"))
dir.create(BASECACHE, showWarnings = FALSE, recursive = TRUE)
basefile <- file.path(BASECACHE, sprintf("%s-%d.rds", MODEL, SEED))
if (file.exists(basefile)) {
  base <- readRDS(basefile)
} else {
  base <- optimise()
  h <- timed(ctsem:::.ctBackendHessianAt(spec, base$x))
  Hb <- h$value
  gb <- JuliaConnectoR::juliaGet(SO("evalgrad")(
    JuliaConnectoR::juliaCall("StochOpt.Ledger", 1L), obj, jv(base$x)))
  gapb <- if (!is.null(Hb)) ctsem:::.ctBackendOptimGap(Hb, as.numeric(gb[[2]])) else NULL
  base$hess_secs <- h$secs
  base$final_gain <- if (!is.null(gapb)) gapb$gap else NA
  # the same optimiser without its diagonal preconditioner
  base$noP <- optimise(precondition = FALSE)
  saveRDS(base, basefile)
}
out$baseline <- base
cat(sprintf("BASE it=%d %.1fs + hess %.1fs f=%.6f gain=%.2g | noP it=%d %.1fs f=%.6f\n",
  base$iterations, base$secs, base$hess_secs, base$f, base$final_gain,
  base$noP$iterations, base$noP$secs, base$noP$f)); flush(stdout())
metric <- ctsem:::.ctJuliaParameterScale(spec, at = start, npar = npar)
metric <- if (is.null(metric)) rep(1, npar) else metric^2
metric[!is.finite(metric) | metric <= 0] <- 1

out$prim <- JuliaConnectoR::juliaGet(SO("primitive_times")(obj, jv(base$x), 3L))
cat("PRIM", paste(names(out$prim), signif(unlist(out$prim), 3)), "\n"); flush(stdout())

if (PHASE == "endgame") {
  runs <- list()
  N <- out$nsubjects
  variants <- list(
    list(curv = "exact"), list(curv = "chord"),
    list(curv = "subset", submax = max(20L, ceiling(N / 8))),
    list(curv = "bhhh"))
  for (tau in c(1, 1e-1)) {
    s1 <- optimise(innergaptol = tau)
    cat(sprintf("STAGE1 tau=%g it=%d %.1fs f=%.6f\n", tau, s1$iterations,
      s1$secs, s1$f)); flush(stdout())
    for (v in variants) {
      if (v$curv == "subset" && N < 40) next
      if (v$curv == "bhhh" && N < 2 * npar) next
      eg <- timed(JuliaConnectoR::juliaGet(SO("endgame")(obj, jv(s1$x),
        curvature = JuliaConnectoR::juliaEval(paste0(":", v$curv)),
        submax = if (is.null(v$submax)) 0L else as.integer(v$submax))))
      e <- eg$value
      rec <- list(tau = tau, curvature = v$curv, submax = .ctsemOr(v$submax),
        stage1_secs = s1$secs, stage1_iterations = s1$iterations,
        stage1_f = s1$f, endgame_secs = eg$secs, status = e$status,
        f = e$f, iterations = e$iterations, refreshes = e$refreshes,
        final_gain = e$final_gain, extra_step = e$extra_step,
        grad = e$grad, value = e$value, hess = e$hess, scores = e$scores,
        engine_secs = e$seconds, x = unlist(e$x))
      runs[[length(runs) + 1L]] <- rec
      cat(sprintf("EG tau=%g %-6s %s it=%d refresh=%d df=%.2e gain=%.2g %.1f+%.1fs (base %.1f+%.1f)\n",
        tau, v$curv, e$status, e$iterations, e$refreshes, e$f - base$f,
        e$final_gain, s1$secs, eg$secs, base$secs, base$hess_secs)); flush(stdout())
    }
  }
  out$endgame <- runs
}

if (PHASE == "pbatch") {
  runs <- list()
  sym <- function(s) JuliaConnectoR::juliaEval(paste0(":", s))
  for (usemetric in c(FALSE, TRUE))
  for (grow in c(TRUE, FALSE)) for (theta in if (grow) c(0.25, 0.5) else 1) {
    pb <- timed(JuliaConnectoR::juliaGet(SO("pbatch")(obj, jv(start),
      theta = theta, grow = grow, tol_switch = 0.1,
      metric = jv(if (usemetric) metric else rep(1, npar)))))
    p <- pb$value
    eg <- timed(JuliaConnectoR::juliaGet(SO("endgame")(obj, jv(p$x),
      curvature = sym("chord"))))
    e <- eg$value
    rec <- list(metric = usemetric, grow = grow, theta = theta, pb_status = p$status,
      pb_iterations = p$iterations, final_batch = p$final_batch,
      growth = unlist(p$growth)[-1], sizes = unlist(p$sizes)[-1],
      pb_grad = p$grad, pb_value = p$value, pb_scores = p$scores,
      pb_secs = pb$secs, eg_status = e$status, eg_iterations = e$iterations,
      f = e$f, final_gain = e$final_gain, eg_grad = e$grad,
      eg_value = e$value, eg_hess = e$hess, eg_secs = eg$secs, x = unlist(e$x))
    runs[[length(runs) + 1L]] <- rec
    cat(sprintf("PB metric=%s grow=%s theta=%g %s it=%d batch=%d sizes=[%s] scores=%.1f val=%.1f %.1fs | EG %s it=%d df=%.2e gain=%.2g %.1fs (base %.1f+%.1f)\n",
      usemetric, grow, theta, p$status, p$iterations, p$final_batch,
      paste(unlist(p$sizes)[-1], collapse = ","), p$scores, p$value, pb$secs,
      e$status, e$iterations, e$f - base$f, e$final_gain, eg$secs,
      base$secs, base$hess_secs)); flush(stdout())
  }
  out$pbatch <- runs
}

if (PHASE == "lbfgsonly") {
  # No Newton endgame: the prototype L-BFGS run to a tight predicted gain, then
  # the one certification Hessian every fit pays. For routes where a Hessian
  # is many gradients (laplace: finite differences), this is the candidate.
  runs <- list()
  for (usemetric in c(FALSE, TRUE)) for (grow in c(TRUE, FALSE)) {
    pb <- timed(JuliaConnectoR::juliaGet(SO("pbatch")(obj, jv(start),
      theta = 0.25, grow = grow, tol_switch = 1e-6, maxit = 3000L,
      metric = jv(if (usemetric) metric else rep(1, npar)))))
    p <- pb$value
    h <- timed(ctsem:::.ctBackendHessianAt(spec, unlist(p$x)))
    g <- JuliaConnectoR::juliaGet(SO("evalgrad")(
      JuliaConnectoR::juliaCall("StochOpt.Ledger", 1L), obj, jv(unlist(p$x))))
    gap <- ctsem:::.ctBackendOptimGap(h$value, as.numeric(g[[2]]))
    f <- as.numeric(g[[1]])
    rec <- list(metric = usemetric, grow = grow, status = p$status,
      iterations = p$iterations, sizes = unlist(p$sizes)[-1], secs = pb$secs,
      hess_secs = h$secs, f = f, gap = gap$gap, x = unlist(p$x))
    runs[[length(runs) + 1L]] <- rec
    cat(sprintf("LO metric=%s grow=%s %s it=%d %.1fs + hess %.1fs df=%.2e gap=%.2g (base %.1f+%.1f)\n",
      usemetric, grow, p$status, p$iterations, pb$secs, h$secs, f - base$f,
      gap$gap, base$secs, base$hess_secs)); flush(stdout())
  }
  out$lbfgsonly <- runs
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, file.path(OUT, sprintf("%s-%s-%d.rds", PHASE, MODEL, SEED)))
cat("SAVED\n")
