# What the quadrature correction costs when a unit carries four or five random
# effects, against the step correction and against no correction.
#
#   Rscript dev/lapcontinue/highdim.R <tree> <outdir> <cores> [k ...]
#
# The model: two latents, two indicators, irregular continuous-time intervals,
# 150 subjects; random DRIFT[1,1] (through -log1p_exp, so the integrand is not
# Gaussian in it), random CINT on both latents, and random T0MEANS on one
# latent (k = 4) or both (k = 5). The data are simulated here, by
# Euler-Maruyama on a fine grid, not by ctGenerate.
#
# Per k: one Laplace fit (laplace_correct = FALSE); the step correction and the
# quadrature correction applied to it post hoc, timed; the engine's own counts
# for the quadrature correction; the time of one 5-node quadrature value at the
# Laplace optimum, which is the unit of the step correction's work; and the
# log posterior at each estimate by importance sampling per unit (t4 proposal
# at 1.5 times the Laplace scale, 50000 draws, seed 1 -- the bench's
# `bench_probe_is`, as verify-cases.R uses it for units wider than three).
# Timings are wall seconds on whatever machine runs it, at `cores` threads.
args <- commandArgs(TRUE)
TREE <- normalizePath(args[1], winslash = "/")
OUT <- args[2]
CORES <- as.integer(args[3])
KS <- if (length(args) > 3) as.integer(args[-(1:3)]) else c(4L, 5L)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(NOT_CRAN = "true", JULIA_NUM_THREADS = as.character(CORES))
if (!nzchar(Sys.getenv("CTSEM_JULIA_AGREE"))) Sys.setenv(CTSEM_JULIA_AGREE = "yes")
suppressMessages(devtools::load_all(TREE, compile = FALSE, quiet = TRUE))
stamp <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush(stdout()) }
stamp("LOADED", TREE, "cores", CORES)
invisible(.ctJuliaModule())
JuliaConnectoR::juliaEval(sprintf('Core.eval(ContinuousTimeSEM, :(include(%s))); nothing',
  deparse(file.path(TREE, "dev/optimbench/reference.jl"))))
jc <- function(fn, ...) JuliaConnectoR::juliaCall(paste0("ContinuousTimeSEM.", fn), ...)
jg <- function(x) if (inherits(x, "JuliaProxy")) JuliaConnectoR::juliaGet(x) else x
jv <- function(x) .ctJuliaNumericVector(as.numeric(x))
now <- function() proc.time()[["elapsed"]]

simulate <- function(k, nsub = 150L, seed = 1L) {
  set.seed(seed)
  rows <- lapply(seq_len(nsub), function(i) {
    raw <- stats::rnorm(1, 0.5, 1.2)
    A <- matrix(c(-log1p(exp(-raw)), 0.1, 0.2, -0.6), 2, 2)
    b <- c(0.5, -0.3) + 0.6 * c(stats::rnorm(1), 0.3 * stats::rnorm(1) +
      sqrt(1 - 0.09) * stats::rnorm(1))
    t0 <- c(1, 0) + c(stats::rnorm(1), if (k >= 5L) stats::rnorm(1) else 0)
    n <- sample(10:14, 1)
    times <- cumsum(c(0, stats::rexp(n - 1L, 1 / 0.8) + 0.05))
    eta <- t0 + 0.3 * stats::rnorm(2)
    Y <- matrix(NA_real_, n, 2)
    Y[1, ] <- eta + 0.4 * stats::rnorm(2)
    for (j in seq_len(n)[-1]) {
      h <- (times[j] - times[j - 1L]) / 50
      for (s in 1:50) eta <- eta + (A %*% eta + b) * h + 0.5 * sqrt(h) * stats::rnorm(2)
      Y[j, ] <- eta + 0.4 * stats::rnorm(2)
    }
    data.frame(id = i, time = times, Y1 = Y[, 1], Y2 = Y[, 2])
  })
  do.call(rbind, rows)
}
model <- function(k) {
  m <- suppressWarnings(suppressMessages(ctModel(type = "ct",
    manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
    LAMBDA = diag(2), MANIFESTMEANS = matrix(0, 2, 1),
    DRIFT = matrix(c("drift11|-log1p_exp(-param)", "drift21", "drift12",
      "drift22"), 2, 2),
    CINT = matrix(c("cint1", "cint2"), 2, 1),
    T0MEANS = matrix(c("t0m1", "t0m2"), 2, 1), T0VAR = diag(0.3, 2),
    DIFFUSION = matrix(c("diff1", 0, 0, "diff2"), 2, 2),
    MANIFESTVAR = matrix(c("mvar1", 0, 0, "mvar2"), 2, 2), silent = TRUE)))
  re <- c("drift11", "cint1", "cint2", "t0m1", if (k >= 5L) "t0m2")
  m$pars$indvarying <- m$pars$param %in% re
  m
}
reference <- function(fit, x) {
  obj <- .ctJuliaObjective(fit)
  ev <- jg(jc("ctsem_laplace_evaluate", obj, jv(x), gradient = FALSE))
  prior <- as.numeric(ev$value) - sum(as.numeric(ev$unit_loglik))
  units <- vapply(seq_along(ev$unit_loglik), function(U) as.numeric(jg(jc(
    "bench_probe_is", obj, jv(x), as.integer(U), n = 50000L, seed = 1L,
    inflate = 1.5))$value), 0)
  list(is = sum(units) + prior, laplace = as.numeric(ev$value), units = units)
}
# The rounds' messages, stamped, so the log shows what each round cost.
timed <- function(expr) withCallingHandlers(expr, message = function(m) {
  stamp(trimws(conditionMessage(m))); invokeRestart("muffleMessage") })

for (k in KS) {
  f <- file.path(OUT, sprintf("highdim-k%d.rds", k))
  if (file.exists(f)) { stamp("have", f); next }
  stamp("K", k)
  d <- simulate(k); m <- model(k)
  rec <- list(k = k, cores = CORES, host = Sys.info()[["nodename"]], tree = TREE,
    sha = tryCatch(system2("git", c("-C", TREE, "rev-parse", "HEAD"), stdout = TRUE),
      error = function(e) NA_character_), nsub = length(unique(d$id)), nrow = nrow(d))
  t0 <- now(); set.seed(1)
  off <- suppressWarnings(suppressMessages(ctFit(d, m, backend = "julia",
    intoverpop = "laplace", cores = CORES,
    optimcontrol = list(laplace_correct = FALSE, finishsamples = 100))))
  rec$fit_seconds <- now() - t0
  rec$npar <- length(off$estimate$raw)
  rec$certification <- off$uncertainty$certification[c("status", "gap")]
  stamp(sprintf("k=%d fit %.0fs, npar %d, certification %s", k, rec$fit_seconds,
    rec$npar, rec$certification$status))
  obj <- .ctJuliaObjective(off)
  x0 <- as.numeric(off$estimate$raw)
  t0 <- now(); q <- jg(jc("ctsem_laplace_quadrature", obj, jv(x0), nodes = 5L))
  rec$quadrature_value_seconds <- now() - t0
  t0 <- now(); jc("ctsem_laplace_evaluate", obj, jv(x0), gradient = TRUE)
  rec$laplace_gradient_seconds <- now() - t0
  stamp(sprintf("one 5-node quadrature value %.2fs, one Laplace gradient %.2fs",
    rec$quadrature_value_seconds, rec$laplace_gradient_seconds))
  t0 <- now(); step <- suppressWarnings(.ctLaplaceAutoCorrect(off, cores = CORES))
  rec$step_seconds <- now() - t0
  rec$step_record <- step$laplace$correction
  stamp(sprintf("step %.0fs: %s, %s steps", rec$step_seconds,
    rec$step_record$status, format(rec$step_record$steps)))
  t0 <- now()
  cont <- suppressWarnings(timed(.ctLaplaceContinue(off, cores = CORES, verbose = 1L)))
  rec$quadrature_seconds <- now() - t0
  cr <- cont$laplace$correction
  rec$quadrature_record <- cr
  stamp(sprintf("quadrature %.0fs: %s/%s, %s rounds of %s, flagged %s of %s, soft blocks %s",
    rec$quadrature_seconds, cr$status, cr$continuation %||% "", format(cr$rounds),
    format(cr$attempts), format(cr$flagged), format(cr$units),
    format(cr$soft_blocks)))
  stamp("evaluations:", paste(names(cr$evaluations), format(cr$evaluations),
    sep = "=", collapse = " "))
  points <- list(laplace = x0, step = as.numeric(step$estimate$raw),
    quadrature = as.numeric(cont$estimate$raw))
  rec$estimates <- points
  rec$reference <- list()
  for (nm in names(points)) {
    same <- Filter(function(o) max(abs(points[[o]] - points[[nm]])) == 0,
      names(rec$reference))
    rec$reference[[nm]] <- if (length(same)) rec$reference[[same[1]]] else {
      t0 <- now(); r <- reference(off, points[[nm]]); r$seconds <- now() - t0; r }
    stamp(sprintf("  %-10s IS %.4f  laplace %.4f", nm, rec$reference[[nm]]$is,
      rec$reference[[nm]]$laplace))
  }
  saveRDS(rec, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  stamp("DONE", k)
}
stamp("ALLDONE")
