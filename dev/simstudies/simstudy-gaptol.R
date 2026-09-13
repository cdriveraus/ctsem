# Where should the optimiser's cheap stopping rule sit relative to the bar the
# certification uses?
#
# `.ctBackendInnerGapTol()` puts it at `converge_tol / 100` -- two orders
# inside -- because the line search's `1/2 g'Bg` is a limited-memory estimate
# of the same quantity the exact `1/2 g'H^-1 g` measures, and aiming *at* the
# bar makes every slightly optimistic proxy fail the check and take a
# correction. A correction is a Hessian plus a resumed optimisation. The
# recorded measurement behind the constant is one fixture: 2729 objective calls
# through the optimiser against 5197 through corrections.
#
# One fixture, and taken before the diagonal preconditioner and the short first
# step. This sweeps the constant over four settings on 120 cells, paired on the
# same data AND the same starting values, and scores the thing the constant is
# a trade between: total objective calls against where the fit lands.
#
#   innergaptol = 1e-8   the default, converge_tol/100
#                 1e-7   converge_tol/10
#                 1e-6   at the bar
#                 0      off -- runs to floating-point stagnation, and is the
#                        precision reference every other row is scored against
#
# Certification must be ON for any of this to mean anything: `estonly = TRUE`
# and `certify = FALSE` both set the inner rule to zero, which is the condition
# being measured rather than a setting to pass through.
#
# Run from the package root (Rscript dev/simstudies/<file>), or set CTSEM_TREE
# to the package directory. Not part of the package build or its tests.
# Set CTSEM_LIB instead to use an installed build -- which is what a dev1 run
# does, since nothing here needs an unexported name.
Sys.setenv(NOT_CRAN = "true")
if (nzchar(Sys.getenv("CTSEM_LIB"))) {
  .libPaths(c(Sys.getenv("CTSEM_LIB"), .libPaths()))
  suppressMessages(library(ctsem))
} else {
  suppressMessages(devtools::load_all(Sys.getenv("CTSEM_TREE", "."),
    compile = FALSE, quiet = TRUE))
}
library(parallel)
cat("LOADED", as.character(packageVersion("ctsem")), "\n"); flush(stdout())

# Julia is deliberately not connected in this process: mclapply forks, and a
# fork inherits the parent's open socket. Each worker opens its own on first
# use. `mclapply` does not fork on Windows; this was run on dev1.

TRUE_DRIFT <- -0.3
TRUE_DIFF  <- 0.8
TRUE_CINTSD <- 0.5
TAU <- c(-1.0, 0.4, 1.9)
NSUB <- 60
NOBS <- 10
NREP <- 15
SETTINGS <- c(default = 1e-8, tenth = 1e-7, atbar = 1e-6, off = 0)

invlog <- function(x) 1 / (1 + exp(-x))
drawcat <- function(eta, tau) {
  cum <- sapply(seq_along(tau), function(k) invlog(tau[k] - eta))
  p <- cbind(cum, 1)
  p <- cbind(p[, 1, drop = FALSE], t(apply(p, 1, diff)))
  apply(p, 1, function(pr) sample.int(length(pr), 1, prob = pmax(pr, 0)))
}

MEASURES <- list(
  gaussian = list(names = c("y1", "y2"), type = c(0L, 0L), ncat = c(0L, 0L)),
  binary   = list(names = c("b1", "b2", "b3"), type = c(1L, 1L, 1L),
    ncat = c(0L, 0L, 0L)),
  ordinal  = list(names = c("o1", "o2", "o3"), type = c(2L, 2L, 2L),
    ncat = c(4L, 4L, 4L)),
  mixed    = list(names = c("o1", "b1", "y1"), type = c(2L, 1L, 0L),
    ncat = c(4L, 0L, 0L)))

make_data <- function(seed, measure) {
  set.seed(seed)
  spec <- MEASURES[[measure]]
  cints <- stats::rnorm(NSUB, 0, TRUE_CINTSD)
  d <- do.call(rbind, lapply(seq_len(NSUB), function(i) {
    gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
      manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
      DRIFT = matrix(TRUE_DRIFT), DIFFUSION = matrix(TRUE_DIFF),
      MANIFESTVAR = matrix(1e-6), T0VAR = matrix(1), T0MEANS = matrix(0),
      CINT = matrix(cints[i]), MANIFESTMEANS = matrix(0), Tpoints = NOBS))
    one <- data.frame(ctGenerate(gen, n.subjects = 1, Tpoints = NOBS,
      backend = "r"))
    one$id <- i
    one
  }))
  eta <- d$eta
  for (j in seq_along(spec$names)) {
    d[[spec$names[j]]] <- switch(as.character(spec$type[j]),
      "0" = eta + stats::rnorm(length(eta), 0, 0.5),
      "1" = stats::rbinom(length(eta), 1, invlog(eta)),
      "2" = drawcat(eta, TAU))
  }
  d$eta <- NULL
  d
}

make_model <- function(measure) {
  spec <- MEASURES[[measure]]
  n <- length(spec$names)
  mvar <- diag(0, n)
  for (i in seq_len(n)) if (spec$type[i] == 0L) mvar[i, i] <- "mvar"
  args <- list(type = "ct", n.latent = 1, n.manifest = n,
    manifestNames = spec$names, latentNames = "eta1",
    LAMBDA = matrix(1, n, 1), MANIFESTMEANS = matrix(0, n, 1),
    CINT = matrix("cint"), T0MEANS = matrix(0), MANIFESTVAR = mvar,
    manifesttype = spec$type)
  if (any(spec$type == 2L)) args$ncategories <- spec$ncat
  m <- suppressWarnings(suppressMessages(do.call(ctModel, args)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$param %in% "cint"] <- TRUE
  m
}

number <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x)) return(default)
  v <- suppressWarnings(as.numeric(x)[1L])
  if (is.finite(v)) v else default
}

# One cell is one dataset fitted under every setting, from the same starting
# values: `ctFit(inits = NULL)` draws `rnorm(npar, 0, .01)`, so seeding
# immediately before each call is what makes the four comparable. Passing
# `inits` instead would have turned `carefulfit` off, which is a different fit.
one_cell <- function(job) {
  d <- make_data(job$seed, job$measure)
  m <- make_model(job$measure)
  rows <- lapply(names(SETTINGS), function(setting) {
    set.seed(job$seed * 1000L + 7L)
    started <- Sys.time()
    f <- try(suppressWarnings(suppressMessages(ctFit(datalong = d,
      ctstanmodel = m, cores = 1, backend = "julia", intoverpop = job$method,
      optimcontrol = list(innergaptol = unname(SETTINGS[[setting]]))))),
      silent = TRUE)
    secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    base <- data.frame(seed = job$seed, measure = job$measure,
      method = job$method, setting = setting,
      innergaptol = unname(SETTINGS[[setting]]), secs = secs,
      stringsAsFactors = FALSE)
    if (inherits(f, "try-error")) {
      return(list(row = cbind(base, ll = NA_real_, iterations = NA_integer_,
        f_calls = NA_integer_, g_calls = NA_integer_, corrections = NA_integer_,
        hessians = NA_integer_, stopped_by_gap = NA, gap = NA_real_,
        certified = NA, converged = FALSE), raw = NULL))
    }
    e <- f$estimate
    certification <- f$uncertainty$certification
    list(row = cbind(base,
      ll = number(e$loglik),
      iterations = as.integer(number(e$iterations)),
      # Totals: when a correction ran, ctFit replaces these with optimiser plus
      # correction, which is the quantity the constant trades against. These
      # are counts rather than seconds deliberately -- wall time here is mostly
      # the uncertainty phase, which the stopping rule does not touch.
      f_calls = as.integer(number(e$f_calls)),
      g_calls = as.integer(number(e$g_calls)),
      # A list of attempts, so its length is the count.
      corrections = length(e$corrections),
      # Each correction costs one; the certification itself costs the first.
      hessians = as.integer(number(e$hessians)),
      stopped_by_gap = isTRUE(e$stopped_by_gap),
      gap = number(certification$gap),
      certified = isTRUE(certification$certified),
      converged = isTRUE(e$converged)),
      raw = as.numeric(e$raw))
  })
  out <- do.call(rbind, lapply(rows, `[[`, "row"))
  # Distance from the setting that optimises furthest, on the raw scale the
  # optimiser actually works on.
  reference <- rows[[match("off", names(SETTINGS))]]$raw
  out$maxrawdiff <- vapply(rows, function(r) {
    if (is.null(r$raw) || is.null(reference) ||
        length(r$raw) != length(reference)) return(NA_real_)
    max(abs(r$raw - reference))
  }, numeric(1))
  out$llgap <- out$ll[out$setting == "off"] - out$ll
  out
}

jobs <- list()
for (seed in seq_len(NREP)) for (measure in names(MEASURES))
  for (method in c("augmented", "laplace"))
    jobs[[length(jobs) + 1L]] <- list(seed = seed, measure = measure,
      method = method)
set.seed(1); jobs <- jobs[sample.int(length(jobs))]
cat("cells:", length(jobs), " fits:", length(jobs) * length(SETTINGS), "\n")
flush(stdout())

started <- Sys.time()
results <- mclapply(jobs, function(j) {
  if (!exists(".jl_ready", envir = globalenv())) {
    suppressMessages(ctJuliaSetup(threads = 1L, force = TRUE))
    assign(".jl_ready", TRUE, envir = globalenv())
  }
  try(one_cell(j), silent = TRUE)
}, mc.cores = 16, mc.preschedule = TRUE)

ok <- !vapply(results, function(x) inherits(x, "try-error"), logical(1))
cat("completed:", sum(ok), "of", length(jobs), "cells in",
  round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
  "minutes\n")
if (any(!ok)) print(unique(vapply(results[!ok], conditionMessage, character(1))))
res <- do.call(rbind, results[ok])
saveRDS(res, "simstudy_gaptol.rds")

res$setting <- factor(res$setting, levels = names(SETTINGS))

cat("\n=== cost: objective calls, iterations, wall time (medians) ===\n")
agg <- do.call(rbind, lapply(split(res, list(res$method, res$setting),
  drop = TRUE), function(g) data.frame(
    method = g$method[1], setting = as.character(g$setting[1]),
    n = nrow(g),
    f_calls = stats::median(g$f_calls, na.rm = TRUE),
    g_calls = stats::median(g$g_calls, na.rm = TRUE),
    iters = stats::median(g$iterations, na.rm = TRUE),
    corrected = sum(g$corrections > 0, na.rm = TRUE),
    hessians = stats::median(g$hessians, na.rm = TRUE),
    stoppedgap = sum(g$stopped_by_gap, na.rm = TRUE),
    secs = stats::median(g$secs, na.rm = TRUE),
    certified = sum(g$certified, na.rm = TRUE),
    converged = sum(g$converged, na.rm = TRUE))))
print(agg[order(agg$method, agg$setting), ], row.names = FALSE, digits = 4)

cat("\n=== what it costs: distance from innergaptol = 0 ===\n")
cat("llgap is (ll at 0) - (ll here); positive means this setting landed lower.\n")
agg2 <- do.call(rbind, lapply(split(res, list(res$method, res$setting),
  drop = TRUE), function(g) data.frame(
    method = g$method[1], setting = as.character(g$setting[1]),
    n = sum(is.finite(g$llgap)),
    median_llgap = stats::median(g$llgap, na.rm = TRUE),
    max_llgap = suppressWarnings(max(g$llgap, na.rm = TRUE)),
    worse_by_gt_0.01 = sum(g$llgap > 0.01, na.rm = TRUE),
    median_maxrawdiff = stats::median(g$maxrawdiff, na.rm = TRUE),
    max_maxrawdiff = suppressWarnings(max(g$maxrawdiff, na.rm = TRUE)),
    median_gap = stats::median(g$gap, na.rm = TRUE),
    max_gap = suppressWarnings(max(g$gap, na.rm = TRUE)))))
print(agg2[order(agg2$method, agg2$setting), ], row.names = FALSE, digits = 4)

cat("\n=== per measurement type, calls and landing ===\n")
agg3 <- do.call(rbind, lapply(split(res, list(res$measure, res$setting),
  drop = TRUE), function(g) data.frame(
    measure = g$measure[1], setting = as.character(g$setting[1]),
    f_calls = stats::median(g$f_calls, na.rm = TRUE),
    g_calls = stats::median(g$g_calls, na.rm = TRUE),
    corrected = sum(g$corrections > 0, na.rm = TRUE),
    max_llgap = suppressWarnings(max(g$llgap, na.rm = TRUE)),
    secs = stats::median(g$secs, na.rm = TRUE))))
print(agg3[order(agg3$measure, agg3$setting), ], row.names = FALSE, digits = 4)

cat("\n=== the cells where a setting landed furthest below the reference ===\n")
worst <- res[is.finite(res$llgap), ]
worst <- worst[order(-worst$llgap), c("seed", "measure", "method", "setting",
  "llgap", "maxrawdiff", "f_calls", "corrections", "converged")]
print(utils::head(worst, 15), row.names = FALSE, digits = 4)

cat("SIMDONE\n")
