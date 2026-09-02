# Individual differences through a categorical measurement, larger, with the
# convergence behaviour recorded rather than inferred.
#
# Sixty replications of four measurement conditions, fitted by state
# augmentation and by the Laplace approximation, and by stan where stan can run.
# Each fit records its iteration count, final gradient and verdict, so the
# stalls can be counted and characterised instead of anecdotally chased.
#
# Bias is reported with a Monte Carlo standard error, because a bias of 0.03
# from sixty replications is only a claim if it is larger than the noise in the
# estimate of it.

# Run from the package root (Rscript dev/simstudies/<file>), or set CTSEM_TREE
# to the package directory. Not part of the package build or its tests.
Sys.setenv(NOT_CRAN = "true")
# Set JULIA_BINDIR here if ctsem cannot find Julia on the machine.
suppressMessages(devtools::load_all(Sys.getenv("CTSEM_TREE", "."), compile = FALSE, quiet = TRUE))
library(parallel)

TRUE_DRIFT <- -0.3
TRUE_DIFF  <- 0.8
TRUE_CINTSD <- 0.5
TAU <- c(-1.0, 0.4, 1.9)
NSUB <- 60
NOBS <- 10
NREP <- 60
CORES <- 20

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

grab_sd <- function(s) {
  tab <- s$popsd
  if (is.null(tab) || nrow(as.matrix(tab)) == 0) return(NA_real_)
  mat <- as.matrix(tab)
  hit <- grep("cint", rownames(mat), ignore.case = TRUE)
  unname(mat[if (length(hit)) hit[1] else 1, 1])
}
grab <- function(mat, name) {
  if (is.null(mat) || !name %in% rownames(mat)) return(NA_real_)
  unname(mat[name, "mean"])
}
num <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[1]

one_cell <- function(job) {
  d <- make_data(job$seed, job$measure)
  m <- make_model(job$measure)
  args <- list(datalong = d, ctstanmodel = m, cores = 1,
    optimcontrol = list(estonly = TRUE))
  if (job$method == "stan") args$backend <- "stan" else {
    args$backend <- "julia"; args$intoverpop <- job$method
  }
  started <- Sys.time()
  f <- try(suppressWarnings(suppressMessages(do.call(ctFit, args))),
    silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  base <- data.frame(seed = job$seed, measure = job$measure,
    method = job$method, secs = secs, stringsAsFactors = FALSE)
  blank <- cbind(base, drift = NA_real_, diffusion = NA_real_,
    cintsd = NA_real_, loglik = NA_real_, iterations = NA_real_,
    gradnorm = NA_real_, converged = FALSE, errored = TRUE)
  if (inherits(f, "try-error")) return(blank)
  s <- try(summary(f), silent = TRUE)
  if (inherits(s, "try-error")) return(blank)
  e <- f$estimate
  cbind(base,
    drift = grab(s$popmeans, "drift_eta1"),
    diffusion = grab(s$popmeans, "diff_eta1"),
    cintsd = grab_sd(s),
    loglik = num(e$loglik),
    iterations = num(e$iterations),
    gradnorm = num(e$gradient_norm),
    # stan carries no `converged` field; a finite estimate is the criterion
    # available for it, and asking each backend in its own terms is what stops
    # every stan row being silently dropped.
    converged = if (job$method == "stan")
      is.finite(grab(s$popmeans, "drift_eta1")) else isTRUE(e$converged),
    errored = FALSE)
}

jobs <- list()
for (seed in seq_len(NREP)) for (measure in names(MEASURES)) {
  methods <- c("augmented", "laplace")
  if (!any(MEASURES[[measure]]$type == 2L)) methods <- c(methods, "stan")
  for (method in methods) jobs[[length(jobs) + 1L]] <-
    list(seed = seed, measure = measure, method = method)
}
set.seed(1); jobs <- jobs[sample.int(length(jobs))]
cat("cells:", length(jobs), "on", CORES, "cores\n"); flush(stdout())

started <- Sys.time()
results <- mclapply(jobs, function(j) {
  if (!exists(".jl_ready", envir = globalenv())) {
    suppressMessages(ctJuliaSetup(threads = 1L, force = TRUE))
    assign(".jl_ready", TRUE, envir = globalenv())
  }
  try(one_cell(j), silent = TRUE)
}, mc.cores = CORES, mc.preschedule = TRUE)

ok <- !vapply(results, function(x) inherits(x, "try-error"), logical(1))
cat("completed:", sum(ok), "of", length(jobs), "in",
  round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
  "minutes\n")
res <- do.call(rbind, results[ok])
saveRDS(res, "simstudy_indiv2.rds")

cells <- function(f) do.call(rbind, lapply(
  split(res, list(res$measure, res$method), drop = TRUE), f))

cat("\n=== convergence and stalls ===\n")
cv <- cells(function(g) data.frame(measure = g$measure[1], method = g$method[1],
  n = nrow(g), converged = sum(g$converged), errored = sum(g$errored),
  # A stall is the signature seen before: it stopped almost immediately with a
  # gradient nowhere near zero, rather than arriving anywhere.
  stalled = sum(!g$converged & is.finite(g$iterations) & g$iterations < 15 &
      is.finite(g$gradnorm) & g$gradnorm > 1, na.rm = TRUE),
  medianiter = stats::median(g$iterations, na.rm = TRUE)))
print(cv[order(cv$measure, cv$method), ], row.names = FALSE)

report <- function(what, truth) {
  cat("\n=== ", what, "  (true ", truth, ") ===\n", sep = "")
  out <- cells(function(g) {
    v <- g[[what]][g$converged & is.finite(g[[what]])]
    data.frame(measure = g$measure[1], method = g$method[1], n = length(v),
      bias = mean(v) - truth, mcse = stats::sd(v) / sqrt(length(v)),
      rmse = sqrt(mean((v - truth)^2)),
      medAE = stats::median(abs(v - truth)))
  })
  out$bias_z <- out$bias / out$mcse
  print(out[order(out$measure, out$method), ], row.names = FALSE, digits = 3)
}
report("cintsd", TRUE_CINTSD)
report("drift", TRUE_DRIFT)
report("diffusion", TRUE_DIFF)

cat("\n=== timing (seconds); the minimum is the warm figure ===\n")
tm <- cells(function(g) data.frame(measure = g$measure[1],
  method = g$method[1], min = min(g$secs), median = stats::median(g$secs)))
print(tm[order(tm$measure, tm$method), ], row.names = FALSE, digits = 3)
cat("SIM3DONE\n")
