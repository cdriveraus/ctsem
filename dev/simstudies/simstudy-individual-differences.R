# Recovering individual differences through a non-Gaussian measurement.
#
# The question is narrow: when subjects differ in the level of the process and
# that process is observed only through categorical indicators, does the fit
# recover the *population spread* of those differences -- and does it matter
# whether the random effects are integrated by state augmentation or by a
# Laplace approximation?
#
# The generating values are the answer. Each subject gets its own CINT drawn
# from N(0, sd), the latent path is simulated from it, and the indicators are
# drawn from the latent path through the appropriate link. What is scored is
# bias and RMSE of DRIFT, DIFFUSION and that population sd, which is the
# parameter the whole exercise is about.
#
# Stan is included where it can run at all: it has a binary measurement but no
# ordinal one, so it appears in two of the four measurement conditions and its
# absence from the others is the point rather than an omission.

# Run from the package root (Rscript dev/simstudies/<file>), or set CTSEM_TREE
# to the package directory. Not part of the package build or its tests.
Sys.setenv(NOT_CRAN = "true")
# Set JULIA_BINDIR here if ctsem cannot find Julia on the machine.
suppressMessages(devtools::load_all(Sys.getenv("CTSEM_TREE", "."), compile = FALSE, quiet = TRUE))
library(parallel)

# Julia is deliberately not connected in this process: mclapply forks, and a
# fork inherits the parent's open socket. Several children then write to the
# same one and the protocol desynchronises, which surfaces as
# "Message type not supported (yet)" and names neither the cause nor the
# process. Each worker opens its own on first use.

TRUE_DRIFT <- -0.3
TRUE_DIFF  <- 0.8
TRUE_CINTSD <- 0.5      # population sd of the individual differences
TAU <- c(-1.0, 0.4, 1.9)
NSUB <- 60
NOBS <- 10
NREP <- 30

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
  # One subject at a time, because each has its own CINT.
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

# The population sd of the random effect is reported in its own section rather
# than among the fixed effects; find it without assuming a name.
grab_sd <- function(s) {
  for (nm in c("popsd", "rawpopsd", "popsds")) {
    tab <- s[[nm]]
    if (!is.null(tab) && nrow(as.matrix(tab)) > 0) {
      mat <- as.matrix(tab)
      hit <- grep("cint", rownames(mat), ignore.case = TRUE)
      if (length(hit)) return(unname(mat[hit[1], 1]))
      return(unname(mat[1, 1]))
    }
  }
  NA_real_
}
grab <- function(mat, name) {
  if (is.null(mat) || !name %in% rownames(mat)) return(NA_real_)
  unname(mat[name, "mean"])
}

one_cell <- function(job) {
  d <- make_data(job$seed, job$measure)
  m <- make_model(job$measure)
  args <- list(datalong = d, ctstanmodel = m, cores = 1,
    optimcontrol = list(estonly = TRUE))
  if (job$method == "stan") {
    args$backend <- "stan"
  } else {
    args$backend <- "julia"
    args$intoverpop <- job$method
  }
  started <- Sys.time()
  f <- try(suppressWarnings(suppressMessages(do.call(ctFit, args))),
    silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  base <- data.frame(seed = job$seed, measure = job$measure,
    method = job$method, secs = secs, stringsAsFactors = FALSE)
  if (inherits(f, "try-error")) {
    return(cbind(base, drift = NA_real_, diffusion = NA_real_,
      cintsd = NA_real_, converged = FALSE))
  }
  s <- try(summary(f), silent = TRUE)
  if (inherits(s, "try-error")) {
    return(cbind(base, drift = NA_real_, diffusion = NA_real_,
      cintsd = NA_real_, converged = FALSE))
  }
  pm <- s$popmeans
  cbind(base,
    drift = grab(pm, "drift_eta1"),
    diffusion = grab(pm, "diff_eta1"),
    cintsd = grab_sd(s),
    converged = isTRUE(f$estimate$converged))
}

jobs <- list()
for (seed in seq_len(NREP)) {
  for (measure in names(MEASURES)) {
    methods <- c("augmented", "laplace")
    # Stan has a binary measurement and no ordinal one.
    if (!any(MEASURES[[measure]]$type == 2L)) methods <- c(methods, "stan")
    for (method in methods) {
      jobs[[length(jobs) + 1L]] <- list(seed = seed, measure = measure,
        method = method)
    }
  }
}
set.seed(1); jobs <- jobs[sample.int(length(jobs))]   # balance the chunks
cat("cells:", length(jobs), "\n"); flush(stdout())

started <- Sys.time()
results <- mclapply(jobs, function(j) {
  # One Julia per worker, established inside the fork. `mc.preschedule = TRUE`
  # keeps it that way: a worker per core handling a contiguous chunk, rather
  # than a fresh process -- and a fresh Julia -- per fit.
  if (!exists(".jl_ready", envir = globalenv())) {
    suppressMessages(ctJuliaSetup(threads = 1L, force = TRUE))
    assign(".jl_ready", TRUE, envir = globalenv())
  }
  try(one_cell(j), silent = TRUE)
}, mc.cores = 16, mc.preschedule = TRUE)

ok <- !vapply(results, function(x) inherits(x, "try-error"), logical(1))
cat("completed:", sum(ok), "of", length(jobs), "in",
  round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
  "minutes\n")
res <- do.call(rbind, results[ok])
saveRDS(res, "simstudy_indiv.rds")

report <- function(what, truth) {
  cat("\n=== ", what, " (true ", truth, ") ===\n", sep = "")
  agg <- do.call(rbind, lapply(split(res, list(res$measure, res$method),
    drop = TRUE), function(g) {
      # `estimate$converged` is a julia field; a stan fit does not carry it,
      # and filtering on it drops every stan row without saying so. Each
      # backend is asked in its own terms.
      conv <- if (g$method[1] == "stan") is.finite(g[[what]]) else g$converged
      v <- g[[what]][conv & is.finite(g[[what]])]
      data.frame(measure = g$measure[1], method = g$method[1],
        n = length(v), mean = mean(v), bias = mean(v) - truth,
        rmse = sqrt(mean((v - truth)^2)),
        secs = stats::median(g$secs, na.rm = TRUE))
    }))
  agg <- agg[order(agg$measure, agg$method), ]
  print(agg, row.names = FALSE, digits = 3)
}
report("cintsd", TRUE_CINTSD)
report("drift", TRUE_DRIFT)
report("diffusion", TRUE_DIFF)

cat("\n=== convergence ===\n")
print(with(res, table(measure, method, converged)))
cat("SIMDONE\n")
