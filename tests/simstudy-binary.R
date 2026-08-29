# Does integrating the binary observation beat linearising it?
#
# The claim under test is narrow and mechanical: the extended-Kalman treatment
# of a binary indicator understates the posterior state covariance, which biases
# process noise low, and the bias grows with the number of binary indicators
# loading on a latent. Earlier single-dataset measurements pointed that way and
# a four-dataset comparison did not separate the methods, so the question is
# whether the difference survives replication.
#
# Design: one latent AR process observed only through binary indicators, varying
# the indicator count, fitted by both backends on identical data. The generating
# values are the answer; what is scored is bias and RMSE of DRIFT and DIFFUSION.
#
# Both backends run on the same generated data, so the comparison is paired --
# a seed that happens to be informative helps both.

Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all("/home/ubuntu/dev/ctsem", compile = FALSE, quiet = TRUE))
library(parallel)

# Deliberately NOT connecting Julia in this process.
#
# `mclapply` forks, and a fork inherits the parent's open JuliaConnectoR socket.
# Several children then write to the same socket and the protocol desynchronises
# -- it fails as `Message type not supported (yet): 5a`, which names neither the
# cause nor the process it happened in. Each worker must own its session, so the
# parent never opens one and every child connects on first use.

TRUE_DRIFT <- -0.3
TRUE_DIFF <- 0.8
TRUE_T0VAR <- 1

make_data <- function(seed, nsubjects, nobs, nindicators, ngauss = 0) {
  invlog <- function(x) exp(x) / (1 + exp(x))
  set.seed(seed)
  gen <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = 1,
    manifestNames = "eta", latentNames = "eta1", LAMBDA = matrix(1),
    DRIFT = matrix(TRUE_DRIFT), DIFFUSION = matrix(TRUE_DIFF),
    MANIFESTVAR = matrix(0.001), T0VAR = matrix(TRUE_T0VAR),
    T0MEANS = matrix(0), CINT = matrix(0), MANIFESTMEANS = matrix(0),
    Tpoints = nobs))
  latent <- ctGenerate(gen, n.subjects = nsubjects, Tpoints = nobs,
    backend = "r")
  d <- data.frame(latent)
  eta <- d$eta
  for (i in seq_len(nindicators)) {
    d[[paste0("b", i)]] <- stats::rbinom(nrow(d), 1, invlog(eta))
  }
  # Optional continuous indicators, to test the mixed case.
  for (i in seq_len(ngauss)) {
    d[[paste0("y", i)]] <- eta + stats::rnorm(nrow(d), 0, 0.5)
  }
  d$eta <- NULL
  d
}

make_model <- function(nindicators, ngauss = 0) {
  names <- c(paste0("b", seq_len(nindicators)),
    if (ngauss > 0) paste0("y", seq_len(ngauss)) else character())
  n <- length(names)
  mvar <- diag(0, n)
  if (ngauss > 0) {
    for (i in (nindicators + 1):n) mvar[i, i] <- "mvar"
  }
  m <- suppressMessages(ctModel(type = "ct", n.latent = 1, n.manifest = n,
    manifestNames = names, latentNames = "eta1", LAMBDA = matrix(1, n, 1),
    MANIFESTMEANS = matrix(0, n, 1), CINT = matrix(0), T0MEANS = matrix(0),
    MANIFESTVAR = mvar))
  m$manifesttype <- c(rep(1L, nindicators), rep(0L, ngauss))
  m$pars$indvarying <- FALSE
  m
}

one_cell <- function(job) {
  d <- make_data(job$seed, job$nsubjects, job$nobs, job$nind, job$ngauss)
  m <- make_model(job$nind, job$ngauss)
  out <- list()
  for (backend in c("julia", "stan")) {
    started <- Sys.time()
    fit <- try(suppressWarnings(suppressMessages(
      ctFit(d, m, backend = backend, cores = 1,
        optimcontrol = list(estonly = TRUE)))), silent = TRUE)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (inherits(fit, "try-error")) {
      out[[backend]] <- data.frame(drift = NA_real_, diffusion = NA_real_,
        t0var = NA_real_, converged = FALSE, secs = elapsed)
      next
    }
    est <- try(summary(fit)$popmeans, silent = TRUE)
    if (inherits(est, "try-error") || !"drift_eta1" %in% rownames(est)) {
      out[[backend]] <- data.frame(drift = NA_real_, diffusion = NA_real_,
        t0var = NA_real_, converged = FALSE, secs = elapsed)
      next
    }
    out[[backend]] <- data.frame(
      drift = unname(est["drift_eta1", "mean"]),
      diffusion = unname(est["diff_eta1", "mean"]),
      t0var = if ("T0var_eta1" %in% rownames(est))
        unname(est["T0var_eta1", "mean"]) else NA_real_,
      converged = isTRUE(fit$estimate$converged),
      secs = elapsed)
  }
  cbind(data.frame(seed = job$seed, nind = job$nind, ngauss = job$ngauss,
    nsubjects = job$nsubjects, nobs = job$nobs,
    backend = rep(c("julia", "stan"), each = 1)),
    do.call(rbind, out))
}

jobs <- list()
for (seed in 1:40) {
  for (cell in list(
    list(nind = 3,  ngauss = 0),
    list(nind = 10, ngauss = 0),
    list(nind = 30, ngauss = 0),
    list(nind = 3,  ngauss = 1)   # mixed binary / gaussian
  )) {
    jobs[[length(jobs) + 1L]] <- list(seed = seed, nind = cell$nind,
      ngauss = cell$ngauss, nsubjects = 50, nobs = 12)
  }
}
cat("cells:", length(jobs), "  fits:", 2 * length(jobs), "\n")

started <- Sys.time()
results <- mclapply(jobs, function(j) {
  # One Julia per worker, established inside the fork.
  if (!exists(".jl_ready", envir = globalenv())) {
    suppressMessages(ctJuliaSetup(threads = 1L, force = TRUE))
    assign(".jl_ready", TRUE, envir = globalenv())
  }
  try(one_cell(j), silent = TRUE)
}, mc.cores = 16, mc.preschedule = FALSE)
ok <- !vapply(results, function(x) inherits(x, "try-error"), logical(1))
cat("cells completed:", sum(ok), "of", length(jobs),
    "  in", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
    "minutes\n")
res <- do.call(rbind, results[ok])
saveRDS(res, "/tmp/simstudy.rds")

summarise <- function(df, what, truth) {
  agg <- do.call(rbind, lapply(split(df, list(df$nind, df$ngauss, df$backend),
    drop = TRUE), function(g) {
      v <- g[[what]][g$converged & is.finite(g[[what]])]
      data.frame(nind = g$nind[1], ngauss = g$ngauss[1],
        backend = as.character(g$backend[1]), n = length(v),
        mean = mean(v), bias = mean(v) - truth,
        rmse = sqrt(mean((v - truth)^2)),
        secs = median(g$secs, na.rm = TRUE))
    }))
  agg[order(agg$ngauss, agg$nind, agg$backend), ]
}
cat("\n=== DIFFUSION (true", TRUE_DIFF, ") ===\n")
print(summarise(res, "diffusion", TRUE_DIFF), row.names = FALSE, digits = 3)
cat("\n=== DRIFT (true", TRUE_DRIFT, ") ===\n")
print(summarise(res, "drift", TRUE_DRIFT), row.names = FALSE, digits = 3)
cat("\n=== convergence ===\n")
print(table(res$backend, res$converged, res$nind))
