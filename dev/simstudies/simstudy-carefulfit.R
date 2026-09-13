# Is carefulfit (prior warm-up) still earning its place now that the julia
# optimiser has a diagonal preconditioner and initial_alpha = 0.1?
#
# The design is the one that produced the evidence recorded beside `careful`
# in ctJuliaBackend.R -- dev/simstudies/simstudy-individual-differences.R --
# with stan dropped and `carefulfit` added as a factor. That study is the
# right one to repeat because the cases carefulfit rescued were basin
# selection, not conditioning: a mixed-indicator fit that *converged* to a
# random-effect sd of 7.23 against a truth of 0.5.
#
# What is scored is the population sd of the individual differences, against a
# known truth of 0.5, and the log likelihood reached. A higher likelihood in
# the wrong basin is still the wrong answer, so the sd is the primary outcome.

# Run from the package root (Rscript dev/simstudies/<file>), or set CTSEM_TREE
# to the package directory. Not part of the package build or its tests.
Sys.setenv(NOT_CRAN = "true")
# Set JULIA_BINDIR here if ctsem cannot find Julia on the machine.
suppressMessages(devtools::load_all(Sys.getenv("CTSEM_TREE", "."),
  compile = FALSE, quiet = TRUE))
library(parallel)
cat("LOADED", as.character(packageVersion("ctsem")), "\n"); flush(stdout())

# Julia is deliberately not connected in this process: mclapply forks, and a
# fork inherits the parent's open socket. Several children then write to the
# same one and the protocol desynchronises. Each worker opens its own on first
# use. `mclapply` does not fork on Windows; this was run on dev1.

TRUE_DRIFT <- -0.3
TRUE_DIFF  <- 0.8
TRUE_CINTSD <- 0.5
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
  args <- list(datalong = d, ctstanmodel = m, cores = 1, backend = "julia",
    intoverpop = job$method,
    optimcontrol = list(estonly = TRUE, carefulfit = job$careful))
  started <- Sys.time()
  f <- try(suppressWarnings(suppressMessages(do.call(ctFit, args))),
    silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  base <- data.frame(seed = job$seed, measure = job$measure,
    method = job$method, careful = job$careful, secs = secs,
    stringsAsFactors = FALSE)
  fail <- cbind(base, drift = NA_real_, diffusion = NA_real_,
    cintsd = NA_real_, ll = NA_real_, iterations = NA_integer_,
    converged = FALSE)
  if (inherits(f, "try-error")) return(fail)
  s <- try(summary(f), silent = TRUE)
  if (inherits(s, "try-error")) return(fail)
  ll <- suppressWarnings(as.numeric(f$estimate$loglik))
  if (!length(ll)) ll <- NA_real_
  pm <- s$popmeans
  cbind(base,
    drift = grab(pm, "drift_eta1"),
    diffusion = grab(pm, "diff_eta1"),
    cintsd = grab_sd(s),
    ll = ll[1],
    iterations = as.integer(f$estimate$iterations)[1],
    converged = isTRUE(f$estimate$converged))
}

jobs <- list()
for (seed in seq_len(NREP)) for (measure in names(MEASURES))
  for (method in c("augmented", "laplace")) for (careful in c(TRUE, FALSE))
    jobs[[length(jobs) + 1L]] <- list(seed = seed, measure = measure,
      method = method, careful = careful)
set.seed(1); jobs <- jobs[sample.int(length(jobs))]
cat("cells:", length(jobs), "\n"); flush(stdout())

started <- Sys.time()
results <- mclapply(jobs, function(j) {
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
saveRDS(res, "simstudy_carefulfit.rds")

report <- function(what, truth) {
  cat("\n=== ", what, " (true ", truth, ") ===\n", sep = "")
  agg <- do.call(rbind, lapply(split(res, list(res$measure, res$method,
    res$careful), drop = TRUE), function(g) {
      v <- g[[what]][g$converged & is.finite(g[[what]])]
      data.frame(measure = g$measure[1], method = g$method[1],
        careful = g$careful[1], n = length(v), mean = mean(v),
        bias = mean(v) - truth, rmse = sqrt(mean((v - truth)^2)),
        worst = if (length(v)) max(abs(v - truth)) else NA_real_,
        secs = stats::median(g$secs, na.rm = TRUE),
        iters = stats::median(g$iterations, na.rm = TRUE))
    }))
  agg <- agg[order(agg$measure, agg$method, agg$careful), ]
  print(agg, row.names = FALSE, digits = 3)
}
report("cintsd", TRUE_CINTSD)
report("drift", TRUE_DRIFT)
report("diffusion", TRUE_DIFF)

cat("\n=== convergence ===\n")
print(with(res, table(measure, method, careful, converged)))

# Paired on the same data: the whole question is whether turning the warm-up
# off changes where a *particular* fit lands, which an aggregate can hide.
cat("\n=== paired, same seed and cell: careful minus plain ===\n")
key <- paste(res$seed, res$measure, res$method)
on <- res[res$careful, ]; off <- res[!res$careful, ]
m <- merge(on, off, by = c("seed", "measure", "method"),
  suffixes = c(".on", ".off"))
m$dll <- m$ll.on - m$ll.off
m$derr <- abs(m$cintsd.off - TRUE_CINTSD) - abs(m$cintsd.on - TRUE_CINTSD)
pair <- do.call(rbind, lapply(split(m, list(m$measure, m$method), drop = TRUE),
  function(g) data.frame(measure = g$measure[1], method = g$method[1],
    n = nrow(g),
    ll_on_better = sum(g$dll > 0.01, na.rm = TRUE),
    ll_off_better = sum(g$dll < -0.01, na.rm = TRUE),
    max_ll_gain = max(g$dll, na.rm = TRUE),
    min_ll_gain = min(g$dll, na.rm = TRUE),
    sd_on_closer = sum(g$derr > 0.01, na.rm = TRUE),
    sd_off_closer = sum(g$derr < -0.01, na.rm = TRUE),
    worst_off_err = max(abs(g$cintsd.off - TRUE_CINTSD), na.rm = TRUE),
    worst_on_err = max(abs(g$cintsd.on - TRUE_CINTSD), na.rm = TRUE),
    secs_on = stats::median(g$secs.on), secs_off = stats::median(g$secs.off))))
print(pair, row.names = FALSE, digits = 3)

cat("\n=== the individual fits where they disagree most on cintsd ===\n")
m$absdiff <- abs(m$cintsd.on - m$cintsd.off)
worst <- m[order(-m$absdiff), c("seed", "measure", "method", "cintsd.on",
  "cintsd.off", "ll.on", "ll.off", "converged.on", "converged.off")]
print(utils::head(worst, 20), row.names = FALSE, digits = 4)

cat("SIMDONE\n")
