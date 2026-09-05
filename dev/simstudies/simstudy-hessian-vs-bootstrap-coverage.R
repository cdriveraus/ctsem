# Do the two uncertainty routes give intervals with nominal coverage?
#
# ctFit reports parameter intervals two ways. The default draws finishsamples
# from the surrogate/sampled posterior; `optimcontrol$uncertainty = 'hessian'`
# takes them from the Hessian at the optimum. Both are asymptotic in different
# senses, and the question is empirical: over repeated datasets from a known
# generating model, what fraction of the 95% intervals contain the generating
# value, for each route and each parameter?
#
# Design: one latent AR process, 200 subjects, individually varying T0MEANS and
# CINT with those two correlated in the population, 10 occasions each. Each
# replication fits the same model twice -- once per uncertainty route, on the
# same data, so the comparison is paired -- and scores whether the generating
# value falls inside the reported 2.5%-97.5% interval.
#
# PROVENANCE. This was `tests/testthat/test-bootHessian.R`. It was committed
# with its body inside `if(FALSE)`, later an unconditional `testthat::skip()`,
# and it has never been executed: it is a coverage simulation study, not a unit
# test -- 100 replications x 2 fits of a 200-subject model -- and it asserted
# nothing even in principle, since the original computed `coverage` and then
# neither printed nor tested it. Moved here by review J15/R1 rather than
# revived in place. Nothing it covered is lost, because it covered nothing;
# what remains genuinely untested in the suite is interval coverage itself, for
# either route and either backend.
#
# UNVERIFIED. `TRUEPARS` below is copied verbatim from the original, where it
# was hand-entered and, since the script never ran, never checked against the
# row order `summary()` actually produces. Check it against the printed row
# names before believing any coverage number this prints.
#
# Run from the package root (Rscript dev/simstudies/<file>), or set CTSEM_TREE
# to the package directory. Not part of the package build or its tests.
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(Sys.getenv("CTSEM_TREE", "."), compile = FALSE, quiet = TRUE))

NREP <- as.integer(Sys.getenv("NREP", "100"))
NSUBJECTS <- as.integer(Sys.getenv("NSUBJECTS", "200"))
CORES <- as.integer(Sys.getenv("CORES", "1"))

# Hand-entered in the original; see UNVERIFIED above.
TRUEPARS <- c(10, -1, 2, .5, 5, 1, 1.116, .45)

onerep <- function(i) {
  message("replication ", i)
  t0m <- 10 + rnorm(NSUBJECTS)
  cint <- t0m / 2 + rnorm(NSUBJECTS)
  d <- do.call(rbind, lapply(seq_len(NSUBJECTS), function(subi) {
    gm <- suppressMessages(ctModel(Tpoints = 10, LAMBDA = matrix(1),
      DRIFT = -1, T0MEANS = t0m[subi], DIFFUSION = .5, MANIFESTVAR = 0.5,
      T0VAR = 0, MANIFESTMEANS = 0, CINT = cint[subi]))
    dd <- suppressMessages(data.frame(ctGenerate(ctmodelobj = gm,
      n.subjects = 1, burnin = 0, dtmean = 1, logdtsd = 0)))
    dd$id <- subi
    dd
  }))

  m <- ctModel(type = 'ct', LAMBDA = matrix(1), CINT = 'cint',
    MANIFESTMEANS = 0)

  fitboot <- ctFit(datalong = d, model = m, cores = CORES, priors = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE,
      finishsamples = 5000))
  fithess <- ctFit(datalong = d, model = m, cores = CORES, priors = TRUE,
    optimcontrol = list(carefulfit = FALSE, stochastic = FALSE,
      finishsamples = 5000, uncertainty = 'hessian'))

  rows <- function(s) rbind(s$popmeans, s$popsd, s$rawpopcorr[, 1:5])
  sb <- rows(summary(fitboot))
  sh <- rows(summary(fithess))
  list(boot = sb, hessian = sh)
}

res <- lapply(seq_len(NREP), function(i) try(onerep(i), silent = TRUE))
res <- res[!vapply(res, function(x) inherits(x, 'try-error'), TRUE)]
cat("completed replications:", length(res), "of", NREP, "\n")
if (!length(res)) stop("no replication completed")

covered <- function(route) {
  m <- res[[1]][[route]]
  stopifnot(nrow(m) == length(TRUEPARS))
  hit <- vapply(res, function(r) {
    x <- r[[route]]
    TRUEPARS > x[, '2.5%'] & TRUEPARS < x[, '97.5%']
  }, logical(nrow(m)))
  rowMeans(hit)
}

cov <- cbind(bootstrap = covered('boot'), hessian = covered('hessian'))
rownames(cov) <- rownames(res[[1]]$boot)
cat("\n=== 95% interval coverage over", length(res), "replications ===\n")
print(round(cov, 3))
cat("SIMDONE\n")
