# Random effects on a VARIANCE matrix: ctsem's parameterisation against a
# hand-written state augmentation.
#
# This is the study that used to sit inside `if(F)` in
# tests/testthat/test-tdeffectvariation_covtest.R, as `randomEffectsDIFFUSION`
# and as the disabled half of `randomEffectsMANIFESTVAR`. It never ran, so it
# is here rather than there. What it asks is worth asking, and the answer is
# why it is not a test:
#
#   The two parameterisations agree on the LIKELIHOOD exactly -- measured gap 0
#   at n = 400, 200, 100 and 50, for both matrices. They do not reliably agree
#   on the population covariance of the random effect, because a random effect
#   on a variance through log1p_exp leaves the population sd of that effect
#   close to unidentified. The fit says so itself: it warns that "some
#   direction of this model is close to unidentified", and summary() reports
#   that population sd with a 95% interval of 0 to 33. Two runs on the same
#   data and models, differing only in the random starting values, gave a
#   population covariance agreeing to 3e-6 in one and differing by 0.04 in the
#   errsd variance in the other, with everything else agreeing to 2e-5. An
#   assertion on it would pass or fail on the starting values, so the test file
#   keeps the loglik claim and nothing else.
#
#   Recovery is worse and is the reason to run this rather than to guess: at
#   n = 400 the estimated population sd of the MANIFESTVAR random effect missed
#   the sd of the generating draws by a factor of five.
#
# So: run this when you change how random effects on DIFFUSION or MANIFESTVAR
# are integrated, and read the printed table. It is not a pass/fail check.
#
# Run from the package root (Rscript dev/simstudies/<file>), or set CTSEM_TREE
# to the package directory. Not part of the package build or its tests.
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(Sys.getenv("CTSEM_TREE", "."),
  compile = FALSE, quiet = TRUE))

NSUB <- as.numeric(Sys.getenv("NSUB", "400"))
NTIMES <- 50
CORES <- 2

gendata <- function(matrixname, nsubjects, ntimes, seed = 1) {
  set.seed(seed)
  baseline <- rnorm(nsubjects, 2, 2)
  t0m <- rnorm(nsubjects, baseline / 2, 1)
  raweffect <- rnorm(nsubjects, if (matrixname == 'DIFFUSION') -baseline / 3 else
    -baseline / 5, if (matrixname == 'DIFFUSION') .1 else .3)
  effect <- log1p(exp(raweffect))
  for (i in 1:nsubjects) {
    args <- list(silent = TRUE, Tpoints = ntimes, LAMBDA = matrix(1), DRIFT = -1,
      T0MEANS = c(t0m[i]), T0VAR = c(0), CINT = baseline[i], MANIFESTMEANS = 0,
      DIFFUSION = 0.5, MANIFESTVAR = 0.5)
    args[[matrixname]] <- effect[i]
    gm <- suppressMessages(do.call(ctModel, args))
    d <- suppressMessages(data.frame(ctGenerate(ctmodelobj = gm, n.subjects = 1,
      burnin = 0, dtmean = .1, logdtsd = 0)))
    d$id <- i
    if (i == 1) dat <- d else dat <- rbind(dat, d)
  }
  list(dat = dat, t0m = t0m, raweffect = raweffect, effect = effect,
    baseline = baseline)
}

# ctsem's own parameterisation, and the hand-written augmented-state twin. The
# augmented model carries the random effect as latent state 2, with the
# population covariance in T0VAR and the transform written out by hand.
models <- function(matrixname) {
  parname <- if (matrixname == 'DIFFUSION') 'diffusion' else 'errsd'
  margs <- list(silent = TRUE, type = 'ct', T0MEANS = 't0m|param',
    MANIFESTMEANS = 0, CINT = 'cint|param', LAMBDA = matrix(1),
    DIFFUSION = .5, MANIFESTVAR = .5)
  margs[[matrixname]] <- paste0(parname, '|log1p_exp(param)|TRUE')
  m <- do.call(ctModel, margs)

  m2args <- list(silent = TRUE, type = 'ct', Tpoints = 3,
    LAMBDA = matrix(c(1, 0, 0), ncol = 3),
    DRIFT = c('drift', 0, 0, 0, -1e-12, 0, 0, 0, -1e-12),
    DIFFUSION = c(.5, 0, 0, 0, 0, 0, 0, 0, 0), MANIFESTVAR = .5,
    T0MEANS = c('t0m|param', paste0(parname, '|param'), 'cint|param'),
    T0VAR = matrix(c('t0var11 | log1p_exp(2*param-1)', 0, 0,
      't0var21', 't0var22 | log1p_exp(2*param-1)', 0,
      't0var31', 't0var32', 't0var33 | log1p_exp(2*param-1)'), 3, 3, byrow = TRUE),
    CINT = c('state[3]', 0, 0), MANIFESTMEANS = 0)
  m2args[[matrixname]] <- if (matrixname == 'DIFFUSION')
    c('log1p_exp(state[2])', 0, 0, 0, 0, 0, 0, 0, 0) else 'log1p_exp(state[2])'
  m2 <- do.call(ctModel, m2args)
  m2$pars$indvarying <- FALSE
  list(m = m, m2 = m2)
}

for (matrixname in c('DIFFUSION', 'MANIFESTVAR')) {
  cat("\n================ ", matrixname, ", n = ", NSUB, " ================\n", sep = "")
  g <- gendata(matrixname, NSUB, NTIMES)
  mm <- models(matrixname)
  f <- ctFit(datalong = g$dat, model = mm$m, cores = CORES)
  s <- summary(f)
  f2 <- ctFit(datalong = g$dat, model = mm$m2, cores = CORES)
  s2 <- summary(f2)

  a <- f$stanfit$transformedparsfull$rawpopcov[1, , ]
  b <- f2$stanfit$transformedparsfull$pop_T0cov[1, , ]
  cat("\nloglik              :", s$loglik, s2$loglik,
    " gap:", abs(s$loglik - s2$loglik), "\n")
  cat("max |cov gap|       :", max(abs(a - b)), "\n")
  cat("max |cov gap| off the random-effect variance:",
    max(abs((a - b)[-5])), "\n")
  cat("max |corr gap|      :", max(abs(cov2cor(a) - cov2cor(b))), "\n")
  cat("\npopulation sd, estimated against generating (raw scale):\n")
  print(data.frame(
    par = c('t0m', 'randomeffect', 'cint'),
    generating = c(sd(g$t0m), sd(g$raweffect), sd(g$baseline)),
    ctsem = c(f$stanfit$transformedparsfull$rawpopsd),
    augmented = sqrt(diag(b))), row.names = FALSE, digits = 4)
  cat("\npopulation correlation, estimated against generating (raw scale):\n")
  print(data.frame(
    generating = cor(cbind(g$t0m, g$raweffect, g$baseline))[lower.tri(diag(3))],
    ctsem = cov2cor(a)[lower.tri(diag(3))],
    augmented = cov2cor(b)[lower.tri(diag(3))]), row.names = FALSE, digits = 4)
}
cat("\nDONE\n")
