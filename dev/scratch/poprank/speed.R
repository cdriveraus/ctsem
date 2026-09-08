## Speed of the regression form, which has never been measured -- the numbers in
## the note are the auxiliary-factor spelling, which carries m more states.
##
## The mechanism from that earlier measurement was that per-row model
## construction dominates, and the reduced form makes each cell an expression
## over the factor states. Here the basis cells still reference one state each
## and only the v regressed cells become expressions over m states, so the
## prediction is that it wins when m is small relative to k and loses when m is
## large. Three shapes test that.
##
## Interleaved, with full rank measured at both ends of each condition order as
## the contamination check, and minima reported.
Sys.setenv(NOT_CRAN = 'true')
devtools::load_all('.', compile = FALSE, quiet = TRUE)
say <- function(...) { cat('@@', ..., '\n', sep = ''); flush.console() }
say('LOADED')

NSUB <- 100L; NT <- 8L

## nlatent latents; `nmean` of the DRIFT diagonals indvarying, all nlatent
## DIFFUSION diagonals indvarying. So k = nmean + nlatent, m = nmean.
mkmodel <- function(nlatent, nmean) {
  DR <- matrix(0, nlatent, nlatent)
  diag(DR) <- paste0('dr', 1:nlatent, '|-log1p_exp(param)')
  DF <- matrix(0, nlatent, nlatent)
  diag(DF) <- paste0('df', 1:nlatent, '|log1p_exp(param)')
  MV <- matrix(0, nlatent, nlatent)
  diag(MV) <- paste0('mv', 1:nlatent, '|log1p_exp(param)')
  m <- suppressMessages(ctModel(type = 'ct', n.latent = nlatent,
    n.manifest = nlatent, LAMBDA = diag(nlatent),
    manifestNames = paste0('Y', 1:nlatent),
    latentNames = paste0('eta', 1:nlatent), DRIFT = DR, DIFFUSION = DF,
    MANIFESTVAR = MV, CINT = matrix(0, nlatent, 1),
    MANIFESTMEANS = matrix(0, nlatent, 1), T0MEANS = matrix(0, nlatent, 1),
    T0VAR = diag(1, nlatent)))
  m$pars$indvarying <- FALSE
  m$pars$indvarying[m$pars$matrix %in% 'DIFFUSION' & !is.na(m$pars$param)] <- TRUE
  drift <- which(m$pars$matrix %in% 'DRIFT' & !is.na(m$pars$param))
  m$pars$indvarying[drift[seq_len(nmean)]] <- TRUE
  m
}
mkdata <- function(nlatent) {
  set.seed(11)
  o <- vector('list', NSUB)
  for (i in seq_len(NSUB)) {
    a <- -log1p(exp(rnorm(nlatent, 0, .5))); q <- log1p(exp(rnorm(nlatent, 0, .5)))
    eta <- matrix(0, NT, nlatent); eta[1, ] <- rnorm(nlatent)
    e <- exp(a); qd <- sqrt(q^2 / (-2 * a) * (1 - e^2))
    for (t in 2:NT) eta[t, ] <- e * eta[t - 1, ] + rnorm(nlatent, 0, qd)
    o[[i]] <- data.frame(id = i, time = seq_len(NT) - 1,
      eta + matrix(rnorm(NT * nlatent, 0, sqrt(.2)), NT, nlatent))
  }
  d <- do.call(rbind, o); names(d)[3:(2 + nlatent)] <- paste0('Y', 1:nlatent); d
}

prep <- function(dat, model, poprank, label) {
  sp <- try(suppressWarnings(suppressMessages(ctFit(datalong = dat, model = model,
    backend = 'julia', fit = FALSE, intoverpop = 'augmented',
    poprank = poprank, cores = 1L))), silent = TRUE)
  if (inherits(sp, 'try-error')) { say(label, ' prepare failed'); return(NULL) }
  h <- structure(sp, class = c('ctJuliaModel', 'ctFitModel'))
  np <- ctsem:::.ctBackendNpar(sp)
  pt <- as.data.frame(sp$parameter_table)
  at <- rep(0.1, np)
  invisible(try(ctJuliaEvaluate(h, at, gradient = TRUE), silent = TRUE))
  list(label = label, h = h, at = at, np = np,
    augdim = max(pt$row[pt$matrix == 'T0MEANS']))
}

for (shape in list(list(nl = 10L, nmean = 10L), list(nl = 10L, nmean = 2L),
                   list(nl = 10L, nmean = 1L))) {
  nl <- shape$nl; nmean <- shape$nmean
  k <- nl + nmean
  dat <- mkdata(nl); model <- mkmodel(nl, nmean)
  say(''); say('##### nlatent ', nl, '  k ', k, '  m ', nmean, ' #####')
  conds <- list(prep(dat, model, NA, 'full (first)'),
                prep(dat, model, 'auto', 'auto'),
                prep(dat, model, NA, 'full (last)'))
  conds <- conds[!vapply(conds, is.null, logical(1))]
  if (length(conds) < 3) { say('  incomplete, skipping'); next }
  NR <- 25L
  tim <- matrix(NA_real_, NR, length(conds))
  for (round in seq_len(NR)) for (j in seq_along(conds)) {
    tim[round, j] <- system.time(ctJuliaEvaluate(conds[[j]]$h, conds[[j]]$at,
      gradient = TRUE))[['elapsed']]
  }
  for (j in seq_along(conds)) say(sprintf(
    '  %-13s npar %3d  augdim %3d   min %.4f  q25 %.4f  med %.4f  max %.4f',
    conds[[j]]$label, conds[[j]]$np, conds[[j]]$augdim,
    min(tim[, j]), quantile(tim[, j], .25), median(tim[, j]), max(tim[, j])))
  b1 <- min(tim[, 1]); b2 <- min(tim[, 3])
  say(sprintf('  contamination check: %.4f vs %.4f  ratio %.3f', b1, b2, b2 / b1))
  say(sprintf('  auto speedup over full rank: %.2fx', min(b1, b2) / min(tim[, 2])))
}
say('DONE')
