#' Half-lives of the latent processes
#'
#' The lag at which a process's model-implied autocorrelation,
#' \eqn{Cor(x_k(t + u), x_k(t))} in its stationary distribution, falls to one
#' half: how long the process takes to forget half of where it is. Computed for
#' every draw of the population parameters, or of each subject's own, so the
#' spread over draws is the uncertainty of the half-life.
#'
#' @param fit A fit from \code{\link{ctFit}}.
#' @param processes Names of the latent processes, all by default.
#' @param subjects \code{'popmean'} (the default) for the population
#'   parameters, or subjects as ids in the data (see \code{realid}). On a
#'   \code{backend='julia'} fit with random effects, each subject's draws carry
#'   their own uncertainty as well as the population's (plausible values; see
#'   \code{\link{ctExtract}}).
#' @param source Optional latent names: the half-life of the part of each
#'   process these latents' innovations produce, the autocorrelation taken in
#'   the stationary distribution they alone would give. With a slow process
#'   beside a fast one, the whole autocorrelation of the fast one may never
#'   halve, and its own component's half-life is the persistence of interest.
#'   The other latents' innovations must be independent of these.
#' @param nsamples Number of draws to use.
#' @param realid If TRUE (the default), subjects are read as the ids in the data.
#' @return An array of half-lives by draw, subject and process, in the units of
#'   the time column (continuous time) or of measurement occasions (discrete
#'   time, interpolated between lags). \code{Inf} where the autocorrelation
#'   never falls to a half within a thousand of the system's slowest time
#'   constants, \code{NA} where the system has no stationary distribution or
#'   the source does not reach the process.
#' @details For a model whose DRIFT depends on the state, the half-life is that
#'   of the model linearised where the matrices were evaluated, which is said.
#' @examples
#' \donttest{
#' hl <- ctHalfLife(ctstantestfit)
#' apply(hl, 3, quantile, probs = c(.025, .5, .975))
#' }
#' @export
ctHalfLife <- function(fit, processes = NULL, subjects = 'popmean', source = NULL,
  nsamples = 200, realid = TRUE) {
  model <- .ctFitModelObject(fit)
  nlatent <- .ctFitNlatent(fit)
  latents <- model$latentNames[seq_len(nlatent)]
  if (is.null(processes)) processes <- latents
  which <- match(processes, latents)
  if (anyNA(which)) stop("processes names latents this model does not have: ",
    paste(processes[is.na(which)], collapse = ", "), ".", call. = FALSE)
  sourceindex <- if (!is.null(source)) {
    s <- match(source, latents)
    if (anyNA(s)) stop("source names latents this model does not have: ",
      paste(source[is.na(s)], collapse = ", "), ".", call. = FALSE)
    s
  }
  population <- identical(subjects, 'popmean')
  labels <- 'popmean'
  if (!population) {
    resolved <- .ctResolveSubjects(fit, subjects, realid)
    ordering <- order(resolved)
    labels <- attr(resolved, 'ids')[ordering]
    subjects <- as.vector(resolved)[ordering]
  }
  e <- ctExtract(fit, subjectMatrices = !population,
    nsamples = min(nsamples, .ctFitNsamples(fit)),
    subjects = if (population) 'all' else subjects, realid = FALSE)
  pick <- function(name) {
    if (population) {
      a <- e[[paste0('pop_', name)]]
      return(array(a, c(dim(a)[1L], 1L, dim(a)[-1L])))
    }
    a <- e[[paste0('subj_', name)]]
    if (is.null(a)) {
      # Nothing varies by subject: each has the population's.
      p <- e[[paste0('pop_', name)]]
      a <- array(NA_real_, c(dim(p)[1L], length(subjects), dim(p)[-1L]))
      for (i in seq_along(subjects)) a[, i, , ] <- p
      return(a)
    }
    # Every subject's, when the extract did not narrow to those asked for.
    if (dim(a)[2L] != length(subjects)) a <- a[, subjects, , , drop = FALSE]
    a
  }
  drift <- pick('DRIFT')
  diffusion <- pick('DIFFUSIONcov')
  ndraws <- dim(drift)[1L]
  nsub <- dim(drift)[2L]
  .ctContextMessage(fit, if (population) .ctContextPopLabel else
    .ctContextSubjectLabel, paste0("The half-life is therefore that of the ",
      "model linearised there."))
  out <- array(NA_real_, c(ndraws, nsub, length(which)), dimnames = list(
    Sample = seq_len(ndraws), Subject = labels, Process = latents[which]))
  for (d in seq_len(ndraws)) for (i in seq_len(nsub)) {
    A <- matrix(drift[d, i, seq_len(nlatent), seq_len(nlatent)], nlatent)
    G <- matrix(diffusion[d, i, seq_len(nlatent), seq_len(nlatent)], nlatent)
    out[d, i, ] <- .ctHalfLifeOne(A, G, which, sourceindex,
      isTRUE(model$continuoustime))
  }
  out
}

# One system's half-lives: the first lag at which `(e^(A u) Q)[k, k] / Q[k, k]`
# reaches a half, `Q` the stationary covariance of the innovations considered.
.ctHalfLifeOne <- function(A, G, which, source, continuoustime) {
  if (!is.null(source)) {
    other <- setdiff(seq_len(nrow(G)), source)
    if (length(other) && any(abs(G[source, other]) > 1e-8 * max(abs(G))))
      stop("The source's innovations are correlated with the other latents', ",
        "so its part of the variance is not separate.", call. = FALSE)
    masked <- matrix(0, nrow(G), ncol(G))
    masked[source, source] <- G[source, source]
    G <- masked
  }
  Q <- .ctNetworkAsym(A, G, continuoustime)
  if (is.null(Q)) return(rep(NA_real_, length(which)))
  values <- eigen(A, only.values = TRUE)$values
  rate <- if (continuoustime) abs(Re(values)) else -log(pmax(Mod(values), 1e-300))
  slow <- 1 / max(min(rate), 1e-12)
  fast <- 1 / max(max(rate), 1e-12)
  lags <- if (continuoustime) c(0, exp(seq(log(fast * 1e-3), log(slow * 1e3),
    length.out = 400L))) else 0:ceiling(min(slow * 1e3, 1e5))
  propagate <- function(u) if (continuoustime) as.matrix(expm::expm(A * u)) else {
    out <- diag(nrow(A))
    for (s in seq_len(u)) out <- out %*% A
    out
  }
  vapply(which, function(k) {
    if (!(Q[k, k] > 0)) return(NA_real_)
    auto <- function(u) (propagate(u) %*% Q)[k, k] / Q[k, k] - 0.5
    f <- vapply(lags, auto, numeric(1))
    cross <- which(f[-1L] <= 0 & f[-length(f)] > 0)
    if (!length(cross)) return(Inf)
    lo <- lags[cross[1L]]
    hi <- lags[cross[1L] + 1L]
    if (continuoustime) return(stats::uniroot(auto, c(lo, hi), tol = 1e-10)$root)
    # Between two occasions, on the log of the autocorrelation.
    a <- log(f[cross[1L]] + 0.5)
    b <- log(f[cross[1L] + 1L] + 0.5)
    lo + (log(0.5) - a) / (b - a)
  }, numeric(1))
}
