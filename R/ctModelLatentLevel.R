#' Put subjects' levels in the latent processes rather than in the indicators
#'
#' In the usual model a subject's level sits in random manifest intercepts
#' (\code{MANIFESTMEANS}), so the latent state is the subject's deviation from
#' it. That is the same model, for linear dynamics, as one whose latent
#' process carries the level itself through a random latent intercept
#' (\code{CINT}) -- but not once the dynamics depend on the state: a term that
#' reads the state reads the deviation in the first, the level in the second.
#' Persistence that differs near the bottom of a scale, or innovations that
#' grow with how high a subject usually is, are dependence on the level, which
#' state dependence written into the first model cannot express and its tests
#' (\code{\link{ctScoreTest}}) cannot see; \code{ctScoreTest(tests = 'level')}
#' asks about it between subjects.
#'
#' For each latent in \code{latents} with a marker -- an indicator whose loading
#' on it is fixed at 1 and on every other latent at 0 -- whose
#' \code{MANIFESTMEANS} intercept varies over subjects, the intercept's
#' variation moves into the latent: \code{CINT} becomes a free parameter,
#' \code{cint_} and the latent's name, with the intercept's random effect and
#' predictor effects; the marker's intercept is fixed at 0; and the latent's
#' initial mean is its stationary one (\code{T0MEANS} \code{'stationary'}), so
#' a subject starts at their own level. The other indicators' intercepts are
#' left as they are, fixed or varying. A latent without such a marker, or
#' whose \code{CINT} is already free, is left as it is, with a message.
#'
#' The level then rides on the dynamics: \code{CINT} is the level times the
#' rate of return to it, so a level's random effect is no longer one an
#' indicator's filter can integrate on its own, and a fit on the Laplace route
#' carries it as a coordinate.
#'
#' @param ctm A model from \code{ctModel(type = 'ct')}.
#' @param latents Which latents, by name or position; by default every one.
#' @return The model, with the latents whose level moved in attribute
#'   \code{moved}.
#' @seealso \code{\link{ctScoreTest}} for the tests of level dependence.
#' @export
#' @examples
#' m <- ctModel(type = 'ct', manifestNames = c('Y1', 'Y2'),
#'   latentNames = c('eta1', 'eta2'), LAMBDA = diag(2), CINT = matrix(0, 2),
#'   MANIFESTMEANS = matrix(c('mm1', 'mm2'), 2))
#' lm <- ctModelLatentLevel(m)
#' attr(lm, 'moved')
#' p <- lm$pars
#' p[p$matrix %in% c('CINT', 'MANIFESTMEANS', 'T0MEANS'), c('matrix', 'row', 'param',
#'   'value', 'indvarying')]
ctModelLatentLevel <- function(ctm, latents = NULL) {
  if (!inherits(ctm, 'ctStanModel')) stop("ctModelLatentLevel() takes a model ",
    "from ctModel(type = 'ct').", call. = FALSE)
  if (!isTRUE(ctm$continuoustime)) stop("A latent level is for continuous time models.",
    call. = FALSE)
  nl <- ctm$n.latent
  ln <- ctm$latentNames
  idx <- if (is.null(latents)) seq_len(nl) else if (is.character(latents))
    match(latents, ln) else as.integer(latents)
  if (anyNA(idx) || any(idx < 1L | idx > nl)) stop("latents must be latents of the model.",
    call. = FALSE)
  p <- ctm$pars
  varies <- as.logical(.ctResolveLevels(ctm)$pars$indvarying) %in% TRUE
  value <- suppressWarnings(as.numeric(p$value))
  fixedat <- function(k, v) length(k) == 1L && is.na(p$param[k]) && isTRUE(value[k] == v)
  plain <- function(x) !is.na(x) & grepl('^[A-Za-z.][A-Za-z0-9._]*$', x)
  effects <- grep('_effect$', names(p), value = TRUE)
  levelflags <- grep('^indvarying_', names(p), value = TRUE)
  cell <- function(m, r, c) which(p$matrix == m & p$row == r & p$col == c)
  moved <- character()
  for (i in idx) {
    marker <- NA_integer_
    for (j in seq_len(ctm$n.manifest)) {
      if (!fixedat(cell('LAMBDA', j, i), 1)) next
      others <- setdiff(seq_len(nl), i)
      if (all(vapply(others, function(o) fixedat(cell('LAMBDA', j, o), 0), logical(1)))) {
        marker <- j
        break
      }
    }
    km <- if (is.na(marker)) integer() else cell('MANIFESTMEANS', marker, 1L)
    kc <- cell('CINT', i, 1L)
    k0 <- cell('T0MEANS', i, 1L)
    why <- if (is.na(marker)) 'no indicator loads on it alone at 1' else
      if (length(km) != 1L || !plain(p$param[km]) || !varies[km])
        paste0("its marker's intercept (", ctm$manifestNames[marker], ") does not vary") else
      if (length(kc) != 1L || !fixedat(kc, 0)) 'its CINT is not fixed at 0' else
      if (length(k0) != 1L) 'it has no T0MEANS cell' else NA_character_
    if (!is.na(why)) {
      message(ln[i], " left as it is: ", why, ".")
      next
    }
    # The intercept's settings -- random effect, scale, predictor effects --
    # go to the latent's CINT under a name of its own.
    rec <- p[km, , drop = FALSE]
    rec$matrix <- 'CINT'; rec$row <- i; rec$col <- 1L
    rec$param <- paste0('cint_', ln[i]); rec$value <- NA_real_; rec$transform <- 'param'
    p[kc, ] <- rec[, names(p)]
    clear <- function(k, param, v) {
      p$param[k] <<- param; p$value[k] <<- v; p$transform[k] <<- NA_character_
      p$indvarying[k] <<- FALSE
      if ('sdscale' %in% names(p)) p$sdscale[k] <<- NA_real_
      for (e in effects) p[k, e] <<- 'FALSE'
      for (l in levelflags) p[k, l] <<- FALSE
    }
    clear(km, NA_character_, 0)
    clear(k0, 'stationary', NA_real_)
    moved <- c(moved, ln[i])
  }
  ctm$pars <- p
  attr(ctm, 'moved') <- moved
  ctm
}
