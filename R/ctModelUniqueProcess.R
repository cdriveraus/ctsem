#' Give indicators processes of their own
#'
#' Each manifest variable in \code{items} gains a latent, named \code{u_} and
#' the variable's name, that loads on it alone with a fixed unit loading and
#' has its own persistence and innovations. The variable then carries, beside
#' what it shares through the model's other latents and its measurement error,
#' a part that is specific to it and carried over time: an item whose wording
#' picks up something besides the construct, and picks it up persistently.
#'
#' Without it, that part can only go to the measurement error, which is
#' independent over time, so the item's prediction errors stay correlated
#' with its own at nearby times -- what the \code{residuals} family of
#' \code{\link{ctFitAdequacy}} shows for that item alone -- or to the shared
#' latents, whose dynamics it then distorts. The specific process and the
#' measurement error are told apart by persistence: the one carries over
#' time, the other does not.
#'
#' The process has no intercept (the item's own carries its level), a free
#' initial variance, ctsem's default drift transform, so it is mean reverting,
#' and innovations independent of every other latent's.
#'
#' @param ctm A model from \code{ctModel(type = 'ct')}.
#' @param items Which manifest variables, by name or position.
#' @return The extended model.
#' @seealso \code{\link{ctModelImprove}}, which adds one for the items the
#'   adequacy check points to.
#' @export
#' @examples
#' m <- ctModel(type = 'ct', manifestNames = c('Y1', 'Y2', 'Y3'),
#'   latentNames = 'eta', LAMBDA = matrix(c(1, 'l2', 'l3'), 3))
#' m2 <- ctModelUniqueProcess(m, 'Y2')
#' m2$latentNames
#' p <- m2$pars
#' p[p$col %in% 2 & p$matrix %in% c('LAMBDA', 'DRIFT', 'DIFFUSION'),
#'   c('matrix', 'row', 'col', 'param', 'value')]
ctModelUniqueProcess <- function(ctm, items) {
  if (!inherits(ctm, 'ctStanModel')) stop("ctModelUniqueProcess() takes a model ",
    "from ctModel(type = 'ct').", call. = FALSE)
  if (!isTRUE(ctm$continuoustime)) stop("Processes of their own are for continuous ",
    "time models.", call. = FALSE)
  mn <- ctm$manifestNames
  j <- if (is.character(items)) match(items, mn) else as.integer(items)
  if (!length(j) || anyNA(j) || anyDuplicated(j) || any(j < 1L | j > length(mn)))
    stop("items must be distinct manifest variables of the model.", call. = FALSE)
  newnames <- paste0('u_', mn[j])
  taken <- newnames[newnames %in% ctm$latentNames]
  if (length(taken)) stop("The model already has ", paste(taken, collapse = ', '), ".",
    call. = FALSE)
  nl <- ctm$n.latent
  ed <- .ctParsLatentEditor(ctm, newnames)
  allnames <- ed$names
  for (k in seq_along(j)) {
    n <- nl + k
    ed$fixed('LAMBDA', j[k], n, 1)
    ed$free('DRIFT', n, n, paste0('drift_', allnames[n]), tipred = FALSE)
    ed$free('DIFFUSION', n, n, paste0('diff_', allnames[n]), tipred = FALSE)
    ed$fixed('T0MEANS', n, 1, 0)
    ed$free('T0VAR', n, n, paste0('T0var_', allnames[n]), tipred = FALSE)
  }
  ed$finish()
}
