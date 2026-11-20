#' Give indicators processes of their own
#'
#' Each manifest variable in \code{items} gets a latent, named \code{u_} and
#' the variable's name, that loads only on that variable with a fixed loading
#' of one, and has its own persistence and innovations. The variable then has,
#' besides what it shares with others through the model's latents and its
#' measurement error, a specific part that persists over time. This suits an
#' item whose wording picks up something besides the construct, consistently
#' over time.
#'
#' Without such a process, the specific part can only go into the measurement
#' error, which is independent over time, or into the shared latents, whose
#' dynamics it then distorts. In the first case the item's prediction errors
#' stay correlated over nearby times, which the \code{residuals} family of
#' \code{\link{ctFitAdequacy}} shows. The specific process and the measurement
#' error are distinguished by persistence over time.
#'
#' The process has no intercept (the item's own intercept carries its level),
#' ctsem's default drift transform, so it is mean reverting, and innovations
#' independent of the other latents. It starts from its stationary
#' distribution, with mean zero and the variance its persistence and
#' innovations imply (the \code{T0VAR} diagonal \code{'stationary'}, see
#' \code{\link{ctFit}}).
#'
#' With \code{shared = TRUE} the items share one process, with a loading of one
#' on the first item and free loadings on the others, representing what moves
#' them together beyond the model's other latents. Measurement errors that are
#' correlated within an occasion are often better modelled this way, as a
#' process faster than the occasions, unless the design gives a reason for
#' correlated error, such as a method shared by the items.
#'
#' @param ctm A model from \code{ctModel(type = 'ct')}.
#' @param items Which manifest variables, by name or position.
#' @param shared If TRUE, one process for all of \code{items}, named \code{u_}
#'   and their names joined by \code{_}, instead of one each.
#' @param name Names for the new latents, replacing the defaults: one for a
#'   shared process, otherwise one per item. Useful for adding a second process
#'   on the same items at another time scale.
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
#' ctModelUniqueProcess(m, c('Y2', 'Y3'), shared = TRUE)$latentNames
ctModelUniqueProcess <- function(ctm, items, shared = FALSE, name = NULL) {
  if (!inherits(ctm, 'ctStanModel')) stop("ctModelUniqueProcess() takes a model ",
    "from ctModel(type = 'ct').", call. = FALSE)
  if (!isTRUE(ctm$continuoustime)) stop("Processes of their own are for continuous ",
    "time models.", call. = FALSE)
  mn <- ctm$manifestNames
  j <- if (is.character(items)) match(items, mn) else as.integer(items)
  if (!length(j) || anyNA(j) || anyDuplicated(j) || any(j < 1L | j > length(mn)))
    stop("items must be distinct manifest variables of the model.", call. = FALSE)
  if (isTRUE(shared) && length(j) < 2L) stop("A shared process needs two items or more.",
    call. = FALSE)
  newnames <- if (isTRUE(shared)) paste0('u_', paste(mn[j], collapse = '_')) else
    paste0('u_', mn[j])
  if (!is.null(name)) {
    if (length(name) != length(newnames)) stop("name gives the new latents' names: ",
      length(newnames), " here.", call. = FALSE)
    newnames <- as.character(name)
  }
  taken <- newnames[newnames %in% ctm$latentNames]
  if (length(taken)) stop("The model already has ", paste(taken, collapse = ', '), ".",
    call. = FALSE)
  nl <- ctm$n.latent
  ed <- .ctParsLatentEditor(ctm, newnames)
  allnames <- ed$names
  for (k in seq_along(newnames)) {
    n <- nl + k
    if (isTRUE(shared)) {
      ed$fixed('LAMBDA', j[1], n, 1)
      for (jj in j[-1]) ed$free('LAMBDA', jj, n, paste0('lambda_', mn[jj], '_', allnames[n]),
        tipred = FALSE)
    } else ed$fixed('LAMBDA', j[k], n, 1)
    ed$free('DRIFT', n, n, paste0('drift_', allnames[n]), tipred = FALSE)
    ed$free('DIFFUSION', n, n, paste0('diff_', allnames[n]), tipred = FALSE)
    ed$fixed('T0MEANS', n, 1, 0)
    ed$stationary('T0VAR', n, n)
  }
  ed$finish()
}
