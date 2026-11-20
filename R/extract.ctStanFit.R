#' Extract samples from a ctsem fit
#'
#' Returns the posterior (or approximate posterior) draws of a fit's
#' parameters, as arrays named for the model matrices.
#'
#' @param object A fit from \code{\link{ctFit}}, on either backend: sampled, or
#' optimized with its approximate posterior draws.
#' @param subjectMatrices Calculate subject specific system matrices? On a
#' \code{backend='julia'} fit with \code{intoverpop='laplace'}, or a sampled
#' fit, each draw's subject matrices combine the population draw with a draw
#' of the subject's own effects from their posterior given its data. Their
#' spread over draws then shows the subject's full uncertainty, as
#' \code{\link{ctSubjectPars}} reports by default.
#' @param cores Only used if subjectMatrices = TRUE . For faster computation use more cores.
#' @param nsamples either 'all' or an integer denoting number of random samples to extract.
#' @param subjects either 'all', or the subjects to extract, as ids in the data
#' (see \code{realid}).
#' @param realid If TRUE (the default), subjects are read as the ids in the data,
#' and integers that are not ids there are read as ctsem's numbering of the
#' subjects, 1 to N, with a message. If FALSE, they are read as that numbering.
#' @param state The latent state at which state dependent cells of the
#' \code{pop_*} matrices are evaluated: \code{'T0MEANS'} (the default),
#' \code{'mean'}, \code{'asymptotic'}, or a numeric state vector. Values other
#' than the default need \code{backend='julia'} and are an error for a stan
#' fit. No effect for a linear model. Only the \code{pop_*} arrays are
#' affected. The subject matrices from \code{subjectMatrices=TRUE} come from
#' the Kalman filter and are at each subject's last observed row.
#' @param ... not used.
#' @return A list of posterior sample arrays. The point at which the
#' \code{pop_*} matrices were evaluated is recorded in the
#' \code{'evaluatedAt'} attribute, without a message.
#' \code{\link{ctModelIsNonlinear}} shows whether the point matters for the
#' model.
#'
#' Where parameters vary over subjects, the list also holds, per draw,
#' \code{popsd} (draws by parameters), \code{rawpopcov} and \code{rawpopcorr}
#' (the covariance and correlations on the raw scale, draws by parameters by
#' parameters) and \code{popcov} (the covariance on the natural scale). On a
#' \code{backend='julia'} fit they are named by parameter, and a model with
#' several grouping levels gives each level's as \code{rawpopcov.<level>} and
#' so on. There, \code{popcov} combines \code{popsd} with the raw
#' correlations, which is exact where a parameter's transform is affine.
#' @aliases extract
#' @examples
#' \donttest{
#' e = ctExtract(ctstantestfit)
#' attr(e, 'evaluatedAt')
#' }
#' @export
ctExtract <- function(object, subjectMatrices=FALSE, cores=getOption("mc.cores", 2L), nsamples='all', subjects='all',
  state=NULL, realid=TRUE, ...) UseMethod("ctExtract")

#' @export
ctExtract.ctStanFit <- function(object,subjectMatrices=FALSE,cores=getOption("mc.cores", 2L),nsamples='all', subjects='all',
  state=NULL, realid=TRUE, ...){
  # inherits(), not class() %in%: a fit carries both 'ctStanFit' and 'ctFit'
  # since the ctFit rename, and `if` on a length-2 condition is an error in
  # R >= 4.2 -- so this guard used to reject every fit it was given.
  if(!inherits(object, c('ctStanFit', 'stanfit'))) stop('Not a ctStanFit or stanfit object')

  # Stan's generated quantities computed these matrices once, during sampling,
  # at the point the filter was at; there is no way to re-materialise them
  # somewhere else afterwards. `state` used to reach here through `...` and be
  # dropped, while ctExtract.ctJuliaFit honoured it -- so the same call returned
  # a different linearisation on each backend and said nothing on either.
  .ctContextRequireStateSupport(object, state)
  if(!'all' %in% subjects && inherits(object, 'ctStanFit'))
    subjects <- sort(as.vector(.ctResolveSubjects(object, subjects, realid)))

  
  
  if(length(object$stanfit$stanfit@sim)==0){
    samps <- object$stanfit$rawposterior
    if(!nsamples %in% 'all') samps <- samps[sample(1:nrow(samps),nsamples),,drop=FALSE]
    if(subjectMatrices && object$standata$savesubjectmatrices==0){
      if(!'all' %in% subjects) object$standata<- standatact_specificsubjects(standata = object$standata,subjects = subjects)
      out = stan_constrainsamples(sm = object$stanmodel,standata = object$standata,
        samples = samps,
        cores = cores,savescores = FALSE,savesubjectmatrices = subjectMatrices,
        dokalman = TRUE,onlyfirstrow = FALSE)
    } else out <- object$stanfit$transformedpars
  }
  
  if(length(object$stanfit$stanfit@sim)>0){
    if(subjectMatrices & object$standata$savesubjectmatrices!=1){
      samps <- t(stan_unconstrainsamples(object$stanfit$stanfit,standata=object$standata))
      if(!nsamples %in% 'all') samps <- samps[sample(1:nrow(samps),nsamples),,drop=FALSE]
      out = stan_constrainsamples(sm = object$stanmodel,standata = object$standata,
        samples = samps,
        cores = cores,savescores = FALSE,savesubjectmatrices = subjectMatrices)
    } else  out <- rstan::extract(object$stanfit$stanfit)
  } 
  
  out$Ygen[out$Ygen==99999] <- NA
  
  # if(!is.null(out$rawpopc)){
  #   out$rawpopcov <- array(out$rawpopc[,4,,],dim=dim(out$rawpopc)[-2])
  #   out$rawpopcorr <-  array(out$rawpopc[,3,,],dim=dim(out$rawpopc)[-2])
  #   out$rawpopcovchol <-  array(out$rawpopc[,2,,],dim=dim(out$rawpopc)[-2])
  #   out$rawpopcovbase <-  array(out$rawpopc[,1,,],dim=dim(out$rawpopc)[-2])
  # }

  # Record rather than message: ctExtract() is called repeatedly by summary(),
  # ctSummaryMatrices() and the plot helpers, and each of those says it once
  # itself through .ctContextMessage(). The label is free; .ctContextAttach()
  # is not -- it re-derives the cell table each time, which measured ~7 ms
  # locally against ~0 for the rest of this function, and ctSummaryMatrices()
  # attaches the same cells to its own output anyway.
  attr(out, 'evaluatedAt') <- .ctContextPopLabel

  return(out)
}
