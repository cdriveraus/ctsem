#' Extract samples from a ctStanFit object
#'
#' @param object ctStanFit object, samples may be from Stan's HMC, or the importance sampling approach of ctsem.
#' @param subjectMatrices Calculate subject specific system matrices?
#' @param cores Only used if subjectMatrices = TRUE . For faster computation use more cores.
#' @param nsamples either 'all' or an integer denoting number of random samples to extract.
#' @param subjects either 'all', or an integer vector denoting subjects to extract.
#' @param state The latent state at which state dependent cells of the
#' \code{pop_*} matrices are evaluated: \code{'T0MEANS'} (the default),
#' \code{'mean'}, \code{'asymptotic'}, or a numeric state vector. Only
#' \code{backend='julia'} can materialise the matrices anywhere other than the
#' default, so anything else is an error on a stan fit rather than being
#' silently ignored. Irrelevant to a linear model, whose matrices are the same
#' everywhere. Note that it moves the \code{pop_*} arrays only: the subject
#' matrices \code{subjectMatrices=TRUE} adds come from the filter, and are
#' always at each subject's own last observed row.
#' @param ... not used.
#' @return A list of posterior sample arrays. The point the \code{pop_*}
#' matrices were evaluated at is recorded in the \code{'evaluatedAt'} attribute;
#' nothing is messaged, because the summaries call this repeatedly and say it
#' once themselves. \code{\link{ctModelIsNonlinear}} answers whether the point
#' matters for this model at all.
#' @aliases extract
#' @examples
#' \donttest{
#' e = ctExtract(ctstantestfit)
#' attr(e, 'evaluatedAt')
#' }
#' @export
ctExtract <- function(object, subjectMatrices=FALSE, cores=1, nsamples='all', subjects='all',
  state=NULL, ...) UseMethod("ctExtract")

#' @export
ctExtract.ctStanFit <- function(object,subjectMatrices=FALSE,cores=1,nsamples='all', subjects='all',
  state=NULL, ...){
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
