#' Convert samples from a stanfit object to the unconstrained scale
#'
#' @param fit stanfit object.
#' @param standata only necessary if R session has been restarted since fitting model -- used to reinitialize 
#' stanfit object.
#'
#' @return Matrix containing columns of unconstrained parameters for each post-warmup iteration.
#' @export
#'
#' @examples
#' \donttest{
#' #get data
#' sunspots<-sunspot.year
#' sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#' id <- 1
#' time <- 1749:1924
#' datalong <- cbind(id, time, sunspots)
#' 
#' #setup model
#' ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
#'   manifestNames='sunspots', 
#'   latentNames=c('ss_level', 'ss_velocity'),
#'   LAMBDA=matrix(c( 1, 'ma1| log(1+(exp(param)))' ), nrow=1, ncol=2),
#'   DRIFT=matrix(c(0, 'a21 | -log(1+exp(param))', 1, 'a22'), nrow=2, ncol=2),
#'   MANIFESTMEANS=matrix(c('m1|param * 10 + 44'), nrow=1, ncol=1),
#'   MANIFESTVAR=diag(0,1), #As per original spec
#'   CINT=matrix(c(0, 0), nrow=2, ncol=1),
#'   DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
#' 
#' #fit
#' ssfit <- ctStanFit(datalong, ssmodel, 
#'   iter=200, chains=2,optimize=FALSE, priors=TRUE,control=list(max_treedepth=4))
#' umat <- stan_unconstrainsamples(ssfit$stanfit$stanfit)
#' }
stan_unconstrainsamples <- function(fit, standata=NA){
  if(!'stanfit' %in% class(fit)) stop('not a stanfit object')
  npars <- try(get_num_upars(fit),silent=TRUE) #$stanmodel)
  
  if(class(npars)[1]=='try-error'){ #in case R has been restarted or similar
    if(any(!is.na(standata))){
      newfit <- stan_reinitsf(fit@stanmodel,standata) 
    } 
    else stop('stanfit object must be reinitialized but no data is provided')
  } else newfit <- fit #no need for reinit
  
  cmat=as.matrix(fit)
  # clist=apply(cmat,1,function(x) relist(flesh = x,skeleton = fit@inits[[1]]))
  
  if(is.null(names(fit@inits[[1]]))) skel = fit@inits else skel=fit@inits[[1]]
  clist=apply(cmat,1,function(x) relistarrays(flesh=x,skeleton=skel))

  ulist=matrix(unlist(lapply(clist,function(x) unconstrain_pars(newfit,x))),ncol=length(clist))
  return(ulist)
}
