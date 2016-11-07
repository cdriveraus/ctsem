#'ctStanGetPars
#'
#'Returns the continuous time parameter matrices for specified subjects of a fit object from ctStanFit
#'
#'@param ctstanfitobj fit object from \link{\code{ctStanFit}}
#'@param subjects Either 'all', or integers denoting which subjects to perform the calculation over. 
#'When multiple subjects are specified, the returned matrices will be a mean over subjects.
#'@param iter Either character string 'all' which will then use all post-warmup iterations, 
#'or an integer specifying which iteration/s to use.
#'@param calcfunc Function to apply over samples, must return a single value. 
#'By default the mean over all samples is returned, but one might also be interested in
#'the \link{\code{sd}} or \link{\code{quantile}} functions.
#'... additional parameters to pass to calcfunc. For instance, with calcfunc = quantile, 
#'the probs argument is needed to ensure only a single value is returned.
#'@examples
#'#posterior mean over all subjects
#'ctStanGetPars(ctstantestfit)
#'
#'#posterior 95% quantiles for subject 2
#'ctStanGetPars(ctstantestfit, subjects=2, calcfunc=quantile, probs = .95)
#'@export
ctStanGetPars <- function(ctstanfitobj,subjects='all',iter='all',
  calcfunc=mean,...){
  
  if(subjects[1] != 'all' && !is.integer(as.integer(subjects))) stop('
  subjects argument must be either "all" or an integer denoting specific subjects')
  
  if(class(ctstanfitobj)!='stanfit') stop('Not an object of class stanfit')
  
  e<-extract(ctstanfitobj) #first dim of subobjects is iter, 2nd subjects
  niter=dim(e$DRIFT)[1]
  
  if(iter!='all') {
    if(!any(iter %in% 1:niter)) stop('Invalid iteration specified!')
    e=lapply(e, function(x) {
      xdims=dim(x)
      out=array(eval(parse(text=
          paste0('x[iter',if(length(xdims)>1) paste0(rep(',',length(xdims)-1),collapse=''),']')
      )),dim=c(length(iter),xdims[-1]))
      return(out)
    }
    )
  }

    nsubjects <- dim(e$indparams)[2]
  if(is.null(nsubjects)) nsubjects=1
  
  if(subjects[1]=='all') subjects=1:nsubjects
    
    collapsemargin<-c(1,2)
    # if(collapseIterations) collapsemargin=1
    # if(collapseSubjects) collapsemargin=c(collapsemargin,2)
  
  for(matname in c('DRIFT','DIFFUSION','CINT','T0MEANS', 
    'T0VAR','MANIFESTMEANS','MANIFESTVAR','LAMBDA', if(!is.null(e$TDPREDEFFECT)) 'TDPREDEFFECT')){
    
    if(dim(e[[matname]])[2] > 1) subselection <- subjects else subselection <- 1
    
    vector <- FALSE
    
    if(matname %in% c('T0MEANS','CINT', 'MANIFESTMEANS')) vector <- TRUE
    
    if(!vector) assign(matname, 
      array(ctCollapse(inarray = e[[matname]][,subselection,,,drop=FALSE],
        collapsemargin = collapsemargin, 
        collapsefunc = calcfunc,...),dim=dim(e[[matname]])[-c(1,2)])
    )
    
    if(vector) assign(matname,
      array(ctCollapse(inarray = e[[matname]][,subselection,,drop=FALSE],
        collapsemargin = collapsemargin, 
        collapsefunc = calcfunc, 
        ...),dim=dim(e[[matname]])[-c(1,2)])
    )
    
  }
  
  model<-list(DRIFT=DRIFT,T0VAR=T0VAR,DIFFUSION=DIFFUSION,CINT=CINT,T0MEANS=T0MEANS,
    MANIFESTMEANS=MANIFESTMEANS,MANIFESTVAR=MANIFESTVAR, LAMBDA=LAMBDA)
  
  if(!is.null(e$TDPREDEFFECT)) model$TDPREDEFFECT<-TDPREDEFFECT
  
  return(model)
}





