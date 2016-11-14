# helper function to generate an index matrix, or return unique elements of a matrix
indexMatrix<-function(dimension,symmetrical=FALSE,upper=FALSE,lowerTriangular=FALSE, sep=NULL,starttext=NULL,endtext=NULL,
  unique=FALSE,rowoffset=0,coloffset=0,indices=FALSE,diagonal=TRUE,namesvector=NULL){
  if(is.null(namesvector)) namesvector=1:9999
  if(indices==T) sep<-c(",")
  tempmatrix<-matrix(paste0(starttext,namesvector[1:dimension+rowoffset],sep,rep(namesvector[1:dimension+coloffset],each=dimension),endtext),nrow=dimension,ncol=dimension)
  if(upper==TRUE) tempmatrix<-t(tempmatrix)
  if(symmetrical==TRUE) tempmatrix[col(tempmatrix)>row(tempmatrix)] <-t(tempmatrix)[col(tempmatrix)>row(tempmatrix)]
  if(unique==TRUE && symmetrical==TRUE) tempmatrix<-tempmatrix[lower.tri(tempmatrix,diag=diagonal)]
  if(lowerTriangular==TRUE) tempmatrix[col(tempmatrix) > row(tempmatrix)] <- 0
  if(indices==TRUE){
    tempmatrix<-matrix(c(unlist(strsplit(tempmatrix,","))[seq(1,length(tempmatrix)*2,2)],
      unlist(strsplit(tempmatrix,","))[seq(2,length(tempmatrix)*2,2)]),ncol=2)
  }
  return(tempmatrix)
}

#' ctCollapse
#' Easily collapse an array margin using a specified function.
#' @param inarray Input array of more than one dimension.
#' @param collapsemargin Integers denoting which margins to collapse.
#' @param collapsefunc function to use over the collapsing margin.
#' @param ... additional parameters to pass to collapsefunc.
#' @examples
#' testarray <- array(rnorm(900,2,1),dim=c(100,3,3))
#' ctCollapse(testarray,1,mean)
#' @export
ctCollapse<-function(inarray,collapsemargin,collapsefunc,...){
  indims<-dim(inarray)
  out<-array(plyr::aaply(inarray,(1:length(indims))[-collapsemargin],collapsefunc,...,
    .drop=TRUE),dim=indims[-collapsemargin])
  return(out)
}

rl<-function(x) { #robust logical - wrap checks likely to return NA's in this
  if(is.na(x)) return(FALSE) else return(x)
}



#' ctWideNames
#' sets default column names for wide ctsem datasets. Primarily intended for internal ctsem usage.
#' @param n.manifest number of manifest variables per time point in the data.
#' @param Tpoints Maximum number of discrete time points (waves of data, or measurement occasions) 
#' for an individual in the input data structure.
#' @param n.TDpred number of time dependent predictors in the data structure.
#' @param n.TIpred number of time independent predictors in the data structure.
#' @param manifestNames vector of character strings giving column names of manifest indicator variables
#' @param TDpredNames vector of character strings giving column names of time dependent predictor variables
#' @param TIpredNames vector of character strings giving column names of time independent predictor variables
#' @export

ctWideNames<-function(n.manifest,Tpoints,n.TDpred=0,n.TIpred=0,manifestNames='auto',TDpredNames='auto',TIpredNames='auto'){
  
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  
  if(length(manifestNames) != n.manifest) stop("Length of manifestNames does not equal n.manifest!") 
  
  if(n.TDpred > 0){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
    if(length(TDpredNames) != n.TDpred) stop("Length of TDpredNames does not equal n.TDpred!") 
  }
  
  if(n.TIpred > 0){
    if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
    if(length(TIpredNames) != n.TIpred) stop("Length of TIpredNames does not equal n.TIpred!") 
  }
  
  manifestnames<-paste0(manifestNames,"_T",rep(0:(Tpoints-1),each=n.manifest))
  if(n.TDpred > 0 && Tpoints > 1) {
      TDprednames<-paste0(TDpredNames,"_T",rep(0:(Tpoints-2),each=n.TDpred))
  } else {
      TDprednames<-NULL
  }
  if (Tpoints > 1) {
      intervalnames<-paste0("dT",1:(Tpoints-1))
  } else {
      intervalnames <- NULL
  }
  if(n.TIpred>0) TIprednames <- paste0(TIpredNames) else TIprednames <- NULL
  return(c(manifestnames,TDprednames,intervalnames,TIprednames))
}

# generates more complex sequences than seq
cseq <- function(from, to, by){
  temp<-c()
  for(i in from){
    temp<-c(temp,seq(i,to,by))
  }
  temp<-sort(temp)
  return(temp)
}

get_stan_params <- function(object) {
  stopifnot(is(object, "stanfit"))
  params <- grep("context__.vals_r", fixed = TRUE, value = TRUE,
    x = strsplit(get_cppcode(get_stanmodel(object)), "\n")[[1]])
  params <- sapply(strsplit(params, "\""), FUN = function(x) x[[2]])
  params <- intersect(params, object@model_pars)
  return(params)
}


get_stan_massmat<-function(fit){
  
  spars<-get_stan_params(fit)
  spars2<-c()
  for(pari in spars){
    spars2<-c(spars2,grep(paste0(pari,'['),names(fit@sim$samples[[1]]),fixed=TRUE))
  }
  
  massmat<-list()
  for(chaini in 1:fit@sim$chains){
    temp<-c()
    for(pari in spars2){
      newval<-cov(cbind(fit@sim$samples[[chaini]][[pari]][(fit@sim$warmup - fit@stan_args[[1]]$control$adapt_term_buffer):fit@sim$warmup]))
      names(newval)<-names(fit@sim$samples[[chaini]])[pari]
      temp<-c(temp,newval)
    }
    massmat[[chaini]]<-temp
  }
  return(massmat)
}


