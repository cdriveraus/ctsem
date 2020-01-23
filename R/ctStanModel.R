ctModelUnlist<-function(ctmodelobj,matnames=c('T0MEANS','LAMBDA','DRIFT','DIFFUSION','MANIFESTVAR','MANIFESTMEANS', 'CINT', 'TDPREDEFFECT', 'T0VAR','PARS')){
  out<-matrix(NA,nrow=0,ncol=5)
  out<-as.data.frame(out,stringsAsFactors =FALSE)
  for(obji in matnames){
    if(!is.null(dim(ctmodelobj[[obji]])) && !is.null(ctmodelobj[[obji]]) && !is.na(ctmodelobj[[obji]][1])){
      for(rowi in 1:nrow(ctmodelobj[[obji]])){
        for(coli in 1:ncol(ctmodelobj[[obji]])){
          out<-rbind(out,data.frame(obji,rowi,coli,
            ifelse(is.na(suppressWarnings(as.numeric(ctmodelobj[[obji]][rowi,coli]))), #ifelse element is character string
              ctmodelobj[[obji]][rowi,coli],
              NA),
            ifelse(!is.na(suppressWarnings(as.numeric(ctmodelobj[[obji]][rowi,coli]))), #ifelse element is numeric
              as.numeric(ctmodelobj[[obji]][rowi,coli]),
              NA),
            stringsAsFactors =FALSE
          ))
        }
      }
    }
  }
  colnames(out)<-c('matrix','row','col','param','value')
  return(out)
}

#' Convert a frequentist (omx) ctsem model specification to Bayesian (Stan).
#'
#' @param ctmodelobj ctsem model object of type 'omx' (default)
#' @param type either 'stanct' for continuous time, or 'standt' for discrete time.
#' @param tipredDefault Logical. TRUE sets any parameters with unspecified time independent 
#' predictor effects to have effects estimated, FALSE fixes the effect to zero unless individually specified.
#'
#' @return List object of class ctStanModel, with random effects specified for any intercept type parameters
#' (T0MEANS, MANIFESTMEANS, and or CINT), and time independent predictor effects for all parameters. Adjust these
#' after initial specification by directly editing the \code{pars} subobject, so \code{model$pars} . 
#' @export
#'
#' @examples
#' model <- ctModel(type='omx', Tpoints=50,
#' n.latent=2, n.manifest=1, 
#' manifestNames='sunspots', 
#' latentNames=c('ss_level', 'ss_velocity'),
#' LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
#' DRIFT=matrix(c(0, 1,   'a21', 'a22'), nrow=2, ncol=2, byrow=TRUE),
#' MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
#' # MANIFESTVAR=matrix(0, nrow=1, ncol=1),
#' CINT=matrix(c(0, 0), nrow=2, ncol=1),
#' DIFFUSION=matrix(c(
#'   0, 0,
#'   0, "diffusion"), ncol=2, nrow=2, byrow=TRUE))
#' 
#' stanmodel=ctStanModel(model)
#' 
#' 
ctStanModel<-function(ctmodelobj, type='stanct',tipredDefault=TRUE){
  if(type=='stanct') continuoustime<-TRUE
  if(type=='standt') continuoustime<-FALSE
  
  ctm <- ctmodelobj
  
  if(!is.null(ctm$timeVarying)) stop('Time varying parameters not allowed for ctsem Stan model at present! Correct ctModel spec')
  
  
  
  
  #read in ctmodel values
  n.latent<-ctm$n.latent
  n.manifest<-ctm$n.manifest
  Tpoints<-ctm$Tpoints
  n.TDpred<-ctm$n.TDpred
  n.TIpred<-ctm$n.TIpred
  
  manifestNames<-ctm$manifestNames
  latentNames<-ctm$latentNames
  TDpredNames<-ctm$TDpredNames
  TIpredNames<-ctm$TIpredNames
  
  ctspec<-ctModelUnlist(ctm)
  
  freeparams<-is.na(ctspec[,'value'])
  
  ctspec$transform<- NA
  ctspec$multiplier <- 1
  ctspec$meanscale <- 1
  ctspec$transform <- 0
  ctspec$offset <- 0
  ctspec$inneroffset <- 0
  
  
  ######### STAN parameter transforms
  for(pi in 1:length(ctspec$matrix)){
    if(freeparams[pi]){
      if(ctspec$matrix[pi] %in% c('T0MEANS','MANIFESTMEANS','TDPREDEFFECT','CINT')) {
        ctspec$meanscale[pi] <-10
      }
      if(ctspec$matrix[pi] %in% c('LAMBDA')) {
        ctspec$offset[pi] <- 0.5
        ctspec$meanscale[pi] <- 5
      }
      
      if(ctspec$matrix[pi] %in% c('DIFFUSION','MANIFESTVAR', 'T0VAR')) {
        if(ctspec$row[pi] != ctspec$col[pi]){
          ctspec$transform[pi] <- 0
          # if(ctspec$matrix[pi] %in% c('DIFFUSION')) ctspec$meanscale[pi] <-2
        }
        if(ctspec$row[pi] == ctspec$col[pi]){
          ctspec$transform[pi] <- 1
          ctspec$meanscale[pi] <- 2
          ctspec$multiplier[pi] <- 5
          ctspec$offset[pi] <- 1e-6
          if(ctspec$matrix[pi] %in% c('DIFFUSION')) ctspec$multiplier[pi] <-10
        }
      }
      if(ctspec$matrix[pi] %in% c('DRIFT')) {
        if(ctspec$row[pi] == ctspec$col[pi]){
          if(continuoustime==TRUE) {
            ctspec$transform[pi] <- 1
            ctspec$meanscale[pi] <- -2
            ctspec$multiplier[pi] <- -2
            ctspec$offset[pi] <- -1e-6
          }
          if(continuoustime==FALSE) {
            ctspec$transform[pi] <- 0
            ctspec$meanscale[pi] <- 1
            ctspec$offset[pi] <- .5
          }
        }
        if(ctspec$row[pi] != ctspec$col[pi]){
          ctspec$transform[pi] <- 0
          ctspec$meanscale[pi] <- 1
        }
      }
    }
  }
  
  ctspec[!is.na(ctspec$value),c('transform','multiplier','meanscale','offset','inneroffset')] <- NA
  
  
  for(ri in 1:nrow(ctspec)){ #convert back to text for new approach
    if(!is.na(as.integer(ctspec$transform[ri]))){
      ctspec$transform[ri] <- Simplify(tform(param = 'param',
        transform = as.integer(ctspec$transform[ri]),
        multiplier = ctspec$multiplier[ri],
        meanscale = ctspec$meanscale[ri],
        offset = ctspec$offset[ri],
        inneroffset =ctspec$inneroffset[ri],
        singletext = TRUE))
    }
  }
  ctspec$multiplier <- NULL
  ctspec$meanscale <- NULL
  ctspec$offset <- NULL
  ctspec$inneroffset <- NULL
  
  
  
  nparams<-sum(freeparams)
  
  # indvarying <- 'all'
  # if(all(indvarying=='all'))  {
  #   indvarying<-rep(TRUE,nparams)
  #   # indvarying[ctspec$matrix[is.na(ctspec$value)] %in% 'T0VAR'] <- FALSE
  # }
  # 
  # if(length(indvarying) != nparams) stop('indvarying must be ', nparams,' long!')
  # nindvarying <- sum(indvarying)
  # 
  
  
  ctspec$indvarying<-FALSE
  ctspec$indvarying[!is.na(ctspec$transform) & ctspec$matrix %in% c('T0MEANS','MANIFESTMEANS','CINT')] <- TRUE
  
  ctspec$sdscale<-NA
  ctspec$sdscale[is.na(ctspec$value)]<-1
  
  
  
  
  if(n.TIpred > 0) {
    tipredspec<-matrix(TRUE,ncol=n.TIpred,nrow=1)
    colnames(tipredspec)<-paste0(TIpredNames,'_effect')
    ctspec<-cbind(ctspec,tipredspec,stringsAsFactors=FALSE)
    ctspec[,paste0(TIpredNames,'_effect')]<-tipredDefault
    for(predi in TIpredNames){
      class(ctspec[,paste0(predi,'_effect')])<-'logical'
    }
  }
  
  
  
  getwords <- function(x) gsub('\\W+',',',x)
  
  
  for(pi in 1:nrow(ctspec)){ #check for complex / split specifications
    if(grepl('|',ctspec$param[pi],fixed=TRUE)){
      ctspec$param[pi] <- gsub('$',' ',ctspec$param[pi])
      split = strsplit(ctspec$param[pi],split = '|',fixed=TRUE)[[1]]
      split=sapply(split,function(x) gsub(' ','',x))
      
      tisplit <- NA
      if(length(split) > 4){ #check for ti pred spec in splits
        if(n.TIpred < 1 || length(split) > 5) stop(paste0('Param spec has too many separators!  ', ctspec$param[pi]))
        tisplit <- split[5]
        split <- split[1:4]
      }
      
      if(grepl('\\W',split[1]) || any(sapply(latentNames,function(x){ #if symbols or latent states
        grepl(paste0('\\b(',x,')\\b'),split[1])
      }))) {
        if(!simpleStateCheck(split[1])) stop(paste0(split[1],' invalid -- Matrix elements involving multiple parameters / latent states cannot have | separators -- transformations should be specified as part of the first element, indvarying and tipredeffects must be specified in the corresponding singular PARS matrix elements.'))
      }
      
      nonzero <- which(!split %in% '')
      ctspec[pi,c('param','transform','indvarying','sdscale')] <- 
        list(NA,ctspec$transform[pi],ctspec$indvarying[pi],1) #base values
      ctspec[pi,c('param','transform','indvarying','sdscale')[1:(length(split))]][nonzero] <- 
        split[1:(length(split))][nonzero] #update with non zero length split elements
      
      whichtipreds <- c()
      timessage <- c()
      if(!is.na(tisplit)){
        tisplit <- strsplit(getwords(tisplit),split = ',')[[1]]
        ctspec[pi,paste0(TIpredNames,'_effect')] <- FALSE #first set all FALSE
        if(!tisplit[1] %in% ''){
          for(ti in TIpredNames){ #check which tipreds were included
            for(spliti in tisplit){
              if(!spliti %in% TIpredNames) stop (spliti,' is not a time independent predictor!')
              if(grepl(paste0('\\b(',ti,')\\b'),spliti)) whichtipreds <- c(whichtipreds,ti)
            }
          }
        }
        
        if(!is.null(whichtipreds))  ctspec[pi,paste0(whichtipreds,'_effect')] <- TRUE #set those effects TRUE
        if(is.null(whichtipreds)) whichtipreds <- 'NULL'
        timessage <- paste0(ctspec$param[pi],' tipred effects from: ', paste0(whichtipreds,collapse=', '))
        
      }
      
      
      message('Custom par ',ctspec$param[pi],' set as: ',paste0(
        colnames(ctspec[pi,c('param','transform','indvarying','sdscale')]),' = ',
        ctspec[pi,c('param','transform','indvarying','sdscale')],'; '),timessage)
    }
  }
  
  if(n.TIpred > 0 && sum(unlist(ctspec[,paste0(TIpredNames,'_effect')]))==0) warning('TI predictors included but no effects specified!')
  
  for(ri in 1:nrow(ctspec)){ #set NA's on complex params
    
    if(grepl('\\W',gsub('.','',ctspec$param[ri],fixed=TRUE)) || ctspec$param[ri] %in% latentNames){
      ctspec$value[ri] <- NA
      if(!simpleStateCheck(ctspec$param[ri])) ctspec$transform[ri] <-NA #allow transforms for simplestates
    }
  }
  
  ctspec$indvarying <- as.logical(ctspec$indvarying)
  
  out<-list(pars=ctspec,n.latent=n.latent,n.manifest=n.manifest,n.TIpred=n.TIpred,n.TDpred=n.TDpred,
    latentNames=latentNames,manifestNames=manifestNames,
    TIpredNames=TIpredNames,TDpredNames=TDpredNames,
    subjectIDname=ctm$id,
    timeName=ctm$time,
    continuoustime=continuoustime)
  class(out)<-'ctStanModel'
  
  out$tipredeffectscale <- 1
  out$tipredsimputedscale <- 1
  
  # out$popsdpriorscale <- 1
  out$rawpopsdbase <- 'normal(0,1)' #'cauchy(0,1)'
  out$rawpopsdbaselowerbound <- NA
  out$rawpopsdtransform <- 'log1p(exp(2*rawpopsdbase-1)) .* sdscale' #'log(1+exp(2*rawpopsdbase)) .* sdscale' #'exp(rawpopsdbase * 2 -2) .* sdscale' # 'rawpopsdbase .* sdscale' #
  # out$stationarymeanprior <- NA
  # out$stationaryvarprior <- NA
  out$manifesttype <- rep(0,n.manifest)
  out$covmattransform <- 'unconstrainedcorr'
  
  return(out)
}

