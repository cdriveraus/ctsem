#' Adds free parameters to a ctsem model in a step-wise fashion.
#' 
#' Takes a base ctmodelobj, and frees parameters in the specified matrices in a stepwise fashion.  
#' At each step, the optimizer is run and rerun for each free parameter up to a maximum number of tries, 
#' until the likelihood from the more constrained model is beaten.
#' Once new estimates are obtained for each candidate free parameter, that parameter which most improves the AIC 
#' is selected, and the process starts again.
#' @seealso \code{\link{ctModel}},\code{\link{ctRefineTo}}
#' @param datawide Data in ctsem wide format
#' @param ctmodelobj A continuous time model specified via the \code{\link{ctModel}} function.
#' @param matrices A vector of continuous time matrix names to iteratively free
#' @param existingfit ctsemFit object already fit from the ctmodelobj with which to start the process --
#' not necessary, but can be helpful to ensure the base model is optimized well to begin with.
#' @param retryattempts Maximum number of times to retry each fit.
#' @param tolerance Numeric. Defaults to 0. Added to best AIC value to determine acceptable threshold for 
#' accepting a best new parameter.
#' @fastrefine If TRUE, instead of only accepting best parameter before beginning next major iteration, 
#' all parameters with improved fit value are accepted.
#' @param fitWithoutInits if TRUE, runs additional fit every iteration using default ctsem inits, rather
#' than best current inits. May help avoid path dependence, but may be more difficult to optimize...
#' @param ... additional parameters to pass to \code{\link{ctFit}}.
#' @export

ctRefineFrom<-function(datawide,ctmodelobj,matrices,existingfit=NULL,retryattempts=0,
  tolerance=0,fastrefine=FALSE,fitWithoutInits=FALSE,...){
  
  if(fastrefine==TRUE & any(matrices !='DRIFT')) stop('fastrefine only appropriate for refining over DRIFT matrix at present')
  
  #internal plot function
  ctRefinePlot<-function(refined){
    best<-which(refined$testsummary[,'AIC'] == min(refined$testsummary[,'AIC'],na.rm=T))[1]
    params<-names(refined$fits[[best]]$mxobj$output$estimate)
    
    
    par(mfrow=c(3,3))
    for(p in params){
      out<-c()
      for( i in 1:length(refined$fits)){
        out<-c(out,refined$fits[[i]]$mxobj$output$estimate[p])
        
      }
      #     out[is.na(out)]<-0
      #     out[out^2>3*sd(out^2)]<-NA
      bplot<-boxplot.stats(out)
      
      out[out>bplot$stats[5] | out < bplot$stats[1]] <- NA
      plot(out,main=p,ylabel='Step',xlabel='Value')
    }
  }
  
  
  
  testcount<-1 #set number of model tests run so far
  
  if(is.null(existingfit)) basefit <- ctFit(datawide,ctmodelobj,carefulFit=T,confidenceintervals='DRIFT',...) #run base fit
  if(!is.null(existingfit)) basefit <- existingfit
  startValues<-omxGetParameters(basefit$mxobj)
  AIC<-summary(basefit$mxobj)$AIC.Mx #extract AIC
  freedParam<-paste0("none") #list the parameter that we freed
  estParams<-summary(basefit$mxobj)$estimatedParameters #and the number of estimated params
  Minus2LogLikelihood <- summary(basefit$mxobj)$Minus2LogLikelihood #and the likelihood
  fit2beat<-Minus2LogLikelihood #update fit to beat
  message('Base -2LL:', fit2beat)
  fits<-list(basefit)
  models<-list(ctmodelobj)
  
  for(m in 1:length(matrices)){ #for every matrix specified to refine over
    symmetricfree<-FALSE
    
    freematrix<-eval(parse(text=(paste0("ctmodelobj$",matrices[m])))) #extract the matrix from the model spec
    
    if(matrices[m]=="DRIFT") freeindex <- which(diag(nrow(freematrix))!=1,arr.ind=TRUE) #index off diagonals to be refined
    
    if(matrices[m]=="CINT" | matrices[m]== 'TDPRED'| matrices[m]== 'T0MEANS'| 
        matrices[m]== 'T0MEANS'| matrices[m]== 'T0MEANS') freeindex <- 
      which(matrix(1,nrow=nrow(freematrix),ncol=ncol(freematrix))<999,arr.ind=TRUE) #index every cell
    
    if(matrices[m]=="MANIFESTVAR") freeindex <- which(diag(ncol(freematrix))==1,arr.ind=TRUE) #index diagonals
    
    if(matrices[m]=="TRAITVAR" | matrices[m]== 'DIFFUSION' | matrices[m]== 'MANIFESTTRAITVAR' | 
        matrices[m]== 'T0VAR') {
      symmetricfree<-TRUE
      blank<-diag(nrow(freematrix))
      blank[lower.tri(blank)]<-T
      freeindex <- which(blank==1,arr.ind=TRUE) #index every cell of the matrix
    }
    for(i in 1:nrow(freeindex)){ #for every parameter to be freed
      freematrix<-eval(parse(text=(paste0("ctmodelobj$",matrices[m]))))
     
      if(!is.na(as.numeric(freematrix[freeindex[i,,drop=F]])) & #if the parameter of interest is not already free
          as.numeric(freematrix[freeindex[i,,drop=F]]) == 0){  #and if it is currently set to 0
        
        testcount<-testcount+1
        newmodel<-ctmodelobj
        freedparamlabel <- paste0('freed',matrices[m],'_',freeindex[i,1],'_',freeindex[i,2])
        freedparam<-0
        names(freedparam)<-freedparamlabel
        startValues<-c(startValues,freedparam)
        newmodel$startValues<-startValues #add startValues from best accepted fit so far
        
        freematrix[freeindex[i,1],freeindex[i,2]] <- freedparamlabel #add free params to freematrix
        if(symmetricfree==TRUE) freematrix[freeindex[i,2],freeindex[i,1]] <- freedparamlabel #add symmetric free params
        message(paste0('Freeing parameter ',matrices[m],freeindex[i,1],freeindex[i,2]))
        eval(parse(text=paste0("newmodel$",matrices[m]," <- freematrix"))) #update the newmodel spec with the freed matrix of params
        
        newfit <- ctFit(datawide,newmodel,carefulFit=F,retryattempts=retryattempts,
          fit2beat=fit2beat,iterationSummary=TRUE,...)
        
        if(fitWithoutInits==TRUE) { #runs additional fit using default ctsem inits rather than best inits
        newmodelnoinits<-newmodel
        newmodelnoinits$startValues<-NULL
        newfitnoinits <- ctFit(datawide,newmodelnoinits,carefulFit=F,retryattempts=retryattempts,
          fit2beat=fit2beat,iterationSummary=TRUE,...)
          if(newfitnoinits[1]!="error" & class(newfitnoinits) != 'try-error'){
            if(newfit[1]=="error" | class(newfit)=='try-error') {
              newfit<-newfitnoinits #if bestinitfit was error and noinit wasn't, use noinit
            } else  #also use noinit fit if its fit was better
              if(summary(newfitnoinits$mxobj)$Minus2LogLikelihood < 
                  summary(newfit$mxobj)$Minus2LogLikelihood)  newfit <- newfitnoinits 
          }
        }
              
              
        if(newfit[1]=="error" | class(newfit)=='try-error'){ #if an error occurs, output as text
          message('Fit error')
          AIC[testcount]<-"error"
          freedParam[testcount]<-"error"
          estParams[testcount]<-"error"
          Minus2LogLikelihood[testcount] <- "error"
          fits[testcount]<-"error"
          models[testcount]<-"error"
        }
        
        if(newfit[1]!="error" & class(newfit) != 'try-error'){
          AIC[testcount]<-summary(newfit$mxobj)$AIC.Mx
          freedParam[testcount]<-paste0(matrices[m],'_',freeindex[i,1],'_',freeindex[i,2])
          estParams[testcount]<-summary(newfit$mxobj)$estimatedParameters
          Minus2LogLikelihood[testcount] <- summary(newfit$mxobj)$Minus2LogLikelihood
          message('refined -2LL = ', Minus2LogLikelihood[testcount], ', old -2LL = ',Minus2LogLikelihood[1])
          message('refined AIC = ', AIC[testcount], ', old AIC = ', AIC[1], '\n')
          fits[[testcount]]<-newfit
          models[[testcount]]<-newmodel
          
        }      
      }
    } 
    
  summary<-cbind(AIC,freedParam,estParams,Minus2LogLikelihood)
  out<-list()
  out$summary<-summary
  out$fits<-fits
    out$models<-models
  }
    
  if(min(AIC,na.rm=TRUE) < AIC[1] + tolerance){ #if one of the constrained models fits better than tolerance
    fit2beat=min(Minus2LogLikelihood,na.rm=T)
    
    if(fastrefine==TRUE){
      newparams<-freedParam[which(AIC < AIC[1] +tolerance)]
      newparams<-newparams[newparams!='none']
      part1<-gsub('_[0-9]+','',newparams)
      index<-gsub('[A-Za-z]+_','',newparams)
    index1<-gsub('_[0-9_]+','',index)
      index2<-gsub('[0-9_]+_','',index)
      
      newmodel<-models[[1]]
      for(k in 1:length(part1)){
        newmodel[[part1[k]]][as.numeric(index1[k]),as.numeric(index2[k])] <-
          paste0(part1[k],'_',index1[k],'_',index2[k])
    }
      message('Accepted new params:  ', paste0(newparams,'  '))
    }
    
    if(fastrefine==FALSE) { newmodel<-models[[which(AIC==min(AIC,na.rm=T))[1]]]    
    message('Accepted new param:  ', freedParam[which(AIC==min(AIC,na.rm=T))], '.  New -2LL = ', fit2beat)
    }
    if(nrow(out$summary > 1)) ctRefinePlot(out) #plot current param history
    smallout<-ctRefineFrom(datawide,
      newmodel,
      matrices=matrices,
      existingfit=fits[[which(AIC==min(AIC,na.rm=T))[1]]],
      retryattempts=retryattempts,...)
    out$summary<-rbind(out$summary,smallout$summary)
    out$models<-c(out$models,smallout$models)
    out$fits<-c(out$fits,smallout$fits)
  }
  
  return(out)
}
