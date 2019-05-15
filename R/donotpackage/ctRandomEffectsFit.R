#' Fits a random effects continuous time model.
#' 
#' Fits a single continuous time structural equation models to multiple groups (where each group contains 1 or more subjects),
#' by default, all parameters are free across groups.  Can also be used to easily estimate seperate models for each group.
#' 
#' @param datawide Wide format data, as used in \code{\link{ctFit}}.  See \code{\link{ctLongToWide}} to
#' easily convert long format data.
#' @param groupings Vector of character labels designating group membership for each row of datawide.  
#' These will be prefixed on relevant parameter estimates in the summary.
#' @param ctmodelobj Continuous time model to fit, specified via \code{\link{ctModel}} function.
#' @param fixedmodel Modified version of ctmodelobj, wherein any parameters you wish to keep 
#' fixed over groups should be given the value 'groupfixed'.  
#' If specified, all other parameters will be free across groups.
#' @param freemodel Modified version of ctmodelobj, wherein any parameters you wish to free across groups
#' should be given the label 'groupfree'.  
#' If specified, all other parameters will be fixed across groups.  
#' If left NULL, the default, all parameters are free across groups.
#' @param showInits Displays start values prior to optimization
#' @param carefulFit if TRUE, first fits the specified model with a penalised likelihood function 
#' to encourage parameters to remain closer to 0, then
#' fits the specified model normally, using these estimates as starting values. 
#' Can help with optimization in some cases, though results in user specified inits being ignored for the final fit.
#' @param retryattempts Number of fit retries to make.
#' @param ... additional arguments to pass to \code{\link{ctFit}}.
#' @return Returns an OpenMx fit object.
#' @details Additional \code{\link{ctFit}} parameters may be specified as required. Confidence intervals for any matrices and or parameters 
#' may be estimated afer fitting using \code{\link{ctCI}}.
#' 
#' @examples 
#' \donttest{
#' 
#' #Two group model, all parameters except LAMBDA[3,1] constrained across groups.
#' data(ctExample4)
#' basemodel<-ctModel(n.latent=1, n.manifest=3, Tpoints=20,
#'                    LAMBDA=matrix(c(1, 'lambda2', 'lambda3'), nrow=3, ncol=1),
#'                    MANIFESTMEANS=matrix(c(0, 'manifestmean2', 'manifestmean3'), 
#'                    nrow=3, ncol=1), TRAITVAR = 'auto')
#' 
#' freemodel<-basemodel
#' freemodel$LAMBDA[3,1]<-'groupfree'
#' groups<-paste0('g',rep(1:2, each=10),'_')
#' 
#' multif<-ctMultigroupFit(datawide=ctExample4, groupings=groups,
#'                        ctmodelobj=basemodel, freemodel=freemodel)
#' summary(multif)
#' 
#' 
#' 
#' #fixed model approach
#' fixedmodel<-basemodel
#' fixedmodel$LAMBDA[2,1]<-'groupfixed'
#' groups<-paste0('g',rep(1:2, each=10),'_')
#' 
#' multif<-ctMultigroupFit(datawide=ctExample4, groupings=groups,
#'                        ctmodelobj=basemodel, fixedmodel=fixedmodel)
#' summary(multif) 
#'}
#' 
#' 
#' @seealso \code{\link{ctFit}} and \code{\link{ctModel}}
#' @export


ctRandomEffectsFit<-function(datawide,groupings,ctmodelobj,fixedmodel=NA,freemodel=NA,
 carefulFit=FALSE,
  retryattempts=5,showInits=FALSE,...){

  if(any(suppressWarnings(!is.na(as.numeric(groupings))))) stop("grouping variable must not contain purely numeric items")
  if(length(groupings)!= nrow(datawide)) stop('length of groupings does not equal number of rows of datawide')
  
  if(all(is.na(fixedmodel))) fixedmodel<-ctmodelobj #so that it is not null or na
  
  startparams<-c() #to fill as needed
  
  
  
  omxmodels<-list() #blank list preparing for model input
  for(i in unique(groupings)){ #for every specified group
    #     
    singlegroup<-datawide[which(groupings == i),,drop=F] #data for the group
    
    singlectspec<-ctmodelobj
    
    if(all(is.na(freemodel))) freemodel <- lapply(ctmodelobj,function(x) { x<-rep('groupfree',length(x))}) #if no freemodel specified then free all params at this point
    
    for(m in 1:length(ctmodelobj)) { #for every element of the ctmodelobj list
      if(is.matrix(ctmodelobj[[m]])){ #if the element is a matrix
      for(j in 1:length(ctmodelobj[[m]])){ #for every slot in the matrix
        if(freemodel[[m]][j]=="groupfree"){ #if the slot is free in freemodel and not fixed in fixedmodel
        jnum<-suppressWarnings(as.numeric(ctmodelobj[[m]][j])) #check if it is numeric
        if(!is.na(ctmodelobj[[m]][j]) && is.na(jnum)) { #if the slot is neither NA or fixed to a value, then
          singlectspec[[m]][j] <- paste0(i,'_',ctmodelobj[[m]][j]) #give the label a group specific prefix
          
          
        }
      }
        if(any(!is.na(fixedmodel))){
        if(fixedmodel[[m]][j]=="groupfixed"){
          jnum<-suppressWarnings(as.numeric(ctmodelobj[[m]][j])) #check if it is numeric
          if(!is.null(ctmodelobj[[m]][j]) && is.na(jnum)) { #if the slot is neither null or fixed to a value, then
            singlectspec[[m]][j] <- paste0(ctmodelobj[[m]][j]) #use the global label
          }
        }
      }
    }
      }
    }
   

    
    
    
    if(carefulFit==TRUE) message('Begin carefulFit start value estimation for group ', i)
    
    omxmodel<-ctFit(singlegroup,singlectspec,nofit=TRUE, carefulFit=carefulFit,...) #omxmodel for group i
    ctfitargs<-omxmodel$ctfitargs
    omxmodel<-omxmodel$mxobj
    
    if(carefulFit==TRUE) {
      startparams<-c( startparams[ !( names(startparams) %in%  #get inits
          names(OpenMx::omxGetParameters(omxmodel))) ], #that are not found in the new fits inits
        OpenMx::omxGetParameters(omxmodel) ) #and combine the two vectors
    }
    
    omxmodel<- OpenMx::mxRename(omxmodel, newname=i) #change name of omxmodel for group i
      
  omxmodels[[i]]<-omxmodel #if fitting single multigroup model, add omxmodel for group i to list of omxmodels for all groups
    
  } #end loop over groups
  
  
  
    
    fullmodel <- OpenMx::mxModel('ctsem multigroup', #output multigroup omxmodel
      mxFitFunctionMultigroup(c(paste0(unique(groupings)))),
#       mxComputeSequence(list(
#         mxComputeGradientDescent(gradientAlgo="central", nudgeZeroStarts=FALSE, 
#           maxMajorIter=1000, gradientIterations = 1),
        # mxComputeReportDeriv(),
      omxmodels)
    

#     hyperpars<-NULL
#     if(!is.null(hyperpars)){ 
      
            params<-OpenMx::omxLocateParameters(fullmodel) #get list of parameters from base model
            
            baseparamindices<-grep(paste0(groupings[1],'_'),params$label) #get base parameters by removing group 1 id from any params with group 1 id
            baseparams <- params$value[baseparamindices]
            names(baseparams) <- gsub(paste0(groupings[1],'_'),'',params$label[baseparamindices])
            baseparams<-baseparams[!duplicated(baseparams)]
            
            meanParams<-lapply(1:nrow(datawide),function(x) {
              basemodel<-omxmodels[[x]]
              OpenMx::omxGetParameters(basemodel)[paste0(groupings[x],'_',names(baseparams)) %in% names(OpenMx::omxGetParameters(basemodel))]
            })
            meanParams<-colMeans(matrix(unlist(meanParams),ncol=length(meanParams[[1]])))
      
      penalisedmodels<-lapply(1:nrow(datawide), function(x) {
        
        basemodel<-omxmodels[[x]]
        indParams<-OpenMx::omxGetParameters(basemodel)[paste0(groupings[x],'_',names(baseparams)) %in% names(OpenMx::omxGetParameters(basemodel))]
        
        algstring<-paste0("mxAlgebra(name='algFit', ",groupings[x],".objective + sum(
          (indParams - meanParams) * (indParams - meanParams)))")
        fitAlg<-eval(parse(text=algstring))
        
        penalisedmodel<-OpenMx::mxModel(paste0('penalised_',groupings[x]), 
          
          basemodel,
          
          mxMatrix(name='meanParams', 
            values=meanParams, 
            labels=paste0('mean_',names(baseparams)),
            free=T,nrow=1,ncol=length(meanParams)),
          
#           mxMatrix(name='hyperpars',
#             values=hyperpars,
#             labels=paste0('penalty_',names(baseparams)),
#             free=F, nrow=1, ncol=length(meanParams)),
          
          mxMatrix(name='indParams', 
            labels=names(indParams), #paste0(groupings[x],'_',names(meanParams)),
            values=indParams,
            free=T,nrow=1,ncol=length(meanParams)),
          
          fitAlg,
          
          mxFitFunctionAlgebra('algFit')
        )
      })
      
      fullmodel<-OpenMx::mxModel('multigroup_ctsem',
        penalisedmodels,
        
        mxMatrix(name='meanParams', 
          values=meanParams, 
          labels=paste0('mean_',names(baseparams)),
          free=T,nrow=1,ncol=length(baseparams)),
        
#         mxMatrix(name='hyperpars',
#           values=hyperpars,
#           labels=paste0('penalty_',names(baseparams)),
#           free=F, nrow=1, ncol=length(baseparams)),
        
        mxFitFunctionMultigroup(paste0('penalised_',groupings))
      )
    
    ###old approach to hypervars
#     if(!is.null(hyperpars)){ #extracts global parameter names for params with variance limits
#       params<-omxLocateParameters(fullmodel) #get list of parameters from base model
#       
#       baseparamindices<-grep(paste0(groupings[1],'_'),params$label) #get base parameters by removing group 1 id from any params with group 1 id
#       baseparams <- params$value[baseparamindices]
#       names(baseparams) <- gsub(paste0(groupings[1],'_'),'',params$label[baseparamindices])
#       baseparams<-baseparams[!duplicated(baseparams)]
# 
#       indparams<-c()
#       for(i in 1:length(baseparams)){ #get individual parameter labels
#         tempparams<-params [grep( #param list where
#           names(baseparams)[i], #main param string i is found
#           params$label),] #in param labels
#         sub1<-gsub('^[[:alpha:]]','',tempparams$label)
#         sub2<-as.numeric(gsub( '[_](.*)' ,'', sub1))
#         indparams<-rbind(indparams,tempparams[order(sub2),])
#       }
# 
#       indparamlabels<-list()
#       for(i in 1:length(baseparams)){ #get individual parameter labels
#         templabels<-indparams$label [grep( #param list where
#           names(baseparams)[i], #main param string i is found
#           indparams$label)] #in param labels
#         
#         templabels<-templabels[!duplicated(templabels)]
# 
#         sub1<-gsub('^[[:alpha:]]','',templabels)
#         sub2<-as.numeric(gsub( '[_](.*)' ,'', sub1))
#         indparamlabels[[i]]<-templabels[order(sub2)]
#         
#       }
#       indparamlabels<-matrix(unlist(indparamlabels),nrow=length(unique(groupings)),ncol=length(baseparams))
#       
#       indparamvalues<-list()
#       for(i in 1:length(baseparams)){ #get individual parameter labels
#         tempvalues<-indparams$value [grep( #param list where
#           names(baseparams)[i], #main param string i is found
#           indparams$label)] #in param labels
#         
#         tempvalues<-tempvalues[!duplicated(tempvalues)]
#         
#         sub1<-gsub('^[[:alpha:]]','',templabels)
#         sub2<-as.numeric(gsub( '[_](.*)' ,'', sub1))
#         indparamvalues[[i]]<-tempvalues[order(sub2)]
#         
#       }
#       indparamvalues<-matrix(unlist(indparamvalues),nrow=length(unique(groupings)),ncol=length(baseparams))
#       
#       
#       #for basing parameters off first parameter (to help optimization)
#       
# #       indparamswithoutfirst<-indparams[-grep(paste0(groupings[1],'_'),indparams$label),]
# #       indparamlabelswithoutfirst<-indparamlabels[-1,]
# #       
# #       indparammodifiermatrix <- OpenMx::mxMatrix(name='indparammodifiermatrix', 
# #         free=c(FALSE, rep(TRUE,length(unique(groupings))-1)),
# #         nrow=length(unique(groupings)), 
# #         values=c(0,stats::rnorm(length(unique(groupings))-1,0,.01)),
# #         ncol=length(hyperpars),
# #         labels=rbind(NA,indparamlabelswithoutfirst))
# #       
# #       indfirstparammatrix<-OpenMx::mxMatrix(name='indfirstparammatrix', 
# #         free=TRUE,
# #         nrow=length(unique(groupings)),
# #         values=rep(indparams$value[which(indparams$label %in% rep(indparamlabels[1,],each=length(unique(groupings))))],
# #           ,each=length(unique(groupings))),
# #         ncol=length(hyperpars),
# #         labels=rep(indparamlabels[1,],each=length(unique(groupings))))
# #       
# #       indparammodifieralg <- OpenMx::mxAlgebra(name='indparammodifieralg', indfirstparammatrix + indparammodifiermatrix)
# 
#   
#       
#       
#       
#       algstring<-c(paste0(unique(groupings),'.objective',collapse=' + '))
#       subfits<-eval(substitute(OpenMx::mxAlgebra(theexpression, name='subfits'),list(theexpression=parse(text=algstring)[[1]])))
# 
# 
#       fullmodel<-OpenMx::mxModel('hypervar',
# #         indparammodifiermatrix,
# #         indfirstparammatrix,
# #         indparammodifieralg,
# #         
#         mxMatrix(name='nparameters', 
#           ncol=1,
#           nrow=1,
#           values=length(baseparams), 
#           free=F),
#         
#         mxMatrix(type='Full',
#           values=length(unique(groupings)), 
#           ncol=1,
#           nrow=1,
#           free=F,
#           name='groupcount'),
# 
#         mxMatrix(name='sumMatrix',
#           values=1,
#           free=F,
#           nrow=1,
#           ncol=length(unique(groupings))),
#         
#         mxAlgebra(name='m', sumMatrix %*% (indParameters) %x% (1/groupcount),dimnames=list(NULL,names(baseparams))), #observed means vector
#         
#         mxMatrix(name='bigMeans',
#           labels=paste0(
#           'm[1,',
#           rep(1:length(baseparams), each= length(unique(groupings))),
#           ']'),
#           ncol=length(baseparams),
#           nrow=length(unique(groupings))),
#         
#         mxAlgebra(name='indparamDeviance', indParameters - bigMeans),
#         mxAlgebra(name='S',t(indparamDeviance) %*% indparamDeviance,dimnames=list(names(baseparams),names(baseparams))), #observed covariance matrix
# 
#         mxAlgebra(subfits + tr(abs (vec2diag(diag2vec(S)) * (vec2diag(hyperpars * hyperpars))) ), name='fitAlgebra'),
#         
#         subfits,
#         
#         omxmodels,
#         
#         mxMatrix(name='hyperpars', 
#           values=hyperpars, 
#           labels=names(hyperpars),
#           nrow=1,
#           ncol=length(hyperpars),
#           free=FALSE),
# 
#         mxFitFunctionAlgebra('fitAlgebra')
#         
#       )      #end L3 model spec
#       
# #       tm<-indparamswithoutfirst #adjust parameter matrices to reference algebra modifier
# #       for(i in 1:length(indparamswithoutfirst$label)){
# #         fullmodel[[ tm[i,'model'] ]] [[  tm[i,'matrix'] ]]$labels[tm[i,'row'], tm[i,'col'] ] <-
# #           paste0('hypervar.indparammodifieralg[', which(indparammodifiermatrix$labels == tm[i,'label'], arr.ind=TRUE)[1], 
# #             ',',
# #             which(indparammodifiermatrix$labels == tm[i,'label'], arr.ind=TRUE)[2],
# #             ']')
# #         
# #         fullmodel[[ tm[i,'model'] ]] [[  tm[i,'matrix'] ]]$free[tm[i,'row'], tm[i,'col'] ] <- FALSE
# #         fullmodel[[ tm[i,'model'] ]] [[  tm[i,'matrix'] ]]$values[tm[i,'row'], tm[i,'col'] ] <- stats::rnorm(1,0,.01)
# #       }
# 
#       fullmodel<-mxModel(fullmodel, 
#         mxMatrix(name='indParameters',
#         nrow=nrow(indparamlabels), 
#           ncol=length(baseparams),
#           # values=rep(baseparams,each=length(unique(groupings))) + stats::rnorm(length(unique(groupings))*length(hyperpars),0,.00),
#           values=indparamvalues,
#         labels=indparamlabels, 
#           free=TRUE))
#       
#       # fullmodel<-omxSetParameters(fullmodel,values=params,labels=names(params))
#       
#       
# 
# #       ####reparameterise in terms of deviations from first individual parameter
# # 
# #       indparammodifiers<-list()
# #       meanparammatrix<-list()
# #       for(i in 1:length(baseparams)){
# #         #these replace the original individual param labels, to reference the individual param algebra
# #         indparammodifierrefs<- paste0('hypervar.indparammodifieralg','[', 1:(length(unique(groupings))),',',i, ']')
# # 
# #         modifyindex<- grep(names(baseparams)[i],indparams$label,perl=T) #for every line related to the base param 
# #         for(x in 1:length(modifyindex)){  
# #           fullmodel[[indparams$model[modifyindex[x]]]][[indparams$matrix[modifyindex[x]]]]$labels[indparams$row[modifyindex[x]],indparams$col[modifyindex[x]]] <- #the related part of the model
# #             indparammodifierrefs[x] #is replaced by the relevent indparammodifierref
# #           fullmodel[[indparams$model[modifyindex[x]]]][[indparams$matrix[modifyindex[x]]]]$free[indparams$row[modifyindex[x]],indparams$col[modifyindex[x]]] <- #the related part of the model
# #             FALSE #is set to fixed
# #         }
# #       }
# # 
# # #this is a matrix of deviations from first individual param
# # indparammodifiers<-mxMatrix(name=paste0('indparammodifiers'),
# #   labels=unlist(indparamlabels),
# #   free=TRUE,
# #   values=rep(baseparams,each=length(unique(groupings))) + stats::rnorm(length(unique(groupings))*length(hyperpars),0,.01),
# #     nrow=length(unique(groupings)),
# #   ncol=length(hyperpars))      
# # 
# #       meanparammatrix<-mxMatrix(name=paste0('meanparammatrix'),
# #         labels=rep(names(baseparams),each=length(unique(groupings))),
# #         nrow=length(unique(groupings)),
# #         ncol=length(hyperpars),
# #         values=rep(baseparams,each=length(unique(groupings))), 
# #         free=TRUE)
# #         
# #       theexpression<-paste0('meanparammatrix + indparammodifiers')
# #       #This algebra adds first individual param to every individual param modifier
# #       indparammodifieralg<- eval(substitute(mxAlgebra(name=paste0('indparammodifieralg'), 
# #         expression=theexpression),
# #         list(theexpression=parse(text=theexpression)[[1]])))
# #       
# #       
# #       fullmodel<-mxModel(fullmodel,indparammodifieralg,indparammodifiers, meanparammatrix)
#       
      
      
      
      
      
    # }  #end variance constraint section


    fullmodel<-OpenMx::omxAssignFirstParameters(fullmodel)

    # if(!is.null(confidenceintervals)) fullmodel <- OpenMx::mxModel(fullmodel, mxCI(confidenceintervals,interval = 0.95,type = "both")) #if 95% confidence intervals are to be calculated

      fullmodel<-OpenMx::mxTryHard(fullmodel,initialTolerance=1e-14,
      showInits=showInits,
      bestInitsOutput=FALSE,
      extraTries=retryattempts,loc=1,scale=.2,paste=FALSE,...) 
    
      fullmodel<-list(mxobj=fullmodel, ctfitargs=ctfitargs, ctmodelobj=ctmodelobj, groups=unique(groupings))
      class(fullmodel)<-'ctsemMultigroupFit'
      
      
    return(fullmodel)
  
}


