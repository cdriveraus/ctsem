ctCompare <- function(fits,leaveOutNseq=1:2){
  #consider check for same data
  if(is.null(names(fits))) names <- paste0('fit',1:length(fits)) else names <- names(fits)
  
  stats <- lapply(fits,function(f){
    out <- data.frame(
      ll=f$stanfit$transformedparsfull$ll,
      lp=f$stanfit$optimfit$value,
      np=length(f$stanfit$rawest)
    )
    out[['aic']]<- 2* out$np - 2*out$lp
    return(out)
  })
  names(stats) <- names
  stats <- rbindlist(stats,idcol = 'Model')
  
  
  llrows <- do.call(what = cbind,args = lapply(fits,function(f) f$stanfit$transformedparsfull$llrow[1,]))
  colnames(llrows) <- names
  llrows<-data.table(row=1:nrow(llrows),llrows)
  meltids <- 'row' 
  
  if(!all(is.na(leaveOutNseq))){ #compute llrows with left out rows
    llrowsCV <- list()
    for(leaveN in leaveOutNseq){
      llrowsCV[[length(llrowsCV)+1]] <- lapply(fits,function(f){
        try(suppressMessages(ctLOO(fit = f,cores = 1,leaveOutN = leaveN,refit = FALSE)))
      })
      llrowsCV[[length(llrowsCV)]] <- llrowsCV[[length(llrowsCV)]][!sapply(llrowsCV[[length(llrowsCV)]],class) %in% 'try-error']
    }
    names(llrowsCV) <- leaveOutNseq
    
    lpCV <- rbindlist(idcol = 'leaveOutN',fill=TRUE, #get left out log probs
      lapply(llrowsCV,function(x) data.table(Model=names(x),sapply(x,function(xx) xx$insampleLogProb)))
    )
    
    llrowsCV <- lapply(llrowsCV,function(leftn) lapply(leftn, function(modeli) modeli$outsampleLogLikRow))
    
    llrowsCV <- rbindlist(idcol='leaveOutN',fill=TRUE, #get all the llrow columns and organise them
      lapply(llrowsCV,function(x){
        # names <- names[!sapply(x,class) %in% 'try-error']
        # x=x[!sapply(x,class) %in% 'try-error']
        lltemp <- do.call(what = cbind,args = x)
        colnames(lltemp) <- names(x)
        lltemp <- data.table(row=1:nrow(llrows),lltemp)
      })
    )

    llrows <- rbind(cbind(leaveOutN=0,llrows),llrowsCV,fill=TRUE)
    meltids <- c(meltids,'leaveOutN')
    
    #compute stats for each left out N
    leftstats <- rbindlist(fill=TRUE,lapply(leaveOutNseq,function(leftn){
      rbindlist(fill=TRUE,
        lapply(names,function(namei){
          message(namei)
          message(leftn)
          if(length(lpCV[leaveOutN==leftn & Model==namei,V2]) == 0){
            return(NULL)
          } else {
            out <- data.frame(Model=namei,leaveOutN=leftn,
              ll=sum(llrows[leaveOutN==leftn,(namei),with=F],na.rm=TRUE),
              lp=lpCV[leaveOutN==leftn & Model==namei,V2],
              np=length(fits[[which(names(fits) %in% namei)]]$stanfit$rawest)
            )
            return(out)
          }
        }))
    }))
    leftstats[,aic:=2* np - 2*lp]
    stats <- rbind(stats[,leaveOutN:=0],leftstats)
  }
  
  mllrows <- melt.data.table(data.table(llrows),id.vars = meltids)
  return(list(stats=stats,llrows=mllrows))
}


#' Chi Square test wrapper for ctStanFit objects.
#'
#' @param fit1 One of the fits to be compared (better fit is assumed as base for comparison)
#' @param fit2 Second fit to be compared
#'
#' @return Numeric probability
#' @export
#'
#' @examples
#' \donttest{
#'     df <- data.frame(id=1, time=1:length(sunspot.year), Y1=sunspot.year)
#'     
#'     m1 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0)
#'     m2 <- ctModel(type='dt', LAMBDA=diag(1),MANIFESTVAR=0,DRIFT = .9)
#'     
#'     f1 <- ctStanFit(df,m1,cores=1)
#'     f2 <- ctStanFit(df,m2,cores=1)
#'     
#'     ctChisqTest(f1,f2)
#' }

ctChisqTest<-function(fit1,fit2){
  d <- data.frame(
    ll=c(fit1$stanfit$optimfit$value,
      fit2$stanfit$optimfit$value),
    npars=c(length(fit1$stanfit$rawest),
      length(fit2$stanfit$rawest))
  )
  d <- d[order(d$npars,decreasing = FALSE),]
  stats::pchisq(q =  diff(2*d$ll),df = diff(d$npars),lower.tail = FALSE)
}
