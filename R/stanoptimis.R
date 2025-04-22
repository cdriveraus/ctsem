

# # states from sigpoints
# gensigstates <- function(mu, ch){
#   ndim=nrow(ch)
#   asquared=1
#   ukfspread = sqrt(1) / (2/sqrt(ndim*2)) #ensuring asquared is 1
#   k=.5
#   b=0
#   lambda=asquared * ( ndim +k)-ndim*2
#   sqrtukfadjust = sqrt(ndim*2 +lambda)
#   
#   sigpoints <- t(chol(as.matrix(Matrix::bdiag((ch%*%t(ch))))))*sqrtukfadjust
#   
#   t0states <- matrix(mu,byrow=TRUE,ndim*2+2,ndim)
#   t0states[3:(2+ndim),] =  t0states[3:(2+ndim),] + t(sigpoints)
#   t0states[(3+ndim):(ndim*2+2),] = t0states[(3+ndim):(ndim*2+2),] - t(sigpoints)
#   return(t0states)
# }


# hesscalc <- function(sf,step){ #input stanfit, get hessian
#   smf <- stan_reinitsf(sf$stanmodel,sf$standata)
#   fgfunc <- function(x) rstan::log_prob(smf,x,gradient=TRUE)
#   # step=findstepsize(est = sf$stanfit$rawest,func = fgfunc,eps = 1e-2,
#   #   maxlpd = 2*length(sf$stanfit$rawest),minlpd=.1*length(sf$stanfit$rawest))
#   H=jac(pars = sf$stanfit$rawest,fgfunc = fgfunc,step = rep(step,length(sf$stanfit$rawest)))
# }

# findstepsize <- function(est,func,eps=1e-3,maxlpd=3,minlpd=1e-3){
#   
#   fmax=func(est)
#   epsbase=.01
#   i=0
#   n<-100
#   devs <- matrix(NA,nrow=n,ncol=length(est))
#   lp<-c()
#   epscount <- 0
#   goodcount <- 0
#   
#   while(i < n && goodcount < 5){
#     i=i+1
#     accepted <- FALSE
#     devs[i,] <- rnorm(length(est))*eps
#     count <- 0
#     while(!accepted){
#       count = count + 1
#       if(count==20) stop('Couldnt determine step size!')
#       fg<-try(func(devs[i,]+est))
#       if('try-error' %in% class(fg)) eps <- eps * .5+1e-6 else accepted <- TRUE
#     }
#     
#     bad=(lp<quantile(lp,.2,na.rm = TRUE) & (fmax-lp) > maxlpd) | 
#       (lp>quantile(lp, .95,na.rm = TRUE) & (fmax-lp) < minlpd*.001)
#     
#     good <- c(1:i)[!bad]
#     lp <- c(lp,fg[1])
#     fl=lp-min(lp)
#     fl <- sqrt(fl)
#     
#     frange <- fmax[1]-lp[i]
#     print(frange)
#     if(length(good) <=5){
#       if(frange < minlpd) eps <- eps * 2
#       if(frange > maxlpd) eps <- eps /2
#     }
#     # message(i)
#     # print(eps)
#     if(frange > minlpd && frange < maxlpd) goodcount <- goodcount + 1 else goodcount <- 0
#     # 
#     if(length(good)>5){
#       if(frange < minlpd) epsbase <- epsbase * 2
#       if(frange > maxlpd) epsbase <- epsbase / 2
#       eps=epsbase/(sqrt(abs(apply(devs,2,function(x) coefficients(lm(fl[good]~abs(x[good])))[-1])))+1e-8)
#       epscount <- epscount + 1
#       print(epsbase)
#     }
#   }
#   return(eps)
# }

# gradcheck <- function(sf,min=1e-5,seqlength=10,whichpars=NA,
#   offdiags=FALSE,scale=TRUE,finite=FALSE){
#   p <- sf$stanfit$rawest
#   smf <- stan_reinitsf(sf$stanmodel,sf$standata)
#   if(is.na(whichpars[1])) whichpars <- 1:length(p)
#   
#   func<-function(x) log_prob(smf,x)
#   b=func(p)
#   
#   for(pi in whichpars){
#     g <- list()
#     gf=c()
#     devs <- min*2^(1:seqlength)
#     devs=sort(c(devs,-devs,0))
#     for(i in seq_along(devs)){
#       px <- p
#       px[pi] <- p[pi] + devs[i]
#       g[[i]] <- log_prob(smf,px,gradient=TRUE)
#       if(finite) gf[i]<-(func(px)-b)/devs[i]
#     }
#     lp=sapply(g,function(x) x[1])
#     plot(devs,lp,type='l',lwd=2,main=pi)
#     abline(h=lp[devs==0],lty=2)
#     abline(v=0,lty=2)
#     gmat<- sapply(g,function(x) attributes(x)$gradient)
#     if(finite) gfmat=sapply(gf,function(x) x*2)
#     if(scale) gmat <- t(t(gmat)/devs)
#     if(finite && scale) gfmat=t(t(gfmat)/devs)
#     plot(devs,gmat[pi,],type='l',lwd=2,main=pi,ylab='grad')
#     if(finite)  points(devs,gfmat,type='l',lwd=2,lty=2,col=2,main=pi,ylab='grad')
#     abline(v=0,lty=2)
#     abline(h=0,lty=2)
#     if(offdiags){
#       for(i in 1:length(p)){
#         if(i != pi) plot(devs,gmat[i,],type='l',lwd=2,main=paste0(pi,'  ',i),ylab='grad')
#         abline(v=0,lty=2)
#         abline(h=0,lty=2)
#       }
#     }
#   }
# }


# jacrandom <- function(grfunc, est, eps=1e-4,
#   maxlpd=.5, minlpd=1e-5){
#   
#   fmax=grfunc(est)
#   n <- max(length(est)*2,200)
#   mcholtest=c()
#   epsbase=.0001
#   mcholtest <- FALSE
#   i=0
#   epsadapted <- FALSE
#   pd<-rep(FALSE,n)
#   H <- NA
#   
#   while(i < n && suppressWarnings(min(which(pd))) > (i-20)){
#     # print(i)
#     if(i==0){
#       devs <- matrix(NA,nrow=n,ncol=length(est))
#       covstore <- devs
#       fg<-list()
#       f<-c()
#       grm <- devs
#     }
#     i=i+1
#     # 
#     accepted <- FALSE
#     trycount <-0
#     while(!accepted){
#       trycount <- trycount + 1
#       if(trycount > 20) stop('Errors encountered estimating Hessian')
#       devs[i,] <- rnorm(length(est),0,eps)
#       fg[[i]] <- try(grfunc((devs[i,])+est))
#       if('try-error' %in% class(fg[[i]])) eps <- eps * .5 else accepted <- TRUE
#     }
#     
#     
#     grm[i,] <- attributes(fg[[i]])$gradient
#     f[i] <- fg[[i]][1]
#     
#     # plot(f)
#     bad=(f<quantile(f,.2,na.rm = TRUE) & (fmax-f) > maxlpd) | 
#       (f>quantile(f, .95,na.rm = TRUE) & (fmax-f) < minlpd*.001)
#     # if(i > length(est)+5) bad <- unique(c(1:(i-length(est)-5)),(bad))
#     
#     good <- c(1:i)[!bad]
#     fl=f-min(f)
#     fl <- sqrt(fl)
#     
#     if(length(good)>5){
#       frange <- fmax-f[i]
#       if(frange < minlpd) epsbase <- epsbase * 1.1
#       if(frange > maxlpd) epsbase <- epsbase / 1.1
#       eps=epsbase/abs(apply(devs,2,function(x) coefficients(lm(fl[good]~abs(x[good])))[-1]))
#       # if(length(good)>length(est)){
#       #   eps=epsbase/abs(coefficients(lm(fl[good]~abs(devs[good,])))[-1])
#       #   }
#       # print(epsbase)
#       
#       if(i==10 && !epsadapted){
#         epsadapted <- TRUE
#         i <- 0
#       }
#     }
#     # print(round(epsbase,4))
#     
#     if(length(good) > length(est) && i > length(est)) {
#       lf=t(apply(grm[good,],2,function(y) coefficients(lm((y) ~ (devs[good,])))[-1])) #without intercept
#       #     lf=t(apply(grm[good,],2,function(y){
#       #       apply(devs[good,],2,function(x) coefficients(lm((y) ~ x-1)))
#       #       }))
#       # print(i)
#       
#       lf=(lf+t(lf))/2
#       # mc=MASS::ginv(-lf)
#       # mc=as.matrix(Matrix::nearPD(mc)$mat)
#       mchol = try(t(chol(solve(-lf))),silent=TRUE)
#       pd[i] <- !'try-error' %in% class(mchol)
#       if(pd[i]) H <- lf
#       # covstore[i,] <- diag(mc)
#     }
#   }
#   if(is.na(H[1])) H <- lf #return last non pos def if no good ones found
#   return(H)
# }


#' Sample more values from an optimized ctstanfit object
#'
#' @param fit fit object
#' @param nsamples number of samples desired
#' @param cores number of cores to use
#'
#' @return fit object with extra samples
#' @export
#'
#' @examples
#' \dontrun{
#' newfit <- ctAddSamples(ctstantestfit, 10, 1)
#' }
ctAddSamples <- function(fit,nsamples,cores=2){
  
  if(length(fit$stanfit$stanfit@sim) > 0) stop('ctStanFit object was sampled and not optimized, cannot add samples!')
  
  mchol <- t(chol(fit$stanfit$cov))
  resamples <- matrix(unlist(lapply(1:nsamples,function(x){
    fit$stanfit$rawest + (mchol) %*% t(matrix(rnorm(length(fit$stanfit$rawest)),nrow=1))
  } )),byrow=TRUE,ncol=length(fit$stanfit$rawest))
  
  fit$stanfit$rawposterior <- rbind(fit$stanfit$rawposterior,resamples)
  
  fit$stanfit$transformedpars=stan_constrainsamples(sm = fit$stanmodel,
    standata = fit$standata,samples=fit$stanfit$rawposterior,
    savescores = fit$standata$savescores,
    savesubjectmatrices=as.logical(fit$standata$savesubjectmatrices),
    dokalman=as.logical(fit$standata$savesubjectmatrices),
    cores=cores)
  return(fit)
}

parallelStanSetup <- function(cl, standata,split=TRUE,nsubsets=1){
  cores <- length(cl)
  if(split) stanindices <- split(unique(standata$subject),(unique(standata$subject) %% min(standata$nsubjects,cores))) #disabled sorting so subset works parallel
  if(!split) stanindices <- lapply(1:cores,function(x) unique(standata$subject))
  if(!split && length(stanindices) < cores){
    for(i in (length(stanindices)+1):cores){
      stanindices[[i]] <- NA
    }
  }
  
  standata$nsubsets <- as.integer(nsubsets)
  if(!split) cores <- 1 #for prior mod
  # clusterIDexport(cl,c('standata','stanindices','cores'))
  benv <- new.env(parent = globalenv())
  
  sapply(c('standata','stanindices','cores','cl'),function(x){
    assign(x,get(x),pos = benv)
    NULL
  })
  
  
  benv$commands <- list(
    "g = eval(parse(text=paste0('gl','obalenv()')))", #avoid spurious cran check -- assigning to global environment only on created parallel workers.
    "environment(parlptext) <- g",
    'if(standata$recompile > 0) load(file=smfile) else sm <- ctsem:::stanmodels$ctsm',
    'eval(parse(text=parlptext))',
    'assign("parlp",parlp,pos=g)',
    "if(length(stanindices[[nodeid]]) < length(unique(standata$subject))) standata <- ctsem:::standatact_specificsubjects(standata,stanindices[[nodeid]])",
    "standata$priormod <- 1/cores",
    "if(FALSE) sm=99",
    "smf=ctsem:::stan_reinitsf(sm,standata)"
  )
  
  eval(parse(text=
      "parallel::clusterExport(cl,c('standata','stanindices','cores','commands'),envir=environment())"),
    envir = benv)
  eval(parse(text="parallel::clusterEvalQ(cl,eval(parse(text=paste0(commands,collapse=';'))))"),envir = benv)
  eval(parse(text="parallel::clusterEvalQ(cl,sapply(commands,function(x) eval(parse(text=x))))"),envir = benv)
  
  # eval(parse(text="parallel::clusterEvalQ(cl,ls(envir=globalenv()))"),envir = benv)
  
  NULL
}

#create as text because of parallel communication weirdness
parlptext <- 
  'parlp <- function(parm){
     a=Sys.time()
          out <- try(rstan::log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
          

        
        if("try-error" %in% class(out) || any(is.nan(attributes(out)$gradient))) {
          outerr <- out
          out <- -1e100
          attributes(out)$gradient <- rep(NaN, length(parm))
          attributes(out)$err <- outerr
        }
        
        attributes(out)$time <- Sys.time()-a
        
        if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
        # attributes(out)$gradient[is.nan(attributes(out)$gradient)] <-
          # rnorm(length(attributes(out)$gradient[is.nan(attributes(out)$gradient)]),0,100)
          
        return(out)
        }'



#based on rstan function, very cut down, may fail in some cases...
#' @importFrom Rcpp cpp_object_initializer
getcxxfun <- function(object) {
  if (length(object@dso_saved) == 0){
    return(function(...) stop("this function should not be called"))
  }  else  return(object@.CXXDSOMISC$cxxfun)
}


#' Quickly initialise stanfit object from model and data
#'
#' @param model stanmodel
#' @param data standata
#' @param fast Use cut down form for speed
#'
#' @return stanfit object
#' @export
#'
#' @examples
#' sf <- stan_reinitsf(ctstantestfit$stanmodel,ctstantestfit$standata)
stan_reinitsf <- function(model, data,fast=FALSE){
  if(fast) sf <- new(model@mk_cppmodule(model),data,0L,getcxxfun(model@dso))
  
  if(!fast) suppressMessages(suppressWarnings(suppressOutput(sf<- 
      rstan::sampling(model,iter=0,chains=0,init=0,data=data,check_data=FALSE,
        control=list(max_treedepth=0),save_warmup=FALSE,test_grad=FALSE))))
  
  return(sf)
}

flexsapply <- function(cl, X, fn,cores=1){
  if(cores > 1) parallel::parSapply(cl,X,fn) else sapply(X, fn)
}

flexlapply <- function(cl, X, fn,cores=1,...){
  if(cores > 1) parallel::parLapply(cl,X,fn,...) else lapply(X, fn,...)
}

flexlapplytext <- function(cl, X, fn,cores=1,...){
  if(cores > 1) {
    nodeindices <- split(1:length(X), sort((1:length(X))%%cores))
    nodeindices<-nodeindices[1:cores]
    clusterIDexport(cl,c('nodeindices'))
    
    out <-unlist(clusterIDeval(cl,paste0('lapply(nodeindices[[nodeid]],',fn,')')),recursive = FALSE)
    # out2<-parallel::parLapply(cl,X,tparfunc,...) 
  } else out <- lapply(X, eval(parse(text=fn),envir =parent.frame()),...)
  return(out)
}


#' Adjust standata from ctsem to only use specific subjects
#'
#' @param standata standata
#' @param subjects vector of subjects
#' @param timestep ignored at present
#'
#' @return list of updated structure
#' @export
#'
#' @examples
#' d <- standatact_specificsubjects(ctstantestfit$standata, 1:2)
standatact_specificsubjects <- function(standata, subjects,timestep=NA){
  long <- standatatolong(standata)
  long <- long[long$subject %in% subjects,]
  standatamerged <- standatalongremerge(long=long, standata=standata)
  standatamerged$ndatapoints <- as.integer(nrow(long))
  if(standata$ntipred > 0) standatamerged$tipredsdata <- standatamerged$tipredsdata[unique(standatamerged$subject),,drop=FALSE]
  standatamerged$nsubjects <- as.integer(length(unique(standatamerged$subject)))
  standatamerged$subject <- array(as.integer(factor(standatamerged$subject)))
  standatamerged$idmap <- standata$idmap[standata$idmap$new %in% subjects,]
  return(standatamerged)
}  


standatalongobjects <- function() {
  longobjects <- c('subject','time','dokalmanrows','nobs_y','ncont_y','nbinary_y',#'nordinal_y','whichordinal_y',
    'Y','tdpreds', 'whichobs_y','whichbinary_y','whichcont_y')
  return(longobjects)
}

standatatolong <- function(standata, origstructure=FALSE,ctm=NA){
  long <- lapply(standatalongobjects(),function(x) as.matrix(standata[[x]]))
  names(long) <- standatalongobjects()
  
  if(origstructure){
    if(is.na(ctm[1])) stop('Missing ctm arg in standatatolong()')
    colnames(long[['Y']]) <- ctm$manifestNames#colnames(standata$Y)
    long[['Y']][long[['Y']] %in% 99999] <- NA
    colnames(long[['subject']]) <- ctm$subjectIDname
    colnames(long[['time']]) <- ctm$timeName
    longout <- data.frame(long[['subject']],long[['time']],long[['Y']])
    if(standata$ntdpred > 0){
      colnames(long[['tdpreds']]) <- colnames(standata$tdpreds)
      if(!is.na(ctm[1])) colnames(long[['tdpreds']]) <- ctm$TDpredNames
      longout <- cbind(longout,long[['tdpreds']])
    }
    if(standata$ntipred > 0){
      tipreds <- standata$tipredsdata[longout[[ctm$subjectIDname]],,drop=FALSE]
      tipreds[tipreds %in% 99999] <- NA
      if(!is.na(ctm[1])) colnames(tipreds) <- ctm$TIpredNames
      longout <- cbind(longout,tipreds)
    }
  } else longout <- data.frame(long) 
  
  #,simplify=data.frame(subject=standata$subject, time=standata$time
  # colnames(long)[colnames(long) %in% 'Y'] <- paste0('Y.1'
  # colnames(long)[colnames(long) %in% 'tdpreds'] <- 'tdpreds.1'
  return(longout)
}

# stanlongtostandatasml <- function(long){
#   standatasml <- sapply( standatalongobjects(), function(x) long[,grep(paste0('^',x),colnames(long)),drop=FALSE])
#   names(standatasml) <- standatalongobjects()
#   return(standatasml)
# }

standatalongremerge <- function(long, standata){ #merge an updated long portion of standata into original standata
  n = names(standata)
  standatamerged <- lapply(names(standata), function(x) {
    if(x %in% standatalongobjects()){
      objdims <- dim(standata[[x]])
      if(is.null(objdims)) objdims <- c()
      objdims[1] <- nrow(long)
      xdat <- unlist(long[,grep(paste0('^',x),colnames(long)),drop=FALSE])
      if(is.null(xdat)) xdat <- NA
      return(array(xdat, dim=objdims))
    } else return(standata[[x]])
  })
  names(standatamerged) <- n
  return(standatamerged)
}

standataFillTime <- function(standata, times, subject, maintainT0=FALSE){
  long <- standatatolong(standata)
  
  if(any(!times %in% long$time)){ #if missing any times, add empty rows
    nlong <- do.call(rbind,
      lapply(subject, function(si){
        mintime <- min(long$time[long$subject==si])
        originaltimes <- round(long$time[long$subject==si],10)
        stimes <- times[(!times %in% originaltimes)]
        if(maintainT0) stimes <-stimes[stimes > mintime]
        data.frame(subject=si,time=stimes)
      })
    )
    
    nlong <- suppressWarnings(data.frame(nlong,long[1,!colnames(long) %in% c('subject','time')]))
    nlong[,grep('(^nobs)|(^which)|(^ncont)|(^nbin)',colnames(nlong))] <- 0L
    nlong[,grep('^dokalman',colnames(nlong))] <- 1L
    nlong[,grep('^Y',colnames(nlong))] <- 99999
    nlong[,grep('^tdpreds',colnames(nlong))] <- 0
    
    long <- rbind(long,nlong)
  } #end empty rows addition
  
  long <- long[order(long$subject,long$time),]
  standatamerged <- standatalongremerge(long=long, standata=standata)
  standatamerged$ndatapoints <- as.integer(nrow(long))
  return(standatamerged)
}




stan_constrainsamples<-function(sm,standata, samples,cores=2, cl=NA,
  savescores=FALSE,
  savesubjectmatrices=TRUE,
  dokalman=TRUE,
  onlyfirstrow=FALSE, #ifelse(any(savesubjectmatrices,savescores),FALSE,TRUE),
  pcovn=2000,
  quiet=FALSE){
  if(savesubjectmatrices && !dokalman){
    dokalman <- TRUE
    warning('savesubjectmatrices = TRUE requires dokalman=TRUE also!')
  }
  standata$savescores <- as.integer(savescores)
  standata$dokalman <- as.integer(dokalman)
  standata$savesubjectmatrices<-as.integer(savesubjectmatrices)
  if(onlyfirstrow) standata$dokalmanrows <- as.integer(c(1,diff(standata$subject)))
  
  if(!quiet) message('Computing quantities for ', nrow(samples),' samples...')
  if(nrow(samples)==1) cores <- 1
  
  if(cores > 1) {
    if(all(is.na(cl))){
      cl <- makeClusterID(cores)
      on.exit(try(parallel::stopCluster(cl),silent=TRUE),add = TRUE)
    }
    clusterIDexport(cl, c('sm','standata','samples'))
    clusterIDeval(cl,list(
      'require(data.table)',
      'smf <- ctsem::stan_reinitsf(sm,standata)',
      'tparfunc <- function(x){ 
         out <- try(data.table::as.data.table(lapply(1:length(x),function(li){
        unlist(rstan::constrain_pars(smf, upars=samples[x[li],]))
      })))
      if(!"try-error" %in% class(out)) return(out)
  }'))
    
  }
  
  if(cores ==1){
    # browser()
    smf <- stan_reinitsf(sm,standata) 
    tparfunc <- function(x){ 
      out <- try(data.table::as.data.table(lapply(1:length(x),function(li){
        unlist(rstan::constrain_pars(smf, upars=samples[x[li],]))
      })))
      if(!'try-error' %in% class(out)) return(out)
    }
  }
  transformedpars <- try(flexlapplytext(cl, 
    1:nrow(samples),
    'tparfunc',cores=cores))
  nulls <- unlist(lapply(transformedpars,is.null))
  if(any(nulls==FALSE)) transformedpars <- transformedpars[!nulls] else stop('No admissable samples!?')
  if(sum(nulls)>0) message(paste0(sum(nulls)/length(nulls)*100,'% of samples inadmissable'))
  
  if(cores >1) smf <- stan_reinitsf(sm,standata) #needs to be after, weird parallel stuff...
  skel= rstan::constrain_pars(smf, upars=samples[which(!nulls)[1],,drop=FALSE]) 
  transformedpars <- t(data.table::as.data.table(transformedpars))
  
  nasampscount <- nrow(transformedpars)-nrow(samples) 
  
  
  if(nasampscount > 0) {
    message(paste0(nasampscount,' NAs generated during final sampling of ', nrow(samples), '. Biased estimates may result -- consider importance sampling, respecification, or full HMC sampling'))
  }
  if(nasampscount < nrow(samples)){ 
    nresamples <- nrow(samples) - nasampscount
  } else{
    message('All samples contain NAs -- returning anyway')
    nresamples <- nrow(samples) 
  }
  transformedpars=tostanarray(flesh=transformedpars, skeleton = skel)
  
  return(transformedpars)
}



tostanarray <- function(flesh, skeleton){
  skelnames <- names(skeleton)
  skelstruc <- lapply(skeleton,dim)
  count=1
  npars <- ncol(flesh)
  niter=nrow(flesh)
  out <- list()
  for(ni in skelnames){
    if(prod(skelstruc[[ni]])>0){
      if(!is.null(skelstruc[[ni]])){
        out[[ni]] <- array(flesh[,count:(count+prod(skelstruc[[ni]])-1)],dim = c(niter,skelstruc[[ni]]))
        count <- count + prod(skelstruc[[ni]])
      } else {
        out[[ni]] <- array(flesh[,count],dim = c(niter))
        count <- count + 1
      }
    }
  }
  return(out)
}


# parsetup <- parsetup <- function(){
#   cl <- parallel::makeCluster(12,type='PSOCK')
#   parallel::clusterCall(cl,function() 1+1)
# }

makeClusterID <- function(cores){
  # benv <- new.env(parent=globalenv())
  # benv$cores <- cores
  # benv$
  cl <- parallel::makeCluster(spec = cores,type = "PSOCK",useXDR=FALSE,outfile='',user=NULL)
  eval(parse(text="parallel::parLapply(cl = cl,X = 1:cores,function(x) assign('nodeid',x,envir=globalenv()))"))
  return(cl)
}

clusterIDexport <- function(cl, vars){
  # benv <- new.env(parent=globalenv())
  # benv$cl <- cl
  # lookframe <- parent.frame()
  # tmp<-lapply(vars,function(x) benv[[x]] <<- eval(parse(text=x),envir =lookframe))
  # eval(
  parallel::clusterExport(cl,vars,envir = parent.frame())
  # ,envir=globalenv())
}

clusterIDeval <- function(cl,commands){
  # benv <- new.env(parent=globalenv())
  # benv$cl <- cl
  clusterIDexport(cl,'commands')
  # unlist(eval(
  unlist(parallel::clusterEvalQ(cl = cl, 
    lapply(commands,function(x){
      eval(parse(text=x),envir = globalenv())
    })),
    #,envir =globalenv())
    recursive = FALSE)
}


#' Optimize / importance sample a stan or ctStan model.
#'
#' @param standata list object conforming to rstan data standards.
#' @param sm compiled stan model object.
#' @param init vector of unconstrained parameter values, or character string 'random' to initialise with
#' random values very close to zero.
#' @param initsd positive numeric specifying sd of normal distribution governing random sample of init parameters,
#' if init='random' .
#' @param sampleinit either NA, or an niterations * nparams matrix of samples to initialise importance sampling.
#' @param deoptim Do first pass optimization using differential evolution? Slower, but better for cases with multiple
#' minima / difficult optimization.
#' @param stochastic Logical. Use stochastic gradient descent instead of mize (bfgs) optimizer.
#' Still experimental, worth trying for either robustness checks or problematic, high dimensional, nonlinear, problems.
#' @param plot Logical. If TRUE, plot iteration details. Probably slower.
#' @param estonly if TRUE,just return point estimates under $rawest subobject.
#' @param verbose Integer from 0 to 2. Higher values print more information during model fit -- for debugging.
#' @param decontrol List of control parameters for differential evolution step, to pass to \code{DEoptim.control}.
#' @param tol objective tolerance.
#' @param priors logical. If TRUE, a priors integer is set to 1 (TRUE) in the standata object -- only has an effect if 
#' the stan model uses this value. 
#' @param carefulfit Logical. If TRUE, priors are always used for a rough first pass to obtain starting values when priors=FALSE
#' @param subsamplesize value between 0 and 1 representing proportion of subjects to include in first pass fit. 
#' @param cores Number of cpu cores to use, should be at least 2.
#' @param bootstrapUncertainty Logical. If TRUE, subject wise gradient contributions are resampled to estimate the hessian, 
#' for computing standard errors or initializing importance sampling.
#' @param is Logical. Use importance sampling, or just return map estimates?
#' @param isloopsize Number of samples of approximating distribution per iteration of importance sampling.
#' @param finishsamples Number of samples to draw (either from hessian
#' based covariance or posterior distribution) for final results computation.
#' @param finishmultiply Importance sampling stops once available samples reach \code{finishsamples * finishmultiply} , then the final samples are drawn
#' without replacement from this set.
#' @param tdf degrees of freedom of multivariate t distribution. Higher (more normal) generally gives more efficent
#' importance sampling, at risk of truncating tails.
#' @param parsteps ordered list of vectors of integers denoting which parameters should begin fixed
#' at zero, and freed sequentially (by list order). Useful for complex models, e.g. keep all cross couplings fixed to zero 
#' as a first step, free them in second step. 
#' @param parstepsAutoModel if TRUE, determines model structure for the parameters specified in parsteps automatically. If 'group', determines this on a group level first and then a subject level. Primarily for internal ctsem use, see \code{?ctFitAuto}.
#' @param groupFreeThreshold threshold for determining whether a parameter is free in a group level model. If the proportion of subjects with a non-zero parameter is above this threshold, the parameter is considered free. Only used with parstepsAutoModel = 'group'.
#' @param chancethreshold drop iterations of importance sampling where any samples are chancethreshold times more likely to be drawn than expected.
#' @param matsetup subobject of ctStanFit output. If provided, parameter names instead of numbers are output for any problem indications.
#' @param nsubsets number of subsets for stochastic optimizer. Subsets are further split across cores, 
#' but each subjects data remains whole -- processed by one core in one subset.
#' @param stochasticTolAdjust Multiplier for stochastic optimizer tolerance. 
#' @param lproughnesstarget target log posterior roughness for stochastic optimizer (suggest between .05 and .4).
#' @return list containing fit elements
#' @importFrom mize mize
#' @importFrom utils head tail
#' @importFrom Rcpp evalCpp

stanoptimis <- function(standata, sm, init='random',initsd=.01,sampleinit=NA,
  deoptim=FALSE, estonly=FALSE,tol=1e-8,
  decontrol=list(),
  stochastic = TRUE,
  priors=TRUE,carefulfit=TRUE,
  bootstrapUncertainty=FALSE,
  subsamplesize=1,
  parsteps=c(),parstepsAutoModel=FALSE,groupFreeThreshold=.5,
  plot=FALSE,
  is=FALSE, isloopsize=1000, finishsamples=1000, tdf=10,chancethreshold=100,finishmultiply=5,
  lproughnesstarget=.2,
  verbose=0,cores=2,matsetup=NA,nsubsets=10, stochasticTolAdjust=1000){
  
  if(!is.null(standata$verbose)) {
    if(verbose > 1) standata$verbose=as.integer(verbose) else standata$verbose=0L
  }
  standata$priors=as.integer(priors)
  
  if(nsubsets > (standata$nsubjects/10)) nsubsets <- ceiling(standata$nsubjects/10)
  if(nsubsets > (standata$nsubjects/cores)) nsubsets <- max(1,ceiling(standata$nsubjects/cores))
  
  if(is.null(decontrol$steptol)) decontrol$steptol=3
  if(is.null(decontrol$reltol)) decontrol$reltol=1e-2
  if(is.null(decontrol$NP)) decontrol$NP='auto'
  if(is.null(decontrol$CR)) decontrol$CR=.9
  if(is.null(decontrol$trace)) decontrol$trace =ifelse(verbose>0,1,0)
  
  if(is.null(init)) init <- 'random'
  if(init[1] !='random') carefulfit <- FALSE
  
  benv <- new.env(parent=globalenv())
  benv$clctsem <- NA #placeholder for flexsapply usage
  optimfit <- c() #placeholder
  
  notipredsfirstpass <- FALSE
  if(init[1] =='random'){# #remove tipreds for first pass
    notipredsfirstpass <- TRUE
    TIPREDEFFECTsetup <- standata$TIPREDEFFECTsetup
    standata$TIPREDEFFECTsetup[,] <- 0L
    standata$ntipredeffects <- 0L
  }
  
  savesubjectmatrices <- standata$savesubjectmatrices
  standata$savesubjectmatrices <- 0L #reinsert when saving samples
  
  #disabled attempt at stochastic grad with mcmc subject pars
  # if(standata$nindvarying > 0 && standata$intoverpop==0){ #detect subject level pars
  #   stochastic <- TRUE
  #   standata$doonesubject <- 1L
  #   message('Using combined stochastic gradient descent + mcmc')
  #   a1=standata$nparams+standata$nindvarying+
  #     (standata$nindvarying^2-standata$nindvarying)/2
  #   whichmcmcpars <- (a1+1):(a1+standata$nindvarying) #*standata$nsubjects)
  # } else 
  # whichmcmcpars <- NA #leftover necessity
  
  
  smf <- stan_reinitsf(sm,standata)
  npars=rstan::get_num_upars(smf)
  if(cores > 1 & !deoptim) rm(smf)
  
  if(stochastic=='auto' && npars > 50){
    message('> 50 parameters and stochastic="auto" so stochastic gradient descent used -- try disabling if slow!')
    stochastic <- TRUE
  } else if(stochastic=='auto') stochastic <- FALSE
  if(length(parsteps)>0 && !stochastic){
    stochastic=TRUE
    message('Stochastic optimizer used for data driven parameter inclusion') 
  }
  
  if(cores<2 && parstepsAutoModel %in% 'group') stop('parstepsAutoModel = "group" requires cores > 1')
  
  # parsets <- 1
  optimcores <- ifelse(length(unique(standata$subject)) < cores, length(unique(standata$subject)),cores)
  
  betterfit<-TRUE
  bestfit <- -9999999999
  try2 <- FALSE
  
  message('Using ',cores,'/', parallel::detectCores(),' CPU cores')
  
  while(betterfit){ #repeat loop if importance sampling improves on optimized max
    betterfit <- FALSE
    
    # if(standata$TIpredAuto >0) parsteps <- c(parsteps,(npars:(npars-max(standata$TIPREDEFFECTsetup)+1)))
    
    if(all(init %in% 'random')){
      init <- rnorm(npars, 0, initsd)
      if(length(parsteps)>0) init[unlist(parsteps)] <- 0 
    }
    if(all(init == 0)) init <- rep(0,npars)
    if(length(init) != npars) init=c(init[1:min(length(init),npars)],rep(0,abs(npars-length(init))))
    init[is.na(init)] <- 0
    if(notipredsfirstpass && standata$ntipredeffects > 0) init[length(init):(length(init)+1-standata$ntipredeffects)] <- 0
    
    
    if(is.na(sampleinit[1])){
      
      storedPars <- as.numeric(c())#matrix(0,nrow=npars,ncol=0)
      gradstore <- rep(0,npars)
      gradmem <- .9
      storedLp <- c()
      
      optimfinished <- FALSE
      on.exit({
        if(!optimfinished){
          message('Optimization cancelled -- restart from current point by including this argument:')
          message((paste0(c('inits = c(',   paste0(round(storedPars,5),collapse=', '), ')'    ))))
          # message('Return inits? Y/N')
          # if(readline() %in% c('Y','y')) returnValue(storedPars)
        }},add=TRUE)
      
      
      
      singletarget<-function(parm,gradnoise=FALSE) {
        a=Sys.time()
        out<- try(log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
        b=Sys.time()
        if('try-error' %in% class(out) || is.nan(out)) {
          out=-1e100
          attributes(out) <- list(gradient=rep(0,length(parm)))
        }
        storedPars <<- parm
        # storedLp <<- c(storedLp,out[1])
        
        evaltime <- b-a
        if(verbose > 0) print(paste('lp= ',out,' ,    iter time = ',round(evaltime,2)),digits=14)
        if(gradnoise) attributes(out)$gradient <-attributes(out)$gradient *
          exp(rnorm(length(parm),0,1e-3))
        return(out)
      }
      
      if(plot > 0 && .Platform$OS.type=="windows") {
        dev.new(noRStudioGD = TRUE)
        on.exit(expr = {try({dev.off()})},add = TRUE)
      }
      
      if(optimcores==1) {
        target = singletarget #we use this for importance sampling
        eval(parse(text=parlptext))
      }
      
      if(cores > 1){ #for parallelised computation after fitting, if only single subject
        
        assign(x = 'clctsem',
          parallel::makeCluster(spec = cores,type = "PSOCK",useXDR=FALSE,outfile='',user=NULL),
          envir = benv)
        benv$cores=cores
        
        eval(parse(text=
            "parallel::parLapply(cl = clctsem,X = 1:cores,function(x){
         assign('nodeid',x,envir=globalenv())
        })"),envir=benv)
        
        
        # eval(clctsem <- makeClusterID(cores),envir = globalenv())
        on.exit(try({parallel::stopCluster(benv$clctsem)},silent=TRUE),add=TRUE)
        
        if(standata$recompile > 0){
          smfile <- file.path(tempdir(),paste0('ctsem_sm_',ceiling(runif(1,0,100000)),'.rda'))
          save(sm,file=smfile,eval.promises = FALSE,precheck = FALSE)
          on.exit(add = TRUE,expr = {file.remove(smfile)})
        } else smfile <- ''
        
        sapply(c('cores','parlptext','smfile'),function(x){
          assign(x,get(x),envir = benv)
          NULL
        })
        
        eval(parse(text="parallel::clusterExport(clctsem,c('cores','parlptext','smfile'),envir=environment())"),envir = benv)
        # system.time(clusterIDexport(benv$clctsem,c('cores','parlptext','sm')))
        # clusterIDeval(clctsem,c(
        #   'eval(parse(text=parlptext))'
        # ))
        
      }
      
      
      if(optimcores > 1) {
        # #crazy trickery to avoid parallel communication pauses
        # env <- new.env(parent = globalenv(),hash = TRUE)
        # environment(parlp) <- env
        # clusterIDexport(cl = clctsem,vars='parlp')
        iter <-0
        
        target<-function(parm,gradnoise=FALSE) {
          # whichframe <- which(sapply(lapply(sys.frames(),ls),function(x){ 'clctsem' %in% x}))
          a=Sys.time()
          # 
          clusterIDexport(benv$clctsem,'parm')
          if('list' %in% class(parm)){
            out2 <- parallel::clusterEvalQ(cl = benv$clctsem,parlp(parm))
            bestset <- which(unlist(out2)==max(unlist(out2)))
            bestset <- bestset[1]
            parm <- parm[[bestset]]
            out <- out2[[bestset]]
            attributes(out)$bestset <- bestset
          } else {
            # browser()
            out2<-  parallel::clusterEvalQ(cl = benv$clctsem,parlp(parm))
            error <- FALSE
            tmp<-sapply(1:length(out2),function(x) {
              if(!is.null(attributes(out2[[x]])$err)){
                if(!error & length(out2) > 1 && as.logical(verbose)){
                  message('Error on core ', x,' but continuing:')
                  error <<- TRUE
                  message(attributes(out2[[x]])$err)
                }
              }
            })
            # 
            out <- try(sum(unlist(out2)),silent=TRUE)
            coretimes <- sapply(out2,function(x) round(attributes(x)$time,3))
            # attributes(out)$gradient <- try(apply(sapply(out2,function(x) attributes(x)$gradient,simplify='matrix'),1,sum))
            
            for(i in seq_along(out2)){
              if(i==1) attributes(out)$gradient <- attributes(out2[[1]])$gradient
              if(i>1) attributes(out)$gradient <- attributes(out)$gradient+attributes(out2[[i]])$gradient
            }
          }
          b=Sys.time()
          b-a
          
          if('try-error' %in% class(out) || is.nan(out)) {
            out=-1e100
            attributes(out) <- list(gradient=rep(0,length(parm)))
          } 
          
          if(plot > 0 && ( (!stochastic &&!carefulfit && nsubsets ==1))){
            if(out[1] > (-1e99)) storedLp <<- c(storedLp,out[1])
            iter <<- iter+1
            # attributes(out)$gradient <- (1-gradmem)*attributes(out)$gradient + gradmem*gradstore
            # gradstore <<- cbind(gradstore,attributes(out)$gradient)
            # gstore <- log(abs(gradstore))*sign(gradstore)
            g=log(abs(attributes(out)$gradient))*sign(attributes(out)$gradient)
            # if(runif(1,0,1) > .9) {
            if(iter %% plot == 0){
              par(mfrow=c(1,3))
              plot(parm,xlab='param',ylab='par value',col=1:length(parm))
              plot(log(1+tail(-storedLp,500)-min(tail(-storedLp,500))),ylab='target',type='l')
              plot(g,type='p',col=1:length(parm),ylab='gradient',xlab='param')
            }
            if(verbose==0) message(paste('\rlp= ',out,' ,    iter time = ',round(b-a,3), '; core times = ',paste0(coretimes,collapse=', ')),appendLF = FALSE) #if not verbose, print lp when plotting
          }
          storedPars <<- parm
          # storedLp <<- c(storedLp,out[1])
          if(verbose > 0) print(paste('lp= ',out,' ,    iter time = ',round(b-a,3), '; core times = ',paste0(coretimes,collapse=', '))) #if not verbose, print lp when plotting
          if(gradnoise) attributes(out)$gradient <-attributes(out)$gradient * 
            exp( rnorm(length(attributes(out)$gradient),0,1e-3))
          return(out)
        }
        
      } #end multicore setup
      
      
      
      if(deoptim){ #init with DE
        message('Using differential evolution for initialization')
        if(requireNamespace('DEoptim',quietly = TRUE)) {
          if(decontrol$NP=='auto') NP=min(c(40,10*npars)) else NP = decontrol$NP
          
          decontrollist <- c(decontrol,DEoptim::DEoptim.control())
          decontrollist <- decontrollist[unique(names(decontrollist))]
          
          lp2 = function(parm) {
            out<- -target(parm)
            attributes(out)$gradient <- -attributes(out)$gradient
            return(out)
          }
          
          deinit <- matrix(rnorm(npars*NP,0,2),nrow = NP)
          deinit[2,] <- rnorm(npars,0,.0002)
          if(length(init)>1 & try2) { #if better fit found during IS
            deinit[1,] <- unconstrain_pars(smf,init)
            if(NP > 10) deinit[3:9,] =  matrix( rnorm(npars*(7),rep(deinit[1,],each=7),.1), nrow = 7)
          }
          decontrollist$initialpop=deinit
          decontrollist$NP = NP
          
          if(optimcores > 1) parallelStanSetup(cl = benv$clctsem,standata = standata,split=TRUE)
          if(optimcores==1) smf<-stan_reinitsf(sm,standata)
          
          optimfitde <- suppressWarnings(DEoptim::DEoptim(fn = lp2,lower = rep(-1e10, npars), upper=rep(1e10, npars),
            control = decontrollist))
          # init=constrain_pars(object = smf,optimfitde$optim$bestmem)
          init=optimfitde$optim$bestmem
          bestfit = -optimfitde$optim$bestval
        } else stop(paste0('use install.packages(\"DEoptim\") to use deoptim')) #end require deoptim
      }
      
      
      
      
      mizelpg=list( #single core mize functions
        fg=function(pars){
          r=-target(pars)
          r=list(fn=r[1],gr= -attributes(r)$gradient)
          return(r)
        },
        fn=function(x) -target(x),
        gr=function(pars) -attributes(target(pars))$gradient
      )
      
      priorsbak <- standata$priors
      # taylorheun <- standata$taylorheun
      finished <- FALSE
      
      if(carefulfit && !deoptim){ #init using priors
        message('1st pass optimization (carefulfit)...')
        
        # standata$taylorheun <- 1L
        # sdscale <- standata$sdscale
        # standata$sdscale <- sdscale * 1e-6
        # tipredeffectscale <- standata$tipredeffectscale
        # standata$tipredeffectscale <- tipredeffectscale* 1e-5
        if(subsamplesize < 1){
          smlnsub <- min(standata$nsubjects,max(min(30,cores*2),ceiling(standata$nsubjects * subsamplesize)))
          standatasml <- standatact_specificsubjects(standata,
            sample(unique(standata$subject),smlnsub))
        } else standatasml <- standata
        standatasml$priors <- 1L
        standatasml$nsubsets <- as.integer(nsubsets)
        
        
        # smlndat <- min(standatasml$ndatapoints,ceiling(max(standatasml$nsubjects * 10, standatasml$ndatapoints*.5)))
        # standatasml$dokalmanrows[sample(1:standatasml$ndatapoints,smlndat)] <- 0L
        # standatasml$dokalmanrows[match(unique(standatasml$subject),standatasml$subject)] <- 1L #ensure first obs is included for t0var consistency
        if(optimcores > 1) parallelStanSetup(cl = benv$clctsem,standata = standatasml,split=TRUE,nsubsets = nsubsets)
        if(optimcores==1) smf<-stan_reinitsf(sm,standatasml)
        
        
        if(stochastic) {
          
          optimfit <- try(sgd(init, fitfunc = function(x) target(x),
            # parsets=parsets,
            nsubsets = nsubsets,
            whichignore = unlist(parsteps),nconvergeiter = 30,
            plot=plot, 
            itertol=tol*1000*stochasticTolAdjust,
            worsecountconverge = 20,
            maxiter=ifelse(standata$ntipred > 0 && notipredsfirstpass, 500,5000)))
          
        }
        
        if(!stochastic || 'try-error' %in% class(optimfit)) {
          optimfit <- mize(init, fg=mizelpg, max_iter=99999,
            method="L-BFGS",memory=100,
            line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
            abs_tol=1e-2,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
          optimfit$value = optimfit$f
        }
      } #end carefulfit  
      
      carefulfit <- FALSE
      iter <- 0
      storedLp <- c()
      # if(length(parsteps) <1 && (standata$ntipred ==0 || notipredsfirstpass ==FALSE)) finished <- TRUE
      
      
      
      
      
      
      if(standata$ntipred > 0 && notipredsfirstpass){
        message('Including TI predictors...')
        
        # standata$tipredeffectscale <- tipredeffectscale
        standata$dokalmanrows[] <- 1L
        # init = optimfit$par
        optimbase = init#[-(length(optimfit$par):(length(optimfit$par)-max(TIPREDEFFECTsetup)+1))]
        finished <- FALSE
        found <- 0
        tia=standata$TIPREDEFFECTsetup
        while(!finished){
          if(found==0){
            if(!priors){ #then we need to reinitialise model
              if(optimcores > 1) parallelStanSetup(cl = benv$clctsem,standata = standata,split=TRUE,nsubsets = nsubsets)
              if(optimcores==1) smf<-stan_reinitsf(sm,standata)
            }
          }
          if(!is.null(standata$TIpredAuto) && standata$TIpredAuto){
            message ('Looking for tipred effects...')
            oldtia=tia
            fit <- list(stanfit=list(rawest=init),standata=standata,stanmodel=sm)
            tia=ctTIauto(fit)
            tia[tia > .05] = 0L
            a <- 0
            for(ci in 1:ncol(tia)){
              for(ri in 1:nrow(tia)){
                if(tia[ri,ci] +oldtia[ri,ci] > 0){
                  a <- a+1
                  tia[ri,ci] <- as.integer(a)
                }
              }
            }
            
            
            if(max(tia) > found){
              tiinit=rep(0,max(tia))
              oldtiinit = tail(init,found)[oldtia[oldtia>0]]
              tiinit[tia[tia>0&oldtia>0]] <- oldtiinit
              standata$TIPREDEFFECTsetup <- array(as.integer(tia),dim=dim(tia))
              found <- max(tia)
              message('Found ',max(tia),' viable TIpred effects')
              standata$ntipredeffects <- as.integer(max(tia))
              init = head(init,length(optimbase))
              init = c(init,tiinit)
            } else {
              finished <- TRUE
              message('No further predictors found, finishing optimization...')
              # standata$taylorheun <- as.integer(taylorheun)
            }
          }
          # if(!finished){
          
          if(!standata$TIpredAuto){
            finished <- TRUE
            init = c(optimbase,rep(0,max(TIPREDEFFECTsetup)))
            standata$TIPREDEFFECTsetup[,] <- TIPREDEFFECTsetup
            standata$ntipredeffects <- as.integer(max(TIPREDEFFECTsetup))
          }
          
          standata$nsubsets <- as.integer(nsubsets)
          if(optimcores > 1) parallelStanSetup(cl = benv$clctsem,standata = standata,split=TRUE,nsubsets = nsubsets)
          if(optimcores==1) smf<-stan_reinitsf(sm,standata)
          
          if(stochastic || nsubsets > 1) optimfit <- try(sgd(init, fitfunc = target,
            lproughnesstarget=lproughnesstarget,
            # parsets=parsets,
            nsubsets = nsubsets,
            whichignore = unlist(parsteps),
            plot=plot,
            maxiter=5000,
            itertol=tol*1000*stochasticTolAdjust,worsecountconverge = 20))
          
          
          if((!stochastic && nsubsets ==1) || 'try-error' %in% class(optimfit)) {
            optimfit <- mize(init, fg=mizelpg, max_iter=99999,
              method="L-BFGS",memory=100,
              line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
              abs_tol=tol*ifelse(finished,1,10000),grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
            optimfit$value = -optimfit$f
          }
          
          
          
        }
        if(length(parsteps)>0) init[-unlist(parsteps)] = optimfit$par else init=optimfit$par
      } #end ti pred auto total loop
      
      
      finished <- TRUE
      standata$nsubsets <- nsubsets <- 1L

      if(optimcores > 1) parallelStanSetup(cl = benv$clctsem,standata = standata,split=TRUE)
      if(optimcores==1) smf<-stan_reinitsf(sm,standata)
      
      
      ##parameter stepwise / selection
      if(length(parsteps) > 0){
        if(parstepsAutoModel %in% FALSE){
          message('Freeing parameters...')
          finished <- FALSE
          while(!finished && length(parsteps)>0){
            if(length(parsteps)>1) parsteps <- parsteps[-1] else parsteps <- c()#parsteps[pstat> pcheck]
            
            optimfit <- sgd(init, fitfunc = target,
              lproughnesstarget=lproughnesstarget,
              # parsets=parsets,
              nsubsets = nsubsets,
              itertol = tol*1000*stochasticTolAdjust,
              maxiter=5000,
              whichignore = unlist(parsteps),
              plot=plot,worsecountconverge = 20)
            
            
            if(length(parsteps)>0) init[-unlist(parsteps)] = optimfit$par else{
              finished <- TRUE
              init = optimfit$par
            }
          }
        }
        if(parstepsAutoModel %in% TRUE){
          # -----------------------------
          # Assume the following are available:
          #   - parFreeList: a list of length 2
          #         parFreeList[[1]]: vector of indices for basic parameters (always estimated)
          #         parFreeList[[2]]: vector of candidate parameter indices that are initially fixed
          #   - init: the current parameter vector (from the previous optimization stage)
          #   - target: a function that returns the log probability, with an attribute "gradient"
          #   - jac: a function to compute a finite-difference approximation of a parameter's Hessian (diagonal element)
          #   - sgd: your stochastic optimizer that accepts an argument 'whichignore' (the parameters to keep fixed)
          #   - tol, stochasticTolAdjust, nsubsets, plot: parameters as in your code.
          #
          # Set a Wald threshold:
          wald_threshold <- 1.96  
          
          jacPars <- function(pars, step = 1e-3, whichpars) {
            # Initialize a vector to store the Hessian diagonal estimates for the specified parameters.
            hess_diag <- numeric(length(whichpars))
            # Loop over each requested parameter index.
            for (i in seq_along(whichpars)) {
              idx <- whichpars[i]
              # Create perturbed parameter vectors: one for the forward difference and one for the backward difference.
              pars_forward <- pars
              pars_backward <- pars
              pars_forward[idx] <- pars_forward[idx] + step
              pars_backward[idx] <- pars_backward[idx] - step
              # Evaluate the target function at both perturbed vectors.
              # It is assumed the 'target' function returns an object with the gradient as an attribute "gradient".
              forward_val <- target(pars_forward)
              backward_val <- target(pars_backward)
              # Extract the gradient for the current parameter.
              grad_forward <- attributes(forward_val)$gradient[idx]
              grad_backward <- attributes(backward_val)$gradient[idx]
              # Estimate the second derivative via the central difference formula.
              hess_diag[i] <- (grad_forward - grad_backward) / (2 * step)
            }
            return(hess_diag)
          }
          
          # Start with all candidate parameters fixed:
          currentFixed <- parsteps[[1]]
          # The permanently free parameter indices are:
          freePars <- (1:length(init))[!parsteps[[1]] %in% (1:length(init))]
          # Track which candidate parameters have been freed (initially, none)
          freed_candidates <- c()
          
          continueFreeing <- TRUE
          improvement_threshold <- 1.96
          while (continueFreeing && length(currentFixed) > 0) {
            
            # Evaluate the current log probability and obtain the gradient.
            obj_val <- target(init)
            grad_vec <- attributes(obj_val)$gradient
            
            # For each candidate parameter currently fixed, compute expected improvement.
            improvement_est <- numeric(length(currentFixed))
            for (i in seq_along(currentFixed)) {
              idx <- currentFixed[i]
              
              # Use our simple jacPars to compute the approximate second derivative for this parameter.
              hess_est <- jacPars(init, step = 1e-6, whichpars = idx)
              
              # For stability, if the second derivative is NA or non-negative (which is
              # unexpected at a maximum), force a small negative value.
              if (is.na(hess_est) || hess_est >= 0) {
                hess_val <- -1e-6
              } else {
                hess_val <- hess_est
              }
              
              # Estimate expected improvement: 0.5 * (grad^2 / |H_ii|)
              improvement_est[i] <- 0.5 * (grad_vec[idx]^2 / abs(hess_val))
            }
            
            # Find the candidate with the highest expected improvement.
            best_candidate_index <- which.max(improvement_est)
            best_improvement <- improvement_est[best_candidate_index]
            
            if (best_improvement >= improvement_threshold) {
              # browser()
              best_param <- currentFixed[best_candidate_index]
              ms=data.frame(standata$matsetup)
              bestpar_ms <- ms[ms$param == best_param,,drop=FALSE][1,]
              message(paste0("Freeing parameter ",names(sort(ctStanMatricesList()$all))[bestpar_ms$matrix],'[',
                bestpar_ms$row,',',bestpar_ms$col,'] with expected improvement ',round(best_improvement,3)))
              
              # Mark the best candidate as freed.
              freed_candidates <- c(freed_candidates, best_param)
              # Remove it from the list of currently fixed candidate parameters.
              currentFixed <- currentFixed[-best_candidate_index]
              
              # Create the list of indices to ignore in the next optimization step.
              # If no candidates remain, use an empty vector.
              ignore_indices <- if (length(currentFixed) > 0) currentFixed else integer(0)
              
              
              # Now, re-run the optimization so that the newly freed parameter is estimated.
              # Note that 'freePars' (the basic parameters) and the already freed candidates remain free.
              optimfit <- sgd(init,
                fitfunc = target,
                itertol = tol * stochasticTolAdjust,
                maxiter = 5000,
                whichignore = ignore_indices,
                plot = plot,
                worsecountconverge = 20)
              
              # Update the current parameter vector with the new (optimized) estimates.
              init[-ignore_indices] <- optimfit$par
            } else {
              message("No candidate parameter meets the improvement threshold. Terminating freeing sequence.")
              continueFreeing <- FALSE
              parsteps <- currentFixed
            }
          }
          
          # At this point, all parameters in `freePars` (the basic ones) and those in `freed_candidates`
          # are free (and have been re–optimized), while those remaining in `currentFixed` are kept fixed.
          # You can now proceed with the rest of your model estimation using 'init' as the final parameter estimate.
          parsteps <- currentFixed
        }
      }
      
      if (parstepsAutoModel %in% 'group') {
        # --------------------------------------------------
        # Group‐level stepwise freeing, with per‐subject re‐fit and init averaging
        
        # thresholds
        improvement_threshold <- 1.96   # per‐subject ΔLL must exceed this
        
        # initial fixed candidates and subjects
        currentFixed <- parsteps[[1]]
        subject_ids  <- unique(standata$subject)
        
        
        parallel::clusterExport(benv$clctsem, c(
          "tol", "stochasticTolAdjust",  "standatact_specificsubjects","subject_ids"
        ), envir = environment())
        
        continueFreeing <- TRUE
        subj_init <- matrix(rep(init, length(subject_ids)), nrow = length(subject_ids), byrow = TRUE)
        while (continueFreeing && length(currentFixed) > 0) {
          # export updated init and fixed set to workers
          parallel::clusterExport(benv$clctsem, c("currentFixed",'subj_init'), envir = environment())
          # 1) fit each subject with currentFixed held fixed
          subj_pars <- parallel::parLapply(benv$clctsem, subject_ids, function(sid) {
            sd    <- standatact_specificsubjects(standata, sid)
            smf   <- stan_reinitsf(sm, sd)
            parlp <- function(parm){
              out <- try(rstan::log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
              if("try-error" %in% class(out) || any(is.nan(attributes(out)$gradient))) {
                outerr <- out
                out <- -1e100
                attributes(out)$gradient <- rep(NaN, length(parm))
                attributes(out)$err <- outerr
              }
              if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
              return(out)
            }
            fit   <- sgd(
              init=subj_init[sid,],
              fitfunc     = parlp,
              itertol     = tol * stochasticTolAdjust,
              maxiter     = 5000,
              whichignore = currentFixed,
              worsecountconverge = 20
            )
            as.numeric(fit$par)
          })
          subj_mat <- do.call(rbind, subj_pars)
          subj_init[,-currentFixed] <- subj_mat # update init for next iteration
          
          # 2) average into init
          init[-currentFixed] <- colMeans(subj_mat)
          
          # 3) compute per‐subject expected ΔLL for each candidate
          subj_imp <- parallel::parLapply(benv$clctsem,seq_along(subject_ids), function(i) {
            sid  <- subject_ids[i]
            pvec <- init
            pvec[-currentFixed] <- subj_mat[i, ]
            sd   <- standatact_specificsubjects(standata, sid)
            smf  <- stan_reinitsf(sm, sd)
            tgt  <- function(p) log_prob(smf, upars = p, adjust_transform = TRUE, gradient = TRUE)
            grad <- attributes(tgt(pvec))$gradient
            # helper: finite‐difference Hessian diagonal
            jacPars <- function(pars, step = 1e-3, whichpars) {
              sapply(whichpars, function(idx) {
                pf <- pars; pb <- pars
                pf[idx] <- pf[idx] + step
                pb[idx] <- pb[idx] - step
                gf <- attributes(tgt(pf))$gradient[idx]
                gb <- attributes(tgt(pb))$gradient[idx]
                (gf - gb) / (2 * step)
              })
            }
            sapply(currentFixed, function(idx) {
              h <- jacPars(pvec, step = 1e-6, whichpars = idx)
              if (is.na(h) || h >= 0) h <- -1e-6
              0.5 * (grad[idx]^2 / abs(h))
            })
          })
          # coerce to matrix if needed
          imp_vecs <- lapply(subj_imp, as.numeric)
          imp_mat  <- do.call(rbind, imp_vecs)
          if (is.null(dim(imp_mat))) {
            imp_mat <- matrix(imp_mat, nrow = length(imp_vecs), byrow = TRUE)
          }
          
          # 4) identify group‐level candidates
          prop_above  <- colMeans(imp_mat >= improvement_threshold)
          group_cands <- currentFixed[prop_above > groupFreeThreshold]
          
          if (length(group_cands) == 0) {
            message("No group‐level parameters exceed threshold; stopping.")
            break
          }
          
          # pick group candidate with highest mean ΔLL
          means      <- colMeans(imp_mat[, currentFixed %in% group_cands, drop = FALSE])
          best_param <- group_cands[which.max(means)]
          
          message(sprintf(
            "Freeing parameter %d (%.0f%% subjects ΔLL ≥ %.2f)",
            best_param,
            100 * prop_above[currentFixed == best_param],
            improvement_threshold
          ))
          
          # update fixed set only
          currentFixed <- setdiff(currentFixed, best_param)
        }
        groupFixed <- currentFixed
        parallel::clusterExport(benv$clctsem, c("groupFixed", "subj_init"), envir = environment())
  
        subj_res <- parallel::parLapply(benv$clctsem, seq_along(subject_ids), function(i) {
        # for(i in 1:length(subject_ids)){
        # lapply(seq_along(subject_ids), function(i) {
          sid    <- subject_ids[i]
          sd     <- standatact_specificsubjects(standata, sid)
          smf    <- stan_reinitsf(sm, sd)
          p_i    <- subj_init[sid,]
          free_i <- setdiff(seq_along(init), groupFixed)
          subjFixed <- groupFixed
          freed_i   <- integer(0)
          parlp <- function(parm){
            out <- try(rstan::log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
            if("try-error" %in% class(out) || any(is.nan(attributes(out)$gradient))) {
              outerr <- out
              out <- -1e100
              attributes(out)$gradient <- rep(NaN, length(parm))
              attributes(out)$err <- outerr
            }
            if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
            return(out)
          }
          repeat {
            lpinit <- parlp(p_i)
            grad <- attributes(lpinit)$gradient
            impr <- sapply(subjFixed, function(idx) {
              pf <- p_i; pb <- p_i
              pf[idx] <- pf[idx] + 1e-6
              pb[idx] <- pb[idx] - 1e-6
              gf <- attributes(parlp(pf))$gradient[idx]
              gb <- attributes(parlp(pb))$gradient[idx]
              h  <- (gf - gb) / (2 * 1e-6)
              if (is.na(h) || h >= 0) h <- -1e-6
              0.5 * (grad[idx]^2 / abs(h))
            })
            if (all(impr < improvement_threshold)) break
            best_idx  <- subjFixed[which.max(impr)]
            freed_i   <- c(freed_i, best_idx)
            subjFixed <- setdiff(subjFixed, freed_i)
            if(length(subjFixed) == 0) break
            fit <- sgd(
              p_i+rnorm(length(p_i),0,.01), #init away from old max
              fitfunc = parlp,
              itertol     = tol * stochasticTolAdjust,
              maxiter     = 5000,
              whichignore = subjFixed,
              worsecountconverge = 20
            )
            p_i[sort(c(free_i,freed_i))] <- fit$par
          }
          list(par = p_i, freed = freed_i)
        })
        
        # build output matrices
        subjPars  <- t(sapply(subj_res, `[[`, "par"))
        subjFreed <- t(sapply(subj_res, function(x) groupFixed %in% x$freed))
        rownames(subjFreed) <- subject_ids
        colnames(subjFreed) <- as.character(groupFixed)
        
        # finalize
        parsteps  <- groupFixed
        optimfit  <- list(
          par       = init[-parsteps],
          subjPars  = subjPars,
          subjFreed = subjFreed
        )
      }
      
      
      
      
      if(parstepsAutoModel %in% FALSE){
        if(stochastic){
          message('Optimizing...')
          optimfit <- try(sgd(init, fitfunc = target,
            lproughnesstarget=lproughnesstarget,
            # parsets=parsets,
            nsubsets = 1,
            itertol = tol*stochasticTolAdjust,
            parrangetol=tol*100,
            maxiter=5000,
            whichignore = unlist(parsteps),
            plot=plot))
        }
        
        if(!'try-error' %in% class(optimfit) & !'NULL' %in% class(optimfit)){
          if(length(parsteps)>0) init[-unlist(parsteps)] = optimfit$par else init=optimfit$par
        }
        
        #use bfgs to double check stochastic fit (or just use bfgs if requested)... 
        message('Finishing optimization...')
        optimfit <- mize(init, fg=mizelpg, max_iter=99999,
          method="L-BFGS",memory=100,
          line_search='Schmidt',c1=1e-4,c2=.9,step0='schmidt',ls_max_fn=999,
          abs_tol=NULL,grad_tol=NULL,
          rel_tol=tol,
          step_tol=NULL,ginf_tol=NULL)
        # if(verbose==0 && as.logical(plot)) message('')
        optimfit$value = -optimfit$f
        init = optimfit$par
      } #end if not auto model parsteps
      
      
      bestfit <-optimfit$value
      est2=init #because init contains the fixed values #unconstrain_pars(smf, est1)
      if(length(parsteps)>0) est2[-parsteps] = optimfit$par else est2=optimfit$par
      
      # if(standata$nindvarying > 0 && standata$intoverpop==0){ #recompose into single model
      #   # standata$doonesubject <- 0L
      #   a1=standata$nparams+standata$nindvarying+
      #     (standata$nindvarying^2-standata$nindvarying)/2
      #   whichmcmcpars <- (a1+1):(a1+standata$nindvarying*standata$nsubjects)
      #   est3=est2[-whichmcmcpars]
      #   est3=c(est3,(optimfit$mcmcpars))
      #   est2=est3
      #   if(standata$ntipredeffects > 0) est3 <- c(est3,est2[(a1+1):length(a1)])
      # }
      npars = length(est2)
    }
    
    if(!estonly){
      # if(cores > 1){
      #   suppressWarnings(rm(smf))
      #   # parallelStanSetup(cl = clctsem,standata = standata,split=TRUE)
      #   hesscl <- NA
      # }
      # if(cores==1){
      #   hesscl <- NA
      #   # smf<-stan_reinitsf(sm,standata)
      # }
      
      # standata$nsubsets <- 1L
      # if(standata$nsubjects > cores){
      #   hesscl <- NA
      # } else {
      #   hesscl <- benv$clctsem
      #   if(cores > 1) {
      #     suppressWarnings(rm(smf))
      #     parallelStanSetup(cl = benv$clctsem,standata = standata,split=TRUE,nsubsets = 1)#,split=TRUE)
      #   }
      #   if(cores==1) smf<-stan_reinitsf(sm,standata)
      # }
      
      
      
      
      
      if(is.na(sampleinit[1])){
        
        message('Estimating Hessian',appendLF = bootstrapUncertainty)
        plot=FALSE
        
        
        if(length(parsteps)>0) grinit= est2[-parsteps] else grinit = est2
        
        if(bootstrapUncertainty %in% 'auto'){
          if(standata$nsubjects < 5){
            bootstrapUncertainty <- FALSE
            message('Too few subjects (<5) for bootstrap uncertainty, classical hessian used')
          } else bootstrapUncertainty=TRUE
        }
        
        if(bootstrapUncertainty){
          # a=Sys.time()
          scores <- scorecalc(standata = standata,est = grinit,stanmodel = sm,
            subjectsonly = standata$nsubjects > 5,returnsubjectlist = F,cores=cores)
          
          # message(paste0('score time = '))
          # print(Sys.time()-a)
          # b=Sys.time()
          num_bootstrap_samples <- max(c(finishsamples,1000))
          alpha_max = 100 # Maximum bootstrap sample size factor
          alpha_min = 1 # Minimum bootstrap sample size factor
          n_threshold=1000 # Threshold n for alpha correction
          alpha <- alpha_max - (alpha_max - alpha_min) * (min(1000,nrow(scores))  / n_threshold)  # Bootstrap sample size factor
          num_bootstrap_samples  # Total number of bootstrap samples
          n <- nrow(scores)  # Number of observations
          p <- ncol(scores)  # Number of parameters
          
          # Step 1: Create a bootstrap resampling matrix
          resample_matrix <- matrix(sample(1:n, size = round(alpha * n) * num_bootstrap_samples, replace = TRUE),
            nrow = num_bootstrap_samples, ncol = round(alpha * n))
          
          # Step 2: Generate random weights for smoothing
          weights <- matrix(runif(length(resample_matrix), min = 0.1, max = 2),
            nrow = num_bootstrap_samples)
          
          # Step 3: Aggregate gradients using matrix multiplication
          gradsamples <- matrix(0, nrow = num_bootstrap_samples, ncol = p)  # Initialize gradsamples
          
          # Compute the weighted sum of gradients for each bootstrap sample
          for (i in 1:num_bootstrap_samples) {
            gradsamples[i, ] <- colSums(scores[resample_matrix[i, ], , drop = FALSE] * weights[i, ])
          }
          
          # Step 4: Trim outliers
          trim_percent <- 0.05  # Proportion of outliers to trim
          if (trim_percent > 0) {
            lower <- apply(gradsamples, 2, quantile, probs = trim_percent, na.rm = TRUE)
            upper <- apply(gradsamples, 2, quantile, probs = 1 - trim_percent, na.rm = TRUE)
            gradsamples <- pmax(pmin(gradsamples, upper), lower)
          }
          
          # Step 5: Compute the Hessian with alpha correction
          hess <- -cov(gradsamples) / alpha
          # message(paste0('boothess time = '))
          # print(Sys.time()-b)
        }
        
        
        
        
        
        
        if(!bootstrapUncertainty){
          jac<-function(pars,step=1e-3,whichpars='all',
            lpdifmin=1e-5,lpdifmax=5, cl=NA,verbose=1,directions=c(-1,1),parsteps=c()){
            if('all' %in% whichpars) whichpars <- 1:length(pars)
            base <- optimfit$value
            
            
            hessout <- sapply( whichpars, function(i){
              
              
              # for(i in whichpars){
              message(paste0("\rEstimating Hessian, par ",i,',', 
                as.integer(i/length(pars)*50+ifelse(directions[1]==1,0,50)),
                '%'),appendLF = FALSE)
              if(verbose) message('### Par ',i,'###')
              stepsize = step
              uppars<-rep(0,length(pars))
              uppars[i]<-1
              accepted <- FALSE
              count <- 0
              lp <- list()
              steplist <- list()
              for(di in 1:length(directions)){
                count <- 0
                accepted <- FALSE
                stepchange = 0
                stepchangemultiplier = 1
                while(!accepted && (count==0 || 
                    ( 
                      (count < 30 && any(is.na(attributes(lp[[di]])$gradient))) || #if NA gradient, try for 30 attempts
                        (count < 15 && all(!is.na(attributes(lp[[di]])$gradient))))#if gradient ok, stop after 15
                )){ 
                  # if(count>8) stepsize=stepsize*-1 #is this good?
                  stepchangemultiplier <- max(stepchangemultiplier,.11)
                  count <- count + 1
                  lp[[di]] <-  suppressMessages(suppressWarnings(target(pars+uppars*stepsize*directions[di])))
                  accepted <- !'try-error' %in% class(lp[[di]]) && all(!is.na(attributes(lp[[di]])$gradient))
                  if(accepted){
                    lpdiff <- base[1] - lp[[di]][1]
                    # if(lpdiff > 1e100) 
                    if(lpdiff < lpdifmin) {
                      if(verbose) message('Increasing step')
                      if(stepchange == -1) stepchangemultiplier = stepchangemultiplier*.5
                      stepchange <- 1
                      stepsize <- stepsize*(1-stepchangemultiplier)+ (stepsize*10)*stepchangemultiplier
                    }
                    if(lpdiff > lpdifmax){
                      if(verbose) message('Decreasing step')
                      
                      
                      if(stepchange == 1) stepchangemultiplier = stepchangemultiplier * .5
                      stepchange <- -1
                      stepsize <- stepsize*(1-stepchangemultiplier)+ (stepsize*.1)*stepchangemultiplier
                    }
                    if(lpdiff > lpdifmin && lpdiff < lpdifmax && lpdiff > 0) accepted <- TRUE else accepted <- FALSE
                    if(lpdiff < 0){
                      base <- lp[[di]]
                      if(verbose) message('Better log probability found during Hessian estimation...')
                      accepted <- FALSE
                      stepchangemultiplier <- 1
                      stepchange=0
                      count <- 0
                      di <- 1
                    }
                  } else stepsize <- stepsize * 1e-3
                  # if(count > 1) print(paste0(count,'___',stepsize,'___',lpdiff))
                }
                if(stepsize < step) step <<- step *.1
                if(stepsize > step) step <<- step *10
                steplist[[di]] <- stepsize
              }
              
              grad<- attributes(lp[[1]])$gradient / steplist[[di]] * directions[di]
              if(any(is.na(grad))){
                warning('NA gradient encountered at param ',i,immediate. =TRUE)
              }
              if(length(directions) > 1) grad <- (grad + attributes(lp[[2]])$gradient / (steplist[[di]]*-1))/2
              return(grad)
            }
            ) #end sapply
            
            out=(hessout+t(hessout))/2
            return(out)
          }
          
          hess1 <- jac(pars = grinit,parsteps=parsteps,
            step = 1e-3,cl=NA,verbose=verbose,directions=1)
          hess2 <- jac(pars = grinit,parsteps=parsteps,#fgfunc = fgfunc,
            step = 1e-3,cl=NA,verbose=verbose,directions=-1)
          message('') #to create new line due to overwriting progress bar
          
          
          
          probpars <- c()
          onesided <- c()
          
          hess <- hess1
          hess[is.na(hess)] <- 0 #set hess1 NA's to 0
          hess[!is.na(hess2)] <- hess[!is.na(hess2)] + hess2[!is.na(hess2)] #add hess2 non NA's
          hess[!is.na(hess1) & !is.na(hess2)] <- hess[!is.na(hess1) & !is.na(hess2)] /2 #divide items where both hess1 and 2 used by 2
          hess[is.na(hess1) & is.na(hess2)] <- NA #set NA when both hess1 and 2 NA
          
          if(any(is.na(c(diag(hess1),diag(hess2))))){
            if(any(is.na(hess))) message ('Problems computing Hessian...')
            onesided <- which(sum(is.na(diag(hess1) & is.na(diag(hess2)))) %in% 1)
          }
          
          hess <- (t(hess)+hess)/2
          
          
          if(length(c(probpars,onesided)) > 0){
            if('data.frame' %in% class(matsetup) || !all(is.na(matsetup[1]))){
              ms=matsetup
              ms=ms[ms$param > 0 & ms$when == 0,]
              ms=ms[!duplicated(ms$param),]
            }
            if(length(onesided) > 0){
              onesided=paste0(ms$parname[ms$param %in% onesided],collapse=', ')
              message ('One sided Hessian used for params: ', onesided)
            }
            if(length(probpars) > 0){
              probpars=paste0(ms$parname[ms$param %in% c(probpars)],collapse=', ')
              message('***These params "may" be not identified: ', probpars)
            }
          }
        } #end classical hessian
        
        # cholcov = try(suppressWarnings(t(chol(solve(-hess)))),silent = TRUE)
        
        # if('try-error' %in% class(cholcov) && !is) message('Approximate hessian used for std error estimation.')
        
        mcov=try(solve(-hess),silent=TRUE)
        if('try-error' %in% class(mcov)){
          mcov=MASS::ginv(-hess) 
          warning('***Generalized inverse required for Hessian inversion -- interpret standard errors with caution. Consider simplification, priors, or alternative uncertainty estimators',call. = FALSE,immediate. = TRUE)
          probpars=which(diag(hess) > -1e-6)
        }
        
        
        
        
        
        # if(robust){
        #   message('Getting scores for robust std errors...')
        #   mcovnonrobust <-mcov
        #   fit <- list(standata=standata,stanmodel=sm) #this is hacky...
        #   fit$stanfit$rawest=est2
        #   
        #   score=as.matrix(scorecalc(fit,subjectsonly = F,returnsubjectlist = F,cores=cores))
        #   mcov=tcrossprod(mcov,MASS::ginv(crossprod(score)))
        # }
        # 
        
        # mcov <- mcov+diag(1e-20,nrow(mcov))
        mcovtmp=try({as.matrix(Matrix::nearPD(mcov,conv.norm.type = 'F')$mat)})
        if(any(class(mcovtmp) %in% 'try-error')) stop('Hessian could not be computed')
        mcov <- diag(1e-10,npars)
        if(length(parsteps)>0) mcov[-parsteps,-parsteps] <- mcovtmp else mcov <- mcovtmp
        mchol = t(chol(mcov))
        
        
      }
      
      
      if(!is.na(sampleinit[1])){
        mcov = cov(sampleinit)*1.5+diag(1e-6,ncol(sampleinit))
        est2 = apply(sampleinit,2,mean)
        bestfit = 9e100
        optimfit <- suppressWarnings(list(par=rstan::sampling(sm,standata,iter=2,control=list(max_treedepth=1),chains=1,show_messages = FALSE,refresh=0)@inits[[1]]))
      }
      mcovl<-list()
      mcovl[[1]]=mcov
      delta=list()
      delta[[1]]=est2
      samples <-matrix(NA)
      resamples <- c()
      prop_dens <-c()
      target_dens<-c()
      sample_prob<-c()
      ess <- 0
      qdiag<-0
      
      if(!is) {
        nresamples = finishsamples
        resamples <- matrix(unlist(lapply(1:nresamples,function(x){
          delta[[1]] + (mchol) %*% t(matrix(rnorm(length(delta[[1]])),nrow=1))
        } )),byrow=TRUE,ncol=length(delta[[1]]))
        # resamples<- rbind(resamples,gensigstates(delta[[1]],mchol))
        
      }
      message('')
      if(is){
        message('Importance sampling...')
        
        log_sum_exp <- function(x) {
          xmax <- which.max(x)
          log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
        }
        
        #configure each node with full dataset for adaptive sampling
        # browser()
        if(cores > 1)   parallelStanSetup(cl=benv$clctsem,standata,split=FALSE)
        targetsamples <- finishsamples * finishmultiply
        # message('Adaptive importance sampling, loop:')
        j <- 0
        while(nrow(samples) < targetsamples){
          j<- j+1
          if(j==1){
            samples <- newsamples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)
          } else {
            delta[[j]]=colMeans(resamples)
            mcovl[[j]] = as.matrix(Matrix::nearPD(cov(resamples))$mat) #+diag(1e-12,ncol(samples))
            
            for(iteri in 1:j){
              if(iteri==1) mcov=mcovl[[j]] else mcov = mcov*.5+mcovl[[j]]*.5 #smoothed covariance estimator
              if(iteri==1) mu=delta[[j]] else mu = mu*.5+delta[[j]]*.5 #smoothed means estimator
            }
            newsamples <- mvtnorm::rmvt(isloopsize, delta = mu, sigma = mcov,   df = tdf)
            samples <- rbind(samples, newsamples)
          }
          
          prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize), delta[[j]], mcovl[[j]], df = tdf,log = TRUE)
          
          if(cores > 1) parallel::clusterExport(benv$clctsem,c('samples'),envir = environment())
          
          # if(cores==1) parlp <- function(parm){ #remove this duplication somehow
          #   out <- try(log_prob(smf,upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
          #   if('try-error' %in% class(out)) {
          #     out[1] <- -1e100
          #     attributes(out)$gradient <- rep(NaN, length(parm))
          #   }
          #   if(is.null(attributes(out)$gradient)) attributes(out)$gradient <- rep(NaN, length(parm))
          #   attributes(out)$gradient[is.nan(attributes(out)$gradient)] <- 
          #     rnorm(length(attributes(out)$gradient[is.nan(attributes(out)$gradient)]),0,100)
          #   return(out)
          # }
          
          # 
          # clusterIDexport(clctsem, c('isloopsize'))
          target_dens[[j]] <- unlist(flexlapplytext(cl = benv$clctsem, 
            X = 1:isloopsize, 
            fn = "function(x){parlp(samples[x,])}",cores=cores))
          
          
          
          target_dens[[j]][is.na(target_dens[[j]])] <- -1e200
          if(all(target_dens[[j]] < -1e100)) stop('Could not sample from optimum! Try reparamaterizing?')
          
          ### this was used to restart optimization in case a better logprob was found during sampling, but might have caused cran problems with break statement.
          # if(any(target_dens[[j]] > bestfit)){
          #   oldfit <- bestfit
          #   try2 <- TRUE
          #   bestfit<-max(target_dens[[j]],na.rm=TRUE)
          #   betterfit<-TRUE
          #   init = samples[which(unlist(target_dens) == bestfit),]
          #   message('Improved fit found - ', bestfit,' vs ', oldfit,' - restarting optimization')
          #   break
          # }
          
          targetvec <- unlist(target_dens)
          
          target_dens2 <- target_dens[[j]] #-max(target_dens[[j]],na.rm=TRUE) + max(prop_dens) #adjustment to get in decent range, doesnt change to prob
          target_dens2[!is.finite(target_dens[[j]])] <- -1e30
          weighted_dens <- target_dens2 - prop_dens
          # psis_dens <- psis(matrix(target_dens2,ncol=length(target_dens2)),r_eff=NA)
          # sample_prob <- weights(psis_dens,normalize = TRUE,log=FALSE)
          # plot(target_dens2,prop_dens)
          
          newsampleprob <- exp((weighted_dens - log_sum_exp(weighted_dens)))
          counter <- 1
          isfinished <- FALSE
          
          # while(counter < 30 && length(unique(sample(x=1:length(newsampleprob),size=length(newsampleprob),prob=newsampleprob,replace=TRUE))
          # ) < 20) {
          #   if(counter==1) message ('Sampling problematic -- trying to recover... ')
          #   counter = counter + 1
          #   if(counter == 30) stop('Importance sampling failed -- either posterior mode not found, or mode is inadequate starting point for sampling')
          #   weighted_dens <- weighted_dens /2
          #   newsampleprob <- exp((weighted_dens - log_sum_exp(weighted_dens)))
          #   # plot(newsampleprob)
          # }
          sample_prob <- c(sample_prob,newsampleprob) #sum to 1 for each iteration, normalise later
          
          sample_prob[!is.finite(sample_prob)] <- 0
          sample_prob[is.na(sample_prob)] <- 0
          
          
          # #check max resample probability and drop earlier samples if too high
          # dropuntil <- ceiling(max(c(0,which(sample_prob > (chancethreshold / isloopsize)) / isloopsize),na.rm=TRUE))*isloopsize
          # if((isloopsize - dropuntil) > isloopsize) dropuntil <- dropuntil -isloopsize
          # if(nrow(samples)-dropuntil < isloopsize*2) dropuntil <- nrow(samples)-isloopsize*2
          # if(nrow(samples) <= isloopsize *2) dropuntil <- 0
          # 
          # if(dropuntil > 0){
          #   targetvec <- targetvec[-(0:dropuntil)]
          #   sample_prob <- sample_prob[-(0:dropuntil)]
          #   samples <- samples[-(0:dropuntil),,drop=FALSE]
          # }
          
          
          # if(plot&runif(1,0,1) > .98){
          #   par(mfrow=c(1,1))
          #   lprob <- Var2 <-value <- NULL
          #   msamps = cbind(data.table(lprob=rep(log(sample_prob),ncol(samples))),
          #     data.table::melt(data.table(cbind(samples))))
          #   ggplot2::ggplot(dat=msamps,aes(y=lprob,x=value,alpha=exp(msamps$lprob))) + 
          #     geom_point(colour='red') +
          #     scale_alpha_continuous(range=c(0,1))+
          #     theme_minimal()+
          #     facet_wrap(vars(Var2),scales='free') 
          # }
          
          
          if(nrow(samples) >= targetsamples && (max(sample_prob)/ sum(sample_prob)) < chancethreshold) {
            message ('Finishing importance sampling...')
            nresamples <- finishsamples 
          }else nresamples = max(5000,nrow(samples)/5)
          
          
          # points(target_dens2[sample_prob> (1/isloopsize * 10)], prop_dens[sample_prob> (1/isloopsize * 10)],col='red')
          resample_i <- sample(1:nrow(samples), size = nresamples, replace = ifelse(nrow(samples) > targetsamples,FALSE,TRUE),
            prob = sample_prob / sum(sample_prob))
          # resample_i <- sample(tail(1:nrow(samples),isloopsize), size = nresamples, replace = ifelse(j == isloops+1,FALSE,TRUE),
          #   prob = tail(sample_prob,isloopsize) / sum(tail(sample_prob,isloopsize) ))
          
          message(paste0('Importance sample loop ',j,', ',length(unique(resample_i)), ' unique samples, from ', nresamples,' resamples of ', nrow(samples),' actual, prob sd = ', round(sd(sample_prob),4),
            ', max chance = ',max(sample_prob) * isloopsize))
          if(length(unique(resample_i)) < 100) {
            message('Sampling ineffective, unique samples < 100 -- try increasing samples per step (isloopsize), or use HMC (non optimizing) approach.')
            # return(est)
          }
          
          resamples <- samples[resample_i, , drop = FALSE]
          
          ess[j] <- (sum(sample_prob[resample_i]))^2 / sum(sample_prob[resample_i]^2)
          qdiag[j]<-mean(unlist(lapply(sample(x = 1:length(sample_prob),size = 500,replace = TRUE),function(i){
            (max(sample_prob[resample_i][1:i])) / (sum(sample_prob[resample_i][1:i]) )
          })))
          
        }
      }
    }
  }#end while no better fit
  if(!estonly){
    if(!is) lpsamples <- NA else lpsamples <- unlist(target_dens)[resample_i]
    
    message('Computing posterior with ',nrow(resamples),' samples')
    standata$savesubjectmatrices=savesubjectmatrices
    
    if(!savesubjectmatrices) sdat=standatact_specificsubjects(standata,1) #only use 1 subject
    if(savesubjectmatrices) sdat=standata
    
    transformedpars=stan_constrainsamples(sm = sm,standata = sdat,
      savesubjectmatrices = savesubjectmatrices, savescores = standata$savescores,
      dokalman=as.logical(standata$savesubjectmatrices),
      samples=resamples,cores=cores, cl=benv$clctsem, quiet=TRUE)
    
    if(cores > 1) {
      parallel::stopCluster(benv$clctsem)
      smf <- stan_reinitsf(sm,standata)
    }
    
    
    
    # transformedparsfull=stan_constrainsamples(sm = sm,standata = standata,
    #   savesubjectmatrices = TRUE, dokalman=TRUE,savescores = TRUE,
    #   samples=matrix(est2,nrow=1),cores=1, quiet = TRUE)
    
    
    
    sds=try(suppressWarnings(sqrt(diag(mcov))))  #try(sqrt(diag(solve(optimfit$hessian))))
    if('try-error' %in% class(sds)[1]) sds <- rep(NA,length(est2))
    lest= est2 - 1.96 * sds
    uest= est2 + 1.96 * sds
    
    transformedpars_old=NA
    try(transformedpars_old<-cbind(unlist(constrain_pars(smf, upars=lest)),
      unlist(constrain_pars(smf, upars= est2)),
      unlist(constrain_pars(smf, upars= uest))),silent=TRUE)
    try(colnames(transformedpars_old)<-c('2.5%','mean','97.5%'),silent=TRUE)
    stanfit=list(optimfit=optimfit,stanfit=stan_reinitsf(sm,standata), rawest=est2, rawposterior = resamples, cov=mcov,
      transformedpars=transformedpars,transformedpars_old=transformedpars_old,
      # transformedparsfull=transformedparsfull,
      standata=list(TIPREDEFFECTsetup=standata$TIPREDEFFECTsetup,ntipredeffects = standata$ntipredeffects),
      isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
    if(bootstrapUncertainty) stanfit$subjectscores <- scores #subjectwise gradient contributions
  }
  
  if(estonly) {
    smf <- stan_reinitsf(sm,standata)
    stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2,parsteps=parsteps)
  }
  optimfinished <- TRUE #disable exit message re pars
  return(stanfit)
}



