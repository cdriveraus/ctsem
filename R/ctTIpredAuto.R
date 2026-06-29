ctRawParnamesRaw <- function(fit){
  if(!class(fit) %in% 'ctStanFit' || !length(fit$stanfit$stanfit@sim)==0) stop('Not an optimized ctStanFit model!')
  ms <- fit$setup$matsetup
  ms <- ms[ms$param > 0 & ms$when %in% 0 & ms$copyrow < 1,]
  ms <- ms[order(ms$param),]
  rawpopmeans<-ms$parname
  iv <- ms$parname[ms$indvarying >0]
  rawpopsd <-c()
  rawpopucorr <- c()
  if(length(iv) > 0){
    rawpopsd <- paste0('rawpopsd_',iv)
    rawpopucorr <- matrix(paste0('rawpopucorr_',rep(iv,times=length(iv)),'_',rep(iv,each=length(iv))),length(iv),length(iv))
    rawpopucorr <- rawpopucorr[lower.tri(rawpopucorr)]
  } 
  tipredeffects <- c()
  if(fit$standata$ntipred > 0){
    tipredeffects <- rep(NA,max(fit$standata$TIPREDEFFECTsetup))
    counter <- 0
    for(j in 1:ncol(fit$standata$TIPREDEFFECTsetup)){
      for(i in 1:nrow(fit$standata$TIPREDEFFECTsetup)){
        if(fit$standata$TIPREDEFFECTsetup[i,j] > 0){
          counter <- counter + 1
          tipredeffects[counter] <- paste0(rawpopmeans[i],'_', fit$ctstanmodel$TIpredNames[j])
        }
      }
    }
  }
  return(list(rawpopmeans=rawpopmeans,rawpopsd=rawpopsd,rawpopucorr=rawpopucorr,tipredeffects=tipredeffects))
}

TIPREDEFFECTnames <- function(fit){
  if(fit$standata$ntipred > 0){
    ms <- fit$setup$matsetup
    tie <- fit$standata$TIPREDEFFECTsetup
    colnames(tie) <- paste0(fit$ctstanmodel$TIpredNames)
    pars <- ms[ms$param > 0 & ms$when %in% 0 & ms$copyrow < 1,]
    pars <- pars[order(pars$param),]
    rownames(tie) <- pars$parname
    return(tie)
  } else return(c())
}

popcovnames <- function(fit){
  if(fit$standata$nindvarying > 0){
    ms <- fit$setup$matsetup
    pars <- ms[ms$param > 0 & ms$when %in% 0 & ms$copyrow < 1 & ms$indvarying > 0,]
    pars <- pars[order(pars$param),]
    popcov <- matrix(NA,fit$standata$nindvarying ,fit$standata$nindvarying,
      dimnames = list( pars$parname, pars$parname) )
    return(popcov)
  } else return(c())
}




checkTIauto <- function(){
  Tpoints=30
  n.manifest=1
  n.TDpred=0
  n.TIpred=1
  n.latent=1
  n.subjects=100
  TI1 <- rnorm(n.subjects)
  gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
    n.TDpred=n.TDpred,n.manifest=n.manifest,
    MANIFESTVAR=diag(0.5,1),
    LAMBDA=diag(1,1),
    DRIFT=matrix(c(-.3),nrow=1),
    DIFFUSION=matrix(c(2),1),
    T0VAR=diag(10,1))
  
  for(i in 1:n.subjects){
    gm$CINT[1,1] <- TI1[i]*.5+rnorm(1,0,.6)
    ndat<-ctGenerate(gm,n.subjects=1,burnin=30,logdtsd=.4)
    ndat <- cbind(ndat,TI1[i])
    ndat[,1] <- i
    if(i>1) tdat <- rbind(tdat,ndat) else tdat <- ndat
  }
  colnames(tdat)[4] <- 'TI1'
  tdat <- data.frame(tdat)
  tdat$TI2 <- rnorm(nrow(tdat))
  # colnames(tdat)[5] <- 'TI2'
  
  tdat[2,'Y1'] <- NA
  tdat[tdat[,'id']==2,'TI1'] <- NA
  
  checkm<-ctModel(type='ct',Tpoints=Tpoints,
    MANIFESTVAR=diag(0.5,1),
    # DRIFT=matrix(c(-.3),nrow=1),
    # DIFFUSION=matrix(c(2),1),
    n.latent=n.latent,n.TDpred=n.TDpred,
    n.TIpred=2,
    MANIFESTMEANS=matrix(0,nrow=n.manifest),
    CINT=matrix(c('cint1'),ncol=1),
    n.manifest=n.manifest,LAMBDA=diag(1))
  
  # checkm$pars$indvarying[!checkm$pars$matrix %in% 'T0MEANS'] <- FALSE
  
  checkm$TIpredAuto <- 1L
  
  fit1<-ctFit(tdat,checkm,chains=1,optimize=TRUE,cores=1,verbose=0,
    # intoverpop=F,
    plot=10,
    # savesubjectmatrices = F,plot=F,
    # init=init,
    # fit=F,
    priors=T)
  summary(fit1)
}

whichsubjectpars <- function(standata,subjects=NA){
  a1=standata$nparams+standata$nindvarying+
    (standata$nindvarying^2-standata$nindvarying)/2
  whichbase <- 1:a1
  if(standata$intoverpop ==0 && standata$nindvarying > 0){ #then there are subject pars
    whichsubjects <- a1+cseq(from=subjects,to=standata$nindvarying*standata$nsubjects,
      by=standata$nsubjects)
    whichbase <- c(whichbase,whichsubjects)
  }
  if(standata$ntipredeffects > 0) {
    tipredstart <- (a1+
        ifelse(standata$intoverpop,0,standata$nindvarying*standata$nsubjects)+1)
    whichbase <- c(whichbase,tipredstart:(tipredstart+standata$ntipredeffects -1
        # ifelse(standata$doonesubject >0,0,-1)
        #disabled the doonesubject thing
      ))
  }
  return(whichbase)
}


scorecalc <- function(standata, est, stanmodel, subjectsonly = TRUE, 
  returnsubjectlist = TRUE, cores = 2) {
  # Set prior modification factor based on whether we use subject-only data
  standata$priormod <- ifelse(subjectsonly, 1 / standata$nsubjects, 1 / standata$ndatapoints)
  
  # Try to initialize fast Stan function
  sf <- suppressMessages(try(stan_reinitsf(stanmodel, standata, fast = TRUE)))
  fast <- !inherits(sf, "try-error")
  
  # Function to compute gradients for a single subject
  compute_subject_gradients <- function(subject_index) {
    whichpars <- whichsubjectpars(standata, subject_index)  # Parameters for the subject
    scores_subject <- matrix(NA, nrow = length(whichpars), 
      ncol = ifelse(subjectsonly, 1, sum(standata$subject == subject_index)))
    
    # Create subject-specific data
    standata1 <- standatact_specificsubjects(standata, subject_index)
    
    # Compute gradients for each data row (or overall for the subject)
    for (j in seq_len(ncol(scores_subject))) {
      standata1$llsinglerow <- as.integer(ifelse(subjectsonly, 0, j))
      sf <- stan_reinitsf(stanmodel, standata1, fast = fast)
      
      if (fast) {
        scores_subject[, j] <- sf$grad_log_prob(
          upars = est[whichpars],
          adjust_transform = TRUE
        )
      } else {
        scores_subject[, j] <- rstan::grad_log_prob(
          sf,
          upars = est[whichpars],
          adjust_transform = TRUE
        )
      }
    }
    return(scores_subject)
  }
  
  # Parallel processing: compute gradients for all subjects
  cl <- parallelly::makeClusterPSOCK(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("standata", "est", "stanmodel", "subjectsonly", "compute_subject_gradients", 
    "stan_reinitsf", "whichsubjectpars", "standatact_specificsubjects"),envir=environment())

  group_size <- ceiling(standata$nsubjects / length(cl))
  
  # Create groups
  task_chunks <- split(1:standata$nsubjects, ceiling(seq_along(1:standata$nsubjects) / group_size))
  # Parallel execution
  scores <- do.call(c, parallel::parLapply(cl, task_chunks, function(chunk) {
    lapply(chunk, compute_subject_gradients)
  }))
  
  
  # Combine results if subjectsonly is TRUE
  if (subjectsonly) {
    scores <- matrix(unlist(scores), nrow = length(scores[[1]]))
  }
  
  # Optionally return as data.table instead of list
  if (!returnsubjectlist) {
    if (is.list(scores)) {
      scores <- lapply(scores, function(x) data.table(t(x)))
      scores <- rbindlist(scores)
    } else {
      scores <- t(scores)
    }
  }
  
  return(scores)
}


ctTIauto <- function(fit,tipreds=NA,cores=2){
  if(is.na(tipreds[1])) tipreds <- fit$standata$tipredsdata
  # colnames(tipreds) <- paste0('ti',1:ncol(tipreds))
  
  if(is.null(fit$stanfit$subjectscores)) scores <- scorecalc(standata = fit$standata,
    est = fit$stanfit$rawest,stanmodel = fit$stanmodel,cores=cores)
  else scores <- fit$stanfit$subjectscores
  scores <- scores[1:fit$standata$nparams,,drop=FALSE]
  rownames(scores) <- paste0('p',1:nrow(scores))
  # matchindex <- match(1:fit$standata$nparams,fit$setup$matsetup$param)
  # rownames(scores)[1:fit$standata$nparams] <- fit$standata$matsetup$parname[match(1:fit$standata$nparams,fit$setup$matsetup$param)]
  sc <- list()
  for(i in 1:nrow(scores)){
    # for(j in 1:ncol(tipreds)){
    # plot(sort(tipreds[,j]),scores[i,][order(tipreds[,j])],ylab=rownames(scores)[i],xlab=colnames(tipreds)[j])
    # }
    sc[[i]]=summary(lm(scores[i,] ~ tipreds))$coefficients
  }
  names(sc)[1:fit$standata$nparams]<-paste0('p',1:nrow(scores)) #fit$setup$matsetup$parname[match(1:fit$standata$nparams,fit$setup$matsetup$param)]
  
  s2=lapply(sc,function(x) {
    x=x[-1,,drop=FALSE]
    rownames(x) <- gsub('^tipreds','',rownames(x))
    rownames(x) <- paste0('ti',1:nrow(x))
    return(x)
  })
  
  TIPREDEFFECTsetup = matrix(NA,length(s2),nrow(s2[[1]]))
  for(i in 1:length(s2)){
    TIPREDEFFECTsetup[i,] <- s2[[i]][,4]
  }
  
  if(fit$standata$nindvarying > 0 && fit$standata$intoverpop > 0){
    fit$setup$matsetup <- data.frame(fit$standata$matsetup)
    e=stan_constrainsamples(sm = fit$stanmodel,standata = fit$standata,
      samples = matrix(fit$stanfit$rawest,nrow=1),savescores=TRUE,quiet=TRUE,pcovn=2)
    p=sort(unique(fit$setup$matsetup$row[fit$setup$matsetup$indvarying>0]))# | fit$setup$matsetup$tipred]))
    firstsub <- rep(TRUE,fit$standata$ndatapoints) #which rows represent first rows per subject
    for(i in 2:fit$standata$ndatapoints){
      if(fit$standata$subject[i] == fit$standata$subject[i-1]) firstsub[i] <- FALSE
    }
    e$etasmooth <-  array(e$etaa[,3,,,drop=FALSE],dim=dim(e$etaa)[-2])
    states <- ctCollapse(e$etasmooth[,firstsub,p,drop=FALSE],1,mean)
    sc=list()
    for(i in 1:ncol(states)){
      sc[[i]]=summary(lm(states[,i] ~ tipreds))$coefficients[-1,,drop=FALSE]
    }
    for(i in 1:length(sc)){
      TIPREDEFFECTsetup[fit$setup$matsetup$param[fit$setup$matsetup$indvarying %in% i],] <- sc[[i]][,4]
    }
  }
  
  if(any(is.na(TIPREDEFFECTsetup))) warning('NA found, probably unused parameters?')
  TIPREDEFFECTsetup[is.na(TIPREDEFFECTsetup)] <- 1
  return(TIPREDEFFECTsetup)
}
