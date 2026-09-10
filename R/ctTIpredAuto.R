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
  # rscript_libs: without it a worker searches its own default .libPaths(),
  # not the caller's, so the exported ctsem internals below (stan_reinitsf
  # etc.) can resolve against a different install than the one running this
  # code -- the same failure mode makeClusterID() (stanoptimis.R) was fixed
  # for. It narrows that gap rather than closing it, so
  # .ctClusterCheckBuild() below reports a build the workers do not share.
  # library(ctsem) is loaded explicitly too, since the closures exported just
  # below reference the ctsem namespace directly rather than requesting it.
  cl <- parallelly::makeClusterPSOCK(cores, rscript_libs = .libPaths())
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(ctsem)))
  .ctClusterCheckBuild(cl)
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
