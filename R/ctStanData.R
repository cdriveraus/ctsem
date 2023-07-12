ctStanData <- function(ctm, datalong,optimize,derrind='all',sameInitialTimes=FALSE){
  
  
  nsubjects <- length(unique(datalong[, ctm$subjectIDname])) 
  
  mats <- ctStanMatricesList()
  
  #simply exponential?
  driftdiagonly <- ifelse(all(!is.na(ctm$pars$value[ctm$pars$matrix == 'DRIFT' & ctm$pars$row != ctm$pars$col]) &
      all(ctm$pars$value[ctm$pars$matrix == 'DRIFT' & ctm$pars$row != ctm$pars$col] == 0) ), 1, 0)
  
  
  ###data checks
  
  if(any(!c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames) %in% colnames(datalong))) stop(paste0('
      variables: ', paste0(c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames)[
        which(!c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames) %in% colnames(datalong))], ', '),' not in data'))
  
  if(!(ctm$subjectIDname %in% colnames(datalong))) stop(paste('id column', (ctm$subjectIDname), "not found in data"))
  
  
  if(ctm$priors){
    if(ctm$n.TIpred > 1 && any(abs(colMeans(datalong[,c(ctm$TIpredNames),drop=FALSE],na.rm=TRUE)) > .3)){
      message('Uncentered TI predictors noted -- interpretability may be hindered and default priors may not be appropriate')
    }
    
    # #scale check
    # if(naf(any(abs(colMeans(datalong[,c(ctm$manifestNames,ctm$TDpredNames),drop=FALSE],na.rm=TRUE)) > 5))){
    #   message('Uncentered data noted -- default priors *may* not be appropriate')
    # }
    # 
    # 
    # if(naf(any(abs(apply(datalong[,c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames),drop=FALSE],2,sd,na.rm=TRUE)) > 3))){
    #   message('Unscaled data noted -- default priors may not be appropriate')
    # }
  }
  
  
  #id mapping
  original <- unique(datalong[,ctm$subjectIDname])
  datalong <- makeNumericIDs(datalong,ctm$subjectIDname,ctm$timeName)
  new <- unique(datalong[,ctm$subjectIDname])
  idmap <- data.frame(original, new)
  
  
  
  
  if( (nrow(ctm$t0varstationary) + nrow(ctm$t0meansstationary)) >0 && 
      length(c(ctm$modelmats$calcs$driftcint, ctm$modelmats$calcs$diffusion)) > 0) message('Stationarity assumptions based on initial states when using non-linear dynamics')
  
  
  if(ctm$intoverstates==FALSE || all(derrind=='all') ) derrind = 1:ctm$n.latent
  # if(all(derrind=='all')) derrind = sort(unique(ctm$pars$col[
  #   ctm$pars$matrix=='DIFFUSION' & (!is.na(ctm$pars$param) | ctm$pars$value!=0)]))
  derrind = as.integer(derrind)
  if(any(derrind > ctm$n.latent)) stop('derrind > than n.latent found!')
  if(length(derrind) > ctm$n.latent) stop('derrind vector cannot be longer than n.latent!')
  if(length(unique(derrind)) < length(derrind)) stop('derrind vector cannot contain duplicates or!')
  ndiffusion=length(derrind)
  
  
  
  nindvarying <- max(ctm$modelmats$matsetup$indvarying)
  nparams <- max(ctm$modelmats$matsetup$param[ctm$modelmats$matsetup$when %in% c(0,-1)])
  nmatrices <- length(mats$base)
  ctm$modelmats$matsetup[which(ctm$modelmats$matsetup$indvarying > 0),]
  indvaryingindex <- ctm$modelmats$matsetup$param[which(ctm$modelmats$matsetup$indvarying > 0)]
  indvaryingindex <- array(indvaryingindex[!duplicated(indvaryingindex)])
  
  sdscale <- array(ctm$modelmats$matvalues$sdscale[match(indvaryingindex,ctm$modelmats$matsetup$param)])
  
  
  
  
  #tipred data
  if(ctm$n.TIpred > 0) {
    tipreds <- datalong[match(unique(datalong[,ctm$subjectIDname]),datalong[,ctm$subjectIDname]),ctm$TIpredNames,drop=FALSE]
    if(any(is.na(tipreds))) {
      if(!optimize){
        message(paste0("NA's in TIpreds - sampling ", sum(is.na(tipreds)),' values'))
        tipreds[is.na(tipreds)] = 99999
      }
      if(optimize){
        message(paste0("NA's in TIpreds - imputing ", sum(is.na(tipreds)),'  NA\'s to allow optimization -- TIpred effect estimates may be overly confident.'))
        # tipreds[is.na(tipreds)] = 0
        timu <- apply(tipreds,2,mean,na.rm=TRUE)
        tisd <- apply(tipreds,2,sd,na.rm=TRUE)
        
        meandat <- data.table((datalong))[ , lapply(.SD, function(x) 
          sum(x,na.rm=TRUE)/sum(!is.na(x))) ,
          # x[1]),
          by=c(ctm$subjectIDname)]
        # sddat <- data.table((datalong))[ , lapply(.SD, function(x) 
        #   sd(x,na.rm=TRUE)) , 
        #   by=c(ctm$subjectIDname)]
        # sddat<-sddat[,!colnames(sddat) %in% ctm$subjectIDname,with=FALSE]
        meandat <- data.frame(scale(meandat[,apply(meandat,2,sd,na.rm=TRUE) > 1e-4,with=FALSE]))
        meandat[is.na(meandat)] <- 0
        # sddat <- sddat[,apply(sddat,2,sd,na.rm=TRUE) > 1e-4,with=FALSE]
        # colnames(sddat) <- paste0('sd_',colnames(sddat))
        # meandat <- cbind(meandat,sddat)
        
        
        # cml <- covml(meandat,reg = TRUE)
        
        for(i in 1:ctm$n.TIpred){
          # cnames <- c(colnames(meandat)[!colnames(meandat) %in% ctm$TIpredNames[i]],ctm$TIpredNames[i])
          # 
          # 
          # chol <- cml$cp$covm
          # dimnames(chol) <- list(colnames(meandat),colnames(meandat))
          # chol <- t(chol(chol[cnames,cnames]))
          # 
          #  tipreds[is.na(tipreds[,i]),i] <- 
          #    apply(meandat[is.na(tipreds[,i]),head(cnames,length(cnames)-1)],1,function(x){
          #      sum(x * c(tail(chol,1))[-ncol(chol)])
          #    }) *tisd[i] + timu[i]
          
          lmform = formula(paste0(ctm$TIpredNames[i],' ~ 1 + ',
            paste0(colnames(meandat)[-which(colnames(meandat) %in% ctm$TIpredNames[i])],
              collapse=' + ')))
          # lmr <- MASS::lm.ridge(formula = lmform, data = meandat,lambda=0,na.action=na.exclude)
          lmf <- lm(formula = lmform,data = meandat,na.action=na.exclude)
          # lmf$coefficients[-1] <- lmr$coef
          # lmf$coefficients[1] <- lmr$ym
          # lmr$coef * 
          #   apply(meandat[,-which(colnames(meandat) %in% ctm$TIpredNames[i]),with=FALSE],
          #     2,sd,na.rm=TRUE) - 
          #   lmf$coefficients[-1]
          # 
          # plot(meandat[,ctm$TIpredNames[i]],predict(lmf),main=ctm$TIpredNames[i],col=as.numeric(is.na(tipreds[,i]))+1)
          
          tipreds[is.na(tipreds[,i]),i] <- predict(lmf)[is.na(tipreds[,i])] * tisd[i] + timu[i]
        }
      }
    }
  }
  
  if(sameInitialTimes){
    datanew <- datalong[!duplicated(datalong[[ctm$subjectIDname]]),]
  datanew[[ctm$timeName]] <- min(datanew[[ctm$timeName]])
  datalong <- merge.data.frame(x = datalong,y = datanew[,c(ctm$subjectIDname,ctm$timeName)],
    by = c(ctm$subjectIDname,ctm$timeName),all = TRUE)
  }
  
  datalong[,ctm$manifestNames][is.na(datalong[,ctm$manifestNames])]<-99999 #missing data
  
  
  standata <- list(
    Y=cbind(as.matrix(datalong[,ctm$manifestNames])),
    subject=array(as.integer(datalong[,ctm$subjectIDname])),
    time=array(datalong[,ctm$timeName]), 
    ndatapoints=as.integer(nrow(datalong)),
    nobs_y=array(as.integer(apply(datalong[,ctm$manifestNames,drop=FALSE],1,function(x) length(x[x!=99999]))),dim=nrow(datalong)),
    whichobs_y=matrix(as.integer(t(apply(datalong[,ctm$manifestNames,drop=FALSE],1,function(x) {
      out<-as.numeric(which(x!=99999))
      if(length(out)==0) out<-rep(0,ctm$n.manifest)
      if(length(out)<ctm$n.manifest) out<-c(out,rep(0,ctm$n.manifest-length(out)))
      out
    }) )),nrow=c(nrow(datalong),ncol=ctm$n.manifest)),
    nbinary_y=array(as.integer(apply(datalong[,ctm$manifestNames,drop=FALSE],1,function(x) 
      length(x[ctm$manifesttype==1 & x!=99999]))),dim=nrow(datalong)),
    whichbinary_y=matrix(as.integer(t(apply(datalong[,ctm$manifestNames,drop=FALSE],1,function(x) {
      out<-as.numeric(which(ctm$manifesttype==1 & x!=99999)) #conditional on whichobs_y
      if(length(out)==0) out<-rep(0,ctm$n.manifest)
      if(length(out)<ctm$n.manifest) out<-c(out,rep(0,ctm$n.manifest-length(out)))
      out
    }) )),nrow=c(nrow(datalong),ncol=ctm$n.manifest)),
    ncont_y=array(as.integer(apply(datalong[,ctm$manifestNames,drop=FALSE],1,function(x) 
      length(x[(ctm$manifesttype==0 | ctm$manifesttype==2) & x!=99999]))),dim=nrow(datalong)),
    whichcont_y=matrix(as.integer(t(apply(datalong[,ctm$manifestNames,drop=FALSE],1,function(x) {
      out<-as.numeric(which( (ctm$manifesttype==0 | ctm$manifesttype==2) & x!=99999)) #conditional on whichobs_y
      if(length(out)==0) out<-rep(0,ctm$n.manifest)
      if(length(out)<ctm$n.manifest) out<-c(out,rep(0,ctm$n.manifest-length(out)))
      out
    }) )),nrow=c(nrow(datalong),ncol=ctm$n.manifest))
  )
  
  if(ctm$n.TDpred > 0) {
    tdpreds <- datalong[,ctm$TDpredNames,drop=FALSE]
    if(any(is.na(tdpreds))) {
      # if(NAtdpreds == 'error'
      message("NA's in TDpreds! Replaced by zeroes, consider appropriateness...")
      tdpreds[is.na(tdpreds)] <-0 ## rough fix for missingness
    }
  }
  if(ctm$n.TDpred ==0) tdpreds <- matrix(0,standata$ndatapoints,0) #standata$ndatapoints,
  standata$tdpreds=array(as.matrix(tdpreds),dim=c(nrow(tdpreds),ncol(tdpreds)))
  
  #subset selection
  if(is.null(ctm$dokalmanrows)) standata$dokalmanrows <- 
    array(rep(1L, standata$ndatapoints)) else standata$dokalmanrows <- array(as.integer(ctm$dokalmanrows))
  standata$priormod = 1L#sum(standata$dokalmanrows)/standata$ndatapoints
  
  standata<-c(standata, 
    list(
      nsubjects=as.integer(nsubjects),
      nmanifest=as.integer(ctm$n.manifest),
      nlatentpop = ctm$nlatentpop,
      Jstep = ctm$nlcontrol$Jstep,
      maxtimestep = ctm$nlcontrol$maxtimestep,
      dokalman=as.integer(is.null(ctm$dokalman)),
      intoverstates=as.integer(ctm$intoverstates),
      verbose=0L,
      manifesttype=array(as.integer(ctm$manifesttype),dim=length(ctm$manifesttype)),
      indvaryingindex=array(as.integer(indvaryingindex)),
      intoverpopindvaryingindex=array(as.integer(ctm$intoverpopindvaryingindex)),
      notindvaryingindex=array(as.integer(which(!(1:nparams) %in% indvaryingindex))),
      continuoustime=as.integer(sum(ctm$continuoustime)),
      nlatent=as.integer(ctm$n.latent),
      ntipred=as.integer(ctm$n.TIpred),
      ntdpred=as.integer(ctm$n.TDpred),
      nparams=as.integer(nparams),
      gendata=0L,
      nindvarying=as.integer(nindvarying),
      nindvaryingoffdiagonals=as.integer((nindvarying^2-nindvarying)/2),
      nt0varstationary=as.integer(nrow(ctm$t0varstationary)),
      nt0meansstationary=as.integer(nrow(ctm$t0meansstationary)),
      t0varstationary=matrix(as.integer(ctm$t0varstationary),ncol=2),
      t0meansstationary=matrix(as.integer(ctm$t0meansstationary),ncol=2),
      derrind=array(as.integer(derrind),dim=ndiffusion),
      ndiffusion=as.integer(ndiffusion),
      driftdiagonly = as.integer(driftdiagonly),
      intoverpop=as.integer(ctm$intoverpop),
      # nlmeasurement=as.integer(nlmeasurement),
      priors=as.integer(ctm$priors),
      savescores=0L,
      savesubjectmatrices=0L,
      nJAxfinite = ifelse(ctm$recompile,0L,length(ctm$JAxfinite)),
      nJyfinite = ifelse(ctm$recompile,0L,length(ctm$Jyfinite))
    ))
  
  if(ctm$n.TIpred == 0) tipreds <- array(0,c(0,0))
  standata$tipredsdata <- as.matrix(tipreds)
  standata$nmissingtipreds <- as.integer(length(tipreds[tipreds== 99999]))
  
  standata$ntipredeffects <- as.integer(ifelse(ctm$n.TIpred > 0, as.integer(max(ctm$modelmats$TIPREDEFFECTsetup)), 0))
  standata$TIPREDEFFECTsetup <- apply(ctm$modelmats$TIPREDEFFECTsetup,c(1,2),as.integer,.drop=FALSE)
  standata$tipredsimputedscale <- ctm$tipredsimputedscale
  standata$tipredeffectscale <- ctm$tipredeffectscale
  
  
  # #drift, jacobian off diagonal check # seems to work ok internally with expm2
  # mx=listOfMatrices(ctm$pars)
  # driftcint <- rbind(cbind(mx$DRIFT[1:ctm$n.latent,1:ctm$n.latent,drop=FALSE],mx$CINT),0)
  # z=c()
  # for(ri in 1:nrow(driftcint)){
  #   if(all(driftcint[ri,-ri] %in% 0) && all(driftcint[-ri,ri] %in% 0)) z[ri]=0L else z[ri]=1L
  # }
  # standata$drcintoffdiag <- array(as.integer(c(z,1)),dim=nrow(driftcint))
  # 
  # jac <- mx$JAx
  # z=c()
  # for(ri in 1:nrow(jac)){
  #   if(all(jac[ri,-ri] %in% 0) && all(jac[-ri,ri] %in% 0)) z[ri]=0L else z[ri]=1L
  # }
  # standata$jacoffdiag <- array(as.integer(c(z,1)),dim=nrow(jac))
  # standata$jacoffdiagindex <- array(as.integer(sort(unique(c(1:ctm$n.latent,which(standata$jacoffdiag ==1))))))
  # standata$njacoffdiagindex <- as.integer(length(standata$jacoffdiagindex))
  
  
  standata$JAxfinite <- ctm$JAxfinite
  standata$Jyfinite <- ctm$Jyfinite
  
  standata$Jycolindexsize <- 1L
  standata$Jycolindex <- array(1L)
  
  standata$difftype <- 0L;
  standata$dotipred <- 1L;
  
  
  if(!ctm$timeName %in% colnames(datalong) && !ctm$continuoustime) {
    datalong[[ctm$timeName]] <- 1:nrow(datalong)
  }
  
  # #t0 index
  # T0check<-rep(1,nrow(datalong))
  # for(i in 2:nrow(datalong)){
  #   T0check[i]<- ifelse(datalong[i,ctm$subjectIDname] != datalong[i-1,ctm$subjectIDname], 1, 0)
  # }
  if (!(ctm$timeName %in% colnames(datalong))) stop(paste('time column', (ctm$timeName), "not found in data"))
  if(any(is.na(datalong[,ctm$timeName]))) stop('Missings in time column!')
  
  oldsubi<-datalong[1,ctm$subjectIDname]-1
  dT<-rep(-1,length(datalong[,ctm$timeName]))
  
  for(rowi in 1:length(datalong[,ctm$timeName])) {
    subi<-datalong[rowi,ctm$subjectIDname]
    # if(rowi==1 && subi!=1) stop('subject id column must ascend from 1 to total subjects without gaps')
    # if(oldsubi!=subi && subi-oldsubi!=1) stop('subject id column must ascend from 1 to total subjects without gaps')
    if(subi - oldsubi == 1) {
      dT[rowi]<-0
      subistartrow<-rowi
    }
    if(subi - oldsubi == 0) {
      dT[rowi]<-datalong[rowi,ctm$timeName] - datalong[rowi-1,ctm$timeName]
      if(dT[rowi] < 0) stop(paste0('A time interval of ', dT[rowi],' was found at row ',rowi))
      if(dT[rowi] == 0) warning(paste0('A time interval of ', dT[rowi],' was found at row ',rowi))
      
    }
    oldsubi<-subi
  }
  
  # if(mean(dT) > 3) message('Average time interval > 3 -- in typical settings this can be too large, consider time scaling...')
  
  
  
  subindices <- lapply(mats$base,function(x) 0)
  for(mati in mats$base){
    if( (!ctm$intoverpop && any(ctm$pars$indvarying[ctm$pars$matrix %in% names(mats$base)[mati]])) || 
        (ctm$n.TIpred >0 && (
          any(unlist(ctm$pars[ctm$pars$matrix %in%  names(mats$base)[mati],paste0(ctm$TIpredNames,'_effect')])) || 
            any(ctm$pars$matrix %in%  names(mats$base)[mati] & grepl('[',ctm$pars$param,fixed=TRUE))  )
        )) subindex <- 1 else subindex <- 0
        subindices[[names(mats$base)[mati]]] <- subindex
  }
  
  if(ctm$stationary || nrow(ctm$t0varstationary) > 0) subindices$T0VAR  <- 
    max(c(subindices$T0VAR,subindices$DRIFT,subindices$DIFFUSION))
  
  if(ctm$stationary || nrow(ctm$t0meansstationary) > 0) subindices$T0MEANS <- 
    max(c(subindices$T0MEANS,subindices$DRIFT,subindices$CINT))
  
  
  
  standata$subindices <- as.integer(unlist(subindices))[order(mats$base)]
  
  
  #state dependence
  statedep <- rep(0L,as.integer(max(mats$all)))
  lhscalcs <- sapply(unique(unlist(ctm$modelmats$calcs)),function(x) gsub('=.*','',x))
  for(i in 1:length(statedep)){
    matname <- try(names(mats$all[mats$all %in% i]),silent=TRUE)
    if(length(matname)==0) next
    statedep[i] <- ifelse(any(
      ctm$modelmats$matsetup$when >0 &
        ctm$modelmats$matsetup$matrix %in% i),1L,0L)
    if(any(sapply(lhscalcs,function(calci) 
      grepl(matname,calci)))) statedep[i] <- 1L #if any calcs modify this matrix, make state dependent
  }
  # 
  # statedep=rep(0L,4)
  # v = 0L
  # if(any(ctm$modelmats$matsetup$when == 2) ||
  #     any(grepl('state[',ctm$jacobian$JAx,fixed=TRUE)) ||
  #     any(grepl('state[',ctm$modelmats$calcs$driftcint,fixed=TRUE)) ||
  #     any(grepl('state[',ctm$modelmats$calcs$diffusion,fixed=TRUE)) 
  # ) statedep[2] = 1L
  # 
  # if(any(ctm$modelmats$matsetup$when == 2 & ctm$modelmats$matsetup$matrix == 4) ||
  #     any(grepl('state[',ctm$modelmats$calcs$diffusion,fixed=TRUE)) 
  # ) multiplicativenoise = 1L
  
  standata$statedep <- statedep
  # standata$nstatedep <- as.integer(length(statedep))
  # standata$multiplicativenoise <- multiplicativenoise
  standata$choleskymats<- 0L
  if(!is.null(ctm$covmattransform)){
    if(ctm$covmattransform=='rawcorr_indep') standata$choleskymats<- -1L
    if(ctm$covmattransform=='cholesky') standata$choleskymats<- 1L
    if(!ctm$covmattransform %in% c('rawcorr','rawcorr_indep','cholesky')) stop(
      'covtransform must be either "rawcorr", "rawcorr_indep", or "cholesky"')
  }
  
  standata$matsetup <- apply(ctm$modelmats$matsetup[,-1],c(1,2),as.integer,.drop=FALSE) #remove parname and convert to int
  standata$matvalues <- apply(ctm$modelmats$matvalues,c(1,2),as.numeric)
  standata$nmatrices <- as.integer(nmatrices)
  standata$matrixdims <- ctm$modelmats$matrixdims
  standata$nrowmatsetup <- as.integer(nrow(ctm$modelmats$matsetup))
  
  standata$sdscale <- array(as.numeric(sdscale),dim=length(sdscale))
  
  standata$approxct <- 0L
  if(!is.null(ctm$approxct)) standata$approxct <- as.integer(ctm$approxct)
  
  standata$taylorheun <- 0L
  if(!is.null(ctm$taylorheun)) standata$taylorheun <- as.integer(ctm$taylorheun)
  
  #fixed hyper pars #ow disabled
  # if(!is.null(ctm$fixedrawpopchol)) {
  #   standata$fixedrawpopmeans = array(ctm$fixedrawpopmeans)
  #   standata$fixedrawpopchol= ctm$fixedrawpopchol
  # }
  # standata$fixedhyper <- as.integer(ifelse(is.null(ctm$fixedrawpopchol),0,1))
  # if(is.null(ctm$fixedrawpopchol)) {
  #   standata$fixedrawpopmeans = array(0,dim = 0)
  #   standata$fixedrawpopchol = matrix(0,0,0)
  # }
  # standata$fixedsubpars <- as.integer(!is.null(ctm$fixedsubpars))
  # if(!is.null(ctm$fixedsubpars)) standata$fixedindparams <- 
  #   ctm$fixedsubpars else standata$fixedindparams <-array(0,dim=c(0,0))
  
  
  standata$idmap <- idmap
  standata$popcovn=1000L
  standata$llsinglerow=0L
  # standata$doonesubject=0L
  if(any(ctm$manifesttype==1)){
    standata$nJyfinite <- as.integer(ctm$n.latent)
    standata$Jyfinite <- array(as.integer(1:ctm$n.latent))
  }
  # } else{
  #   standata$nJyfinite <- 0L
  #   standata$Jyfinite <- integer()
  # }
  
  if(!is.null(ctm$TIpredAuto) && ctm$TIpredAuto %in% c(1L,TRUE)) standata$TIpredAuto <- 1L else standata$TIpredAuto <- 0L
  
  mc=c(ctStanMatricesList()$all)#base,ctStanMatricesList()$jacobian)
  ms=data.frame(standata$matsetup)
  ms=ms[order(ms$param),]
  standata$whenmat <- array(0L,dim=c(max(mc),5)) #whenmat contains 0's when matrix isn't computed, 1's when it is. 'when 5' is indvaryig.
  
  for(mi in mc){
    mrows=which(ms$matrix==mi)
    for(wheni in 1:4){
      standata$whenmat[mi,wheni] = as.integer(any(
        ms$when[mrows] %in% c(wheni,100))) #100 for PARS, find a better approach to work out exactly when PARS are needed
    }
    standata$whenmat[mi,5] = as.integer(any(
      ms$param[mrows] > 0 & ms$when[mrows] %in% c(0,100) & (ms$indvarying[mrows] > 1 | ms$tipred[mrows] > 0) ))
  }
  rownames(standata$whenmat)[mc] <- names(mc)

  if(optimize){ #indvarying elements need to be updated dynamically
    for(ri in 1:nrow(standata$whenmat)){
      if(standata$whenmat[ri,5]>0){ #if any of this matrices pars are indvarying
        standata$whenmat[mi,
          switch(rownames(standata$whenmat)[ri],PARS=5,T0MEANS=5,LAMBDA=4,DRIFT=2,DIFFUSION=2,MANIFESTVAR=4,MANIFESTMEANS=4,CINT=2,T0VAR=5,TDPREDEFFECT=3,
            JAx=2,Jtd=3,Jy=4,asymCINT=2,asymDIFFUSIONcov=2,DIFFUSIONcov=2,MANIFESTcov=4,T0cov=5)] <- 1L
        # standata$whenmat[ri,5] <- 0 #set indvarying to zero because updating dynamically -- don't do this not all are dynamically updated.
      }
    }
  }

  
  #this PARS when = 100 thing is annoyinh, improve it...
  standata$whenvecp <- array(0L, c(2,standata$nparams)) #whenvecp contains 0's for unchanging pars
  standata$whenvecp[1,] <- as.integer(1:standata$nparams) #base parameters
  standata$whenvecp[2,ms$param[ms$when %in% c(0,100) & ms$copyrow <1 & (ms$tipred > 0 | ms$indvarying > 0) & ms$param > 0]] <- 
    as.integer(ms$param[ms$when %in% c(0,100) & ms$copyrow <1 & (ms$tipred > 0 | ms$indvarying > 0) & ms$param > 0])
  # standata$whenvecp[3,] <- as.integer(1:ncol(standata$whenvecp))
  
  standata$whenvecs <- array(0L,dim=c(6,standata$nlatentpop)) #when do we need to compute transformed states?
  for(wheni in 1:4){ #whenvecs specifies array of when values for the state transform vector -- 1 to compute, 0 not
    standata$whenvecs[wheni,ms$param[ms$when %in% c(wheni,100) & ms$copyrow <=0 & ms$param > 0]] <- 
      as.integer(ms$param[ms$when  %in% c(wheni,100) & ms$copyrow <=0 & ms$param > 0]) #100 is for PARS - needed everywhere.
  }
  
  # why was this in the code? when do we need the below line... ind varying based on states?
  # if(standata$intoverpop==1 && standata$nlatentpop > standata$nlatent){
  #   standata$whenvecs[5,ms$param[ms$when==0 & ms$param > 0 & ms$copyrow < 1 & (ms$indvarying > 0 | ms$tipred > 0)]] <- 
  #     as.integer(ms$param[ms$when==0 & ms$param > 0 & ms$copyrow < 1 & (ms$indvarying > 0 | ms$tipred > 0)])
  # }
  # standata$whenvecs[6,] <- as.integer(1:ncol(standata$whenvecs))
  
  #special matrix adjustments
  standata$whenmat[mc[names(mc) == 'asymDIFFUSIONcov'],] <- 
    apply(standata$whenmat[ mc[names(mc) %in% c('DIFFUSION','DRIFT')],],2,max)
  standata$whenmat[mc[names(mc) == 'asymCINT'],] <- 
    apply(standata$whenmat[mc[names(mc) %in% c('CINT','DRIFT')],],2,max)
  standata$whenmat[mc[names(mc) == 'DIFFUSIONcov'],] <- standata$whenmat[mc[names(mc) == 'DIFFUSION'],]
  standata$whenmat[mc[names(mc) == 'T0cov'],] <- standata$whenmat[mc[names(mc) == 'T0VAR'],]
  standata$whenmat[mc[names(mc) == 'MANIFESTcov'],] <- standata$whenmat[mc[names(mc) == 'MANIFESTVAR'],]
  
  standata$statedep[31:33] <- standata$statedep[c(4,5,8)]
  
  #laplace priors
  standata$laplaceprior <- array(rep(0L,standata$nparams))
  standata$laplacetipreds <- 0L
  if(!is.null(ctm$laplaceprior)){
    if('tipreds' %in% ctm$laplaceprior) standata$laplacetipreds <- 1L
    ms <- data.frame(standata$matsetup)
    standata$laplaceprior[
      ms$param[
        ms$matrix %in% ctStanMatricesList()$all[names(ctStanMatricesList()$all) %in% ctm$laplaceprior] & 
          ms$param > 0 & 
          ms$row!=ms$col & 
          ms$when==0 & 
          ms$copyrow<1]
    ] <- 1L
  }
  standata$laplaceprioronly <- ifelse(is.null(ctm$laplaceprioronly),0L,as.integer(ctm$laplaceprioronly))
  
  
  #CINT non zero
  ms <- data.frame(standata$matsetup)
  CINTnonzero <- c()#1:standata$nlatent
  for(i in 1:standata$nlatent){
    ri=which(ms$matrix %in% ctStanMatricesList()$all[names(ctStanMatricesList()$all) %in% 'CINT'] & ms$row %in% i)
    if(ms$param[ri]==0 && standata$matvalues[ri,'value']==0) next else CINTnonzero <- c(CINTnonzero,i)
  }
  standata$CINTnonzero <- array(as.integer(CINTnonzero))
  standata$CINTnonzerosize <- length(CINTnonzero)
  
  #JAx different from DRIFT?
  standata$JAxDRIFTequiv <- 1L
  ml <- listOfMatrices(ctm$pars)
  for(i in 1:nrow(ml$JAx)){
    for(j in 1:ncol(ml$JAx)){
      if(i <= nrow(ml$DRIFT) && j <= nrow(ml$DRIFT)){ #check drift equivalence
        # if(i==j && !ctm$continuoustime && ml$JAx[i,j] %in% 1) ml$JAx[i,j] <- 
        if(ml$JAx[i,j] != ml$DRIFT[i,j])   standata$JAxDRIFTequiv <- 0L
      }
      if(i > nrow(ml$DRIFT) || j > nrow(ml$DRIFT)){ #check drift equivalence
        if(i != j && !ml$JAx[i,j] %in% 0)   standata$JAxDRIFTequiv <- 0L
        if(i == j && !ml$JAx[i,j] %in% ifelse(ctm$continuoustime,0,1))   standata$JAxDRIFTequiv <- 0L
      }
    }
  }
  
  standata$recompile <- as.integer(ctm$recompile)
  standata$nsubsets <- 1L
  
  return(standata)
}
