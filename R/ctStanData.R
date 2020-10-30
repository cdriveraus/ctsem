ctStanData <- function(ctm, datalong,optimize,derrind='all'){
  
  if(is.null(datalong[[ctm$timeName]]) && ctm$continuoustime == FALSE) {
    datalong <- data.frame(datalong)
    datalong[ctm$timeName] <- 1:nrow(datalong)
  }
  datalong <- datalong[,c(ctm$timeName,ctm$subjectIDname,
    ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames)]
  #start data section
  # if('data.table' %in% class(datalong) || 'tbl' %in% class(datalong)) 
  
  
  nsubjects <- length(unique(datalong[, ctm$subjectIDname])) 
  
  #create random effects indices for each matrix
  
  mats <- ctStanMatricesList()
  
  #simply exponential?
  driftdiagonly <- ifelse(all(!is.na(ctm$pars$value[ctm$pars$matrix == 'DRIFT' & ctm$pars$row != ctm$pars$col]) &
      all(ctm$pars$value[ctm$pars$matrix == 'DRIFT' & ctm$pars$row != ctm$pars$col] == 0) ), 1, 0)
  
  
  ###data checks
  
  if(any(!c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames) %in% colnames(datalong))) stop(paste0('
      variables: ', paste0(c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames)[
        which(!c(ctm$manifestNames,ctm$TDpredNames,ctm$TIpredNames) %in% colnames(datalong))], ', '),' not in data'))
  
  if (!(ctm$subjectIDname %in% colnames(datalong))) stop(paste('id column', (ctm$subjectIDname), "not found in data"))
  
  
  if(ctm$nopriors==FALSE){
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
  # browser()
  original <- unique(datalong[,ctm$subjectIDname])
  datalong <- makeNumericIDs(datalong,ctm$subjectIDname,ctm$timeName)
  new <- unique(datalong[,ctm$subjectIDname])
  idmap <- data.frame(original, new)
  
  
  
  
  
  #linearity checks
  
  #force nonlinear due to cran cuts:
  ctm$nlcontrol$nldynamics=1L
  
  # if(sum(sapply(ctm$modelmats$calcs[!names(ctm$modelmats$calcs) %in% c('jacobian','measurement')],length)) > 0 || 
  #     any(ctm$modelmats$matsetup$when %in% c(1,2,3)) ||
  #     length(ctm$modelmats$calcs$jacobian) - sum(grepl('sJy[',unlist(ctm$modelmats$calcs$jacobian),fixed=TRUE)) > 0 #non measurement jacobians
  #     ){
  #   if(ctm$nlcontrol$nldynamics == FALSE) warning('Linear model requested but nonlinear model specified! May be a poor approximation') else ctm$nlcontrol$nldynamics <- TRUE 
  # }
  
  
  if( (nrow(ctm$t0varstationary) + nrow(ctm$t0meansstationary)) >0 && 
      length(c(ctm$modelmats$calcs$driftcint, ctm$modelmats$calcs$diffusion)) > 0) message('Stationarity assumptions based on initial states when using non-linear dynamics')
  
  
  
  if(ctm$nlcontrol$nldynamics == 'auto') ctm$nlcontrol$nldynamics <- FALSE
  
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
        message(paste0('Missingness in TIpreds - sampling ', sum(is.na(tipreds)),' values'))
        tipreds[is.na(tipreds)] = 99999
      }
      if(optimize){
        message(paste0('Missingness in TIpreds - single imputing ', sum(is.na(tipreds)),'  NA\'s to allow optimization -- TI predictor effect estimates will be overly confident.'))
        tipreds[is.na(tipreds)] = 0
        
        meandat <- data.table((datalong))[ , lapply(.SD, function(x) 
          mean(x,na.rm=TRUE)) , 
          by=c("id")]
        sddat <- data.table((datalong))[ , lapply(.SD, function(x) 
          sd(x,na.rm=TRUE)) , 
          by=c("id")]
        sddat<-sddat[,!colnames(sddat) %in% ctm$subjectIDname,with=FALSE]
        meandat <- meandat[,apply(meandat,2,sd,na.rm=TRUE) > 1e-4,with=FALSE]
        sddat <- sddat[,apply(sddat,2,sd,na.rm=TRUE) > 1e-4,with=FALSE]
        colnames(sddat) <- paste0('sd_',colnames(sddat))
        meandat <- cbind(meandat,sddat)
        
        for(i in 1:ctm$n.TIpred){
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
          # plot(c(meandat[,ctm$TIpredNames[i],with=FALSE])[[1]],predict(lmf),main=ctm$TIpredNames[i])
          tipreds[is.na(tipreds[,1]),1] <- predict(lmf)[is.na(tipreds[,1])]
        }
      }
    }
  }
  datalong[,ctm$manifestNames][is.na(datalong[,ctm$manifestNames])]<-99999 #missing data
  
  
  standata <- list(
    Y=cbind(as.matrix(datalong[,ctm$manifestNames])),
    subject=as.integer(datalong[,ctm$subjectIDname]),
    time=datalong[,ctm$timeName], #not used in model but used elsewhere
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
      message('Missingness in TDpreds! Replaced by zeroes...')
      tdpreds[is.na(tdpreds)] <-0 ## rough fix for missingness
    }
  }
  if(ctm$n.TDpred ==0) tdpreds <- matrix(0,standata$ndatapoints,0) #standata$ndatapoints,
  standata$tdpreds=array(as.matrix(tdpreds),dim=c(nrow(tdpreds),ncol(tdpreds)))
  
  #subset selection
  if(is.null(ctm$dokalmanrows)) standata$dokalmanrowsdata <- 
    rep(1L, standata$ndatapoints) else standata$dokalmanrowsdata <- as.integer(ctm$dokalmanrows)
  standata$dokalmanpriormodifier = sum(standata$dokalmanrows)/standata$ndatapoints
  
  standata<-c(standata, 
    list(
      nsubjects=as.integer(nsubjects),
      nmanifest=as.integer(ctm$n.manifest),
      nlatentpop = ctm$nlatentpop,
      nldynamics=as.integer(ctm$nlcontrol$nldynamics),
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
      nopriors=as.integer(ctm$nopriors),
      savescores=0L,
      savesubjectmatrices=0L,
      nsJAxfinite = ifelse(ctm$recompile,0L,length(ctm$sJAxfinite))
    ))
  
  if(ctm$n.TIpred == 0) tipreds <- array(0,c(0,0))
  standata$tipredsdata <- as.matrix(tipreds)
  standata$nmissingtipreds <- as.integer(length(tipreds[tipreds== 99999]))
  
  # browser()
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
  
  
  standata$sJAxfinite <- ctm$sJAxfinite
  
  standata$sJycolindexsize <- 1L
  standata$sJycolindex <- array(1L)
  
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
  for(mati in names(mats$base)){
    if( (!ctm$intoverpop && any(ctm$pars$indvarying[ctm$pars$matrix==mati])) || 
        (ctm$n.TIpred >0 && (
          any(unlist(ctm$pars[ctm$pars$matrix==mati,paste0(ctm$TIpredNames,'_effect')])) || 
            any(ctm$pars$matrix==mati & grepl('[',ctm$pars$param,fixed=TRUE))  )
        )) subindex <- 1 else subindex <- 0
        subindices[[mati]] <- subindex
  }
  
  if(ctm$stationary || nrow(ctm$t0varstationary) > 0) subindices$T0VAR  <- 
    max(c(subindices$T0VAR,subindices$DRIFT,subindices$DIFFUSION))
  
  if(ctm$stationary || nrow(ctm$t0meansstationary) > 0) subindices$T0MEANS <- 
    max(c(subindices$T0MEANS,subindices$DRIFT,subindices$CINT))
  
  
  
  standata$subindices <- as.integer(unlist(subindices))
  
  #state dependence
  statedep <- rep(0L,10)
  lhscalcs <- sapply(unique(unlist(ctm$modelmats$calcs)),function(x) gsub('=.*','',x))
  for(i in 1:length(statedep)){
    matname <- try(names(ctStanMatricesList()$all[ctStanMatricesList()$all %in% i]),silent=TRUE)
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
  standata$choleskymats<- ifelse(ctm$covmattransform=='unconstrainedcorr',0L,1L)
  if(!ctm$covmattransform %in% c('unconstrainedcorr','cholesky')) stop('covtransform must be either "unconstrainedcorr" or "cholesky"')
  
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
  
  #fixed hyper pars
  if(!is.null(ctm$fixedrawpopchol)) {
    standata$fixedrawpopmeans = array(ctm$fixedrawpopmeans)
    standata$fixedrawpopchol= ctm$fixedrawpopchol
  }
  standata$fixedhyper <- as.integer(ifelse(is.null(ctm$fixedrawpopchol),0,1))
  if(is.null(ctm$fixedrawpopchol)) {
    standata$fixedrawpopmeans = array(0,dim = 0)
    standata$fixedrawpopchol = matrix(0,0,0)
  }
  standata$fixedsubpars <- as.integer(!is.null(ctm$fixedsubpars))
  if(!is.null(ctm$fixedsubpars)) standata$fixedindparams <- 
    ctm$fixedsubpars else standata$fixedindparams <-array(0,dim=c(0,0))
  
  
  standata$idmap <- idmap
  standata$popcovn=1000L
  standata$llsinglerow=0L
  standata$doonesubject=0L
  if(any(ctm$manifesttype==1)){
    standata$nJyfinite <- as.integer(ctm$n.latent)
    standata$sJyfinite <- array(as.integer(1:ctm$n.latent))
  } else{
    standata$nJyfinite <- 0L
    standata$sJyfinite <- integer()
  }
  
  if(!is.null(ctm$TIpredAuto) && ctm$TIpredAuto %in% c(1L,TRUE)) standata$TIpredAuto <- 1L else standata$TIpredAuto <- 0L
  
  mc=c(ctStanMatricesList()$all)#base,ctStanMatricesList()$jacobian)
  ms=data.frame(standata$matsetup)
  ms=ms[order(ms$param),]
  standata$whenmat <- array(0L,dim=c(max(mc),5)) #whenmat contains 0's when matrix isn't computed, 1's when it is. 'when 5' is indvaryig.
  for(mi in mc){
    mrows=which(ms$matrix==mi)
    for(wheni in 1:4){
      standata$whenmat[mi,wheni] = as.integer(any(
        ms$when[mrows]==wheni))
    }
    standata$whenmat[mi,5] = as.integer(any(
      ms$param[mrows] > 0 & ms$when[mrows] <=0 & (ms$indvarying[mrows] == 1 | ms$tipred[mrows] > 0) ))
  }
  rownames(standata$whenmat)[mc] <- names(mc)
  
  standata$whenvecp <- array(0L, c(3,standata$nparams)) #whenvecp contains 0's for unchanging pars, 1's for changing pars
  standata$whenvecp[1,] <- as.integer(1:standata$nparams) #base parameters
  standata$whenvecp[2,ms$param[ms$when == 0 & ms$copyrow <1 & (ms$tipred > 0 | ms$indvarying > 0) & ms$param > 0]] <- 
    as.integer(ms$param[ms$when == 0 & ms$copyrow <1 & (ms$tipred > 0 | ms$indvarying > 0) & ms$param > 0])
  standata$whenvecp[3,] <- as.integer(1:ncol(standata$whenvecp))
  
  standata$whenvecs <- array(0L,dim=c(6,standata$nlatentpop)) #when do we need to compute transformed states?
  for(wheni in 1:4){ #whenvecs specifies array of when values for the state transform vector -- 1 to compute, 0 not
    standata$whenvecs[wheni,ms$param[ms$when == wheni & ms$copyrow <=0 & ms$param > 0]] <- 
      as.integer(ms$param[ms$when == wheni & ms$copyrow <=0 & ms$param > 0])
  }
  # browser()
  #when do we need the below line... ind varying based on states?
  if(standata$intoverpop==1 && standata$nlatentpop > standata$nlatent){
    standata$whenvecs[5,ms$param[ms$when==0 & ms$param > 0 & ms$copyrow < 1 & (ms$indvarying > 0 | ms$tipred > 0)]] <- 
      as.integer(ms$param[ms$when==0 & ms$param > 0 & ms$copyrow < 1 & (ms$indvarying > 0 | ms$tipred > 0)])
  }
  standata$whenvecs[6,] <- as.integer(1:ncol(standata$whenvecs))
  
  #special matrix adjustments
  standata$whenmat[mc[names(mc) == 'asymDIFFUSION'],] <- 
    apply(standata$whenmat[ mc[names(mc) %in% c('DIFFUSION','DRIFT')],],2,max)
  standata$whenmat[mc[names(mc) == 'asymCINT'],] <- 
    apply(standata$whenmat[mc[names(mc) %in% c('CINT','DRIFT')],],2,max)
  standata$whenmat[mc[names(mc) == 'DIFFUSIONcov'],] <- standata$whenmat[mc[names(mc) == 'DIFFUSION'],]
  
  return(standata)
}
