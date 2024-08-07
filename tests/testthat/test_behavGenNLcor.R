if(identical(Sys.getenv("NOT_CRAN"), "true")& !(.Platform$OS.type=="windows" && 
    R.version$major %in% 4 && as.numeric(R.version$minor) >= 2 &&
    unlist(utils::packageVersion('rstan'))[2] < 25) ){
  library(ctsem)
  library(testthat)
  library(data.table)
  cores=2
  
  
  test_that("behavGenNLcor", {
    
    
    # data gen ----------------------------------------------------------------
    
    dt <- .1
    maxtime <- 20
    steps=maxtime/dt
    nsubjects <- 50
    nmanifest <- 2
    drift <- -.2 #relaxation of changes in slope
    slopecrosseffect <- -.02 #effect of level on slope
    diffusion <- -1
    t0corz <- 2
    zyg_t0corz <- .3
    age_corz <- -1
    agesq_corz <-.1
    zyg_age_corz <- .2
    initialslopesd <- .1
    slopeint <- 1
    diffusionage <- -.1
    ressd <- .2
    ressdByEarlyTest <- .3
    
    y <- matrix(NA,ncol=nmanifest,nrow=steps*nsubjects)
    zyg <- sample(0:1,nsubjects,replace=TRUE)
    for(subi in 1:nsubjects){
      corz <- (t0corz+zyg[subi]*zyg_t0corz)
      cholm <- t(chol(matrix(c(1,inv_logit(corz)*2-1,inv_logit(corz)*2-1,1),2,2)))
      slope <- slopeint + cholm %*% rnorm(2,0,initialslopesd)
      age <- 0
      for(stepi in 2:steps){
        age <- age + dt
        i <- (subi-1)*steps+stepi
        if(stepi==2) y[i-1,] <- rnorm(1,0,.01) #first time point
        corz <- (t0corz+zyg[subi]*zyg_t0corz) + age_corz*age +
          zyg_age_corz * age * zyg[subi]+
          agesq_corz*age^2
        # if(subi==1) print(paste0('age = ',age,'  corz= ', corz))
        cholm <- t(chol(matrix(c(1,(inv_logit(corz)*2-1)*.99,(.99*inv_logit(corz)*2-1),1),2,2)))
        slope <- slope + (slopecrosseffect*y[i-1,] + slope * drift) * dt +
          cholm %*% rnorm(2,0,exp(diffusion+diffusionage*age)) * dt
        y[i,] <- y[i-1,] + slope * dt
      }
    }
    
    #set zygosity coding to mz = -1 and dz = +1
    zyg = -(zyg*2-1)
    
    dat <- data.table(time=rep(seq(dt,maxtime,dt),nsubjects),
      id=factor(rep(1:nsubjects,each=steps)),adjValue=y,zyg=rep(zyg,each=steps))
    
    colnames(dat)[3:4] <- c('adjValue1r','adjValue2r')
    
    dat[, c('adjValue1','adjValue2'):= list(c(scale(adjValue1r)),c(scale(adjValue2r))), by=time] #scaled vars
    # dat[, c('adjValue1','adjValue2'):= list(c((adjValue1r)),c((adjValue2r))), by=time] #unscaled vars
    
    # raw cor -----------------------------------------------------------------
    
    # dat<-data.table(dat)
    dat[, "latentcor":= list(cor(adjValue1r,adjValue2r)), by=list(time,zyg)] #correlation by time
    dat[, "latentcors":= list(cor(adjValue1,adjValue2)), by=list(time,zyg)] #correlation on scaled (as check)
    dat[, "latentsds":= list(sd(c(adjValue1,adjValue2))), by=list(time,zyg)] #correlation on scaled (as check)
    
    # plot(unique(dat[zyg == 1,time]),unique(dat[zyg == 1,latentcor]),type='l',ylim=range(dat$latentcor))
    # points(unique(dat[zyg == -1,time]),unique(dat[zyg == -1,latentcor]),type='l',col='blue')
    
    # points(dat[id==1,time],dat[id==1,latentcors],type='l',lty=2,col='red')
    
    
    Nsubjects <- 50
    Nobs <- 10
    dat<-dat[as.numeric(id)<=Nsubjects,]
    vars <- c('adjValue1','adjValue2')
    dat$level1 <- dat$adjValue1
    dat$level2 <- dat$adjValue2
    
    #measurement error
    dat$earlyTest <- as.numeric(dat$time < 4)
    dat$adjValue1 <- rnorm(nrow(dat), dat$adjValue1, ressd+ressdByEarlyTest * dat$earlyTest)
    dat$adjValue2 <- rnorm(nrow(dat), dat$adjValue2, ressd+ressdByEarlyTest * dat$earlyTest)
    
    #compute observed correlations after measurement error
    dat[, "obscor":= list(cor(adjValue1r,adjValue2r)), by=list(time,zyg)] #observed correlation by time
    dat[, "obscors":= list(cor(adjValue1,adjValue2)), by=list(time,zyg)] #observed correlation on scaled (as check)
    dat[, "obssds":= list(sd(c(adjValue1,adjValue2))), by=list(time,zyg)] #observed sd on scaled (as check)
    
    #select random data points but include missing first obs (so t0var is coherent)
    samp <- sort(unique(c(sample((1:nrow(dat))[dat$time > .5],Nsubjects*Nobs), which(dat$time==min(dat$time)))))
    ltsData <- dat[samp,]
    ltsData[time %in% min(time),c('adjValue1','adjValue2')] <- NA #na initial obs
    ltsData[time %in% min(time),time:= 0.5] #set first obs to occur just before earliest observation for everyone
    
    
    #create centered age
    ltsData[,centeredAge:=time-mean(time)]
    
    
    # modelling ---------------------------------------------------------------
    agevarying=2
    
    latentNames=c('level1',"level2",'slope1','slope2')
    
    PARS <- c('baselevel','baseslope' ) #no extra pars unless added lower
    
    basecor <- '2/(1 + exp(-basecor)) - 1'
    
    levelonslope <- 'levelonslope'
    
    if(agevarying  %in% c(0,1)){
      diffsd <- "diffsd"
      drift <- "driftbase"
      ressd <- 'ressd'
    }
    
    if(agevarying == 0){  #no change in correlation
      diffcor <- basecor
    }
    
    if(agevarying == 1){ #simplest change in correlation over time
      diffcor <- "2/(1 + exp(-(basecor + diffcor))) - 1"
    }
    
    if(agevarying %in% c(0,1,2))   PARS = c(
      PARS,
      c('basecor ||||zyg',
        'diffcor||||zyg'))
    
    if(agevarying == 2){ #more flexibly changing correlation and slope variance over time
      diffcor <- "2/(1 + exp(-(basecor + diffcor+ diffcorByAge*tdpreds[rowi,2]+ diffcorByAgeSq*(tdpreds[rowi,2])^2))) - 1" #
      diffsd <- "log1p(exp(diffsd+ diffsdByAge*tdpreds[rowi,2] + diffsdByAgesq*tdpreds[rowi,2]^2))+1e-5" #
      drift <- "-log1p(exp(-(drift+driftByAge*tdpreds[rowi,2]+driftByAgesq*tdpreds[rowi,2]^2)))-1e-6" #
      ressd <- 'log1p_exp(ressd+ressdByAge*tdpreds[rowi,2]+ressdByEarlyTest*tdpreds[rowi,1])'
      levelonslope <- '-log1p_exp(levelonslope+levelonslopeByAge*tdpreds[rowi,2]+levelonslopeByAgeSq*tdpreds[rowi,2]^2)'
      
      # latentNames <- c(latentNames,"agetrend")   #add trend
    }
    
    #base matrices
    T0VAR = matrix(c(
      "t0levelsd",0,0,0,
      basecor,"t0levelsd", 0,0,
      0,0,'baseslopesd',0,
      0,0,diffcor,'baseslopesd'),byrow=TRUE,4,4)
    
    DRIFT=matrix(c(
      0,0,1,0,
      0,0,0,1,
      levelonslope,0,drift,0,
      0,levelonslope,0,drift),byrow=TRUE,4,4)
    
    
    DIFFUSION=matrix(c(
      0,0,0,0,
      0,0,0,0,
      0,0,diffsd,0,
      0,0,diffcor,diffsd),byrow=TRUE,4,4)
    
    CINT=c(0,0,'cint1||FALSE',"cint1||FALSE") #disable individual variation in par
    
    T0MEANS=c('baselevel',"baselevel",'baseslope','baseslope') #disable individual variation in pars
    
    LAMBDA = matrix(c(
      1,0,0,0,
      0,1,0,0),byrow=TRUE,2,4)
    
    
    if(agevarying==2){
      
      PARS=c(PARS,
        "diffcorByAge||||zyg","diffcorByAgeSq||||zyg",
        "diffsd","diffsdByAge",'diffsdByAgesq',
        'drift','driftByAge','driftByAgesq',
        "ressd","ressdByAge",'ressdByEarlyTest',
        'levelonslope','levelonslopeByAge','levelonslopeByAgeSq'
      )
    }
    
    m1 <- ctModel(type = 'stanct',tipredDefault = FALSE,
      manifestNames = c('adjValue1',"adjValue2"),
      latentNames=latentNames,
      TIpredNames = 'zyg',
      TDpredNames = c('earlyTest','centeredAge'),#c("test1","test2","test3","test4","test5"),
      T0VAR = T0VAR,
      DRIFT=DRIFT,
      DIFFUSION=DIFFUSION,
      CINT=CINT,
      T0MEANS=T0MEANS,
      TDPREDEFFECT = 0,
      LAMBDA = LAMBDA,
      MANIFESTMEANS=0,
      # MANIFESTMEANS=matrix(c('td1*tdpreds[rowi,1]+td2*tdpreds[rowi,2]+td3*tdpreds[rowi,3]+td4*tdpreds[rowi,4]+td5*tdpreds[rowi,5]',
      # 'td1*tdpreds[rowi,1]+td2*tdpreds[rowi,2]+td3*tdpreds[rowi,3]+td4*tdpreds[rowi,4]+td5*tdpreds[rowi,5]'),2,1),
      MANIFESTVAR=c(ressd,0,
        'ressdcor',ressd),
      PARS = PARS
    )
    
    m1$pars$sdscale[!grepl('drift',m1$pars$param)] <- 10
    
    # fit ---------------------------------------------------------------------
    f <- ctStanFit(datalong = ltsData, ctstanmodel = m1,cores=cores,saveComplexPars = T,fit=T)#,optimcontrol=list(stochastic=F,carefulfit=F),init=rep(0,30)
    
    testthat::expect_equivalent(class(f),'ctStanFit')
    
    
    
  })
}





