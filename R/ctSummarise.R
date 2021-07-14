ctSummarise<-function(sf,name='ctSummary',cores=2, times=seq(0,10,.1),quantiles=c(.025,.5,.975),
  folder = 'ctSummary/',nsamples=200,lags=1:3,
  ctCheckFit = TRUE, ctStanPlotPost=TRUE, ctStanTIpredEffects = TRUE,IndDifCorrelations=TRUE,
  latents=1:sf$ctstanmodelbase$n.latent, manifests=1:sf$ctstanmodelbase$n.manifest){
  
  Sig. <- NULL
  
  if(class(sf)=='ctStanFit'){ #avoid overwriting plots if error!
    # library(data.table)
    oldwd = getwd()
    on.exit(add = TRUE,expr = setwd(oldwd))
    dir.create(file.path(oldwd,folder))
    setwd(file.path(oldwd,folder))
    
    mp <- options()$max.print
    on.exit(expr = options(max.print=mp),add=TRUE)
    options(max.print=99999)
    
    
    ydat=standatatolong(sf$standata,origstructure = TRUE,ctm=sf$ctstanmodelbase)
    
    sm <- sf$ctstanmodelbase 
    summ=summary(sf)
    if(is.null(sf$generated)) sf <- ctStanGenerateFromFit(fit = sf,nsamples = nsamples,fullposterior = FALSE,cores = cores)
    cp <- ctStanContinuousPars(sf)
    
    n<-sapply(unique(ydat[[sm$subjectIDname]]),function(x) sum(ydat[[sm$subjectIDname]] %in% x))
    whichsubfull <- unique(ydat[[sm$subjectIDname]])[order(n,decreasing = TRUE)][1:(min(30,length(n)))]
    # whichsubfull <- match(whichsubfull,sf$setup$idmap[,1]) #set to numeric
    
    nl <- sm$n.latent
    nm <- sm$n.manifest
    
    viscovf <- function(mat){
      matname <- mat
      if(mat %in% 'MANIFESTcov') {
        matname <- 'Residual'
        varnames <- sf$ctstanmodelbase$manifestNames
      } else varnames <- sf$ctstanmodelbase$latentNames
      #correlation ci's
      diffusion=sf$stanfit$transformedpars[[paste0('pop_',mat)]]
      if(mat %in% 'T0cov' && sf$standata$nindvarying > 0 && sf$standata$intoverpop){
        diffusion = diffusion[,latents,latents,drop=FALSE]
      }
      diffindices <- unique(c(which(diffusion[1,,] != 0,arr.ind = TRUE)))
      if(length(diffindices) > 0){
      diffusion <- plyr::aaply(diffusion,1,function(x) matrix(x[diffindices,diffindices,drop=FALSE],length(diffindices),length(diffindices)))
      dimnames(diffusion) <- list(NULL,varnames[diffindices],varnames[diffindices])
      dcor=plyr::aaply(diffusion,1, function(x) cov2cor(x))
      dcorq=plyr::aaply(dcor,c(2,3),quantile,probs=c(.025,.5,.975))
      dcorm <- as.data.table(dcorq)
      dcorm <- dcast(dcorm,'V1+V2 ~ V3')
      print(corplotmelt(meltcov(dcorq[,,1]),title = paste0(gsub('cov','',matname),' corr. 2.5% quantile')))
      print(corplotmelt(meltcov(dcorq[,,2]),title = paste0(gsub('cov','',matname),' corr. 50% quantile')))
      print(corplotmelt(meltcov(dcorq[,,3]),title = paste0(gsub('cov','',matname),' corr. 97.5% quantile')))
      options(max.print=100000)
      cat(gsub('cov','',matname),' correlation matrix confidence intervals')
      print(dcorm,nrows = 10000)
      }
    }
    sink(paste0(name,'_corrCI.txt'))
    pdf(paste0(name,'_VisualParMats.pdf'))
    lapply(c('MANIFESTcov','DIFFUSIONcov','asymDIFFUSIONcov','T0cov'),function(x) try(viscovf(x)))
    sink()
    # print(corplotmelt(meltcov(cov2cor(cp$MANIFESTcov[manifests,manifests,drop=FALSE])),title =  'Residual correlations'))
    # print(corplotmelt(meltcov(cov2cor(cp$DIFFUSIONcov[latents,latents,drop=FALSE])),title =  'Random latent change correlations'))
    # print(corplotmelt(meltcov(cov2cor(cp$asymDIFFUSION[latents,latents,drop=FALSE])),title =  'Asymptotic latent correlations'))
    # print(corplotmelt(meltcov(cov2cor(cp$T0cov[latents,latents,drop=FALSE])),title =  'Initial latent correlations'))
    dr=cp$DRIFT[latents,latents,drop=FALSE]
    dr[diag(nrow(dr))==1] <- -dr[diag(nrow(dr))==1]
    print(corplotmelt(meltcov(cov2cor(dr %*% t(dr))),title =  'Std. deterministic change relations'))
    dev.off()
    
    ### GENERATING OUTPUT ### 
    
    # generate summary of output in textfile
    sink(file = paste0(name,'_summary.txt'))
    print(summ)
    cat("\n")
    print("Z > 1.96 Correlations:")
    if(!is.null(summ$rawpopcorr)) print(summ$rawpopcorr[which(abs(summ$rawpopcorr[, "z"]) > 1.96), ])
    if (sf$ctstanmodel$n.TIpred > 0) {
      cat("\n")
      print("Z > 1.96 Covariates (TIPs):")
      print(summ$tipreds[which(abs(summ$tipreds[, "z"]) > 1.96), ])
    }
    cat("\n")
    print("popmeans with 0 NOT within range 2.5% - 97.5%")
    print(
      subset(summ$popmeans, summ$popmeans$`2.5%` > 0 & summ$popmeans$`97.5%` > 0 | summ$popmeans$`2.5%` < 0 & summ$popmeans$`97.5%` < 0)
    )
    cat("\n")
    print('Par matrices')
    print(cp)
    print('T0VAR correlation')
    print(cov2cor(cp$T0cov))
    print('DIFFUSION correlation')
    print(cov2cor(cp$DIFFUSIONcov))
    print('asymptotic DIFFUSION correlation')
    print(cov2cor(cp$asymDIFFUSION))
    print('Residual correlation')
    print(cov2cor(cp$MANIFESTcov))
    sink()
    
    # #Tables output -> makes word document
    # if(!requireNamespace('flextable')) {
    #   warning('flextable package required for word table output')
    # } else {
    #   if(!requireNamespace('officer')) {
    #     warning('officer package required for word table output')
    #   } else{ #do tables
    sitems=c('popmeans','popsd','tipreds','rawpopcorr')
    ftlist <- lapply(sitems,
      function(x){
        if(!is.null(summ[[x]])){
          fti=data.frame(round(summ[[x]] ,2))
          fti=data.frame(Parameter=rownames(fti),fti[,colnames(fti)[!colnames(fti) %in% 'X50.']])
          colnames(fti)[1:5] <- c('Parameter','Mean','SD','2.5%','97.5%')
          return(fti)
        } else return(NULL)
      })
    names(ftlist)=sitems
    
    ftshort <- lapply(c('tipreds','rawpopcorr'), function(x){
      if(!is.null(ftlist[[x]])){
        fts <- ftlist[[x]][apply(ftlist[[x]][,c('2.5%','97.5%')],1,function(x) abs(sum(sign(x)))==2),]
      } else return(NULL)})
    names(ftshort) = c('tipreds','rawpopcorr')
    
    ftlist <- c(ftlist,ftshort)
    ftlist <- ftlist[c(1,2,5,6,3,4)]
    names(ftlist) <- c('Pop. mean parameters','Pop. SD parameters','Significant covariate effects',
      'Significant random effect correlations','All covariate effects','All random effect correlations')
    
    sink(file = paste0(name,'_tables.txt'))
    print(ftlist,row.names=FALSE)
    sink()
    
    #   ftlist <- lapply(ftlist,function(x){
    #     if(!is.null(x)){
    #       x=flextable::autofit(flextable::font(
    #         flextable::flextable(as.data.frame(x))
    #         ,fontname = 'Times New Roman'))
    #     }
    #     return(x)})
    #   flextable::save_as_docx(values=ftlist,path=paste0(name,'_tables.docx'))
    # }}
    
    
    
    # #how many subjects have max observed number of time points
    # tldatNArm <- data.frame(ydat)
    # tldatNArm <- tldatNArm[apply(ydat[,sm$manifestNames],1,function(x) any(!is.na(x))),]
    # n<-sapply(unique(tldatNArm$id),function(x) length(tldatNArm$id[tldatNArm$id %in% x]))
    # whichsubfull <- unique(data$id)[n %in% c(max(n):(max(n)-30))]
    # whichsubfull <- match(whichsubfull,sf$setup$idmap[,1]) #set to numeric
    # sink(file = paste0('WhichSubFull_',name,'.txt'))
    # print("WHICH SUBJECTS OBSERVED AT ALL TIME POINTS")
    # cat("\n")
    # print("number of subject with one or more time points")
    # print(length(n))
    # print("Most time points for a subject")
    # print(max(n))
    # print("IDs of subjects with all time points")
    # print(whichsubfull)
    # print("Average Timepoints per ID")
    # print(sum(n)/(length(n)))
    # cat("\n")
    # print("Timepoints per ID")
    # print(n)
    # sink()
    
    #print Latex equation 
    try(ctModelLatex(sf,folder = './',filename = paste0(name,'_tex_fit'),open = FALSE))
    try(ctModelLatex(sm,folder = './',filename = paste0(name,'_tex_model'),open = FALSE))
    
    # #TI Correlations matrix plot
    #   sub <- ctStanSubjectPars(sf)
    #   pdf('TI_correlations.pdf')
    #   corm=meltcov(cov2cor(cov(
    #     cbind(sf$standata$tipredsdata,sub[1,,order(dimnames(sub)[[3]])]))))
    #   a=corplotmelt(corm, label = 'Corr.')
    #   print(a+ggplot2::labs(title='Individual difference correlations'))
    #   dev.off()
    #   
    
    
    
    #individual difference / ti pred correlations
    
    #subject par correlation quantiles
    # library(plyr)
    # library(ggplot2)
    # browser()
    if(IndDifCorrelations){
      if(sf$standata$ntipred > 0 || sf$standata$nindvarying > 0){
        try({
          sp=ctStanSubjectPars(sf,pointest=FALSE,cores=cores,nsamples = nsamples)
          
          if(sf$standata$ntipred > 0){
            tip <- array(rep(sf$standata$tipredsdata,each=dim(sp)[[1]]),dim=c(dim(sp)[1:2],ncol(sf$standata$tipredsdata)))
            spti=array(c(sp, tip), dim = c(dim(sp)[1], dim(sp)[2], dim(sp)[3]+dim(tip)[3]))
          } else spti <- sp
          
          dn=dimnames(sp)
          dn[[3]] <- c(dn[[3]],sf$ctstanmodelbase$TIpredNames)
          dimnames(spti) <- dn
          
          cornames <- matrix(
            paste0(dimnames(spti)[[3]],'_',rep(dimnames(spti)[[3]],each=length(dimnames(spti)[[3]]))),
            nrow=length(dimnames(spti)[[3]]))
          
          cornames <- cornames[lower.tri(cornames)]
          cm=aaply(spti,1,cor,.drop=FALSE)
          qc=aaply(cm,c(3,2),function(x) round(c(Mean=mean(x),SD=sd(x),
            `2.5%`=quantile(x,probs=c(.025)),`50%`=quantile(x,probs=c(.5)),
            `97.5%`=quantile(x,probs=c(.975)),z=mean(x)/sd(x)),2),.drop=FALSE)
          
          qcsig <- apply(qc[,,c("2.5%.2.5%","97.5%.97.5%"),drop=FALSE],c(1,2),function(x) abs(sum(sign(x)))==2)
          cors <- as.data.table(qc)
          cors <- data.frame(dcast(cors,'V1+V2~V3'))
          colnames(cors) <- c('Par1','Par2','2.5%','50%','97.5%','Mean','SD','z')
          cors <- cors[,c('Par1','Par2','Mean','SD','2.5%','50%','97.5%','z')]
          cors <- cors[paste0(cors$Par1,'_',cors$Par2) %in% cornames,]
          corssig <- cors[apply(cors,1,function(x) sum(sign(as.numeric(x[3:5])))==3),]
          
          #now do similar to point estimate
          spp=ctStanSubjectPars(sf,pointest=TRUE,cores=1)
          
          if(sf$standata$ntipred > 0){
            tip <- array(rep(sf$standata$tipredsdata,each=dim(spp)[[1]]),dim=c(dim(spp)[1:2],ncol(sf$standata$tipredsdata)))
            spti=array(c(spp, tip), dim = c(dim(spp)[1], dim(spp)[2], dim(spp)[3]+dim(tip)[3]))
          } else spti <- spp
          dn=dimnames(spp)
          dn[[3]] <- c(dn[[3]],sf$ctstanmodelbase$TIpredNames)
          dimnames(spti) <- dn
          cm=cor(spti[1,,])
          
          mcor=meltcov(cm)
          msig=data.frame(matrix(as.numeric(qcsig),nrow=nrow(qcsig)))
          dimnames(msig) <- dimnames(cm)
          mcorsig=meltcov(msig)
          mcorsig$sig <- mcorsig$value
          mcorsig$value <- NULL
          mcor <- merge(mcor,mcorsig)
          mcor$`Sig.` <- as.logical(mcor$sig)
          pdf(paste0(name,'_IndDifCorrelations.pdf'))
          print(ggplot(data=(mcor),
            aes_string(x='Var1',y='Var2',fill=('value')))+ #ifelse(groups,NULL,'Group')))+
              geom_tile( width=1,height=1,colour='black')+
              scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                midpoint = 0, limits = c(-1,1), space = "Lab", 
                name='Corr')  +
              geom_point(mapping=aes(alpha=Sig.),stroke=0,size=.6)+
              theme_minimal()+ theme(axis.text.y=element_text(size=8),
                axis.text.x = element_text(angle = 90,size=8)))
          dev.off()
        }) #end try
      }
      
      sink(paste0(name,'_IndDifCorrelations.txt'))
      print(cm,digits=3)
      sink()
      
      
      if(sf$ctstanmodelbase$n.TIpred > 0){
        
        if(sf$ctstanmodelbase$n.TIpred > 1){
          #TIP Effect expectations 2
          ktifull=ctKalmanTIP(sf,kalmanvec='yprior',subject = 1,
            #elementNames=c(sf$ctstanmodelbase$latentNames),
            plot=TRUE,polygonsteps=0)+
            scale_colour_manual(name='Covariate',values=c(1,rep(2:((sm$n.TIpred+1)),each=2)))+
            ggtitle('Expectations')
          
          lktifull=ctKalmanTIP(sf,kalmanvec='etaprior',subject = 1,
            plot=TRUE,polygonsteps=0)+
            scale_colour_manual(name='Covariate',values=c(1,rep(2:((sm$n.TIpred+1)),each=2)))+
            ggtitle('Latent Expectations')
        }
        
        pdf(paste0(name,'_tieffects_expectations.pdf'))
        if(sf$ctstanmodelbase$n.TIpred > 1){
          print(lktifull)
          print(ktifull)
          rm(lktifull);rm(ktifull);
        }
        for(ti in sm$TIpredNames){
          print(plot(ctKalmanTIP(sf,tipreds = ti,subject = whichsubfull[1]),
            kalmanvec='etaprior', plot=FALSE,polygonsteps=0)+scale_color_discrete(name='Covariate')+ggtitle('Latent Expectations'))
          print( plot(ctKalmanTIP(sf,tipreds = ti,subject = whichsubfull[1]),kalmanvec='yprior',
            plot=FALSE,polygonsteps=0)+scale_color_discrete(name='Covariate')+ggtitle('Expectations'))
        }
        dev.off()
      }
    }
    
    
    #subject expectation plots 
    k<-ctKalman(sf,subjects = whichsubfull,realid = TRUE)
    krem<-ctKalman(sf, subjects = whichsubfull,removeObs = TRUE,realid = TRUE)
    pdf(paste0(name,'_subjectexpectations.pdf'))
    print(plot(krem,polygonsteps=FALSE, kalmanvec='etaprior',plot=FALSE)+ggtitle('Latent Expectations Conditional on Covariates'))
    print(plot(k,polygonsteps=FALSE, kalmanvec='etasmooth',plot=FALSE) +ggtitle('Latent Expectations Conditional on All Data'))
    print(plot(krem,polygonsteps=FALSE, kalmanvec='yprior',plot=FALSE)+ggtitle('Expectations Conditional on Covariates'))
    print(plot(k,polygonsteps=FALSE, kalmanvec='ysmooth',plot=FALSE)+ggtitle('Expectations Conditional on All Data'))
    dev.off()
    rm(k);rm(krem)
    
    #Discrete pars plots
    pdf(paste0(name,'_discretepars.pdf'))
    dtpars=ctStanDiscretePars(sf,plot=FALSE,times=times) #,indices=cbind(1:9,i))
    dtparso=ctStanDiscretePars(sf,plot=FALSE,observational = TRUE,times=times) #
    # dtparsc=ctStanDiscretePars(sf,plot=FALSE,observational = TRUE,cov = TRUE,nsamples = 500,times=times) #
    # dtparsoc=ctStanDiscretePars(sf,plot=FALSE,observational = FALSE,cov = TRUE,nsamples = 500,times=times) #
    for(i in 1:sf$ctstanmodelbase$n.latent){
      print(ctStanDiscreteParsPlot(dtpars,indices=cbind(latents,i),quantiles = quantiles))
      print(ctStanDiscreteParsPlot(dtparso,indices=cbind(latents,i),quantiles=quantiles))
      # print(ctStanDiscreteParsPlot(dtparsc,indices=cbind(1:9,i),quantiles = c(.4,.5,.6)))
      # print(ctStanDiscreteParsPlot(dtparsoc,indices=cbind(1:9,i),quantiles = c(.4,.5,.6)))
    }
    dev.off()
    
    #tipred effect plots
    if(ctStanTIpredEffects){
      try({
        if(sf$ctstanmodelbase$n.TIpred >0){
          pdf(paste0(name,'_TIpredEffects.pdf'))
          for(tip in 1:sf$ctstanmodelbase$n.TIpred){
            tieffects <- ctStanTIpredeffects(sf,nsamples = 200,whichTIpreds = tip,timeinterval = 1)
            tieffects$y <- tieffects$y[, #drop unchanging parameters
              apply(tieffects$y[,,2,drop=FALSE],2,function(x) length(unique(x))>1),,drop=FALSE]
            matnames <- unique(unlist(sapply(colnames(tieffects$y),function(x) gsub(pattern = '\\[.*','',x))))
            for(m in matnames){
              tmp <- tieffects
              tmp$y <- tmp$y[,grep(paste0('^',m,'\\['),colnames(tmp$y) ),,drop=FALSE]
              print(ctPlotArrayGG(input = tmp))
            }
          }
          dev.off()
        }
      })
    }
    
    ### CHECKING MODELFIT ###
    # Checkfit function 
    if(ctCheckFit){
      pdf(paste0(name,'_checkfit.pdf'))
      by=sm$timeName
      try(ctCheckFit(fit = sf,data = TRUE,postpred = TRUE,statepred = FALSE,by = by,breaks=2,covplot = TRUE,smooth=TRUE,reg=TRUE))
      try(ctCheckFit(fit = sf,data = TRUE,postpred = TRUE,statepred = TRUE,by = by,breaks=4,covplot = FALSE,smooth=TRUE,reg=TRUE))
      try(ctCheckFit(fit = sf,data = TRUE,postpred = TRUE,statepred = TRUE,by = by,fastcov = TRUE,
        breaks=1,covplot = TRUE,smooth=TRUE,reg=TRUE,lagcovplot = TRUE,lag = lags,groupbysplit = TRUE))
      try(ctCheckFit(fit = sf,data = TRUE,postpred = TRUE,statepred = FALSE,by = by,breaks=4,covplot = TRUE,smooth=TRUE,reg=TRUE))
      try(ctCheckFit(fit = sf,data = TRUE,postpred = TRUE,statepred = FALSE,by = by,breaks=2,covplot = TRUE,lag=1,smooth=TRUE,reg=TRUE))
      for(mani in sm$manifestNames){
        try(ctCheckFit(fit = sf,data = TRUE,postpred = TRUE,statepred = FALSE,by = mani,
          breaks=2,covplot = TRUE,smooth=TRUE,reg=TRUE))
      }
      
      dev.off()
    }
    
    if(ctStanPlotPost){
      #paramter posterior plots
      if(sf$standata$nopriors==0){
        pdf(paste0(name,'_parameter_posterior.pdf'))
        ctStanPlotPost(sf,priorwidth = TRUE,cores=cores)
        dev.off()
      }
    }
    
    # #K-fold Cross Validation - OOS entropy score
    #   loo=ctLOO(fit = sf,parallelFolds = FALSE,folds = 10,cores = cores,subjectwise = TRUE,keepfirstobs = FALSE)
    #   save(loo,file='loo.rda')
    #   sink('loo.txt')
    #   print(loo)
    #   sink()  
    
  }
}

