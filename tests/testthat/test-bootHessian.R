if(FALSE){
  if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
    
    
    # library(future)
    # library(data.table)
    # plan(strategy = 'multisession',workers=10)
    
    testfunc <- function(i){
      message(i)
      nsubjects=200
      t0m <- 10+rnorm(nsubjects)
      cint <- t0m/2+rnorm(nsubjects)
      # cor(t0m,cint)
      for(subi in 1:nsubjects){
        gm <- suppressMessages(ctModel(Tpoints=10,
          LAMBDA=matrix(1), 
          DRIFT= -1,
          T0MEANS = t0m[subi],
          DIFFUSION=.5,
          MANIFESTVAR = 0.5,
          T0VAR = 0,
          MANIFESTMEANS = 0,
          CINT=cint[subi]))
        
        dd <- suppressMessages(data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 1,
          burnin = 0,dtmean = 1,logdtsd = 0)))
        dd$id <- subi
        if(subi==1) d <- dd else d <- rbind(d,dd)
      }
      
      m <- ctModel(type='stanct',LAMBDA=matrix(1),CINT='cint',MANIFESTMEANS=0)
      # m$pars$indvarying=F
      
      f <- ctStanFit(datalong = d,ctstanmodel = m,cores=1,nopriors=F,optimcontrol=list(carefulfit=F,stochastic=F,finishsamples=5000))
      s=summary(f)
      # print(s)
      
      scores=t(ctsem:::scorecalc(standata = f$standata,est = f$stanfit$rawest,stanmodel = f$stanmodel,subjectsonly = T,cores=1))
      
      fc <- ctStanFit(datalong = d,ctstanmodel = m,cores=1,nopriors=F,optimcontrol=list(carefulfit=F,stochastic=F,finishsamples=5000,bootstrapUncertainty=F))
      sc=summary(fc)
      # print(sc)
      
      # print(sc$popmeans)
      # print(s$popmeans)
      colnames(sc$popmeans) <- paste0(colnames(sc$popmeans),'c')
      colnames(sc$rawpopcorr) <- paste0(colnames(sc$rawpopcorr),'c')
      colnames(sc$popsd) <- paste0(colnames(sc$popsd),'c')
      dt=data.frame(rbind(s$popmeans,s$popsd,s$rawpopcorr[,1:5]),rbind(sc$popmeans,sc$popsd,sc$rawpopcorr[,1:5]))
      rownames(dt)=c(rownames(s$popmeans),paste0('popsd_',rownames(s$popsd)),rownames(s$rawpopcorr))
      return(dt)
    }
    
    out <- list()
    for(i in 1:100){
      out[[i]] <- future(testfunc(i))
    }
    out <- value(out)
    
    
    truepars <- out[[1]][,'mean',drop=FALSE]
    truepars[] <- c(10,-1,2,.5,5,1,1.116,.45)
    
    out <- lapply(1:length(out),function(o) data.frame(run=o,par=rownames(out[[o]]),out[[o]]))
    outdt <- rbindlist(out)
    outdtb <- data.frame(outdt)[,colnames(outdt)[!grepl('c$',colnames(outdt))]]
    outdtc <- data.frame(outdt)[,c('run','par',colnames(outdt)[grepl('c$',colnames(outdt))])]
    colnames(outdtc)=gsub('c$','',colnames(outdtc))
    outdt <- rbind(data.table(type='bootstrap',outdtb),data.table(type='classic',outdtc))
    outdt <- melt(outdt,id.vars = c('type','par','run'))
    require(ggplot2)
    ggplot(outdt[variable %in% c('X97.5.','X2.5.','X50.'),],aes(x=variable,y=value,col=type,group=interaction(run,type)))+
      geom_jitter(alpha=.5,width=.2,height=0)+geom_line(alpha=.2)+theme_bw()+
      geom_hline(aes(yintercept=value),data=data.frame(par=rownames(truepars),value=truepars$mean))+
      facet_wrap(vars(par),scales = 'free')
    
    
    covered <- lapply(out,function(x) data.table(
      b= truepars > x$`X2.5.` & truepars < x$`X97.5.`, 
      c= truepars > x$`X2.5.c` & truepars < x$`X97.5.c`))
    
    coverage <- as.matrix(covered[[1]])
    for(i in 1:nrow(coverage)){
      for(j in 1:ncol(coverage)){
        coverage[i,j]<-sum(sapply(covered,function(x) as.matrix(x)[i,j]))/length(out)
      }
    }
    rownames(coverage) <- rownames(truepars)
    # print(coverage)
    
    
    
  }
}
