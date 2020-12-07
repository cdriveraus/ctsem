if(1==99 && identical(Sys.getenv("NOT_CRAN"), "true") & .Machine$sizeof.pointer != 4){
  # Sys.setenv(NOT_CRAN='true')
  
  set.seed(1)
  library(ctsem)
  library(testthat)
  
  test_that("discretepars", {
    
    n <- 300
    dri <- rnorm(n,0,2)
    dr11=sort(-log1p(exp(dri+rnorm(n,0,.2))))
    dr12=(dri+1+rnorm(1,0,.3))/10
    
    for(i in 1:n){
      gm <- ctModel(LAMBDA=diag(2), DRIFT = c(dr11[i],dr12[i],0,-1),
        DIFFUSION = c(2,2,0,.1),
        MANIFESTVAR = diag(.2,2),Tpoints=30,T0MEANS = 10)
      d <- data.frame(ctGenerateOld(ctmodelobj = gm,n.subjects = 1,dtmean = .3),dri=dri[i])
      d$id <- i
      if(i==1) dat <- d else dat <- rbind(dat,d)
    }
    plot(dat$time,dat$Y1)
    plot(dat$Y2,dat$Y1)
    
    m <- ctModel(LAMBDA=diag(2),type='stanct',
      DRIFT = c('dr1||||dri','dr12||||dri','dr21','dr22'),
      TIpredNames = 'dri',tipredDefault = F)
    m$pars$indvarying<-FALSE
    
    f=ctStanFit(datalong = dat,ctstanmodel = m,plot=50,cores=6,
      savesubjectmatrices = T,verbose=0)
    
    
    ctStanDiscretePars(f,subjects=1:3,plot=T,nsamples = 50,standardise = TRUE)
    ctStanDiscretePars(f,subjects=1:3,plot=T,nsamples = 50,observational = T)
    ctStanDiscretePars(f,subjects=1:3,plot=T,nsamples = 50,observational = T,cov=T)
    ctStanDiscretePars(f,subjects=5:8,plot=T,nsamples = 50,facets='Subject')
    
    g=ctStanDiscretePars(f,subjects=1:3,plot=T,nsamples = 50,ggcode=T)
    cat(g$ggcode)
    
    ggplot2::ggplot(data = g$dt,mapping=aes(y=value,x=`Time interval`,
      colour=Effect,
      fill=Effect))+
      theme_bw()+ylab('Coef.')+
      ggplot2::labs(title = 'My plot')+  
      stat_summary( #ribbon
        fun.data = function(x) list(
          y=quantile(x,.5),
          ymin=quantile(x,.1), 
          ymax=quantile(x,.9)
        ),
        geom = "ribbon",
        alpha= .1,
        linetype=3)+
      stat_summary( #center line
        fun.data = function(x) list(
          y=quantile(x,.5)
        ),
        geom = "line",
        linetype=1)
    
  })
}
