# vgrep <- function(patterns,x){
#   unique(unlist(sapply(patterns,function(pattern) grep(pattern,x))))
# }

# mupd <- function(m, nr=NA,nc=NA, upd=NA, sr=NA,sc=NA){
#   if(is.na(nr)) nr <- max(nrow(upd)+sr-1,nrow(m))
#   if(is.na(nc)) nc <- max(ncol(upd)+sc-1,nrow(m))
#   mo <- rbind(
#     cbind(m, matrix(0,nrow(m),nc-ncol(m))),
#     matrix(0,nr-nrow(m),nc))
#   if(!is.na(upd)) mo[sr:(sr+nrow(upd)-1),sc:(ncol(upd)+sc-1)] <- upd
#   return(mo)
# }

# diagt <- function(d){
#   m <- diag(0,length(d))
#   m[diag(1,length(d))==1] <- d
#   return(m)
# }

# findmatrixslots <- function(pars,l){
#   p<-list()
#   for(pi in pars){
#     for(mi in 1:length(l)){
#       if(any(l[[mi]] %in% pi)){
#         arrind <- arrayInd(which(l[[mi]] %in% pi),dim(l[[mi]]))
#         p[[pi]] <- paste0(names(l)[mi],'[',arrind[1,1],',',arrind[1,2],']')
#         next
#       }
#     }
#   }
#   return(p)
# }



ctstantestfitfunc<-function(){
  checkm<-ctModel(
    type='ct',
    n.latent=2,n.TDpred=1,n.TIpred=1,n.manifest=2,
    MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
    MANIFESTMEANS=0,
    DIFFUSION=c('diff11',0,'diff21','diff22||||TI1'),
    CINT=matrix(c('cint1||||TI1','cint2||||TI1'),ncol=1),
    LAMBDA=diag(2),tipredDefault=FALSE)  
  
  ctstantestfit<-ctStanFit(ctstantestdat,checkm,cores=1,
    inits = c(0.748310681869536,0.945659953796114,0.0964592332562144,
      0.029153487981562,0.651471066485501,0.0314778013950629,
      0.217818608752396,1.10441297459423,-0.801320300354595,
      0.647010811111734,-0.7344068376597,-1.04150782976995,
      0.0558819480347101,-0.108435212373754,-0.225029736388403,
      -0.203457959897841,-0.736264486213394,-0.687369939087293,
      0.576641002392084,0.248625561427667,-0.0683189297539777,
      -0.230342395895042,0.205299380670756,-0.34522281922735,
      0.0829819407118698,0.0137367678089216,-0.0611697475527028),
    optimize = TRUE,optimcontrol=list(finishsamples=20),priors=TRUE)
  
  ctstantestfit <- ctStanGenerateFromFit(ctstantestfit,nsamples = 20,fullposterior = TRUE)
  
  return(ctstantestfit)
}




# removeOutliers <- function(dat,multiplier,by=2){
#   dat2 <- array(apply(dat,by,function(x){
#     s=sd(x,na.rm=TRUE)
#     m=mean(x,na.rm=TRUE)
#     message("Removed ", sum(abs(x-m) > (multiplier*s),na.rm=TRUE)," outliers...")
#     x[abs(x-m) > (multiplier*s)] <- NA
#     return(x)
#   }),dim=dim(dat))
# }


testall<- function(cores=4,folder = '/tests/testthat',examples=TRUE){
  requireNamespace('testthat')
  Sys.setenv(NOT_CRAN='true')
  pdf(NULL)
  tests <- dir(paste0('.',folder))
  tests <- tests[grepl('^test',tests)]
  runex <- grep('runExamples',tests)
  tests <- c(tests[runex],tests[-runex]) #do examples first
  if(!examples) tests <- tests[-grep('runExamples',tests)]
  a=Sys.time()

  if(cores > 1){
    cl <- parallelly::makeClusterPSOCK(cores)
    on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
    out <- parallel::parLapplyLB(cl,paste0(getwd(),folder,'/',tests),function(x){
      Sys.setenv(NOT_CRAN='true')
      pdf(NULL)
    out<-testthat::test_file(x, reporter = "minimal")
    dev.off()
    return(out)
  })
  }
  if(cores==1){
    out <- lapply(paste0(getwd(),folder,'/',tests),function(x){
      cat(x)
      out<-testthat::test_file(x, reporter = "minimal")
      print(out)
      return(out)
    })
  }
  out2 <- do.call(what = rbind,lapply(out,utils::getS3method('as.data.frame','testthat_results')))
  dev.off()
  print(out2[,colnames(out2)!='result'])
  print(Sys.time()-a)
  if(cores > 1) parallel::stopCluster(cl)
  return(invisible(out2))
}
  


suppressOutput <- function(...,verbose=0){
  if(verbose > 0) return(eval(...)) else return(capture.output(eval(...)))
}

  openPDF <- function(f) {
    os <- .Platform$OS.type
    if (os=="windows")
      shell.exec(normalizePath(f))
    else {
      pdf <- getOption("pdfviewer", default='')
      if (nchar(pdf)==0)
        stop("The 'pdfviewer' option is not set. Use options(pdfviewer=...)")
      system2(pdf, args=c(f))
    }
  }

naf <-function(x){
  x[is.na(x)] <- FALSE
  return(x)
}



meltkalman <- function(l){
if(1==99) Row <- Col <- NULL
  Time <- data.table(l$time,Obs=1:length(l$time))
  Subject <- data.table(Subject=factor(l$id),Obs=seq_along(l$time))
  l$id <- NULL
  l$time <- NULL
  TimeSubject <- merge(Time,Subject)
  
  dout <- NULL
  i=0
  while(i < length(l)){
    i=i+1
    
    while(length(dim(l[[i]])) < 3){
      dn = c(dimnames(l[[i]]),list(NULL))
      l[[i]] <- array(l[[i]],dim = c(dim(l[[i]]),1))
      dimnames(l[[i]]) <- dn
      
    }
      x <- melt(as.data.table(l[[i]],keep.rownames = TRUE,na.rm = FALSE),measure.vars = 'value')
      x$Obs <- as.integer(x$Obs)
      x <- merge(cbind(x,Element = names(l)[i]),TimeSubject,by='Obs')
      if(grepl('(prior)|(upd)|(smooth)',names(l)[i]) &!grepl('(cov$)|(^err)',names(l)[i])){
        # browser()
        xsd<-melt(as.data.table(l[[ paste0(names(l)[i],'cov') ]],
          keep.rownames = TRUE,na.rm = FALSE),measure.vars='value')
        xsd <-subset(xsd,Row==Col)
        xsd$value <- sqrt(xsd$value+1e-8)
        xsd$Obs <- as.integer(xsd$Obs)
        setnames(xsd,'value','sd')
        x<-merge(x,xsd,by=colnames(xsd)[colnames(xsd) %in% colnames(x)]) #data.table(sd=(sqrt(xsd$value))))
      }
      if(is.null(dout)) dout <- x else dout <- rbind(dout,x,fill=TRUE)
  }
  
  dout <- data.frame(dout) #because of weird temporal data.table behaviour
  class(dout) <- c('ctKalmanDF',class(dout))
return(dout)
}

# gridplot <- function(m, maxdim=c(3,3),...){
#   d=n2mfrow(dim(m)[length(dim(m))])
#   d[d>maxdim] <-maxdim[d>maxdim]
#   oldpar<-par(no.readonly=TRUE)
#   par(mfrow=d,mar=c(1.1,1.1,1.1,0),mgp=c(.1,.1,0))
#   for(i in 1:dim(m)[length(dim(m))]){
#     n=colnames(m)[i]
#     if('matrix' %in% class(m)) plot(m[,i],main=ifelse(is.null(n),i,n),col='red',xlab='',ylab='',...)
#     if('array' %in% class(m)) matplot(m[,,i],main=ifelse(is.null(n),i,n),type='l',xlab='',ylab='',...)
#   }
#   suppressWarnings(do.call(par,oldpar))
# }

# perm <- function(v) {
#   n <- length(v)
#   if (n == 1) v
#   else {
#     X <- NULL
#     for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
#     X
#   }
# }

# helper function to generate an index matrix, or return unique elements of a matrix
indexMatrix<-function(dimension,symmetrical=FALSE,upper=FALSE,lowerTriangular=FALSE, sep=NULL,starttext=NULL,endtext=NULL,
  unique=FALSE,rowoffset=0,coloffset=0,indices=FALSE,diagonal=TRUE,namesvector=NULL,shortdiag=FALSE){
  if(is.null(namesvector)) namesvector=1:9999
  if(indices==TRUE) sep<-c(",")
  tempmatrix<-matrix(paste0(starttext,rep(namesvector[1:dimension+coloffset],dimension)),nrow=dimension,ncol=dimension)
  for(i in 1:nrow(tempmatrix)){ #append step by step for shortdiag
    for(j in 1:ncol(tempmatrix)){
        # tempmatrix[i,j] <- paste0(tempmatrix[i,j],sep,namesvector[1:dimension+rowoffset][i])
        # print(tempmatrix[i,j])
        # print(namesvector[1:dimension+rowoffset][i])
        if(i!=j || !shortdiag) tempmatrix[i,j] <- paste0(tempmatrix[i,j],sep,namesvector[1:dimension+rowoffset][j])
    }
  }
  tempmatrix[,]<-paste0(tempmatrix[,],endtext)
      
  if(upper==TRUE) tempmatrix<-t(tempmatrix)
  if(symmetrical==TRUE) tempmatrix[col(tempmatrix)>row(tempmatrix)] <-t(tempmatrix)[col(tempmatrix)>row(tempmatrix)]
  if(unique==TRUE && symmetrical==TRUE) tempmatrix<-tempmatrix[lower.tri(tempmatrix,diag=diagonal)]
  if(lowerTriangular==TRUE) tempmatrix[col(tempmatrix) > row(tempmatrix)] <- 0
  if(indices==TRUE){
    tempmatrix<-matrix(c(unlist(strsplit(tempmatrix,","))[seq(1,length(tempmatrix)*2,2)],
      unlist(strsplit(tempmatrix,","))[seq(2,length(tempmatrix)*2,2)]),ncol=2)
  }
  return(tempmatrix)
}


relistarrays <- function(flesh, skeleton){
  skelnames <- names(skeleton)
  skelstruc <- lapply(skeleton,dim)
  count=1
  npars <- length(flesh)
  out <- list()
  for(ni in skelnames){
    if(!is.null(skelstruc[[ni]])){
      out[[ni]] <- array(flesh[count:(count+prod(skelstruc[[ni]]))],dim = skelstruc[[ni]])
      count <- count + prod(skelstruc[[ni]])
    } else {
      out[[ni]] <- flesh[count]
      count <- count + 1
    }
  }
  return(out)
}

#convert id's to ascending numeric 
makeNumericIDs <- function(datalong,idName='id',timeName='time'){
  originalid <- unique(datalong[,idName])
  datalong[,idName] <- match(datalong[,idName],originalid)
  
  datalong <- datalong[order(datalong[,idName],datalong[,timeName]),] #sort by id then time
  
  if(any(is.na(as.numeric(datalong[,idName])))) stop('id column may not contain NA\'s or character strings!')
  return(datalong)
}


crosscov <- function(a,b){
  da <- a-matrix(colMeans(a),nrow=nrow(a),ncol=ncol(a),byrow=TRUE)
  db <- b-matrix(colMeans(b),nrow=nrow(b),ncol=ncol(b),byrow=TRUE)
  t(da) %*% db / (nrow(a)-1)
  
  # cc <- matrix(NA,nrow=nrow(a))
  # for(i in 1:nrow(a)){
  #   cc[i,] <- da[i,] * db[i,]
  
}



#' ctCollapse
#' Easily collapse an array margin using a specified function.
#' @param inarray Input array of more than one dimension.
#' @param collapsemargin Integers denoting which margins to collapse.
#' @param collapsefunc function to use over the collapsing margin.
#' @param plyr Whether to use plyr.
#' @param ... additional parameters to pass to collapsefunc.
#' @examples
#' testarray <- array(rnorm(900,2,1),dim=c(100,3,3))
#' ctCollapse(testarray,1,mean)
#' @export
ctCollapse<-function(inarray,collapsemargin,collapsefunc,plyr=TRUE,...){
  indims<-dim(inarray)
  if(plyr) out<-array(plyr::aaply(inarray,(1:length(indims))[-collapsemargin],collapsefunc,...,
    .drop=TRUE),dim=indims[-collapsemargin])
  if(!plyr) out<-array(apply(inarray,(1:length(indims))[-collapsemargin],collapsefunc,...,
    .drop=TRUE),dim=indims[-collapsemargin])
  dimnames(out)=dimnames(inarray)[-collapsemargin]
  return(out)
}

rl<-function(x) { #robust logical - wrap checks likely to return NA's in this
  x[is.na(x)] <- FALSE
  return(x)
}




#' Inverse logit
#' 
#' Maps the stan function so the same code works in R.
#'
#' @param x value to calculate the inverse logit for. 
#'
#' @examples
#' inv_logit(-3)
#' @export
inv_logit<-function(x) {
  1/(1+exp(-x))
}

#' log1p_exp
#' 
#' Maps the stan function so the same code works in R.
#'
#' @param x value to use. 
#'
#' @examples
#' log1p_exp(-3)
#' @export
log1p_exp <- function(x) log1p(exp(x))

#' ctDensity
#'
#' Wrapper for base R density function that removes outliers and computes 'reasonable' bandwidth and x and y limits.
#' Used for ctsem density plots.
#' 
#' @param x numeric vector on which to compute density.
#' @param bw either 'auto' or a numeric indicating bandwidth.
#' @param plot logical to indicate whether or not to plot the output.
#' @param ... Further args to density.
#' 
#' @examples
#' y <- ctDensity(exp(rnorm(80)))
#' plot(y$density,xlim=y$xlim,ylim=y$ylim)
#' 
#' #### Compare to base defaults:
#' par(mfrow=c(1,2))
#' y=exp(rnorm(10000))
#' ctdens<-ctDensity(y)
#' plot(ctdens$density, ylim=ctdens$ylim,xlim=ctdens$xlim)
#' plot(density(y))
#' @export

ctDensity<-function(x,bw='auto',plot=FALSE,...){
  xlims=stats::quantile(x,probs=c(.05,.95),na.rm=TRUE)
  sd=sd(xlims)
  xlims[1] = xlims[1] - sd
  xlims[2] = xlims[2] + sd
  # x=x[x>xlims[1]*1.2 & x<xlims[2]*1.2]
  # bw=(max(x)-min(x))^1.2 / length(x)^.4 *.4
  if(bw=='auto') bw=min(sd/100,1e-4)
  
  # xlims=stats::quantile(x,probs=c(.01,.99))
  # mid=mean(c(xlims[2],xlims[1]))
  # xlims[1] = xlims[1] - (mid-xlims[1])/8
  # xlims[2] = xlims[2] + (xlims[2]-mid)/8
  
  out1<-stats::density(x,bw=bw,n=5000,from=xlims[1]-sd,to=xlims[2]+sd)
  ylims=c(0,max(out1$y)*1.1)
  
  if(plot) plot(out1$x, out1$y,type='l', xlim=xlims,ylim=ylims,ylab='Density',xlab='Par. Value',...)
  
  return(list(density=out1,xlim=xlims,ylim=ylims))
}



ctDensityList<-function(x,xlimsindex='all',ylimsindex='all',cut=FALSE,plot=FALSE,
  grouplabels=names(x),
  ylab='Density',
  xlab='Par. Value',probs=c(.05,.95),main='',colvec=NA){
  
  if(all(xlimsindex=='all')) xlimsindex <- 1:length(x)

  for(i in xlimsindex){
    newxlims=stats::quantile(x[[i]],probs=probs,na.rm=TRUE)
    if(i==1) {
      xlims=newxlims
    }
    else {
      xlims <- range(c(xlims,newxlims))
    }
  }
  sd=sd(xlims)
  xlims[1] = xlims[1] - sd/2
  xlims[2] = xlims[2] + sd/2
# browser()
  bw=sapply(x,function(d) {
    out<-try(bw.SJ(na.omit(c(d))),silent=TRUE)
    if('try-error' %in% class(out)) out <- bw.nrd0(na.omit(c(d)))
    return(out)
    })
  # browser()
  # logbwmean = mean(logbw)# + ifelse(length(logbw) > 1,sd(logbw),0))
# browser()
  denslist<-lapply(1:length(x),function(xi) {
    bw=exp(log(bw[xi]) + (mean(log(bw))-log(bw[xi]))*.3)
    # print(bw)
    d=stats::density(x[[xi]],bw=bw,n=5000,from=xlims[1],to=xlims[2],na.rm=TRUE)
    # d$y=d$y/ sum(d$y)/range(d$x)[2]*length(d$y)
    return(d)
  })
if(xlimsindex[1] %in% 'all') xlimsindex <- 1:length(denslist)
if(ylimsindex[1] %in% 'all') ylimsindex <- 1:length(denslist)
  xlims=range(sapply(denslist[xlimsindex],function(d) d$x[d$y> (.1*max(d$y))]))
  xlims <- xlims +c(-1,1)*sd(xlims)
  ylims=c(0,max(unlist(lapply(denslist[ylimsindex],function(li) max(li$y))))*1.1)
  
  if(cut){
    denslist<-lapply(denslist,function(d){
      d$x[d$x < min(xlims)] <- NA
      d$x[d$x > max(xlims)] <- NA
      d$y<-d$y[!is.na(d$x)] 
      d$x<-d$x[!is.na(d$x)] 
      return(d)
    })
  }

  denslist <- lapply(denslist,function(d){
    o<-list()
    o$y <- d$y[between(d$x,xlims[1],xlims[2]) & d$y > .001*max(d$y)]
    o$x <- d$x[between(d$x,xlims[1],xlims[2])& d$y > .001*max(d$y)]
    return(o)
  })
  
  if(plot) {
   if(is.null(names(x))) names(x) <- paste0('var ', 1:length(x))
    
    pd <- data.table(Source=names(x)[1],denslist[[1]]$x,denslist[[1]]$y)
    if(length(denslist)>1){
      for(di in 2:length(denslist)){
      pd <- rbind(pd,data.table(Source=names(x)[di],denslist[[di]]$x,denslist[[di]]$y))
      }
    }
    names(pd)[3] <- 'y'
    names(pd)[2] <- 'x'
    # pd$Source <- factor(pd$Source)
    
    if(1==99) x <- y <- Source <- NULL
    # plot(denslist[[1]]$x, denslist[[1]]$y,type='l', xlim=xlims,ylim=ylims,ylab=ylab,xlab=xlab,col=colvec[1],lty=ltyvec[1],...)

    plots<- ggplot(pd,aes(x=x,y=y,group=Source,colour=Source)) +
      geom_line() +
      coord_cartesian(xlim = xlims, ylim = ylims) +
      theme_minimal()+
      labs(y=ylab,x=xlab)+
      theme(legend.title = element_blank())
    if(!all(is.na(colvec))) plots <- plots + scale_colour_manual(values=setNames(colvec, names(x)))

  #     
  #   if(length(denslist)>1){
  #     for(ci in 2:length(denslist)){
  #       points(denslist[[ci]]$x, denslist[[ci]]$y,type='l', col=colvec[ci],lty=ltyvec[ci],...)
  #     }
  #   }
  #   if(all(legend!=FALSE)) {
  #     if(is.null(legendargs$col)) legendargs$col = colvec
  #     if(is.null(legendargs$text.col)) legendargs$text.col = colvec
  #     if(is.null(legendargs$lty)) legendargs$lty = ltyvec
  #     if(is.null(legendargs$x)) legendargs$x='topright'
  #     if(is.null(legendargs$bty)) legendargs$bty='n'
  #     legendargs$legend = legend
  #     do.call(graphics::legend,legendargs)
  #   }
    
  }
  
  if(plot) return(plots) else return(list(density=denslist,xlim=xlims,ylim=ylims))
}


#' Plots uncertainty bands with shading
#'
#' @param x x values
#' @param y y values
#' @param ylow lower limits of y
#' @param yhigh upper limits of y
#' @param steps number of polygons to overlay - higher integers lead to 
#' smoother changes in transparency between y and yhigh / ylow.
#' @param ... arguments to pass to polygon()
#'
#' @return Nothing. Adds a polygon to existing plot.
#' @export
#'
#' @examples
#' plot(0:100,sqrt(0:100),type='l')
#' ctPoly(x=0:100, y=sqrt(0:100), 
#' yhigh=sqrt(0:100) - runif(101), 
#' ylow=sqrt(0:100) + runif(101),
#' col=adjustcolor('red',alpha.f=.1))
ctPoly <- function(x,y,ylow,yhigh,steps=20,...){
  for(i in 1:steps){
    tylow= y + (ylow-y)*i/steps
    tyhigh= y + (yhigh-y)*i/steps
    xf <- c(x,x[length(x):1])
    yf <- c(tylow,tyhigh[length(tyhigh):1])
    polygon(xf,yf,border=NA,...)
  }
}



#' ctWideNames
#' sets default column names for wide ctsem datasets. Primarily intended for internal ctsem usage.
#' @param n.manifest number of manifest variables per time point in the data.
#' @param Tpoints Maximum number of discrete time points (waves of data, or measurement occasions) 
#' for an individual in the input data structure.
#' @param n.TDpred number of time dependent predictors in the data structure.
#' @param n.TIpred number of time independent predictors in the data structure.
#' @param manifestNames vector of character strings giving column names of manifest indicator variables
#' @param TDpredNames vector of character strings giving column names of time dependent predictor variables
#' @param TIpredNames vector of character strings giving column names of time independent predictor variables
#' @export

ctWideNames<-function(n.manifest,Tpoints,n.TDpred=0,n.TIpred=0,manifestNames='auto',TDpredNames='auto',TIpredNames='auto'){
  
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  
  if(length(manifestNames) != n.manifest) stop("Length of manifestNames does not equal n.manifest!") 
  
  if(n.TDpred > 0){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
    if(length(TDpredNames) != n.TDpred) stop("Length of TDpredNames does not equal n.TDpred!") 
  }
  
  if(n.TIpred > 0){
    if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
    if(length(TIpredNames) != n.TIpred) stop("Length of TIpredNames does not equal n.TIpred!") 
  }
  
  manifestnames<-paste0(manifestNames,"_T",rep(0:(Tpoints-1),each=n.manifest))
  if(n.TDpred > 0 && Tpoints > 1) {
    TDprednames<-paste0(TDpredNames,"_T",rep(0:(Tpoints-1),each=n.TDpred))
  } else {
    TDprednames<-NULL
  }
  if (Tpoints > 1) {
    intervalnames<-paste0("dT",1:(Tpoints-1))
  } else {
    intervalnames <- NULL
  }
  if(n.TIpred>0) TIprednames <- paste0(TIpredNames) else TIprednames <- NULL
  return(c(manifestnames,TDprednames,intervalnames,TIprednames))
}

# generates more complex sequences than seq
cseq <- function(from, to, by){
  temp<-c()
  for(i in from){
    temp<-c(temp,seq(i,to,by))
  }
  temp<-sort(temp)
  return(temp)
}

# get_stan_params <- function(object) {
#   stopifnot(methods::is(object, "stanfit"))
#   params <- grep("context__.vals_r", fixed = TRUE, value = TRUE,
#     x = strsplit(rstan::get_cppcode(rstan::get_stanmodel(object)), "\n")[[1]])
#   params <- sapply(strsplit(params, "\""), FUN = function(x) x[[2]])
#   params <- intersect(params, object@model_pars)
#   return(params)
# }


# get_stan_massmat<-function(fit){
#   
#   spars<-get_stan_params(fit)
#   spars2<-c()
#   for(pari in spars){
#     spars2<-c(spars2,grep(paste0(pari,'['),names(fit@sim$samples[[1]]),fixed=TRUE))
#   }
#   
#   massmat<-list()
#   for(chaini in 1:fit@sim$chains){
#     temp<-c()
#     for(pari in spars2){
#       newval<-stats::cov(cbind(fit@sim$samples[[chaini]][[pari]][(fit@sim$warmup - fit@stan_args[[1]]$control$adapt_term_buffer):fit@sim$warmup]))
#       names(newval)<-names(fit@sim$samples[[chaini]])[pari]
#       temp<-c(temp,newval)
#     }
#     massmat[[chaini]]<-temp
#   }
#   return(massmat)
# }


