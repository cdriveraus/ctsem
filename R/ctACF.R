

#' Continuous Time Autocorrelation Function (ctACF)
#'
#' This function computes an approximate continuous time autocorrelation function (ACF) for data
#' containing multiple subjects and/or variables.
#'
#' @param dat The input data in data frame or data table format.
#' @param varnames Character vector of variable names in the data to compute the ACF for. 'auto' uses all columns that are not time / id. 
#' @param ccfnames Character vector of variable names in the data to compute cross correlation for. 'all' uses all variables in varnames, NA uses none.  
#' @param idcol The name of the column containing subject IDs (default is 'id').
#' @param timecol The name of the column containing time values (default is 'time').
#' @param plot A logical value indicating whether to create a plot (default is TRUE).
#' @param timestep The time step for discretizing data. 'auto' to automatically determine
#'        the timestep based on data distribution (default is 'auto'). 
#'        In this case the timestep is computed as half of the median for time intervals in the data.
#' @param time.max The maximum time lag to compute the ACF (default is 10). If 'auto', is set to 10 times the 90th percentile interval in the data.
#' @param nboot The number of bootstrap samples for confidence interva1l estimation (default is 100).
#' @param scale if TRUE, scale variables based on within-subject standard deviation.
#' @param center if TRUE, center variables based on within-subject mean.
#' @param ... additional arguments (such as demean=FALSE) to pass to the \code{stats::acf} function.
#'
#' @return If 'plot' is TRUE, the function returns a ggplot object of the ACF plot. If 'plot' is
#'         FALSE, it returns a data table with ACF estimates and confidence intervals.
#'
#' @details This function computes the continuous time ACF by discretizing the data and then
#' performing bootstrapped ACF calculations to estimate the confidence intervals. It can create
#' ACF plots with confidence intervals if 'plot' is set to TRUE.
#'
#' @seealso \code{\link{ctDiscretiseData}}
#'
#' @examples
#' data.table::setDTthreads(1) #ignore this line
#' # Example usage:
#' head(ctstantestdat)
#' ctACF(ctstantestdat,varnames=c('Y1'),idcol='id',timecol='time',nboot=5)
#'
#' @export
ctACF <- function(dat, varnames='auto',ccfnames='all',idcol='id', timecol='time',
  plot=TRUE,timestep='auto',time.max='auto',nboot=100,scale=FALSE, center=FALSE,...){
  
  if(requireNamespace('collapse')){
  } else stop('collapse package needed for ACF!')
  
  if(F) .timediff = ACFhigh= ACFlow= Element= Estimate= SignificanceLevel=  TimeInterval= Variable= boot= ci=NULL
  dat=copy(data.table(dat))
  
  if(length(varnames)==1 && varnames[1] == 'auto') varnames <- colnames(dat)[!colnames(dat) %in% c(idcol,timecol,'.timediff')]
  if(length(ccfnames )==1 && ccfnames %in% 'all') ccfnames=varnames
  if(length(ccfnames)==1 && is.na(ccfnames)) ccfnames=NULL
  acfnames <- varnames
  varnames <- unique(c(acfnames,ccfnames))
  
  if(timestep == 'auto'){
    timestep = 0.5 * quantile(dat[,.timediff:= c(NA,diff(get(timecol))),by=idcol][['.timediff']],
      probs=.5,na.rm=TRUE)
    message('Timestep used: ',timestep)
  }
  if(time.max == 'auto')time.max = 2/3*mean(dat[,.timerange:= diff(range(get(timecol))),by=idcol][['.timerange']])
  
  
  # Function to remove rows with missing values before and after non-missing values
  # remove_rows_with_missing <- function(dt, col_name) {
  #   dt[, {
  #     if (all(is.na(.SD[[col_name]]))) {
  #       .SD[.N]  # Return a single row with NA values if all values are missing
  #     } else {
  #       first_non_missing <- min(which(!is.na(.SD[[col_name]])))
  #       last_non_missing <- max(which(!is.na(.SD[[col_name]])))
  #       .SD[first_non_missing:last_non_missing, ]
  #     }
  #   }, by = idcol]
  # }
  
  # remove_rows_with_missing <- function(dt, col_names) {
  #   dt[, {
  #     row_has_data <- rowSums(!is.na(.SD[, col_names,with=FALSE])) > 0
  #     first_non_missing <- min(which(row_has_data))
  #     last_non_missing <- max(which(row_has_data))
  #     if (is.na(first_non_missing) || is.na(last_non_missing)) {
  #       .SD[1:.N]
  #     } else {
  #       .SD[first_non_missing:last_non_missing, ]
  #     }
  #   }, by = idcol]
  # }
  
  if(scale || center){
    dat[[idcol]] <- factor(dat[[idcol]])
    for(vi in varnames){ #scale each subject (skip demeaning so we can see trait variance)
      model <- lme4::lmer(formula=paste0(vi,' ~ 0+(1 | ',idcol,')'), data = dat)
      randomint=as.data.table(lme4::ranef(model))
      setnames(randomint,old = 'grp',new=idcol)
      dat <- merge(dat,randomint[,c('id','condval'),with=FALSE],by = idcol)
      dat[, (vi) := (get(vi) - condval) / sd(get(vi),na.rm=TRUE)+condval/sd(get(vi),na.rm=TRUE),by=idcol]
      dat[['condval']] <- NULL
    }
    dat[[idcol]] <- as.character(dat[[idcol]])
  }
  
  counter <- 0
  for(vari in acfnames){
    for(varj in unique(c(vari,ccfnames))){
      counter <- counter + 1
      varji <- unique(c(vari,varj))
      
      vdat <- dat[,c(idcol,timecol,varji),with=FALSE]
      vdat <- vdat[apply(vdat[,(varji),with=F],1,function(x) !all(is.na(x))),] #remove rows with total missingness first
      for(col in varji) set(vdat, j = col, value = as.numeric(vdat[[col]])) #ensure cols are numeric!
      # vdat[apply(.SD,1,function(x) !(all(is.na(x)))), ,.SDcols=varji]
      # vdat<-remove_rows_with_missing(vdat,varji)
      # for(varjik in varji) vdat[,(varjik):=scale(get(varjik)),by=idcol]
      vdat = data.table(ctDiscretiseData(vdat,timecol = timecol,idcol = idcol,timestep=timestep)) #collapse to specified time step
      
      # datt=vdat[,tail(.SD,1),by=idcol] #adding distant time point (outside acf range) to ensure sufficient padding with NA's for stacking (avoiding attenuation of low lag corr)
      # datt=datt[,list(.newtime=seq(get(timecol),get(timecol)+time.max+timestep,timestep)),by=idcol]
      # setnames(datt,old = '.newtime',timecol)
      # vdat=merge(vdat,datt,by = c(idcol,timecol),all=TRUE)
      # 
      # acf=vdat[,{
      #   acftemp=acf(get(varji),na.action = na.pass,plot=FALSE)$acf[,1,1];
      #   data.table(TimeInterval=seq(0,length.out=length(acftemp),by=timestep),ACF=acftemp)
      # }, by=idcol]
      
      # psacf(vdat[[vari]][rows], g=vdat[[idcol]][rows],lag.max=ceiling(time.max/timestep), 
      # plot=F, na.action=na.pass,...)
      # browser()
      ACF=lapply(0:nboot,function(x) {
        message(sprintf(paste0("\r ",paste(varji,collapse=' ')," %3d%%"), as.integer(x/nboot*100)),appendLF = FALSE)
        
        if(x==0){
          subjects <- vdat[[idcol]]
          rows <- 1:nrow(vdat)
        } else{
          subjects <- sample(unique(vdat[[idcol]]),replace=TRUE)
          rows <- unlist(lapply(subjects,function(si) which(vdat[[idcol]] %in% si)))
          subjects <- unlist(lapply(subjects,function(si){
            r=which(vdat[[idcol]] %in% si)
            paste0('s',round(rnorm(1),5),vdat[[idcol]][r])
          }))
        }
        if(vari==varj) acfout <- collapse:::psacf(vdat[[vari]][rows], 
          g=subjects,gscale=FALSE,lag.max=ceiling(time.max/timestep), 
          plot=F,...)$acf[,1,1] #
        if(vari!=varj) acfout <- collapse:::psccf(x = vdat[[vari]][rows],y = vdat[[varj]][rows],
          g=subjects,gscale=FALSE,
          lag.max =ceiling(time.max/timestep),plot=FALSE,...)$acf[,1,1]
        
        out <- data.table(ACF=acfout)
        out[,TimeInterval:=seq(0,length.out=.N,by=timestep)]
        if(vari!=varj) out[,TimeInterval:=TimeInterval-median(TimeInterval)]
        return(out)
      })
      message('')
      ACF=rbindlist(ACF,idcol = 'Sample')
      ACF[,Sample:=Sample-1]
      ACF=ACF[TimeInterval!=0,]
      ACF[,Variable:=ifelse(vari==varj,vari,paste0(vari,'_',varj))]
      ACF[,SignificanceLevel:= ifelse(vari==varj,
        qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(dat[[vari]]))),
        NA)]
      # if(nboot >0)  ACF <- ctACFquantiles(ACF,c(.025,.975)) else ACF[,Sample:=NULL]
      
      if(counter>1) fullACF <- rbind(fullACF,ACF) else fullACF <- ACF
    }
  }
  
  
  if(plot){
    return(plotctACF(fullACF))
  } else return(fullACF)
}

ctACFquantiles<-function(ctacfobj,quantiles=c(.025,.5,.975),separateLearnRates=FALSE,df='auto'){
  
  if(df =='auto'){
    ndt=length(unique(ctacfobj$TimeInterval))
    dfmax=ndt-1
    df=min(c(dfmax,
      ceiling(ndt/3)))
  }
  ctacfobj <- data.table(data.frame(ctacfobj))
  
  if(!requireNamespace('qgam')) stop("qgam package required for ACF plots: install.packages('qgam')")
  
  learnrate <- NA #estimate the learning rates based on one variable combination, but if that fails, try the next until one succeeds. 
  for(vari in unique(ctacfobj$Variable)){
    if(is.na(learnrate[1])){
      try({
        learnrate=qgam::tuneLearnFast(data = ctacfobj[Variable %in% vari,],qu = quantiles,
          form = ACF ~ s(TimeInterval,bs='cc',k=df),
          control=list(tol=.Machine$double.eps^0.25,progress=FALSE),
          argGam=list(select=TRUE))
      })
    }
  }
  # browser()
  for(vari in unique(ctacfobj$Variable)){
    if(separateLearnRates) try({
      learnrate=qgam::tuneLearnFast(data = ctacfobj[Variable %in% vari,],qu = quantiles,
        form = ACF ~ s(TimeInterval,bs='tp',k=df),
        control=list(tol=.Machine$double.eps^0.1,progress=FALSE),
        argGam=list(select=TRUE,optimizer=c('outer','bfgs')))
    })
    
    qg=qgam::mqgam(data = ctacfobj[Variable %in% vari,],qu = quantiles,
      form = ACF ~ s(TimeInterval,bs='tp',k=df),
      # control=list(tol=.Machine$double.eps^0.1,progress=T),
      err=learnrate$err,
      lsig = learnrate$lsig,
      argGam=list(select=TRUE,gamma=2))
    pred=qgam::qdo(obj = qg,qu = quantiles, fun = predict,newdata=ctacfobj[Variable %in% vari,])
    ctacfobj[Variable %in% vari,paste0('Q',quantiles*100,'%'):=pred]
    # })
  }
  return(ctacfobj)
}

#' Plot an approximate continuous-time ACF object from ctACF
#'
#' @param ctacfobj object
#' @param df df for the basis spline.
#' @param separateLearnRates if TRUE, estimate the learning rate for the quantile splines for each combination of variables. Slower but theoretically more accurate. 
#' @param reducedXlim if non-zero, n timesteps are removed from the upper and lower end of the x range 
#' where the spline estimates are less likely to be reasonable. 
#' @param estimateSpline if TRUE, quantile spline regression is used, otherwise the samples are simply plotted as lines and the other arguments here are not used.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' data.table::setDTthreads(1) #ignore this line
#' # Example usage:
#' head(ctstantestdat)
#' ac=ctACF(ctstantestdat,varnames=c('Y1'),idcol='id',timecol='time',timestep=.5,nboot=5,plot=FALSE)
#' plotctACF(ac, reducedXlim=0)
plotctACF <- function(ctacfobj,df='auto',quantiles=c(.025,.5,.975),
  separateLearnRates=FALSE, reducedXlim=1,estimateSpline=TRUE){
  
  
  # browser()
  # ctacfobj=melt(ctacfobj,measure.vars = paste0('Q',quantiles*100,'%'),variable.name = 'Quantile',value.name = 'Corr')
  if(estimateSpline){
    ctacfobj=copy(ctACFquantiles(ctacfobj,quantiles=quantiles,
      separateLearnRates=separateLearnRates,df=df))
    ctacfobj <- ctacfobj[abs(TimeInterval) < min(tail(sort(unique(TimeInterval)),reducedXlim)),] #drop the ends of the spline
    gg <- ggplot(ctacfobj,aes(y=!!sym(paste0('Q',quantiles[2]*100,'%')),
      x=TimeInterval))+
      geom_ribbon(aes(ymin=!!sym(paste0('Q',quantiles[1]*100,'%')),
        ymax=!!sym(paste0('Q',quantiles[3]*100,'%'))),alpha=.2,linetype='dashed',linewidth=.5)+
      geom_line(linewidth=1)+
      # guides(linewidth='none',linetype='none')+
      geom_hline(yintercept=0,linetype='dotted')+
      geom_vline(xintercept=0,linetype='dotted')+
      ylab('ACF')+
      theme_bw()
  } else {
    gg <- ggplot(ctacfobj,aes(y=ACF,x=TimeInterval,group=Sample))+
      geom_line(alpha=.2)+
      geom_hline(yintercept=0,linetype='dotted')+
      geom_vline(xintercept=0,linetype='dotted')+
      ylab('ACF')+
      theme_bw()
  }
  
  if(length(unique(ctacfobj[['Variable']])) > 1) gg <- gg+facet_wrap(vars(Variable),scales = 'free')
  return(gg)
}

#' Calculate Continuous Time Autocorrelation Function (ACF) for Standardized Residuals of ctsem fit.
#'
#' This function takes a fit object from ctsem and computes the continuous time autocorrelation
#' function (ACF) on the standardized residuals.
#'
#' @param fit A fitted model object generated by the ctsem package.
#' @param ... Additional arguments to be passed to the \code{\link{ctACF}} function.
#'
#' @return A data table containing the continuous time ACF estimates for standardized residuals.
#'
#' @details This function first extracts the standardized residuals from the fit object using
#' the \code{\link{ctStanKalman}} function. Then, it calculates the continuous time ACF for these residuals
#' and returns the results as a data table.
#'
#' @seealso \code{\link{ctStanKalman}}
#'
#' @examples
#' data.table::setDTthreads(1) #ignore this line
#' # Example usage:
#' ctACFresiduals(ctstantestfit, varnames='Y1',nboot=5)
#'
#' @export
ctACFresiduals <- function(fit,...){
  k=ctResiduals(fit)
  ctACF(k,timecol='Time',idcol='Subject',...)
}

#' Extract Standardized Residuals from a ctsem Fit
#'
#' This function takes a fit object from the ctsem package and extracts the standardized residuals.
#'
#' @param fit A fitted model object generated by the ctsem package.
#'
#' @return A data table containing the standardized residuals for each subject and time point.
#'
#' @details This function uses the \code{\link{ctStanKalman}} function to calculate the standardized residuals
#' and then extracts and formats them as a data table. The standardized residuals represent the differences
#' between the observed and predicted values, divided by the standard errors of the observations.
#'
#' @seealso \code{\link{ctStanKalman}}
#'
#' @examples
#' data.table::setDTthreads(1) #ignore this line
#' # Example usage:
#' residuals <- ctResiduals(ctstantestfit)
#'
#' @export
ctResiduals <- function(fit){
  if(F) Element=NULL
  k = data.table(meltkalman(suppressMessages(ctStanKalman(fit, standardisederrors = TRUE))))
  k = k[Element %in% 'errstdprior' & !is.na(value), ]
  k = dcast(k, formula = formula(Subject + Time ~ Row))
}


