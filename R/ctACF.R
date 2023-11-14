

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
#' @param nboot The number of bootstrap samples for confidence interval estimation (default is 100).
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
  plot=TRUE,timestep='auto',time.max='auto',nboot=100,...){
  
  if(requireNamespace('collapse')){
  } else stop('collapse package needed for ACF!')
  
  if(F) .timediff = ACFhigh= ACFlow= Element= Estimate= SignificanceLevel=  TimeInterval= Variable= boot= ci=NULL
  dat=copy(data.table(dat))
  if(timestep == 'auto'){
    timestep = 0.5 * quantile(dat[,.timediff:= c(NA,diff(get(timecol))),by=idcol][['.timediff']],
      probs=.5,na.rm=TRUE)
    message('Timestep used: ',timestep)
  }
  if(time.max == 'auto')time.max = 10 * quantile(dat[,.timediff:= c(NA,diff(get(timecol))),by=idcol][['.timediff']],probs=.9,na.rm=TRUE)
  
  if(length(varnames)==1 && varnames[1] == 'auto') varnames <- colnames(dat)[!colnames(dat) %in% c(idcol,timecol,'.timediff')]
  if(length(ccfnames )==1 && ccfnames %in% 'all') ccfnames=varnames
  if(length(ccfnames)==1 && is.na(ccfnames)) ccfnames=NULL
  acfnames <- varnames
  varnames <- unique(c(acfnames,ccfnames))
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
  
  # for(vi in v){ #demean via random effects
  # browser()
  #   model <- lmer(formula=paste0(vi,' ~ 1 + (1 | id)'), data = dat)
  #   
  #   # Extract estimated random intercept and SD
  #   random_intercept <- as.numeric(ranef(model)$id[, "(Intercept)"])
  #   
  #   datsd <- dat[,sdvi:=sd(get(vi),na.rm=TRUE),by=id]
  #   datsd <- datsd[,weight:=sum(!is.na(get(vi)))/.N,by=id]
  #   weights <- datsd$weight
  #   model <- lmer(formula='sdvi ~ 1 + (1 | id)', data = datsd,weights=weights)
  #   random_sd <- as.numeric(ranef(model)$id[, "(Intercept)"])
  #   
  #   
  #   # Scale the 'Variable' column using the estimated parameters
  #   dat[, (vi) := (Variable - random_intercept) / random_sd]
  #   }
  # 
  
  
  counter <- 0
  for(vari in acfnames){
    for(varj in unique(c(vari,ccfnames))){
      counter <- counter + 1
      varji <- unique(c(vari,varj))
      # browser()
      vdat <- dat[,c(idcol,timecol,varji),with=FALSE]
      vdat <- vdat[apply(vdat[,(varji),with=F],1,function(x) !all(is.na(x))),] #remove rows with total missingness first
      for(col in varji) set(vdat, j = col, value = as.numeric(vdat[[col]])) #ensure cols are numeric!
      # vdat[apply(.SD,1,function(x) !(all(is.na(x)))), ,.SDcols=varji]
      # vdat<-remove_rows_with_missing(vdat,varji)
      # browser()
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
          g=subjects,lag.max=ceiling(time.max/timestep), 
          plot=F,...)$acf[,1,1] #na.action=na.pass
        if(vari!=varj) acfout <- collapse:::psccf(x = vdat[[vari]][rows],y = vdat[[varj]][rows],
          g=subjects,
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

ctACFquantiles<-function(ACF,probs=c(.025,.975)){
  
  ACF[,ACFlow:=quantile(ACF[Sample!=0],probs = probs[1],na.rm = TRUE),by=TimeInterval]
  ACF[,ACFhigh:=quantile(ACF[Sample!=0],probs = probs[2],na.rm = TRUE),by=TimeInterval]
  ACF<- melt(ACF[Sample==0,colnames(ACF)[colnames(ACF)!='Sample'],with=FALSE],
    measure.vars = c('ACF','ACFlow','ACFhigh'),value.name = 'ACF',variable.name='Estimate')
  ACF[,ci:=grepl('(high)|(low)',Estimate)]
  return(ACF)
}

#' Plot an approximate continuous-time ACF object from ctACF
#'
#' @param ctacfobj object
#' @param df df for the basis spline.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' data.table::setDTthreads(1) #ignore this line
#' # Example usage:
#' head(ctstantestdat)
#' ac=ctACF(ctstantestdat,varnames=c('Y1'),idcol='id',timecol='time',timestep=.5,nboot=5,plot=FALSE)
#' plotctACF(ac)
plotctACF <- function(ctacfobj,df='auto'){
  if(df =='auto'){
    ndt=length(unique(ctacfobj$TimeInterval))
    dfmax=ndt-1
    df=min(c(dfmax,
      ceiling(sqrt(ndt))+5))
  }
      
  gg <- ggplot((ctacfobj),aes(y=ACF,
    # colour=Estimate,
    # linewidth=ci,
    # linetype=ci,
    # group=Estimate,
    x=TimeInterval))+
    scale_linewidth_manual(values = c(1,.1))+
    scale_linetype_manual(values = c('solid','dashed'))+
    # geom_smooth(se=T,method='gam', formula = y ~ s(x, bs = "tp"))+
    stat_quantile(method='rqss', formula=formula(paste0('y ~ bs(x,df=',df,')')), 
      quantiles = .5, linewidth=2)+
    stat_quantile(method='rqss', formula=formula(paste0('y ~ bs(x,df=',df,')')), 
      quantiles = c(.025,.975), linetype='dashed')+
    # geom_point(alpha=.1)+
    guides(linewidth='none',linetype='none')+
    # geom_hline(aes(yintercept=SignificanceLevel),colour='blue',linetype='dotted')+
    # geom_hline(aes(yintercept=-SignificanceLevel),colour='blue',linetype='dotted')+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    ylab('ACF')+
    # coord_cartesian(ylim=c(-1,1))+
    theme_bw()
  
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


