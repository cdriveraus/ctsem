#' Continuous Time Autocorrelation Function (ctACF)
#'
#' This function computes an approximate continuous time autocorrelation function (ACF) for data
#' containing multiple subjects and/or variables.
#'
#' @param dat The input data in data frame or data table format.
#' @param varnames Character vector of variable names in the data to compute the ACF for. 'auto' uses all columns that are not time / id. 
#' @param idcol The name of the column containing subject IDs (default is 'id').
#' @param timecol The name of the column containing time values (default is 'time').
#' @param plot A logical value indicating whether to create a plot (default is TRUE).
#' @param timestep The time step for discretizing data. 'auto' to automatically determine
#'        the timestep based on data distribution (default is 'auto'). 
#'        In this case the timestep is computed as 1/10th of the 10th percentile for time intervals in the data.
#' @param time.max The maximum time lag to compute the ACF (default is 10). If 'auto', is set to 10 times the 90th percentile interval in the data.
#' @param nboot The number of bootstrap samples for confidence interval estimation (default is 100).
#' @param ... additional arguments (such as demean=FALSE) to pass to the \code{\link{stats::acf}} function.
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
#' # Example usage:
#' head(ctstantestdat)
#' ctACF(ctstantestdat,varnames=c('Y1','Y2'),idcol='id',timecol='time',nboot=10)
#'
#' @export
ctACF <- function(dat, varnames='auto',idcol='id', timecol='time',plot=TRUE,timestep='auto',time.max='auto',nboot=100,...){

  dat=copy(data.table(dat))
  if(timestep == 'auto') timestep = 0.1 * quantile(dat[,.timediff:= c(NA,diff(get(timecol))),by=idcol][['.timediff']],probs=.1,na.rm=TRUE)
  if(time.max == 'auto') time.max = 10 * quantile(dat[,.timediff:= c(NA,diff(get(timecol))),by=idcol][['.timediff']],probs=.9,na.rm=TRUE)
  

  
  if(length(varnames)==1 && varnames[1] == 'auto') varnames <- colnames(dat)[!colnames(dat) %in% c(idcol,timecol,'.timediff')]
  
  # Function to remove rows with missing values before and after non-missing values
  remove_rows_with_missing <- function(dt, col_name) {
    dt[, {
      if (all(is.na(.SD[[col_name]]))) {
        .SD[.N]  # Return a single row with NA values if all values are missing
      } else {
        first_non_missing <- min(which(!is.na(.SD[[col_name]])))
        last_non_missing <- max(which(!is.na(.SD[[col_name]])))
        .SD[first_non_missing:last_non_missing, ]
      }
    }, by = idcol]
  }

  
  for(vari in varnames){
    
    vdat<-remove_rows_with_missing(dat[,c(idcol,timecol,vari),with=FALSE],vari)
    datt=vdat[,tail(.SD,1),by=idcol] #adding distant time point (outside acf range) to ensure sufficient padding with NA's for stacking (avoiding attenuation of low lag corr)
    datt[,(timecol):=get(timecol)+time.max]
    datt <- datt[,c(idcol,timecol),with=FALSE]
    vdat=merge(vdat,datt,by = c(idcol,timecol),all=TRUE)
    
    vdat = data.table(ctDiscretiseData(vdat,timecol = timecol,idcol = idcol,timestep=timestep))
    

    ACF=lapply(0:nboot,function(x) {
      message(sprintf(paste0("\r ",vari," %3d%%"), as.integer(x/nboot*100)),appendLF = FALSE)
      subjects <- sample(unique(vdat[[idcol]]),replace=TRUE)
      rows <- unlist(sapply(subjects,function(x) which(vdat[[idcol]] %in% x)))
      if(x==0) acfout <- acf(vdat[[vari]], lag.max=ceiling(time.max/timestep), plot=F, na.action=na.pass,...)$acf[,1,1]
      if(x > 0) acfout <- acf(vdat[[vari]][rows], lag.max=ceiling(time.max/timestep), plot=F, na.action=na.pass,...)$acf[,1,1]
      out <- data.table(ACF=acfout)
      out[,TimeInterval:=seq(0,length.out=.N,by=timestep)]
    })
    message('')
    ACF=rbindlist(ACF,idcol = 'boot')
    ACF[,boot:=boot-1]
    ACF[,ACFlow:=quantile(ACF[boot!=0],probs = .025,na.rm = TRUE),by=TimeInterval]
    ACF[,ACFhigh:=quantile(ACF[boot!=0],probs = .975,na.rm = TRUE),by=TimeInterval]
    ACF<- melt(ACF[boot==0,colnames(ACF)[colnames(ACF)!='boot'],with=FALSE],
      id.vars = c('TimeInterval'),value.name = 'ACF',variable.name='Estimate')
    ACF[,SignificanceLevel:= qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(dat[[vari]])))]
    ACF[,ci:=grepl('(high)|(low)',Estimate)]
    ACF=ACF[TimeInterval!=0,]
    ACF[,Variable:=vari]
    if(vari!=varnames[1]) fullACF <- rbind(fullACF,ACF) else fullACF <- ACF
  }
  
  
  if(plot){
    gg <- ggplot(fullACF,aes(y=ACF,
      colour=Estimate,
      linewidth=ci,
      linetype=ci,
      group=Estimate,
      x=TimeInterval))+
      scale_linewidth_manual(values = c(1,.1))+
      scale_linetype_manual(values = c('solid','dashed'))+
      geom_smooth(se=F,colour='red',method='gam', formula = y ~ s(x, bs = "tp"))+
      # geom_point()+
      guides(linewidth='none',linetype='none')+
      geom_hline(aes(yintercept=SignificanceLevel),colour='blue',linetype='dotted')+
      geom_hline(aes(yintercept=-SignificanceLevel),colour='blue',linetype='dotted')+
      geom_hline(yintercept=0)+
      geom_vline(xintercept=0)+
      ylab('ACF')+
      coord_cartesian(ylim=c(-1,1))+
      theme_bw()
    
    if(length(varnames) > 1) gg <- gg+facet_wrap(vars(Variable))
    return(gg)
  } else return(fullACF)
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
#' # Example usage:
#' ctACFresiduals(ctstantestfit,nboot=10)
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
#' # Example usage:
#' residuals <- ctResiduals(ctstantestfit)
#'
#' @export
ctResiduals <- function(fit){
  k = data.table(meltkalman(suppressMessages(ctStanKalman(fit, standardisederrors = TRUE))))
  k = k[Element %in% 'errstdprior' & !is.na(value), ]
  k = dcast(k, formula = formula(Subject + Time ~ Row))
}

  
