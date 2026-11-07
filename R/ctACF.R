#' Continuous Time Autocorrelation Function (ctACF)
#'
#' Computes autocorrelations and cross correlations by time lag, pooled over subjects, for data
#' observed at irregular times.
#'
#' @param dat The input data in data frame or data table format, in long format.
#' @param varnames Character vector of variable names in the data to compute the ACF for. 'auto' uses all columns that are not time / id.
#' @param ccfnames Character vector of variable names in the data to compute cross correlation for. 'all' uses all variables in varnames, NA uses none.
#' @param idcol The name of the column containing subject IDs (default is 'id').
#' @param timecol The name of the column containing time values (default is 'time').
#' @param plot A logical value indicating whether to create a plot (default is TRUE).
#' @param timestep The width of the time lag bins. 'auto' (default) uses half the median
#'        time interval between consecutive rows of a subject.
#' @param time.max The largest time lag to compute. 'auto' (default) uses two thirds of the
#'        mean, over rows, of the time span of the row's subject.
#' @param nboot The number of bootstrap samples, each resampling subjects with replacement (default is 100).
#' @param scale if TRUE, divide each variable by its within-subject standard deviation.
#' @param center if TRUE, subtract each variable's within-subject mean.
#'
#' @return If 'plot' is TRUE, the function returns a ggplot object of the ACF plot. If 'plot' is
#'         FALSE, it returns a data table with a row for each sample, variable and time lag:
#'         \code{Sample} (0 for the data, then the bootstrap samples), \code{ACF},
#'         \code{TimeInterval}, \code{Variable}, \code{SignificanceLevel} (the 95\% bound of a
#'         correlation of zero computed from that many independent pairs), and \code{Pairs}, the
#'         number of pairs the correlation is computed from.
#'
#' @details Within each subject, every pair of observations is put in the time lag bin nearest
#' to the time between them, at multiples of \code{timestep} up to \code{time.max}. The value at a
#' lag is the correlation of the pairs in its bin, pooled over subjects: of a variable at time
#' t + lag with the same variable at time t for an autocorrelation, and for the cross correlation
#' \code{'x_y'} of x at time t + lag with y at time t, so that at a positive lag y is the earlier.
#' Lags with fewer than three pairs are left out. Autocorrelations are returned from the first
#' non-zero lag, as the lag 0 autocorrelation is always 1. Cross correlations include lag 0, the
#' contemporaneous correlation. Bootstrap samples resample whole subjects, and the plot shows a
#' spline through the estimates with a band from their bootstrap spread (see \code{\link{plotctACF}}).
#'
#' @examples
#' data.table::setDTthreads(1) #ignore this line
#' # Example usage:
#' head(ctstantestdat)
#' ctACF(ctstantestdat,varnames=c('Y1'),idcol='id',timecol='time',nboot=5)
#'
#' @export
ctACF <- function(dat, varnames='auto',ccfnames='all',idcol='id', timecol='time',
  plot=TRUE,timestep='auto',time.max='auto',nboot=100,scale=FALSE, center=FALSE){

  if(F) Sample = ACF = TimeInterval = Variable = SignificanceLevel = Pairs = NULL
  dat <- as.data.table(dat)

  if(length(varnames)==1 && varnames[1] == 'auto') varnames <- colnames(dat)[!colnames(dat) %in% c(idcol,timecol)]
  if(length(ccfnames )==1 && ccfnames %in% 'all') ccfnames=varnames
  if(length(ccfnames)==1 && is.na(ccfnames)) ccfnames=NULL
  acfnames <- varnames
  varnames <- unique(c(acfnames,ccfnames))
  absent <- setdiff(c(idcol,timecol,varnames), colnames(dat))
  if(length(absent)) stop('Not found in the data: ', paste(absent, collapse=', '))
  if(anyNA(dat[[timecol]])) stop('Cannot compute the ACF with missing time data!')
  if(anyNA(dat[[idcol]])) stop('Cannot compute the ACF with missing id data!')

  dat <- copy(dat[, c(idcol,timecol,varnames), with=FALSE])
  for(v in varnames) set(dat, j = v, value = as.numeric(dat[[v]])) #ensure cols are numeric!
  setorderv(dat, c(idcol,timecol))

  # Lags are measured within subjects, from a subject's consecutive rows
  subject <- match(as.character(dat[[idcol]]), unique(as.character(dat[[idcol]])))
  time <- as.numeric(dat[[timecol]])
  within <- subject[-1] == subject[-length(subject)]
  if(identical(timestep, 'auto')){
    timestep = 0.5 * stats::median(diff(time)[within])
    message('Timestep used: ',timestep)
  }
  if(!is.numeric(timestep) || length(timestep) != 1 || !isTRUE(timestep > 0) || !is.finite(timestep))
    stop('timestep must be a single positive number, or auto from subjects with at least two distinct times')
  if(identical(time.max, 'auto')){
    span <- tapply(time, subject, function(x) diff(range(x)))
    time.max = 2/3 * mean(span[subject])
  }
  nlag <- ceiling(time.max/timestep)
  if(!isTRUE(nlag >= 0)) stop('time.max must be a positive number')

  if(scale || center){ #within each subject
    for(v in varnames) set(dat, j = v, value = dat[, {
      x <- get(v)
      m <- if(center) mean(x, na.rm = TRUE) else 0
      s <- if(scale) stats::sd(x, na.rm = TRUE) else 1
      (x - m) / s
    }, by = idcol]$V1) #dat is sorted by id, so the groups come back in row order
  }

  # Rows with no value cannot pair; drop them before pairing
  hasdata <- rowSums(!is.na(as.matrix(dat[, varnames, with = FALSE]))) > 0
  dat <- dat[hasdata,]
  time <- time[hasdata]
  subject <- match(subject[hasdata], unique(subject[hasdata]))
  nsubj <- max(c(0L, subject))

  pairs <- .ctACFPairs(subject = subject, time = time, timestep = timestep, nlag = nlag)

  # Bootstrap samples as counts of each subject drawn; sample 0 is the data
  W <- matrix(1, 1, nsubj)
  if(nboot > 0) W <- rbind(W, t(matrix(vapply(seq_len(nboot), function(b)
    tabulate(sample.int(nsubj, nsubj, replace = TRUE), nsubj), numeric(nsubj)), nrow = nsubj)))

  # Values centred on their overall means: the correlations do not change, and the
  # sums of squares below lose nothing to a large mean
  vals <- lapply(stats::setNames(varnames, varnames), function(v) dat[[v]] - mean(dat[[v]], na.rm = TRUE))

  fullACF <- list()
  for(vari in acfnames){
    for(varj in unique(c(vari,ccfnames))){
      if(vari == varj){ #each pair once, earlier and later
        sel <- pairs$k >= 1
        x <- vals[[vari]][pairs$j[sel]]; y <- vals[[vari]][pairs$i[sel]]
        bin <- pairs$k[sel]; s <- subject[pairs$i[sel]]
        bins <- seq_len(nlag)
      } else { #vari later at positive lags, earlier at negative ones; a row with itself once, at lag 0
        d <- pairs$i != pairs$j
        x <- c(vals[[vari]][pairs$j], vals[[vari]][pairs$i[d]])
        y <- c(vals[[varj]][pairs$i], vals[[varj]][pairs$j[d]])
        bin <- c(pairs$k, -pairs$k[d])
        s <- subject[c(pairs$i, pairs$i[d])]
        bins <- -nlag:nlag
      }
      r <- .ctACFPooledCor(x = x, y = y, subject = s, bin = bin, bins = bins, W = W)
      out <- data.table(Sample = rep(0:nboot, times = length(bins)), ACF = as.vector(r$cor),
        TimeInterval = rep(bins * timestep, each = nboot + 1),
        Variable = ifelse(vari==varj, vari, paste0(vari,'_',varj)),
        Pairs = as.vector(r$n))
      out[, SignificanceLevel := stats::qnorm((1 + 0.95)/2)/sqrt(Pairs)]
      setcolorder(out, c('Sample','ACF','TimeInterval','Variable','SignificanceLevel','Pairs'))
      setorderv(out, c('Sample','TimeInterval'))
      fullACF[[length(fullACF) + 1]] <- out[!is.na(ACF),]
    }
  }
  fullACF <- rbindlist(fullACF)

  if(plot){
    return(plotctACF(fullACF))
  } else return(fullACF)
}

# Every pair of rows of one subject whose time difference rounds to at most nlag
# steps: row indices i (earlier, or the same row) and j, and k, the difference in
# steps rounded to the nearest. Rows are sorted by subject and time, so pairs are
# gathered by row offset, and an offset at which no pair qualifies ends the search.
.ctACFPairs <- function(subject, time, timestep, nlag){
  n <- length(time)
  out <- list()
  for(m in seq_len(n) - 1L){
    i <- seq_len(n - m)
    j <- i + m
    k <- floor((time[j] - time[i]) / timestep + 0.5)
    keep <- subject[i] == subject[j] & k <= nlag
    if(!any(keep)) break
    out[[m + 1L]] <- list(i = i[keep], j = j[keep], k = as.integer(k[keep]))
  }
  if(!length(out)) return(data.table(i = integer(0), j = integer(0), k = integer(0)))
  rbindlist(out)
}

# Pearson correlation of the pairs (x, y) in each bin, pooled over subjects, for
# each row of W, which weights the subjects (the data, then bootstrap draws).
# Returns matrices with a row for each row of W and a column for each bin: the
# correlation, NA from fewer than three pairs or without variation, and the
# number of pairs.
.ctACFPooledCor <- function(x, y, subject, bin, bins, W){
  ok <- !is.na(x) & !is.na(y) & bin %in% bins
  d <- data.table(s = subject[ok], b = match(bin[ok], bins), x = x[ok], y = y[ok])
  d <- d[, list(n = .N, sx = sum(x), sy = sum(y), sxx = sum(x*x), syy = sum(y*y), sxy = sum(x*y)),
    by = c('s','b')]
  tot <- lapply(c(n='n', sx='sx', sy='sy', sxx='sxx', syy='syy', sxy='sxy'), function(st){
    M <- matrix(0, ncol(W), length(bins))
    M[cbind(d$s, d$b)] <- d[[st]]
    W %*% M
  })
  n <- tot$n
  vx <- tot$sxx - tot$sx^2 / n
  vy <- tot$syy - tot$sy^2 / n
  cv <- tot$sxy - tot$sx * tot$sy / n
  r <- cv / sqrt(vx * vy)
  r[n < 3 | !(vx > 0) | !(vy > 0)] <- NA
  list(cor = r, n = n)
}

ctACFquantiles<-function(ctacfobj,quantiles=c(.025,.5,.975),separateLearnRates=FALSE,df='auto'){
  
  if(df =='auto'){
    ndt=length(unique(ctacfobj$TimeInterval))
    dfmax=ndt-1 #max df is limited by the number of time intervals
    df=min(c(dfmax,
      ceiling(ndt/10))) #and 10% of the total observations
  }
  ctacfobj <- data.table(data.frame(ctacfobj))
  
  if(!requireNamespace('qgam')) stop("qgam package required for ACF plots: install.packages('qgam')")
  
  qcols <- paste0('Q',quantiles*100,'%')

  # Estimate the learning rate on one variable combination; if that fails, try the
  # next until one succeeds. A short panel has too few distinct time intervals for
  # the spline basis and every attempt fails -- so keep NULL rather than NA as the
  # "nothing yet" marker, because the failure path below has to be distinguishable
  # from a fitted rate rather than indexed as one.
  learnrate <- NULL
  for(vari in unique(ctacfobj$Variable)){
    if(is.null(learnrate)){
      lr <- try(qgam::tuneLearnFast(data = ctacfobj[Variable %in% vari,],qu = quantiles,
        form = ACF ~ s(TimeInterval,bs='cc',k=df),
        control=list(tol=.Machine$double.eps^0.25,progress=FALSE),
        argGam=list(select=TRUE)), silent=TRUE)
      if(!inherits(lr,'try-error')) learnrate <- lr
    }
  }

  # No learning rate anywhere means no quantile band. Substituting a default rate
  # would return a plausible band with nothing behind it, so return the samples
  # unchanged and let the caller notice the absent columns.
  if(is.null(learnrate)){
    message('Quantile splines could not be fitted (too few distinct time intervals?); ACF samples returned without quantile estimates.')
    return(ctacfobj)
  }

  failed <- character(0)
  for(vari in unique(ctacfobj$Variable)){
    lrv <- learnrate
    if(separateLearnRates){
      lr <- try(qgam::tuneLearnFast(data = ctacfobj[Variable %in% vari,],qu = quantiles,
        form = ACF ~ s(TimeInterval,bs='tp',k=df),
        control=list(tol=.Machine$double.eps^0.1,progress=FALSE),
        argGam=list(select=TRUE,optimizer=c('outer','bfgs'))), silent=TRUE)
      if(!inherits(lr,'try-error')) lrv <- lr #otherwise fall back on the shared rate
    }

    pred <- try({
      qg=qgam::mqgam(data = ctacfobj[Variable %in% vari,],qu = quantiles,
        form = ACF ~ s(TimeInterval,bs='tp',k=df),
        err=lrv$err,
        lsig = lrv$lsig,
        argGam=list(select=TRUE,gamma=2))
      qgam::qdo(obj = qg,qu = quantiles, fun = predict,newdata=ctacfobj[Variable %in% vari,])
    }, silent=TRUE)

    if(inherits(pred,'try-error')){
      failed <- c(failed,vari)
      next #leaves this variable's quantile cells NA
    }
    ctacfobj[Variable %in% vari,(qcols):=pred]
  }

  if(!any(qcols %in% names(ctacfobj))){ #every variable failed
    message('Quantile splines could not be fitted (too few distinct time intervals?); ACF samples returned without quantile estimates.')
    return(ctacfobj)
  }
  if(length(failed)) message('Quantile splines could not be fitted for: ',
    paste(failed,collapse=', '),'.')

  return(ctacfobj)
}

# Precision weighted spline through the full data ACF estimates (Sample 0), each
# lag weighted by the inverse of its bootstrap variance, which is treated as
# known (scale=1). The smoothness is then chosen by whether the estimates depart
# from a smooth curve by more than their own standard errors: a precisely
# estimated spike is followed, while scatter consistent with its standard errors
# is smoothed away. An unweighted spline can only choose one smoothness for
# both, and either flattens the spike or follows the noise. Returns the same
# quantile columns as ctACFquantiles, from the spline estimate and its standard
# error, so plotctACF can draw either.
ctACFweightedSpline <- function(ctacfobj,quantiles=c(.025,.5,.975),df='auto'){
  if(F) Sample = ACF = TimeInterval = Variable = Estimate = SE = NULL
  if(!requireNamespace('mgcv')) stop("mgcv package required for ACF plots: install.packages('mgcv')")
  ctacfobj <- data.table(data.frame(ctacfobj))
  qcols <- paste0('Q',quantiles*100,'%')

  failed <- character(0)
  for(vari in unique(ctacfobj$Variable)){
    vdat <- ctacfobj[Variable %in% vari & !is.na(ACF),]
    se <- vdat[Sample > 0, list(SE=sd(ACF)), by=TimeInterval] #needs bootstrap samples, nboot > 1
    est <- merge(vdat[Sample == 0, list(TimeInterval, Estimate=ACF)], se, by='TimeInterval')
    est <- est[is.finite(SE) & SE > 0,]
    k <- min(nrow(est)-1, if(df == 'auto') 40 else df) #basis can be generous, the penalty does the smoothing

    fit <- NULL
    if(k >= 3) fit <- try(mgcv::gam(Estimate ~ s(TimeInterval, k=k), data=est,
      weights=1/SE^2, scale=1, method='REML'), silent=TRUE)
    if(is.null(fit) || inherits(fit,'try-error')){
      failed <- c(failed,vari)
      next #leaves this variable's quantile cells NA
    }
    newdata <- data.frame(TimeInterval=sort(unique(ctacfobj[Variable %in% vari,TimeInterval])))
    pred <- stats::predict(fit, newdata=newdata, se.fit=TRUE)
    qdat <- data.table(TimeInterval=newdata$TimeInterval)
    for(qi in seq_along(quantiles)) set(qdat, j=qcols[qi],
      value=as.numeric(pred$fit + stats::qnorm(quantiles[qi]) * pred$se.fit))
    rows <- which(ctacfobj$Variable %in% vari)
    for(qc in qcols) set(ctacfobj, i=rows, j=qc,
      value=qdat[[qc]][match(ctacfobj$TimeInterval[rows], qdat$TimeInterval)])
  }

  if(!any(qcols %in% names(ctacfobj))){ #every variable failed
    message('Weighted splines could not be fitted (too few distinct time intervals, or fewer than 2 bootstrap samples?); ACF samples returned without quantile estimates.')
    return(ctacfobj)
  }
  if(length(failed)) message('Weighted splines could not be fitted for: ',
    paste(failed,collapse=', '),'.')

  return(ctacfobj)
}

#' Plot an approximate continuous-time ACF object from ctACF
#'
#' @param ctacfobj object
#' @param df basis dimension for the spline. 'auto' uses up to 40 for \code{method='weighted'},
#' and 10\% of the number of time intervals for \code{method='quantile'}.
#' @param quantiles quantiles to plot.
#' @param separateLearnRates if TRUE, estimate the learning rate for the quantile splines for each combination of variables. Slower but theoretically more accurate.
#' Only used with \code{method='quantile'}.
#' @param reducedXlim if non-zero, n timesteps are removed from the upper and lower end of the x range
#' where the spline estimates are less likely to be reasonable.
#' @param estimateSpline if TRUE, a spline is estimated according to \code{method}, otherwise the samples are simply plotted as lines and the other arguments here are not used.
#' @param method \code{'weighted'} (default) fits a spline through the full data estimates, weighting each
#' time interval by the inverse of its bootstrap variance, so that precisely estimated features such as a
#' sharp peak are followed while noisy estimates are smoothed. Needs \code{nboot > 1} in \code{ctACF}.
#' \code{'quantile'} fits quantile regression splines to the pooled bootstrap samples, which smooths
#' all time intervals alike and can flatten sharp peaks.
#'
#' @details With a spline, the full data estimate at each time interval is shown as a point, the spline
#' as a line, and the band spans the outer quantiles.
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
  separateLearnRates=FALSE, reducedXlim=1,estimateSpline=TRUE,method=c('weighted','quantile')){

  if(F) Sample = ACF = TimeInterval = NULL
  method <- match.arg(method)
  # browser()
  # ctacfobj=melt(ctacfobj,measure.vars = paste0('Q',quantiles*100,'%'),variable.name = 'Quantile',value.name = 'Corr')
  if(estimateSpline){
    if(method == 'weighted') ctacfobj=copy(ctACFweightedSpline(ctacfobj,quantiles=quantiles,df=df))
    if(method == 'quantile') ctacfobj=copy(ctACFquantiles(ctacfobj,quantiles=quantiles,
      separateLearnRates=separateLearnRates,df=df))
    #both spline functions return the samples without the quantile columns when the
    #spline cannot be fitted at all -- then there is no band to draw.
    estimateSpline <- all(paste0('Q',quantiles*100,'%') %in% names(ctacfobj))
    if(!estimateSpline) message('No quantile estimates available; plotting ACF samples instead.')
  }
  if(estimateSpline){
    #drop the ends of the spline; reducedXlim = 0 drops nothing (min() of an empty vector warned)
    if(reducedXlim > 0) ctacfobj <- ctacfobj[abs(TimeInterval) < min(tail(sort(unique(TimeInterval)),reducedXlim)),]
    gg <- ggplot(ctacfobj,aes(y=!!sym(paste0('Q',quantiles[2]*100,'%')),
      x=TimeInterval))+
      geom_ribbon(aes(ymin=!!sym(paste0('Q',quantiles[1]*100,'%')),
        ymax=!!sym(paste0('Q',quantiles[3]*100,'%'))),alpha=.2,linetype='dashed',linewidth=.5)+
      geom_line(linewidth=1)+
      geom_point(data=ctacfobj[Sample == 0 & !is.na(ACF),],aes(y=ACF),size=1)+ #full data estimates
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
#' the \code{\link{ctKalmanArray}} function. Then, it calculates the continuous time ACF for these residuals
#' and returns the results as a data table.
#'
#' @seealso \code{\link{ctKalmanArray}}
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
#' @return A data table containing the standardized residuals for each subject and time point,
#' with subjects labelled by their ids in the data.
#'
#' @details This function uses the \code{\link{ctKalmanArray}} function to calculate the standardized residuals
#' and then extracts and formats them as a data table. The standardized residuals represent the differences
#' between the observed and predicted values, divided by the standard errors of the observations.
#'
#' @seealso \code{\link{ctKalmanArray}}
#'
#' @examples
#' data.table::setDTthreads(1) #ignore this line
#' # Example usage:
#' residuals <- ctResiduals(ctstantestfit)
#'
#' @export
ctResiduals <- function(fit){
  if(F) Element=NULL
  k = data.table(meltkalman(suppressMessages(ctKalmanArray(fit, standardisederrors = TRUE))))
  k = k[Element %in% 'errstdprior' & !is.na(value), ]
  k = dcast(k, formula = formula(Subject + Time ~ Row))
  # ctKalmanArray numbers subjects by position; report them by their ids.
  map <- try(.ctFitIdMap(fit), silent = TRUE)
  if(!inherits(map, 'try-error')){
    ids <- as.character(map$original)
    k$Subject <- factor(ids[match(as.integer(as.character(k$Subject)),
      as.integer(map$new))], levels = unique(ids))
  }
  k
}


