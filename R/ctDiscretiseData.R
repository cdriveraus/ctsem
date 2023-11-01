#' Discretise long format continuous time (ctsem) data to specific timestep.
#'
#' Extends and rounds timing information so equal intervals, according to specified
#' timestep, are achieved. NA's are inserted in other columns as necessary,
#' any columns specified by TDpredNames or TIpredNames have zeroes rather than NA's
#' inserted (because some estimation routines do not tolerate NA's in covariates).
#'
#' @param dlong Long format data
#' @param timestep Positive real value to discretise
#' @param timecol Name of column containing absolute (not intervals) time information.
#' @param idcol Name of column containing subject id variable.
#' @param TDpredNames Vector of column names of any time dependent predictors
#' @param TIpredNames Vector of column names of any time independent predictors
#'
#' @return long format ctsem data.
#' @export
#'
#' @examples
#' long <- ctDiscretiseData(dlong=ctstantestdat, timestep = .1,
#' TDpredNames=c('TD1'),TIpredNames=c('TI1','TI2','TI3'))

ctDiscretiseData <- function(dlong, timestep, timecol = 'time', idcol = 'id', TDpredNames = NULL, TIpredNames = NULL) {
  
  dlong <- data.table(dlong)
  
  if (any(is.na(dlong[[timecol]]))) stop('Cannot discretise with missing time data!')
  if (any(is.na(dlong[[idcol]]))) stop('Cannot discretise with missing id data!')
  
  # Calculate the time offset that minimizes information loss
  offset <- mean(diff(unique(dlong[[timecol]]))) %% timestep
  
  originalrows <- sum(apply(dlong, 1, function(x) sum(!is.na(x)) - 2))
  
  dlong[[timecol]] <- dlong[[timecol]] - offset
  dlong <- melt(dlong, id.vars = c(idcol, timecol))
  dlong <- dlong[!is.na(value),]
  dlong <- dcast(dlong, formula = formula(paste0(idcol, '+', timecol, '~variable')), fun.aggregate = mean, na.rm = TRUE)
  
  dnew <- dlong
  dnew <- dnew[, .(newtime = seq(min(get(timecol)), max(get(timecol)), timestep)), by = idcol]
  setnames(dnew, old = 'newtime', timecol)
  dlong <- merge(dlong, dnew, all = TRUE, by = c(idcol, timecol))
  setorderv(dlong, cols = c(idcol, timecol))
  
  newrows <- sum(apply(dlong, 1, function(x) sum(!is.na(x)) - 2))
  
  if (newrows != originalrows) warning(paste0(originalrows - newrows, ' cells of data removed due to time overlap -- reduce timestep if problematic'))
  
  dlong <- data.frame(dlong)
  dlong[, TDpredNames][is.na(dlong[, TDpredNames])] <- 0 
  dlong[, TIpredNames][is.na(dlong[, TIpredNames])] <- NA
  
  return(dlong)
}


  
