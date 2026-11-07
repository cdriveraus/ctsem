#' Discretise long format continuous time (ctsem) data to specific timestep.
#'
#' Rounds timing information so equal intervals, according to specified
#' timestep, are achieved. Each subject's times are moved to the nearest point of a
#' grid that starts at the subject's first time; values that land on the same point
#' are averaged, with a warning. NA's are inserted in other columns as necessary.
#' Any columns specified by TDpredNames have zeroes rather than NA's inserted
#' (because some estimation routines do not tolerate NA's in covariates), while
#' columns specified by TIpredNames are left as NA's -- since a time independent
#' predictor is constant within a subject, inserting zeroes for it would be incorrect.
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

  if(F) .origin = NULL
  dlong <- copy(data.table(dlong))

  if (any(is.na(dlong[[timecol]]))) stop('Cannot discretise with missing time data!')
  if (any(is.na(dlong[[idcol]]))) stop('Cannot discretise with missing id data!')
  if (!is.numeric(timestep) || length(timestep) != 1 || !isTRUE(timestep > 0) || !is.finite(timestep))
    stop('timestep must be a single positive number')

  valuecols <- setdiff(colnames(dlong), c(idcol, timecol))
  originalcells <- sum(!is.na(as.matrix(dlong[, valuecols, with = FALSE])))

  # Times are counted in whole steps from each subject's first time until the end.
  # A grid built as a sequence of doubles does not meet the observed times
  # exactly, so merging on them left observations beside grid rows a rounding
  # error away, and the intervals unequal.
  dlong[, .origin := min(get(timecol)), by = idcol]
  origins <- unique(dlong[, c(idcol, '.origin'), with = FALSE])
  dlong[, (timecol) := floor((get(timecol) - .origin) / timestep + 0.5)]
  dlong[, .origin := NULL]
  dlong <- melt(dlong, id.vars = c(idcol, timecol))
  dlong <- dlong[!is.na(value),]
  dlong <- dcast(dlong, formula = formula(paste0(idcol, '+', timecol, '~variable')), fun.aggregate = mean, na.rm = TRUE, fill = NA)

  dnew <- dlong[, list(.step = seq(min(get(timecol)), max(get(timecol)))), by = idcol]
  setnames(dnew, old = '.step', timecol)
  dlong <- merge(dlong, dnew, all = TRUE, by = c(idcol, timecol))
  dlong <- merge(dlong, origins, by = idcol)
  dlong[, (timecol) := .origin + get(timecol) * timestep]
  dlong[, .origin := NULL]
  setorderv(dlong, cols = c(idcol, timecol))

  newcells <- sum(!is.na(as.matrix(dlong[, setdiff(colnames(dlong), c(idcol, timecol)), with = FALSE])))

  if (newcells != originalcells) warning(paste0(originalcells - newcells, ' cells of data removed due to time overlap -- reduce timestep if problematic'))

  dlong <- data.frame(dlong)
  dlong[, TDpredNames][is.na(dlong[, TDpredNames])] <- 0
  dlong[, TIpredNames][is.na(dlong[, TIpredNames])] <- NA

  return(dlong)
}
