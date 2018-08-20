#' Update an already compiled and fit ctStanFit object
#' 
#' Allows one to change model elements that don't require recompiling, then re fit.
#'
#' @param fit 
#' @param datalong 
#' @param ctstanmodel 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
ctStanUpdModel <- function(fit, datalong, ctstanmodel,...){
  
  new <-ctStanFit(datalong = datalong, ctstanmodel = ctstanmodel,fit=FALSE,...)
  
  fit$standata <- new$standata
  fit$data <- new$data
  fit$setup <- new$setup
  fit$args <- match.call
  return(fit)
}
