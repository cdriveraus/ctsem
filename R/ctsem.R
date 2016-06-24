#' ctsem
#' 
#' ctsem is an R package for continuous time structural equation modelling of panel (N > 1) 
#' and time series (N = 1) data, using either full information maximum likelihood (FIML) or the Kalman filter. Most 
#' dynamic models for longitudinal data in the social and behavioural sciences are discrete time models. An assumption of 
#' discrete time models is that time intervals between measurements are equal, and that all subjects were assessed at the same 
#' intervals. Violations of this assumption are regularly ignored due to the difficulty of accounting for varying time intervals, 
#' therefore parameter estimates can be severely biased. By using stochastic differential equations and estimating an underlying 
#' continuous process, continuous time models allow for any pattern of measurement occasions. By interfacing to a general purpose 
#' SEM package (OpenMx), ctsem combines the flexible specification of structural equation models with the enhanced 
#' data gathering opportunities and improved estimation of continuous time models. ctsem can estimate relationships over 
#' time for multiple latent processes, measured by multiple noisy indicators with varying time intervals between observations. 
#' Within and between effects are estimated simultaneously by modelling both observed covariates and unobserved heterogeneity. 
#' Exogenous shocks with different shapes, group differences, higher order diffusion effects and oscillating processes can all 
#' be simply modelled. To use ctsem, one first specifies a model using \code{\link{ctModel}}, fits this to wide format data using \code{\link{ctFit}}, then 
#' \code{plot} (\code{\link{plot.ctsemFit}}) and \code{summary} (\code{\link{summary.ctsemFit}}) methods are available to analyse the fitted object.  
#' \code{\link{ctMultigroupFit}} may be used in place of \code{\link{ctFit}} to specify a multi group model.
#' For examples, see \code{\link{ctFit}}. For more detailed information, see the vignette by running: \code{vignette('ctsem')}
#'  
#' @docType package
#' @name ctsem
NULL