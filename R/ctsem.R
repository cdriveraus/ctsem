utils::globalVariables(c("invDRIFT","II","DRIFTexp","vec2diag","diag2vec",
  "mxData","mxMatrix","mxAlgebra","MANIFESTVARbase","MANIFESTVARcholdiag",
  "MANIFESTVARchol","T0VARbase","T0VARcholdiag","T0VARchol","DIFFUSIONbase",
  "DIFFUSIONcholdiag","DIFFUSIONchol","invDRIFTHATCH","cvectorize","DRIFTHATCH",
  "TRAITVARbase","TRAITVARcholdiag","TRAITVARchol","MANIFESTTRAITVARbase",
  "MANIFESTTRAITVARcholdiag","MANIFESTTRAITVARchol","mxComputeSequence",
  "mxComputeGradientDescent","mxComputeReportDeriv","TDPREDVARbase",
  "TDPREDVARcholdiag","TDPREDVARchol","TIPREDVARbase","TIPREDVARcholdiag",
  "TIPREDVARchol","mxExpectationRAM","mxFitFunctionML","Ilatent","Alatent",
  "Amanifestcov","invIminusAlatent","Smanifest","Amanifest","Mmanifest",
  "mxExpectationNormal","omxSelectRowsAndCols","expCov","existenceVector",
  "omxSelectCols","expMean","log2pi","numVar_i","filteredExpCov","%&%",
  "filteredDataRow","filteredExpMean","firstHalfCalc","secondHalfCalc",
  "rowResults","mxFitFunctionRow","TdpredNames","discreteCINT_T1","discreteDRIFT_T1",
  "discreteDIFFUSION_T1","mxExpectationStateSpace","mxExpectationSSCT","ctsem.fitfunction",
  "ctsem.penalties","FIMLpenaltyweight","ctsem.simpleDynPenalty","ieigenval",
  "mxFitFunctionAlgebra","mxCI","mxComputeConfidenceInterval","DRIFT",
  "n.latent","DIFFUSION","TRAITVAR","n.TDpred","TDPREDEFFECT","TDPREDMEANS",
  "TDPREDVAR","TRAITTDPREDCOV","n.TIpred","TIPREDEFFECT","TIPREDMEANS",
  "TIPREDVAR","CINT","n.manifest","LAMBDA","MANIFESTMEANS","MANIFESTVAR",
  "mxFitFunctionMultigroup", "asymDIFFUSION", 'data.id',
  'filteredExpCovchol','filteredExpCovcholinv',
  'A','M','testd',
  'T0VAR','T0MEANS',
  'TDpredNames', 'TIpredNames', 'Tpoints', 'extract', 'latentNames', 'manifestNames',
  'plot', 'points','T0TRAITEFFECT'))

#' ctsem
#' 
#' ctsem is an R package for continuous time structural equation modelling of panel (N > 1) 
#' and time series (N = 1) data, using either a frequentist or Bayesian approach. 
#' The frequentist approach is faster but can only estimate random-effects on the intercepts, 
#' while the Bayesian approach allows for random-effects across all model parameters. 
#' 
#' The general workflow begins by specifying a model using the \code{\link{ctModel}} function, 
#' in which the \code{type} of model is also specified. Then the model is fit to data using 
#' either \code{\link{ctFit}} if an 'omx' (OpenMx, frequentist) model is specified or 
#' \code{\link{ctStanFit}} if a 'stanct' or 'standt' (Stan, continuous / discrete time, Bayesian)
#' model is specified. 
#' For examples, see either \code{\link{ctFit}} or \code{\link{ctStanFit}}. 
#' For more detailed information, see the frequentist vignette by running: \code{vignette('ctsem')} 
#' For citation info, please run \code{citation('ctsem')} .
#'  
#' @docType package
#' @name ctsem
#' @import grDevices methods stats rstan OpenMx graphics plyr
NULL
