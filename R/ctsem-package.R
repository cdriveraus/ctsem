
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
  'A','M','testd','ctstantestdat','smfnode',
  'T0VAR','T0MEANS', 'MANIFESTTRAITVAR',
  'TDpredNames', 'TIpredNames', 'Tpoints', 'extract', 'latentNames', 'manifestNames',
  'plot', 'points','T0TRAITEFFECT',
  'T0VARsubindex','DRIFTsubindex','DIFFUSIONsubindex','CINTsubindex','.',
  'Var1',
  'Var2',
  'value',
  'WhichObs',
  'variable','Datasource','id','verbose','param',
  'manifest',
  'Original','sysnoise','starts', 'obsNames',
  'parlp',
  "ACF", "ACFhigh", "ACFlow", "Element", "Estimate", "Model", "NobsDT", "Obs", 
  "ObsVsGenID", "ObsVsGenRow", "OutOf95Row", "Row", "Sample", "Time", 
  "TimeInterval", "V2", "Variable", "aic", "ci", "highdat", "leaveOutN", 
  ".ObsCount", ".splitmedian", ".timerange", "DataType", "Iter", "Split", "Time.interval",
  "Type", "condval", "ll", "var1",
  "lowdat", "lp", "mediandat", "mediandatRank", "np", "obsValue"))

if(1==99){
  `:=` = NULL
  `.` =NULL
  .N = .id = id= . = grp = NULL # due to NSE notes in R CMD check
  tibble() #due to weird hidden tibble requirement in plyr or ggplot?
}






#' ctsem
#' 
#' ctsem is an R package for continuous time structural equation modelling of panel (N > 1) 
#' and time series (N = 1) data, using either a frequentist or Bayesian approach, or middle
#' ground forms like maximum a posteriori. 
#'  
#' The general workflow begins by specifying a model using the \code{\link{ctModel}} function, 
#' in which the \code{type} of model is also specified. Then the model is fit to data using 
#' \code{\link{ctStanFit}}. The ctFit function which allows for fitting using the OpenMx / SEM form,
#' as described in the original JSS ctsem paper, can now be found in the ctsemOMX package.  
#' The omx forms are no longer in 
#' development and for most purposes, the newer stan based forms are more robust and flexible.
#' For examples, see  \code{\link{ctStanFit}}. 
#' For citation info, please run \code{citation('ctsem')} .
#'  
#' @import grDevices methods stats graphics data.table ggplot2
#' @import Rcpp
#' @importFrom RcppParallel CxxFlags RcppParallelLibs
#' @importFrom tibble tibble
#' @import expm expm
#' @importFrom rstan constrain_pars sampling unconstrain_pars stan_model log_prob monitor get_num_upars stanc get_sampler_params As.mcmc.list monitor
#' @importFrom plyr aaply alply round_any
#' @importFrom utils relist as.relistable tail capture.output
#' @importFrom Deriv Simplify 
#' @importFrom cOde jacobianSymb prodSymb
#' @importFrom splines bs
#' @useDynLib ctsem, .registration = TRUE
#' 
#' @references 
#' https://www.jstatsoft.org/article/view/v077i05
#' 
#' Driver, C. C., & Voelkle, M. C. (2018). Hierarchical Bayesian continuous time dynamic modeling. 
#' Psychological Methods. Advance online publication.http://dx.doi.org/10.1037/met0000168
#' 
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.3. http://mc-stan.org
#' 
#' #' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end

NULL

.onAttach <- function(libname, pkgname) {
  # to show a startup message
  packageStartupMessage("ctsem also changes in time, for manual run ctDocs(), for blog see https://cdriver.netlify.app/, for citation info run citation('ctsem'), for original OpenMx functionality install.packages('ctsemOMX'), and for discussion https://github.com/cdriveraus/ctsem/discussions")
  
  try({
    a=sapply(c('rstan','ctsem'),utils::packageVersion)
    apkgs = data.frame(utils::available.packages(repos = 'https://cloud.r-project.org'))
    apkgs = apkgs[apkgs$Package %in% names(a),]
    apkgs = sapply(names(a), function(x){
      v=paste0(a[[which(names(a) %in% x)]],collapse='.')
      utils::compareVersion(a=apkgs$Version[apkgs$Package %in% x],b=v)
    })
    if(any(apkgs > 0)) warning('The following important packages for ctsem are out of date: ', paste0(names(a)[apkgs > 0],collapse=', '))
    
  })
  
}

#' Get documentation pdf for ctsem
#'
#' @return Nothing. Opens a pdf when run interactively.
#' @export
#'
#' @examples
#' ctDocs()
ctDocs <- function(){
  if(interactive()){
    r=runif(1,0,9999999)
    pdfpath=file.path(tempdir(),paste0('/ctsemManual_',r,'.pdf'))
    utils::download.file(url="https://github.com/cdriveraus/ctsem/raw/master/vignettes/hierarchicalmanual.pdf",
      destfile=pdfpath,mode='wb')
    try(openPDF(pdfpath))
  }
}

#' ctFit function placeholder
#' 
#' For the original ctsem OpenMx functionality, the package ctsemOMX should be loaded.
#'
#' @param ... arguments to pass to ctFit, if ctsemOMX is loaded.
#'
#' @return message or fit object.
#' @export
#'
#' @examples
#' \donttest{
#' data(AnomAuth) 
#' AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
#'   Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL) 
#' AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)
#' }
ctFit <- function(...){
  if('ctsemOMX' %in% utils::installed.packages()[,1]){
    if(!'ctsemOMX' %in% (.packages())){
      message('for original ctsem functionality using OpenMx, please use: library(ctsemOMX)')
    } else message('call ctFit from ctsemOMX package, e.g. ctsemOMX::ctFit(...)')
  } else message('For original ctsem functionality using OpenMx, install.packages("ctsemOMX")')  
}

#' Tests if 2 values are close to each other
#'
#' @param ... values to compare
#' @param tol tolerance
#'
#' @return Logical or testthat output.
#' @export
#'
#' @examples
#' test_isclose(1,1.0000001, tol=1e-4)
test_isclose <- function(..., tol = 1e-8) {
  values <- list(...)
  n <- length(values)
  out <- TRUE
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (any(abs(values[[i]] - values[[j]]) >= tol)) {
        out <- FALSE
      }
    }
  }
  
  if(!requireNamespace("testthat", quietly = TRUE)){
    message("testthat not installed, returning logical")
    return(out)
  } else return(testthat::expect_true(out))
}
