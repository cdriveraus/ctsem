
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
utils::globalVariables(c(".ObsCol", ".ObsRow", ".splitgroup",
  ".splitscore", "Sig", "colvar", "empirical", "n", "n_empirical",
  "q025", "q50", "q975", "rowvar"))

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
#' then fitting it to data using \code{\link{ctFit}} (alias \code{\link{ctStanFit}}).
#' For original OpenMx / SEM functionality from the first ctsem versions, use the
#' \pkg{ctsemOMX} package.
#' For most purposes, the Stan-based forms in \pkg{ctsem} are more robust and flexible.
#' For examples, see \code{\link{ctFit}}.
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
  if(interactive()){
    # to show a startup message
    packageStartupMessage("ctsem also changes in time, for manual run ctDocs(), see https://github.com/cdriveraus/ctsem/ for quick start / resources, for citation info run citation('ctsem'), for discussion https://github.com/cdriveraus/ctsem/discussions")
    
    try({
      a=sapply(c('rstan','ctsem'),utils::packageVersion)
      apkgs = data.frame(utils::available.packages(repos = getOption('repos')))
      apkgs = apkgs[apkgs$Package %in% names(a),]
      apkgs = sapply(names(a), function(x){
        v=paste0(a[[which(names(a) %in% x)]],collapse='.')
        utils::compareVersion(a=apkgs$Version[apkgs$Package %in% x],b=v)
      })
      if(any(apkgs > 0)) warning('The following important packages for ctsem are out of date: ', paste0(names(a)[apkgs > 0],collapse=', '))
      
    })
  }
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
