
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
  'A','M','testd','ctsemtestdat','smfnode',
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
  "ACF", "ACFhigh", "ACFlow", "Element", "Estimate", "Model", "Obs", 
  "Row", "Sample", "Time", 
  "TimeInterval", "V2", "Variable", "aic", "ci", "leaveOutN", 
  ".ObsCount", ".splitmedian", ".timerange", "DataType", "Iter", "Split", "Time.interval",
  "Type", "condval", "ll", "var1",
  "lp", "np", "obsValue"))
utils::globalVariables(c(".ObsCol", ".ObsRow", ".splitgroup",
  ".splitscore", "Sig", "colvar", "empirical", "n", "n_empirical",
  "ntime", "q025", "q50", "q975", "rowvar", ".order", "dx", "dy", "dxx",
  "dxy", "dyy", "sxx", "sxy", "syy", "x", "y", ".xgroup", ".xvalue", ".xpos"))

if(1==99){
  `:=` = NULL
  `.` =NULL
  .N = .id = id= . = grp = NULL # due to NSE notes in R CMD check
  tibble() #due to weird hidden tibble requirement in plyr or ggplot?
}






#' ctsem
#'
#' ctsem fits continuous and discrete time dynamic models to longitudinal
#' data: panel data, intensive longitudinal data or a single time series, with
#' observations at any times. The processes can interact, be of higher order
#' or nonlinear, and be measured by continuous, binary, ordinal, count or
#' censored variables. Any parameter can differ between subjects, at one or
#' more levels, and depend on covariates. Estimation is by maximum likelihood,
#' or Bayesian by optimisation or sampling, with an engine written in Julia
#' (the default) or with Stan.
#'
#' A model is specified with \code{\link{ctModel}}, or in lavaan, Mplus or
#' formula syntax with \code{\link{ctModelLavaan}}, \code{\link{ctModelMplus}}
#' and \code{\link{ctModelFormula}}, and fitted with \code{\link{ctFit}}.
#' \code{summary()}, \code{\link{ctDiscretePars}}, \code{\link{ctPredict}} and
#' \code{\link{ctSubjectPars}} describe a fit, and \code{\link{ctScoreTest}},
#' \code{\link{ctFitAdequacy}} and \code{\link{ctModelImprove}} check and
#' improve it. \code{\link{ctJuliaInstall}} sets up the Julia engine.
#'
#' \code{vignette('ctsem')} works through a first analysis, and the website,
#' \url{https://cdriveraus.github.io/ctsem/}, has the vignettes and the help
#' pages grouped by task. For the original OpenMx based ctsem, see the
#' \pkg{ctsemOMX} package. \code{citation('ctsem')} gives the references to
#' cite.
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
#' Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous time
#' structural equation modeling with R package ctsem. Journal of Statistical
#' Software, 77(5), 1-35. \doi{10.18637/jss.v077.i05}
#'
#' Driver, C. C., & Voelkle, M. C. (2018). Hierarchical Bayesian continuous
#' time dynamic modeling. Psychological Methods, 23(4), 774-799.
#' \doi{10.1037/met0000168}
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end

NULL

.onAttach <- function(libname, pkgname) {
  if(interactive()){
    # to show a startup message
    packageStartupMessage("ctsem also changes in time. Getting started: vignette('ctsem') or https://cdriveraus.github.io/ctsem/, citation: citation('ctsem'), questions: https://github.com/cdriveraus/ctsem/discussions")
    
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
#' @param ... values to compare. At least two, and none of them empty.
#'   Comparing a zero-length or NULL value is an error.
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
  # A zero-length or NULL argument used to make every comparison against it
  # vacuously true, because `any(numeric(0) >= tol)` is FALSE. That reported
  # success for a comparison never made, and it is how a test comparing a
  # summary field that does not exist passed for years.
  if (n < 2L) stop('test_isclose needs at least two values to compare.')
  empty <- vapply(values, function(v) length(v) == 0L, logical(1))
  if (any(empty)) stop('test_isclose was given a zero-length or NULL value ',
    'in position(s) ', paste(which(empty), collapse = ', '),
    ': there is nothing to compare.')
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
