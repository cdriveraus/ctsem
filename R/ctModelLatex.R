ctModelBuildPopCov <- function(ctm,linearise){ #for latex
  ctm <- T0VARredundancies(ctm)
  ctm$pars <- ctStanModelCleanctspec(ctm$pars)
  freepars <- !is.na(ctm$pars$param) &
    !grepl('\\W', gsub('.', '', ctm$pars$param, fixed=TRUE)) &
    !ctm$pars$param %in% ctm$latentNames
  pars <- unique(ctm$pars$param[ctm$pars$indvarying & freepars])
  d=length(pars)
  m <- matrix(paste0(ifelse(linearise,'','raw'),'PCov_',rep(1:d,d),'_',rep(1:d,each=d)),d,d,dimnames = list(pars,pars))
  m[upper.tri(m)]=t(m)[upper.tri(m)]
  return(m)
}

# getPopEffectsFromFit <- function(x,linearise=TRUE,digits=3){
#   # browser()
#   ms=ctMatsetupFreePars(x$setup$matsetup)
#   e=x$stanfit$transformedparsfull#ctExtract(x)
#   if(x$standata$ntipred > 0){
#     if(linearise) timat <- e$linearTIPREDEFFECT
#     if(!linearise) timat <- e$TIPREDEFFECT
#     dimnames(timat) <- list(iter=1:(dim(timat)[1]),param= ms$parname,
#       TIpred= x$ctstanmodel$TIpredNames)
#     timat <- timat[,apply(x$standata$TIPREDEFFECTsetup,1,function(x) any(x!=0)),,drop=FALSE]
#   } else timat <- diag(0,0)
#   
#   if(!is.null(e$rawpopc) && !is.null(e$popcov)){ 
#     if(!linearise) popcov <- matrix(e$rawpopc[1,4,,],dim(e$rawpopc)[3])
#     if(linearise) {
#       popcov <- matrix(e$popcov[1,,],dim(e$popcov)[3])
#       if(x$standata$intoverpop==1){
#         t0index <- ms$indvarying[ms$param > 0 & ms$row <= x$standata$nlatent & ms$matrix %in% 1 & ms$indvarying > 0]
#         popcov[t0index,t0index] <- e$pop_T0VAR[1,
#           t0index,t0index] #is this correct...?
#       }
#     }
#     dimnames(popcov) <- list(#iter=1:(dim(popcov)[1]),
#       ms$parname[as.logical(ms$indvarying)],
#       ms$parname[as.logical(ms$indvarying)] )
#     
#   }
#   return(list(popcov=popcov,tieffects=timat))
# }

ctModelBuildTIeffects <- function(ctm){ #for latex
  ctm$pars <- ctStanModelCleanctspec(ctm$pars)
  tieffects <- unique(colnames(ctm$pars)[grep('_effect',colnames(ctm$pars),fixed=TRUE)])
  freepars <- !is.na(ctm$pars$param) &
    !grepl('\\W', gsub('.', '', ctm$pars$param, fixed=TRUE)) &
    !ctm$pars$param %in% ctm$latentNames
  pars <- unique(ctm$pars$param[
    freepars & apply(ctm$pars[,tieffects,drop=FALSE],1,any)])
  timat <- matrix(0,length(pars),length(tieffects),dimnames = list(pars,gsub('_effect','',tieffects)))
  if(length(tieffects)){
    for(p in seq_along(pars)){
      timat[p,] <- unlist(ctm$pars[match(x = pars[p],ctm$pars$param),tieffects,drop=FALSE])
    }
  }
  for(i in seq_len(nrow(timat))){
    for(j in seq_len(ncol(timat))){
      if(timat[i,j] !=0) timat[i,j] = paste0(pars[i],'_',gsub('_effect','',tieffects[j],fixed=TRUE))
    }}
  return(timat)
}

ctModelLatexT0Labels <- function(ctm){
  nlatent <- ctm$n.latent
  labels <- paste0('T0MEANS_', ctm$latentNames[seq_len(nlatent)])
  t0pars <- ctm$pars[ctm$pars$matrix %in% 'T0MEANS' &
      ctm$pars$row <= nlatent & ctm$pars$col == 1,,drop=FALSE]
  t0pars <- t0pars[match(seq_len(nlatent), t0pars$row),,drop=FALSE]
  named <- !is.na(t0pars$param) & nzchar(as.character(t0pars$param))
  labels[named] <- as.character(t0pars$param[named])
  labels
}

ctModelLatexT0DisplayCov <- function(t0var, digits=3){
  t0var <- as.matrix(t0var)
  d <- nrow(t0var)
  if(d < 1) return(matrix(numeric(0), 0, 0))
  num <- suppressWarnings(matrix(as.numeric(t0var), nrow=d, ncol=d))
  if(!any(is.na(num))) return(round(sdpcor2cov(num), digits))
  
  out <- matrix(0, d, d)
  dimnames(out) <- dimnames(t0var)
  for(ri in seq_len(d)){
    for(ci in seq_len(d)){
      if(ri == ci){
        out[ri,ci] <- if(!is.na(num[ri,ri]) && num[ri,ri] == 0) 0 else
          paste0('T0cov_', ri, '_', ci)
      } else {
        lr <- max(ri, ci)
        lc <- min(ri, ci)
        out[ri,ci] <- if(!is.na(num[lr,lc]) && num[lr,lc] == 0) 0 else
          paste0('T0cov_', ri, '_', ci)
      }
    }
  }
  out[upper.tri(out)] <- t(out)[upper.tri(out)]
  out
}

ctModelLatexAugmentT0 <- function(popmeans, popcov, timat, ctm, digits=3,
  t0means=NULL, t0cov=NULL, latentPopNames=NULL, t0covIsTotal=FALSE){
  
  nlatent <- ctm$n.latent
  if(nlatent < 1) return(list(popmeans=popmeans, popcov=popcov, timat=timat))
  labels <- ctModelLatexT0Labels(ctm)
  if(is.null(t0means)) t0means <- as.matrix(ctm$T0MEANS)[seq_len(nlatent),1]
  if(is.null(t0cov)) t0cov <- ctModelLatexT0DisplayCov(
    as.matrix(ctm$T0VAR)[seq_len(nlatent), seq_len(nlatent), drop=FALSE],
    digits=digits)
  t0cov <- as.matrix(t0cov)
  
  labelsInPopcov <- labels %in% rownames(popcov)
  
  pars <- unique(c(rownames(popcov), labels))
  newpopcov <- matrix(0, length(pars), length(pars), dimnames=list(pars, pars))
  if(nrow(popcov) > 0) {
    ii <- match(rownames(popcov), rownames(newpopcov))
    newpopcov[ii, ii] <- popcov
  }
  if(t0covIsTotal) {
    newpopcov[labels, labels] <- t0cov[seq_len(nlatent), seq_len(nlatent), drop=FALSE]
  } else if(any(!labelsInPopcov)){
    missingLabels <- labels[!labelsInPopcov]
    missingIndex <- which(!labelsInPopcov)
    newpopcov[missingLabels, missingLabels] <- t0cov[missingIndex, missingIndex, drop=FALSE]
  }
  
  if(!is.null(latentPopNames) && length(latentPopNames) > 0 &&
      nrow(t0cov) >= nlatent + length(latentPopNames)){
    for(pi in seq_along(latentPopNames)){
      pname <- latentPopNames[pi]
      if(pname %in% rownames(newpopcov)){
        newpopcov[labels, pname] <- t0cov[seq_len(nlatent), nlatent + pi]
        newpopcov[pname, labels] <- t0cov[nlatent + pi, seq_len(nlatent)]
      }
    }
  }
  popcov <- newpopcov
  
  if(length(popmeans) > 0){
    newpopmeans <- rep(NA, nrow(popcov))
    names(newpopmeans) <- rownames(popcov)
    oldnames <- names(popmeans)
    if(is.null(oldnames)) oldnames <- rownames(popcov)[seq_along(popmeans)]
    oldmatch <- match(oldnames, names(newpopmeans))
    newpopmeans[oldmatch[!is.na(oldmatch)]] <- popmeans[!is.na(oldmatch)]
    newpopmeans[labels] <- t0means
    popmeans <- newpopmeans
  }
  
  if(nrow(timat) > 0){
    newtimat <- matrix(0, nrow(popcov), ncol(timat),
      dimnames=list(rownames(popcov), colnames(timat)))
    newtimat[rownames(timat),] <- timat
    timat <- newtimat
  }
  
  list(popmeans=popmeans, popcov=popcov, timat=timat)
}

ctMatsetupFreePars <- function(m,intoverpop){
  m=m[m$when %in% c(0,-1) & m$param > 0,,drop=FALSE]
  m=m[match(unique(m$param),m$param),,drop=FALSE]
  m = m[order(m$param),,drop=FALSE]
}

texPrep <- function(x){ #replaces certain characters with tex safe versions
  for(i in 1:length(x)){
    x[i]=gsub('_', '\\_',x[i],fixed=TRUE)
    x[i]=gsub('^', '\\textasciicircum',x[i],fixed=TRUE)
  }
  return(x)
}

# Rendering a fitted model's context-dependent cells ---------------------------
#
# For an ordinary cell, substituting the fitted number is exactly right. For a
# cell written as an expression referencing a latent process or a time dependent
# predictor it is not: collapsing `-log1p(exp(dr11)) * eta2` to a single number
# produces the equation of a *linear* model, and a reader documenting a
# nonlinear fit from this output would publish the wrong model.
#
# So such a cell keeps its structure and loses only its labels: each free
# parameter is replaced by its estimated value, latent references become
# subscripted etas and TD predictor references keep their names. The result
# reads as the model that was fitted -- numbers where the estimates are,
# structure where the structure is.

# The estimated raw value of every free parameter, by name. Raw, not
# transformed: inside a cell expression a parameter appears before its
# transform has been applied, so the transformed value would be wrong there.
.ctLatexRawEstimates <- function(fit, e = NULL) {
  if (inherits(fit, 'ctJuliaFit')) {
    cells <- .ctBackendFreeParameterCells(fit)
    values <- as.numeric(fit$estimate$raw)[cells$parnumber]
    names(values) <- .ctBackendParameterNames(cells)
    return(values[!is.na(names(values))])
  }
  if (is.null(e)) e <- ctExtract(fit)
  values <- as.numeric(ctCollapse(e$rawpopmeans, 1, mean))
  names(values) <- getparnames(fit)
  values[!is.na(names(values))]
}

# The expression a cell was actually written as.
#
# Two spec syntaxes reach the same place. A bare expression lives in `param`
# and `listOfMatrices` already returns it. A `name | transform` cell keeps only
# the name there, with the expression in `transform` and the name written as
# `param` inside it -- so the expression has to be reassembled before it can be
# rendered, or the cell collapses to a single estimate again.
.ctLatexCellExpression <- function(pars, mi, i, j) {
  row <- which(pars$matrix %in% mi & pars$row %in% i & pars$col %in% j)
  if (!length(row)) return(NA_character_)
  row <- row[1L]
  transform <- pars$transform[row]
  written <- !is.na(transform) && nzchar(as.character(transform)) &&
    is.na(suppressWarnings(as.numeric(transform)))
  if (written) return(gsub('\\bparam\\b', as.character(pars$param[row]),
    as.character(transform), perl = TRUE))
  as.character(pars$param[row])
}

# One cell's expression, with estimates for labels and math for references.
.ctLatexRenderExpression <- function(expression, estimates, latentNames,
  TDpredNames, digits = 3) {

  if (length(expression) != 1 || is.na(expression) || !nzchar(expression)) return(expression)
  quoted <- function(x) paste0('\\b\\Q', x, '\\E\\b')

  # Longest names first, so that a parameter called `a1` cannot eat part of
  # `a12` before `a12` has had its turn.
  named <- names(estimates)[order(nchar(names(estimates)), decreasing = TRUE)]
  for (nm in named) {
    expression <- gsub(quoted(nm), format(round(estimates[[nm]], digits), trim = TRUE),
      expression, perl = TRUE)
  }

  for (k in seq_along(latentNames)) {
    expression <- gsub(quoted(latentNames[k]), paste0('\\\\eta_{', k, '}'),
      expression, perl = TRUE)
  }
  expression <- gsub('\\bstate\\s*\\[\\s*([0-9]+)\\s*\\]', '\\\\eta_{\\1}',
    expression, perl = TRUE)

  for (k in seq_along(TDpredNames)) {
    expression <- gsub(quoted(TDpredNames[k]),
      paste0('\\\\text{', texPrep(TDpredNames[k]), '}'), expression, perl = TRUE)
  }
  expression <- gsub('tdpreds\\s*\\[\\s*rowi\\s*,\\s*([0-9]+)\\s*\\]',
    '\\\\text{TD}_{\\1}', expression, perl = TRUE)

  expression <- gsub('*', ' \\cdot ', expression, fixed = TRUE)

  # Function names set upright, so `log1p(...)` does not render as a product of
  # the letters l, o, g and the number 1p. Done name by name rather than with
  # one replacement pattern because the name needs tex-escaping on the way in,
  # and an underscore must be escaped there while the one in `\eta_{2}` must
  # not. `\text` and `\mathrm` are never caught: they are followed by a brace,
  # not a parenthesis.
  # The lookbehind keeps the `\cdot` just inserted above from being read as a
  # function name when the next token happens to be a parenthesis.
  notmacro <- '(?<![\\\\A-Za-z0-9_.])'
  functions <- unique(regmatches(expression,
    gregexpr(paste0(notmacro, '[A-Za-z][A-Za-z0-9_.]*(?=\\s*\\()'),
      expression, perl = TRUE))[[1L]])
  for (fn in functions) {
    replacement <- gsub('\\', '\\\\', paste0('\\mathrm{', texPrep(fn), '}'), fixed = TRUE)
    expression <- gsub(paste0(notmacro, '\\Q', fn, '\\E\\b(\\s*\\()'),
      paste0(replacement, '\\1'), expression, perl = TRUE)
  }
  trimws(gsub(' {2,}', ' ', expression))
}

ctModelLatexMathElement <- function(x){
  for(i in seq_along(x)){
    if(is.na(suppressWarnings(as.numeric(x[i]))) &&
        grepl('\\',x[i],fixed=TRUE) == FALSE) {
      x[i] <- paste0('\\text{',texPrep(x[i]),'}')
    }
  }
  x
}

bmatrix = function(x, digits=NULL,nottext=FALSE, ...) {
  if(!is.null(x)){
    if(!nottext){
      for(i in 1:length(x)){
        if(is.na(suppressWarnings(as.numeric(x[i]))) & #if x[i] cannot be numeric and 
            grepl('\\',x[i],fixed=TRUE) == FALSE) {
          x[i] = texPrep(x[i])
          x[i] = paste0('\\text{',x[i],'}')
        }
      }
    }
    x=as.matrix(x)
    
    out=c()
    for(i in 1:nrow(x)){
      for(j in 1:ncol(x)){
        out=c(out,x[i,j])
        if(j!=ncol(x)) out=c(out,paste0(' & '))
        if(j==ncol(x) & i!=nrow(x)) out=c(out,paste0('\\\\ \n'))
        if(j==ncol(x) & i==nrow(x)) out=c(out, '\n')
      }
    }
    
    out = paste0("\\begin{bmatrix}\n",
      paste0(out,collapse=''),
      "\\end{bmatrix}",collapse='')
  } else out=""
  return(out)
}

ctModelLatexDynamicsBlock <- function(ctmodel, showd, continuoustime, matrixnames=TRUE,
  splitDynamics=TRUE){
  drift <- paste0("\\underbrace{
        ",bmatrix(ctmodel$DRIFT),"
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{A}}",ifelse(!matrixnames,"}","_\\textrm{DRIFT}}")," \\underbrace{
        ",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
        \\big(t\\big)
      }_{\\vect{\\eta} (t",ifelse(continuoustime,"","-1"),")}	+ \\underbrace{
        ",bmatrix(ctmodel$CINT),"
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{b}}",ifelse(!matrixnames,"}","_\\textrm{CINT}}"),
    if(ctmodel$n.TDpred > 0) paste0( "+ \\underbrace{
        ",bmatrix(ctmodel$TDPREDEFFECT),"
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{M}}",ifelse(!matrixnames,"}","_\\textrm{TDPREDEFFECT}}"),"
      \\underbrace{
        ",bmatrix(matrix(ctmodel$TDpredNames))," 
      }_{\\vect{\\chi} (t)}"))
  diffusion <- paste0("\\underbrace{UcorSDtoChol\\left\\{
      ",bmatrix(ctmodel$DIFFUSION),"\\right\\}
    ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{G}}",ifelse(!matrixnames,"}","_\\textrm{DIFFUSION}}"),"
    \\underbrace{",showd,"
      ",bmatrix(matrix(paste0('W_{',1:ctmodel$n.latent,'}')),nottext=TRUE)," 
      (t)}_{",showd," \\vect{W}(t)}")
  
  if(splitDynamics) return(paste0("\\parbox{10em}{\\centering{Deterministic\\linebreak change:}}
  &\\underbrace{",showd,"
    ",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
    \\big(t\\big)}_{",showd," \\vect{\\eta} (t)}	=  \\left(
      ",drift,"
      \\right) ",ifelse(continuoustime,"\\mathrm{d}t","")," \\quad + \\nonumber \\\\ \\\\
    \\parbox{10em}{\\centering{Random\\linebreak change:}}
    & \\qquad \\qquad \\quad ",diffusion," \\\\ \\\\"))
  
  paste0("\\parbox{10em}{\\centering{Dynamics:}}
  &\\underbrace{",showd,"
    ",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
    \\big(t\\big)}_{",showd," \\vect{\\eta} (t)}	=  \\left(
      ",drift,"
      \\right) ",ifelse(continuoustime,"\\mathrm{d}t","")," \\quad + ",diffusion," \\\\ \\\\")
}

ctModelLatexMeasurementBlock <- function(ctmodel, matrixnames=TRUE,
  splitMeasurement=TRUE){
  manifesttype <- rep(0, ctmodel$n.manifest)
  if(!is.null(ctmodel$manifesttype)) manifesttype <- as.integer(ctmodel$manifesttype)
  binary <- manifesttype == 1
  ordinal <- manifesttype == 2
  count <- manifesttype == 3
  censored <- manifesttype == 4
  manifestNames <- ctmodel$manifestNames
  manifestIndex <- seq_along(manifestNames)
  binaryIndex <- which(binary)
  ordinalIndex <- which(ordinal)
  countIndex <- which(count)
  censoredIndex <- which(censored)
  
  if(!any(manifesttype != 0)){
    observation <- paste0("\\underbrace{
          ",bmatrix(ctmodel$LAMBDA)," 
        ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\Lambda}}",ifelse(!matrixnames,"}","_\\textrm{LAMBDA}}")," \\underbrace{
          ",bmatrix(matrix(ctmodel$latentNames))," 
          (t)}_{\\vect{\\eta}(t)} +
        \\underbrace{
          ",bmatrix(ctmodel$MANIFESTMEANS)," 
        ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\tau}}",ifelse(!matrixnames,"}","_\\textrm{MANIFESTMEANS}}"))
    noise <- paste0("\\underbrace{UcorSDtoChol \\left\\{
                ",bmatrix(ctmodel$MANIFESTVAR),"\\right\\}  
              ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\Theta}}",ifelse(!matrixnames,"}","_\\textrm{MANIFESTVAR}}"),"
              \\underbrace{
          ",bmatrix(matrix(paste0('\\epsilon_{',1:ctmodel$n.manifest,'}')))," 
          (t)}_{\\vect{\\epsilon}(t)}")
    
    if(splitMeasurement) return(paste0("\\parbox{10em}{\\centering{Observations:}}
&\\underbrace{
      ",bmatrix(matrix(ctmodel$manifestNames),nottext=FALSE),"  
      (t)}_{\\vect{Y}(t)} = 
        ",observation," + \\nonumber \\\\ \\\\
    \\parbox{10em}{\\centering{Observation\\linebreak noise:}}
    & \\qquad \\qquad \\quad  ",noise," \\\\ \\\\"))
    
    return(paste0("\\parbox{10em}{\\centering{Measurement:}}
&\\underbrace{
      ",bmatrix(matrix(ctmodel$manifestNames),nottext=FALSE),"  
      (t)}_{\\vect{Y}(t)} = 
        ",observation," + ",noise," \\\\ \\\\"))
  }
  
  linearPredictor <- paste0("\\underbrace{
          ",bmatrix(ctmodel$LAMBDA)," 
        ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\Lambda}}",ifelse(!matrixnames,"}","_\\textrm{LAMBDA}}")," \\underbrace{
          ",bmatrix(matrix(ctmodel$latentNames))," 
          (t)}_{\\vect{\\eta}(t)} +
        \\underbrace{
          ",bmatrix(ctmodel$MANIFESTMEANS)," 
        ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\tau}}",ifelse(!matrixnames,"}","_\\textrm{MANIFESTMEANS}}"))
  predictorLine <- paste0("\\parbox{10em}{\\centering{Linear\\linebreak predictor:}}
&\\underbrace{",bmatrix(matrix(paste0('\\nu_{',manifestIndex,'}(t)')),nottext=TRUE),"
      }_{\\vect{\\nu}(t)} =
        ",linearPredictor," \\\\ \\\\")
  yhat <- matrix(paste0('\\nu_{',manifestIndex,'}(t)'))
  yhat[binary] <- paste0('\\operatorname{logit}^{-1}\\left(\\nu_{',binaryIndex,'}(t)\\right)')
  # An ordinal row's predicted quantity is a set of cumulative probabilities
  # rather than a single number, so the cell carries the cumulative form and
  # the reader differences adjacent entries to get a category probability.
  yhat[ordinal] <- paste0('\\Pr\\left(Y_{',ordinalIndex,
    '}(t)\\leq k\\right)=\\operatorname{logit}^{-1}\\left(\\tau_{',
    ordinalIndex,',k}-\\nu_{',ordinalIndex,'}(t)\\right)')
  # A count's predicted quantity is its rate, and the link is the exponential:
  # the linear predictor is the log rate.
  yhat[count] <- paste0('\\exp\\left(\\nu_{',countIndex,'}(t)\\right)')
  # A censored variable's predicted quantity is its uncensored mean; what the
  # censoring changes is which values can be recorded, which the note explains
  # rather than the cell.
  yhat[censored] <- paste0('\\nu_{',censoredIndex,'}(t)')
  predictedLine <- paste0("\\parbox{10em}{\\centering{Predicted\\linebreak observations:}}
&\\underbrace{",bmatrix(matrix(paste0('\\hat{Y}_{',manifestIndex,'}(t)')),nottext=TRUE),"
      }_{\\widehat{\\vect{Y}}(t)} =
        ",bmatrix(yhat,nottext=TRUE)," \\\\ \\\\")
  observationLine <- paste0("\\parbox{10em}{\\centering{Observations:}}
&\\underbrace{
      ",bmatrix(matrix(manifestNames),nottext=FALSE),"  
      (t)}_{\\vect{Y}(t)} = 
        \\widehat{\\vect{Y}}(t) + \\vect{\\epsilon}(t) \\\\ \\\\")
  
  errorCov <- as.matrix(ctmodel$MANIFESTVAR)
  errorCovNumeric <- suppressWarnings(matrix(as.numeric(errorCov),
    nrow=nrow(errorCov), ncol=ncol(errorCov)))
  diagonalErrorCov <- !any(is.na(errorCovNumeric[row(errorCovNumeric) != col(errorCovNumeric)])) &&
    all(errorCovNumeric[row(errorCovNumeric) != col(errorCovNumeric)] == 0)
  errorCov <- ctModelLatexMathElement(errorCov)
  # Censored keeps the Gaussian form: it is the one non-Gaussian type with a
  # measurement error of its own, and that error is what MANIFESTVAR holds.
  gaussianerror <- !binary & !ordinal & !count
  diag(errorCov)[gaussianerror] <-
    paste0('\\left[',diag(errorCov)[gaussianerror],'\\right]^2')
  diag(errorCov)[binary] <- paste0('\\hat{Y}_{',binaryIndex,'}(t)\\left(1-\\hat{Y}_{',
    binaryIndex,'}(t)\\right)')
  # An ordinal observation has no Gaussian error term at all: the julia filter
  # integrates the categorical likelihood over the latent rather than matching
  # it to a normal, so there is nothing to put in this cell.
  diag(errorCov)[ordinal] <- '\\textrm{--}'
  # A Poisson's variance is its mean, so the cell is the predicted rate itself
  # rather than a free parameter -- there is nothing in MANIFESTVAR to show.
  diag(errorCov)[count] <- paste0('\\hat{Y}_{',countIndex,'}(t)')
  errorCovNote <- paste0(
    if(!diagonalErrorCov) paste0(" \\\\ 
&\\textrm{Note: off-diagonal entries in }\\vect{\\Theta}\\textrm{ are shown as specified; binary diagonal entries are conditional approximations.}") else "",
    if(any(ordinal)) paste0(" \\\\ 
&\\textrm{Note: ordinal observations carry no measurement error term -- the likelihood is integrated over }\\vect{\\eta}\\textrm{ directly.}") else "",
    if(any(count)) paste0(" \\\\ 
&\\textrm{Note: count observations are Poisson with a log link -- the variance shown is the predicted rate rather than a free parameter, and the rate is the exponential of the linear predictor in }\\vect{\\eta}\\textrm{ .}") else "",
    if(any(censored)) paste0(" \\\\ 
&\\textrm{Note: a censored observation is recorded as the value above clamped to its limits, so the equation gives the uncensored mean and observations pile up at censormin and censormax rather than passing them; the error term shown is real for this kind, unlike the others, and applies to }\\vect{\\eta}\\textrm{ .}") else "")
  errorLine <- paste0("\\parbox{10em}{\\centering{Observation\\linebreak error:}}
& \\qquad \\qquad \\quad \\vect{\\epsilon}(t) \\sim \\mathrm{N}\\left(\\mathbf{0},
      \\underbrace{",bmatrix(errorCov,nottext=TRUE),"
              ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\Theta}}",ifelse(!matrixnames,"}","_\\textrm{effective measurement covariance}}"),"\\right)",errorCovNote," \\\\ \\\\")
  
  if(splitMeasurement) return(paste0(
    predictorLine,predictedLine,observationLine,errorLine))
  
  paste0("\\parbox{10em}{\\centering{Measurement:}}
& \\underbrace{",bmatrix(matrix(paste0('\\nu_{',manifestIndex,'}(t)')),nottext=TRUE),"
      }_{\\vect{\\nu}(t)} = ",linearPredictor," \\\\ \\\\
& \\underbrace{",bmatrix(matrix(manifestNames),nottext=FALSE),"
      (t)}_{\\vect{Y}(t)} =
      \\widehat{\\vect{Y}}(t) + \\vect{\\epsilon}(t),\\quad
      \\widehat{\\vect{Y}}(t) = ",bmatrix(yhat,nottext=TRUE),
    " \\\\ \\\\
& \\vect{\\epsilon}(t) \\sim \\mathrm{N}\\left(\\mathbf{0},
      \\underbrace{",bmatrix(errorCov,nottext=TRUE),"
              ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\Theta}}",ifelse(!matrixnames,"}","_\\textrm{effective measurement covariance}}"),"\\right)",errorCovNote," \\\\ \\\\")
}


#' Generate and optionally compile latex equation of subject level ctsem model.
#'
#' @param x ctsem model object or ctStanFit object.
#' @param matrixnames Logical. If TRUE, includes ctsem matrix names such as DRIFT and DIFFUSION under the matrices.
#' @param digits Precision of decimals for numeric values.
#' @param linearise Logical. Show the linearised normal approximation for subject parameters and 
#' covariate effects, or the raw parameters?
#' @param textsize Standard latex text sizes -- 
#' tiny scriptsize footnotesize small normalsize large Large LARGE huge Huge. 
#' Useful if output overflows page. 
#' @param filename filename, without suffix, to output .tex and .pdf files too.
#' @param folder Character string specifying folder to save to, defaults to temporary directory, use "./" for working directory.
#' @param tex Save .tex file? Otherwise latex is simply returned within R as a string.
#' @param equationonly Logical. If TRUE, output is only the latex relevant to the equation, not a compileable document.
#' @param minimal if TRUE, outputs reduced form version displaying matrix dimensions and equation structure only.
#' @param splitDynamics Logical. If TRUE, split latent process dynamics across deterministic and random change lines.
#' If FALSE, show the full dynamics equation on one line.
#' @param splitMeasurement Logical. If TRUE, split measurement equations across observation and observation noise lines.
#' If FALSE, show the full measurement equation on one line.
#' @param compile Compile to .pdf? (Depends on \code{tex = TRUE}) 
#' @param open Open after compiling? (Depends on \code{compile = TRUE})
#' @param includeNote Include text describing matrix transformations and subject notation?
#' triangular matrices (which results in a covariance or Cholesky matrix) is shown -- 
#' the latter is a more direct representation of the model, while the former is often simpler to convey.
#' @param savepng Logical. If TRUE, renders the equation as a png file instead of a pdf, viewing the png in RStudio viewer when available.
#'
#' @return character string of latex code. Side effects include saving a .tex, .pdf, and displaying the pdf. 
#' @export
#' @importFrom tools texi2pdf
#'
#' @examples
#' ctmodel <- ctModel(type='ct', 
#' n.latent=2, n.manifest=1, 
#' manifestNames='sunspots', 
#' latentNames=c('ss_level', 'ss_velocity'),
#' LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
#' DRIFT=matrix(c(0, 1,   'a21', 'a22'), nrow=2, ncol=2, byrow=TRUE),
#' MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
#' CINT=matrix(c(0, 0), nrow=2, ncol=1),
#' DIFFUSION=matrix(c(
#'   0, 0,
#'   0, "diffusion"), ncol=2, nrow=2, byrow=TRUE))
#'   
#' l=ctModelLatex(ctmodel,compile=FALSE, open=FALSE)
#' cat(l)
ctModelLatex<- function(x,matrixnames=TRUE,digits=3,linearise=class(x) %in% 'ctStanFit',textsize='normalsize',folder=tempdir(),
  filename=paste0('ctsemTex',as.numeric(Sys.time())),tex=TRUE, equationonly=FALSE, compile=TRUE, open=TRUE, includeNote=TRUE,
  minimal=FALSE, splitDynamics=TRUE, splitMeasurement=TRUE, savepng=FALSE){
  #library(ctsem)
  dopopcov <- FALSE
  t0cov <- NULL
  t0means <- NULL
  latentPopNames <- NULL
  # Only a fit can have these; an unfitted model shows its expressions anyway.
  contextcells <- NULL
  
  # When savepng is TRUE, force compilation settings
  if(savepng) {
    equationonly <- FALSE
    compile <- TRUE
  }
  
  if('ctStanFit' %in% class(x)){
    ms=ctMatsetupFreePars(x$setup$matsetup)
    e=ctExtract(x)
    if(x$standata$ntipred > 0){
      if(linearise) timat <- round(ctCollapse(e$linearTIPREDEFFECT,1,mean),digits)
      if(!linearise) timat <- round(ctCollapse(e$TIPREDEFFECT,1,mean),digits)
      rownames(timat) <- ms$parname
      colnames(timat) <- x$ctstanmodel$TIpredNames
      timat <- timat[apply(x$standata$TIPREDEFFECTsetup,1,function(x) any(x!=0)),,drop=FALSE]
    } else timat <- diag(0,0)
    
    if(!is.null(e$rawpopcov)){ 
      
      if(!linearise) popcov <- round(ctCollapse(e$rawpopcov,1,mean),digits)
      if(linearise) {
        popcov <- stan_constrainsamples(x$stanmodel,x$standata,matrix(x$stanfit$rawest,nrow=1),
          cores=1,pcovn =1000,dokalman=FALSE,savesubjectmatrices = FALSE)$popcov
        popcov <- round(ctCollapse(e$popcov,1,mean),digits=digits)
        if(!is.null(e$pop_T0cov)) t0cov <- round(ctCollapse(e$pop_T0cov,1,mean),digits=digits)
        if(!is.null(e$pop_T0MEANS)) t0means <- round(ctCollapse(e$pop_T0MEANS,1,mean),digits=digits)
        if(!is.null(x$ctstanmodel$latentPopNames)) latentPopNames <- x$ctstanmodel$latentPopNames
        if(x$standata$intoverpop==1){
          t0index <- ms$indvarying[ms$param > 0 & ms$row <= x$standata$nlatent & ms$matrix %in% 1 & ms$indvarying > 0]
          popcov[t0index,t0index] <- round(ctCollapse(e$pop_T0cov,1,mean),digits=digits)[
            t0index,t0index] 
        }
        rownames(popcov) <- ms$parname[as.logical(ms$indvarying)]
        colnames(popcov) <- ms$parname[as.logical(ms$indvarying)]
      }
    } else popcov <- diag(0,0)
    
    if(!linearise) popmeans <- round(ctCollapse(e$rawpopmeans,1,mean),digits)[
      as.logical(ms$indvarying + ms$tipred),drop=FALSE]
    if(linearise) {
      popmeans <- round(ctCollapse(e$popmeans,1,mean),digits)[
        as.logical(ms$indvarying + ms$tipred),drop=FALSE]
      if(x$standata$intoverpop==1){
        popmeans[t0index]<- round(ctCollapse(e$pop_T0MEANS,1,mean),digits=digits)[
          t0index,1]
      }
    }
    
    # parmats <- summary(x,residualcov=FALSE,priorcheck=FALSE,digits=digits)
    # parmats <- data.frame(parmats$parmatrices,matrix=rownames(parmats$parmatrices))
    ctmodelmats <- listOfMatrices((x$ctstanmodelbase$pars))
    ctmodel <- x$ctstanmodelbase
    ####################################################################
    # Cells written as expressions over the state or the TD predictors keep
    # their structure; see .ctLatexRenderExpression above. Detected from the
    # *base* model, so the coordinates match ctmodelmats and so that an
    # intoverpop carrier state -- which is a parameter, not a reference to the
    # dynamics -- is not caught up in it.
    contextcells <- try(.ctFitContextDependentCells(x$ctstanmodelbase), silent = TRUE)
    if(inherits(contextcells,'try-error')) contextcells <- NULL
    estimates <- if(!is.null(contextcells) && nrow(contextcells)){
      try(.ctLatexRawEstimates(x, e), silent = TRUE)
    } else NULL
    if(inherits(estimates,'try-error')) estimates <- numeric()
    isContextCell <- function(mi,i,j) !is.null(contextcells) &&
      any(contextcells$matrix %in% mi & contextcells$row %in% i & contextcells$col %in% j)

    for(mi in names(ctmodelmats)){
      mimean <- ctCollapse(e[[paste0('pop_',mi)]],1,mean)
      for(i in 1:nrow(ctmodelmats[[mi]])){
        for(j in 1:ncol(ctmodelmats[[mi]])){
          if(isContextCell(mi,i,j)){
            ctmodelmats[[mi]][i,j] <- .ctLatexRenderExpression(
              .ctLatexCellExpression(x$ctstanmodelbase$pars,mi,i,j), estimates,
              ctmodel$latentNames, ctmodel$TDpredNames, digits)
          } else ctmodelmats[[mi]][i,j] <- round(mimean[i,j],digits)
        }
      }
    }
    ctmodel <- c(ctmodel,ctmodelmats)
    class(ctmodel) <- 'ctStanModel'
  } else ctmodel <- x
  
  if('ctStanModel' %in% class(ctmodel)) {
    
    if(!'ctStanFit' %in% class(x)){ #construct pop effects
      popcov <- ctModelBuildPopCov(ctmodel,linearise=linearise)
      if(ctmodel$n.TIpred > 0) timat <- ctModelBuildTIeffects(ctmodel) else timat <- diag(0,0)
      if(!linearise) timat[,] <- paste0('raw_',timat)
      popmeans<-paste0(ifelse(linearise,'','raw_'),unique(c(rownames(popcov),rownames(timat))))
      ctmodel<-T0VARredundancies(ctmodel)
      ctmodel <- c(ctmodel,listOfMatrices(ctmodel$pars)) 
    }
    
    dopopcov <- as.logical(nrow(popcov))
    doti <- as.logical(nrow(timat))
    
    if(doti){
      # both <- rownames(timat) %in% rownames(popcov)
      pars <- unique(c(rownames(popcov),rownames(timat)))
      newpopcov <- matrix(0,length(pars),length(pars),dimnames=list(pars,pars))
      newtimat <- matrix(0,length(pars),ncol(timat),dimnames=list(pars,colnames(timat)))
      newtimat[na.omit(match(rownames(timat),rownames(newtimat))),] <- timat
      newpopcov[na.omit(match(rownames(popcov),rownames(newpopcov))), 
        na.omit(match(rownames(popcov),rownames(newpopcov)))] <- popcov
      popcov <- newpopcov
      timat <- newtimat
    }
    
    if(dopopcov || doti){
      t0covIsTotal <- !is.null(t0cov)
      if(exists('t0cov') && is.null(t0cov) && exists('ctmodel')) {
        t0cov <- ctModelLatexT0DisplayCov(
          as.matrix(ctmodel$T0VAR)[seq_len(ctmodel$n.latent), seq_len(ctmodel$n.latent), drop=FALSE],
          digits=digits)
      }
      if(exists('t0means') && !is.null(t0means)) {
        t0means <- as.matrix(t0means)[seq_len(ctmodel$n.latent),1]
      }
      t0aug <- ctModelLatexAugmentT0(popmeans=popmeans, popcov=popcov,
        timat=timat, ctm=ctmodel, digits=digits, t0means=t0means,
        t0cov=t0cov, latentPopNames=latentPopNames, t0covIsTotal=t0covIsTotal)
      popmeans <- t0aug$popmeans
      popcov <- t0aug$popcov
      timat <- t0aug$timat
      dopopcov <- as.logical(nrow(popcov))
      doti <- as.logical(nrow(timat))
    }
    
    dopop <- doti||dopopcov
    
    ### this section replaced t0var fixed values with params, seemed broken...
    # if(!'ctStanFit' %in% class(x)){ #if a model
    #   
    #   t0index <- unique(ctmodel$pars$row[ctmodel$pars$matrix %in% 'T0MEANS' & 
    #       ctmodel$pars$indvarying & is.na(ctmodel$pars$value)]) #which t0means are indvarying
    #   if(length(t0index)){
    #     t0varpopcov <- matrix(
    #       paste0('Pcorsqrt_',t0index,'_',rep(t0index,each=length(t0index))),
    #       nrow=length(t0index),ncol=length(t0index))
    #     t0varpopcov[upper.tri(t0varpopcov)] <- 0
    #     # ms$indvarying[ms$param > 0 & ms$row <= x$standata$nlatent & ms$matrix %in% 1 & ms$indvarying > 0]
    #     ctmodel$T0VAR[t0index,t0index] <- t0varpopcov
    #   }
    # }
    
    continuoustime <- ctmodel$continuoustime
  } else {
    dopop <- FALSE
    if(! 'ctsemInit' %in% class(ctmodel)) stop('not a ctsem model!')
    continuoustime <- TRUE
  }
  
  if(equationonly) compile <- FALSE
  
  
  
  
  
  
  W <- diag(1,1)
  if(continuoustime) diag(W) <- 't-u'
  
  #out = 'Hello' 
  
  
  out <- ifelse(equationonly,"",paste0("
\\documentclass[a4paper]{article}
\\usepackage{geometry}
\\geometry{paperwidth=\\maxdimen,paperheight=\\maxdimen,margin=1cm}

\\usepackage[fleqn]{amsmath} %for multiple line equations
\\usepackage[active,tightpage,displaymath]{preview}
\\usepackage{bm}
\\newcommand{\\vect}[1]{\\boldsymbol{\\mathbf{#1}}}

\\begin{document}
\\pagenumbering{gobble}
\\begin{",textsize,"}
"))
  
  if (minimal){
    dict = list('A' = 'DRIFT','b'='CINT','M'='TDPREDEFFECT','G'='DIFFUSION','tau'='MANIFESTMEANS')
    
    for (name in names(dict)) {
      chmat = ctmodel[[dict[[name]]]]
      if (!is.numeric(chmat)){
        dict[[name]] = TRUE
      } else {
        #print('Recognized as numeric')
        if (max(abs(chmat)) < 1e-3) {
          dict[[name]] = FALSE
        } else dict[[name]] = TRUE
      }
    }
    
    nu = ctmodel$n.latent
    c = ctmodel$n.manifest
    l = ctmodel$n.TDpred
    
    tablestring <- '\\begin{center}
\\begin{tabular}'#{c|c|c|c} missing
    equationstring <- '\\begin{align*} \n'
    
    tabledim = '{c'
    tablecont1 = '$\\eta(t)$'
    tablecont2 = paste0('$',nu,'$')
    noisestring = ''
    equationcont = 'd\\eta(t) &= '
    
    if (!dict[['A']] & !dict[['b']] & !dict[['M']] ){ #if all of these are not in the equation, we only have noise.
      if (!dict[['G']]) equationcont = paste0(equationcont,'\\mathbf{0}') 
    } else {
      equationcont = paste0(equationcont, '\\left(')
      if (dict[['A']]){
        tabledim = paste0(tabledim,'|c')
        tablecont1 = paste0(tablecont1,'& $\\mathbf{A}$')
        tablecont2 = paste0(tablecont2,'& $',nu,'\\times', nu,'$')
        equationcont = paste0(equationcont,'\\mathbf{A} \\eta(t)')
      }
      if (dict[['b']]){
        tabledim = paste0(tabledim,'|c')
        tablecont1 = paste0(tablecont1,'& $\\mathbf{b}$')
        tablecont2 = paste0(tablecont2,'& $',nu,'$')
        if(dict[['A']]) equationcont = paste0(equationcont,'+')
        equationcont = paste0(equationcont,'\\mathbf{b}')
      }
      if (dict[['M']] && !l==0){
        tabledim = paste0(tabledim,'|c|c')
        tablecont1 = paste0(tablecont1,'& $\\mathbf{M}$ & $\\chi(t)$')
        tablecont2 = paste0(tablecont2,'& $',nu,'\\times', l,'$','& $',l,'$')
        if(dict[['A']]||dict[['b']]) equationcont = paste0(equationcont,'+')
        equationcont = paste0(equationcont,'\\mathbf{M} \\chi(t)')
      }
    }
    
    if (dict[['A']] || dict[['b']] || dict[['M']] ) equationcont = paste0(equationcont,'\\right) dt')
    
    
    if (dict[['G']]){
      tabledim = paste0(tabledim,'|c|c')
      tablecont1 = paste0(tablecont1,'& $\\mathbf{G}$ & $d\\mathbf{W}(t) $')
      tablecont2 = paste0(tablecont2,'& $',nu,'\\times', nu,'$','& $',nu,'$')
      noisestring = paste0(noisestring,'\\mathbf{W}(t+\\Delta t)-\\mathbf{W}(t) &\\sim N(\\mathbf{0},\\mathrm{diag}(\\Delta t)) \\\\','\n' )
      if (dict[['A']] || dict[['b']] || dict[['M']] ) equationcont = paste0(equationcont,'+')
      equationcont = paste0(equationcont,'\\mathbf{G} d\\mathbf{W}(t)')
    }
    
    manifesttype <- rep(0, ctmodel$n.manifest)
    if(!is.null(ctmodel$manifesttype)) manifesttype <- as.integer(ctmodel$manifesttype)
    nbinary <- sum(manifesttype == 1)
    ncount <- sum(manifesttype == 3)
    
    if(!any(manifesttype != 0)){
      tabledim = tabledim = paste0(tabledim,'|c|c')
      tablecont1 = paste0(tablecont1,'& $\\mathbf{y}$ & $\\Lambda $')
      tablecont2 = paste0(tablecont2,'& $',c,'$','& $',nu,'\\times', c,'$')
      equationcont = paste0(equationcont,'\\\\','\n', '\\mathbf{y}(t) &= \\Lambda \\eta(t)')
      
      if (dict[['tau']]){
        tabledim = paste0(tabledim,'|c')
        tablecont1 = paste0(tablecont1,'& $\\tau$')
        tablecont2 = paste0(tablecont2,'& $',c,'$')
        equationcont = paste0(equationcont,'+ \\tau')
      }
      
      tabledim = tabledim = paste0(tabledim,'|c|c}')
      tablecont1 = paste0(tablecont1,'& $\\epsilon(t)$ & $\\Theta$ \\\\','\n', '\\hline')
      tablecont2 = paste0(tablecont2,'& $',c,'$','& $',c,'\\times', c,'$')
      noisestring = paste0(noisestring,'\\epsilon(t) &\\sim N(\\mathbf{0},\\Theta) \\\\','\n' )
      equationcont = paste0(equationcont,'+ \\epsilon(t)')
    } else {
      tabledim = tabledim = paste0(tabledim,'|c|c|c')
      tablecont1 = paste0(tablecont1,'& $\\nu(t)$ & $\\hat{Y}(t)$ & $\\Lambda $')
      tablecont2 = paste0(tablecont2,'& $',c,'$','& $',c,'$','& $',nu,'\\times', c,'$')
      equationcont = paste0(equationcont,'\\\\','\n',
        '\\nu(t) &= \\Lambda \\eta(t)')
      
      if (dict[['tau']]){
        tabledim = paste0(tabledim,'|c')
        tablecont1 = paste0(tablecont1,'& $\\tau$')
        tablecont2 = paste0(tablecont2,'& $',c,'$')
        equationcont = paste0(equationcont,'+ \\tau')
      }
      
      binaryIndex <- which(manifesttype == 1)
      ordinalIndex <- which(manifesttype == 2)
      countIndex <- which(manifesttype == 3)
      yhat <- paste0('\\nu_{',seq_len(c),'}(t)')
      yhat[binaryIndex] <- paste0('\\operatorname{logit}^{-1}(\\nu_{',binaryIndex,'}(t))')
      yhat[ordinalIndex] <- paste0('\\operatorname{logit}^{-1}(\\tau_{',
        ordinalIndex,',k}-\\nu_{',ordinalIndex,'}(t))')
      yhat[countIndex] <- paste0('\\exp(\\nu_{',countIndex,'}(t))')
      yhatString <- paste0('\\begin{bmatrix}',paste(yhat,collapse=' \\\\ '),'\\end{bmatrix}')
      equationcont = paste0(equationcont,'\\\\','\n',
        '\\hat{Y}(t) &= ',yhatString)
      equationcont = paste0(equationcont,'\\\\','\n',
        '\\mathbf{y}(t) &= \\hat{Y}(t) + \\epsilon(t)')
      
      tabledim = paste0(tabledim,'|c|c')
      tablecont1 = paste0(tablecont1,'& $\\epsilon(t)$ & $\\Theta$')
      tablecont2 = paste0(tablecont2,'& $',c,'$','& $',c,'\\times', c,'$')
      noisestring = paste0(noisestring,'\\epsilon(t) &\\sim N(\\mathbf{0},\\Theta) \\\\','\n' )
      if(nbinary > 0){
        noisestring = paste0(noisestring,'\\Theta_{',binaryIndex,',',binaryIndex,'} &= \\hat{Y}_{',
          binaryIndex,'}(t)(1-\\hat{Y}_{',binaryIndex,'}(t)) \\\\','\n',
          collapse='')
      }
      if(ncount > 0){
        noisestring = paste0(noisestring,'\\Theta_{',countIndex,',',countIndex,
          '} &= \\hat{Y}_{',countIndex,'}(t) \\\\','\n', collapse='')
      }
      
      tabledim = paste0(tabledim,'}')
      tablecont1 = paste0(tablecont1,' \\\\','\n', '\\hline')
    }
    
    
    tablestring = paste0(tablestring,tabledim,'\n',tablecont1,'\n',tablecont2,'\n','\\end{tabular}','\n','\\end{center}')
    
    equationstring = paste0(equationstring,noisestring,'\\\\',equationcont,'\n','\\end{align*}')
    
    out= paste0(out,tablestring,'\n',equationstring)
    
    
  } else { #end minimal
    showd <- ifelse(continuoustime,"\\mathrm{d}","") #for continuous or discrete system
    
    # if(covMatrices){
    #   if('ctStanFit' %in% class(m)){
    #     cp <- ctSummaryMatrices(m)
    #     ctmodel$T0VAR <- cp$T0COV
    #     ctmodel$DIFFUSION <- cp$DIFFUSIONcov
    #     } else {
    #       ctmodel$T0VAR[upper.tri(ctmodel$T0VAR)] <- t(ctmodel$T0VAR)[upper.tri(ctmodel$T0VAR)]
    #       ctmodel$DIFFUSION[upper.tri(ctmodel$DIFFUSION)] <- t(ctmodel$DIFFUSION)[upper.tri(ctmodel$DIFFUSION)]
    #     }
    # }
    
    
    
    
    initialStateLatex <- if(!dopop) paste0("\\parbox{10em}{\\centering{Initial\\linebreak latent\\linebreak state:}}
  &\\underbrace{",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
    \\big(t_0\\big)}_{\\vect{\\eta} (t_0)}	\\sim \\mathrm{N} \\left(
              \\underbrace{
        ",bmatrix(ctmodel$T0MEANS),"
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{}}",ifelse(!matrixnames,"}","_\\textrm{T0MEANS}}"),",
      \\underbrace{UcorSDtoCov \\left\\{","
        ",bmatrix(ctmodel$T0VAR),"\\right\\}"," 
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{Q^{*}}_{t0}}",ifelse(!matrixnames,"}","_\\textrm{T0VAR}}"),"
      \\right) \\\\
") else ""
    
    out <- paste0(out, "
 \\setcounter{MaxMatrixCols}{200}
 \\begin{flalign*}
  &\\begin{aligned}
  ",if(dopop) paste0("\\parbox{10em}{\\centering{Initial\\linebreak and subject\\linebreak parameter\\linebreak distribution:}}
             &\\underbrace{",bmatrix(matrix(paste0('\\text{',
               texPrep(colnames(popcov)),'}_i')),nottext=TRUE)," 
            }_{\\vect{\\phi}(i)} ",ifelse(linearise,"\\approx","\\sim"),
    ifelse(linearise,"","\\textrm{tform}\\left\\{"),
    " \\mathrm{N} \\left(
              ",bmatrix(popmeans),", ", bmatrix(popcov)," \\right) ",
    if(doti) paste0(" + \\underbrace{",bmatrix(timat),"}_{\\vect{",ifelse(linearise,"\\hat",""),"\\beta}}","
  \\underbrace{
    ",bmatrix(matrix(colnames(timat))),"}_{\\vect{z}}"),
    ifelse(linearise,"","\\right\\}")," \\\\"), initialStateLatex,
      ctModelLatexDynamicsBlock(ctmodel=ctmodel, showd=showd,
        continuoustime=continuoustime, matrixnames=matrixnames,
        splitDynamics=splitDynamics),
      ctModelLatexMeasurementBlock(ctmodel=ctmodel, matrixnames=matrixnames,
        splitMeasurement=splitMeasurement),
      "\\parbox{10em}{\\centering{System noise\\linebreak distribution per time step:}}
          &",ifelse(continuoustime,'\\Delta ',''),"\\big[W_{j \\in [1,",ctmodel$n.latent,"]}\\big](t",
      ifelse(continuoustime,'-u',''),")   \\sim  \\mathrm{N}(0,",W,") \\\\ \\\\
      \\end{aligned} \\\\",
      if(includeNote) paste0("&\\textrm{Note: } UcorSDtoChol\\textrm{ converts lower tri matrix of standard deviations and unconstrained correlations to Cholesky factor,} \\\\
&UcorSDtoCov =\\textrm{ transposed cross product of UcorSDtoChol, to give covariance, See Driver \\& Voelkle (2018) p11.} \\\\",
        if(dopop) paste0("&\\textrm{Individual specific notation (subscript i) only shown for subject parameter distribution -- pop. means shown elsewhere.} \\\\
",if(linearise) "&\\textrm{Linearised approximation of subject parameter distribution shown.} \\\\"),
        # Only for a model that has such cells: for every other model this line
        # would be a caveat about something the reader cannot see.
        if(isTRUE(nrow(contextcells) > 0)) "&\\textrm{Cells depending on the latent state or a time dependent predictor keep their expression, with estimates substituted for parameter labels.} \\\\"),
      "\\end{flalign*}
      ")
  }
  
  if(!equationonly) out <- paste0(out, 
    "  \\end{",textsize,"}
\\end{document}")
  
  
  if(tex) {
    oldwd <- getwd()
    setwd(dir = folder)
    on.exit(setwd(oldwd))
    write(x = out,file = paste0(filename,'.tex'))
    if(compile){
      hastex <- !Sys.which('pdflatex') %in% ''
      a=try(tools::texi2pdf(file=paste0(filename,'.tex'),quiet=FALSE, clean=TRUE))
      if('try-error' %in% class(a)) {
        
        if(!grepl('SunOS',Sys.info()['sysname']) && requireNamespace('tinytex',quietly=TRUE)){
          a=try(tinytex::pdflatex(file=paste0(filename,'.tex'), clean=TRUE))
          if('try-error' %in% class(a)) 'Error - Perhaps tinytex needs to be installed via: tinytex::install_tinytex()' 
        } else {
          open <- FALSE
          message('Tex compiler not found -- you could install the tinytex package using:\ninstall.packages("tinytex")\ntinytex::install_tinytex()')
        }
      }
      
      # Handle PNG rendering
      if(!'try-error' %in% class(a)) {
        if(requireNamespace('pdftools', quietly=TRUE)) {
          bitmap <- pdftools::pdf_render_page(paste0(filename,'.pdf'), dpi=300)
          
          # Display PNG in RStudio viewer if available
          if(interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
            # Write bitmap to temp PNG file for viewer
            temp_png <- tempfile(tmpdir = tempdir(), fileext = ".png")
            png::writePNG(bitmap, temp_png)
            # Use relative path from HTML file
            temp_html <- tempfile(tmpdir = dirname(temp_png), fileext = ".html")
            png_name <- basename(temp_png)
            writeLines(paste0('<img src="', png_name, '" style="max-width:100%; height:auto;">'), temp_html)
            rstudioapi::viewer(temp_html)
          } else if(interactive()) {
            # Fallback to plot device - convert bitmap to raster first
            temp_png <- tempfile(fileext = ".png")
            png::writePNG(bitmap, temp_png)
            img <- png::readPNG(temp_png)
            old_par <- par(mar=c(0,0,0,0))  # Save old settings
            plot(as.raster(img), asp=1)
            par(old_par)  # Restore original settings
            unlink(temp_png)
          }
        } else {
          message('pdftools package required for viewing output within R. Install with: install.packages("pdftools")')
        }
        
        if(savepng){
          if(requireNamespace('pdftools', quietly=TRUE)) {
            # Save the bitmap as a PNG file
            png_filename <- paste0(filename, '.png')
            png::writePNG(bitmap, png_filename)
            message(paste("PNG saved to", png_filename))
          } else {
            stop("pdftools package required for render_png=TRUE. Install with: install.packages('pdftools')")
          }
        }
        
        if(interactive() && open) try(openPDF(paste0(filename,'.pdf')))
      } #end if succcessful compile
    } #end if compile
  }# end if tex
  return(invisible(out))
}

