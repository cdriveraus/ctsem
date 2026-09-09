ctModelUnlist<-function(ctmodelobj,
  matnames=c('T0MEANS','LAMBDA','DRIFT','DIFFUSION','MANIFESTVAR','MANIFESTMEANS', 'CINT', 'TDPREDEFFECT', 'T0VAR','PARS','THRESHOLDS')){
  out<-data.frame(matrix=as.character(NA), row=as.integer(NA), col=as.integer(NA), param=as.character(NA), value=as.numeric(NA),
    stringsAsFactors =FALSE) 
  out[1:sum(sapply(ctmodelobj[names(ctmodelobj) %in% matnames],length)),]=out

  rowcount <- 0
  for(obji in matnames){
    if(!is.null(dim(ctmodelobj[[obji]])) && !is.null(ctmodelobj[[obji]]) && !is.na(ctmodelobj[[obji]][1])){
      for(rowi in 1:nrow(ctmodelobj[[obji]])){
        for(coli in 1:ncol(ctmodelobj[[obji]])){
          rowcount <- rowcount + 1
          out[rowcount,]<-list(obji,rowi,coli,
            ifelse(is.na(suppressWarnings(as.numeric(ctmodelobj[[obji]][rowi,coli]))), #ifelse element is character string
              ctmodelobj[[obji]][rowi,coli],
              NA),
            ifelse(!is.na(suppressWarnings(as.numeric(ctmodelobj[[obji]][rowi,coli]))), #ifelse element is numeric
              as.numeric(ctmodelobj[[obji]][rowi,coli]),
              NA)
          )
        }
      }
    }
  }
  return(out)
}

.ctModelDefaultFreePar <- function(matrix, row, col, continuoustime){
  transform <- 0
  multiplier <- 1
  meanscale <- 1
  offset <- 0
  inneroffset <- 0

  if(matrix %in% c('T0MEANS','MANIFESTMEANS','TDPREDEFFECT','CINT')) {
    meanscale <- 10
  }
  if(matrix %in% c('LAMBDA')) {
    offset <- 0.5
    meanscale <- 5
  }
  # Column 1 is the first threshold, free on the real line like any intercept.
  # Later columns are gaps to the previous threshold and must be positive, so
  # they take the same positive transform the variances use. At a raw value of
  # zero a gap is log(2)*2 = 1.39, which is a sane spacing on a logit scale.
  if(matrix %in% c('THRESHOLDS')) {
    if(col == 1) meanscale <- 10
    if(col > 1) {
      transform <- 1
      meanscale <- 2
      multiplier <- 2
    }
  }

  if(matrix %in% c('DIFFUSION','MANIFESTVAR', 'T0VAR')) {
    if(row != col){
      transform <- 3
      multiplier <- 2
      offset <- -1
      meanscale <- 1
    }
    if(row == col){
      transform <- 1
      meanscale <- 2
      multiplier <- 5
      offset <- 1e-10
      if(matrix %in% c('DIFFUSION')) multiplier <- 10
    }
  }
  if(matrix %in% c('DRIFT')) {
    if(row == col){
      if(continuoustime==TRUE) {
        transform <- 1
        meanscale <- -2
        multiplier <- -2
        offset <- -1e-6
      }
      if(continuoustime==FALSE) {
        transform <- 3
        meanscale <- 2
        offset <- 0
      }
    }
    if(row != col){
      transform <- 0
      meanscale <- 1
    }
  }

  list(
    transform = Simplify(tform(parin = 'param',
      transform = as.integer(transform),
      multiplier = multiplier,
      meanscale = meanscale,
      offset = offset,
      inneroffset = inneroffset,
      singletext = TRUE)),
    indvarying = matrix %in% c('T0MEANS','MANIFESTMEANS','CINT'),
    sdscale = 1
  )
}

.ctModelMatrixValue <- function(value){
  if(length(value) != 1) stop('Matrix elements must have length 1')
  if(is.factor(value)) value <- as.character(value)
  if(is.na(value)) return(list(param=NA_character_, value=NA_real_))

  numericvalue <- suppressWarnings(as.numeric(value))
  if(!is.na(numericvalue)) {
    return(list(param=NA_character_, value=numericvalue))
  }
  list(param=as.character(value), value=NA_real_)
}

.ctModelMatricesPlaceholder <- function(){
  'pars-backed matrix view -- use model$matrices or ctModelMatrices(model)'
}

.ctModelUpdateParsFromMatrices <- function(ctm, matrices){
  if(!'ctStanModel' %in% class(ctm)) stop('x must be a ctStanModel object')
  if(!is.list(matrices) || is.null(names(matrices))) {
    stop('matrices must be a named list of matrices')
  }

  # RAWPOPVAR is not a system matrix and has no rows in `pars`; it is stored on
  # the model and taken out here before the loop below, which requires every
  # matrix it sees to be present in `pars`.
  if('RAWPOPVAR' %in% names(matrices)){
    ctm <- .ctModelRawPopVarAssign(ctm, matrices[['RAWPOPVAR']])
    matrices[['RAWPOPVAR']] <- NULL
  }

  pars <- ctm[['pars']]
  tieffects <- colnames(pars)[grep('_effect', colnames(pars), fixed=TRUE)]

  if(!length(matrices)) return(ctm)

  for(matrixname in names(matrices)){
    mat <- matrices[[matrixname]]
    if(!is.matrix(mat)) stop(matrixname, ' must be a matrix')

    matrixrows <- pars$matrix %in% matrixname
    if(!any(matrixrows)) stop(matrixname, ' is not present in x$pars')

    expecteddim <- c(max(pars$row[matrixrows]), max(pars$col[matrixrows]))
    if(!identical(dim(mat), expecteddim)) {
      stop(matrixname, ' must have dimensions ', paste(expecteddim, collapse=' x '))
    }

    for(rowi in seq_len(nrow(mat))){
      for(coli in seq_len(ncol(mat))){
        parrow <- which(matrixrows & pars$row == rowi & pars$col == coli)
        if(length(parrow) != 1) {
          stop('Could not match one pars row for ', matrixname, '[', rowi, ',', coli, ']')
        }

        wasfixed <- !is.na(pars$value[parrow])
        parsed <- .ctModelMatrixValue(mat[rowi,coli])
        # A cell assigned here can carry the same fields a cell given to
        # ctModel() can, written either way (R/ctParSpec.R). Without this the
        # separators or `key=value` text would be stored as the parameter name.
        spec <- .ctCellSpecFields(parsed$param)
        parsed$param <- spec$param
        pars$param[parrow] <- parsed$param
        pars$value[parrow] <- parsed$value

        if(!is.na(parsed$value)){
          pars$transform[parrow] <- NA
          pars$indvarying[parrow] <- FALSE
          pars$sdscale[parrow] <- NA_real_
          if(length(tieffects) > 0) pars[parrow,tieffects] <- 'FALSE'
        } else if(!is.na(parsed$param)){
          defaults <- .ctModelDefaultFreePar(
            matrix=matrixname,
            row=rowi,
            col=coli,
            continuoustime=ctm[['continuoustime']])
          if(wasfixed || is.na(pars$transform[parrow])) pars$transform[parrow] <- defaults$transform
          if(wasfixed || is.na(pars$sdscale[parrow])) pars$sdscale[parrow] <- defaults$sdscale
          pars$indvarying[parrow] <- as.logical(pars$indvarying[parrow])
          if(wasfixed || is.na(pars$indvarying[parrow])) pars$indvarying[parrow] <- defaults$indvarying

          # Anything the cell states explicitly beats the default for that cell.
          if(!is.na(spec$transform)) pars$transform[parrow] <- spec$transform
          if(!is.na(spec$indvarying)) pars$indvarying[parrow] <- spec$indvarying
          if(!is.na(spec$sdscale)) pars$sdscale[parrow] <- spec$sdscale
          if(!is.null(spec$tipreds)){
            wanted <- paste0(spec$tipreds, '_effect')
            unknown <- setdiff(wanted, tieffects)
            if(length(unknown)) {
              stop(paste(sub('_effect$', '', unknown), collapse=', '),
                ' is not a time independent predictor of this model')
            }
            pars[parrow,tieffects] <- 'FALSE'
            pars[parrow,wanted] <- 'TRUE'
          }
        }
      }
    }
  }

  ctm[['pars']] <- pars
  if(is.null(ctm[['matrices']])) ctm[['matrices']] <- .ctModelMatricesPlaceholder()
  ctm
}

#' Matrix view for ctStanModel objects
#'
#' Access or replace the matrix representation of a \code{ctStanModel} while
#' keeping the \code{pars} data frame as the canonical model specification.
#'
#' @param x A \code{ctStanModel} object.
#' @param value A named list of matrices, typically copied from
#' \code{ctModelMatrices(x)} or \code{x$matrices}.
#'
#' @details
#' \code{ctModelMatrices(x)} reconstructs matrices from \code{x$pars}. The
#' replacement form updates matching \code{matrix}, \code{row}, and \code{col}
#' entries in \code{x$pars}. Numeric matrix entries become fixed values;
#' non-numeric entries become parameter labels. Metadata such as transforms and
#' individual variation settings is preserved when possible and reset to
#' fit-ready defaults when fixed/free status changes.
#'
#' The same view is available as \code{x$matrices}. Direct replacements such as
#' \code{x$matrices$DRIFT[1, 2] <- "cross"} update \code{x$pars}; detached copies
#' must be assigned back with \code{x$matrices <- mats}. A placeholder
#' \code{matrices} element is stored in the object so the view is visible in
#' \code{names(x)} and printed objects, but the matrix values are always derived
#' from \code{x$pars}.
#'
#' @return \code{ctModelMatrices()} returns a named list of matrices. The
#' replacement form returns the updated \code{ctStanModel}.
#' @export
ctModelMatrices <- function(x){
  if(!'ctStanModel' %in% class(x)) stop('x must be a ctStanModel object')
  out <- listOfMatrices(x[['pars']])
  # Rebuilt on read rather than trusted: `indvarying` is routinely set directly
  # after the model is built, and a RAWPOPVAR describing a different set of random
  # effects than the model currently has would be worse than none.
  synced <- .ctModelRawPopVarSync(x)
  if(!is.null(synced[['RAWPOPVAR']])) out$RAWPOPVAR <- synced[['RAWPOPVAR']]
  out
}

#' @rdname ctModelMatrices
#' @export
`ctModelMatrices<-` <- function(x, value){
  .ctModelUpdateParsFromMatrices(x, value)
}

#' @export
`$.ctStanModel` <- function(x, name){
  nameindex <- pmatch(name, names(x), duplicates.ok=FALSE)
  if(!is.na(nameindex) && identical(names(x)[nameindex], 'matrices')) return(ctModelMatrices(x))
  x[[name, exact=FALSE]]
}

#' @export
`$<-.ctStanModel` <- function(x, name, value){
  if(identical(name, 'matrices')) return(.ctModelUpdateParsFromMatrices(x, value))
  x[[name]] <- value
  x
}

#' Convert an old OpenMx-style ctsem model to the modern ctsem model format.
#'
#' Convert an old OpenMx-style ctsem model to the modern ctsem model format.
#' \code{ctStanModel} is maintained as a backward-compatible alias.
#'
#' @param ctmodelobj ctsem model object created by \code{\link{ctModel}} with
#' \code{type='omx'}.
#' @param type Either \code{'ct'} for continuous time, or \code{'dt'} for
#' discrete time.
#' @param tipredDefault Logical. TRUE sets any parameters with unspecified time independent 
#' predictor effects to have effects estimated, FALSE fixes the effect to zero unless individually specified.
#'
#' @return List object of class ctStanModel in the modern ctsem model format,
#' with a \code{pars} data frame used by \code{\link{ctFit}}. Random effects are
#' specified for any intercept type parameters (T0MEANS, MANIFESTMEANS, and or
#' CINT), and time independent predictor effects are specified for all
#' parameters. Adjust these after initial specification by directly editing the
#' \code{pars} subobject, so \code{model$pars}, or by editing the pars-backed
#' matrix view, so \code{model$matrices}.
#' @details
#' \code{ctModelConvertOMX()} is primarily a compatibility bridge for older
#' workflows that first build a matrix-list model with \code{ctModel(type='omx')}.
#' New Stan-based models can usually be created directly with
#' \code{ctModel(type='ct')} or \code{ctModel(type='dt')}, which already return
#' the modern format.
#' @importFrom rstantools rstan_config
#' @aliases ctStanModel
#' @export
#'
#' @examples
#' model <- ctModel(type='omx', Tpoints=50,
#' n.latent=2, n.manifest=1, 
#' manifestNames='sunspots', 
#' latentNames=c('ss_level', 'ss_velocity'),
#' LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
#' DRIFT=matrix(c(0, 1,   'a21', 'a22'), nrow=2, ncol=2, byrow=TRUE),
#' MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
#' # MANIFESTVAR=matrix(0, nrow=1, ncol=1),
#' CINT=matrix(c(0, 0), nrow=2, ncol=1),
#' DIFFUSION=matrix(c(
#'   0, 0,
#'   0, "diffusion"), ncol=2, nrow=2, byrow=TRUE))
#' 
#' modernmodel=ctModelConvertOMX(model)
#' 
#' 
ctModelConvertOMX<-function(ctmodelobj, type='ct',tipredDefault=TRUE){
  if(FALSE) rstanconfig() #placeholder to use rstantools
  if(type=='stanct' | type=='standt'){
    warning('type should now be specified as simply ct or dt, without the stan prefix')
    type <- gsub('stan','',type)
  }
  
  if(type=='ct') continuoustime<-TRUE
  if(type=='dt') continuoustime<-FALSE
  
  if(!type %in% c('ct','dt')) stop('type must be either ct or dt!')
  
  ctm <- ctmodelobj
  
  if(!is.null(ctm$timeVarying)) stop('Time varying parameters not allowed for ctsem Stan model at present! Correct ctModel spec')
  
  
  
  
  #read in ctmodel values
  n.latent<-ctm$n.latent
  n.manifest<-ctm$n.manifest
  Tpoints<-ctm$Tpoints
  n.TDpred<-ctm$n.TDpred
  n.TIpred<-ctm$n.TIpred
  
  manifestNames<-ctm$manifestNames
  latentNames<-ctm$latentNames
  TDpredNames<-ctm$TDpredNames
  TIpredNames<-ctm$TIpredNames
  
  ctspec<-ctModelUnlist(ctm)

  # A cell may state its fields with `|` separators ('mm||TRUE|0.5') or by name
  # ('mm, indvarying=TRUE, sdscale=0.5'). Normalise the named form to the `|`
  # form here so the parser below has exactly one thing to read; see
  # .ctCellSpecToPipe in R/ctParSpec.R for how the two are told apart.
  ctspec$param <- vapply(ctspec$param, .ctCellSpecToPipe, character(1),
    USE.NAMES = FALSE)

  freeparams<-is.na(ctspec[,'value'])
  
  ctspec$transform<- NA
  ctspec$multiplier <- 1
  ctspec$meanscale <- 1
  ctspec$transform <- 0
  ctspec$offset <- 0
  ctspec$inneroffset <- 0
  
  
  ######### STAN parameter transforms
  for(pi in 1:length(ctspec$matrix)){
    if(freeparams[pi]){
      if(ctspec$matrix[pi] %in% c('T0MEANS','MANIFESTMEANS','TDPREDEFFECT','CINT')) {
        ctspec$meanscale[pi] <-10
      }
      if(ctspec$matrix[pi] %in% c('LAMBDA')) {
        ctspec$offset[pi] <- 0.5
        ctspec$meanscale[pi] <- 5
      }
      if(ctspec$matrix[pi] %in% c('THRESHOLDS')) {
        if(ctspec$col[pi] == 1) ctspec$meanscale[pi] <- 10
        if(ctspec$col[pi] > 1) {
          ctspec$transform[pi] <- 1
          ctspec$meanscale[pi] <- 2
          ctspec$multiplier[pi] <- 2
        }
      }
      
      if(ctspec$matrix[pi] %in% c('DIFFUSION','MANIFESTVAR', 'T0VAR')) {
        if(ctspec$row[pi] != ctspec$col[pi]){
          ctspec$transform[pi] <- 3
          ctspec$multiplier[pi] <- 2
          ctspec$offset[pi] <- -1
          ctspec$meanscale[pi] <-1
        }
        if(ctspec$row[pi] == ctspec$col[pi]){
          ctspec$transform[pi] <- 1
          ctspec$meanscale[pi] <- 2
          ctspec$multiplier[pi] <- 5
          ctspec$offset[pi] <- 1e-10
          if(ctspec$matrix[pi] %in% c('DIFFUSION')) ctspec$multiplier[pi] <-10
        }
      }
      if(ctspec$matrix[pi] %in% c('DRIFT')) {
        if(ctspec$row[pi] == ctspec$col[pi]){
          if(continuoustime==TRUE) {
            ctspec$transform[pi] <- 1
            ctspec$meanscale[pi] <- -2
            ctspec$multiplier[pi] <- -2
            ctspec$offset[pi] <- -1e-6
          }
          if(continuoustime==FALSE) {
            ctspec$transform[pi] <- 3
            ctspec$meanscale[pi] <- 2
            ctspec$offset[pi] <- 0
          }
        }
        if(ctspec$row[pi] != ctspec$col[pi]){
          ctspec$transform[pi] <- 0
          ctspec$meanscale[pi] <- 1
        }
      }
    }
  }
  
  ctspec[!is.na(ctspec$value),c('transform','multiplier','meanscale','offset','inneroffset')] <- NA
  
  
  for(ri in 1:nrow(ctspec)){ #convert back to text for new approach
    if(!is.na(as.integer(ctspec$transform[ri]))){
      ctspec$transform[ri] <- Simplify(tform(parin = 'param',
        transform = as.integer(ctspec$transform[ri]),
        multiplier = ctspec$multiplier[ri],
        meanscale = ctspec$meanscale[ri],
        offset = ctspec$offset[ri],
        inneroffset =ctspec$inneroffset[ri],
        singletext = TRUE))
    }
  }
  ctspec$multiplier <- NULL
  ctspec$meanscale <- NULL
  ctspec$offset <- NULL
  ctspec$inneroffset <- NULL
  
  
  
  nparams<-sum(freeparams)

  ctspec$indvarying<-FALSE
  ctspec$indvarying[!is.na(ctspec$transform) & ctspec$matrix %in% c('T0MEANS','MANIFESTMEANS','CINT')] <- TRUE
  
  ctspec$sdscale<-NA
  ctspec$sdscale[is.na(ctspec$value)]<-1


  # One sdscale per id element, defaulting to 1. `sdscale` is the subject
  # level's; each grouping level above it gets `sdscale_<idname>`, the same
  # rectangular encoding `indvarying_<idname>` uses. A genuine list column
  # would be a more literal "one vector per parameter" but would break every
  # consumer of `pars`, which indexes it as an ordinary data frame.
  #
  # `indvarying_<idname>` is created here too, defaulting to FALSE. Defaulting
  # an outer level the way `indvarying` defaults -- TRUE for T0MEANS,
  # MANIFESTMEANS and CINT -- would silently make every model with a grouping
  # id enormously parameterised, so an outer level varies only where it is
  # asked to.
  for(nm in if(length(ctm$id) > 1) ctm$id[-1] else character()){
    ctspec[[paste0('sdscale_', nm)]] <- ctspec$sdscale
    ctspec[[paste0('indvarying_', nm)]] <- FALSE
  }
  
  
  
  
  if(n.TIpred > 0) {
    tipredspec<-matrix(TRUE,ncol=n.TIpred,nrow=1)
    colnames(tipredspec)<-paste0(TIpredNames,'_effect')
    ctspec<-cbind(ctspec,tipredspec,stringsAsFactors=FALSE)
    # Character, not logical. An effect can now be free ('TRUE'), free and named
    # ('myeffect', so two parameters can share one), fixed ('4.3') or absent
    # ('FALSE'); see R/ctTipredEffect.R. A logical could say only the first and
    # last of those, and a numeric could not tell 'free' from 'fixed at 1'.
    for(predi in TIpredNames){
      ctspec[[paste0(predi,'_effect')]] <- as.character(tipredDefault)
    }
  }
  
  
  
  getwords <- function(x) gsub('\\W+',',',x)
  
  
  for(pi in 1:nrow(ctspec)){ #check for complex / split specifications
    if(grepl('|',ctspec$param[pi],fixed=TRUE)){
      ctspec$param[pi] <- gsub('$',' ',ctspec$param[pi])
      split = strsplit(ctspec$param[pi],split = '|',fixed=TRUE)[[1]]
      split=sapply(split,function(x) gsub(' ','',x))
      
      tisplit <- NA
      if(length(split) > 4){ #check for ti pred spec in splits
        if(n.TIpred < 1 || length(split) > 5) stop(paste0('Param spec has too many separators!  ', ctspec$param[pi]))
        tisplit <- split[5]
        split <- split[1:4]
      }
      
      # A number in the first element fixes the cell, and nothing after it can
      # apply.
      #
      # A fixed value cannot individually differ, cannot carry a covariate and
      # cannot have a transform -- there is no parameter left for any of those
      # to act on. So the value is heeded and the rest dropped, with a warning
      # naming exactly what was dropped, because the alternative is what this
      # branch used to do: `'2|param|TRUE|3'` became a *free* parameter named
      # "2" with value NA and indvarying TRUE, and `'2||||TI1=4.3'` became a
      # free parameter named "2" carrying a covariate of 4.3. Both silently.
      #
      # Only non-negative integers even got that far. `'-2|...'` and
      # `'-0.5|...'` both hit the simpleStateCheck error below, because a minus
      # sign and a decimal point are non-word characters -- so one
      # specification errored or silently misbehaved depending on which number
      # was written in it.
      fixedvalue <- suppressWarnings(as.numeric(split[1]))
      if(!is.na(fixedvalue)){
        ignored <- c(
          if(length(split) > 1 && nzchar(split[2])) paste0('transform ',split[2]),
          if(length(split) > 2 && nzchar(split[3])) paste0('indvarying ',split[3]),
          if(length(split) > 3 && nzchar(split[4])) paste0('sdscale ',split[4]),
          if(!is.na(tisplit) && nzchar(tisplit)) paste0('tipred effects ',tisplit))
        ctspec$value[pi] <- fixedvalue
        ctspec$param[pi] <- NA
        # Cleared here rather than left to the earlier pass that blanks the
        # transform of a valued cell: that pass runs before this loop, when the
        # cell still looked like a free parameter whose name happened to
        # contain separators.
        ctspec$transform[pi] <- NA
        for(cl in grep('^indvarying', names(ctspec), value=TRUE)) ctspec[[cl]][pi] <- FALSE
        for(cl in grep('^sdscale', names(ctspec), value=TRUE)) ctspec[[cl]][pi] <- NA
        for(cl in grep('_effect$', names(ctspec), value=TRUE)) ctspec[[cl]][pi] <- 'FALSE'
        if(length(ignored)) warning(paste0(ctspec$matrix[pi],'[',ctspec$row[pi],',',
          ctspec$col[pi],'] is fixed to ',format(fixedvalue),', so ',
          paste0(ignored,collapse=', '),
          if(length(ignored) > 1) ' were' else ' was',
          ' ignored: a fixed value cannot be transformed, cannot vary between ',
          'individuals and cannot carry a predictor effect. Give the cell a ',
          'parameter name if any of those were meant.'), call.=FALSE)
        next
      }

      if(grepl('\\W',split[1]) || 
          any(sapply(latentNames,function(x){ #if symbols or latent states
        grepl(paste0('\\b(',x,')\\b'),split[1])
      }))) {
        if(!simpleStateCheck(split[1])) stop(paste0(split[1],' invalid -- Matrix elements involving multiple parameters / latent states cannot have | separators -- transformations should be specified as part of the first element, indvarying and tipredeffects must be specified in the corresponding singular PARS matrix elements.'))
      }
      
      nonzero <- which(!split %in% '')
      ctspec[pi,c('param','transform','indvarying','sdscale')] <- 
        list(NA,ctspec$transform[pi],ctspec$indvarying[pi],1) #base values
      ctspec[pi,c('param','transform','indvarying','sdscale')[1:(length(split))]][nonzero] <- 
        split[1:(length(split))][nonzero] #update with non zero length split elements
      
      whichtipreds <- c()
      timessage <- c()
      if(!is.na(tisplit)){
        # `TI1, TI2=4.3, TI3=myeffect` -- a bare name is a free effect as it
        # always was, `=` gives it a value or a name of its own.
        #
        # Deliberately not `getwords()`, which every other field here uses:
        # it replaces each non-word character with a comma, so `TI1=4.3` came
        # out as `TI1,4,3` and there was no way to write a value at all.
        tisplit <- trimws(strsplit(tisplit, ',', fixed=TRUE)[[1]])
        ctspec[pi,paste0(TIpredNames,'_effect')] <- 'FALSE' #first set all off
        for(entry in tisplit){
          if(!nzchar(entry)) next
          parts <- trimws(strsplit(entry, '=', fixed=TRUE)[[1]])
          ti <- parts[1]
          if(!ti %in% TIpredNames) stop(ti,' is not a time independent predictor!')
          spec <- if(length(parts) > 1 && nzchar(parts[2])) parts[2] else 'TRUE'
          ctspec[pi,paste0(ti,'_effect')] <- spec
          whichtipreds <- c(whichtipreds,
            if(identical(spec,'TRUE')) ti else paste0(ti,'=',spec))
        }
        if(is.null(whichtipreds)) whichtipreds <- 'NULL'
        timessage <- paste0(ctspec$param[pi],' tipred effects from: ', paste0(whichtipreds,collapse=', '))

      }
      
      
      message('Custom par ',ctspec$param[pi],' set as: ',paste0(
        colnames(ctspec[pi,c('param','transform','indvarying','sdscale')]),' = ',
        ctspec[pi,c('param','transform','indvarying','sdscale')],'; '),timessage)
    }
  }
  
  if(n.TIpred > 0 && !any(.ctTipredEffectActive(ctspec[,paste0(TIpredNames,'_effect')]))) warning('TI predictors included but no effects specified!')
  
  for(ri in 1:nrow(ctspec)){ #set NA's on complex params
    
    if(grepl('\\W',gsub('.','',ctspec$param[ri],fixed=TRUE)) || ctspec$param[ri] %in% latentNames){ #if non word chars or state ref
      ctspec$value[ri] <- NA
      if(!simpleStateCheck(ctspec$param[ri])) ctspec$transform[ri] <-NA #allow transforms for simplestates
      if(simpleStateCheck(ctspec$param[ri])) ctspec$transform[ri] <-'param' #allow transforms for simplestates
    }
  }
  

  # Refuse a circular dependency here rather than later: the generated Stan
  # program breaks the T0MEANS / state / PARS loop by evaluation order, so a
  # cell caught in a real loop is read before anything writes it and returns a
  # plausible number instead of an error. See R/ctModelCycleCheck.R.
  .ctCheckModelCycles(ctspec, latentNames, n.latent)

  ctspec$indvarying <- as.logical(ctspec$indvarying)
  
  out<-list(pars=ctspec,n.latent=n.latent,n.manifest=n.manifest,n.TIpred=n.TIpred,n.TDpred=n.TDpred,
    latentNames=latentNames,manifestNames=manifestNames,
    TIpredNames=TIpredNames,TDpredNames=TDpredNames,
    # `id` may name a hierarchy rather than a single column: innermost
    # (subject) first, then any grouping levels above it -- e.g.
    # id = c('subject','study'). `subjectIDname` stays the scalar every
    # existing code path expects, and the rest is carried alongside so nothing
    # downstream has to learn about vectors it does not use.
    subjectIDname=ctm$id[1],
    groupIDnames=if(length(ctm$id) > 1) ctm$id[-1] else character(),
    timeName=ctm$time,
    continuoustime=continuoustime,
    manifesttype=ctmodelobj$manifesttype,
    # Zero for every non-ordinal variable, so `any(manifesttype %in% 2)` is the
    # only thing that ever has to be tested before reading it.
    ncategories=if(is.null(ctmodelobj$ncategories))
      rep(0L, n.manifest) else as.integer(ctmodelobj$ncategories),
    # Infinite for every non-censored variable, so a limit cannot apply where
    # it was not asked for.
    censormin=if(is.null(ctmodelobj$censormin))
      rep(-Inf, n.manifest) else as.numeric(ctmodelobj$censormin),
    censormax=if(is.null(ctmodelobj$censormax))
      rep(Inf, n.manifest) else as.numeric(ctmodelobj$censormax))
  class(out)<-'ctStanModel'
  
  out$tipredeffectscale <- 1
  out$tipredsimputedscale <- 1
  
  # out$popsdpriorscale <- 1
  out$rawpopsdbase <- 'normal(0,1)' #'cauchy(0,1)'
  out$rawpopsdbaselowerbound <- NA
  out$rawpopsdtransform <- 'log1p_exp(2*rawpopsdbase-1) .* sdscale' #'log(1+exp(2*rawpopsdbase)) .* sdscale' #'exp(rawpopsdbase * 2 -2) .* sdscale' # 'rawpopsdbase .* sdscale' #
  # out$stationarymeanprior <- NA
  # out$stationaryvarprior <- NA
  out$covmattransform <- 'rawcorr'
  # The population covariance, one row and column per varying parameter. Held
  # beside `pars` rather than in it, deliberately: everything that enumerates
  # free parameters walks `pars`, and the population covariance is not one of
  # those -- both backends build it from their own parameterisation. See
  # R/ctModelRawPopVar.R.
  out <- .ctModelRawPopVarSync(out)
  out[['matrices']] <- .ctModelMatricesPlaceholder()
  # out$NOrdinalIntegrationPoints <- 9L
  
  return(out)
}

#' @export
ctStanModel <- ctModelConvertOMX


# --- which parameters vary, at which level -----------------------------------

# The subject level is `indvarying`; each grouping level above it is
# `indvarying_<idname>`, created in `ctModelConvertOMX()` above.
#
# Everything that asks "does this model have random effects" has to ask about
# every level, not just the subject's. A study effect with no subject effect is
# an ordinary multilevel model -- one intercept per study, subjects treated as
# exchangeable within a study -- and testing `indvarying` alone refused it as
# having nothing to integrate over.
.ctVaryingColumns <- function(model){
  c('indvarying',
    if(length(model$groupIDnames)) paste0('indvarying_', model$groupIDnames))
}

# Which rows are free and varying at the given levels. A cell with a value is
# fixed and never varies, whatever its flag says. `%in% TRUE` rather than
# `isTRUE` per element because these columns carry NA on rows added after the
# model was built.
.ctVaryingRows <- function(model, columns = .ctVaryingColumns(model)){
  columns <- columns[columns %in% names(model$pars)]
  free <- is.na(model$pars$value)
  out <- rep(FALSE, nrow(model$pars))
  for(cl in columns) out <- out | (model$pars[[cl]] %in% TRUE & free)
  out
}

# Does anything vary at any of the given levels?
.ctAnyVarying <- function(model, columns = .ctVaryingColumns(model))
  any(.ctVaryingRows(model, columns))

# The same, named. For messages that have to say which parameters a route
# cannot carry; an unnamed row cannot be reported and is dropped here only.
.ctVaryingParams <- function(model, columns = .ctVaryingColumns(model)){
  rows <- .ctVaryingRows(model, columns)
  unique(model$pars$param[rows & !is.na(model$pars$param)])
}

# The grouping levels above the subject, alone. Used where a route handles the
# subject level and nothing else, so it can say which parameters it would have
# had to drop rather than dropping them.
.ctOuterVaryingColumns <- function(model){
  cols <- .ctVaryingColumns(model)
  cols[-1L]
}
