ctModeltoNumeric <- function(ctmodelobj){
  ###read in model
  #set any matrices to numeric elements
  #
  # Free parameters need *a* value before anything can be simulated. Zero,
  # matching ctsem's behaviour before commit 437bfbc0 so existing scripts keep
  # producing the data they always have -- except DRIFT, which cannot be a
  # literal zero: a zero DRIFT is singular, so `fQinf()` cannot solve for the
  # asymptotic covariance and generation fails outright. DRIFT diagonals use
  # the same near-zero convention `ctModel0DRIFT()` applies to a fixed DRIFT
  # diagonal of exactly zero.
  #
  # `.ctGenerateDefaults()` holds the per-matrix choices, shared with the julia
  # generation path so both produce the same kind of data from an underspecified
  # model. They are defaults for *simulation*, not estimates: the point is that
  # generation does not error, and anything a user cares about they should set.
  defaults <- .ctGenerateDefaults(continuoustime =
    if(!is.null(ctmodelobj$continuoustime)) ctmodelobj$continuoustime else TRUE)
  filled <- character()
  sapply(names(ctmodelobj), function(x){
    if(is.matrix(ctmodelobj[[x]])){
      m <- ctmodelobj[[x]]
      free <- suppressWarnings(is.na(as.numeric(m)))
      dim(free) <- dim(m)
      if(any(free)){
        spec <- defaults[[x]]
        idx <- which(free, arr.ind=TRUE)
        for(k in seq_len(nrow(idx))){
          i <- idx[k,1]; j <- idx[k,2]
          value <- if(is.null(spec)) 0 else
            if(i == j) spec$diagonal else spec$offdiagonal
          m[i,j] <- value
          # Every filled cell is named, including the ones filled with zero.
          # Reporting only the non-zero ones made the common case silent: a
          # model leaving only MANIFESTVAR or T0MEANS free said nothing at all,
          # so a script -- or a test under a set.seed() -- generated from
          # values it never stated and had no way to notice when those values
          # changed underneath it. That is how test-binary-binary-mix.R came to
          # characterise a dataset that no longer existed.
          filled <<- c(filled, sprintf('%s[%d,%d]=%s', x, i, j, format(value)))
        }
      }
      ctmodelobj[[x]] <<- matrix(as.numeric(m),nrow=nrow(m), ncol=ncol(m))
    }
  })
  if(length(filled)) message(length(filled),
    ' free parameters were given generating values: ',
    paste(utils::head(filled, 8), collapse=', '),
    if(length(filled) > 8) ', ...' else '',
    '. Set them in the model if they matter.')
  
  return(ctmodelobj)
}

#' Generate data from a ctstanmodel object
#'
#' Generate data from a ctstanmodel object.
#' \code{ctStanGenerate} is maintained as a backward-compatible alias.
#'
#' @param cts \code{\link{ctModelConvertOMX}}, \code{\link{ctModel}}, or
#' \code{\link{ctStanFit}} object.
#' @param datastruct long format data structure as used by ctsem. 
#' Not used if cts is a ctStanFit object.
#' @param is Deprecated and ignored, with a warning if set. Importance
#' sampling reweights draws taken from a gaussian approximation onto a target
#' that cannot be sampled directly. Fitted to an empty dataset there is no
#' such target: the likelihood contributes nothing and what remains is the
#' prior, which ctsem holds as independent univariate densities in the raw
#' space these draws are taken in. There is no approximation to correct.
#' @param fullposterior Generate from the full prior, or from its mode (the raw
#' origin)?
#' @param nsamples How many samples to generate?
#' @param parsonly If TRUE, only return samples of raw parameters, don't generate data.
#' @param cores Number of cpu cores to use.
#' @param backend Which engine turns a parameter draw into data: \code{'auto'}
#' (the default -- \code{'julia'} when a julia session is available,
#' \code{'stan'} otherwise), \code{'julia'} or \code{'stan'}. The draws
#' themselves do not depend on this: they come from the model's prior, which is
#' the same object either way. \code{parsonly=TRUE} returns those draws and
#' generates nothing, so it does not use this argument.
#'
#' @return List containing Y, an array of nsamples by data rows by manifest
#' variables, and llrow, an array of nsamples by data rows log likelihoods.
#' With \code{parsonly=TRUE}, the prepared model carrying the prior draws
#' under \code{$stanfit$rawposterior} and \code{$stanfit$transformedpars}.
#' @aliases ctStanGenerate
#' @export
#'
#' @examples
#' \donttest{
#' #generate and plot samples from prior predictive
#' priorpred <- ctGenerateFromPriors(cts = ctstantestfit,cores=2,nsamples = 50)
#'}
ctGenerateFromPriors <- function(cts,datastruct=NA, is=FALSE,
  fullposterior=TRUE, nsamples=200, parsonly=FALSE,cores=2,
  backend=c('auto','julia','stan')){

  backend <- match.arg(backend)

  # Named rather than asserted. "Not a ctStanModel object" (from ctFit(), further
  # downstream) is opaque when the caller is plainly holding a fit; what it means
  # is that this function reads stan fit structures (ctstanmodelbase, standata,
  # args) that a julia fit does not carry.
  if(inherits(cts, 'ctJuliaFit')) stop(
    'This function is not available for julia backend fits yet: it reads the ',
    'stan fit structures (ctstanmodelbase, standata, args) that a julia fit ',
    'does not carry. ctFitCheck(), ctFitCheckCov(), ctACFresiduals() and ',
    'ctPostPredPlots() do work on a julia fit -- ctFitCheck() simply omits its ',
    'prior predictive panel, which is the one thing that needs this function. ',
    "Note that backend='julia' here selects the generator, not the kind of fit ",
    'this accepts.',
    call.=FALSE)

  # `is` selected stan's optimize-then-importance-sample route, back when
  # ctFit() had an `optimcontrol$is` to select it with. It is not rewired,
  # because nothing below approximates anything for it to correct.
  if(!identical(is, FALSE)) .Deprecated(msg = paste0(
    'The `is` argument of ctGenerateFromPriors() is deprecated and ignored. ',
    'Importance sampling corrects a gaussian approximation to a posterior; ',
    'this function draws from the prior directly, so there is no ',
    'approximation to correct.'))

  if('ctStanFit' %in% class(cts)){
    # Three places, in order, because a fit made before `$args$resolved`
    # existed has only the other two -- and `ctstantestfit`, the fit every
    # example on this page uses, is one of them.
    priors <- cts$args$resolved$priors
    if(is.null(priors)) priors <- cts$args$priors
    if(is.null(priors) && !is.null(cts$standata$priors)) priors <- as.logical(cts$standata$priors)
    datastruct <- standatatolong(cts$standata, origstructure=TRUE, ctm=cts$ctstanmodelbase)

    cts <- cts$ctstanmodelbase

  } else priors<-TRUE

  if(!is.null(priors) && !as.logical(priors)) stop('Priors disabled, cannot sample from prior!')

  # -99, not NA. The placeholder says "generate a value for this row"; NA says
  # the row is missing, and a missing row is one the generator steps over.
  datastruct[,cts$manifestNames] <- -99

  cts$TIpredAuto <- 0L

  # `parsonly` asks about parameters and generates nothing, so there is no
  # generator to choose; the prior draws are the same object either way. It
  # takes the stan preparation because the constrained draws it returns are
  # stan shaped -- ctPlotPosterior() reads `$stanfit$transformedpars`.
  if(backend == 'auto') backend <-
    if(isTRUE(tryCatch(ctJuliaStatus()$available, error=function(e) FALSE)))
      'julia' else 'stan'
  if(parsonly) backend <- 'stan'

  # Prepared, not fitted.
  #
  # This used to fit the model to a one-row-per-subject dataset with every
  # manifest set to NA, so that the "posterior" it optimised would be the
  # prior, and then draw from the hessian covariance at the mode. That was a
  # way to get a stan model instance to draw through, not a statement about
  # the prior, and it cost an optimisation and a hessian to arrive at an
  # answer that is written down in the model: ctsem's raw priors are
  # independent normals, so the draws are one rnorm() each. It also went
  # wrong in the ways an unnecessary optimisation does -- an empty dataset
  # with priors off is a completely flat objective, and the optimiser wanders.
  #
  # `fit=FALSE` gives the same prepared data and model this needs, over the
  # real row and missingness structure rather than a dummy one, without
  # running anything.
  args <- cts$args
  args$model <- cts
  args$ctstanmodel <- NULL
  args$datalong <- datastruct
  args$fit <- FALSE
  args$priors <- TRUE
  args$optimize <- TRUE
  args$intoverstates <- TRUE
  args$intoverpop <- TRUE
  args$cores <- cores
  args$backend <- backend

  prepared <- do.call(ctFit, args)

  if(backend == 'julia'){
    spec <- .ctBackendSpec(prepared)
    if(is.null(spec$priors) || !length(spec$priors$index)) stop(
      'The prepared model carries no prior specification, so there is nothing ',
      "to draw from. Use backend='stan'.", call.=FALSE)
    # The spec accounts for every free parameter or refuses to be built, so its
    # index is exactly seq_len(npar).
    npar <- length(spec$priors$index)
    draws <- .ctPriorRawDraws(prepared$standata, npar, nsamples)
    prepared$estimate <- list(raw = rep(0, npar), rawposterior = draws)
    message('Generating ', nsamples, ' datasets from the prior, backend julia.')
    ppf <- .ctBackendGenerateFromFit(prepared, nsamples=nsamples,
      fullposterior=fullposterior, cores=cores)
    return(list(Y = ppf$generated$Y, llrow = ppf$generated$llrow))
  }

  # stan. `sm` follows ctFit()'s own rule: the built-in model unless this
  # model's text differs from it, in which case there is nothing compiled to
  # reuse.
  standata <- prepared$standata
  if(isTRUE(prepared$setup$recompile)){
    message('Compiling model -- usually ~ 1 min.')
    sm <- rstan::stan_model(model_code = prepared$stanmodeltext)
  } else sm <- stanmodels$ctsm

  npar <- rstan::get_num_upars(stan_reinitsf(sm, standata))
  draws <- .ctPriorRawDraws(standata, npar, nsamples)

  prepared$stanmodel <- sm
  prepared$stanfit <- list(rawest = rep(0, npar), rawposterior = draws)
  class(prepared) <- c('ctStanFit','ctFit')

  if(parsonly){
    # dokalman=FALSE matches what the fitted route computed here: it took
    # `dokalman` from `savescores`, which is off for this data.
    prepared$stanfit$transformedpars <- suppressMessages(stan_constrainsamples(
      sm = sm, standata = standata, samples = draws, cores = cores,
      savescores = FALSE, savesubjectmatrices = FALSE, dokalman = FALSE,
      pcovn = FALSE))
    return(prepared)
  }

  ppf <- ctGenerateFromFit(fit = prepared, nsamples = nsamples,
    fullposterior = fullposterior, cores = cores)

  # Named by ctGenerateFromFit() and by .ctBackendGenerateFromFit() alike, and
  # named correctly: dim 1 is the sample. This used to relabel dim 1
  # `datapoints` and dim 2 `samples`, which is the wrong way round and
  # contradicted the @return text directly above.
  list(Y = ppf$generated$Y, llrow = ppf$generated$llrow)
}

# ctsem's prior over the raw parameters, sampled directly.
#
# `.ctBackendPriorSpec()` reads the index and scale of each raw parameter's
# prior out of the prepared data. It is the same spec the julia engine is
# handed and the same layout the generated stan model uses, and every family in
# it is normal: the density the engines evaluate is
# `weight * sum(dnorm(value / scale, log = TRUE))`, which as a distribution over
# `value` is `normal(0, scale / sqrt(weight))`. `weight` is
# `priormod / nsubsets`, and is 1 for anything but a subsetted fit.
#
# The spec accounts for every free parameter or refuses to be built, so there
# is no unnamed remainder; the zero-filled matrix is there to make that
# assumption visible rather than to be relied on.
.ctPriorRawDraws <- function(standata, npar, nsamples){
  .ctPriorRejectUnsamplable(standata)
  spec <- .ctBackendPriorSpec(standata, npar)
  sdvec <- spec$scale / sqrt(spec$weight)
  draws <- matrix(0, nrow = nsamples, ncol = npar)
  draws[, spec$index] <- stats::rnorm(nsamples * length(spec$index)) *
    rep(sdvec, each = nsamples)
  draws
}

# The prior families this cannot draw from, refused where the caller can act.
#
# `.ctBackendRejectLaplacePriors()` refuses the same models but advises
# `backend='stan'`, which is wrong here: the obstacle is the density, not the
# engine, and no backend samples it.
#
# Worth saying plainly what changed, because it reads like a lost capability
# and is not one. The fit-to-empty-data route this replaced did not sample
# these either -- it optimised to the mode and drew from the hessian covariance
# there, so a laplaceprior parameter's "prior" draws came back gaussian
# whatever laplaceprior said, and nothing reported it.
.ctPriorRejectUnsamplable <- function(standata){
  laplaceprior <- as.integer(standata$laplaceprior)
  if((length(laplaceprior) && any(laplaceprior == 1L)) ||
      isTRUE(as.integer(standata$laplacetipreds)[1L] == 1L)) stop(
    'ctGenerateFromPriors() draws from the prior directly, and the smoothed ',
    'double exponential ctsem uses for laplaceprior parameters is not a ',
    'density it can draw from. No backend changes that -- the obstacle is the ',
    'prior rather than the engine. The route this replaced did not sample them ',
    'either: it drew from a gaussian approximation, so those draws were normal ',
    'whatever laplaceprior said. Drop laplaceprior for the affected matrices ',
    'to generate from priors.', call.=FALSE)
  invisible(TRUE)
}

#' @export
ctStanGenerate <- ctGenerateFromPriors



#' ctGenerate
#' 
#' This function generates data according to the specified ctsem model object. 
#' 
#' @param ctmodelobj ctsem model object from \code{\link{ctModel}}.
#' @param n.subjects Number of subjects to output.
#' @param burnin Number of initial time points to discard (to simulate stationary data)
#' @param dtmean Positive numeric. Median time interval (delta T) to use.
#' Intervals are drawn as \code{exp(rnorm(n, log(dtmean), logdtsd))}, so
#' \code{dtmean} is the median rather than the mean whenever \code{logdtsd}
#' is non-zero: the mean interval is \code{dtmean * exp(logdtsd^2/2)}, which
#' at \code{logdtsd = 0.6} is about 20 percent longer than \code{dtmean}.
#' @param logdtsd Numeric. Standard deviation of the log time interval. Zero
#' gives an equal-interval design.
#' @param dtmat Either NA, or numeric matrix of n.subjects rows and Tpoints-1 columns, 
#' containing positive numeric values for all time intervals between measurements. 
#' If not NA, dtmean and logdtsd are ignored.
#' @param Tpoints Optional number of time points to generate. If supplied, this overrides
#' any \code{Tpoints} stored in \code{ctmodelobj}. If not supplied, \code{ctGenerate}
#' uses \code{ctmodelobj$Tpoints} when available.
#' @param wide Logical. Output in wide format?
#' @param backend Which generator to use: \code{'r'}, \code{'julia'}, or
#' \code{'auto'} (the default). \code{'r'} integrates the linear system with a
#' matrix exponential -- fast, and exact for a linear Gaussian model, but it has
#' no measurement link, so a model with non-Gaussian indicators would silently
#' get continuous values, and a nonlinear (state-dependent) model is refused
#' outright. \code{'julia'} generates through the same filter that fits the
#' model, so it handles nonlinear dynamics and non-Gaussian (binary, ordinal,
#' count) indicators, and requires \code{ctmodelobj} to be a
#' \code{ctStanModel} (as returned by \code{ctModel(type='ct'/'dt')}) rather
#' than the matrix-list form. \code{'auto'} picks \code{'julia'} exactly when
#' the model is nonlinear or declares a non-Gaussian indicator, and \code{'r'}
#' otherwise -- the split is by capability, not preference, so linear Gaussian
#' models keep the seed-for-seed output every existing caller already gets.
#' @param intoverstates For \code{backend='julia'}: \code{'auto'} (the
#' default), \code{TRUE} or \code{FALSE}, choosing how the latent states are
#' handled while generating.
#'
#' \code{TRUE} draws each row from the filter's own one-step-ahead predictive
#' and lets the filter condition on the draw, so the dataset is an exact draw
#' from the density a fit maximises. \code{FALSE} samples the latent
#' trajectory from the process and then each observation from its conditional
#' distribution given the state at its row -- a draw from the model rather than
#' from the filter's approximation of it.
#'
#' \code{'auto'} picks \code{TRUE} exactly when the filter's predictive *is*
#' the model's: linear dynamics and Gaussian indicators. There the two agree in
#' distribution, and \code{TRUE} keeps the output every existing caller gets
#' for a given seed. Otherwise it picks \code{FALSE}. A categorical indicator
#' makes the measurement update an assumed-density projection, which moves the
#' state it conditions on -- and with an unbounded indicator (a count) an
#' improbable draw can move it far enough that the following rows are drawn
#' from a rate that has already run away. A state-dependent drift makes the
#' prediction a moment approximation in the same way.
#'
#' \code{backend='r'} already generates this way and ignores the argument.
#' @details Covariance related matrices are treated as Cholesky factors.
#' TRAITTDPREDCOV and TIPREDCOV matrices are not accounted for, at present. 
#' The first 1:n.TDpred rows and columns of TDPREDVAR are used for generating
#' tdpreds at each time point. 
#' @examples 
#' #generate data for 2 process model, each process measured by noisy indicator, 
#' #stable individual differences in process levels.
#' 
#' generatingModel<-ctModel(Tpoints=8,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'  MANIFESTVAR=diag(.1,2),
#'  LAMBDA=diag(1,2),
#'  DRIFT=matrix(c(-.2,-.05,-.1,-.1),nrow=2),
#'  DIFFUSION=matrix(c(1,.2,0,4),2),
#'  CINT=matrix(c(1,0),nrow=2),
#'  T0MEANS=matrix(0,ncol=1,nrow=2),
#'  T0VAR=diag(1,2))
#'
#' nsubjects <- 15
#' traitChol <- matrix(c(.5,.2,0,.8),nrow=2)
#' subjectCint <- t(replicate(nsubjects, as.numeric(traitChol %*% rnorm(2))))
#' datalist <- vector("list", nsubjects)
#' for(i in seq_len(nsubjects)){
#'   subjectModel <- generatingModel
#'   #through $matrices: a ctStanModel's specification is $pars, which every
#'   #top level matrix is rebuilt from below, so `subjectModel$CINT <- ` would
#'   #be discarded and every subject generated with the same CINT.
#'   subjectModel$matrices$CINT <- matrix(subjectCint[i,], ncol = 1)
#'   d <- ctGenerate(subjectModel,n.subjects=1,burnin=10)
#'   d[,'id'] <- i
#'   datalist[[i]] <- d
#' }
#' data <- do.call(rbind, datalist)
#' @export

ctGenerate<-function(ctmodelobj,n.subjects=100,burnin=0,dtmean=1,logdtsd=0,dtmat=NA,
  Tpoints=NULL, wide=FALSE, backend=c('auto','r','julia'), intoverstates='auto'){
  backend <- match.arg(backend)
  # `auto` routes to julia only what the generator below cannot do. That
  # generator integrates the linear system with a matrix exponential, which is
  # exact for a linear model and simply inapplicable to a state-dependent one;
  # the engine filters and generates both. Defaulting to julia for everything
  # would change the numbers under every existing caller for no gain on the
  # models they use, so the split is by capability rather than by preference.
  nonlinear <- isTRUE(try(ctModelIsNonlinear(ctmodelobj), silent=TRUE))
  # A categorical indicator is the same situation as a nonlinear one: the
  # generator below integrates a linear Gaussian system and has no notion of a
  # link, so it produces continuous values for a manifest the model declares
  # binary or ordinal -- silently, which is the worst of both. Measured before
  # this line existed: a model with `manifesttype = 1` generated values with a
  # mean of 0.063 and no zeros or ones among them.
  categorical <- !is.null(ctmodelobj$manifesttype) &&
    any(ctmodelobj$manifesttype > 0)
  if(backend == 'auto') backend <- if(nonlinear || categorical) 'julia' else 'r'
  # `intoverstates='auto'` asks the same question the backend choice asks, one
  # level down: is the filter's one-step-ahead predictive the model's own?
  #
  # It is exactly when the dynamics are linear and every indicator Gaussian.
  # There the filter's predictive is exact, the two routes agree in
  # distribution, and TRUE is kept -- it preserves the seed-for-seed output
  # every existing caller gets, and the round-trip identity that a generated
  # dataset's likelihood is the one reported while generating it.
  #
  # It is not, for anything else. A categorical indicator makes the update an
  # assumed-density projection, which moves the state it conditions on and
  # can run away on an unbounded one -- the count case this route exists for.
  # A state-dependent drift makes the *prediction* a moment approximation in
  # the same way. In both, sampling the trajectory and then the observations
  # given it is a draw from the model where the filter route is a draw from
  # the filter's approximation of it.
  if(identical(intoverstates,'auto')) intoverstates <- !(nonlinear || categorical)
  intoverstates <- isTRUE(as.logical(intoverstates)[1])
  if(backend == 'r' && categorical){
    warning('This model declares non-Gaussian indicators (manifesttype ',
      "1, 2 or 3) and backend='r' generates continuous values for them: the R ",
      'generator has no measurement link. Use backend="julia" for ',
      'binary, ordinal or count data.', call.=FALSE)
  }
  if(backend == 'r' && nonlinear) {
    stop("This model is nonlinear, and ctGenerate's own generator integrates a ",
      "linear system: it has no way to apply a state-dependent specification. ",
      "Use backend='julia' (the default for such models), which generates ",
      "through the same filter that fits them.", call.=FALSE)
  }
  if(backend == 'julia'){
    if(!'ctStanModel' %in% class(ctmodelobj)) {
      stop("backend='julia' generation needs a ctModel(type='ct'/'dt') object, ",
        "which carries the parameter specification the engine reads. The ",
        "matrix-list form does not.", call.=FALSE)
    }
    modelTpoints <- if(!is.null(Tpoints) && !is.na(Tpoints[1])) Tpoints[1] else
      if(!is.null(ctmodelobj$Tpoints) && !is.na(ctmodelobj$Tpoints[1]))
        ctmodelobj$Tpoints[1] else
          stop('Tpoints not found in ctmodelobj and no Tpoints argument supplied. Provide Tpoints explicitly.')
    fullTpoints <- burnin + as.integer(modelTpoints)
    times <- lapply(seq_len(n.subjects), function(si){
      dtvec <- if(is.na(dtmat[1])) exp(rnorm(fullTpoints,log(dtmean),logdtsd)) else
        c(rep(1,burnin), dtmat[si,,drop=TRUE])
      tv <- numeric(fullTpoints)
      for(t in 2:fullTpoints) tv[t] <- round(tv[t-1] + dtvec[t-1], 6)
      tv
    })
    out <- .ctGenerateJulia(ctmodelobj, n.subjects, times,
      intoverstates = intoverstates)
    if(burnin > 0){
      keep <- unlist(lapply(seq_len(n.subjects), function(si)
        (si-1)*fullTpoints + (burnin+1):fullTpoints))
      out <- out[keep, , drop=FALSE]
      # Time restarts at zero for each subject once the burnin is dropped, as
      # the generator below does. Leaving it running from the burnin would make
      # the first observed interval look like the whole burnin period.
      for(si in seq_len(n.subjects)){
        rows <- (si-1)*(fullTpoints-burnin) + seq_len(fullTpoints-burnin)
        out[rows,'time'] <- out[rows,'time'] - out[rows[1],'time']
      }
    }
    if(wide) return(ctLongToWide(out, id='id', time='time',
      manifestNames=ctmodelobj$manifestNames,
      TDpredNames=ctmodelobj$TDpredNames, TIpredNames=ctmodelobj$TIpredNames))
    return(out)
  }
  if('ctStanModel' %in% class(ctmodelobj)){
    # Reconstruct matrix-style slots when a ctStanModel is supplied.
    mlist <- listOfMatrices(ctmodelobj$pars)
    for(nm in names(mlist)){
      ctmodelobj[[nm]] <- mlist[[nm]]
    }
  }
  
  ctmodelobj <- ctModeltoNumeric(ctmodelobj)
  
  m <- ctmodelobj
  
  modelTpoints <- NULL
  if(!is.null(m$Tpoints) && !is.na(m$Tpoints[1])) modelTpoints <- m$Tpoints[1]
  if(!is.null(Tpoints) && !is.na(Tpoints[1])) modelTpoints <- Tpoints[1]
  
  if(is.null(modelTpoints)) stop('Tpoints not found in ctmodelobj and no Tpoints argument supplied. Provide Tpoints explicitly.')
  if(length(modelTpoints) != 1 || !is.finite(modelTpoints) || modelTpoints < 1) stop('Tpoints must be a single finite value >= 1.')
  
  modelTpoints <- as.integer(modelTpoints)
  m$Tpoints <- modelTpoints
  fullTpoints<-burnin+m$Tpoints

  for(si in 1:n.subjects){
    
    if(is.na(dtmat[1])) dtvec<- exp(rnorm(fullTpoints,log(dtmean),logdtsd))
    if(!is.na(dtmat[1])) dtvec <- c(rep(1,burnin),dtmat[si,,drop=FALSE])
    time=rep(0,fullTpoints)
    for(t in 2:fullTpoints) time[t] = round(time[t-1] + dtvec[t-1],6)
    
    if(m$n.TDpred > 0) {
      tdpreds <- rbind(matrix(0,nrow=1+(burnin),ncol=m$n.TDpred)[-1,,drop=FALSE], #additional row added then removed in case no burnin
        matrix(m$TDPREDMEANS,ncol=m$n.TDpred))
    }
    
    # #convert to triangular...
    # m$T0VAR <- t(chol(m$T0VAR))
    # m$MANIFESTVAR <- t(chol(m$MANIFESTVAR))
    
    sm=m
    if(any(m$TRAITVAR != 0)) {
      traits = m$TRAITVAR %*% rnorm(m$n.latent,0,1)
      sm$CINT = m$CINT +  traits
      sm$T0MEANS = m$T0MEANS + m$T0TRAITEFFECT %*% traits
    }
    
    if(any(m$MANIFESTTRAITVAR != 0)) {
      sm$MANIFESTMEANS = m$MANIFESTMEANS + m$MANIFESTTRAITVAR %*% rnorm(m$n.manifest,0,1)
    }
    
    if(m$n.TIpred > 0) {
      tipreds <- m$TIPREDMEANS + m$TIPREDVAR %*% rnorm(m$n.TIpred,0,1)
      sm$CINT = m$CINT + m$TIPREDEFFECT %*% tipreds
    }
    
    manifests<-matrix(NA,fullTpoints,m$n.manifest)
    latents<-matrix(NA,fullTpoints,m$n.latent)
    
    sdat <- cbind(si,time,manifests,
      if(m$n.TDpred > 0) tdpreds,
      if(m$n.TIpred > 0) matrix(tipreds,byrow=TRUE,nrow=fullTpoints,ncol=m$n.TIpred))
    
    colnames(sdat) <- c('id','time',m$manifestNames,m$TDpredNames,m$TIpredNames)
    sdat <- data.frame(sdat)

    
    latents[1,] <- sm$T0MEANS+m$T0VAR %*% rnorm(m$n.latent)
    Qinf <- fQinf(sm$DRIFT,sm$DIFFUSION)
    
    for(i in 2:nrow(latents)){
      dtA=expm::expm(sm$DRIFT * (sdat$time[i]-sdat$time[i-1]))
      # message('dtA')
      # print(dtA)
      # message('dtCINT')
      # print(solve(sm$DRIFT,(dtA - diag(m$n.latent))) %*% sm$CINT)
      # message('dtG')
      # print(t(chol(fdtQ(Qinf,dtA))))
      latents[i,] <- dtA %*% latents[i-1,] +
        solve(sm$DRIFT,(dtA - diag(m$n.latent))) %*% sm$CINT + 
        t(chol(fdtQ(Qinf,dtA))) %*% rnorm(m$n.latent)
      # browser()
      if(m$n.TDpred > 0) latents[i,] <- latents[i,] + sm$TDPREDEFFECT %*% 
        t(as.matrix(sdat[i,m$TDpredNames,drop=FALSE]))
    }
    
    for(i in seq_len(nrow(sdat))){
      sdat[i,m$manifestNames] <- sm$LAMBDA %*% latents[i,] + sm$MANIFESTMEANS + 
        sm$MANIFESTVAR %*% rnorm(m$n.manifest)
    }
        
    
    
    
    sdat=sdat[(burnin+1):fullTpoints,]
    
    sdat[,'time'] = sdat[,'time'] - sdat[1,'time'] 
    
    if(si==1) datalong <- sdat else datalong <- rbind(datalong,sdat)
  }
  
  datalong<-as.matrix(datalong)
  
  
  if(wide==FALSE) return(datalong) else {
    datawide <- ctLongToWide(datalong = datalong,id = 'id',time = 'time',
      manifestNames = m$manifestNames, TDpredNames = m$TDpredNames,TIpredNames = m$TIpredNames)
    datawide <- ctIntervalise(datawide = datawide,Tpoints = m$Tpoints,n.manifest = m$n.manifest,n.TDpred = m$n.TDpred,n.TIpred = m$n.TIpred,
      manifestNames=m$manifestNames,TDpredNames=m$TDpredNames,TIpredNames=m$TIpredNames)
    return(datawide)
  }
}
