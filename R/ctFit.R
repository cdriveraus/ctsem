# One control list, one vocabulary.
#
# `optimcontrol` holds every optimiser setting for both backends, and a name in
# it means the same thing on each. Stan's names are on CRAN and fixed, so the
# julia side -- unreleased, and free -- was moved onto them. The correspondence
# is exact, because the two optimisers expose the same stopping rules:
#
#   optimcontrol   stan (mize)   julia engine   stops on
#   tol            abs_tol       f_tol          change in the objective
#   g_tol          grad_tol      g_tol          l2 norm of the gradient
#   x_tol          step_tol      x_tol          size of the parameter step
#   maxiter        max_iter      maxiter        iteration count
#   lbfgs_memory   memory        lbfgs_memory   curvature pairs L-BFGS keeps
#   initsd         initsd        (the same)     sd of the random start
#
# Each of those was previously honoured on one side and ignored on the other:
# `tol` and `initsd` did nothing on julia, and `g_tol`, `x_tol`, `maxiter` and
# `lbfgs_memory` existed only as `backendcontrol` names that reached nothing on
# stan. mize's gradient and step tolerances were pinned at zero by ctsem rather
# than absent from mize, which is what made the split look irreducible.
#
# The *defaults* stay each backend's own, and they are mirror images: stan stops
# on the objective (tol = 1e-8, g_tol off), julia on the gradient (g_tol = 1e-8,
# tol off). Leaving a name unset keeps that backend's default; setting it is
# honoured by both.
#
# `backendcontrol` is gone. It existed because these settings had nowhere to go
# and `optimcontrol` was stan-shaped, which made a second list the mismatch
# rather than the fix. Its two session controls were already ctJuliaSetup()'s
# own arguments -- `julia_project` is `project`, `restart_session` is
# `force = TRUE` -- and the rest are the rows above.

# Names ctFit() sets on optimcontrol itself on the way into stanoptimis() (see
# the `optimcontrol$cores <- cores` block below). A caller who passes one of
# these is overridden on stan and ignored on julia -- the same nothing on both
# sides -- so they are accepted in silence. `is` is here because the stan path
# refuses it by name a few lines below and .ctJuliaUnsupported() refuses it on
# julia; both messages are better than anything this function would say.
.ctOptimcontrolInert <- c('init','priors','plot','verbose','cores','matsetup',
  'standata','sm','is')

# Honoured by both backends, with the same meaning. `stallretries` is how many
# times a fit that stopped short of a maximum is tried again, with each
# backend's own retry: stan restarts from fresh random values, julia pulls the
# flat coordinates back and refits. It used to be listed below twice, stan-only
# and julia-only, and a list lookup returns the first -- so julia refused it
# with stan's message and stan dropped it with the julia-only names.
.ctOptimcontrolShared <- c('tol','g_tol','x_tol','maxiter','lbfgs_memory',
  'initsd','carefulfit','estonly','finishsamples','uncertainty',
  'uncertaintyDraws','uncertaintyControl','stallretries','stochastic')

# What is left after the vocabulary above is a genuine capability difference: a
# phase or a hook one backend has and the other does not. Each such name is
# listed here with the backend that has it and, where there is one, the value
# that asks for nothing.
#
# The rule is on the value, not the name. A value that asks for a capability the
# chosen backend lacks is refused; a value that merely *describes* what that
# backend already does is accepted, because it is true -- `gradient =
# 'adjoint'` on stan names the reverse-mode gradient Stan's autodiff already
# takes. (`stochastic` was the other example until julia gained the same
# phase; it is a shared name now.)
.ctOptimcontrolSplit <- function() list(

  # -- stan only: the stochastic-gradient family's tuning. `stochastic` itself
  # is shared (julia's sgd phase, `_ctsem_sgd`); these tune stan's own sgd(),
  # whose subsets and roughness targets the julia phase does not have.
  nsubsets = list(only = 'stan',
    inert = function(v) isTRUE(all(v == 1)),
    msg = paste0("splits the data for stan's stochastic optimizer, which the ",
      "julia backend does not have. Drop it")),
  subsamplesize = list(only = 'stan',
    inert = function(v) isTRUE(all(v >= 1)),
    msg = paste0("gives stan's first pass a proportion of the subjects; the ",
      "julia warm-up uses the priors over all of them instead. Use ",
      "optimcontrol$carefulfit")),
  lproughnesstarget = list(only = 'stan',
    inert = function(v) FALSE,
    msg = paste0("tunes stan's stochastic optimizer, which the julia backend ",
      "does not have. Drop it")),
  stochasticTolAdjust = list(only = 'stan',
    inert = function(v) FALSE,
    msg = paste0("tunes stan's stochastic optimizer, which the julia backend ",
      "does not have. Drop it")),
  parsteps = list(only = 'stan',
    inert = function(v) length(v) == 0L,
    msg = paste0("holds parameters at zero during a stepwise stan ",
      "optimisation, which the julia optimiser has no step for. Drop it")),
  stalltol = list(only = 'stan',
    inert = function(v) FALSE,
    msg = paste0("is the gradient at which stan calls a fit stalled and ",
      "restarts it; the julia backend judges a stall by its progress instead. ",
      "Use optimcontrol$stallwindow")),

  # -- julia only.
  gradient = list(only = 'julia',
    # Stan's autodiff is reverse mode, so 'adjoint' is a true description of it
    # and costs nothing to accept. 'forward' names a second gradient that only
    # the julia engine has.
    inert = function(v) identical(as.character(v)[1L], 'adjoint'),
    msg = paste0("selects the julia engine's gradient; stan has one, and it ",
      "is already reverse mode. Drop it, or use gradient='adjoint'")),
  datastart = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("reads starting values off the data through the julia ",
      "parameter table, which the stan path does not build. Pass inits, or ",
      "optimcontrol$initsd for the scale of the random start")),
  precondition = list(only = 'julia',
    inert = function(v) !isFALSE(v),
    msg = paste0("switches off the diagonal metric the julia optimiser takes ",
      "from the model's own transforms, so a raw unit means the same amount of ",
      "model in every coordinate. The stan path has no equivalent")),
  batch = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("switches off the julia optimiser's progressive batching, ",
      "which starts on a subset of subjects and grows it as the fit needs more ",
      "data. The stan path has no equivalent")),
  newton = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("switches off the julia optimiser's Newton finish on the exact ",
      "Hessian. The stan path has no equivalent")),
  restarts = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 0)),
    msg = paste0("sets how many random restarts a julia fit tries when it has ",
      "not converged. The stan path has no equivalent")),
  restartsd = list(only = 'julia',
    inert = function(v) TRUE,
    msg = paste0("sets the spread of the julia route's random restarts. The ",
      "stan path has no equivalent")),
  explosive_forward = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("takes forward-mode gradients for subjects with explosive ",
      "dynamics on the julia route. The stan path has no equivalent")),
  lbfgs_diagonal = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("learns the julia optimiser's initial inverse Hessian per ",
      "coordinate. The stan path has no equivalent")),
  lbfgs_gll = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 0)),
    msg = paste0("makes the julia optimiser's line search judge a step ",
      "against the worst of its recent values. The stan path has no equivalent")),
  lbfgs_nonmonotone = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 0)),
    msg = paste0("makes the julia optimiser's line search non-monotone. The ",
      "stan path has no equivalent")),
  initial_alpha = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 1)),
    msg = paste0("sets the length of the julia optimiser's first trial step ",
      "(0.1 raw units by default). The stan path has no equivalent")),
  certify = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("switches off the curvature check that certifies a julia fit ",
      "has converged, which needs the exact Hessian the engine differentiates ",
      "out of its own gradient. The stan path has no equivalent")),
  innergaptol = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 0)),
    msg = paste0("stops the julia optimiser once a step is predicted to gain ",
      "less objective than this, which needs the line search's own directional ",
      "derivative. The stan path has no equivalent")),
  gaptol = list(only = 'julia',
    inert = function(v) FALSE,
    msg = paste0("sets how much objective a julia fit may still have available ",
      "before it is continued, and the stan path certifies nothing to apply ",
      "it to. Use optimcontrol$tol for stan's own objective criterion")),
  overshoot = list(only = 'julia',
    inert = function(v) identical(as.character(v)[1L], 'magnitude'),
    msg = paste0("selects which directions the julia engine pulls back after ",
      "a fit to test whether it stopped at a maximum, and the stan path has ",
      "no such probe. Drop it")),
  stallwindow = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 80)),
    msg = paste0("sets how many iterations the julia optimiser may make ",
      "negligible progress over before it stops, which the stan path does ",
      "not do. Drop it")),
  stallfraction = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 1e-2)),
    msg = paste0("sets what share of a julia fit's own progress counts as ",
      "negligible for the stall check, and the stan path has no such check. ",
      "Drop it")),
  stallratio = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 1e-3)),
    msg = paste0("sets how much of its live responsiveness a transform must ",
      "have lost before it counts as the reason a julia fit stopped. The ",
      "stan path has no such check. Drop it")),
  escapesaturated = list(only = 'julia',
    inert = function(v) !isTRUE(v),
    msg = paste0("asks the julia route to refit a saturated fit from a zeroed ",
      "boundary and keep whichever wins, and the stan path has no such ",
      "mechanism. Drop it")),
  escapepin = list(only = 'julia',
    inert = function(v) !isTRUE(v),
    msg = paste0("asks the julia route to hold the escaped coordinates still ",
      "while the rest of the model re-optimises around them before releasing ",
      "them, and the stan path has no such mechanism. Drop it")),
  stallcooldown = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 30)),
    msg = paste0("sets how long the julia stall check waits after finding a ",
      "fit slow rather than stuck, and the stan path has no such check. ",
      "Drop it")),
  stalltighten = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 0.1)),
    msg = paste0("sets how much the julia stall check tightens its bar after ",
      "finding a fit slow rather than stuck. The stan path has no such ",
      "check. Drop it")),
  stalltightenings = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 2)),
    msg = paste0("caps how often the julia stall check tightens its bar, and ",
      "the stan path has no such check. Drop it")),
  gapretries = list(only = 'julia',
    inert = function(v) isTRUE(all(v == 0)),
    msg = paste0("caps how many curvature-based corrections a julia fit may ",
      "take, which the stan path does not do. Drop it")),
  callback = list(only = 'julia',
    inert = function(v) is.null(v),
    msg = paste0("is called by the julia engine while the fit runs, and the ",
      "stan optimizer has no such hook. Use verbose=1 for stan's own ",
      "iteration output")),
  saveEffects = list(only = 'julia',
    inert = function(v) !isTRUE(v),
    msg = paste0("keeps every draw of every random effect from the julia ",
      "sampler. Use ctSubjectPars() on a stan fit for per-subject draws")),
  progress = list(only = 'julia',
    inert = function(v) !isTRUE(v),
    msg = paste0("forces the julia engine's progress line on or off, and the ",
      "stan optimizer does not draw one. Use verbose=1")),
  tipredMissingIncludeOutcome = list(only = 'julia',
    inert = function(v) isTRUE(v),
    msg = paste0("chooses whether the outcome informs julia's imputation of ",
      "missing TI predictors; stan imputes them as parameters of the joint ",
      "posterior, so on stan it always does. Drop it")),

  # -- julia only, and within julia only for intoverpop='laplace'. Prefixed
  # because they tune the *inner* solve -- one Newton iteration per unit per
  # outer evaluation -- and every other stopping name in this list is about the
  # outer optimiser. There is no inner solve to tune on the augmented route,
  # where the random effects are latent states, nor on stan, which has neither.
  laplace_inner_maxiter = list(only = 'julia',
    inert = function(v) FALSE,
    msg = paste0("caps the Newton iterations each unit's random-effect mode ",
      "is solved with, and stan has no such inner solve -- it augments the ",
      "latent state instead. Drop it")),
  laplace_inner_tol = list(only = 'julia',
    inert = function(v) FALSE,
    msg = paste0("is the gradient at which a unit's random-effect mode counts ",
      "as found, and stan has no such inner solve -- it augments the latent ",
      "state instead. Drop it")),
  # The Laplace term's prior floor, julia only: intoverpop='laplace', and
  # sampling with intoverpop=FALSE, whose sampler is placed on the Laplace fit.
  # 'total' used to be accepted anywhere because it was what every fit did;
  # since the default became 'gated' neither value describes a stan fit, so
  # both are refused there.
  laplace_floor = list(only = 'julia',
    inert = function(v) FALSE,
    msg = paste0("chooses how julia's Laplace term treats a unit whose ",
      "likelihood is convex in its random effects, and stan has no Laplace ",
      "term -- it augments the latent state instead. Drop it")),
  # SAEM, julia only: the random effects julia's Laplace route integrates are
  # sampled instead (saem.jl). FALSE describes what stan does.
  saem_proposal = list(only = 'julia',
    inert = function(v) FALSE,
    msg = paste0("chooses the proposal of the julia engine's SAEM phase, and ",
      "stan has no SAEM. Drop it")),
  saem = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("runs the julia engine's SAEM phase, which samples the random ",
      "effects julia's Laplace route integrates; stan augments them as latent ",
      "states instead. Drop it")),
  # The quadrature correction every julia Laplace fit gets. FALSE describes
  # what stan does, so it is accepted there, as `stochastic=FALSE` is on julia.
  laplace_correct = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("corrects a julia Laplace fit by quadrature over its random ",
      "effects, and stan has no Laplace term to correct. Drop it"))
)

# Refuse a control-list name the chosen backend cannot honour, before any data
# preparation happens. Called from ctFit() for both backends.
#
# Most of what this used to refuse now simply works, so what is left is the
# capability split above and a misspelling check. stanoptimis() would report an
# unknown name itself as "unused argument", but only when optimize=TRUE: an
# optimcontrol passed with optimize=FALSE never reaches it at all, so a
# misspelled name was silent on exactly the route with no other feedback.
.ctFitCheckControls <- function(optimcontrol, backend){
  supplied <- names(optimcontrol)
  if(is.null(supplied)) supplied <- character()
  supplied <- supplied[nzchar(supplied)]
  split <- .ctOptimcontrolSplit()

  refused <- character()
  for(nm in intersect(supplied, names(split))){
    entry <- split[[nm]]
    if(identical(entry$only, backend)) next
    inert <- try(entry$inert(optimcontrol[[nm]]), silent=TRUE)
    if(isTRUE(inert)) next
    refused <- c(refused, paste0("optimcontrol$", nm, " ", entry$msg,
      ", or refit with backend='", entry$only, "'."))
  }
  if(length(refused)) stop(paste(refused, collapse='\n'), call.=FALSE)

  onlyhere <- names(split)[vapply(split, function(x) identical(x$only, backend),
    logical(1))]
  known <- c(.ctOptimcontrolShared, .ctOptimcontrolInert, names(split),
    if(identical(backend,'stan')) names(formals(stanoptimis)))
  unknown <- setdiff(supplied, known)
  if(length(unknown)) stop(
    "Unrecognised optimcontrol name(s): ", paste(unknown, collapse=', '),
    ". Both backends take ", paste(sort(.ctOptimcontrolShared), collapse=', '),
    "; backend='", backend, "' also takes ", paste(sort(onlyhere), collapse=', '),
    ".", call.=FALSE)
  invisible(TRUE)
}

# The stan path hands its optimcontrol straight to stanoptimis(), which has no
# `...`, so the julia-only names have to come off first. Each has already been
# checked above and is inert here by definition -- `gradient='adjoint'` is what
# Stan's autodiff does, `progress=FALSE` is what its optimizer draws -- so
# dropping them removes nothing the caller asked for.
.ctOptimcontrolForStan <- function(optimcontrol){
  split <- .ctOptimcontrolSplit()
  drop <- names(split)[vapply(split, function(x) identical(x$only, 'julia'),
    logical(1))]
  optimcontrol[setdiff(names(optimcontrol), drop)]
}

#' Update a fit to new data, or to the current version of ctsem
#'
#' Either to include different data, or because you have upgraded ctsem and
#' the internal data structure has changed. Works on a fit from either backend.
#'
#' The fit is rebuilt from its own call: the arguments it was made with are
#' passed to \code{\link{ctFit}} again, with the model it was given and, unless
#' \code{data} is supplied, the data it was fitted to. Arguments in \code{...}
#' replace the fit's own.
#'
#' @param oldfit fit object to be updated, from \code{\link{ctFit}} with either
#' backend.
#' @param data replacement long format data object. If not supplied, the data
#' the fit was made with.
#' @param recompile whether to force a recompile of the Stan model -- safer but
#' slower and usually unnecessary. Stan fits only.
#' @param refit if TRUE, refits the model using the old estimates as a starting
#' point. Only applicable for optimized fits, not sampling. If FALSE, the old
#' fit is returned with its data and model specification replaced; its
#' estimates, and everything computed from them when it was fitted, are left
#' as they were.
#' @param ... extra arguments to pass to ctFit
#'
#' @return updated fit object, of the same class as \code{oldfit}.
#' @aliases ctStanFitUpdate
#' @export
#'
#' @examples
#' newfit <- ctFitUpdate(ctstantestfit,refit=FALSE)

ctFitUpdate <- function(oldfit, data=NA, recompile=FALSE,refit=FALSE,...){

  julia <- .ctFitIsJulia(oldfit)
  if(julia && isTRUE(recompile)) stop("recompile applies to the compiled ",
    "program of a stan fit, and a julia fit has none. Drop it.", call.=FALSE)
  # The model as the caller wrote it. Not `.ctFitModelObject()`, which is that
  # model after preparation and cannot be prepared a second time.
  basemodel <- .ctFitBaseModel(oldfit)
  if(is.null(basemodel)) stop("This fit does not carry the model ",
    "it was built from, which julia fits made before ctFitUpdate() supported ",
    "them do not, so it cannot be rebuilt. Fit again from its estimate with ",
    "the model you passed to ctFit(): ctFit(datalong, model, inits = ",
    "fit$estimate$raw, backend = 'julia').", call.=FALSE)

  dots <- list(...)
  sampled <- if(julia) !is.null(oldfit$sample) else
    length(oldfit$stanfit$stanfit@sim) > 0
  if(sampled && refit){
    message('A sampled fit is not refitted; updating it with refit=FALSE')
    refit <- FALSE
  }
  if(!refit) message('Trying to do a quick update -- if there are problems, try with refit=TRUE for more robustness')
  # Kept, the estimates are coordinates of the backend that made them.
  if(!refit && !is.null(dots$backend) &&
      !identical(as.character(dots$backend)[1L], if(julia) 'julia' else 'stan'))
    stop("A fit updated with refit=FALSE keeps its estimates, which belong to ",
      "the backend that made them. Use refit=TRUE to fit on another backend.",
      call.=FALSE)

  # `$args$input` -- the literal call, still carrying 'auto'/'maxneeded' and
  # whatever else was unresolved -- not `$args$resolved`, which would freeze
  # this refit at whatever a previous 'auto' happened to route to instead of
  # letting it re-route against the new data or overrides in `...`. A fit made
  # by 3.11.1 or earlier keeps its call in `$args` itself, and reading only
  # `$input` replayed every argument of such a fit at its default -- the
  # priors it was estimated with among them -- when bringing an old fit up to
  # date is what this function is for.
  args <- as.list(if(!is.null(oldfit$args$input)) oldfit$args$input else
    oldfit$args)
  # That capture is the calling environment rather than the literal call, so a
  # defaulted argument is in it and `do.call()` below hands it back looking as
  # though the caller had typed it. For `poprank` that is the difference
  # between a no-op and an error: `ctFit()` reads `!missing(poprank)` as "a
  # rank was asked for" and refuses one on the stan backend, so replaying the
  # default 'auto' made this function fail on every stan fit, the documented
  # example included. The capture recorded the distinction next to the value;
  # use it, before the `...` overrides, so a rank passed here still counts as
  # asked for.
  if(!isTRUE(args$poprankexplicit)) args$poprank <- NULL
  # `priors` is captured after ctFit() has reduced it to a logical, so the julia
  # default, 'randomCorr', reads TRUE there -- and TRUE replayed is a prior on
  # every coordinate, a different estimator. `priorscope` beside it still says
  # which. 'all' and 'none' are what TRUE and FALSE mean, and a stan fit never
  # records 'randomCorr'.
  if(identical(args$priorscope, 'randomCorr')) args$priors <- 'randomCorr'
  # `iter`, `chains` and `control` are captured as what they resolved to, so
  # replayed they read as the deprecated spellings and draw that warning at a
  # caller who never used them. When `sampleControl` alone resolves to the
  # same settings they add nothing, and are left out.
  sampling <- c('iter', 'chains', 'control')
  if(isTRUE(all.equal(.ctSampleControlResolve(args$sampleControl),
      args[sampling]))) args[sampling] <- NULL
  # Only ctFit()'s arguments are replayed. The capture also holds its locals
  # (`datavars`, `priorscope`, `poprankexplicit`), `nopriors`, whose effect
  # `priors` already carries, and arguments since removed -- `vb`, which every
  # 3.11.1 fit carries and ctFit() now refuses by name. Anything else was the
  # caller's own `...`, which only stan's sampler reads, and a fit refitted
  # here was optimised, so nothing ever read it.
  args <- args[intersect(names(args),
    setdiff(names(formals(ctFit)), c('...', 'nopriors')))]
  for(n in names(dots)){
    args[[n]] <- dots[[n]]
  }
  args$fit <- refit
  args$inits <- .ctFitRawEstimate(oldfit)
  args$model <- basemodel
  args$ctstanmodel <- NULL
  # The data the fit was made with, unless replacement data is given. (This was
  # `length(data==1)`, which is TRUE for any data set, so the old data was
  # rebuilt on every call and then discarded whenever new data was supplied.)
  newdata <- !is.null(data) && !(is.atomic(data) && length(data) == 1L &&
    is.na(data))
  args$datalong <- if(newdata) data else .ctFitLongData(oldfit)
  newfit <- do.call(ctFit,args)
  if(refit) return(newfit)

  if(julia){
    # A prepared julia model is its own specification, carrying the call, the
    # prepared data and the model beside it; a fit keeps those apart.
    spec <- unclass(newfit)
    spec[c('args', 'standata', 'modelbase')] <- NULL
    # `nlcontrol$nsubsteps = 'auto'` leaves the fit with a mesh: one substep
    # count per row of the data it was chosen for, at the estimate. The same
    # rows keep it; other rows get the one the fit would have chosen for them
    # at the same estimate.
    if(is.integer(oldfit$model_spec$max_timestep) && !is.null(spec$substeps)){
      if(newdata) spec <- .ctJuliaAutoSubsteps(spec,
        .ctFitRawEstimate(oldfit)[seq_len(.ctBackendNpar(spec))])$spec else
        spec$max_timestep <- oldfit$model_spec$max_timestep
    }
    oldfit$model_spec <- spec
    oldfit$model <- spec$model
  } else {
    if(oldfit$ctstanmodel$recompile || recompile) oldfit$stanmodel <- rstan::stan_model(model_code = newfit$stanmodeltext) else
      oldfit$stanmodel <- stanmodels$ctsm
  }
  oldfit$standata <- newfit$standata
  oldfit$data <- .ctStandataNA(newfit$standata)
  oldfit
}

# The prepared data as a fit reports it in `$data`: `standata` with the 99999
# missing-value sentinel replaced by NA. The `$tipreds` line is a no-op --
# that is not a field of the list, and the time-invariant predictors in
# `$tipredsdata` keep their sentinel in both copies -- and is kept only because
# every fit's `$data` has been built with it.
.ctStandataNA <- function(standata){
  standata$Y[standata$Y==99999] <- NA
  standata$tipreds[standata$tipreds==99999] <- NA
  standata
}

#' @export
ctStanFitUpdate <- ctFitUpdate


# Disable the T0VAR rows and columns that RAWPOPVAR already accounts for.
#
# An individually varying T0MEANS gets no carrier state of its own: that latent
# state *is* the carrier (see `.ctModelIntOverPop`, which skips T0MEANS when
# appending states). So that state's initial covariance is the population
# covariance of its random effect, and RAWPOPVAR is what states it. T0VAR
# stating it as well would be two matrices specifying one quantity, which is
# how the two came to be entangled in the first place. This is the whole of the
# remaining relationship between them.
#
# The values these cells are fixed to do no work on the stan path. That row
# and column of T0cov are dropped where it is assembled, and the entries
# RAWPOPVAR spans are written back there, so neither the diagonal nor the
# off-diagonals reach the fit. The zero covariance between an indvarying
# T0MEANS and a non-indvarying one is the absence of any statement about that
# pair rather than T0VAR stating zero -- which is why the disabling is done at
# assembly and not left to these values surviving a construction.
#
# Level-specific, and it has to be. `indvarying` is the subject level, and
# that is the only flag this reads. A T0MEANS that varies over studies but not
# over subjects keeps its T0VAR: within a study every subject has the same
# T0MEANS, so T0VAR's dispersion really is that subject's initial covariance,
# and the study effect shifts the value for the whole study rather than adding
# to any individual's initial spread. Generalising this to "varies at any
# level" would delete a parameter the data can identify, quietly. There is a
# test for the three cases in test-t0var-redundancies.R.
#
# The composition a reader might expect -- T0VAR plus the subject effects plus
# the study effects -- is the *marginal* initial covariance over all
# individuals, and it is not what T0cov holds. T0cov is conditional on the
# levels above it: T0VAR for the latents that keep it, and the subject level
# population covariance for the rest. The study level reaches an individual by
# moving that individual's parameter values, so adding a study covariance into
# T0cov as well would count it twice.
#
# They are still *fixed* rather than left free, because taking them out of the
# parameter vector is the whole job of this function. Neither backend now reads
# the values: both drop these rows and columns where T0cov is assembled and
# write RAWPOPVAR's entries over them. A small positive diagonal rather than
# zero all the same, because the R-side paths that build a covariance straight
# from T0VAR -- ctGenerate(backend='r'), ctGraph, ctModelLatex -- would meet a
# singular matrix, and they are not what this function is about.
T0VARredundancies <- function(ctm) {
  whichT0VAR_T0MEANSindvarying <- ctm$pars$matrix %in% 'T0VAR'  &
    is.na(ctm$pars$value) &
    (ctm$pars$row %in% ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying] |
        ctm$pars$col %in% ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying])
  if(any(whichT0VAR_T0MEANSindvarying)){
    message('T0VAR rows/columns for latents with individually varying T0MEANS disabled: RAWPOPVAR gives their covariance.')
    ctm$pars$value[whichT0VAR_T0MEANSindvarying & ctm$pars$col == ctm$pars$row ] <- 1e-6
    ctm$pars$value[whichT0VAR_T0MEANSindvarying & ctm$pars$col != ctm$pars$row ] <- 0
    ctm$pars$param[whichT0VAR_T0MEANSindvarying] <- NA
    ctm$pars$transform[whichT0VAR_T0MEANSindvarying] <- NA
    ctm$pars$indvarying[whichT0VAR_T0MEANSindvarying] <- FALSE
    #these cells are no longer parameters, so they carry no TI predictor
    #effects either -- a stale effect here makes ctStanData compute a
    #per subject T0VAR that cannot vary by subject. 'FALSE' as a string:
    #the column holds a character spec, not a logical (R/ctTipredEffect.R),
    #and this has to clear a fixed size or an effect name as well as TRUE.
    if(ctm$n.TIpred > 0) ctm$pars[whichT0VAR_T0MEANSindvarying,
      paste0(ctm$TIpredNames,rep('_effect',ctm$n.TIpred))] <- 'FALSE'
  }
  return(ctm)
}

# A latent whose DRIFT row is fixed at zero never settles anywhere -- a trait
# held constant, or a random walk -- so stationary = TRUE has nothing to start
# it from. Refused by name here, because the engine would otherwise meet a
# singular DRIFT at every trial point and the fit would fail at its starting
# values with nothing said about why.
.ctStationaryCheckDrift <- function(ctm) {
  drift <- ctm$pars[ctm$pars$matrix %in% 'DRIFT', , drop = FALSE]
  fixedzero <- is.na(drift$param) & !is.na(drift$value) & drift$value == 0
  static <- vapply(seq_len(ctm$n.latent), function(i)
    any(drift$row == i) && all(fixedzero[drift$row == i]), logical(1))
  if(any(static)) stop('stationary = TRUE needs every latent process to ',
    'settle somewhere, and DRIFT is fixed at zero for ',
    paste(ctm$latentNames[static], collapse = ', '),
    '. Estimate T0MEANS and T0VAR instead.', call. = FALSE)
  invisible(TRUE)
}



#' Fit a ctsem model
#'
#' Fits a ctsem model specified via \code{\link{ctModel}} with type either 'ct' or 'dt'.
#' \code{ctStanFit} is maintained as a backward-compatible alias.
#'
#' A julia fit with random effects reports, at its estimate, how much of each
#' random effect each subject's own data determine, and how many subjects' worth
#' of information each population standard deviation rests on
#' (\code{fit$identifiability$effects}). An effect resting on fewer than two is
#' named in a message, in \code{summary()} and in \code{\link{ctReport}}, with a
#' lower \code{poprank} suggested where its level can take one and
#' \code{indvarying = FALSE} after it.
#'
#' @aliases ctStanFit
#' @seealso \code{\link{ctOptimControl}} for every optimiser setting;
#' \code{\link{ctFitUncertainty}} for uncertainty after a fit;
#' \code{\link{ctIdentify}} reports which parameters the data can inform,
#' before a fit is spent finding out; \code{\link{ctTracePlot}} draws the
#' optimisation trace a julia fit records.
#' @param datalong long format data containing columns for subject id, time,
#' manifest variables, any time dependent (i.e. varying within subject)
#' predictors, and any time independent (not varying within subject)
#' predictors.
#' @param model model object as generated by \code{\link{ctModel}} with type='ct' or 'dt', for continuous or discrete time
#' models respectively.
#' @param ctstanmodel Deprecated. Use \code{model}.
#' @param stanmodeltext already specified Stan model character string, generally leave NA unless modifying Stan model directly.
#' (Possible after modification of output from fitting with argument fit=FALSE.) \code{backend='stan'} only.
#' @param intoverstates Logical. TRUE (the default) integrates over the latent
#' states with a Kalman filter. FALSE makes the states part of what is
#' estimated -- the target becomes the joint density of the data and the
#' states, and the fit reports the trajectory in \code{fit$estimate$states} --
#' and is for sampling: with \code{backend='julia'} that route makes no
#' Gaussian assumption about the state anywhere, where the filter's update for
#' a binary, ordinal or count indicator is an assumed-density projection.
#' \code{optimize=TRUE} is refused, because the joint mode is degenerate rather
#' than an estimate; \code{optimcontrol$estonly=TRUE} returns it anyway, with
#' no standard errors.
#' @param binomial Deprecated. Logical indicating the use of binary rather than Gaussian data, as with IRT analyses.
#' This now sets the \code{manifesttype} of every indicator to 1, for binary.
#' @param fit If FALSE, prepare the model and data without fitting: stan
#' returns its fit object unfitted, julia the specification its engine reads.
#' @param poprank Rank of the population covariance of the individually
#' varying parameters. \code{'auto'}, the default, uses the number of varying
#' parameters that reach the observation mean: under
#' \code{intoverpop='augmented'} that is all the route can identify, so it
#' costs no likelihood, it does nothing on a model whose varying parameters all
#' reach the mean, and it says so when it reduces. \code{NA} leaves the
#' covariance unrestricted. A smaller whole number is an approximation -- fewer
#' dimensions than the data support, which lowers the likelihood and distorts
#' the parameters retained -- useful for parsimony or speed.
#'
#' Requires \code{backend='julia'}, and the default applies to
#' \code{'augmented'} only: under \code{'laplace'} and \code{'none'} the
#' dropped coordinates are identified, so a rank has to be asked for. There it
#' may be given per grouping level by name, as \code{poprank=c(study=2)}, with
#' \code{'auto'} resolved per level, and it may not fall below the number of
#' individually varying T0MEANS with a free T0VAR, whose latents then take their
#' initial variance from this covariance. Under \code{'augmented'} a varying T0MEANS cannot be
#' reduced at all, and an explicit rank is refused. A restricted covariance is
#' \code{L \%*\% t(L)} for a loading matrix \code{L}, reported as
#' \code{poploading_<parameter>_dim<j>.<level>}; a loading's sign is
#' arbitrary, so read the standard deviations and correlations
#' \code{summary()} reports. A stated \code{RAWPOPVAR} cell turns a defaulted
#' reduction off, and is refused alongside a rank asked for. May also be set on
#' the model, as \code{model$poprank <- 2}.
#' @param intoverpop How declared individual differences (random effects) are
#' handled: integrated out of the likelihood, or sampled with everything else.
#'
#' \code{'auto'}, the default, integrates them out when optimizing, as
#' \code{TRUE} does, and samples them (\code{'none'}) when
#' \code{optimize=FALSE} -- except where TI predictor values are missing and
#' the augmented route is exact for the model, since only it samples missing
#' predictors.
#'
#' \code{TRUE} integrates them out by whichever method suits the model:
#' \code{'augmented'} where every individual difference shifts a mean
#' affinely and every indicator is Gaussian, where it is exact and cheapest,
#' and \code{'laplace'} otherwise -- a varying parameter in DRIFT, DIFFUSION,
#' MANIFESTVAR or LAMBDA, a non-affine transform, a parameter entering another
#' cell's expression, a non-Gaussian indicator, or a grouping level above the
#' subject (see \code{id} in \code{\link{ctModel}}). A choice of
#' \code{'laplace'} is announced, and \code{fit$args$resolved$intoverpopreason}
#' records why. Before 3.12.0 \code{TRUE} always meant \code{'augmented'}.
#'
#' \code{'augmented'} makes each random effect a state of the filter; variation
#' on a DIFFUSION or MANIFESTVAR parameter is then only partly identified, and
#' ctFit warns. \code{'laplace'} integrates each subject's effects by a Laplace
#' approximation at their mode (\code{backend='julia'}): exact when they enter
#' the state mean linearly, an approximation elsewhere, which \code{summary()}
#' says. \code{FALSE}, or \code{'none'}, samples them, so it needs
#' \code{optimize=FALSE}: maximising over every subject's effects would drive
#' the population variance to zero. With \code{optimize=FALSE} on julia,
#' \code{'none'} (and so \code{'auto'}) samples the joint posterior of
#' parameters and random effects; \code{'laplace'} and \code{'augmented'}
#' sample the population parameters with the effects integrated out, by the
#' Laplace approximation or the filter, and \code{TRUE} by whichever of the two
#' it resolves to. \code{sampleControl$target} overrides this.
#' @param sameInitialTimes if TRUE, include an empty observation for every subject that has no observation
#' at the earliest observation time of the dataset. This ensures that the T0MEANS occurs for every subject at the same time,
#' rather than just at the earliest observation for that subject. Important when modelling trends over time, age, etc.
#' @param plot if TRUE, for sampling, a Shiny program is launched upon fitting to interactively plot samples.
#' May struggle with many (e.g., > 5000) parameters. For optimizing, the optimisation trace is plotted --
#' with \code{backend='julia'} once the fit returns, since the fit is one call
#' into the engine; for live output use \code{optimcontrol$callback}.
#' @param derrind deprecated, latents involved in dynamic error calculations are determined automatically now.
#' @param optimize if TRUE, estimate by optimisation -- maximum likelihood, or
#' maximum a posteriori with \code{priors=TRUE} -- with uncertainty from the
#' curvature at the estimate. If FALSE, sample the posterior: with Stan's HMC
#' sampler for \code{backend='stan'}, and for \code{backend='julia'} with
#' SAEM's kernel on the joint posterior of parameters and random effects, or
#' NUTS when there are none or \code{intoverpop} integrates them out -- slower,
#' but exact rather than a normal approximation at a mode. Other uncertainty methods for an optimised fit,
#' importance sampling among them, are in \code{\link{ctFitUncertainty}}. A
#' sampled fit's point estimate is the per-parameter median of the draws on
#' stan (\code{fit$stanfit$rawest}) and their mean on julia
#' (\code{fit$estimate$raw}).
#' @param optimcontrol List of optimiser settings, read by both backends with
#' the same meaning: stopping rules (\code{tol}, \code{g_tol},
#' \code{maxiter}), the start (\code{initsd}, \code{carefulfit}), the
#' uncertainty computed at the end (\code{uncertainty}, \code{estonly}), and
#' optional phases such as \code{stochastic} and, under
#' \code{intoverpop='laplace'}, \code{saem} and \code{laplace_correct}.
#' \code{\link{ctOptimControl}} lists every name with its backend and default.
#' A misspelled name, or one the chosen backend cannot honour, is an error.
#' @param nopriors deprecated, use priors argument. logical. If TRUE, any priors are disabled -- sometimes desirable for optimization.
#' @param iter \strong{Deprecated} -- use \code{sampleControl$iter}. Still
#' honoured, with a warning.
#' @param inits \code{NULL} (the default) for ctsem's own starting values, or
#' a vector of starting values on the raw (unconstrained) scale, one per free
#' parameter, as a fit holds its estimate: \code{fit$estimate$raw} on julia,
#' \code{fit$stanfit$rawest} on stan. With \code{backend='stan'} and
#' \code{optimize=FALSE}, the string \code{'optimize'} optimises first and
#' starts the chains there.
#' @param priors \code{'randomCorr'}, \code{TRUE} or \code{FALSE}. \code{TRUE}
#' adds ctsem's \code{normal(0,1)} raw-scale prior to every parameter, making
#' an optimised fit maximum a posteriori; \code{FALSE} uses none.
#' \code{'randomCorr'}, the default, applies it to the
#' random-effect correlations only, at every level: those are where unbounded
#' maximum likelihood is ill-posed, since a bounded correlation on an
#' unbounded coordinate walks along any direction the data do not determine.
#' Sampling (\code{optimize=FALSE}) without a prior on every parameter warns,
#' because a parameter the data leave flat has an improper posterior; use
#' \code{priors=TRUE} there.
#' \code{'randomCorr'} needs \code{backend='julia'}: asked for on stan it is an
#' error, and left at the default on stan it means \code{FALSE}.
#' @param cores number of cpu cores to use: a positive whole number, or
#' \code{'maxneeded'} for all but one (capped at the number of chains on
#' stan). Defaults to \code{getOption("mc.cores", 2)}. It is the most CPU
#' cores the call uses at once, counting every process it starts: on julia
#' the threads a fit uses, BLAS included, and a sampling run's worker
#' processes, which share it between them. A julia fit at \code{cores > 1} is not
#' reproducible to the last decimal, because the chunk tuner times candidate
#' splits; use \code{cores = 1} for a before-and-after comparison. The Julia
#' session starts as wide as the call that starts it, and a later call asking
#' for more runs at what it has and says so: \code{ctJuliaSetup(threads = n,
#' force = TRUE)} restarts it wider, and
#' \code{options(ctsem.julia.restart = TRUE)} has a fit do so itself.
#' @param backend Either 'stan' (the default) or 'julia'. The julia backend is a
#' separate engine with the same model definitions and the same summaries. It
#' optimises and samples (\code{optimize=FALSE}), and adds
#' \code{intoverpop='laplace'}, grouping levels above the subject, and
#' ordinal, count and censored indicators, and it integrates a binary one rather
#' than linearising it. It
#' needs a working Julia -- see
#' \code{\link{ctJuliaSetup}} and \code{\link{ctJuliaInstall}}.
#' @param sampleControl Used when \code{optimize=FALSE}: a list of sampler
#' settings. \code{chains} (4 on julia, 2 on stan); \code{iter} (1000, warmup
#' and sampling together); \code{warmup} (200, or half of \code{iter} if that
#' is less); \code{draws}, the post-warmup count, which wins over \code{iter}
#' when given; \code{seed}; \code{saveEffects}, whether every random-effect
#' draw is kept; and \code{processes} (TRUE), whether the chains run in worker
#' processes -- at most \code{cores} at once, so four chains on two cores run
#' two at a time.
#'
#' For \code{backend='julia'} also \code{target}, which posterior:
#' \code{'auto'} (the one \code{intoverpop} names), \code{'joint'} (parameters
#' and random effects) or \code{'marginal'} (parameters, effects integrated
#' out); \code{sampler} (\code{'saem'}, SAEM's
#' kernel, on the joint posterior, \code{'nuts'} otherwise or by choice),
#' \code{placement} (\code{'saem'} on the joint posterior: SAEM's state places
#' the chains, so no Laplace optimum is needed; \code{'fit'} places them
#' around an optimised fit), the stopping target \code{minESS} (200, of the
#' worst parameter) and \code{rhatTarget} (1.01), with \code{meanESS},
#' \code{maxDraws} and \code{settleTol}, and NUTS's \code{maxdepth},
#' \code{target_accept}, \code{maxdelta}, \code{init_scale},
#' \code{adapt_metric} and \code{adapt_effects} -- all described under
#' \code{uncertainty = 'sample'} in \code{\link{ctFitUncertainty}}, with how
#' lower targets give a quicker approximate posterior. A run stops once the
#' target is met, and may run to four times the post-warmup count to reach it.
#' A name the sampler does not read is an error. For \code{backend='stan'} the
#' settings are rstan's, passed to \code{\link[rstan]{stan}}'s \code{control}.
#' @param chains \strong{Deprecated} -- use \code{sampleControl$chains}. Still
#' honoured, with a warning.
#' @param control \strong{Deprecated} -- use \code{sampleControl}, which is the
#' same list under a name that says what it controls. Still honoured, with a
#' warning; where both name the same setting, \code{sampleControl} wins.
#' @param nlcontrol List of non-linear control parameters.
#' \code{maxtimestep} must be a positive numeric,  specifying the largest time
#' span covered by the numerical integration. The large default ensures that for each observation time interval,
#' only a single step of exponential integration is used. When \code{maxtimestep} is smaller than the observation time interval,
#' the integration is nested within an Euler like loop.
#' Smaller values may offer greater accuracy, but are slower and not always necessary. Given the exponential integration,
#' linear model elements are fit exactly with only a single step.
#' \code{nsubsteps = 'auto'} (julia backend only) instead measures, at the starting values and again at the
#' optimum, how nonlinear each observation interval is and refines only the intervals that need it
#' (under \code{intoverpop = 'laplace'}, with each subject at its random-effect modes);
#' \code{substeptol} (default 0.01) is the largest acceptable linearisation error as a fraction of the
#' predicted state standard deviation, and \code{maxsubsteps} (default 64) caps an interval.
#' \code{maxtimestep} remains a ceiling on the step. The choice is reported in \code{fit$substeps}.
#' With \code{optimize = FALSE} that optimum is the one that places the sampler, and the chains keep the
#' mesh chosen there.
#' \code{transition = 'euler'} (julia backend, \code{intoverstates = FALSE} only) replaces the exponential
#' step between substeps of the state-explicit path with plain Euler-Maruyama, for a reference that
#' shares no approximation with the filter; it needs a fine \code{maxtimestep}. See also
#' \code{\link{ctParticleLik}}.
#' @param verbose Integer from 0 to 2. 1 reports progress while the model
#'   fits; 2 additionally keeps every progress line rather than overwriting one
#'   in place, and prints more for debugging. Whether overwriting is possible is
#'   detected from where the output is going; set
#'   \code{options(ctsem.progress.overwrite = FALSE)} if that detection is wrong
#'   for your front end -- a Shiny app capturing stdout, for instance -- or
#'   \code{TRUE} to force it on.
#' @param stationary Logical. If TRUE, each subject's latent processes start
#' from the distribution DRIFT, CINT and DIFFUSION settle into -- the mean
#' \code{-DRIFT^-1 CINT} and the asymptotic covariance -- rather than from
#' T0MEANS and T0VAR, which are then not estimated. Suits processes observed
#' long after they began, and saves estimating the initial state. Requires
#' \code{backend='julia'}, a continuous time model, and a DRIFT, CINT and
#' DIFFUSION that are fixed within each subject: no cell may depend on the
#' latent states or on time dependent predictors. Individual differences in
#' those matrices are allowed with \code{intoverpop='laplace'}, which
#' \code{intoverpop='auto'} then chooses; each subject starts from its own
#' stationary distribution. A latent process whose DRIFT row is fixed at zero
#' never settles, and is refused.
#' @param forcerecompile logical. For development purposes.
#' If TRUE, stan model is recompiled, regardless of apparent need for compilation.
#' @param saveCompile if TRUE and compilation is needed / requested, writes the stan model to
#' the parent frame as ctsem.compiled (unless that object already exists and is not from ctsem), to avoid unnecessary recompilation.
#' @param savescores Logical. If TRUE, output from the Kalman filter is saved in output. For datasets with many variables
#' or time points, will increase file size substantially. \code{backend='stan'} only;
#' on julia, \code{\link{ctKalman}} computes it from the fit.
#' @param savesubjectmatrices Logical. If TRUE, subject specific matrices are saved --
#' only relevant when either time dependent predictors or individual differences are
#' used. Can increase memory usage dramatically in large models, and can be computed after fitting using ctExtract
#' or ctSubjectPars. \code{backend='stan'} only.
#' @param saveComplexPars Logical. If TRUE, also save rowwise output of any complex parameters specified,
#' i.e. combinations of parameters, functions and states.
#' @param gendata Logical -- If TRUE, uses provided data for only covariates and a time and missingness structure, and
#' generates random data according to the specified model / priors.
#' Generated data is in the $Ygen subobject after running \code{extract} on the fit object.
#' For datasets with many manifest variables or time points, file size may be large.
#' To generate data based on the posterior of a fitted model, see \code{\link{ctGenerateFromFit}}.
#' \code{backend='stan'} only; julia generates with \code{\link{ctGenerate}}.
#' @param compileArgs List of arguments to pass to \code{\link[rstan]{stan_model}} for compilation of the Stan model.
#' @param ... additional arguments to pass to \code{\link[rstan]{stan}} function
#' (\code{backend='stan'} only; refused on julia).
#' @return A fitted object of class \code{ctStanFit} (\code{backend='stan'}) or
#' \code{ctJuliaFit} (\code{backend='julia'}), both also classed \code{ctFit}.
#' Besides backend-specific components, every fit carries \code{$args}, a list
#' with two sublists that mean the same thing on both backends:
#' \describe{
#'  \item{\code{input}}{Exactly what was passed to \code{ctFit()}, or the
#'  formal default when an argument was not supplied -- \code{intoverpop} may
#'  still read \code{'auto'}, and \code{cores} \code{'maxneeded'} if that is
#'  what was passed.}
#'  \item{\code{resolved}}{What the fit actually ran with, after \code{'auto'}
#'  routing, deprecated-argument merges (e.g. \code{nopriors} into
#'  \code{priors}) and backend-specific defaults were settled -- \code{cores}
#'  is a concrete integer, \code{intoverpop} is one of \code{'augmented'},
#'  \code{'laplace'} or \code{'none'}, and \code{priors} is \code{TRUE},
#'  \code{FALSE} or \code{'randomCorr'}. A call that gives the same \code{input}
#'  on both backends gives the same \code{resolved} on both.}
#' }
#' Before this, \code{$args} itself held the raw call on a stan fit and only
#' the resolved settings on a julia fit, so the same field name meant opposite
#' things across backends; code reading \code{fit$args$intoverpop} or
#' \code{fit$args$cores} directly should now read \code{fit$args$resolved} (or
#' \code{$input}, if the literal call is what is wanted, as
#' \code{\link{ctFitUpdate}} does).
#' @export
#' @examples
#' \donttest{
#'
#' #generate a modern ctsem model relying heavily on defaults
#' model<-ctModel(type='ct',
#'   latentNames=c('eta1','eta2'),
#'   manifestNames=c('Y1','Y2'),
#'   MANIFESTVAR=diag(.1,2),
#'   TDpredNames='TD1',
#'   TIpredNames=c('TI1','TI2','TI3'),
#'   LAMBDA=diag(2))
#'
#' fit<-ctFit(ctstantestdat, model,priors=TRUE)
#'
#' summary(fit)
#'
#' plot(fit,wait=FALSE)
#'
#' #### extended examples
#'
#' library(ctsem)
#' set.seed(3)
#'
#' #  Data generation (run this, but no need to understand!) -----------------
#'
#' Tpoints <- 20
#' nmanifest <- 4
#' nlatent <- 2
#' nsubjects<-20
#'
#' #random effects
#' age <- rnorm(nsubjects) #standardised
#' cint1<-rnorm(nsubjects,2,.3)+age*.5
#' cint2 <- cint1*.5+rnorm(nsubjects,1,.2)+age*.5
#' tdpredeffect <- rnorm(nsubjects,5,.3)+age*.5
#'
#' for(i in 1:nsubjects){
#'   #generating model
#'   gm<-ctModel(Tpoints=Tpoints,n.manifest = nmanifest,n.latent = nlatent,n.TDpred = 1,
#'   type='omx',
#'     LAMBDA = matrix(c(1,0,0,0, 0,1,.8,1.3),nrow=nmanifest,ncol=nlatent),
#'     DRIFT=matrix(c(-.3, .2, 0, -.5),nlatent,nlatent),
#'     TDPREDMEANS=matrix(c(rep(0,Tpoints-10),1,rep(0,9)),ncol=1),
#'     TDPREDEFFECT=matrix(c(tdpredeffect[i],0),nrow=nlatent),
#'     DIFFUSION = matrix(c(1, 0, 0, .5),2,2),
#'     CINT = matrix(c(cint1[i],cint2[i]),ncol=1),
#'     T0VAR=diag(2,nlatent,nlatent),
#'     MANIFESTVAR = diag(.5, nmanifest))
#'
#'   #generate data
#'   newdat <- ctGenerate(ctmodelobj = gm,n = 1,burnin = 2,
#'     dtmat<-rbind(c(rep(.5,8),3,rep(.5,Tpoints-9))))
#'   newdat[,'id'] <- i #set id for each subject
#'   newdat <- cbind(newdat,age[i]) #include time independent predictor
#'   if(i ==1) {
#'     dat <- newdat[1:(Tpoints-10),] #pre intervention data
#'     dat2 <- newdat #including post intervention data
#'   }
#'   if(i > 1) {
#'     dat <- rbind(dat, newdat[1:(Tpoints-10),])
#'     dat2 <- rbind(dat2,newdat)
#'   }
#' }
#' colnames(dat)[ncol(dat)] <- 'age'
#' colnames(dat2)[ncol(dat)] <- 'age'
#'
#'
#' #plot generated data for sanity
#' plot(age)
#' matplot(dat[,gm$manifestNames],type='l',pch=1)
#' plotvar <- 'Y1'
#' plot(dat[dat[,'id']==1,'time'],dat[dat[,'id']==1,plotvar],type='l',
#'   ylim=range(dat[,plotvar],na.rm=TRUE))
#' for(i in 2:nsubjects){
#'   points(dat[dat[,'id']==i,'time'],dat[dat[,'id']==i,plotvar],type='l',col=i)
#' }
#'
#'
#' dat2[,gm$manifestNames][sample(1:length(dat2[,gm$manifestNames]),size = 100)] <- NA
#'
#'
#' #data structure
#' head(dat2)
#'
#'
#' # Model fitting -----------------------------------------------------------
#'
#' ##simple univariate default model
#'
#' m <- ctModel(type = 'ct', manifestNames = c('Y1'), LAMBDA = diag(1))
#' ctModelLatex(m)
#'
#' #Specify univariate linear growth curve
#'
#' m1 <- ctModel(type = 'ct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   DRIFT=matrix(-.0001,nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   T0VAR=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror'),nrow=1,ncol=1))
#'
#' ctModelLatex(m1)
#'
#' #fit
#' f1 <- ctFit(datalong = dat2, model = m1, optimize=TRUE, priors=FALSE)
#'
#' summary(f1)
#'
#' #plots of individual subject models v data
#' ctPredict(f1,plot=TRUE,subjects=1,kalmanvec=c('y','yprior'),timestep=.01)
#' ctPredict(f1,plot=TRUE,subjects=1:3,kalmanvec=c('y','ysmooth'),timestep=.01,errorvec=NA)
#'
#' #compare randomly generated data from the posterior to the observed data,
#' #including the lagged covariance structure and the mean trajectory over time
#' ctPostPredict(f1, wait=FALSE)
#'
#'  ### Further example models
#'
#' #Include intervention
#' m2 <- ctModel(type = 'ct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix(-1e-5,nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix(0,nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror'),nrow=1,ncol=1))
#'
#'
#'
#' #Individual differences in intervention, Bayesian estimation, covariates
#' m2i <- ctModel(type = 'ct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   TIpredNames = 'age',
#'   TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect||TRUE'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix(-1e-5,nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix(0,nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror'),nrow=1,ncol=1))
#'
#'
#' #Including covariate effects
#' m2ic <- ctModel(type = 'ct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   n.TIpred = 1, TIpredNames = 'age',
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix(-1e-5,nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix(0,nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror'),nrow=1,ncol=1))
#'
#' m2ic$pars$indvarying[m2ic$pars$matrix %in% 'TDPREDEFFECT'] <- TRUE
#'
#'
#' #Include deterministic dynamics
#' m3 <- ctModel(type = 'ct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix('drift11',nrow=1,ncol=1),
#'   DIFFUSION=matrix(0,nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix('t0var11',nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c('merror1'),nrow=1,ncol=1))
#'
#'
#'
#'
#'
#' #Add system noise to allow for fluctuations that persist in time
#' m3n <- ctModel(type = 'ct',
#'   manifestNames = c('Y1'), latentNames=c('eta1'),
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix('drift11',nrow=1,ncol=1),
#'   DIFFUSION=matrix('diffusion',nrow=1,ncol=1),
#'   CINT=matrix(c('cint1'),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix('t0var11',nrow=1,ncol=1),
#'   LAMBDA = diag(1),
#'   MANIFESTMEANS=matrix(0,ncol=1),
#'   MANIFESTVAR=matrix(c(0),nrow=1,ncol=1))
#'
#'
#'
#' #include 2nd latent process
#'
#' m4 <- ctModel(n.manifest = 2,n.latent = 2, type = 'ct',
#'   manifestNames = c('Y1','Y2'), latentNames=c('L1','L2'),
#'   n.TDpred=1,TDpredNames = 'TD1',
#'   TDPREDEFFECT=matrix(c('tdpredeffect1','tdpredeffect2'),nrow=2,ncol=1),
#'   DRIFT=matrix(c('drift11','drift21','drift12','drift22'),nrow=2,ncol=2),
#'   DIFFUSION=matrix(c('diffusion11','diffusion21',0,'diffusion22'),nrow=2,ncol=2),
#'   CINT=matrix(c('cint1','cint2'),nrow=2,ncol=1),
#'   T0MEANS=matrix(c('t0m1','t0m2'),nrow=2,ncol=1),
#'   T0VAR=matrix(c('t0var11','t0var21',0,'t0var22'),nrow=2,ncol=2),
#'   LAMBDA = matrix(c(1,0,0,1),nrow=2,ncol=2),
#'   MANIFESTMEANS=matrix(c(0,0),nrow=2,ncol=1),
#'   MANIFESTVAR=matrix(c('merror1',0,0,'merror2'),nrow=2,ncol=2))
#'
#' #dynamic factor model -- fixing CINT to 0 and freeing indicator level intercepts
#'
#' m3df <- ctModel(type = 'ct',
#'   manifestNames = c('Y2','Y3'), latentNames=c('eta1'),
#'   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
#'   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
#'   DRIFT=matrix('drift11',nrow=1,ncol=1),
#'   DIFFUSION=matrix('diffusion',nrow=1,ncol=1),
#'   CINT=matrix(c(0),ncol=1),
#'   T0MEANS=matrix(c('t0m1'),ncol=1),
#'   T0VAR=matrix('t0var11',nrow=1,ncol=1),
#'   LAMBDA = matrix(c(1,'Y3loading'),nrow=2,ncol=1),
#'   MANIFESTMEANS=matrix(c('Y2_int','Y3_int'),nrow=2,ncol=1),
#'   MANIFESTVAR=matrix(c('Y2residual',0,0,'Y3residual'),nrow=2,ncol=2))
#'
#' }

ctFit<-function(datalong, model, stanmodeltext=NA, iter=1000, intoverstates=TRUE, binomial=FALSE,
  fit=TRUE, intoverpop='auto', poprank='auto', sameInitialTimes=FALSE, stationary=FALSE,plot=FALSE,  derrind=NA,
  optimize=TRUE,  optimcontrol=list(),
  backend=c('stan','julia'),
  nlcontrol = list(), nopriors=NA, priors='randomCorr', chains=2,
  cores=getOption("mc.cores", 2L),
  inits=NULL,
  compileArgs=list(),
  forcerecompile=FALSE,saveCompile=TRUE,savescores=FALSE,
  savesubjectmatrices=FALSE, saveComplexPars=FALSE,
  gendata=FALSE,
  sampleControl=list(),
  control=list(),verbose=0,..., ctstanmodel){

  # `vb` (Stan's variational Bayes) was removed: it is a stan-only path,
  # crashed on the default fit=TRUE, and stan is being deprecated. Without
  # this check a caller's `vb=TRUE` would silently fall into `...` and be
  # ignored rather than erroring, which would fit MCMC or MLE while the
  # caller believed they had asked for variational inference.
  # `iter`, `chains` and `control` are entries of `sampleControl` now. Folded
  # here, at the top, so that everything below -- on both backends -- goes on
  # receiving exactly what it received before: the rename cannot change what a
  # fit does, and `control` still reaches rstan unchanged on the stan path.
  # `names(match.call())` because a default cannot otherwise be told from a
  # value that happens to equal it, and only an argument the caller actually
  # wrote should draw a deprecation warning.
  # Whether the caller chose the chain count: the julia sampler defaults to 4,
  # as ctFitUncertainty(fit, 'sample') does, and the formal's 2 is stan's.
  chainsgiven <- !is.null(sampleControl$chains) || !is.null(control$chains) ||
    'chains' %in% names(match.call())
  .ctsample_resolved <- .ctSampleControlResolve(sampleControl,
    given = names(match.call()), iter = iter, chains = chains,
    control = control)
  iter <- .ctsample_resolved$iter
  chains <- .ctsample_resolved$chains
  control <- .ctsample_resolved$control
  # Sampling settings on a fit that does not sample were dropped without a
  # word, and the deprecation warning above even said they took effect.
  sampling_given <- c(if(length(sampleControl)) 'sampleControl',
    intersect(c('iter', 'chains', 'control'), names(match.call())))
  if(isTRUE(optimize) && length(sampling_given)) warning(
    paste(sampling_given, collapse=', '), " set how optimize=FALSE samples; ",
    "this fit optimises, so ", if(length(sampling_given) > 1L) "they are" else
      "it is", " not used.", call.=FALSE)

  if('vb' %in% ...names()) stop(
    "the 'vb' (variational Bayes) argument to ctFit() has been removed -- ",
    "stan's variational inference was broken and stan is being deprecated. ",
    "Use optimize=TRUE or optimize=FALSE (sampling) instead.", call.=FALSE)

  # `backendcontrol` was a second, julia-only control list for settings that
  # `optimcontrol` had no room for; they are all in `optimcontrol` now and mean
  # the same thing on both backends. Caught here rather than left to `...`,
  # where it would be accepted in silence -- which is the fault the merge was
  # for.
  if('backendcontrol' %in% ...names()) stop(
    "backendcontrol has been merged into optimcontrol, which now means the ",
    "same thing on both backends: maxiter, g_tol, f_tol (now tol), x_tol and ",
    "lbfgs_memory are optimcontrol names, and so are gradient, progress and ",
    "tipredMissingIncludeOutcome. julia_project and restart_session were ",
    "already ctJuliaSetup()'s project and force=TRUE.", call.=FALSE)

  if(missing(model)){
    if(missing(ctstanmodel)) stop('model must be supplied')
    warning('ctstanmodel argument is deprecated, use model instead')
    model <- ctstanmodel
  } else if(!missing(ctstanmodel)) {
    stop('Use only one of model or deprecated ctstanmodel')
  }
  ctstanmodel <- model

  # A `ctModel(type='omx')` object is a list of matrices, not a built model: it
  # has no `pars` table and none of the fields read below. The first of them
  # reached is `$timeName`, and a NULL on the left of `&&` gave "invalid
  # argument type" -- naming neither the argument nor what to do. Refuse it
  # here, where both are still known.
  if(inherits(model, 'ctsemInit')) stop(
    'model is a ctModel(type=\'omx\') object, which is a list of matrices ',
    'rather than a model ctFit() can use. Convert it first with ',
    'ctModelConvertOMX(model), or build the model with ctModel(type=\'ct\') ',
    'instead.', call.=FALSE)

  backend <- match.arg(backend)
  # Whether `poprank` was asked for or merely defaulted. Taken here because
  # `missing()` has to be evaluated before the argument is touched, and it
  # decides whether an inapplicable rank is an error or a no-op. A rank stated
  # on the model sets it too, once the model is read (`model$poprank` below).
  poprankexplicit <- !missing(poprank)
  # Before any data preparation and before the Julia install prompt: a control
  # name the chosen backend cannot honour is a mistake to report immediately,
  # not after a wait.
  .ctFitCheckControls(optimcontrol, backend)
  # Both are coerced further down without a check: cores = 0 ran a fit and
  # recorded 0 as what it used, and verbose = 'yes' became NA with only a
  # coercion warning.
  if(!(identical(cores, 'maxneeded') || (is.numeric(cores) &&
      length(cores) == 1L && !is.na(cores) && cores >= 1 && cores == round(cores))))
    stop("cores must be a positive whole number, or 'maxneeded'.", call.=FALSE)
  if(!(is.numeric(verbose) || is.logical(verbose)) || length(verbose) != 1L ||
      is.na(verbose) || verbose < 0 || verbose != round(verbose))
    stop("verbose must be 0, 1 or 2.", call.=FALSE)
  if(backend %in% 'julia') {
    dots <- if(...length()) {
      dotnames <- ...names()
      if(is.null(dotnames)) dotnames <- rep('', ...length())
      ifelse(nzchar(dotnames), dotnames, '(unnamed)')
    } else character()
    .ctJuliaUnsupported(ctstanmodel, optimize=optimize, priors=priors,
      intoverpop=intoverpop, gendata=gendata,
      stanmodeltext=stanmodeltext, compileArgs=compileArgs,
      forcerecompile=forcerecompile, optimcontrol=optimcontrol,
      savescores=savescores, savesubjectmatrices=savesubjectmatrices,
      dots=dots)
    # Before any data preparation, so that a first-time user is asked about the
    # setup they need rather than being told about it after a wait. Only when
    # the fit will actually run: preparation is pure R, and stays usable -- and
    # testable -- on a machine with no Julia at all.
    if(isTRUE(fit)) .ctJuliaEnsureInstalled()
  } else if(!is.null(optimcontrol$is)) {
    # `optimcontrol$is` selected Stan's optimization-plus-importance-sampling
    # route before the `uncertainty` argument replaced it. stanoptimis() has
    # no `is` parameter and no `...` catch-all, so leaving this in place would
    # reach `do.call(stanoptimis, optimcontrol)` and fail with a raw "unused
    # argument" error. Caught here, before any work happens.
    stop("optimcontrol$is is no longer supported; use optimcontrol$uncertainty='is' instead.", call.=FALSE)
  }

  # Before either is touched: missing() is unreliable once an argument is
  # assigned.
  priorsdefaulted <- missing(priors) && is.na(nopriors)
  if(!is.na(nopriors)){
    warning('nopriors argument is deprecated, use priors argument in future')
    priors <- !nopriors
  }

  # `priors` is a logical everywhere below this, and on the stan path, because
  # that is what every existing branch and both `as.integer`/`as.logical`
  # coercions expect. The scope rides alongside it rather than replacing it.
  #
  # 'randomCorr' is the default: a prior on the random-effect correlations and
  # nothing else. Those are the coordinates where unbounded maximum likelihood
  # is ill-posed rather than merely uncertain -- a correlation is bounded, its
  # coordinate is not, and along a direction the data do not determine the
  # coordinate walks until the transform saturates. Measured on a 4x4
  # population covariance with 50 subjects: coordinates at 118, a fit
  # reporting convergence because nothing was moving, and every correlation in
  # the summary `NA`, including the ones the data determine. See
  # `.ctBackendRandomCorrPriorSpec`.
  priorscope <- if(is.character(priors)) match.arg(priors, 'randomCorr') else
    if(isTRUE(priors)) 'all' else 'none'
  if(!is.character(priors) && !is.logical(priors)) stop(
    "priors must be TRUE, FALSE, or 'randomCorr'", call.=FALSE)
  # Sampling needs a prior wherever the likelihood can go flat -- under
  # 'randomCorr' a population sd has none, and on the ?ctFit example the chains
  # wandered off along popsd_mm (R-hat 2.1). The argument is not changed for
  # it: the default is one value whatever the route, and sampling without a
  # prior on every parameter warns below, before the run.
  # The generated Stan model builds its priors in, so a subset of coordinates
  # is not expressible there. Asking for one explicitly is refused by name;
  # arriving at one only because it is the default falls back to stan's own
  # previous behaviour rather than failing a call the user did not make.
  if(priorscope %in% 'randomCorr' && !backend %in% 'julia'){
    if(!missing(priors)) stop("priors='randomCorr' applies a prior to a subset ",
      "of coordinates, which the generated Stan model cannot express. Use ",
      "priors=TRUE or priors=FALSE with backend='stan', or backend='julia'.",
      call.=FALSE)
    priorscope <- 'none'
  }
  priors <- !priorscope %in% 'none'

  if(any(!is.na(derrind))) warning('derrind argment is deprecated, computed automatically now')

  datalong <- data.frame(datalong)

  if(!ctstanmodel$timeName %in% colnames(datalong) && !ctstanmodel$continuoustime) {
    dtable <- data.table(datalong)
    dtable[,.ObsCount:=1:.N,by=ctstanmodel$id]
    datalong[[ctstanmodel$timeName]] <- dtable[['.ObsCount']]
    rm(dtable)
  }

  # Every column the model names, checked before anything reads one, and the
  # missing ones reported together with the role each was named for.
  dataroles <- list(time = ctstanmodel$timeName, id = ctstanmodel$subjectIDname,
    manifest = ctstanmodel$manifestNames,
    `time dependent predictor` = ctstanmodel$TDpredNames,
    `time independent predictor` = ctstanmodel$TIpredNames)
  missingcols <- unlist(lapply(names(dataroles), function(role) {
    absent <- setdiff(dataroles[[role]], colnames(datalong))
    if(length(absent)) paste0(paste(absent, collapse = ', '), ' (', role, ')')
  }))
  if(length(missingcols)) stop(call. = FALSE, 'Columns not found in the data: ',
    paste(missingcols, collapse = '; '), '.')
  # as.numeric() turns text into NA with only a warning, which would drop those
  # observations without a word, so a value that does not read as a number is
  # refused. The check this replaces tested is.numeric(as.numeric(x)), which is
  # always TRUE.
  # Columns that pass, numbers held as text or a logical, are used as numbers.
  for(x in setdiff(unique(unlist(dataroles)), ctstanmodel$subjectIDname)){
    values <- datalong[[x]]
    if(is.numeric(values)) next
    if(is.logical(values)) { datalong[[x]] <- as.numeric(values); next }
    numbers <- suppressWarnings(as.numeric(as.character(values)))
    notnumber <- !is.na(values) & is.na(numbers)
    if(any(notnumber)) stop(call. = FALSE, 'Column ', x, ' holds values that are not numbers: ',
      paste0("'", utils::head(unique(as.character(values[notnumber])), 3), "'", collapse = ', '),
      if(length(unique(values[notnumber])) > 3) ', ...', '.')
    datalong[[x]] <- numbers
  }

  datalong <- datalong[order(datalong[[ctstanmodel$subjectIDname]],datalong[[ctstanmodel$timeName]]),] #sort by subject, time.

  if(!'ctStanModel' %in% class(ctstanmodel)) stop('not a ctStanModel object')

  #set nlcontrol defaults
  if(is.null(nlcontrol$maxtimestep)) nlcontrol$maxtimestep = 999999
  if(is.null(nlcontrol$Jstep)) nlcontrol$Jstep = 1e-6
  # nsubsteps = 'auto' lets the julia engine choose the number of prediction
  # substeps per observation interval from how nonlinear that interval turns
  # out to be (ctsem_auto_substeps in the engine); maxtimestep stays a ceiling
  # on the step. substeptol is the largest acceptable linearisation error, as
  # a fraction of the predicted state standard deviation.
  if(!is.null(nlcontrol$nsubsteps)) {
    if(!identical(nlcontrol$nsubsteps, 'auto')) stop("nlcontrol$nsubsteps must be NULL or 'auto'")
    if(!backend %in% 'julia') stop("nlcontrol$nsubsteps = 'auto' requires backend = 'julia'")
  }
  if(is.null(nlcontrol$substeptol)) nlcontrol$substeptol = 0.01
  if(is.null(nlcontrol$maxsubsteps)) nlcontrol$maxsubsteps = 64
  # transition = 'euler' replaces the exponential (locally linearised) step of
  # the state-explicit path (intoverstates = FALSE) with plain Euler-Maruyama,
  # a reference that shares no approximation with the filter. It needs a fine
  # maxtimestep. No effect on the filter itself.
  if(is.null(nlcontrol$transition)) nlcontrol$transition = 'exponential'
  if(!nlcontrol$transition %in% c('exponential', 'euler')) stop("nlcontrol$transition must be 'exponential' or 'euler'")
  if(nlcontrol$transition == 'euler' && !backend %in% 'julia') stop("nlcontrol$transition = 'euler' requires backend = 'julia'")

  args=c(as.list(environment()), list(...)) #as.list((match.call(expand.dots=FALSE)))
  args$datalong <- NULL
  args$model <- NULL
  args$ctstanmodel <- NULL

  ctm <- ctstanmodel

  # A model saved before the correlation squash moved into the covariance
  # construction still carries it in its own parameter table, so it would be
  # applied twice. Both backends come through here.
  .ctCheckLegacyCovTransformModel(ctm$pars)
  # Covariances given as ctCov(), rewritten for the construction in force now.
  ctm <- .ctCovRefresh(ctm)

  if(!is.null(ctm$TIpredAuto) && ctm$TIpredAuto %in% c(1L,TRUE)){ #if auto tipred, set all effects to true
    for(tip in ctm$TIpredNames){
      ctm$pars[[paste0(tip,'_effect')]] <- 'TRUE'
    }
  }

  # The prior has something to act on only where a level has two or more
  # varying parameters; naming it on a model with none read as a prior that
  # was not there.
  hascorrelations <- any(vapply(.ctVaryingColumns(ctm), function(cl)
    length(.ctVaryingParams(ctm, cl)) > 1L, logical(1)))
  if(optimize && priorscope %in% 'randomCorr') message(
    if(hascorrelations) paste0("Maximum likelihood estimation requested, ",
      "with priors on the random-effect correlations. priors=FALSE for none, ",
      "TRUE for all.") else "Maximum likelihood estimation requested")
  if(optimize && !priors) message("Maximum likelihood estimation requested")
  # `optimcontrol$is` is refused above, so it is NULL by the time we get here
  # and the importance-sampling wording this used to choose is unreachable.
  # Importance sampling is now `optimcontrol$uncertainty='is'`, which runs
  # after optimization rather than instead of it, so the estimation this
  # message describes is a posteriori either way.
  if(optimize && priorscope %in% 'all') message(
    "Maximum a posteriori estimation requested")
  # Naming stan here was wrong for half the fits it described: with
  # backend='julia' the engine runs its own sampler -- SAEM's kernel on the
  # joint posterior by default, NUTS on a marginal one -- and a user reading
  # "Stan's NUTS sampler" on a julia fit has no reason to believe the julia
  # sampler ran at all. Which kernel is said when sampling starts.
  if(!optimize) message("Bayesian estimation via ",
    if(identical(backend,'julia')) "the julia engine's sampler" else
      "Stan's NUTS sampler", " requested")


  ###stationarity
  # The latent processes start from the distribution DRIFT, CINT and DIFFUSION
  # settle into, so T0MEANS and T0VAR are no longer parameters: their cells
  # are fixed here and the engine writes the stationary moments over them
  # (`_ctsem_stationary!`). The fixed values are placeholders nothing reads.
  # Set on the model now rather than with the other flags below, because
  # `.ctIntOverPopAuto()` reads it to keep random effects on the dynamics off
  # the augmented route, where they are states and there is no one system to
  # be stationary in.
  if(!isTRUE(stationary) && !isFALSE(stationary)) stop(
    'stationary must be TRUE or FALSE.', call. = FALSE)
  if(any(ctm$pars$param %in% 'stationary')) stop(
    "Naming a T0MEANS or T0VAR cell 'stationary' is no longer supported; ",
    "use stationary = TRUE, which applies to every latent process.",
    call. = FALSE)
  if(stationary) {
    if(!identical(backend, 'julia')) stop(
      "stationary = TRUE needs backend = 'julia'.", call. = FALSE)
    if(!isTRUE(as.logical(ctm$continuoustime))) stop(
      'stationary = TRUE is implemented for continuous time models only.',
      call. = FALSE)
    .ctStationaryCheckDrift(ctm)
    t0 <- ctm$pars$matrix %in% c('T0MEANS', 'T0VAR')
    ctm$pars$param[t0] <- NA
    ctm$pars$value[t0] <- ifelse(ctm$pars$matrix[t0] %in% 'T0VAR' &
        ctm$pars$row[t0] == ctm$pars$col[t0], 1, 0)
    ctm$pars$transform[t0] <- NA
    ctm$pars$indvarying[t0] <- FALSE
    # 'FALSE' as a string, for the reason T0VARredundancies() gives.
    for(effect in intersect(paste0(ctm$TIpredNames, '_effect'),
      names(ctm$pars))) ctm$pars[t0, effect] <- 'FALSE'
  }
  ctm$stationary <- as.integer(stationary)


  if(length(unique(datalong[,ctm$subjectIDname]))==1 && any(ctm$pars$indvarying[is.na(ctm$pars$value)]==TRUE)){
    # is.null(ctm$fixedrawpopmeans) && is.null(ctm$fixedsubpars) & is.null(ctm$forcemultisubject)) {
    ctm$pars$indvarying <- FALSE
    message('Individual variation not possible as only 1 subject! indvarying set to FALSE on all parameters')
  }

  if(length(unique(datalong[,ctm$subjectIDname]))==1 & any(is.na(ctm$pars$value[ctm$pars$matrix %in% 'T0VAR']))){
    # is.null(ctm$fixedrawpopmeans) & is.null(ctm$fixedsubpars) & is.null(ctm$forcemultisubject)) {
    for(ri in 1:nrow(ctm$pars)){
      if(is.na(ctm$pars$value[ri]) && ctm$pars$matrix[ri] %in% 'T0VAR'){
        ctm$pars$value[ri] <- ifelse(ctm$pars$row[ri] == ctm$pars$col[ri], 1e-3, 0)
      }
    }
    message('Free T0VAR parameters fixed to diagonal matrix of 0.001 as only 1 subject - consider appropriateness!')
  }

  if(binomial){
    # A warning, as the other two deprecations are. A message is easy to miss,
    # and this one silently changes `intoverstates` and every indicator's type.
    warning('binomial argument is deprecated -- set manifesttype in the model object to 1 for binary indicators instead. It has set manifesttype=1 for every indicator.', call.=FALSE)
    # It used to set `intoverstates <- FALSE` as well, which is a leftover from
    # when binary data meant sampling the latent states rather than integrating
    # them. The very next check refuses `intoverstates=FALSE` with
    # `optimize=TRUE` -- so the documented shortcut would now put a user
    # straight into an error, under the default `optimize=TRUE`. Setting `manifesttype` directly never did that, and the
    # linearised measurement handles binary indicators with the filter intact:
    # on a three-indicator model it recovers a generating drift of -0.3 as
    # -0.279 and a diffusion of 0.8 as 0.681, both intervals containing the
    # truth.
    ctm$manifesttype[] <- 1
  }

  recompile <- FALSE
  if(!optimize && priorscope %in% c('none', 'randomCorr')) warning(
    if(priorscope %in% 'none') 'Sampling with priors=FALSE' else
      "Sampling with priors='randomCorr', which leaves the population sds without a prior",
    ": where the data leave a parameter flat its posterior is improper. ",
    "priors=TRUE puts a prior on every parameter.", call. = FALSE)
  # `intoverstates=FALSE` with `optimize=TRUE` maximises the joint density of
  # the parameters and the innovations that build the states, and that mode is
  # degenerate rather than merely biased: with an innovation per observation
  # the states re-optimise to absorb almost any change in the parameters, so
  # the profile is nearly flat -- its largest eigenvalue measured 0.05 on 15
  # subjects x 6 rows, against 22.8 for a well-determined count parameter --
  # and a joint fit reports a better log likelihood and a better conditioned
  # Hessian for changes that mean nothing. Sampling the same density is sound,
  # which is what the route is for, so this is refused and points there.
  #
  # `estonly` is the way through for someone who wants the mode regardless;
  # the warnings below still say what it does. Not for `fit=FALSE`, which
  # optimises nothing and returns the prepared model, whose joint density is a
  # fair thing to evaluate. A property of the estimator rather than of a
  # backend, so it stands for both.
  if(isTRUE(optimize) && !isTRUE(intoverstates) && isTRUE(fit) &&
      !isTRUE(optimcontrol$estonly)) stop(
    'optimize=TRUE with intoverstates=FALSE maximises over the latent states, ',
    'and that joint mode is degenerate rather than an estimate. Use ',
    'optimize=FALSE to sample the states, or intoverstates=TRUE to integrate ',
    'them out. optimcontrol$estonly=TRUE returns the joint mode anyway, ',
    'without standard errors.', call.=FALSE)
  # Maximising over the states rather than integrating them out biases the
  # variance parameters downward -- a variance whose own realisations are
  # being chosen at the same time can always be made to look smaller -- so
  # the joint mode is not the maximum likelihood estimate. That is a property
  # of the estimator and not of a backend, so the warning stands for both.
  # Sampling the same density has no such problem, which is why it points
  # there.
  if(optimize && !intoverstates){
    # Two strengths, because the damage is not uniform. Maximising over the
    # states biases every variance parameter downward, which is the general
    # case. But a parameter that governs how tightly the transition prior
    # constrains the freely chosen states is not merely biased: the joint
    # mode buys cheaper innovations by weakening mean reversion, so such a
    # parameter runs to the flat end of its transform and stays there.
    # Measured on the engine's own fixture: a free DRIFT landed at a raw
    # value of -17 to -19, against a saturation guard at 20, at every sample
    # size from 20 to 120 observations and in six of seven seeds. The
    # engine's count-model test asserts the same runaway independently, on
    # data generated from the model. Fix those parameters and the joint mode
    # is well behaved -- the same fixture with DRIFT fixed converges to a
    # gradient norm of 1e-9 or better -- which is the case this route is for.
    jointrunaway <- ctm$pars$matrix %in% c('DRIFT','DIFFUSION') &
      is.na(ctm$pars$value)
    if(any(jointrunaway)){
      warning(
        'intoverstates=FALSE with optimize=TRUE maximises over the latent ',
        'states rather than integrating them out, and this model leaves ',
        paste0(unique(ctm$pars$matrix[jointrunaway]), collapse=' and '),
        ' free. Those parameters set how tightly the transition prior holds ',
        'the states, so the joint mode weakens them without limit to buy ',
        'cheaper innovations: they run to the end of their transform and ',
        'more data does not help. Treat their estimates as unusable. Fix ',
        'them and free only the means, or use optimize=FALSE to sample the ',
        'states, or intoverstates=TRUE to integrate them out.',
        call.=FALSE)
    } else {
      warning(
        'intoverstates=FALSE maximises over the latent states rather than ',
        'integrating them out, which biases variance parameters downward: ',
        'the joint mode is not the maximum likelihood estimate. Use ',
        'intoverstates=TRUE, or optimize=FALSE to sample the states instead.',
        call.=FALSE)
    }
    # And for a Gaussian indicator with free measurement error it is worse
    # than biased: the joint density is *unbounded*. Send the measurement
    # variance to zero and let the trajectory interpolate the data exactly,
    # and the density diverges -- there is no maximum to find, and an
    # optimiser correctly runs off toward the boundary. Measured on a
    # one-indicator model: MANIFESTVAR reached a raw value of -10.6 with a
    # gradient of 3e9, reported as not converged.
    #
    # Said here rather than left to that non-convergence, which describes
    # the symptom and not the cause. A categorical indicator has no such
    # parameter and is unaffected; so is a Gaussian one whose MANIFESTVAR
    # is fixed. Confirmed degenerate rather than merely biased, so this one
    # is a hard error rather than a warning: there is no fit to return.
    freevar <- ctm$pars$matrix %in% 'MANIFESTVAR' &
      ctm$pars$row == ctm$pars$col & is.na(ctm$pars$value)
    if(any(freevar)){
      gaussian <- ctm$manifesttype[ctm$pars$row[freevar]] %in% 0
      if(any(gaussian)) stop(
        'With intoverstates=FALSE the joint density is unbounded for a ',
        'Gaussian indicator whose MANIFESTVAR is free: the measurement ',
        'variance goes to zero and the latent trajectory interpolates the ',
        'data exactly. The optimiser will run to that boundary and report ',
        'not converged. Fix MANIFESTVAR for ',
        paste(ctm$manifestNames[unique(ctm$pars$row[freevar][gaussian])],
          collapse=', '), ', or use optimize=FALSE.', call.=FALSE)
    }
  }

  # `intoverpop` selects how declared individual differences are handled.
  # TRUE, FALSE and 'auto' keep their existing meanings exactly. The two
  # character forms name the two methods: 'augmented' is the existing
  # state-augmentation, which appends a static latent state per varying
  # parameter and lets the ordinary filter integrate them out; 'laplace'
  # integrates them out per subject instead, so each subject keeps the
  # single-subject state space and the cost of a random effect stops being
  # cubic in the augmented dimension.
  #
  # `intoverpopmethod` carries that choice onwards. `intoverpop` itself stays
  # the logical that drives augmentation everywhere downstream, so the Laplace
  # route reaches the backend with an unaugmented model and every existing
  # `if(intoverpop)` keeps meaning what it meant.
  intoverpopmethod <- 'none'
  # Why 'auto' went the way it did, kept for `fit$args$resolved`; NA when the
  # route was named.
  intoverpopreason <- NA_character_
  # TRUE asks for the random effects to be integrated out, as it always has,
  # and since there are now two ways to do that it takes the one the model
  # suits: the rule 'auto' applies when optimizing, whatever `optimize` is.
  # Before 3.12.0 it meant 'augmented', the only way there was.
  if(isTRUE(intoverpop)){
    auto <- .ctIntOverPopAuto(ctm, backend = backend, optimize = TRUE,
      intoverstates = intoverstates)
    intoverpopmethod <- auto$route
    intoverpopreason <- auto$reason
    intoverpop <- identical(auto$route, 'augmented')
    if(isTRUE(auto$announce)) message("intoverpop=TRUE integrates by 'laplace': ",
      auto$reason, ".")
  } else if(is.character(intoverpop)){
    intoverpop <- match.arg(intoverpop[1], c('auto','augmented','laplace','none'))
    if(intoverpop %in% 'auto'){
      # See `.ctIntOverPopAuto()` for the rule and the measurements behind it.
      # An outer level resolves to 'laplace' silently, as it always has: there
      # is one route for that model, not a choice between two.
      auto <- .ctIntOverPopAuto(ctm, backend = backend, optimize = optimize,
        intoverstates = intoverstates,
        timissing = length(ctm$TIpredNames) > 0L &&
          anyNA(datalong[, ctm$TIpredNames, drop = FALSE]))
      intoverpopmethod <- auto$route
      intoverpopreason <- auto$reason
      intoverpop <- identical(auto$route, 'augmented')
      if(isTRUE(auto$announce)) message("intoverpop='auto' chose '",
        auto$route, "': ", auto$reason, ".")
    } else {
      intoverpopmethod <- intoverpop
      intoverpop <- identical(intoverpopmethod,'augmented')
    }
  }
  intoverpop <- isTRUE(intoverpop)
  if(intoverpop) intoverpopmethod <- 'augmented'

  if(identical(intoverpopmethod,'laplace')){
    # Reached through 'auto' too (a grouping level above the subject), where
    # naming 'laplace' alone blamed a choice the caller never made.
    if(!backend %in% 'julia') stop(
      if(!is.na(intoverpopreason)) paste0("This model needs ",
        "intoverpop='laplace' (", intoverpopreason, "), which ") else
        "intoverpop='laplace' ",
      "requires backend='julia'; the generated Stan model ",
      "does not provide the higher-order derivatives it needs.", call.=FALSE)
    if(!.ctAnyVarying(ctm)) stop(
      "intoverpop='laplace' was requested but no free parameters are marked ",
      "indvarying at any level, so there is nothing to integrate over.",
      call.=FALSE)
  }

  # An outer level on a route that cannot carry it.
  #
  # `.ctModelIntOverPop()` reads `indvarying` and nothing else, so a study
  # effect declared in `indvarying_study` would be dropped without trace and
  # the fit would report a single-level model as though that were what was
  # asked for. Refuse by name instead. There is no reason for a lower level to
  # vary before an upper one does -- a study effect with exchangeable subjects
  # inside it is an ordinary model -- so this is about which route can
  # represent the request, not about which requests are meaningful.
  # A fixed value asked to vary between individuals, which it cannot.
  #
  # `.ctVaryingRows()` treats a cell with a value as fixed whatever its
  # `indvarying` flag says, so neither preparation route augments it and the
  # request simply evaporated -- and RAWPOPVAR used to go on offering a
  # population spread for the parameter anyway. Setting `indvarying` directly
  # on `pars` is how nearly every multilevel model here is written, so this is
  # a reachable mistake rather than a hypothetical one.
  #
  # Named cells only. A number written in a matrix has no parameter name, and
  # the common idiom `model$pars$indvarying <- TRUE` sets the flag on every row
  # including those; naming each fixed LAMBDA and T0VAR cell back at someone
  # who wrote that would be noise. A *named* parameter that also carries a
  # value got there by a deliberate assignment, which is the case worth
  # reporting.
  #
  # A warning rather than an error, matching what the model spec parser does
  # with the same contradiction written in one cell: the value is kept, the
  # individual differences are dropped, and the drop is said out loud.
  fixedvarying <- .ctVaryingColumns(ctm)
  fixedvarying <- fixedvarying[fixedvarying %in% names(ctm$pars)]
  if(length(fixedvarying)){
    flagged <- rep(FALSE, nrow(ctm$pars))
    for(cl in fixedvarying) flagged <- flagged | ctm$pars[[cl]] %in% TRUE
    inert <- flagged & !is.na(ctm$pars$value) & !is.na(ctm$pars$param)
    if(any(inert)) warning(
      'Individual differences were requested for ', 
      paste(unique(ctm$pars$param[inert]), collapse=', '),
      ', which ', if(sum(inert) > 1) 'are' else 'is',
      ' fixed to a value and so cannot vary between individuals. Fitted as ',
      'fixed. Clear the value to estimate ', 
      if(sum(inert) > 1) 'them' else 'it', ' with a population spread.',
      call.=FALSE)
  }

  outervarying <- .ctVaryingParams(ctm, .ctOuterVaryingColumns(ctm))
  if(length(outervarying)){
    named <- paste(outervarying, collapse=', ')
    if(!backend %in% 'julia') stop(
      "Random effects above the subject level (", named, ") are represented by ",
      "backend='julia' only; the generated Stan model has one grouping level.",
      call.=FALSE)
    if(intoverpop || (isTRUE(optimize) && !identical(intoverpopmethod,'laplace'))) stop(
      "Random effects above the subject level (", named, ") are integrated out ",
      "by intoverpop='laplace' only -- the augmented route gives carrier states ",
      "to subject level effects and would drop these.", call.=FALSE)
  }

  # Optimizing without integrating over the population distribution is not a
  # combination that means anything. `intoverpop=FALSE` leaves each subject's
  # random effects as free parameters of the objective, so maximizing it
  # maximizes over the effects as well as over the population parameters, and
  # the population variance it lands on is the one that makes those particular
  # effects most likely -- which is zero, or as near as the data allow. Sampling
  # is what handles that model, which is why `intoverpop='auto'` resolves to
  # FALSE exactly when `optimize` is FALSE.
  #
  # It is refused rather than quietly corrected because both readings of the
  # request are plausible -- integrate them, or sample instead -- and guessing
  # would silently answer a different question. Before this it died inside a
  # transform rendering with `invalid format '%.17g'`, which named neither.
  if(isTRUE(optimize) && !intoverpop && identical(intoverpopmethod,'none') &&
      .ctAnyVarying(ctm)) stop(
    "intoverpop=FALSE with optimize=TRUE leaves each subject's random effects ",
    "as free parameters to be maximized over, which drives the population ",
    "variance to zero rather than estimating it. Use intoverpop='augmented' ",
    "or 'laplace' to integrate them out, or optimize=FALSE to sample.",
    call.=FALSE)

  # if(optimize && !intoverpop && any(ctm$pars$indvarying[is.na(ctm$pars$value)]) &&
  #     is.null(ctm$fixedrawpopchol) && is.null(ctm$fixedsubpars)){
  #   intoverpop <- TRUE
  #   message('Setting intoverpop=TRUE to enable optimization of random effects...')
  # }

  # if(intoverpop==TRUE && !any(ctm$pars$indvarying[is.na(ctm$pars$value)])) {
  #   # message('No individual variation -- disabling intoverpop switch');
  #   intoverpop <- FALSE
  # }

  # Individual variation on a variance cell is only partially identified under
  # the augmented route, and it is not backend specific -- the filter is the
  # same on both, so this sits ahead of the stan / julia dispatch deliberately.
  #
  # The augmented route makes a random effect a static latent state with LAMBDA
  # zero on it. A mean-affecting cell (MANIFESTMEANS, CINT, T0MEANS, DRIFT)
  # reaches the observation mean, so the Kalman update moves that state and its
  # population variance is informed. A DIFFUSION or MANIFESTVAR cell is built
  # from the state's mean only: the observation mean function's Jacobian with
  # respect to it is zero, so the update can never move it, and the data learns
  # about the effect solely through its correlation with states the filter can
  # update. The covariance sd_i x sd_j x corr_ij is then identified and its
  # split into an sd and correlations is a ridge -- bounded below, unbounded
  # above. See review/RANDOMEFFECTS-partial-identification-2026-09-07.md.
  #
  # This is deliberately coarse and will occasionally fire where the parameter
  # IS identified, because a parameter also referenced inside a mean-affecting
  # expression -- a DRIFT string, say -- is identified and has no DRIFT row in
  # `$pars` to reveal it. That false positive is accepted rather than chased
  # with a substring scan of the character cells. The accurate distinction is
  # computable, not pattern-matched: `.ctBackendIdentifiability()` returns each
  # flat direction's loadings, so one helper checking that the covariance
  # functionals are orthogonal to the direction would serve both this pre-fit
  # site and the post-fit `.ctBackendIdentifyWarn()`. That is where an exact
  # version belongs.
  #
  # It must run before `.ctModelIntOverPop()` below, which clears `indvarying`
  # on the cells it rewrites into state references, and MANIFESTVAR rows for
  # non-Gaussian indicators are excluded because they are fixed further down
  # (the `errfix` block) and warning about them would contradict that message.
  if(intoverpop){
    revarpars <- ctm$pars$matrix %in% c('DIFFUSION','MANIFESTVAR') &
      ctm$pars$indvarying & is.na(ctm$pars$value)
    if(any(revarpars)) revarpars <- revarpars &
      !(ctm$pars$matrix %in% 'MANIFESTVAR' &
          ctm$pars$row %in% which(ctm$manifesttype > 0 & ctm$manifesttype != 4))
    if(any(revarpars)) warning(
      "Individual variation on DIFFUSION or MANIFESTVAR is only partially ",
      "identified with intoverpop='augmented': the random effect enters only ",
      "the predicted covariance, so the filter never updates it, and the data ",
      "determines its covariance with the other random effects but not the ",
      "split of that covariance into a standard deviation and correlations. ",
      "Affected: ", paste(unique(ctm$pars$param[revarpars]), collapse=', '),
      ". intoverpop='laplace' identifies these separately.", call.=FALSE)
  }

  ctm <- ctModel0DRIFT(ctm, ctm$continuoustime) #offset 0 drift
  ctm$pars <- ctModelStatesAndPARS(ctm$pars,statenames = ctm$latentNames,tdprednames=ctm$TDpredNames) #replace latent states and PARS with state and PAR[] refs, need this early because we rely on [] detection

  # A reduced-rank population covariance, written as a regression of the
  # remaining random effects on a basis of them. Three steps around the
  # augmentation rather than a second path through it: resolve the basis while
  # `indvarying` still says which parameters vary, clear the flag on the
  # regressed ones so `.ctModelIntOverPop()` gives carrier states to the basis
  # alone, then write the regression expressions once those states exist. The
  # rewrite deliberately lands *before* the second `ctModelStatesAndPARS()`
  # call below, so the new mean and coefficient parameters can be introduced as
  # plain labels and turned into `PARS[r,c]` references by the machinery that
  # already does exactly that. See R/ctPopRegression.R.
  # `poprank` defaults to 'auto', so the places it does not apply have to be
  # inapplicable rather than errors: only a rank the user asked for is refused.
  #
  # One structure, a loading matrix on standardised dimensions, reached two
  # ways because the routes build their population covariance in two places.
  # Under `'laplace'` and `'none'` -- which prepare the same engine
  # specification, one integrating the effects and the other sampling them --
  # it is built in the engine, per level, by one piece of code (`laplacerank`
  # below, read by `.ctJuliaLevelRank()`). Under `'augmented'` the effects are
  # carrier states in the filter, so the rank has to be written into the model
  # around `.ctModelIntOverPop()` (R/ctPopRegression.R); that route has only the
  # subject level, so the rewrite's limit to the innermost level costs it
  # nothing. `'none'` used to take a separate rewrite into regression
  # coordinates -- the form `203e223f` measured as reaching the optimum from 1
  # of 12 matched starts against 8 of 12 for loadings, and which tied levels
  # together through one coefficient.
  #
  # What the restriction *means* differs between them, and only the message says
  # so: on the augmented route it removes coordinates the filter cannot see and
  # costs no likelihood, while under laplace those coordinates are identified
  # and removing them is an approximation. The same restriction, different
  # claims.
  # The rank may be stated on the model instead, as `model$poprank`, which is
  # where it belongs for anyone who thinks of it as part of the specification --
  # `indvarying` is set that way in nearly every multilevel model in the tests,
  # so the idiom is already the house one. Not a `ctModel()` argument: that
  # list goes through `ctModelConvertOMX()` and an unknown field there is a
  # risk for no gain, where a plain assignment onto the returned model works and
  # survives.
  #
  # A rank passed to `ctFit()` wins, because an argument at the call site is the
  # more specific statement of the two; the model's value is used only when the
  # argument was left at its default.
  #
  # A rank stated on the model was asked for just as much as one passed here, so
  # from this point it counts as explicit: applied under laplace and 'none',
  # refused by name where it cannot apply. Before this the flag still read "left
  # at its default", and a laplace fit took the model's rank and then ignored
  # it, running at full rank without a word. The `args` capture above already
  # holds the flag as it was, so `ctFitUpdate()` replays the argument only when
  # one was passed, and the model it refits still carries its own rank.
  if(!poprankexplicit && !is.null(ctm[['poprank']])){
    poprank <- ctm[['poprank']]
    poprankexplicit <- TRUE
  }

  # Under `laplace` a rank restricts each level's population covariance where
  # that covariance is actually built -- in the engine, as a loading matrix --
  # rather than by rewriting the model into a basis and a set of regressions.
  # Two reasons it has to be done that way here and not the augmented way.
  #
  # It reaches every level. The rewrite works from `pars$indvarying`, which is
  # the innermost level alone, so on a burst/subject/study model it can restrict
  # the burst covariance and cannot touch the study one -- and the study level,
  # with the fewest groups, is the one that needs it.
  #
  # And it does not entangle the levels. A regressed effect is written into its
  # cell as `(p + beta * b)`, so if the basis effect `b` also varies at another
  # level then `p` inherits *that* level's deviation of `b` through the same
  # `beta`. One coefficient tying two levels together is not the model anyone
  # asked for, and nothing about the resulting fit would look wrong.
  #
  # A loading matrix per level has neither problem: each level's deviation is
  # `L_level * u_level` with its own `u`, and the levels stay independent.
  #
  # Laplace takes this route alone, whether or not it ends up restricting
  # anything: a request that resolves to full rank at every level is answered
  # by full rank, not by handing the rank on to the augmented machinery below,
  # which reads it as a single number for the innermost level.
  #
  # Not by default, and this is the whole reason a request is distinguished
  # from the default. On the augmented route the coordinates `'auto'` removes
  # cannot be identified, so removing them costs nothing and is a good default.
  # Under laplace they *are* identified, and measured on a 250 x 50 design the
  # same restriction costs 48 log likelihood units and takes the fit to the
  # boundary -- basis sd to zero with the coefficient to -612. A default that
  # does that to a user who chose the route precisely because it identifies
  # these things would be indefensible, so here it has to be asked for. The
  # same holds under 'none', below.
  levelnames <- c(ctm$subjectIDname, ctm$groupIDnames)
  laplacerank <- NULL
  enginerank <- intoverpopmethod %in% c('laplace','none') && poprankexplicit &&
    !(length(poprank)==1 && is.na(poprank))
  if(enginerank && !identical(backend,'julia')) stop(
    "poprank requires backend='julia'.", call.=FALSE)
  if(enginerank && !any(ctm$pars$indvarying[is.na(ctm$pars$value)]) &&
      !.ctAnyVarying(ctm, .ctOuterVaryingColumns(ctm))) stop(
    "poprank restricts the population covariance, so it needs a model with ",
    "individually varying parameters.", call.=FALSE)
  if(enginerank){
    levelcolumns <- stats::setNames(c('indvarying', if(length(ctm$groupIDnames))
      paste0('indvarying_', ctm$groupIDnames)), levelnames)
    asked <- .ctPoprankByLevel(poprank, levelnames)
    laplacerank <- vapply(names(asked), function(level){
      value <- asked[[level]]
      if(!identical(value, 'auto')) return(value)
      # `'auto'` keeps the meaning it has on the augmented route -- how many of
      # this level's effects reach the observation mean -- resolved per level
      # rather than once for the model. A level with none keeps its full
      # covariance rather than a rank of zero, which would be no variation.
      roles <- .ctPopEffectRoles(ctm$pars, column=levelcolumns[[level]])
      counted <- if(nrow(roles)) as.integer(sum(roles$mean)) else NA_integer_
      if(is.na(counted) || counted < 1L) NA_integer_ else counted
    }, integer(1L))
    laplacerank <- laplacerank[!is.na(laplacerank)]
    if(!length(laplacerank)) laplacerank <- NULL
    # An innermost-level varying T0MEANS whose latent has a free T0VAR takes
    # its initial variance from this covariance instead (`T0VARredundancies()`,
    # below, disables those T0VAR cells), so a rank below their number fixes
    # the initial state along a direction: on the ?ctFit example poprank 1 cost
    # 4447 log likelihood units and reported convergence. A fixed T0VAR keeps
    # its own variance, and there a lower rank is an approximation like any
    # other. The augmented route refuses the same case.
    subjectrank <- laplacerank[ctm$subjectIDname]
    t0cells <- ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying %in% TRUE &
      is.na(ctm$pars$value)
    freevar <- ctm$pars$matrix %in% 'T0VAR' & is.na(ctm$pars$value) &
      ctm$pars$row == ctm$pars$col
    t0vary <- unique(as.character(ctm$pars$param[t0cells &
      ctm$pars$row %in% ctm$pars$row[freevar]]))
    if(length(subjectrank) && isTRUE(subjectrank < length(t0vary))) stop(
      "poprank ", subjectrank, " is below the ", length(t0vary),
      " individually varying T0MEANS (", paste(t0vary, collapse=', '), "): ",
      "their latents take their initial variance from the population ",
      "covariance, so a lower rank fixes the initial state along a direction. ",
      "Use poprank=", length(t0vary), " or more, or poprank=NA.", call.=FALSE)
    ctm$laplacerank <- laplacerank
    # Said once, for the levels it reduces: on these routes the coordinates a
    # rank drops are identified, so the restriction is an approximation.
    reduced <- vapply(names(laplacerank), function(level){
      k <- nrow(.ctPopEffectRoles(ctm$pars, column=levelcolumns[[level]]))
      if(k > laplacerank[[level]]) paste0(level, ' ', laplacerank[[level]],
        ' of ', k) else NA_character_
    }, character(1L))
    reduced <- reduced[!is.na(reduced)]
    if(length(reduced)) message('poprank: population covariance reduced to rank ',
      paste(reduced, collapse=', '), " -- an approximation under intoverpop='",
      intoverpopmethod, "'; poprank=NA estimates it in full.")
  }

  popregression <- NULL
  if(intoverpop && !(length(poprank)==1 && is.na(poprank))){
    if(!identical(backend,'julia')){
      if(poprankexplicit) stop("poprank requires backend='julia'.", call.=FALSE)
    } else {
      # The augmented route has the subject level only: the rewrite works from
      # `pars$indvarying`. So a per-level rank may name that level and no
      # other, and an unnamed one is a single number for it.
      asked <- .ctPoprankByLevel(poprank,
        if(is.null(names(poprank))) ctm$subjectIDname else levelnames)
      outer <- setdiff(names(asked), ctm$subjectIDname)
      if(length(outer)) stop("poprank names level ",
        paste(outer, collapse=', '), ", but under intoverpop='augmented' ",
        "a rank restricts the subject level ('", ctm$subjectIDname,
        "') alone. A rank per level needs intoverpop='laplace', or ",
        "sampling (optimize=FALSE).", call.=FALSE)
      poprank <- asked[[ctm$subjectIDname]]
      popregression <- .ctPopRegressionSpec(ctm$pars, poprank,
        explicit=poprankexplicit, model=ctm)
      if(!is.null(popregression)) ctm <- .ctPopRegressionDemote(ctm, popregression)
    }
  }
  if(intoverpop)   ctm <- .ctModelIntOverPop(ctm) #extend system matrices for individual differences
  if(!is.null(popregression)){
    ctm <- .ctPopRegressionRewrite(ctm, popregression)
    popregression <- ctm$popregression
    if(!is.null(popregression)) message(.ctPopRegressionMessage(popregression))
  }

#   #check this *after* replacing PARS references as needed
#   if(any(duplicated(ctm$pars$param[ctm$pars$matrix %in% 'T0MEANS' &
#       !grepl('[',ctm$pars$param,fixed=TRUE) &
#       !is.na(ctm$pars$param)]))) stop(paste0(
#         'Unfortunately, duplicate T0MEANS parameters must be specified via inclusion of additional PARS matrix in ctModel: e.g.,
# ctModel(... #regular model code
# PARS=c("t0mPar||TRUE"), #specify an additional parameter called t0mPar, with random effects
# T0MEANS=c("t0mPar","t0mPar"), #insert this parameter into the T0MEANS matrix as many times as needed.
# ... #regular model code)
# '))

  #jacobian addition
  ctm$jacobian <- try(ctJacobian(ctm))
  if('try-error' %in% class(ctm$jacobian)) ctm$jacobian <- ctJacobian(ctm,simplify=FALSE)
  # ctm$jacobian <- unfoldmats( #replaces matrix references with base parameter
  #   c(listOfMatrices(ctm$pars),ctm$jacobian))
  ctm$jacobian <- ctm$jacobian[names(.ctMatricesList()$jacobian)]
  jl <- ctModelUnlist(ctm$jacobian,names(ctm$jacobian))
  jl <- jl[apply(jl,1,function(x) any(!is.na(x))),] #clean up messy leftovers of NA's
  jl2 <- as.data.frame(rbind(data.table(ctm$pars[1,]),data.table(jl),fill=TRUE))[-1,]


  for(i in 1:nrow(jl2)){ #copy base parameter transforms etc to jacobian when needed
    if(!is.na(jl2$param[i]) && jl2$param[i] %in% ctm$pars$param){
      jl2[i, !colnames(jl2) %in% list('matrix','row','col')] <-
        ctm$pars[which(ctm$pars$param %in% jl2$param[i])[1], !colnames(ctm$pars) %in% list('matrix','row','col')]
    }
  }

  ctm$pars <- rbind(ctm$pars,jl2)

  ctm$pars <- ctModelStatesAndPARS(ctm$pars,statenames = ctm$latentNames,tdprednames=ctm$TDpredNames) #replace any new state and par refs with square bracket refs

  # The transform text as the model states it, keyed by cell, kept before
  # `ctModelTransformsToNum` replaces it with four numbers recovered from it by
  # a grid search. That search scores candidates by squared residual, so it
  # cannot see a constant small enough not to move the residual: every variance
  # diagonal carries a `1e-10` floor -- `1e-10 + 5 * log1p_exp(2 * param)` --
  # and the `round(x, 6)` that follows finishes it off. DRIFT's floor is `1e-06`
  # and survives, which is why only the variances lost theirs.
  #
  # A floorless variance reaches *exactly* zero, since `log1p_exp` is exactly
  # zero once `1 + exp(x)` rounds to one, and a zero variance is then divided
  # by and differentiated through. The julia backend reads this rather than
  # re-rendering the numbers; see `.ctJuliaParameterTable`.
  #
  # Keyed by cell and kept *off* `ctm$pars`, because parts of the Stan pipeline
  # read that frame's columns by position -- carrying it there as a character
  # column turned Stan's `pop_CINT` into NaN.
  if(is.null(ctm$transformtext) && is.character(ctm$pars$transform)){
    # A cell whose transform column already parses as a bare number (the
    # output of an earlier ctModelTransformsToNum() call, which
    # ctEBadjustModel() makes before handing its adjusted model to ctFit --
    # see R/ctEmpiricalBayesFit.R) is not descriptive source text: it is a
    # transform *code*, with the real expression already discarded, and this
    # column is character-typed regardless because it mixes such codes with
    # genuine expression strings. Recording a bare code as `text` fed
    # .ctJuliaParameterTable() a rendered transform of literally "1", with no
    # param[] reference left for the engine's adjoint to find -- see the
    # comment there. NA leaves that cell to the numeric multiplier/meanscale
    # reconstruction, which is what a bare code needs regardless of whether it
    # arrived that way from the user or from an earlier numeric reduction.
    bareCode <- !is.na(suppressWarnings(as.numeric(ctm$pars$transform)))
    text <- as.character(ctm$pars$transform)
    text[bareCode] <- NA_character_
    ctm$transformtext <- data.frame(
      matrix = as.character(ctm$pars$matrix),
      row = as.integer(ctm$pars$row), col = as.integer(ctm$pars$col),
      text = text, stringsAsFactors = FALSE)
  }
  ctm <- ctModelTransformsToNum(ctm)

  ctm$pars <- .ctModelCleanctspec(ctm$pars)

  ctm <- T0VARredundancies(ctm)

  if(!all(ctm$pars$transform[!is.na(suppressWarnings(as.integer(ctm$pars$transform)))] %in% c(0,1,2,3,4))) stop('Unknown transform specified -- integers should be 0 to 4')

  #fix binary manifestvariance

  if(any(ctm$manifesttype %in% 2) && !identical(backend, 'julia')){
    stop('Ordinal manifest variables (manifesttype 2) need backend="julia". ',
      'The stan model has no ordinal measurement, so its thresholds would ',
      'never be read; the julia filter integrates the observation over the ',
      'latent instead. Ordinal variable(s): ',
      paste(ctm$manifestNames[ctm$manifesttype %in% 2], collapse=', '), '.',
      call.=FALSE)
  }
  # Same reason as ordinal, and stated separately because the reason a stan fit
  # would be wrong is different: stan has no Poisson measurement at all, so it
  # would treat a count as a Gaussian observation of the log rate and return
  # numbers that look entirely reasonable.
  if(any(ctm$manifesttype %in% 3) && !identical(backend, 'julia')){
    stop('Count manifest variables (manifesttype 3) need backend="julia". ',
      'The stan model has no count measurement, so the observation would be ',
      'treated as Gaussian; the julia filter integrates it against the ',
      'predicted state under a Poisson log link instead. Count variable(s): ',
      paste(ctm$manifestNames[ctm$manifesttype %in% 3], collapse=', '), '.',
      call.=FALSE)
  }
  # Censored is refused on stan for the same reason as the others: stan has no
  # censored measurement, so a value pinned at a limit would be treated as an
  # ordinary observation of that number and the censoring simply ignored.
  if(any(ctm$manifesttype %in% 4) && !identical(backend, 'julia')){
    stop('Censored manifest variables (manifesttype 4) need backend="julia". ',
      'The stan model has no censored measurement, so a value at its limit ',
      'would be read as an ordinary observation; the julia filter integrates ',
      'the censored likelihood against the predicted state instead. Censored ',
      'variable(s): ',
      paste(ctm$manifestNames[ctm$manifesttype %in% 4], collapse=', '), '.',
      call.=FALSE)
  }
  if(any(ctm$manifesttype %in% 1)) .ctDataBinary(datalong, ctm)
  if(any(ctm$manifesttype %in% 2)) .ctDataCategories(datalong, ctm)
  if(any(ctm$manifesttype %in% 3)) .ctDataCounts(datalong, ctm)
  if(any(ctm$manifesttype %in% 4)) .ctDataCensored(datalong, ctm)

  if(any(ctm$manifesttype > 0)){ #if any non continuous variables, (with free parameters)...
    # Binary and ordinal only. For those two a normal term on the linear
    # predictor is not identified rather than merely unwanted: with a probit
    # link E_e Phi(lambda'x + mu + e) is exactly Phi((lambda'x + mu)/sqrt(1+v)),
    # so it is absorbed into the loadings and thresholds, and the logistic link
    # differs only in that the absorption is approximate. Fixing it is the only
    # coherent thing to do.
    #
    # Censored is excluded because that entry is the scale of the Gaussian
    # inside the limits. Count is excluded because the Poisson's variance is
    # locked to its mean, so the term is not absorbable anywhere: it shifts the
    # mean by v/2, which MANIFESTMEANS takes, and multiplies the variance by a
    # factor nothing else in the model can produce. It is the only parameter
    # that moves a count's variance-to-mean ratio, and leaving it out does not
    # make the data equidispersed -- it makes DIFFUSION absorb the difference,
    # which is the parameter the model is usually for.
    deterministic <- which(ctm$manifesttype %in% c(1L, 2L))
    errfix <- which(ctm$pars$matrix %in% 'MANIFESTVAR' &
        (ctm$pars$row %in% deterministic |
            ctm$pars$col %in% deterministic) &
        is.na(suppressWarnings(as.numeric(
          ctm$pars$value))))

    if(length(errfix) > 0){
      message('MANIFESTVAR fixed to 0 for binary / ordinal indicators: the ',
        'measurement link already carries the randomness.')
      ctm$pars$value[errfix] <- 1e-5
      ctm$pars[errfix,c('param','transform','multiplier','offset','meanscale','inneroffset','sdscale')] <- NA
      ctm$pars$indvarying[errfix] <- FALSE
    }

    # A count's dispersion is scalar: each count row is updated on its own,
    # sequentially and conditionally independent of the others given the state,
    # so there is no off-diagonal for it to be correlated through and the
    # engine reads only the diagonal. A free off-diagonal would therefore be
    # accepted here and ignored there, which is the failure this fixes rather
    # than reports -- the cell has one right value and it is zero.
    countrows <- which(ctm$manifesttype %in% 3L)
    offfix <- which(ctm$pars$matrix %in% 'MANIFESTVAR' &
        ctm$pars$row != ctm$pars$col &
        (ctm$pars$row %in% countrows | ctm$pars$col %in% countrows) &
        is.na(suppressWarnings(as.numeric(ctm$pars$value))))
    if(length(offfix) > 0){
      message('Fixing free off-diagonal MANIFESTVAR parameters for count ',
        'indicators to zero -- a count is scored one row at a time, so it has ',
        'no correlated measurement error.')
      ctm$pars$value[offfix] <- 0
      ctm$pars[offfix,c('param','transform','multiplier','offset','meanscale','inneroffset','sdscale')] <- NA
      ctm$pars$indvarying[offfix] <- FALSE
    }

    # A *fixed* non-zero variance on a binary indicator is left alone by the
    # block above, and it is almost certainly a mistake: the Bernoulli link
    # already supplies the randomness, so an extra measurement variance adds
    # noise on top of noise. Free ones are fixed silently because there is a
    # right answer; a deliberate value is the user's, so it is questioned
    # rather than overwritten.
    # What the approximation costs, said once, because it is invisible
    # otherwise and it is large.
    #
    # A binary observation is handled by a moment-matched Gaussian update: the
    # predicted probability from `inv_logit`, a finite-difference Jacobian, and
    # `ycov = Jy etacov Jy' + p(1-p)`. The covariance update that follows,
    # `etacov -= K Jy etacov`, is the standard EKF one and ignores the link's
    # curvature, so posterior uncertainty is understated and the understatement
    # compounds over time steps. The filter ends up believing the latent is
    # pinned down, its innovations shrink, and less process noise is needed to
    # explain them.
    #
    # Measured on a one-latent process with true DIFFUSIONcov 0.16, 50
    # subjects, 50 timepoints, observed only through binary indicators:
    #
    #   5 indicators   0.114   (0.71 of truth)
    #  10 indicators   0.107   (0.67)
    #  30 indicators   0.073   (0.46)
    #
    # DRIFT and CINT survive as long as the latent has few binary indicators;
    # the bias gets *worse* with more of them, which is the signature of a
    # systematic error rather than sampling noise -- more data makes it more
    # confidently wrong. A latent that also has a continuous indicator is
    # unaffected in the same fit.
    #
    # The sharpest evidence is `test-ctRaschExampleTest.R`, which fits the same
    # data twice -- linearised, and sampling the states exactly by HMC -- and
    # compares. The item parameters agree to within 0.012. The dynamics do not:
    #
    #   diff  (process noise)   0.254  linearised   0.373  exact
    #   drift_eta1             -0.095  linearised  -0.935  exact
    #
    # A drift of -0.095 against -0.935 is not a small bias: the linearised fit
    # describes a near random walk where the exact one finds strong mean
    # reversion. That test has been failing, and it is this.
    # Backend specific, because they no longer do the same thing. The julia
    # engine integrates the Bernoulli observation against the predicted state
    # (a mode-centred Gauss-Hermite rule on the scalar linear predictor); the
    # stan model moment-matches it to a Gaussian and linearises. On 40
    # replications with 50 subjects and 12 timepoints, RMSE for a true
    # DIFFUSION of 0.8:
    #
    #   indicators   julia   stan
    #            3   0.134   0.449
    #           10   0.075   0.297
    #           30   0.045   0.229
    #
    # Paired on the same data, every cell favours julia (p <= 0.014). Note the
    # shape of it: stan's *bias* is modest, its RMSE is three to five times
    # worse -- the linearised estimate is unstable rather than uniformly low,
    # which is why single-dataset comparisons looked so different from each
    # other.
    # Nothing is said for the julia backend. That it integrates rather than
    # linearises is a property of the engine the user cannot act on, and the
    # threshold parameterisation is now reported readably by summary()'s
    # Ordinal thresholds section rather than explained at every fit. The stan
    # warning below stays because it is a caution about estimates being
    # unreliable, which is actionable.
    if(!backend %in% 'julia'){
      message('Binary indicators use a linearised (moment-matched Gaussian) ',
        'measurement update on the stan backend, which makes DRIFT and ',
        'especially DIFFUSION unreliable for a latent seen only through them. ',
        'Over 40 replications with a true diffusion of 0.8, RMSE was 0.45 with ',
        '3 indicators and 0.23 with 30, against 0.13 and 0.04 for ',
        "backend='julia', which integrates the observation instead. Prefer ",
        'the julia backend for binary data, or treat process noise as ',
        'indicative.')
    }

    # `deterministic` rather than every non-Gaussian type, so that censored is
    # excluded here too: its MANIFESTVAR entry is the standard deviation of the
    # Gaussian inside the limits, so a fixed non-zero value there is the model
    # working as specified rather than a mistake to warn about.
    linkrows <- which(ctm$pars$matrix %in% 'MANIFESTVAR' &
        ctm$pars$row %in% deterministic &
        ctm$pars$row == ctm$pars$col)
    stated <- linkrows[!is.na(ctm$pars$value[linkrows]) &
        abs(ctm$pars$value[linkrows]) > 1e-4]
    if(length(stated)){
      warning('MANIFESTVAR is fixed to a non-zero value for indicator',
        if(length(stated) > 1) 's ' else ' ',
        paste(ctm$manifestNames[ctm$pars$row[stated]], collapse=', '),
        '. A binary or ordinal indicator gets its randomness from its ',
        'measurement link -- the Bernoulli link for binary, the cumulative ',
        'logit for ordinal -- so this adds measurement noise on top of it, and ',
        'is not separately identified from the loadings and thresholds. Set it ',
        'to 0 unless that is meant.', call.=FALSE)
    }}

  ctm$modelmats <- .ctModelMatSetup(ctm) #slow!
  ctm <- .ctCalcsList(ctm,save=saveComplexPars) #get extra calculations and adjust model spec as needed???

  #store values in ctm
  ctm$intoverpop <- as.integer(intoverpop)
  ctm$nlatentpop <- as.integer(ifelse(ctm$intoverpop ==1, max(ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS']),  ctm$n.latent))
  ctm$intoverstates <- as.integer(intoverstates)
  ctm$priors <- as.integer(priors)
  ctm$stationary <- as.integer(stationary)
  ctm$nlcontrol <- nlcontrol



  #recompile checks
  if(forcerecompile) recompile <- TRUE
  if(naf(!is.na(ctm$rawpopsdbaselowerbound))) recompile <- TRUE
  if(ctm$rawpopsdbase != 'normal(0,1)') recompile <- TRUE
  if(ctm$rawpopsdtransform != 'log1p_exp(2*rawpopsdbase-1) .* sdscale') recompile <- TRUE
  if(any(ctm$modelmats$matsetup[,'transform'] < -10)) recompile <- TRUE #if custom transforms needed

  ncalcsNoJ<- length(unlist(ctm$modelmats$calcs)[!grepl('JAx[',unlist(ctm$modelmats$calcs),fixed=TRUE)])
  if(ncalcsNoJ > 0) recompile <- TRUE

  # Fires when the only calcs are JAx entries: those no longer force a fresh
  # program, so the model reuses the precompiled one with finite difference
  # jacobians. Tested against the TOTAL calc count -- testing ncalcsNoJ, as this
  # did from the commit that introduced it, is the same condition that has just
  # set recompile, so the message could never appear.
  #
  # Stan only, and said so because it was not: `JAxfinite` and `Jyfinite` are
  # read by `ctData.R` into standata and from there by the generated Stan
  # program, and nothing under `inst/julia` reads either. The julia engine
  # differentiates its own model text, so it has analytic jacobians whatever
  # `recompile` says, and `forcerecompile=TRUE` is advice it cannot act on --
  # a message offering a fix for a problem the reader does not have.
  if(!recompile && length(unlist(ctm$modelmats$calcs)) > 0 &&
      !identical(backend, 'julia'))
    message('Finite difference jacobian used to avoid recompiling -- use forcerecompile=TRUE for analytic jacobians')


  #further model adjustments conditional on recompile
  if(!recompile){ #then use finite diffs for some elements

    #collect row and column of complicated jacobian elements into vector
    ctm$JAxfinite <- array(as.integer(unique(
      unlist(ctm$modelmats$matsetup[ctm$modelmats$matsetup$matrix %in% 52 &
          ctm$modelmats$matsetup$when == -999, 'col']))))# &
    # ctm$modelmats$matsetup$copyrow < 1,c('row','col')]))))

    ctm$Jyfinite <- array(as.integer(unique(
      unlist(ctm$modelmats$matsetup[ctm$modelmats$matsetup$matrix %in% 54 &
          ctm$modelmats$matsetup$when == -999, 'col']))))

  }
  if(recompile) ctm$JAxfinite <- ctm$Jyfinite <- array(as.integer(c()))
  ctm$recompile <- recompile


  standata <- .ctPrepareData(ctm,datalong,optimize=optimize, sameInitialTimes=sameInitialTimes) #bit slow
  standata$verbose=as.integer(verbose)
  standata$savesubjectmatrices=as.integer(savesubjectmatrices)
  standata$gendata=as.integer(gendata)

  if(standata$savesubjectmatrices==1L) savescores = TRUE
  standata$savescores=as.integer(savescores)

  # Reached only when the caller asks for 'maxneeded' by name. It is not a
  # default and must not become one: it resolves to the whole machine, and CRAN
  # allows at most two cores unless the user has asked for more. It was the
  # `optimize=FALSE` default until 3.12.0.
  #
  # Resolved here rather than after the julia branch below, which returns before
  # ever reaching the old resolution site: the string arrived at
  # .ctFitJuliaBackend(), `as.integer()` made it NA, and the NA guard there
  # turned it into 1 -- so asking for every core silently got one.
  #
  # `chains` does not cap the julia backend -- it has no MCMC chains. Its
  # parallelism is over subject chunks, and `ctsem_tune_chunks!` measures within
  # whatever ceiling it is given, so the ceiling is simply the machine.
  if(is.character(cores) && cores=='maxneeded') {
    cores <- if(backend %in% 'julia') max(1, parallel::detectCores()-1) else
      max(1,min(c(chains,parallel::detectCores()-1)))
  }

  # `args` (captured above, at the point nothing had been routed yet) is what
  # the caller literally passed, or the formal default -- 'auto', 'maxneeded'
  # and all. `argsresolved` is what the fit actually runs with, once 'auto'
  # routing and backend-independent defaults are settled. Before this, a stan
  # fit's `$args` was the raw call and a julia fit's `$args` was a curated
  # resolved list, so `fit$args$intoverpop` and `fit$args$cores` meant opposite
  # things across backends for the same call (review/jobs/J10). Both sublists
  # are attached below, under `$args$input` and `$args$resolved`, and mean the
  # same thing on both backends -- the fields overridden here are exactly the
  # ones that get routed or cast between capture and use.
  argsresolved <- args
  argsresolved$backend <- backend
  argsresolved$cores <- cores
  argsresolved$intoverpop <- intoverpopmethod
  argsresolved$intoverpopreason <- intoverpopreason
  # The rank actually used, not the argument: 'auto' resolves to a number, and a
  # model the restriction did not apply to reports NA whatever was asked for.
  argsresolved$poprank <- if(!is.null(laplacerank)) laplacerank else
    if(is.null(popregression)) NA_integer_ else as.integer(popregression$rank)
  # The scope, not the logical `priors` was reduced to: the julia default reads
  # TRUE there, and TRUE is a prior on every coordinate, a different estimator.
  argsresolved$priors <- switch(priorscope, all = TRUE, none = FALSE,
    randomCorr = 'randomCorr')
  argsresolved$optimize <- isTRUE(optimize)
  argsresolved$intoverstates <- isTRUE(intoverstates)

  if(backend %in% 'julia') {
    # The resolved route rather than the logical `intoverpop`, which cannot
    # say 'laplace': 'auto' can resolve there, and the refusal of
    # intoverstates=FALSE with laplace reads this argument.
    .ctJuliaUnsupported(ctm, optimize=optimize, priors=priors,
      intoverpop=if(identical(intoverpopmethod, 'laplace')) 'laplace' else
        intoverpop, gendata=gendata,
      stanmodeltext=stanmodeltext, compileArgs=compileArgs,
      forcerecompile=forcerecompile, intoverstates=intoverstates,
      optimcontrol=optimcontrol)
    # `optimize` and `intoverpop` are orthogonal here. `intoverpop` says which
    # random effects are integrated out and how; `optimize` says whether the
    # remaining parameters are maximised or sampled. Every combination is
    # meaningful for this backend:
    #
    #   optimize  route          what runs                 sampled dimension
    #   TRUE      'laplace'      Laplace ML                --
    #   TRUE      'augmented'    augmented ML              --
    #   FALSE     'laplace'      NUTS, Laplace marginal    npar
    #   FALSE     'augmented'    NUTS, filter marginal     npar
    #   FALSE     'none'         sampled, joint            npar + effects
    #
    # The route is what `intoverpop` resolved to: TRUE takes 'laplace' or
    # 'augmented' as the model suits, and 'auto' under optimize=FALSE takes
    # 'none'. The joint row is drawn by SAEM's kernel and placed by SAEM's
    # state unless `sampleControl` says otherwise (`sampler`, `placement`);
    # `sampleControl$target` overrides the target, so 'joint' on 'laplace'
    # samples the joint posterior with the Laplace fit placing the chains.
    #
    # Only 'laplace' is julia-only; the guard for that is above. `'none'` is
    # the third route the engine needs, and it prepares the Laplace structure
    # -- which is what says *which* parameters vary -- without integrating.
    juliaintoverpop <- if(identical(intoverpopmethod,'laplace')) 'laplace' else
      if(intoverpop) 'augmented' else
        if(!optimize && .ctAnyVarying(ctm)) 'none' else
          'augmented'
    # Which posterior a sampled fit draws, settled here because only here is
    # it still known whether `intoverpop` named the route or 'auto' chose it
    # (`.ctBackendSampleMarginal()`).
    if(!isTRUE(optimize)) control$target <- if(.ctBackendSampleMarginal(
      control$target, juliaintoverpop, requested = args[['intoverpop']]))
      'marginal' else 'joint'
    juliafit <- .ctFitJuliaBackend(datalong=datalong, model=ctm, prepared_data=standata, inits=inits,
      cores=cores, optimcontrol=optimcontrol,
      verbose=verbose, fit=fit, priors=priors, priorscope=priorscope,
      optimize=optimize,
      chains=if(chainsgiven) chains else 4L, iter=iter, control=control,
      intoverpop=juliaintoverpop, intoverstates=intoverstates)
    # Replaces whatever narrower `$args` the julia backend built internally
    # (it only ever had the resolved settings, and not all of them) with the
    # two-sublist form every fit now carries -- see the `argsresolved` comment
    # above. This also reaches the `fit=FALSE` case: `.ctFitJuliaBackend()`
    # returns the prepared model spec then, unclassed by `$args` before, and
    # assigning a list element to it here does not disturb its class.
    juliafit$args <- list(input = args, resolved = argsresolved)
    # The model as the caller wrote it, which a stan fit carries as
    # `$ctstanmodelbase`. `$model` is not that: it is what the engine ran,
    # after `.ctModelIntOverPop()` and the rest of the preparation above, and
    # handing it back to ctFit() prepares it a second time -- which errors on
    # the augmented route. ctFitUpdate() rebuilds a fit from this one. Read it
    # through `.ctFitBaseModel()`, which says why the name differs.
    juliafit$modelbase <- ctstanmodel
    # `$data`/`$standata` mean the same thing on both backends: `$standata` is
    # the prepared data with the 99999 missing-value sentinel intact, `$data`
    # is the same thing with that sentinel replaced by `NA` (see
    # `.ctStandataNA()`). `standata` was already computed above,
    # unconditionally, before backend dispatch -- .ctPrepareData() runs for julia
    # too, purely to prepare `prepared_data` for `.ctFitJuliaBackend()` -- so
    # attaching it here costs nothing further and is not a second computation.
    juliafit$standata <- standata
    # Not on a prepared model (`fit = FALSE`). That object *is* the
    # specification, and its own `$data` is the long data frame every julia
    # accessor reads through `.ctBackendSpec()`; replacing it with the prepared
    # list left ctKalman() on a subset of subjects finding no subjects at all.
    if(isTRUE(fit)) juliafit$data <- .ctStandataNA(standata)
    # `plot` draws the trace *after* the fit here, not during it.
    #
    # The Stan path can plot live because it writes sample files a second
    # process reads. A julia fit is one blocking call into the engine, so R
    # cannot draw anything until it returns -- there is no point in this
    # function where a live plot could be made. What is drawn is the same
    # information, recorded every iteration and handed back on the fit; for
    # genuinely live output, `optimcontrol$callback` is called while the fit
    # runs and can draw whatever it likes.
    if(isTRUE(fit) && !identical(plot, FALSE) && !is.null(juliafit$optim$trace)) {
      try(ctTracePlot(juliafit), silent=TRUE)
    }
    return(juliafit)
  }

  # print(standata$savesubjectmatrices)

  #####post model / data checks
  # (`cores='maxneeded'` is resolved above, before the julia backend returns.)

  if(is.logical(stanmodeltext)) {
    stanmodeltext<- ctStanModelWriter(ctm, gendata, ctm$modelmats$extratforms,ctm$modelmats$matsetup)
  }




  if(fit){
    # if(gendata && stanmodels$ctsmgen@model_code != stanmodeltext) recompile <- TRUE
    # if(!gendata && paste0(stanmodels$ctsm@model_code) != paste0(stanmodeltext)) recompile <- TRUE

    # STAN_NUM_THREADS <- Sys.getenv('STAN_NUM_THREADS',unset=NA)
    # Sys.setenv(STAN_NUM_THREADS=cores)

    if(recompile || forcerecompile) {
      message('Compiling model...')

      #r4.2 / windows / old rstan check
      if(.Platform$OS.type=="windows" && R.version$major %in% 4 && as.numeric(R.version$minor) >= 2 &&
          unlist(utils::packageVersion('rstan'))[2] < 25){
        stop('
*****
Compiling not possible with R version 4.2+ on Windows with Rstan version < 2.26
To upgrade rstan, close and restart all R sessions then run:

install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
*****')
        # a=readline('Attempt to upgrade Rstan from Stan repository? y/n')
        # if(a %in% c('yes','y','Y','Yes','YES')){
        #   message(
        #   install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
        #   install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
        # }
      }
      compileArgs$model_name='ctsem'
      compileArgs$model_code=stanmodeltext
      compileArgs$auto_write=TRUE
      sm <- do.call(stan_model,compileArgs)
      # ,allow_undefined = TRUE,verbose=TRUE,
      # includes = paste0(
      #   '\n#include "', file.path(getwd(), 'syl2.hpp'),'"',
      #   '\n')
      if(saveCompile){
        if(exists(x = 'ctsem.compiled',envir= parent.frame())
          && !'stanmodel' %in% class(get('ctsem.compiled',envir = parent.frame()))){
          warning('ctsem.compiled object already exists, not saving compile')
        } else  assign(x = 'ctsem.compiled',sm,envir = parent.frame())
      }
    }
    if(!recompile && !forcerecompile) {
      if(!gendata) sm <- stanmodels$ctsm else sm <- stanmodels$ctsmgen
    }


    # configure inits ---------------------------------------------------------
    initOptim <- (length(inits)==1 && (inits=='optimize' || inits=='Optimize' || inits=='optimise' || inits=='Optimise'))
    if(optimize || initOptim){ #then set optimcontrol list
      # The julia-only names come off here: they were checked at the top of
      # ctFit() and are inert on stan, and stanoptimis() has no `...`.
      optimcontrol <- .ctOptimcontrolForStan(optimcontrol)
      optimcontrol$cores <- cores
      optimcontrol$verbose <- verbose
      optimcontrol$priors <- as.logical(priors)
      optimcontrol$standata=standata
      optimcontrol$sm=sm
      optimcontrol$init=inits
      optimcontrol$plot=plot
      optimcontrol$matsetup <- data.frame(ctm$modelmats$matsetup)
    }

    if(!optimize && !is.null(inits)){
      if('list' %in% class(inits)){
        staninits=inits
      } else {
        if(initOptim){ #then first optimize to get inits
          if(!intoverpop && length(unique(datalong[[ctm$subjectIDname]]) > 1) && any(ctm$pars$indvarying)) stop('Cannot optimize to get inits unless intoverpop=TRUE')
          optimcontrol$init <- NULL
          optimcontrol$tol=1e-7
          if(!intoverpop & ! intoverstates) stop('Cannot initialize with optimization unless intoverpop and intoverstates are set to TRUE')
          inits <- do.call(stanoptimis,optimcontrol)$rawest
        }
        sf <- stan_reinitsf(sm,standata)
        staninits <- list()
        if(chains > 1){ #set all chains to same inits
          for(i in 1:chains){
            staninits[[i]]<-constrain_pars(sf,inits+rnorm(length(inits),0,.01))
          }
        }
      }
    }


    if(!optimize){

      #control arguments for rstan
      # if(is.null(control$adapt_term_buffer)) control$adapt_term_buffer <- min(c(iter/10,max(iter-20,75)))
      if(is.null(control$adapt_delta)) control$adapt_delta <- .8
      if(is.null(control$adapt_window)) control$adapt_window <- 5
      if(is.null(control$max_treedepth)) control$max_treedepth <- 10
      if(is.null(control$adapt_init_buffer)) control$adapt_init_buffer=2
      if(is.null(control$stepsize)) control$stepsize=.001
      if(is.null(control$metric)) control$metric='diag_e'


      message('Sampling...')
      #
      stanargs <- list(object = sm,
        init_r=.03,
        save_warmup=as.logical(plot),
        refresh=20,
        iter=iter,
        data = standata, chains = chains, control=control,
        cores=cores,
        ...)
      if(!is.null(inits)) stanargs$init=staninits

      if(plot==TRUE) stanfit <- suppressWarnings(list(stanfit=do.call(stanWplot,stanargs))) else stanfit <- suppressWarnings(list(stanfit=do.call(sampling,stanargs)))

      #find the median sample and compute kalman scores etc for this
      # browser()
      e=rstan::extract(stanfit$stanfit)
      # middle <- which(abs(e$ll-quantile(e$ll,.5)) == min(abs(e$ll-quantile(e$ll,.5) )))
      # middle <- which(e$ll==max(e$ll))
      stanfit$rawposterior <- t(stan_unconstrainsamples(fit = stanfit$stanfit,standata = standata))
      stanfit$rawest <- apply(stanfit$rawposterior,2,median)
    }

    if(optimize==TRUE) {
      stanfit <- do.call(stanoptimis,optimcontrol) #eval(parse(text=opcall))

      #update data that may have changed during optimization
      for(ni in names(stanfit$standata)){
        standata[[ni]] <- stanfit$standata[[ni]]
      }
      if(ctm$n.TIpred>0){
        ctm$modelmats$TIPREDEFFECTsetup <- stanfit$standata$TIPREDEFFECTsetup
        ms <- ctm$modelmats$matsetup
        ms$tipred <- 0L
        defining <- .ctMatsetupFreeRows(ms, defining = TRUE)
        parswithtipreds <- sort(unique(ms$param[defining]))
        parswithtipreds<-parswithtipreds[apply(stanfit$standata$TIPREDEFFECTsetup,1,sum)>0]
        ms$tipred[defining & ms$param %in% parswithtipreds] <- 1L
        ctm$modelmats$matsetup <- ms
      }
    }

    # if(is.na(STAN_NUM_THREADS)) Sys.unsetenv('STAN_NUM_THREADS') else Sys.setenv(STAN_NUM_THREADS = STAN_NUM_THREADS) #reset sys env
  } # end if fit==TRUE
  #convert missings back to NA's for data output
  standataout <- .ctStandataNA(standata)

  setup=list(recompile=recompile,idmap=standata$idmap,matsetup=ctm$modelmats$matsetup,matvalues=ctm$modelmats$matvalues,
    popsetup=ctm$modelmats$matsetup[.ctMatsetupFreeRows(ctm$modelmats$matsetup),],
    popvalues=ctm$modelmats$matvalues[.ctMatsetupFreeRows(ctm$modelmats$matsetup),],
    extratforms=ctm$modelmats$extratforms)
  if(fit) {
    stanfit$transformedparsfull <- suppressMessages(stan_constrainsamples(sm = sm,standata = standata,
      savesubjectmatrices = TRUE, samples = matrix(stanfit$rawest,1),cores=1,savescores=TRUE,pcovn=5000))

    out <- list(args=list(input=args, resolved=argsresolved),
      setup=setup,
      stanmodeltext=stanmodeltext, data=standataout, ctdatastruct=datalong[c(1,nrow(datalong)),],standata=standata,
      ctstanmodelbase=ctstanmodel, ctstanmodel=ctm,stanmodel=sm, stanfit=stanfit)
    out$backend <- 'stan'
    class(out) <- c('ctStanFit','ctFit')
    # Not before this point: naming the raw draws and covariance needs
    # `ctstanmodelbase`, which the fit only has once assembled here.
    out <- .ctFitNameRawUncertainty(out)
    out$stanfit$kalman<-suppressMessages(ctKalmanArray(out,pointest = TRUE))
  }

  if(!fit) out=list(args=list(input=args, resolved=argsresolved),setup=setup,
    stanmodeltext=stanmodeltext,data=standataout,  ctdatastruct=datalong[c(1,nrow(datalong)),],standata=standata,
    ctstanmodelbase=ctstanmodel,  ctstanmodel=ctm)


  return(out)
}

#' @export
ctStanFit <- ctFit

