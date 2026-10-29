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
# backend's own retry: stan restarts from fresh random values, julia resumes
# from the better point its probe measured or, when the certification finds
# the fit short, from where it stopped (`.ctBackendCorrectResult()`). It used
# to be listed below twice, stan-only and julia-only, and a list lookup returns
# the first -- so julia refused it with stan's message and stan dropped it with
# the julia-only names.
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
  coreset = list(only = 'julia',
    inert = function(v) isFALSE(v) || identical(v, 'auto'),
    msg = paste0("runs the julia optimiser's coreset stages, stretches of ",
      "L-BFGS on a weighted few of the units on the Laplace route. The stan ",
      "path has no equivalent")),
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
  lbfgs_scores = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("runs the julia optimiser with a second model of the ",
      "curvature from the per-unit score outer products. The stan path has no ",
      "equivalent")),
  lbfgs_local = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("starts the julia optimiser's curvature from exact Hessian ",
      "columns for parameters few units read. The stan path has no ",
      "equivalent")),
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
  hessianreuse = list(only = 'julia',
    inert = function(v) FALSE,
    msg = paste0("sets how far the julia estimate may move from where its ",
      "Hessian was formed for that Hessian still to give the standard errors. ",
      "The stan path forms its Hessian at the estimate")),
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
  accuracy = list(only = 'julia',
    inert = function(v) is.null(v),
    msg = paste0("asks the julia uncertainty stage to sample to an accuracy ",
      "where its method falls short, and the stan path has no sampler to ",
      "escalate to. Drop it")),
  sampleControl = list(only = 'julia',
    inert = function(v) is.null(v) || !length(v),
    msg = paste0("holds the sampler settings for the julia uncertainty ",
      "stage's escalation under accuracy, which the stan path does not do. ",
      "Drop it")),
  callback = list(only = 'julia',
    inert = function(v) is.null(v),
    msg = paste0("is called by the julia engine while the fit runs, and the ",
      "stan optimizer has no such hook. Use verbose=1 for stan's own ",
      "iteration output")),
  checkpoint = list(only = 'julia',
    inert = function(v) is.null(v),
    msg = paste0("writes the julia optimiser's point and modes to a file while ",
      "the fit runs, and resumes from it. The stan optimizer has no such hook")),
  checkpoint_every = list(only = 'julia',
    inert = function(v) TRUE,
    msg = paste0("sets how often the julia optimiser writes its checkpoint. The ",
      "stan optimizer has no checkpoint")),
  saveEffects = list(only = 'julia',
    inert = function(v) !isTRUE(v),
    msg = paste0("keeps every draw of every random effect from the julia ",
      "sampler. Use ctSubjectPars() on a stan fit for per-subject draws")),
  progress = list(only = 'julia',
    inert = function(v) !isTRUE(v),
    msg = paste0("forces the julia engine's progress line on or off, and the ",
      "stan optimizer does not draw one. Use verbose=1")),
  tipredMissing = list(only = 'julia',
    inert = function(v) identical(v, 'auto'),
    msg = paste0("chooses whether julia imputes or samples missing TI ",
      "predictors; stan imputes them when optimizing and samples them when ",
      "sampling. Drop it")),
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
  # The Laplace check every julia Laplace fit gets, and the integrated refit
  # asked for. FALSE describes what stan does, so it is accepted there, as
  # `stochastic=FALSE` is on julia.
  laplace_check = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("checks a julia Laplace fit's error by integrating over its ",
      "random effects, and stan has no Laplace term to check. Drop it")),
  laplace_refit = list(only = 'julia',
    inert = function(v) isFALSE(v),
    msg = paste0("refits a julia Laplace fit with its random effects integrated, ",
      "and stan has no Laplace term to refit. Drop it"))
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

  args <- .ctFitReplayArgs(oldfit)
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
    # Each group's posterior was solved from the old data; with new data it is
    # another set of groups' posterior, at the same estimate.
    if(newdata && !is.null(oldfit$randomeffects)) {
      at <- tryCatch(.ctEffectEvaluate(spec, .ctFitRawEstimate(oldfit)),
        error = function(e) NULL)
      oldfit$randomeffects <- if(is.null(at)) NULL else
        .ctRandomEffectsStore(oldfit, at$posterior)
    }
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

# The arguments `oldfit` was made with, as `ctFit()` takes them: what
# `ctFitUpdate()` replays, and what `ctScoreTest()` prepares an extended model
# with, so the extension is the fitted model plus its new coordinates and
# nothing else.
.ctFitReplayArgs <- function(oldfit) {
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
  # `priors` is captured after ctFit() has reduced it to a logical, so
  # 'randomCorr' reads TRUE there -- and TRUE replayed is a prior on every
  # coordinate, a different estimator. `priorscope` beside it still says
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
  args
}



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
# levels above it: T0VAR plus the subject level population covariance. The
# study level reaches an individual by moving that individual's parameter
# values, so adding a study covariance into T0cov as well would count it twice.
#
# They are still *fixed* rather than left free, because taking them out of the
# parameter vector is the whole job of this function, and every route reads
# them as the same thing: given a subject's effects its initial state is known,
# and this diagonal is what is left of its variance. Julia's Laplace route and
# a sampled fit read it so, and Stan's and julia's augmented routes add the
# population covariance to it, so all of them fit one model. With MANIFESTVAR
# fixed at zero it is also all of the first observation's innovation variance,
# so it cannot be zero: a floor, at sd 1e-3 (a variance of 1e-6), negligible
# beside any population variance worth estimating. A floor that small puts a
# curvature of 1e6 or more on the Laplace inner solve, which stops on its step
# for that reason (`_laplace_newton_unit_mode`). Positive for the R-side paths
# that build a covariance straight from T0VAR too -- ctGenerate(backend='r'),
# ctGraph, ctModelLatex -- which would meet a singular matrix.
T0VARredundancies <- function(ctm) {
  # Per-level columns, so `indvarying` below is the subject level's logical.
  ctm <- .ctResolveLevels(ctm)
  whichT0VAR_T0MEANSindvarying <- ctm$pars$matrix %in% 'T0VAR'  &
    is.na(ctm$pars$value) &
    (ctm$pars$row %in% ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying] |
        ctm$pars$col %in% ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying])
  if(any(whichT0VAR_T0MEANSindvarying)){
    message('T0VAR rows/columns for latents with individually varying T0MEANS disabled: RAWPOPVAR gives their covariance.')
    ctm$pars$value[whichT0VAR_T0MEANSindvarying & ctm$pars$col == ctm$pars$row ] <- 1e-3
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
#' the model, as \code{model$poprank <- 2}. A rank of 0 is allowed at a level
#' whose every varying parameter has an own variance (\code{popunique}).
#' @param popunique Own variances on a reduced-rank level: the parameters
#' selected get a variance of their own beside the \code{poprank} shared
#' dimensions, so that level's covariance is \code{L \%*\% t(L)} plus a
#' diagonal -- the factor-analytic form. \code{TRUE} selects every parameter
#' varying at a level, \code{FALSE} none, and names select matrices or
#' parameters, e.g. \code{popunique = c(subject = 'MANIFESTMEANS')}. One
#' unnamed value applies to every level; named entries, or a named list, apply
#' per level. Each selected level needs a \code{poprank} below its number of
#' varying parameters; \code{poprank = 0} with every parameter selected gives
#' independent random effects. Reported as
#' \code{popunique_<parameter>.<level>}, on the same raw scale as a full-rank
#' sd. Requires \code{backend = 'julia'} and \code{intoverpop = 'laplace'} or
#' sampling (\code{'none'}); must be passed by name.
#' @param popblocks Shared dimensions in blocks, on a reduced-rank level: a
#' list of blocks, each a list with \code{rank}, how many shared dimensions,
#' and \code{effects}, which parameters load on them -- \code{TRUE} for every
#' parameter varying at the level, or names of matrices or parameters, as for
#' \code{popunique}. A parameter loads only on the dimensions of the blocks
#' that name it, so parameters with no block in common are uncorrelated, and
#' one named in no block needs an own variance (\code{popunique}). A level's
#' rank is the sum of its blocks' ranks, so \code{poprank} need not be given
#' for it; \code{'auto'} there defers to the blocks, and a number must equal
#' the sum. One unnamed list of blocks applies to every level; a list named by
#' level gives each level its own. For example, person-level item intercepts
#' that do not covary with a latent level but may covary with each other and
#' with the dynamics: \code{popblocks = list(level = list(rank = 1, effects =
#' c('lev', 'DRIFT')), items = list(rank = 2, effects = c('MANIFESTMEANS',
#' 'DRIFT')))}. Within a block, its first \code{rank} parameters, in the
#' order the model lists them, form its basis. Loadings are reported as
#' \code{poploading_<parameter>_dim<j>}, the dimensions numbered across the
#' blocks in order. Requires \code{backend = 'julia'} and \code{intoverpop =
#' 'laplace'} or sampling (\code{'none'}); must be passed by name.
#' @param weights Sampling weights, one per subject: the name of a column of
#' \code{datalong} that is constant within subject, or a numeric vector named
#' by subject id (or unnamed, in the order subjects first appear in the data).
#' Each subject's contribution to the likelihood -- its marginal likelihood,
#' random effects integrated out -- is multiplied by its weight (the weighted
#' pseudo-likelihood), and the prior is not; a weight of zero leaves the
#' subject out. Weights are normalised to a mean of one over the weighted
#' subjects, and standard errors come from the sandwich estimator unless
#' \code{optimcontrol$uncertainty} says otherwise. \code{fit$estimate$loglik}
#' stays the data's log likelihood at the estimate;
#' \code{fit$estimate$loglik_weighted} is the weighted one. Requires
#' \code{backend = 'julia'} and optimisation; with nested ids, weights must be
#' constant within each group of the outermost level. Must be passed by name.
#' @param intoverpop How declared individual differences (random effects) are
#' handled: integrated out of the likelihood, or sampled with everything else.
#'
#' \code{'auto'}, the default, integrates them out when optimizing, as
#' \code{TRUE} does, and samples them (\code{'none'}) when
#' \code{optimize=FALSE} -- except on julia where every subject's likelihood is
#' Gaussian in its effects (Gaussian indicators, effects only on means through
#' the affine transform), where they are integrated exactly and the population
#' parameters sampled (\code{'laplace'}).
#'
#' \code{TRUE} integrates them out: by \code{'laplace'} on
#' \code{backend='julia'}, and by \code{'augmented'} on stan, which has only
#' that route, or with \code{intoverstates=FALSE}, which composes with it
#' only. \code{fit$args$resolved$intoverpopreason} records why. Before 4.0.0
#' \code{TRUE} always meant \code{'augmented'}.
#'
#' \code{'laplace'} (\code{backend='julia'}) integrates an effect that shifts
#' a mean in a linear Gaussian model -- a manifest mean, an initial mean, a
#' continuous intercept -- exactly, inside each subject's filter, and every
#' other effect by a Laplace approximation at its mode, which \code{summary()}
#' reports. \code{'augmented'} makes each random effect a state of the filter:
#' exact for effects that shift a mean with Gaussian indicators, a
#' linearisation elsewhere that shrinks random effects on the dynamics, and
#' variation on a DIFFUSION or MANIFESTVAR parameter is then only partly
#' identified, which ctFit warns of. \code{FALSE}, or \code{'none'}, samples them, so it needs
#' \code{optimize=FALSE}: maximising over every subject's effects would drive
#' the population variance to zero. With \code{optimize=FALSE} on julia,
#' \code{'none'} (and so \code{'auto'}, but for the exception above) samples the
#' joint posterior of parameters and random effects; \code{'laplace'} and \code{'augmented'}
#' sample the population parameters with the effects integrated out, by the
#' Laplace approximation or the filter, and \code{TRUE} by whichever of the two
#' it resolves to.
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
#' NUTS through a learnt map on the joint posterior of parameters and random
#' effects, or NUTS or tempered SMC when there are none or \code{intoverpop}
#' integrates them out -- slower,
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
#' \code{intoverpop='laplace'}, \code{saem}, \code{laplace_check} and
#' \code{laplace_refit}.
#' \code{\link{ctOptimControl}} lists every name with its backend and default.
#' A misspelled name, or one the chosen backend cannot honour, is an error.
#' @param nopriors deprecated, use priors argument. logical. If TRUE, any priors are disabled -- sometimes desirable for optimization.
#' @param iter \strong{Deprecated} -- use \code{sampleControl$iter}. Still
#' honoured, with a warning.
#' @param inits \code{NULL} (the default) for ctsem's own starting values, or
#' a vector of starting values on the raw (unconstrained) scale, one per free
#' parameter, as a fit holds its estimate: \code{fit$estimate$raw} on julia,
#' \code{fit$stanfit$rawest} on stan. On julia, a fit itself: its estimates are
#' matched to this model's parameters by name, so a model extended by a
#' predictor or a smooth starts where the fit stopped, with the new
#' coefficients at zero. With \code{backend='stan'} and
#' \code{optimize=FALSE}, the string \code{'optimize'} optimises first and
#' starts the chains there.
#' @param priors \code{FALSE} (the default), \code{TRUE} or \code{'randomCorr'}.
#' \code{FALSE} uses no prior: maximum likelihood. \code{TRUE} adds ctsem's
#' \code{normal(0,1)} raw-scale prior to every parameter, making an optimised
#' fit maximum a posteriori. \code{'randomCorr'} applies it to the
#' random-effect correlations only, at every level: a bounded correlation on an
#' unbounded coordinate can walk along a direction the data do not determine,
#' and the prior holds it. Sampling (\code{optimize=FALSE}) without a prior on
#' every parameter warns, because a parameter the data leave flat has an
#' improper posterior; use \code{priors=TRUE} there. \code{'randomCorr'} needs
#' \code{backend='julia'}.
#' @param cores number of cpu cores to use: a positive whole number, or
#' \code{'maxneeded'} for all but one (capped at the number of chains on
#' stan). Defaults to \code{getOption("mc.cores", 2)}. It is the most CPU
#' cores the call uses at once, counting every process it starts: on julia
#' the threads a fit uses, BLAS included, and a sampling run's worker
#' processes, which share it between them. A julia fit at \code{cores > 1} is
#' not reproducible to the last decimal, because its sums are divided among
#' workers in an order that varies from run to run; use \code{cores = 1} for a
#' before-and-after comparison. The Julia
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
#' draw is kept; and \code{processes} (FALSE), whether the chains run in worker
#' processes -- at most \code{cores} at once, so four chains on two cores run
#' two at a time -- rather than together in this session.
#'
#' For \code{backend='julia'} also \code{sampler} (\code{'map'}, NUTS
#' through a learnt map, on the joint posterior; \code{'saem'}, SAEM's kernel,
#' or \code{'nuts'} by choice), \code{placement} (the map starts from an
#' optimised fit; for \code{'saem'} and \code{'nuts'} on the joint posterior
#' \code{'saem'}: SAEM's state places the chains, so no Laplace optimum is
#' needed; \code{'fit'} places them around an optimised fit), the stopping
#' target \code{minESS} (200, of the
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
#'
#' To start only some initial conditions there, name those cells
#' \code{'stationary'} in the model (\code{backend='julia'}); every other
#' initial condition stays as written. A T0MEANS cell takes its entry of the
#' stationary mean, asymCINT. A T0VAR diagonal cell gives its latent the
#' stationary variance, asymDIFFUSION's diagonal, keeping that latent's
#' T0VAR correlations with the others. A T0VAR off-diagonal cell gives the
#' pair the stationary covariance, and needs both their variances, and every
#' covariance among the latents it joins, named too. The stationary moments
#' are those of the marked latents and every latent that drives them, so a
#' latent driven by others works, and DRIFT, CINT and DIFFUSION need to be
#' fixed within each subject only for those.
#' @param forcerecompile logical. For \code{backend='stan'}, recompiles the
#' stan model regardless of apparent need (for development). For
#' \code{backend='julia'}, compiles the model's state-dependent expressions --
#' matrix cells that read a latent state, and their derivatives -- and any
#' parameter transforms of your own into code for this model, rather than
#' evaluating them as data. The model's first fit in a session then compiles
#' for a minute or more, and its gradients run faster: by up to a third for a
#' heavily nonlinear model (many state-dependent cells, or long formulas in
#' them), and by little otherwise. Worth it only when such a model's fit takes
#' many minutes.
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
#'     dtmat=rbind(c(rep(.5,8),3,rep(.5,Tpoints-9))))
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
  nlcontrol = list(), nopriors=NA, priors=FALSE, chains=2,
  cores=getOption("mc.cores", 2L),
  inits=NULL,
  compileArgs=list(),
  forcerecompile=FALSE,saveCompile=TRUE,savescores=FALSE,
  savesubjectmatrices=FALSE, saveComplexPars=FALSE,
  gendata=FALSE,
  sampleControl=list(),
  control=list(),verbose=0,..., popunique=NULL, popblocks=NULL, weights=NULL,
  ctstanmodel){

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
  # Starting values matched by name from another fit are the julia path's.
  if(inherits(inits, 'ctFit') && !(backend %in% 'julia' && .ctFitIsJulia(inits)))
    stop("inits as a fit takes a backend='julia' fit, for a backend='julia' ",
      "fit; otherwise pass its raw estimate.", call.=FALSE)
  # Before anything prepares, warns or fits.
  # Sampling weights (R/ctWeights.R): a weighted pseudo-likelihood, which the
  # julia engine maximises. Sampling it would not give a posterior, and the
  # state-explicit route optimises a joint density weights do not apply to.
  if(!is.null(weights)){
    if(!backend %in% 'julia') stop("weights requires backend='julia'.", call.=FALSE)
    if(!isTRUE(optimize)) stop("weights are for optimised fits: a weighted ",
      "likelihood is not a posterior to sample. Use optimize=TRUE; standard ",
      "errors then come from the sandwich estimator.", call.=FALSE)
    if(!isTRUE(intoverstates)) stop("weights require intoverstates=TRUE.",
      call.=FALSE)
  }
  # Whether `poprank` was asked for or merely defaulted. Taken here because
  # `missing()` has to be evaluated before the argument is touched, and it
  # decides whether an inapplicable rank is an error or a no-op. A rank stated
  # on the model sets it too, once the model is read (`model$poprank` below).
  # NULL is the argument left out: the default, and not a rank asked for.
  poprankexplicit <- !missing(poprank) && !is.null(poprank)
  if(is.null(poprank)) poprank <- 'auto'
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
  # 'randomCorr': a prior on the random-effect correlations and nothing else.
  # A correlation is bounded and its coordinate is not, so along a direction the
  # data do not determine the coordinate can walk until the transform
  # saturates -- measured on a 4x4 population covariance with 50 subjects:
  # coordinates at 118, and every correlation in the summary `NA`. It was the
  # julia default from 2026-09-22 (4b2dbf05) until 2026-10-19, when the
  # default went back to FALSE, the released one (Charles: the prior was
  # covering another issue). See `.ctBackendRandomCorrPriorSpec`.
  priorscope <- if(is.character(priors)) match.arg(priors, 'randomCorr') else
    if(isTRUE(priors)) 'all' else 'none'
  if(!is.character(priors) && !is.logical(priors)) stop(
    "priors must be TRUE, FALSE, or 'randomCorr'", call.=FALSE)
  # Sampling needs a prior wherever the likelihood can go flat -- with none, or
  # under 'randomCorr', a population sd has none, and on the ?ctFit example the
  # chains wandered off along popsd_mm (R-hat 2.1). The argument is not changed
  # for it: the default is one value whatever the route, and sampling without
  # a prior on every parameter warns below, before the run.
  # The generated Stan model builds its priors in, so a subset of coordinates
  # is not expressible there.
  if(priorscope %in% 'randomCorr' && !backend %in% 'julia') stop(
    "priors='randomCorr' applies a prior to a subset of coordinates, which the ",
    "generated Stan model cannot express. Use priors=TRUE or priors=FALSE with ",
    "backend='stan', or backend='julia'.", call.=FALSE)
  priors <- !priorscope %in% 'none'

  if(any(!is.na(derrind))) warning('derrind argment is deprecated, computed automatically now')

  datalong <- data.frame(datalong)

  # Smooth TI predictor terms: their basis columns, built from these data where
  # the model has none yet, and kept on the model the fit stores
  # (R/ctTIpredSmooth.R).
  if(length(ctstanmodel$TIsmooths)){
    if(!backend %in% 'julia') stop("Smooth TI predictor terms, s(), need ",
      "backend='julia'.", call.=FALSE)
    if(!isTRUE(optimize)) stop("Smooth TI predictor terms are estimated by ",
      "optimisation, which sets their smoothness: use optimize=TRUE, then ",
      "ctFitUncertainty(fit, 'sample') to sample with the smoothness held at ",
      "its estimate.", call.=FALSE)
    smoothed <- .ctTIsmoothData(ctstanmodel, datalong,
      includeOutcome = !isFALSE(optimcontrol$tipredMissingIncludeOutcome))
    ctstanmodel <- smoothed$model
    datalong <- smoothed$data
  }

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

  nlcontrol <- .ctNlcontrolDefaults(nlcontrol, backend)

  args=c(as.list(environment()), list(...)) #as.list((match.call(expand.dots=FALSE)))
  # A fit given as starting values is not kept in the new fit's call: its
  # estimates matched this model's parameters once, and a replay starts afresh.
  if(inherits(args$inits, 'ctFit')) args['inits'] <- list(NULL)
  args$datalong <- NULL
  args$model <- NULL
  args$ctstanmodel <- NULL

  setup <- .ctFitModelSetup(ctstanmodel, backend = backend,
    intoverpop = intoverpop, intoverstates = intoverstates, optimize = optimize,
    priors = priors, priorscope = priorscope, stationary = stationary,
    binomial = binomial, poprank = poprank, poprankexplicit = poprankexplicit,
    popblocks = popblocks, popunique = popunique, nlcontrol = nlcontrol,
    saveComplexPars = saveComplexPars,
    nsubjects = length(unique(datalong[, ctstanmodel$subjectIDname])),
    checks = TRUE, fit = fit, estonly = isTRUE(optimcontrol$estonly))
  ctm <- setup$ctm
  intoverpop <- setup$intoverpop
  intoverpopmethod <- setup$intoverpopmethod
  intoverpopreason <- setup$intoverpopreason
  laplacerank <- setup$laplacerank
  popregression <- setup$popregression
  recompile <- FALSE
  # The data against the indicator types: read here, as they read the data.
  if(any(ctm$manifesttype %in% 1)) .ctDataBinary(datalong, ctm)
  if(any(ctm$manifesttype %in% 2)) .ctDataCategories(datalong, ctm)
  if(any(ctm$manifesttype %in% 3)) .ctDataCounts(datalong, ctm)
  if(any(ctm$manifesttype %in% 4)) .ctDataCensored(datalong, ctm)



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


  standata <- .ctPrepareData(ctm,datalong,optimize=optimize, sameInitialTimes=sameInitialTimes, #bit slow
    tipredMissing = if (identical(backend, 'julia')) 'leave' else if (optimize) 'impute' else 'sample')
  standata$verbose=as.integer(verbose)
  standata$savesubjectmatrices=as.integer(savesubjectmatrices)
  standata$gendata=as.integer(gendata)

  if(standata$savesubjectmatrices==1L) savescores = TRUE
  standata$savescores=as.integer(savescores)

  # Reached only when the caller asks for 'maxneeded' by name. It is not a
  # default and must not become one: it resolves to the whole machine, and CRAN
  # allows at most two cores unless the user has asked for more. It was the
  # `optimize=FALSE` default until 4.0.0.
  #
  # Resolved here rather than after the julia branch below, which returns before
  # ever reaching the old resolution site: the string arrived at
  # .ctFitJuliaBackend(), `as.integer()` made it NA, and the NA guard there
  # turned it into 1 -- so asking for every core silently got one.
  #
  # `chains` does not cap the julia backend: its chains share `cores`
  # (`.ctBackendWorkerThreads()`), and the rest of its parallelism is over
  # subject chunks, run at whatever ceiling it is given, so the ceiling is
  # simply the machine.
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
    #   FALSE     'laplace'      NUTS/SMC, Laplace marginal npar
    #   FALSE     'augmented'    NUTS/SMC, filter marginal  npar
    #   FALSE     'none'         sampled, joint            npar + effects
    #
    # The route is what `intoverpop` resolved to: TRUE takes 'laplace' or
    # 'augmented' as the model suits, and 'auto' under optimize=FALSE takes
    # 'none'. The joint row is drawn through a learnt map from an optimised
    # placement unless `sampleControl` says otherwise (`sampler`, `placement`).
    # Under optimize=FALSE the route is what is sampled, and nothing else
    # says it (`.ctBackendSampleRoute()`): 'auto' resolves to 'laplace' there
    # only where the Laplace term is exact, which the prepared specification
    # says (`exactmarginal`, `.ctFitJuliaBackend()`), and the fit records.
    #
    # Only 'laplace' is julia-only; the guard for that is above. `'none'` is
    # the third route the engine needs, and it prepares the Laplace structure
    # -- which is what says *which* parameters vary -- without integrating.
    juliaintoverpop <- if(identical(intoverpopmethod,'laplace')) 'laplace' else
      if(intoverpop) 'augmented' else
        if(!optimize && .ctAnyVarying(ctm)) 'none' else
          'augmented'
    juliafit <- .ctFitJuliaBackend(datalong=datalong, model=ctm, prepared_data=standata, inits=inits,
      cores=cores, optimcontrol=optimcontrol, forcerecompile=forcerecompile,
      verbose=verbose, fit=fit, priors=priors, priorscope=priorscope,
      optimize=optimize,
      chains=if(chainsgiven) chains else NULL, iter=iter, control=control,
      intoverpop=juliaintoverpop, intoverstates=intoverstates, weights=weights,
      exactmarginal=!isTRUE(optimize) && identical(juliaintoverpop, 'none') &&
        !is.na(intoverpopreason) && isTRUE(intoverstates))
    if(!identical(juliaintoverpop, .ctBackendIntOverPop(.ctBackendSpec(juliafit)))){
      argsresolved$intoverpop <- .ctBackendIntOverPop(.ctBackendSpec(juliafit))
      argsresolved$intoverpopreason <- paste0("optimize=FALSE samples the ",
        "exact marginal: every subject's likelihood is Gaussian in its effects")
    }
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
    # With the smoothing parameters the fit estimated, so a refit or an update
    # starts from them (R/ctTIpredSmooth.R).
    if(length(ctstanmodel$TIsmooths) && is.list(juliafit$model_spec) &&
        !is.null(juliafit$model_spec$penalty))
      juliafit$modelbase$TIsmooths <- .ctTIsmoothKeepLambda(ctstanmodel$TIsmooths,
        juliafit$model_spec$penalty)
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
      try(print(ctTracePlot(juliafit)), silent=TRUE)
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

