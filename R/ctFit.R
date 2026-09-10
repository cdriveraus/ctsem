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

# Honoured by both backends, with the same meaning.
.ctOptimcontrolShared <- c('tol','g_tol','x_tol','maxiter','lbfgs_memory',
  'initsd','carefulfit','estonly','finishsamples','uncertainty',
  'uncertaintyDraws','uncertaintyControl')

# What is left after the vocabulary above is a genuine capability difference: a
# phase or a hook one backend has and the other does not. Each such name is
# listed here with the backend that has it and, where there is one, the value
# that asks for nothing.
#
# The rule is on the value, not the name. A value that asks for a capability the
# chosen backend lacks is refused; a value that merely *describes* what that
# backend already does is accepted, because it is true -- `stochastic = FALSE`
# on julia names the deterministic optimiser julia already runs, and
# `gradient = 'adjoint'` on stan names the reverse-mode gradient Stan's autodiff
# already takes. That is why several honest calls carrying `stochastic = FALSE`
# keep working.
.ctOptimcontrolSplit <- function() list(

  # -- stan only: the stochastic-gradient family. The julia optimiser is L-BFGS
  # over the full data and has no SGD phase for these to tune.
  stochastic = list(only = 'stan',
    inert = function(v) !isTRUE(v),
    msg = paste0("asks for stochastic gradient descent; the julia backend ",
      "optimises with L-BFGS over the full data. Drop it -- stochastic=FALSE ",
      "is what julia does")),
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
  stallretries = list(only = 'stan',
    inert = function(v) isTRUE(all(v == 0)),
    msg = paste0("restarts the stan optimizer from fresh values when it stops ",
      "short of a maximum; the julia backend warms every fit from the priors ",
      "instead. Use optimcontrol$carefulfit")),
  stalltol = list(only = 'stan',
    inert = function(v) FALSE,
    msg = paste0("is the gradient at which stan calls a fit stalled and ",
      "restarts it, and the julia backend has no such retry. Use ",
      "optimcontrol$carefulfit")),

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
      "posterior, so on stan it always does. Drop it"))
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

#' Update a ctStanFit object
#'
#' Either to include different data, or because you have upgraded ctsem and the internal data structure has changed.
#'
#' @param oldfit fit object to be upgraded
#' @param data replacement long format data object
#' @param recompile whether to force a recompile -- safer but slower and usually unnecessary.
#' @param refit if TRUE, refits the model using the old estimates as a starting point. Only applicable for
#' optimized fits, not sampling.
#' @param ... extra arguments to pass to ctFit
#'
#' @return updated ctStanFit object.
#' @aliases ctStanFitUpdate
#' @export
#'
#' @examples
#' newfit <- ctFitUpdate(ctstantestfit,refit=FALSE)

ctFitUpdate <- function(oldfit, data=NA, recompile=FALSE,refit=FALSE,...){

  if(!refit) message('Trying to do a quick update -- if there are problems, try with refit=TRUE for more robustness')

  dots <- list(...)
  # `$args$input` -- the literal call, still carrying 'auto'/'maxneeded' and
  # whatever else was unresolved -- not `$args$resolved`, which would freeze
  # this refit at whatever a previous 'auto' happened to route to instead of
  # letting it re-route against the new data or overrides in `...`.
  args <- as.list(oldfit$args$input)
  for(n in names(dots)){
    args[[n]] <- dots[[n]]
  }
  if(length(oldfit$stanfit$stanfit@sim) > 0) refit=FALSE
  args$fit <- refit
  args$inits <- oldfit$stanfit$rawest
  args$model <- oldfit$ctstanmodelbase
  args$ctstanmodel <- NULL

  newargs <- as.list(args(ctFit))
  for(argi in names(args)){
    if(argi %in% names(args)) newargs[[argi]] <- args[[argi]] else message(argi, ' is no longer a valid argument, dropping...')
  }


  if(length(data==1)) args$datalong <- standatatolong(oldfit$standata,origstructure = TRUE,ctm=oldfit$ctstanmodelbase)
  if(length(data) > 1) args$datalong <- data
  newfit <- do.call(ctFit,args)

  if(!refit){
    oldfit$standata <- newfit$standata
    if(oldfit$ctstanmodel$recompile || recompile) oldfit$stanmodel <- rstan::stan_model(model_code = newfit$stanmodeltext) else
      oldfit$stanmodel <- stanmodels$ctsm
  }
  if(refit) oldfit <- newfit
  return(oldfit)
}

#' @export
ctStanFitUpdate <- ctFitUpdate


T0VARredundancies <- function(ctm) { #check for redundant T0VAR parameters (because indvarying t0means) and disable
  whichT0VAR_T0MEANSindvarying <- ctm$pars$matrix %in% 'T0VAR'  &
    is.na(ctm$pars$value) &
    (ctm$pars$row %in% ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying] |
        ctm$pars$col %in% ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying])
  if(any(whichT0VAR_T0MEANSindvarying)){
    message('Free T0VAR parameters as well as indvarying T0MEANS -- fixing T0VAR pars to diag matrix of 1e-6')
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



#' Fit a ctsem model
#'
#' Fits a ctsem model specified via \code{\link{ctModel}} with type either 'ct' or 'dt'.
#' \code{ctStanFit} is maintained as a backward-compatible alias.
#'
#' @aliases ctStanFit
#' @seealso \code{\link{ctIdentify}} reports which parameters the data can
#' inform, before a fit is spent finding out; \code{\link{ctTracePlot}} draws
#' the optimisation trace a julia fit records.
#' @param datalong long format data containing columns for subject id (numeric values, 1 to max subjects), manifest variables,
#' any time dependent (i.e. varying within subject) predictors,
#' and any time independent (not varying within subject) predictors.
#' @param model model object as generated by \code{\link{ctModel}} with type='ct' or 'dt', for continuous or discrete time
#' models respectively.
#' @param ctstanmodel Deprecated. Use \code{model}.
#' @param stanmodeltext already specified Stan model character string, generally leave NA unless modifying Stan model directly.
#' (Possible after modification of output from fitting with argument fit=FALSE)
#' @param intoverstates logical indicating whether or not to integrate over
#' latent states using a Kalman filter. \code{FALSE} instead makes the
#' latent states part of what is estimated: the target becomes the joint
#' density of the data and the states, and the fit reports the trajectory
#' alongside the parameters in \code{fit$estimate$states}.
#'
#' With \code{backend='julia'} that route is exact -- no Gaussian
#' assumption is made about the state anywhere, where the filter's update
#' for a binary, ordinal or count indicator is an assumed-density
#' projection. Pair it with \code{optimize=FALSE}: sampling the joint
#' density gives the posterior of parameters and states together, while
#' \code{optimize=TRUE} gives its joint mode, whose variance parameters
#' are biased downward. Standard errors for an optimised fit come from the
#' Hessian with the states profiled out, and \code{uncertainty} is
#' restricted to \code{'hessian'} for that reason.
#' Generally recommended to set TRUE unless using non-gaussian measurement model.
#' @param binomial Deprecated. Logical indicating the use of binary rather than Gaussian data, as with IRT analyses.
#' This now sets \code{intoverstates = FALSE} and the \code{manifesttype} of every indicator to 1, for binary.
#' @param fit If TRUE, fit specified model using Stan, if FALSE, return stan model object without fitting.
#' @param poprank Rank of the population covariance of the individually
#' varying parameters.
#'
#' \code{'auto'}, the default, uses the number of varying parameters that reach
#' the observation mean. That is the most the \code{intoverpop='augmented'}
#' filter can identify: under that route a random effect on a variance cell
#' (DIFFUSION, MANIFESTVAR) reaches the likelihood only through the predicted
#' covariance, so the filter never updates its carrier state, and the data
#' determines its covariance with the mean-affecting effects but not the split
#' of that covariance into a standard deviation and correlations. \code{'auto'}
#' estimates exactly the part that is determined, which costs no likelihood, and
#' it is a no-op on any model where every varying parameter reaches the mean.
#' When it does reduce the rank it says so, naming the parameters affected.
#'
#' \code{NA} leaves the covariance unrestricted, as in versions before this
#' argument existed. The extra parameters are then estimated but not identified:
#' their reported values are one arbitrary point on a ridge, and their intervals
#' are not trustworthy in either direction.
#'
#' A whole number below the \code{'auto'} value is an explicit
#' **approximation**. It describes the individual differences with fewer
#' dimensions than the data supports, which lowers the likelihood and
#' \strong{distorts the parameters it retains} -- the fit is joint, so nothing
#' holds the retained covariances fixed while the rest is squeezed into them.
#' Useful for parsimony, or for speed in high dimensions, and not otherwise.
#'
#' The population covariance is \code{Sigma = [[S, S b'], [b S, b S b']]}, for a
#' freely estimated \code{S} over the basis effects and regression coefficients
#' \code{b} for the rest; each regressed effect has no variance independent of
#' the basis. \code{summary()} reports the standard deviations and correlations
#' this implies, with a note saying which of them follow from the structure
#' rather than being estimated.
#'
#' Applies under \code{intoverpop='augmented'}, \code{'laplace'} and
#' \code{'none'}, and requires \code{backend='julia'} -- but the **default only
#' applies to** \code{'augmented'}. What it means differs between them: on the
#' augmented route the coordinates it removes cannot be identified, so removing
#' them costs no likelihood; under \code{'laplace'} they are identified, and
#' removing them is an approximation that on one 250-subject design cost 48 log
#' likelihood units. So off the augmented route it has to be asked for
#' explicitly, and the message then says which of the two it is doing.
#'
#' May also be stated on the model, as \code{model$poprank <- 2}; an argument
#' here wins over that.
#' @param intoverpop how to handle declared individual differences. If 'auto',
#' set to TRUE if optimizing and FALSE if using hmc -- except when a grouping
#' level above the subject varies (see \code{id} in \code{\link{ctModel}}),
#' which only 'laplace' can integrate out, so 'auto' resolves to that.
#' if TRUE, integrates over population distribution of parameters rather than full sampling.
#' Allows for optimization of non-linearities and random effects, via state expansion.
#' 'augmented' names that state-expansion method explicitly. Individual
#' variation on a DIFFUSION or MANIFESTVAR parameter is only partially
#' identified under 'augmented' -- the data determines that effect's covariance
#' with the other random effects but not the split of it into a standard
#' deviation and correlations -- and \code{ctFit} warns when it sees one.
#' 'laplace' instead
#' integrates the random effects out subject by subject with a Laplace
#' approximation, leaving each subject's filtered state space at its
#' single-subject size, so the cost of a random effect stops growing cubically
#' with the number of varying parameters. 'laplace' requires
#' \code{backend='julia'} and works with either \code{optimize=TRUE} (Laplace
#' maximum likelihood) or \code{optimize=FALSE} (NUTS over the Laplace
#' marginal). It is exact whenever the
#' varying parameters enter the state mean linearly; elsewhere it is an
#' approximation, and \code{summary()} says so.
#' \code{FALSE} is the other route, and the one \code{'auto'} chooses when
#' \code{optimize=FALSE}: the individual parameters are sampled rather than
#' integrated over, so HMC targets the joint posterior over the population
#' parameters and every subject's random effects. That is exact whatever the
#' model, and its dimension grows with the number of subjects rather than
#' staying at the parameter count. \code{TRUE} and \code{FALSE} may be given
#' in place of the character forms above.
#' @param sameInitialTimes if TRUE, include an empty observation for every subject that has no observation
#' at the earliest observation time of the dataset. This ensures that the T0MEANS occurs for every subject at the same time,
#' rather than just at the earliest observation for that subject. Important when modelling trends over time, age, etc.
#' @param plot if TRUE, for sampling, a Shiny program is launched upon fitting to interactively plot samples.
#' May struggle with many (e.g., > 5000) parameters. For optimizing, various optimization details are plotted -- in development.
#' With \code{backend='julia'} the trace is plotted once the fit returns
#' rather than during it: a julia fit is a single blocking call into the
#' engine, so there is no point at which R could draw anything while it runs.
#' For genuinely live output use \code{optimcontrol$callback}.
#' @param derrind deprecated, latents involved in dynamic error calculations are determined automatically now.
#' @param optimize if TRUE, use \code{\link{stanoptimis}} function for maximum a posteriori / importance sampling estimates,
#' otherwise use the HMC sampler from Stan, which is (much) slower, but generally more robust for complex individual differences.
#' When \code{optimize=FALSE}, the stored point estimate (\code{stanfit$rawest}) is the per-parameter
#' median of the posterior draws; the julia backend's sampled point estimate (see \code{\link{ctSample}})
#' is the per-parameter mean instead.
#' @param optimcontrol list of parameters sent to \code{\link{stanoptimis}}
#' governing optimization. It is the only optimizer control list: the julia
#' backend's \code{backendcontrol} was merged into it, and every name below
#' means the same thing on both backends.
#'
#' Stopping rules, honoured by both: \code{tol} (the objective stops changing),
#' \code{g_tol} (l2 norm of the gradient), \code{x_tol} (size of the parameter
#' step), \code{maxiter}, and \code{lbfgs_memory} for how many curvature pairs
#' L-BFGS keeps. \code{initsd} is the sd of the random start on both.
#' \code{estonly}, \code{carefulfit}, \code{finishsamples}, \code{uncertainty},
#' \code{uncertaintyDraws} and \code{uncertaintyControl} also work on both.
#'
#' The \emph{defaults} are each backend's own, and are mirror images: stan stops
#' on the objective (\code{tol = 1e-8}, \code{g_tol} off), julia on the gradient
#' (\code{g_tol = 1e-8}, \code{tol} off). Leaving a name unset keeps that
#' backend's default; setting one is honoured by both.
#'
#' What is left is a capability one backend has and the other does not.
#' \code{backend='stan'} alone has the stochastic-gradient optimizer
#' (\code{stochastic}, \code{nsubsets}, \code{subsamplesize},
#' \code{lproughnesstarget}, \code{stochasticTolAdjust}, \code{parsteps}) and
#' the stall-and-restart pass (\code{stallretries}, \code{stalltol});
#' \code{backend='julia'} alone has \code{gradient}, \code{datastart},
#' \code{callback}, \code{saveEffects}, \code{progress} and
#' \code{tipredMissingIncludeOutcome}. A \emph{value} asking for a capability
#' the chosen backend does not have is refused by name before anything else
#' happens; a value that describes what it already does is simply accepted, so
#' \code{stochastic=FALSE} works on julia and \code{gradient='adjoint'} works on
#' stan.
#' With \code{backend='julia'}, \code{optimcontrol$gradient} selects the
#' gradient method: \code{'adjoint'} (reverse mode, the default) or
#' \code{'forward'} (ForwardDiff). Both compute the same gradient; 'adjoint'
#' costs the same regardless of the number of free parameters, so it is
#' dramatically faster for larger models and marginally slower for very small
#' ones.
#' \code{optimcontrol$datastart} (\code{backend='julia'}, default \code{TRUE})
#' takes the starting values for the diagonals of \code{DRIFT},
#' \code{DIFFUSION}, \code{T0VAR} and \code{MANIFESTVAR} from the data
#' instead of from the fixed point every fit used to start at. That point put
#' \code{DIFFUSION} at 6.93 and the variances at 3.47 whatever the data were
#' measured in, and starting far above the data's scale lets the optimiser
#' collapse the latent process to zero variance and call the whole signal
#' measurement error. Each process's scale is read off the indicators loading
#' on it -- on the scale its link puts them, so a count is read through
#' \code{log1p} and a binary or ordinal indicator, whose logit fixes the
#' latent scale by itself, contributes none -- and each cell's own transform is
#' then inverted numerically to get the raw value. Everything is clamped, and a
#' cell that cannot be solved keeps the old default.
#'
#' Measured on a three-indicator factor model with the data rescaled, over
#' four starting-jitter seeds each: at x0.01 the fit converged 0 times out of
#' 4 with this off and 4 out of 4 with it on, reaching 4565.60 against
#' 4203.69-4237.21 and recovering the first free loading at 0.808 against
#' -844 to -1206, for a generating value of 0.8. At x1 the two are identical.
#' At x100 it is neutral: that model is near a basin boundary at that scale
#' and converges twice in four either way, decided by the jitter on the
#' parameters this does not set rather than by the derived ones. Supplied
#' \code{inits} are never overridden, and \code{optimcontrol$datastart =
#' FALSE} restores the fixed start.
#'
#' \code{optimcontrol$carefulfit} works for \code{backend='julia'} as it does
#' for Stan: when \code{priors=FALSE}, a rough first pass is run \emph{with}
#' ctsem's \code{normal(0,1)} raw priors to obtain starting values, and the
#' likelihood is then maximised from there. It defaults to \code{TRUE}, capped
#' at 10 iterations, and is skipped when \code{inits} are supplied. Set
#' \code{optimcontrol$carefulfit = FALSE} to switch it off or to a number to
#' choose the cap.
#'
#' It does not make fits faster: measured over 800 fits, total iterations came
#' to 0.91-1.14 times a plain fit at ten prior iterations, 1.09-1.50 at twenty
#' and 1.47-1.64 at forty. It is on by default for where the fit lands rather
#' than how quickly it gets there. Over 720 fits judged on a random-effect SD
#' with a true value of 0.5, a cap of 10 gave an RMSE of 0.190 and a worst
#' error of 0.547, against 0.478 and 6.73 with the pass switched off. A longer
#' pass is not a safer one -- a cap of 20 scored 0.550 and 7.99, worse than not
#' warming up at all -- because the prior pass pulls the start toward the prior
#' mode and past about ten iterations that is what it hands the likelihood.
#'
#' \code{fit$estimate$carefulfit} records whether the pass ran, and
#' \code{$carefulfit_iterations} how long it was allowed.
#' With \code{backend='julia'}, \code{optimcontrol$callback} is a function
#' called while the fit runs, with \code{(iteration, total, objective,
#' gradient_norm)}. It is for a front end that wants to draw progress live:
#' the engine calls it on a time cadence rather than once per iteration,
#' because a callback costs about half a millisecond through the Julia
#' bridge, and always once more at the end. An error inside it disables it
#' and warns, leaving the fit unaffected. If output after the fit is enough,
#' \code{fit$trace} holds every iteration and \code{\link{ctTracePlot}}
#' draws it.
#' \code{backend='julia'} also finishes by estimating uncertainty, as the stan
#' backend does, and reads the same \code{stanoptimis} control names for it:
#' \code{uncertainty} (default \code{'hessian'}), \code{uncertaintyDraws},
#' \code{finishsamples}, and \code{uncertaintyControl}. Set
#' \code{optimcontrol$estonly = TRUE} for point estimates only.
#'
#' \code{optimcontrol$stallretries} and \code{optimcontrol$stalltol} govern what happens
#' when the stan optimizer stops somewhere that is not a maximum. A rough likelihood can
#' leave it unable to find any improving step, which it reports as convergence -- returning
#' the starting values with the uncertainty computed about them. The gradient where it
#' stopped separates the two cases by orders of magnitude, so the fit is restarted from
#' fresh values (twice by default) when it exceeds \code{stalltol} per data point, and warns
#' if it still ends up there.
#'
#' @param nopriors deprecated, use priors argument. logical. If TRUE, any priors are disabled -- sometimes desirable for optimization.
#' @param priors if TRUE, priors are included in computations, otherwise specified priors are ignored.
#' @param iter \strong{Deprecated} -- use \code{sampleControl$iter}. Still
#' honoured, with a warning.
#' @param inits either character string 'optimize, NULL, or vector of (unconstrained)
#' parameter start values, as returned by the rstan function \code{rstan::unconstrain_pars}, or the parameter values
#' found in a ctsem fit object \code{myfit$stanfit$rawest} (or \code{$rawposterior}) for instance.
#' @param cores number of cpu cores to use. A positive integer, or 'maxneeded' for
#' as many as available minus one (capped at the number of chains on the stan
#' backend, uncapped on julia, whose parallelism is over subject chunks). Defaults
#' to \code{getOption("mc.cores", 2)}. More cores are generally faster when
#' \code{optimize=TRUE}. Note that a julia fit at \code{cores > 1} is not
#' reproducible to the last decimal, because the chunk tuner times candidate
#' splits and the timings vary; use \code{cores = 1} for a before-and-after
#' comparison.
#'
#' On julia, \code{cores} cannot exceed the Julia session's thread count, which
#' Julia fixes when the session starts. A fit asked for more than the session
#' has runs at the count it has and says so;
#' \code{\link{ctJuliaSetup}(threads = n, force = TRUE)} restarts the session
#' wider, \code{ctJuliaStatus()$threads} reports what it currently has, and
#' \code{options(ctsem.julia.restart = TRUE)} has a fit restart it for itself
#' when it is short. Starting the engine before the first fit -- which any
#' script that warms it up front does -- otherwise pins every later fit to one
#' thread.
#' @param backend Either 'stan' (the default) or 'julia'. The julia backend is a
#' separate maximum-likelihood engine with the same model definitions and the
#' same summaries; it takes its own reverse-mode gradient, supports
#' \code{intoverpop='laplace'} for random effects, and can be sampled afterwards
#' with \code{\link{ctSample}}. It needs a working Julia -- see
#' \code{\link{ctJuliaSetup}} and \code{\link{ctJuliaInstall}}.
#' @param sampleControl Used when \code{optimize=FALSE}: a list holding
#' everything about how to sample. \code{iter} (default 1000) counts warmup and
#' sampling together, \code{warmup} (200, or half of \code{iter} if that is
#' less) how much of it is
#' discarded, \code{draws} the post-warmup count directly -- given, it wins and
#' \code{iter} is not consulted -- \code{chains} (2) how many chains,
#' \code{seed} (20260828), \code{saveEffects} whether individual random-effect
#' draws are kept, and \code{processes} (TRUE) whether the chains get their own
#' R processes.
#'
#' For \code{backend='julia'} it also carries the sampler's own settings:
#' \code{maxdepth}/\code{max_treedepth} (default 10),
#' \code{target_accept}/\code{adapt_delta} (0.8), \code{maxdelta} (1000),
#' \code{init_scale} (1), \code{adapt_metric} (FALSE), \code{adapt_effects}
#' (FALSE), and the effective-sample-size target that decides when a run stops:
#' \code{minESS} (200, the size the worst parameter must reach),
#' \code{rhatTarget} (1.01), \code{meanESS}, \code{maxDraws} and
#' \code{settleTol} -- all documented in full under \code{sampleControl} in
#' \code{\link{ctSample}}. A name the sampler does not read is an error rather
#' than ignored, because a name the list drops silently costs a whole run.
#' Given \code{minESS} or \code{meanESS}, \code{iter} becomes the budget the
#' run may take rather than the count it must: it stops as soon as the target
#' is met.
#'
#' For \code{backend='stan'} the sampler settings are rstan's, and are passed
#' to \code{\link[rstan]{stan}}'s own \code{control} argument.
#' @param chains \strong{Deprecated} -- use \code{sampleControl$chains}. Still
#' honoured, with a warning.
#' @param control \strong{Deprecated} -- use \code{sampleControl}, which is the
#' same list under a name that says what it controls. Still honoured, with a
#' warning; where both name the same setting, \code{sampleControl} wins.
#' For \code{backend='stan'}, a list of arguments sent to \code{\link[rstan]{stan}} control argument,
#' regarding warmup / sampling behaviour. Unless specified, values used are:
#' list(adapt_delta = .8, adapt_window=5, max_treedepth=10, adapt_init_buffer=2, stepsize = .001).
#' For \code{backend='julia'}, the same argument instead carries the julia sampler's own settings:
#' \code{maxdepth}/\code{max_treedepth} (default 10), \code{target_accept}/\code{adapt_delta} (0.8),
#' \code{maxdelta} (1000), \code{init_scale} (1), \code{adapt_metric} (FALSE), \code{adapt_effects} (FALSE),
#' and the optional effective-sample-size target \code{minESS}, \code{meanESS}, \code{maxDraws},
#' \code{rhatTarget} (1.01) and \code{settleTol} -- all documented in full under \code{control} in
#' \code{\link{ctSample}} -- plus \code{warmup} (default half of \code{iter}), \code{seed} (default
#' 20260828) and \code{processes} (default TRUE), which \code{\link{ctSample}} instead takes as
#' separate named arguments.
#' @param nlcontrol List of non-linear control parameters.
#' \code{maxtimestep} must be a positive numeric,  specifying the largest time
#' span covered by the numerical integration. The large default ensures that for each observation time interval,
#' only a single step of exponential integration is used. When \code{maxtimestep} is smaller than the observation time interval,
#' the integration is nested within an Euler like loop.
#' Smaller values may offer greater accuracy, but are slower and not always necessary. Given the exponential integration,
#' linear model elements are fit exactly with only a single step.
#' \code{nsubsteps = 'auto'} (julia backend only) instead measures, at the starting values and again at the
#' optimum, how nonlinear each observation interval is and refines only the intervals that need it;
#' \code{substeptol} (default 0.01) is the largest acceptable linearisation error as a fraction of the
#' predicted state standard deviation, and \code{maxsubsteps} (default 64) caps an interval.
#' \code{maxtimestep} remains a ceiling on the step. The choice is reported in \code{fit$estimate$substeps}.
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
#' @param stationary Logical. If TRUE, T0VAR and T0MEANS input matrices are ignored,
#' the parameters are instead fixed to long run expectations. More control over this can be achieved
#' by instead setting parameter names of T0MEANS and T0VAR matrices in the input model to 'stationary', for
#' elements that should be fixed to stationarity.
#' @param forcerecompile logical. For development purposes.
#' If TRUE, stan model is recompiled, regardless of apparent need for compilation.
#' @param saveCompile if TRUE and compilation is needed / requested, writes the stan model to
#' the parent frame as ctsem.compiled (unless that object already exists and is not from ctsem), to avoid unnecessary recompilation.
#' @param savescores Logical. If TRUE, output from the Kalman filter is saved in output. For datasets with many variables
#' or time points, will increase file size substantially.
#' @param savesubjectmatrices Logical. If TRUE, subject specific matrices are saved --
#' only relevant when either time dependent predictors or individual differences are
#' used. Can increase memory usage dramatically in large models, and can be computed after fitting using ctExtract
#' or ctSubjectPars .
#' @param saveComplexPars Logical. If TRUE, also save rowwise output of any complex parameters specified,
#' i.e. combinations of parameters, functions and states.
#' @param gendata Logical -- If TRUE, uses provided data for only covariates and a time and missingness structure, and
#' generates random data according to the specified model / priors.
#' Generated data is in the $Ygen subobject after running \code{extract} on the fit object.
#' For datasets with many manifest variables or time points, file size may be large.
#' To generate data based on the posterior of a fitted model, see \code{\link{ctGenerateFromFit}}.
#' @param compileArgs List of arguments to pass to \code{\link[rstan]{stan_model}} for compilation of the Stan model.
#' @param ... additional arguments to pass to \code{\link[rstan]{stan}} function.
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
#'  \code{'laplace'} or \code{'none'}. A call that gives the same \code{input}
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
#'   newdat <- ctGenerate(ctmodelobj = gm,n.subjects = 1,burnin = 2,
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
  nlcontrol = list(), nopriors=NA, priors=FALSE, chains=2,
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
  .ctsample_resolved <- .ctSampleControlResolve(sampleControl,
    given = names(match.call()), iter = iter, chains = chains,
    control = control)
  iter <- .ctsample_resolved$iter
  chains <- .ctsample_resolved$chains
  control <- .ctsample_resolved$control

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
  backend <- match.arg(backend)
  # Whether `poprank` was asked for or merely defaulted. Taken here because
  # `missing()` has to be evaluated before the argument is touched, and it
  # decides whether an inapplicable rank is an error or a no-op.
  poprankexplicit <- !missing(poprank)
  # Before any data preparation and before the Julia install prompt: a control
  # name the chosen backend cannot honour is a mistake to report immediately,
  # not after a wait.
  .ctFitCheckControls(optimcontrol, backend)
  if(backend %in% 'julia') {
    .ctJuliaUnsupported(ctstanmodel, optimize=optimize, priors=priors,
      intoverpop=intoverpop, gendata=gendata,
      stanmodeltext=stanmodeltext, compileArgs=compileArgs,
      forcerecompile=forcerecompile, optimcontrol=optimcontrol)
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

  if(!is.na(nopriors)){
    warning('nopriors argument is deprecated, use priors argument in future')
    priors <- !nopriors
  }

  if(any(!is.na(derrind))) warning('derrind argment is deprecated, computed automatically now')

  datalong <- data.frame(datalong)

  if(!ctstanmodel$timeName %in% colnames(datalong) && !ctstanmodel$continuoustime) {
    dtable <- data.table(datalong)
    dtable[,.ObsCount:=1:.N,by=ctstanmodel$id]
    datalong[[ctstanmodel$timeName]] <- dtable[['.ObsCount']]
    rm(dtable)
  }

  datalong <- datalong[order(datalong[[ctstanmodel$subjectIDname]],datalong[[ctstanmodel$timeName]]),] #sort by subject, time.

  datavars <- c(ctstanmodel$timeName,ctstanmodel$subjectIDname, ctstanmodel$manifestNames,ctstanmodel$TDpredNames,ctstanmodel$TIpredNames)
  sapply(datavars,function(x){
    if(!x %in% colnames(datalong)) stop(paste0(x,' column not found in data!'))
    if(!x %in% ctstanmodel$subjectIDname){ #if not an id column
      if(any(!is.numeric(as.numeric(datalong[!is.na(datalong[,x]),x])))) stop(x ,' column contains non-numeric data!')
    }
  })

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

  if(!is.null(ctm$TIpredAuto) && ctm$TIpredAuto %in% c(1L,TRUE)){ #if auto tipred, set all effects to true
    for(tip in ctm$TIpredNames){
      ctm$pars[[paste0(tip,'_effect')]] <- 'TRUE'
    }
  }

  if(optimize && !priors) message("Maximum likelihood estimation requested")
  # `optimcontrol$is` is refused above, so it is NULL by the time we get here
  # and the importance-sampling wording this used to choose is unreachable.
  # Importance sampling is now `optimcontrol$uncertainty='is'`, which runs
  # after optimization rather than instead of it, so the estimation this
  # message describes is a posteriori either way.
  if(optimize && priors) message("Maximum a posteriori estimation requested")
  # Naming stan here was wrong for half the fits it described: with
  # backend='julia' the engine runs its own NUTS over the joint posterior of
  # parameters and random effects, and a user reading "Stan's NUTS sampler" on a
  # julia fit has no reason to believe the julia sampler ran at all.
  if(!optimize) message("Bayesian estimation via ",
    if(identical(backend,'julia')) "the julia engine's" else "Stan's",
    " NUTS sampler requested")


  ###stationarity
  if(stationary) {
    stop('Stationary option temporarily unavailable -- reductions needed to pass all CRAN checks')
    ctm$pars$param[ctm$pars$matrix %in% c('T0VAR','T0MEANS')] <- 'stationary'
    ctm$pars$value[ctm$pars$matrix %in% c('T0VAR','T0MEANS')] <- NA
    ctm$pars$indvarying[ctm$pars$matrix %in% c('T0VAR','T0MEANS')] <- FALSE
  }

  #collect individual stationary elements and update ctm$pars
  if(any(ctm$pars$param %in% 'stationary'))  stop('Stationary option temporarily unavailable -- reductions needed to pass all CRAN checks')
    ctm$t0varstationary <- as.matrix(rbind(ctm$pars[which(ctm$pars$param %in% 'stationary' & ctm$pars$matrix %in% 'T0VAR'),c('row','col')]))
    if(nrow(ctm$t0varstationary) > 0){ #ensure upper tri is consistent with lower
      for(i in 1:nrow(ctm$t0varstationary)){
        if(ctm$t0varstationary[i,1] != ctm$t0varstationary[i,2]) ctm$t0varstationary <- rbind(ctm$t0varstationary,ctm$t0varstationary[i,c(2,1)])
      }}
    ctm$t0varstationary = unique(ctm$t0varstationary) #remove any duplicated rows
    ctm$t0meansstationary <- as.matrix(rbind(ctm$pars[which(ctm$pars$param[ctm$pars$matrix %in% 'T0MEANS'] %in% 'stationary'),c('row','col')]))
    ctm$pars$value[ctm$pars$param %in% 'stationary'] <- -99 #does this get inserted?
    ctm$pars$indvarying[ctm$pars$param %in% 'stationary'] <- FALSE
    ctm$pars$transform[ctm$pars$param %in% 'stationary'] <- NA
    ctm$pars$param[ctm$pars$param %in% 'stationary'] <- NA


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
    # them. The very next check warns that `intoverstates=TRUE` is required for
    # sensible optimization -- so the documented shortcut put a user straight
    # into the state the code itself calls unreliable, under the default
    # `optimize=TRUE`. Setting `manifesttype` directly never did that, and the
    # linearised measurement handles binary indicators with the filter intact:
    # on a three-indicator model it recovers a generating drift of -0.3 as
    # -0.279 and a diffusion of 0.8 as 0.681, both intervals containing the
    # truth.
    ctm$manifesttype[] <- 1
  }

  recompile <- FALSE
  if(!optimize && !priors){
    message('HMC sampling requested, but priors disabled -- are you sure? consider setting priors=TRUE')
    # !priors <- FALSE
  }
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
  if(is.character(intoverpop)){
    intoverpop <- match.arg(intoverpop[1], c('auto','augmented','laplace'))
    if(intoverpop %in% 'auto'){
      intoverpop <- isTRUE(optimize) && .ctAnyVarying(ctm)
      # The augmented layout gives a carrier state to every `indvarying` cell
      # and knows nothing about the columns a grouping level uses, so a model
      # with effects above the subject has one route rather than two and 'auto'
      # has to take it. Resolving to 'augmented' here would silently fit a
      # model without the study effect that was asked for.
      if(intoverpop && .ctAnyVarying(ctm, .ctOuterVaryingColumns(ctm))){
        intoverpopmethod <- 'laplace'
        intoverpop <- FALSE
      }
    } else {
      intoverpopmethod <- intoverpop
      intoverpop <- identical(intoverpopmethod,'augmented')
    }
  }
  intoverpop <- isTRUE(intoverpop)
  if(intoverpop) intoverpopmethod <- 'augmented'

  if(identical(intoverpopmethod,'laplace')){
    if(!backend %in% 'julia') stop(
      "intoverpop='laplace' requires backend='julia'; the generated Stan model ",
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
  # Two routes to the same restriction, because the two have different
  # machinery to hang it on. Under `'augmented'` a basis effect has a carrier
  # state and a regressed cell references `state[j]`, so the rewrite has to
  # follow `.ctModelIntOverPop()`. Under `'laplace'` -- and `'none'`, which
  # prepares the same structure without integrating -- there are no carrier
  # states, so the basis effects move into PARS and the regressed cells
  # reference them as parameters. Both land before the second
  # `ctModelStatesAndPARS()` call below, which is what turns the new labels into
  # `PARS[r,c]` references.
  #
  # What the restriction *means* differs between them, and only the message says
  # so: on the augmented route it removes coordinates the filter cannot see and
  # costs no likelihood, while under laplace those coordinates are identified
  # and removing them is an approximation. Same structure, different claim.
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
  if(!poprankexplicit && !is.null(ctm[['poprank']])) poprank <- ctm[['poprank']]

  popregression <- NULL
  if(!(length(poprank)==1 && is.na(poprank))){
    if(!identical(backend,'julia')){
      if(poprankexplicit) stop("poprank requires backend='julia'.", call.=FALSE)
    } else if(!intoverpop && !identical(intoverpopmethod,'laplace') &&
        !any(ctm$pars$indvarying[is.na(ctm$pars$value)])){
      if(poprankexplicit) stop(
        "poprank restricts the population covariance, so it needs a model with ",
        "individually varying parameters.", call.=FALSE)
    } else if(!intoverpop && !poprankexplicit){
      # Not by default off the augmented route, and this is the whole reason
      # the two are distinguished. On the augmented route the coordinates
      # `'auto'` removes cannot be identified, so removing them costs nothing
      # and is a good default. Under laplace they *are* identified, and
      # measured on a 250 x 50 design the same restriction costs 48 log
      # likelihood units and takes the fit to the boundary -- basis sd to zero
      # with the coefficient to -612. A default that does that to a user who
      # chose the route precisely because it identifies these things would be
      # indefensible, so here it has to be asked for.
      popregression <- NULL
    } else {
      popregression <- .ctPopRegressionSpec(ctm$pars, poprank,
        explicit=poprankexplicit, model=ctm)
      if(!is.null(popregression)) ctm <- .ctPopRegressionDemote(ctm, popregression)
    }
  }
  if(intoverpop)   ctm <- .ctModelIntOverPop(ctm) #extend system matrices for individual differences
  if(!is.null(popregression)){
    ctm <- if(intoverpop) .ctPopRegressionRewrite(ctm, popregression) else
      .ctPopRegressionRewriteParameters(ctm, popregression)
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
  if(any(ctm$manifesttype %in% 2)) .ctDataCategories(datalong, ctm)
  if(any(ctm$manifesttype %in% 3)) .ctDataCounts(datalong, ctm)
  if(any(ctm$manifesttype %in% 4)) .ctDataCensored(datalong, ctm)

  if(any(ctm$manifesttype > 0)){ #if any non continuous variables, (with free parameters)...
    # Censored variables are excluded: a censored observation is Gaussian
    # within its limits, so its measurement standard deviation is the scale of
    # the whole thing and has to stay free. Every other non-Gaussian type
    # supplies its own randomness through the link and would be adding noise on
    # top of noise.
    deterministic <- which(ctm$manifesttype > 0 & ctm$manifesttype != 4)
    errfix <- which(ctm$pars$matrix %in% 'MANIFESTVAR' &
        (ctm$pars$row %in% deterministic |
            ctm$pars$col %in% deterministic) &
        is.na(suppressWarnings(as.numeric(
          ctm$pars$value))))

    if(length(errfix) > 0){
      message('Fixing any free MANIFESTVAR parameters for binary / ordinal indicators to deterministic calculation')
      ctm$pars$value[errfix] <- 1e-5
      ctm$pars[errfix,c('param','transform','multiplier','offset','meanscale','inneroffset','sdscale')] <- NA
      ctm$pars$indvarying[errfix] <- FALSE
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
    if(backend %in% 'julia'){
      message('Binary and ordinal indicators are integrated rather than ',
        'linearised on this backend: the observation is taken against the ',
        'predicted state by quadrature, so DRIFT and DIFFUSION are estimated ',
        'without the linearisation bias the stan path carries.')
      if(any(ctm$manifesttype %in% 2)) message(
        'Ordinal thresholds are reported as the first threshold followed by ',
        'the gap to each subsequent one, which is what keeps them ordered; ',
        'cumulate them to read the thresholds themselves.')
    } else {
      message('Binary indicators use a linearised (moment-matched Gaussian) ',
        'measurement update on the stan backend, which makes DRIFT and ',
        'especially DIFFUSION unreliable for a latent seen only through them. ',
        'Over 40 replications with a true diffusion of 0.8, RMSE was 0.45 with ',
        '3 indicators and 0.23 with 30, against 0.13 and 0.04 for ',
        "backend='julia', which integrates the observation instead. Prefer ",
        'the julia backend for binary data, or treat process noise as ',
        'indicative.')
    }

    binaryrows <- which(ctm$pars$matrix %in% 'MANIFESTVAR' &
        ctm$pars$row %in% which(ctm$manifesttype > 0) &
        ctm$pars$row == ctm$pars$col)
    stated <- binaryrows[!is.na(ctm$pars$value[binaryrows]) &
        abs(ctm$pars$value[binaryrows]) > 1e-4]
    if(length(stated)){
      warning('MANIFESTVAR is fixed to a non-zero value for categorical indicator',
        if(length(stated) > 1) 's ' else ' ',
        paste(ctm$manifestNames[ctm$pars$row[stated]], collapse=', '),
        '. A categorical indicator gets its randomness from its measurement ',
        'link -- the Bernoulli link for binary, the cumulative logit for ',
        'ordinal -- so ',
        'this adds measurement noise on top of it. Set it to 0 unless that is ',
        'meant.', call.=FALSE)
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
  # The rank actually used, not the argument: 'auto' resolves to a number, and a
  # model the restriction did not apply to reports NA whatever was asked for.
  argsresolved$poprank <- if(is.null(popregression)) NA_integer_ else
    as.integer(popregression$rank)
  argsresolved$priors <- as.logical(priors)
  argsresolved$optimize <- isTRUE(optimize)
  argsresolved$intoverstates <- isTRUE(intoverstates)

  if(backend %in% 'julia') {
    .ctJuliaUnsupported(ctm, optimize=optimize, priors=priors,
      intoverpop=intoverpop, gendata=gendata,
      stanmodeltext=stanmodeltext, compileArgs=compileArgs,
      forcerecompile=forcerecompile, intoverstates=intoverstates,
      optimcontrol=optimcontrol)
    # `optimize` and `intoverpop` are orthogonal here. `intoverpop` says which
    # random effects are integrated out and how; `optimize` says whether the
    # remaining parameters are maximised or sampled. Every combination is
    # meaningful for this backend:
    #
    #   optimize  intoverpop     what runs                 sampled dimension
    #   TRUE      'laplace'      Laplace ML                --
    #   TRUE      TRUE           augmented ML              --
    #   FALSE     'laplace'      NUTS, Laplace marginal    npar
    #   FALSE     TRUE           NUTS, filter marginal     npar
    #   FALSE     FALSE          NUTS over parameters      npar + effects
    #                            *and* effects
    #
    # Only 'laplace' is julia-only; the guard for that is above. `'none'` is
    # the third route the engine needs, and it prepares the Laplace structure
    # -- which is what says *which* parameters vary -- without integrating.
    juliaintoverpop <- if(identical(intoverpopmethod,'laplace')) 'laplace' else
      if(intoverpop) 'augmented' else
        if(!optimize && .ctAnyVarying(ctm)) 'none' else
          'augmented'
    juliafit <- .ctFitJuliaBackend(datalong=datalong, model=ctm, prepared_data=standata, inits=inits,
      cores=cores, optimcontrol=optimcontrol,
      verbose=verbose, fit=fit, priors=priors, optimize=optimize,
      chains=chains, iter=iter, control=control,
      intoverpop=juliaintoverpop, intoverstates=intoverstates)
    # Replaces whatever narrower `$args` the julia backend built internally
    # (it only ever had the resolved settings, and not all of them) with the
    # two-sublist form every fit now carries -- see the `argsresolved` comment
    # above. This also reaches the `fit=FALSE` case: `.ctFitJuliaBackend()`
    # returns the prepared model spec then, unclassed by `$args` before, and
    # assigning a list element to it here does not disturb its class.
    juliafit$args <- list(input = args, resolved = argsresolved)
    # `$data`/`$standata` mean the same thing on both backends: `$standata` is
    # the prepared data with the 99999 missing-value sentinel intact, `$data`
    # is the same thing with that sentinel replaced by `NA` in `$Y` (and,
    # replicating the stan path exactly -- see the identical two lines below --
    # a no-op attempt at `$tipreds`, which is not a field of this list; the
    # real time-invariant predictor data lives in `$tipredsdata` and keeps its
    # sentinel in both copies). `standata` was already computed above,
    # unconditionally, before backend dispatch -- .ctPrepareData() runs for julia
    # too, purely to prepare `prepared_data` for `.ctFitJuliaBackend()` -- so
    # attaching it here costs nothing further and is not a second computation.
    standataout <- standata
    standataout$Y[standataout$Y==99999] <- NA
    standataout$tipreds[standataout$tipreds==99999] <- NA
    juliafit$standata <- standata
    juliafit$data <- standataout
    # `plot` draws the trace *after* the fit here, not during it.
    #
    # The Stan path can plot live because it writes sample files a second
    # process reads. A julia fit is one blocking call into the engine, so R
    # cannot draw anything until it returns -- there is no point in this
    # function where a live plot could be made. What is drawn is the same
    # information, recorded every iteration and handed back on the fit; for
    # genuinely live output, `optimcontrol$callback` is called while the fit
    # runs and can draw whatever it likes.
    if(isTRUE(fit) && !identical(plot, FALSE) && !is.null(juliafit$trace)) {
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
        parswithtipreds <- sort(unique(ms$param[ms$param >0 & ms$when %in% c(0,-1) & ms$copyrow < 1]))
        parswithtipreds<-parswithtipreds[apply(stanfit$standata$TIPREDEFFECTsetup,1,sum)>0]
        ms$tipred[ms$param >0 & ms$when %in% c(0,-1) & ms$copyrow < 1 & ms$param %in% parswithtipreds] <- 1L
        ctm$modelmats$matsetup <- ms
      }
    }

    # if(is.na(STAN_NUM_THREADS)) Sys.unsetenv('STAN_NUM_THREADS') else Sys.setenv(STAN_NUM_THREADS = STAN_NUM_THREADS) #reset sys env
  } # end if fit==TRUE
  #convert missings back to NA's for data output
  standataout<-standata
  standataout$Y[standataout$Y==99999] <- NA
  standataout$tipreds[standataout$tipreds==99999] <- NA
  # standataout <- utils::relist((standataout),skeleton=standata)

  setup=list(recompile=recompile,idmap=standata$idmap,matsetup=ctm$modelmats$matsetup,matvalues=ctm$modelmats$matvalues,
    popsetup=ctm$modelmats$matsetup[ctm$modelmats$matsetup$when %in% c(0,-1) & ctm$modelmats$matsetup$param > 0,],
    popvalues=ctm$modelmats$matvalues[ctm$modelmats$matsetup$when %in% c(0,-1) & ctm$modelmats$matsetup$param > 0,],
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

