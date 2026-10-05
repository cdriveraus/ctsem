# The optimiser settings, as a help page of their own.
#
# ?ctFit carried all of this in its `optimcontrol` entry, about 380 lines with
# the measurements behind each default, which made one setting hard to find.
# The page below says what each name does, on which backend, and its default.
# Why each default is what it is lives beside the code that reads it: the
# carefulfit studies at `.ctJuliaOptimiseFit()`, the datastart measurements in
# R/ctDataStart.R, the Laplace correction at
# `.ctLaplaceContinue()`, the floor at `.ctJuliaLaplaceFloorDefault`, SAEM in
# `saem.jl` and review/SAEM-*.md. The vocabulary itself -- which names are
# shared and which belong to one backend -- is `.ctOptimcontrolShared` and
# `.ctOptimcontrolSplit()` in ctFit.R, which this page must agree with.

#' Optimiser settings for ctFit
#'
#' The names \code{\link{ctFit}} reads from \code{optimcontrol}. A name means
#' the same thing on both backends. One a backend does not have is refused by
#' name, unless its value asks for nothing that backend does not already do
#' (\code{gradient = 'adjoint'} on stan, \code{nsubsets = 1} on julia), and a
#' misspelled name is an error. Leaving a name unset keeps that backend's own
#' default.
#'
#' @section Stopping:
#' \describe{
#'  \item{\code{tol}}{Both. Stop when the objective changes by less than this.
#'  Default \code{1e-8} on stan, off on julia.}
#'  \item{\code{g_tol}}{Both. Stop when the l2 norm of the gradient falls below
#'  this. Default \code{1e-8} on julia, off on stan.}
#'  \item{\code{x_tol}}{Both. Stop when the parameter step falls below this. Off
#'  by default.}
#'  \item{\code{maxiter}}{Both. Iteration cap. Default 1000 on julia, mize's on
#'  stan.}
#'  \item{\code{lbfgs_memory}}{Both. How many curvature pairs L-BFGS keeps.
#'  Default 20 on julia, mize's on stan.}
#'  \item{\code{stallretries}}{Both. Default 2. How often a fit that stopped
#'  somewhere that is not a maximum is tried again: stan restarts from fresh
#'  random values when the gradient exceeds \code{stalltol} per data point, and
#'  julia pulls the flat coordinates back and refits, keeping the refit only if
#'  it improves the objective. 0 turns it off.}
#' }
#'
#' @section Starting values:
#' \describe{
#'  \item{\code{initsd}}{Both. Default 0.01. The sd of the random start, on the
#'  raw scale.}
#'  \item{\code{datastart}}{Julia. Default TRUE. Starting values for the
#'  diagonals of DRIFT, DIFFUSION, T0VAR and MANIFESTVAR are read off the
#'  data's scale -- through the link, so a count is read through \code{log1p}
#'  and a binary or ordinal indicator contributes none -- rather than from one
#'  fixed point whatever the units. Supplied \code{inits} are never
#'  overridden.}
#'  \item{\code{carefulfit}}{Both. A rough first pass with ctsem's
#'  \code{normal(0,1)} raw priors on every coordinate, for starting values; the
#'  fit's own objective is then maximised from there. Skipped with
#'  \code{priors = TRUE}, which already has those priors, and when \code{inits}
#'  are supplied. Default TRUE on stan; on julia TRUE when any indicator is not
#'  Gaussian and FALSE otherwise, capped at 10 iterations, and a number sets the
#'  cap. A longer pass is not a safer one: it pulls the start toward the prior
#'  mode. \code{fit$optim$carefulfit} records whether it ran.}
#' }
#'
#' @section The optimiser:
#' \describe{
#'  \item{\code{stochastic}}{Both. Start with ctsem's stochastic-gradient
#'  optimiser and finish with L-BFGS; \code{'auto'} does so above 50
#'  parameters. Default TRUE on stan, FALSE on julia. On large, poorly
#'  conditioned models it can be much faster than L-BFGS alone, and it can
#'  carry a fit further toward a degenerate limit, such as a variance running
#'  to zero. \code{fit$optim$sgd_iterations} records the phase.}
#'  \item{\code{gradient}}{Julia. \code{'adjoint'} (reverse mode, the default)
#'  or \code{'forward'}. The same gradient; the adjoint's cost does not grow
#'  with the number of parameters.}
#'  \item{\code{batch}}{Julia. Default TRUE. Start on a random subset of the
#'  subjects and grow it as the fit needs more data; needs at least 80
#'  subjects or top-level groups.}
#'  \item{\code{lbfgs_diagonal}}{Julia. Default TRUE. L-BFGS's initial inverse
#'  Hessian takes a scale per parameter, learned from the curvature pairs
#'  (Gilbert and Lemarechal); FALSE uses one scale for all, which can be much
#'  slower on large models. Like \code{stochastic}, it can carry a fit further
#'  toward a degenerate limit.}
#'  \item{\code{lbfgs_gll}}{Julia. Default 0. A step is accepted against the
#'  worst of the last \code{W} objective values (Grippo, Lampariello and
#'  Lucidi), so it may cross a curved valley.}
#'  \item{\code{lbfgs_nonmonotone}}{Julia. Between 0 (the default) and 1. A step
#'  is accepted against a running average of recent objective values (Zhang
#'  and Hager).}
#'  \item{\code{initial_alpha}}{Julia. Default 0.1. The length of the first
#'  trial step, in raw units.}
#'  \item{\code{precondition}}{Julia. Default TRUE. FALSE drops the diagonal
#'  metric taken from the model's transforms.}
#'  \item{\code{saem}}{Julia, \code{intoverpop = 'laplace'}. Default FALSE.
#'  Start with SAEM: the random effects are sampled from their conditional
#'  distribution rather than integrated by the Laplace approximation, heading
#'  for the exact marginal posterior mode. The usual optimiser then polishes
#'  and certifies its point on the Laplace objective, so the estimate usually
#'  ends where a fit without SAEM does; what it can change is which maximum is
#'  reached when the Laplace objective has more than one. Not a better
#'  approximate estimate on its own -- for that, sample with
#'  \code{optimize = FALSE}. A number caps its iterations (TRUE is 10000).
#'  \code{fit$optim$saem_iterations}, \code{saem_settled}, \code{saem_trend},
#'  \code{saem_chains}, \code{saem_acceptance} and \code{saem_trace} record the
#'  phase.}
#'  \item{\code{saem_proposal}}{Julia. SAEM's proposal, \code{'rw'} or
#'  \code{'laplace'}.}
#'  \item{\code{restarts}}{Julia. Default 0. Random restarts for a fit that has
#'  still not converged after its corrections -- typically a likelihood with
#'  more than one basin -- each from the fit's own start plus normal noise of sd
#'  \code{restartsd} (default 1) on the raw scale, keeping the best if it is
#'  better. They run in worker processes when \code{cores > 1} and the fit took
#'  more than 30 seconds. Not when \code{inits} or \code{maxiter} were supplied.
#'  \code{fit$optim$restarts} records each start.}
#' }
#'
#' @section The finish and its certification (julia):
#' \describe{
#'  \item{\code{newton}}{Default TRUE. Newton steps on the exact Hessian once
#'  close, ending on the Hessian the certification and the standard errors
#'  use; it tries the direction of negative curvature at a saddle and checks
#'  what the directions with no curvature are worth. Under
#'  \code{intoverpop = 'laplace'} it runs only when the fit is certified. May
#'  name the curvature its steps use: \code{'exact'}, \code{'chord'} or
#'  \code{'subset'}. \code{fit$optim$newton_steps} records it.}
#'  \item{\code{certify}}{Default TRUE. The curvature check that certifies a
#'  fit has converged. \code{estonly} skips it too.}
#'  \item{\code{gaptol}}{Default \code{1e-6}. How much objective the curvature
#'  may still predict at the estimate before the fit is continued.}
#'  \item{\code{gapretries}}{Default 2. How often a fit the check finds short of
#'  a maximum is resumed; a resumed stage that stops gaining is stopped.}
#'  \item{\code{innergaptol}}{Stop once a step is predicted to gain less than
#'  this; by default a hundredth of the certification's tolerance.}
#'  \item{\code{overshoot}}{\code{'magnitude'} (the default),
#'  \code{'saturation'} or \code{'off'}: which directions are pulled back after
#'  a fit, to test whether it stopped at a maximum.}
#'  \item{\code{stallwindow}, \code{stallfraction}, \code{stallratio},
#'  \code{stallcooldown}, \code{stalltighten}, \code{stalltightenings}}{The
#'  stall check (defaults 80, 0.01, 0.001, 30, 0.1, 2): how many iterations of
#'  negligible progress end a fit, what share of its own progress counts as
#'  negligible, how much of its responsiveness a transform must have lost to
#'  count as the reason, and how long and how far the bar tightens after a fit
#'  is found slow rather than stuck.}
#'  \item{\code{escapesaturated}}{Default FALSE. Refit a fit stalled on a
#'  saturated transform from a zeroed boundary, keeping whichever wins.}
#'  \item{\code{escapepin}}{Default FALSE. Hold the escaped coordinates still
#'  while the rest re-optimises around them, then release them.}
#' }
#'
#' @section The Laplace route (julia, \code{intoverpop = 'laplace'}):
#' \describe{
#'  \item{\code{laplace_inner_maxiter}, \code{laplace_inner_tol}}{Defaults 200
#'  and \code{1e-10}. The inner solve that finds each unit's random-effect mode.
#'  Raise the first when the fit reports modes that did not converge; the
#'  tolerance bounds the inner gradient, floored at
#'  \code{1e-10 * (1 + abs(value))}.}
#'  \item{\code{laplace_floor}}{\code{'gated'} (the default) or \code{'total'}:
#'  how a unit whose likelihood is convex in its random effects is scored.
#'  \code{'gated'} scores a unit whose curvature has an eigenvalue below 0.7 by
#'  a three-point quadrature along that direction, handing over to
#'  \code{'total'}, which floors the log determinant at zero, between 0.2 and
#'  0.7. Recorded in \code{fit$laplace$floor}, and used by every function that
#'  rebuilds the objective; \code{fit$laplace$conditioning} counts such units.
#'  It applies too when sampling with \code{intoverpop = FALSE} and
#'  \code{sampleControl$placement = 'fit'}.}
#'  \item{\code{laplace_correct}}{\code{'quadrature'} (the default, or TRUE),
#'  \code{'step'} or FALSE. Corrects the estimate for the Laplace
#'  approximation's error by the adaptive Gauss-Hermite quadrature of
#'  \code{\link{ctLaplaceCheck}} (5 nodes per random effect), last, after the
#'  certification and the standard errors. Where the two agree to 0.01 in total
#'  -- every model whose effects enter linearly -- nothing changes. Elsewhere
#'  \code{fit$estimate$loglik} is the quadrature log likelihood (the Laplace one
#'  is kept in \code{loglik_laplace}). \code{'quadrature'} continues the fit on
#'  the quadrature objective for the units whose values differ; \code{'step'}
#'  takes up to three Newton steps toward its optimum. A unit with more than
#'  five effects keeps the Laplace term. Refused with \code{estonly}.
#'  \code{fit$laplace$correction} records what was done, and \code{print(fit)}
#'  says when the estimate moved.}
#' }
#'
#' @section Uncertainty:
#' \describe{
#'  \item{\code{uncertainty}}{Both. Default \code{'hessian'}. The method
#'  computed when the optimiser ends, as in \code{\link{ctFitUncertainty}}.
#'  \code{'sample'} is refused here: use \code{optimize = FALSE}, or
#'  \code{ctFitUncertainty(fit, 'sample')} on the fit.}
#'  \item{\code{uncertaintyDraws}, \code{finishsamples},
#'  \code{uncertaintyControl}}{Both. \code{ctFitUncertainty}'s \code{draws}
#'  (default \code{'auto'}), \code{finishsamples} (1000) and \code{control}.}
#'  \item{\code{estonly}}{Both. Default FALSE. Point estimates only: no Hessian,
#'  so no standard errors and no certification; under
#'  \code{intoverpop = 'laplace'} no Newton finish and no quadrature
#'  correction either.}
#' }
#'
#' @section Progress and output (julia):
#' \describe{
#'  \item{\code{callback}}{A function called while the fit runs, with
#'  \code{(iteration, total, objective, gradient_norm, parameters)}, where
#'  \code{parameters} is the raw vector the objective and gradient describe; one
#'  declaring only the first four is called with four. It is called on a time
#'  cadence and once more at the end. For live progress in a front end, or for
#'  checkpointing: \code{parameters} written from it and passed back as
#'  \code{inits} resume a fit, though the quasi-Newton history restarts. An
#'  error inside it disables it with a warning. \code{fit$optim$trace} holds
#'  every iteration after the fit, and \code{\link{ctTracePlot}} draws it.}
#'  \item{\code{progress}}{Force the progress line on or off.}
#'  \item{\code{tipredMissingIncludeOutcome}}{Default TRUE. Whether the outcome
#'  informs the imputation of missing TI predictors.}
#'  \item{\code{saveEffects}}{A sampling setting; read here for older scripts.
#'  Use \code{sampleControl$saveEffects}.}
#' }
#'
#' @section Stan only:
#' \describe{
#'  \item{\code{nsubsets}, \code{subsamplesize}, \code{lproughnesstarget},
#'  \code{stochasticTolAdjust}, \code{parsteps}}{Tuning of stan's
#'  stochastic-gradient optimiser and its stepwise phase.}
#'  \item{\code{stalltol}}{Default 0.01. The gradient per data point above which
#'  stan calls a fit stalled and restarts it (see \code{stallretries}).}
#' }
#'
#' @seealso \code{\link{ctFit}}, \code{\link{ctFitUncertainty}},
#' \code{\link{ctTracePlot}}.
#' @name ctOptimControl
NULL
