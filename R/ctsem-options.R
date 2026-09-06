# Options, in one place.
#
# Each of these was documented where it took effect, which is the right place to
# find it if you already know it exists and the wrong place to discover it.
# There was no way to answer "what can I configure?" short of reading the
# source.

#' Options ctsem reads
#'
#' Settings that change how ctsem behaves across a session, rather than
#' per-call. Set them with \code{options()}, typically once at the top of a
#' script or in a front end's startup.
#'
#' @section Progress and output:
#' \describe{
#'   \item{\code{ctsem.progress.overwrite}}{Whether progress reporting rewrites
#'     a single line with a carriage return (\code{TRUE}) or prints occasional
#'     separate lines (\code{FALSE}). Unset, ctsem detects it: an active Shiny
#'     session, a \code{sink()}, or a knitr chunk all capture output, where a
#'     carriage return is a character rather than a cursor movement. Set it
#'     explicitly if that detection is wrong for your front end.}
#'   \item{\code{ctsem.cluster.outfile}}{Where parallel workers send their
#'     output. Unset, they are silent. Setting it to \code{""} forwards worker
#'     stdout and stderr to the console, which is useful only when debugging a
#'     worker -- it also restores the startup chatter that made a two-core fit
#'     announce itself six times.}
#' }
#'
#' @section Per-call alternatives:
#' Progress and reporting are also controlled per fit, and those win over these:
#' \code{verbose} in \code{\link{ctFit}} (0 quiet, 1 reports progress, 2 keeps
#' every line rather than overwriting), \code{optimcontrol$callback} for a
#' function called while the fit runs, and \code{optimcontrol$progress} to
#' force progress on or off for one fit.
#'
#' @section Julia:
#' The Julia backend is mostly located through the environment rather than
#' options: \code{JULIA_BINDIR} names the binary directory,
#' \code{JULIA_NUM_THREADS} the thread count of the session that starts next,
#' and \code{ctJuliaSetup(project=)} points the session at a different engine
#' environment. See \code{\link{ctJuliaSetup}}.
#' \describe{
#'   \item{\code{ctsem.julia.restart}}{Whether a fit asking for more
#'     \code{cores} than the running Julia session has threads may restart the
#'     session to get them. \code{FALSE} by default: Julia fixes its thread
#'     count at process start, so the only way to raise it is to end the
#'     process, which discards the engine's compiled model shapes and costs
#'     seconds. Left \code{FALSE}, such a fit runs at the threads it has and
#'     says so.}
#' }
#'
#' @return Nothing; this page documents options rather than defining a
#'   function.
#' @name ctsem-options
#' @aliases ctsem-options
NULL
