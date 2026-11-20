# Options, in one place.
#
# Each of these was documented where it took effect, which is the right place to
# find it if you already know it exists and the wrong place to discover it.
# There was no way to answer "what can I configure?" short of reading the
# source.

#' Options ctsem reads
#'
#' Settings that change how ctsem behaves for a whole session. Set them with
#' \code{options()}, usually once at the top of a script or when a front end
#' starts.
#'
#' @section Progress and output:
#' \describe{
#'   \item{\code{ctsem.progress.overwrite}}{Whether progress reporting rewrites
#'     a single line (\code{TRUE}) or prints separate lines now and then
#'     (\code{FALSE}). If unset, ctsem decides from where the output goes. In
#'     a Shiny session, under \code{sink()} or in a knitr chunk the line is
#'     not rewritten. Set it if that choice is wrong for your front end.}
#'   \item{\code{ctsem.cluster.outfile}}{Where parallel workers send their
#'     output. If unset, they are silent. \code{""} sends their output to the
#'     console, which is useful when debugging a worker.}
#' }
#'
#' @section Per-call alternatives:
#' Progress and reporting can also be set for each fit, and these settings
#' take precedence: \code{verbose} in \code{\link{ctFit}} (0 quiet, 1 reports
#' progress, 2 keeps every progress line), \code{optimcontrol$callback} for a
#' function called while the fit runs, and \code{optimcontrol$progress} to
#' turn progress on or off for one fit.
#'
#' @section Julia:
#' The Julia backend is mostly configured through environment variables:
#' \code{JULIA_BINDIR} names the Julia binary directory, and
#' \code{JULIA_NUM_THREADS} the number of threads of the next session that
#' starts. \code{ctJuliaSetup(project=)} points the session at a different
#' engine. See \code{\link{ctJuliaSetup}}.
#' \describe{
#'   \item{\code{ctsem.julia.restart}}{Whether a fit asking for more
#'     \code{cores} than the running Julia session has threads may restart the
#'     session to get them. \code{FALSE} by default, because Julia fixes its
#'     number of threads when it starts, and a restart loses the code it has
#'     compiled. With \code{FALSE}, such a fit runs with the threads available
#'     and says so.}
#' }
#'
#' @return Nothing. This page documents options.
#' @name ctsem-options
#' @aliases ctsem-options
NULL
