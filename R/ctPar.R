# Writing a parameter specification without counting pipes.
#
# ctsem's matrix cells carry more than a name. A cell can say what the parameter
# is called, how it is transformed, whether it varies over subjects, how its
# prior is scaled, and which time-independent predictors act on it -- and it
# says all of that in one string, with the fields separated by `|` and absent
# fields left empty:
#
#   'mm||TRUE||TI1'
#
# That is a parameter named `mm`, default transform, individually varying,
# default sdscale, with an effect from TI1. It is compact and it is completely
# unreadable: the meaning of `TRUE` depends on counting the pipes before it, and
# a reader has no way to know how many fields there are or what order they come
# in. Writing one means consulting the documentation every time; reading someone
# else's means counting.
#
# `ctPar()` builds the same string from named arguments. It is not a new
# specification language -- the output is exactly the string the parser already
# takes, and hand-written strings keep working -- it just moves the field names
# from the documentation into the call.

#' Write a parameter specification
#'
#' Builds the \code{|}-separated cell specification ctsem's model matrices use,
#' from named arguments rather than positional fields.
#'
#' @param name Parameter name. Two cells given the same name are the same
#'   parameter, which is how equality constraints are written.
#' @param transform Transform expression in terms of \code{param}, e.g.
#'   \code{'exp(param)'}. Omitted leaves the default for that matrix, which
#'   depends on the cell -- a variance, a drift diagonal and a mean all differ.
#' @param indvarying Whether the parameter varies over subjects.
#' @param sdscale Multiplier on the population standard deviation's prior. Only
#'   meaningful when \code{indvarying} is TRUE.
#' @param tipreds Time-independent predictors acting on this parameter: a
#'   character vector of names for effects to estimate.
#'
#' @return A single character string, suitable for a model matrix cell.
#'
#' @details The output is the same string the \code{|} syntax produces, so this
#'   is a way of writing one rather than an alternative to it. Fields left at
#'   their default are emitted empty, which is what the parser reads as
#'   "unspecified".
#'
#' @examples
#' ctPar('mm', indvarying = TRUE)
#' ctPar('drift11', transform = '-exp(param)')
#' ctPar('mm', indvarying = TRUE, tipreds = c('age', 'sex'))
#'
#' \dontrun{
#' model <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
#'   LAMBDA = matrix(1),
#'   MANIFESTMEANS = matrix(ctPar('mm', indvarying = TRUE, sdscale = 0.5)))
#' }
#'
#' @seealso \code{\link{ctModel}} for the matrices these cells go in.
#' @export
ctPar <- function(name, transform = NA, indvarying = NA, sdscale = NA,
  tipreds = NULL) {

  if (missing(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("ctPar() needs a parameter name.", call. = FALSE)
  }
  field <- function(x) {
    if (length(x) != 1L || is.na(x)) return("")
    if (is.logical(x)) return(if (x) "TRUE" else "FALSE")
    as.character(x)
  }
  parts <- c(as.character(name), field(transform), field(indvarying),
    field(sdscale))

  if (!is.null(tipreds) && length(tipreds)) {
    if (!is.character(tipreds)) {
      stop("ctPar(tipreds=) takes predictor names. An effect is estimated; ",
        "to fix one to a value, set it in the model's matrices afterwards.",
        call. = FALSE)
    }
    parts <- c(parts, paste(tipreds, collapse = ","))
  }

  # Trailing empty fields carry no information and only add pipes to read past.
  while (length(parts) > 1L && !nzchar(parts[length(parts)])) {
    parts <- parts[-length(parts)]
  }
  paste(parts, collapse = "|")
}
