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
# `ctParSpec()` builds the same string from named arguments. It is not a new
# specification language -- the output is exactly the string the parser already
# takes, and hand-written strings keep working -- it just moves the field names
# from the documentation into the call.
#
# The same named fields can also be written straight into a cell, so nothing has
# to be called at all. These three specify the same parameter:
#
#   'mm||TRUE|0.5'                            the `|` form
#   ctParSpec('mm', indvarying = TRUE, sdscale = 0.5)
#   'mm, indvarying=TRUE, sdscale=0.5'        the named form, in the cell
#
# `.ctCellSpecToPipe()` below normalises the named form to the `|` form, so
# there stays exactly one thing for the model parser to read.

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
#'   character vector of names for effects to estimate, e.g. \code{'age'} or
#'   \code{c('age', 'sex')}. Naming any predictor switches the rest off for
#'   this cell, so the list is the whole set of effects on it. A name must be
#'   one of the model's \code{TIpredNames}; anything else is an error when the
#'   model is built.
#'
#' @return A single character string, suitable for a model matrix cell.
#'
#' @details The output is the same string the \code{|} syntax produces, so this
#'   is a way of writing one rather than an alternative to it. Fields left at
#'   their default are emitted empty, which is what the parser reads as
#'   "unspecified".
#'
#'   The same named fields can also be written directly into a model matrix
#'   cell, in which case nothing needs to be called. A cell containing
#'   \code{|} is read as a \code{|}-separated specification; otherwise a cell
#'   with a field of the form \code{name=value} is read as a named
#'   specification, and anything else as a plain parameter name. All three of
#'   these specify the same parameter:
#'
#'   \preformatted{
#'   MANIFESTMEANS = matrix('mm||TRUE|0.5')
#'   MANIFESTMEANS = matrix(ctParSpec('mm', indvarying = TRUE, sdscale = 0.5))
#'   MANIFESTMEANS = matrix('mm, indvarying=TRUE, sdscale=0.5')
#'   }
#'
#'   In the named form, fields are separated by commas (commas inside brackets
#'   belong to the value, so \code{transform=pnorm(param, 0, 1)} is one field),
#'   the parameter name comes first, and several predictors are given as
#'   \code{tipreds=c(age, sex)}. A field name that is not one of the above is
#'   an error naming the nearest valid field, so a misspelling such as
#'   \code{topreds=age} is refused rather than becoming part of a parameter
#'   name. Comparisons are unaffected: a state-dependent cell may contain
#'   \code{==}, \code{>=}, \code{<=} or \code{!=}, and a named argument
#'   inside a call belongs to the value it appears in.
#'
#' @examples
#' ctParSpec('mm', indvarying = TRUE)
#' ctParSpec('drift11', transform = '-exp(param)')
#' ctParSpec('mm', indvarying = TRUE, tipreds = c('age', 'sex'))
#'
#' \dontrun{
#' # These two models are identical.
#' model <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
#'   LAMBDA = matrix(1),
#'   MANIFESTMEANS = matrix(ctParSpec('mm', indvarying = TRUE, sdscale = 0.5)))
#'
#' model <- ctModel(type = 'ct', n.latent = 1, n.manifest = 1,
#'   LAMBDA = matrix(1),
#'   MANIFESTMEANS = matrix('mm, indvarying=TRUE, sdscale=0.5'))
#' }
#'
#' @seealso \code{\link{ctModel}} for the matrices these cells go in.
#' @export
ctParSpec <- function(name, transform = NA, indvarying = NA, sdscale = NA,
  tipreds = NULL) {

  if (missing(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("ctParSpec() needs a parameter name.", call. = FALSE)
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
      stop("ctParSpec(tipreds=) takes predictor names. An effect is estimated; ",
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

# --- reading a cell, either way it was written -------------------------------
#
# The fields a cell can carry are exactly ctParSpec()'s arguments, in the order
# the `|` form puts them.
.ctCellSpecKeys <- c('name', 'transform', 'indvarying', 'sdscale', 'tipreds')

# A cell is a named specification when one of its comma-separated fields opens
# with `<name>=`. Deliberately not keyed on the field being one of the *known*
# names: a misspelling then went undetected, and `'mm, topreds=age'` became a
# free parameter literally named "mm, topreds=age" carrying whatever tipred
# effects the model defaults to. Nothing warned, and the eventual failure was
# `Could not retrieve body of '=()'` from the symbolic differentiator, which
# names neither the cell nor the typo. An unrecognised name is now an error
# from .ctCellSpecToPipe(), which can say what was meant.
#
# Two things must still keep meaning what they say, and the shape of the test
# is what protects them:
#
#  * A comparison. `state[1]>=2`, `x==0` and `a!=b` have no `<name>=` at the
#    start of a field -- the character before the `=` is an operator, and
#    `=(?!=)` rules out the `==` case besides.
#  * The ordered form, whose tipreds field may itself contain `name=value`
#    ('mm||||age=4.3,sex=2'). A `|` anywhere says the cell is that form, so it
#    is never read as a named specification.
#
# Fields are split bracket-aware, so a named argument inside a call --
# `transform=pnorm(param, mean=0)` -- is part of a value and not a field of
# its own.
.ctCellSpecIsNamed <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) return(FALSE)
  if (grepl('|', x, fixed = TRUE)) return(FALSE)
  parts <- trimws(.ctCellSpecSplit(x))
  any(grepl('^[A-Za-z.][A-Za-z0-9._]*[[:space:]]*=(?!=)', parts, perl = TRUE))
}

# Split on commas that are not inside brackets, so a value may contain commas:
# `transform=pnorm(param, 0, 1)` is one field, not three.
.ctCellSpecSplit <- function(x) {
  chars <- strsplit(x, '', fixed = TRUE)[[1]]
  depth <- 0L
  parts <- character(0)
  start <- 1L
  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (ch %in% c('(', '[', '{')) depth <- depth + 1L
    else if (ch %in% c(')', ']', '}')) depth <- depth - 1L
    else if (ch == ',' && depth <= 0L) {
      parts <- c(parts, paste(chars[start:(i - 1L)], collapse = ''))
      start <- i + 1L
    }
  }
  c(parts, if (start <= length(chars)) paste(chars[start:length(chars)],
    collapse = '') else '')
}

# Named form -> `|` form. Anything that is not a named form is returned
# untouched, so this is safe to run over every cell of every matrix.
.ctCellSpecToPipe <- function(x) {
  if (!.ctCellSpecIsNamed(x)) return(x)
  parts <- trimws(.ctCellSpecSplit(x))
  parts <- parts[nzchar(parts)]
  args <- list()
  for (i in seq_along(parts)) {
    part <- parts[i]
    keyed <- regmatches(part,
      regexpr('^[A-Za-z.][A-Za-z0-9._]*[[:space:]]*=(?!=)', part, perl = TRUE))
    if (!length(keyed)) {
      if (i != 1L) {
        stop("In the parameter specification '", x, "', the field '", part,
          "' has no name. Only the parameter name is positional, and it comes ",
          'first; write the rest as ',
          paste(paste0(.ctCellSpecKeys[-1], '='), collapse = ', '),
          call. = FALSE)
      }
      args$name <- part
      next
    }
    key <- trimws(sub('=$', '', keyed))
    value <- trimws(substring(part, nchar(keyed) + 1L))
    if (!key %in% .ctCellSpecKeys) {
      stop("In the parameter specification '", x, "', '", key, "' is not a ",
        'field. Valid fields are: ',
        paste(.ctCellSpecKeys, collapse = ', '), '.',
        .ctCellSpecDidYouMean(key), call. = FALSE)
    }
    if (key %in% names(args)) {
      stop("Field '", key, "' given twice in the parameter specification '",
        x, "'.", call. = FALSE)
    }
    args[[key]] <- switch(key,
      indvarying = .ctCellSpecLogical(value, x),
      sdscale = .ctCellSpecNumeric(value, x),
      tipreds = .ctCellSpecNames(value),
      value)
  }
  if (is.null(args$name)) {
    stop("The parameter specification '", x, "' has no parameter name. It ",
      'comes first, before any named field.', call. = FALSE)
  }
  do.call(ctParSpec, args)
}

# A misspelled field name is the common case, so say which field was probably
# meant rather than only listing all of them. `tipred` for `tipreds` and
# `indvarying` for `indvaring` are both one edit away.
.ctCellSpecDidYouMean <- function(key) {
  d <- as.integer(utils::adist(tolower(key), .ctCellSpecKeys)[1, ])
  if (min(d) > max(1L, nchar(key) %/% 3L)) return('')
  paste0(" Did you mean '", .ctCellSpecKeys[which.min(d)], "'?")
}

.ctCellSpecLogical <- function(value, x) {
  out <- switch(toupper(value), 'TRUE' = TRUE, 'T' = TRUE, 'FALSE' = FALSE,
    'F' = FALSE, NULL)
  if (is.null(out)) {
    stop("indvarying=", value, " in '", x, "' must be TRUE or FALSE.",
      call. = FALSE)
  }
  out
}

.ctCellSpecNumeric <- function(value, x) {
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out)) stop("sdscale=", value, " in '", x, "' must be a number.",
    call. = FALSE)
  out
}

# `tipreds=age`, `tipreds=c(age, sex)` and `tipreds=c('age','sex')` all mean the
# same list of names.
.ctCellSpecNames <- function(value) {
  value <- sub('^c\\(', '', sub('\\)$', '', trimws(value)))
  out <- trimws(strsplit(value, ',', fixed = TRUE)[[1]])
  out <- gsub('^["\']|["\']$', '', out)
  out[nzchar(out)]
}

# Read a cell into its fields, whichever way it was written. Fields the cell
# leaves unspecified come back NA (NULL for tipreds). Used by the
# `model$matrices <- ` / `ctModelMatrices() <- ` path, which assigns into an
# existing pars row rather than building one; ctModel() normalises to the `|`
# form instead and lets its own parser read the fields.
.ctCellSpecFields <- function(x) {
  out <- list(param = NA_character_, transform = NA_character_,
    indvarying = NA, sdscale = NA_real_, tipreds = NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x)) return(out)
  x <- .ctCellSpecToPipe(x)
  if (!grepl('|', x, fixed = TRUE)) {
    out$param <- trimws(x)
    return(out)
  }
  split <- gsub(' ', '', strsplit(x, '|', fixed = TRUE)[[1]], fixed = TRUE)
  if (length(split) > length(.ctCellSpecKeys)) {
    stop('Param spec has too many separators!  ', x, call. = FALSE)
  }
  for (i in seq_along(split)) {
    if (!nzchar(split[i])) next
    key <- .ctCellSpecKeys[i]
    if (identical(key, 'name')) out$param <- split[i]
    else if (identical(key, 'indvarying')) out$indvarying <- as.logical(split[i])
    else if (identical(key, 'sdscale')) out$sdscale <- as.numeric(split[i])
    else if (identical(key, 'tipreds')) out$tipreds <- .ctCellSpecNames(split[i])
    else out[[key]] <- split[i]
  }
  out
}
