# Model syntax in lavaan, Mplus and formula (brms) styles.
#
# Three readers and one writer around ctModel(). Each reader turns its text or
# formulas into one short list of statements -- a loading, a regression, a
# (co)variance, an intercept, a statement about a named parameter -- and
# `.ctSynModel()` turns that list into the cell strings ctModel() already
# reads (R/ctPar.R). Nothing downstream knows which syntax a model came from.
# `ctModelSyntax()` goes the other way: a model's `pars` to the same kind of
# statements, printed in the style asked for.
#
# What the three styles share, so that it is decided once:
#
#  * A path that is not written is zero: loadings, cross effects, predictor
#    effects. Everything else keeps ctModel()'s default: auto effects,
#    variances and the covariances of the processes are free, intercepts of
#    the indicators are free and vary between subjects, process intercepts
#    are zero. `.ctSynCells()` is that table.
#  * Every cell has a name, ctModel()'s default one or the user's label, and a
#    statement whose left side is such a name is about that parameter. That is
#    how the styles reach what they have no words for (initial states, a
#    random effect on any parameter, predictors of any parameter).
#  * A modifier may be a number, a label, a quoted ctsem cell or a ctPar()
#    call, so a cell can say anything a ctModel() cell can.

#' Write a model in lavaan, Mplus or formula syntax
#'
#' \code{ctModelLavaan}, \code{ctModelMplus} and \code{ctModelFormula} create
#' the same model object as \code{\link{ctModel}}, from lavaan style syntax,
#' Mplus style syntax, or R formulas in the style of brms. \code{ctModelSyntax}
#' writes an existing model in any of the three styles.
#'
#' @param syntax Character string containing the model.
#' @param ... For \code{ctModelFormula}, the formulas, followed by any named
#'   arguments for \code{\link{ctModel}}. For \code{ctModelLavaan} and
#'   \code{ctModelMplus}, named arguments for \code{\link{ctModel}}. Useful
#'   ones are \code{ncategories}, \code{censormin}, \code{censormax},
#'   \code{TIpredNames}, \code{tipredDefault}, \code{PARS} and
#'   \code{driftScaled}. \code{latentNames}, \code{manifestNames} and
#'   \code{TDpredNames} set the order of the variables. The model matrices
#'   themselves are written by the syntax.
#' @param type 'ct' for a continuous time model, 'dt' for discrete time.
#' @param id Name of the subject id column, or several names for nested
#'   levels, innermost first, as in \code{\link{ctModel}}.
#' @param time Name of the time column.
#' @param manifesttype Types of the indicators that are not Gaussian, as a
#'   named vector such as \code{c(y3 = 'binary', y4 = 'ordinal')}. Types are
#'   'continuous', 'binary', 'ordinal', 'count' and 'censored', or the numbers
#'   0 to 4 that \code{\link{ctModel}} uses. An unnamed vector with one entry
#'   per indicator is also accepted. Ordinal indicators need
#'   \code{ncategories}, and censored ones \code{censormin} or
#'   \code{censormax}, each as a named vector.
#' @param family For \code{ctModelFormula}, a brms style alternative to
#'   \code{manifesttype}. Either one family for every indicator or a named
#'   list, given as a name, a function or a family object.
#'   \code{gaussian}, \code{bernoulli}, \code{cumulative} and \code{poisson}
#'   are recognised.
#' @param model A model from \code{\link{ctModel}} or one of the functions
#'   here, or a fit.
#' @param style 'lavaan', 'mplus' or 'formula'.
#'
#' @details
#' \strong{What is estimated by default.} Loadings, effects between processes
#' and predictor effects are in the model only if they are written. The rest
#' follows \code{\link{ctModel}}:
#' \itemize{
#'   \item the first indicator of each process has a loading of 1, and other
#'     loadings that are written are free;
#'   \item each process has a free auto effect, free noise variance and free
#'     initial mean and variance, and the noise and initial states of
#'     different processes are correlated;
#'   \item each indicator has a free residual variance and a free intercept,
#'     and the intercept varies between subjects;
#'   \item process intercepts are zero, and residuals of different indicators
#'     are uncorrelated.
#' }
#' Fix or free any of these with the statements below.
#'
#' \strong{lavaan style.}
#' \preformatted{
#' eta1 =~ y1 + y2 + y3        # eta1 is measured by y1, y2, y3
#' eta1 ~ eta1 + a*eta2 + x    # effects on eta1: auto, cross (labelled a), predictor x
#' eta1 ~ 1                    # free process intercept
#' eta1 ~~ 0*eta2              # no noise covariance
#' y1 ~~ y2                    # residual covariance
#' y1 ~ 0*1                    # intercept fixed to zero
#' y4 | t1 + t2                # y4 is ordinal with three categories
#' a ~ age                     # age predicts the parameter a
#' a ~~ a                      # a varies between subjects
#' }
#' A number before \code{*} fixes a parameter, a name labels it, and \code{NA}
#' frees it. Two parameters with the same label are equal, as are labels
#' joined by \code{==}. A variable on the right of \code{~} that is not a
#' process is a time dependent predictor. \code{level: 2} or
#' \code{level: between} starts a block in which an indicator's name stands
#' for its intercept and a process's name for its intercept, so that
#' \code{y1 ~ age} and \code{y1 ~~ 0*y1} there are about the intercept of
#' \code{y1}. With nested ids, name the level, as in \code{level: school}.
#'
#' \strong{Mplus style.}
#' \preformatted{
#' eta1 BY y1 y2 y3;
#' eta1 ON eta1 eta2 (a);      ! label a
#' eta1 ON x;
#' [eta1];                     ! free process intercept
#' eta1 WITH eta2@0;
#' y1 WITH y2;
#' [y1@0];
#' [y4$1 y4$2];                ! y4 is ordinal with three categories
#' s | eta2 ON eta1;           ! s varies between subjects
#' \%BETWEEN\%
#' s ON age;
#' y1@0;                       ! the intercept of y1 does not vary
#' }
#' \code{@} fixes a parameter, \code{*} frees it, and a name in parentheses
#' labels the parameters before it on the same line. A lag written as
#' \code{eta1&1} is read as \code{eta1}. If the text has sections, only
#' \code{MODEL:} is read. With nested ids use \code{\%BETWEEN school\%}.
#'
#' \strong{Formula style.}
#' \preformatted{
#' ctModelFormula(
#'   y1 + y2 + y3 ~ eta1,                 # measurement
#'   d(eta1) ~ eta1 + eta2 + x,           # dynamics of eta1
#'   d(eta2) ~ 1 + eta2 + (1 | id),       # with an intercept that varies
#'   drift_eta1_eta2 ~ age + (1 | id),    # a formula for one parameter
#'   mm_y1 ~ s(age) + (0 | id))
#' }
#' A formula with \code{d()} on the left gives the effects on a process. Its
#' intercept is zero unless \code{1} is written. Other formulas with a
#' variable on the left give the processes that the variable measures. A
#' formula with a parameter name on the left describes that parameter:
#' \code{1} estimates it, a number fixes it, predictors on the right are time
#' independent predictors, \code{s()} is a smooth effect and \code{(1 | id)}
#' makes it vary. \code{(0 | id)} stops it varying. Random effect terms may
#' also be written in the other formulas, where \code{(1 + eta2 | id)} lets
#' the intercept and the effect of \code{eta2} vary. Fix or label a term with
#' \code{0.5*eta2}, \code{"a"*eta2} or \code{NA*eta1}.
#'
#' \strong{Parameter names.} Every parameter has a name, either your label or
#' the default shown by \code{print(model)} and \code{ctModelSyntax}:
#' \code{drift_eta1} and \code{drift_eta1_eta2} (effect of eta2 on eta1),
#' \code{diff_eta1} and \code{diff_eta2_eta1}, \code{cint_eta1},
#' \code{T0m_eta1}, \code{T0var_eta1}, \code{lambda_y2_eta1}, \code{mm_y1},
#' \code{mvary1}, \code{mvary2_y1} and \code{td_eta1_x}. In all three styles a
#' statement with a parameter name on the left is about that parameter. For
#' example \code{T0m_eta1 ~ 0*1} (lavaan) or \code{[T0m_eta1@0];} (Mplus)
#' fixes an initial mean, and \code{drift_eta1 ~~ drift_eta1} or
#' \code{drift_eta1;} makes an auto effect vary between subjects. A value
#' given this way is on the scale ctsem stores, which for variances is a
#' standard deviation. Values given with \code{~~}, \code{WITH} or a variance
#' statement are variances.
#'
#' \strong{Other ctsem features.} Wherever a label can be written, a quoted
#' ctsem cell or a \code{\link{ctPar}} call can be written instead, such as
#' \code{ctPar("a", transform = "exp(param)", sdscale = 0.5)*eta2} in lavaan
#' and formula style, or \code{eta1 ON eta2 (a|exp(param));} in Mplus style.
#' This gives transforms, prior scales, stated population means and state
#' dependent expressions. The result is an ordinary model, so
#' \code{model$pars}, \code{\link{ctModelMatrices}},
#' \code{\link{ctModelHigherOrder}} and the other model functions apply to it.
#' Some things belong elsewhere. Start values are given to \code{\link{ctFit}}
#' as \code{inits}, and priors with \code{priors = TRUE}. A stationary initial
#' state is \code{stationary = TRUE} in \code{\link{ctFit}}. Groups are
#' compared with a time independent predictor. Random effects are
#' correlated with each other, and \code{poprank} and \code{popblocks} in
#' \code{\link{ctFit}} give their covariance a factor structure. Higher order
#' lags are added with \code{\link{ctModelHigherOrder}}.
#'
#' @return \code{ctModelLavaan}, \code{ctModelMplus} and \code{ctModelFormula}
#'   return a model of class \code{ctStanModel}, as \code{\link{ctModel}} does.
#'   \code{ctModelSyntax} returns the code that creates the model, as a
#'   character string that prints as code.
#'
#' @examples
#' # The same model four ways
#' m1 <- ctModel(LAMBDA = diag(2), manifestNames = c('y1', 'y2'),
#'   latentNames = c('eta1', 'eta2'), silent = TRUE)
#'
#' m2 <- ctModelLavaan('
#'   eta1 =~ y1
#'   eta2 =~ y2
#'   eta1 ~ eta1 + eta2
#'   eta2 ~ eta1 + eta2')
#'
#' m3 <- ctModelMplus('
#'   eta1 BY y1;
#'   eta2 BY y2;
#'   eta1 eta2 ON eta1 eta2;')
#'
#' m4 <- ctModelFormula(
#'   y1 ~ eta1,
#'   y2 ~ eta2,
#'   d(eta1) ~ eta1 + eta2,
#'   d(eta2) ~ eta1 + eta2)
#'
#' # A cross effect that varies between subjects and depends on age
#' m5 <- ctModelLavaan('
#'   eta1 =~ y1
#'   eta2 =~ y2
#'   eta1 ~ eta1 + a*eta2
#'   eta2 ~ eta2
#'   a ~ age
#'   a ~~ a')
#'
#' ctModelSyntax(m5, 'mplus')
#' ctModelSyntax(m5, 'formula')
#'
#' @seealso \code{\link{ctModel}}, \code{\link{ctPar}}, \code{\link{ctFit}}
#' @export
ctModelSyntax <- function(model, style = c('lavaan', 'mplus', 'formula')) {
  style <- match.arg(style)
  if (!inherits(model, 'ctStanModel')) {
    base <- tryCatch(.ctFitBaseModel(model), error = function(e) NULL)
    if (!inherits(base, 'ctStanModel')) stop('ctModelSyntax() needs a model ',
      'from ctModel() or a fit.', call. = FALSE)
    model <- base
  }
  w <- .ctSynDescribe(model)
  code <- .ctSynWrite(w, style)
  # A transform is written only where it is not the cell's default, and the
  # defaults are whatever the model built from the code has.
  built <- tryCatch(suppressMessages(eval(parse(text = code),
    envir = environment(ctModelSyntax))), error = function(e) NULL)
  differs <- w$cells$kind %in% 'free'
  if (!is.null(built)) {
    at <- match(paste(w$cells$matrix, w$cells$row, w$cells$col),
      paste(built$pars$matrix, built$pars$row, built$pars$col))
    differs <- differs & !(as.character(w$cells$transform) ==
        as.character(built$pars$transform[at])) %in% TRUE
  }
  if (any(differs)) {
    w$cells$tfexplicit <- differs
    code <- .ctSynWrite(w, style)
  }
  structure(code, class = 'ctModelSyntax')
}

#' @export
print.ctModelSyntax <- function(x, ...) {
  cat(unclass(x), '\n', sep = '')
  invisible(x)
}

#' @rdname ctModelSyntax
#' @export
ctModelLavaan <- function(syntax, type = 'ct', id = 'id', time = 'time',
  manifesttype = NULL, ...) {
  .ctSynModel(.ctSynLavaan(syntax), type = type, id = id, time = time,
    manifesttype = manifesttype, args = list(...))
}

#' @rdname ctModelSyntax
#' @export
ctModelMplus <- function(syntax, type = 'ct', id = 'id', time = 'time',
  manifesttype = NULL, ...) {
  .ctSynModel(.ctSynMplus(syntax), type = type, id = id, time = time,
    manifesttype = manifesttype, args = list(...))
}

#' @rdname ctModelSyntax
#' @export
ctModelFormula <- function(..., type = 'ct', id = 'id', time = 'time',
  manifesttype = NULL, family = NULL) {
  dots <- list(...)
  # an unnamed formula is part of the model, a named one an argument
  isformula <- vapply(dots, function(x) inherits(x, 'formula') ||
      (is.list(x) && length(x) > 0L &&
          all(vapply(x, inherits, NA, what = 'formula'))), NA) &
    !nzchar(if (is.null(names(dots))) rep('', length(dots)) else names(dots))
  formulas <- do.call(c, lapply(dots[isformula], function(x)
    if (inherits(x, 'formula')) list(x) else unname(x)))
  args <- dots[!isformula]
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args))))) {
    stop('ctModelFormula() takes formulas and named arguments for ctModel().',
      call. = FALSE)
  }
  .ctSynModel(.ctSynFormula(formulas), type = type, id = id, time = time,
    manifesttype = manifesttype, args = args, family = family)
}

# --- reading text -----------------------------------------------------------

# Which characters of a string sit outside brackets and quotes. Operators and
# separators count only there, so a ctPar() call or a quoted cell can hold
# anything.
.ctSynTop <- function(x) {
  ch <- strsplit(x, '', fixed = TRUE)[[1L]]
  top <- logical(length(ch))
  depth <- 0L
  quote <- ''
  for (i in seq_along(ch)) {
    c <- ch[i]
    if (nzchar(quote)) {
      if (c == quote) quote <- ''
      next
    }
    if (c %in% c('"', "'")) quote <- c
    else if (c %in% c('(', '[', '{')) depth <- depth + 1L
    else if (c %in% c(')', ']', '}')) depth <- depth - 1L
    else top[i] <- depth == 0L
  }
  list(ch = ch, top = top)
}

# Split at a single character separator outside brackets and quotes.
.ctSynSplit <- function(x, sep) {
  s <- .ctSynTop(x)
  at <- which(s$top & s$ch == sep)
  if (sep == '+' && length(at)) {
    # the sign of an exponent, 1e+05, is no separator
    exponent <- vapply(at, function(i) i > 2L && i < length(s$ch) &&
        s$ch[i - 1L] %in% c('e', 'E') && grepl('[0-9.]', s$ch[i - 2L]) &&
        grepl('[0-9]', s$ch[i + 1L]), NA)
    at <- at[!exponent]
  }
  starts <- c(1L, at + 1L)
  ends <- c(at - 1L, length(s$ch))
  vapply(seq_along(starts), function(i) if (starts[i] > ends[i]) '' else
    paste(s$ch[starts[i]:ends[i]], collapse = ''), '')
}

# The first operator outside brackets and quotes, the longest at its position.
.ctSynFindOp <- function(x, ops) {
  s <- .ctSynTop(x)
  ops <- ops[order(-nchar(ops))]
  n <- length(s$ch)
  for (i in which(s$top)) {
    for (op in ops) {
      k <- nchar(op)
      if (i + k - 1L <= n && all(s$top[i:(i + k - 1L)]) &&
          identical(paste(s$ch[i:(i + k - 1L)], collapse = ''), op)) {
        return(list(pos = i, op = op))
      }
    }
  }
  NULL
}

# Text up to a comment character that is not inside a quote.
.ctSynUncomment <- function(line, mark) {
  ch <- strsplit(line, '', fixed = TRUE)[[1L]]
  quote <- ''
  for (i in seq_along(ch)) {
    if (nzchar(quote)) {
      if (ch[i] == quote) quote <- ''
    } else if (ch[i] %in% c('"', "'")) quote <- ch[i]
    else if (ch[i] == mark) return(substr(line, 1L, i - 1L))
  }
  line
}

.ctSynIsName <- function(x) grepl('^[A-Za-z.][A-Za-z0-9._]*$', x)

# A ctPar() call written in the syntax. Bare names are read as text, so
# ctPar(a, tipreds = c(age, sex)) needs no quotes.
.ctSynEvalCtPar <- function(e) {
  asText <- function(a) {
    if (is.name(a) && !as.character(a) %in% c('T', 'F')) return(as.character(a))
    if (is.call(a) && identical(a[[1L]], as.name('c'))) {
      return(vapply(as.list(a)[-1L], function(b)
        if (is.name(b)) as.character(b) else as.character(eval(b, baseenv())),
        ''))
    }
    a
  }
  for (i in seq_along(e)[-1L]) e[[i]] <- asText(e[[i]])
  eval(e, list(ctPar = ctPar), baseenv())
}

# What stands before a term: a number fixes, NA frees, a name labels, and a
# quoted cell or a ctPar() call gives a ctsem cell with any of its fields.
.ctSynMod <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(NULL)
  if (identical(x, 'NA')) return(list(type = 'free'))
  value <- suppressWarnings(as.numeric(x))
  if (!is.na(value)) return(list(type = 'fixed', value = value))
  if (grepl('^".*"$', x) || grepl("^'.*'$", x)) {
    return(.ctSynModCell(substr(x, 2L, nchar(x) - 1L)))
  }
  call <- regmatches(x, regexec(
    '^([A-Za-z.][A-Za-z0-9._]*)[[:space:]]*\\((.*)\\)$', x))[[1L]]
  if (length(call)) {
    fn <- call[2L]
    arg <- gsub('^["\']|["\']$', '', trimws(call[3L]))
    if (fn == 'ctPar') {
      e <- tryCatch(parse(text = x)[[1L]], error = function(err) NULL)
      if (is.null(e)) stop('Could not read ', x, call. = FALSE)
      return(.ctSynModCell(.ctSynEvalCtPar(e)))
    }
    if (fn %in% c('label', 'equal')) return(list(type = 'label', value = arg))
    if (fn == 'start') return(list(type = 'start'))
    if (fn == 'c') stop(x, ': a modifier per group is not available. ',
      'Compare groups with a time independent predictor of the parameter.',
      call. = FALSE)
    if (fn %in% c('prior', 'lower', 'upper')) stop(x, ' is not available. ',
      'Bounds are set by a transform, ctPar(name, transform = ), and priors ',
      'by ctFit(priors = TRUE).', call. = FALSE)
    stop('Modifier ', x, ' is not recognised.', call. = FALSE)
  }
  if (.ctSynIsName(x)) return(list(type = 'label', value = x))
  stop('Modifier ', x, ' is not recognised. A modifier is a number, NA, a ',
    'label, a quoted ctsem cell or a ctPar() call.', call. = FALSE)
}

.ctSynModCell <- function(text) {
  list(type = 'cell', fields = .ctCellSpecFields(as.character(text)))
}

# lavaan style ---------------------------------------------------------------

.ctSynLavaan <- function(x) {
  if (!is.character(x) || !length(x)) stop('syntax must be a character string.',
    call. = FALSE)
  lines <- unlist(strsplit(gsub('\r', '', paste(x, collapse = '\n')), '\n',
    fixed = TRUE))
  lines <- vapply(lines, .ctSynUncomment, '', mark = '#', USE.NAMES = FALSE)
  lines <- trimws(unlist(lapply(lines, .ctSynSplit, sep = ';')))
  lines <- lines[nzchar(lines)]
  ops <- c('=~', '<~', '~*~', '~~', ':=', '==', '~', '|', '<', '>')
  header <- '^(level|group|block)[[:space:]]*:[[:space:]]*([^=].*)$'
  # A line without an operator continues the one before.
  stmts <- character()
  for (l in lines) {
    if (grepl(header, l, ignore.case = TRUE) || !is.null(.ctSynFindOp(l, ops)) ||
        !length(stmts) || grepl(header, stmts[length(stmts)], ignore.case = TRUE)) {
      stmts <- c(stmts, l)
    } else stmts[length(stmts)] <- paste(stmts[length(stmts)], l)
  }
  st <- list()
  add <- function(...) st[[length(st) + 1L]] <<- list(...)
  level <- NULL
  for (s in stmts) {
    if (grepl(header, s, ignore.case = TRUE)) {
      what <- tolower(sub(header, '\\1', s, ignore.case = TRUE))
      if (what != 'level') stop("'", s, "': groups are not part of the model ",
        'syntax. Compare groups with a time independent predictor of the ',
        'parameters that may differ.', call. = FALSE)
      level <- trimws(sub(header, '\\2', s, ignore.case = TRUE))
      next
    }
    op <- .ctSynFindOp(s, ops)
    if (is.null(op)) stop("Could not read '", s, "'.", call. = FALSE)
    lhs <- trimws(.ctSynSplit(trimws(substr(s, 1L, op$pos - 1L)), '+'))
    rhs <- trimws(substr(s, op$pos + nchar(op$op), nchar(s)))
    if (op$op %in% c(':=', '<', '>', '<~', '~*~')) stop("'", s, "': ", switch(
      op$op,
      ':=' = paste('defined parameters are not part of a ctsem model. Compute',
        'them from the fit, for example with ctSummaryMatrices().'),
      '<' = , '>' = paste('bounds are set by a transform, for example',
        'ctPar("a", transform = "exp(param)")*x for a positive parameter.'),
      'this operator is not available.'), call. = FALSE)
    if (op$op == '==') {
      if (length(lhs) != 1L || !.ctSynIsName(lhs) || !.ctSynIsName(rhs)) stop(
        "'", s, "': only two labels can be set equal. Giving two parameters ",
        'the same label does the same.', call. = FALSE)
      add(kind = 'equal', a = lhs, b = rhs, src = s)
      next
    }
    terms <- lapply(.ctSynSplit(rhs, '+'), function(term) {
      p <- trimws(.ctSynSplit(term, '*'))
      mod <- NULL
      for (m in p[-length(p)]) {
        new <- .ctSynMod(m)
        if (!is.null(new) && !(identical(new$type, 'start') && !is.null(mod))) {
          mod <- new
        }
      }
      list(var = p[length(p)], mod = mod)
    })
    if (any(!nzchar(lhs)) || any(!nzchar(vapply(terms, `[[`, '', 'var')))) {
      stop("Could not read '", s, "'.", call. = FALSE)
    }
    if (op$op == '|') {
      add(kind = 'thresh', manifest = lhs[1L], n = length(terms) + 1L, src = s)
      next
    }
    for (l in lhs) for (t in terms) {
      smooth <- grepl('^s[[:space:]]*\\(', t$var)
      if (!smooth && !.ctSynIsName(t$var) && !(op$op == '~' && t$var == '1')) {
        stop("'", s, "': ", t$var, ' is not a variable name.', call. = FALSE)
      }
      switch(op$op,
        '=~' = add(kind = 'load', latent = l, manifest = t$var, mod = t$mod,
          src = s),
        '~' = if (t$var == '1') {
          add(kind = 'int', lhs = l, mod = t$mod, level = level, src = s)
        } else {
          add(kind = 'reg', lhs = l, rhs = if (smooth) NA_character_ else t$var,
            smooth = if (smooth) t$var, mod = t$mod, level = level, src = s)
        },
        '~~' = add(kind = 'cov', lhs = l, rhs = t$var, mod = t$mod,
          level = level, src = s))
    }
  }
  st
}

# Mplus style ----------------------------------------------------------------

# The items of an Mplus statement: names with @value, *value or * behind
# them, and labels in parentheses that apply to the items before them on the
# same line.
.ctSynMplusItems <- function(text, smooth = FALSE, src = text) {
  items <- list()
  for (line in strsplit(text, '\n', fixed = TRUE)[[1L]]) {
    ch <- strsplit(line, '', fixed = TRUE)[[1L]]
    tokens <- list()
    i <- 1L
    cur <- ''
    flush <- function() {
      if (nzchar(cur)) tokens[[length(tokens) + 1L]] <<- list(text = cur)
      cur <<- ''
    }
    while (i <= length(ch)) {
      c <- ch[i]
      if (grepl('[[:space:]]', c)) flush()
      else if (c == '(') {
        depth <- 1L
        j <- i
        quote <- ''
        while (j < length(ch) && depth > 0L) {
          j <- j + 1L
          if (nzchar(quote)) { if (ch[j] == quote) quote <- '' }
          else if (ch[j] %in% c('"', "'")) quote <- ch[j]
          else if (ch[j] == '(') depth <- depth + 1L
          else if (ch[j] == ')') depth <- depth - 1L
        }
        if (depth > 0L) stop("Unmatched '(' in '", src, "'.", call. = FALSE)
        inner <- paste(ch[seq_len(j - i - 1L) + i], collapse = '')
        if (smooth && identical(cur, 's')) {
          cur <- paste0('s(', inner, ')')
        } else {
          flush()
          tokens[[length(tokens) + 1L]] <- list(label = trimws(inner))
        }
        i <- j
      } else cur <- paste0(cur, c)
      i <- i + 1L
    }
    flush()
    pending <- integer()
    for (tok in tokens) {
      if (!is.null(tok$label)) {
        if (!length(pending)) stop("In '", src, "' the label (", tok$label,
          ') follows no parameter on its line.', call. = FALSE)
        lab <- tok$label
        raw <- grepl('^["\']', lab) || grepl('[|=,]', lab)
        if (raw) {
          mod <- .ctSynModCell(gsub('^["\']|["\']$', '', lab))
          for (k in pending) items[[k]]$mod <- mod
        } else {
          labs <- .ctSynMplusRange(strsplit(lab, '[[:space:]]+')[[1L]])
          labs <- ifelse(grepl('^[0-9]', labs), paste0('p', labs), labs)
          if (length(labs) == 1L) labs <- rep(labs, length(pending))
          if (length(labs) != length(pending)) stop("In '", src, "' the ",
            'labels (', lab, ') do not match the ', length(pending),
            ' parameters before them.', call. = FALSE)
          for (k in seq_along(pending)) {
            items[[pending[k]]]$mod <- list(type = 'label', value = labs[k])
          }
        }
        pending <- integer()
        next
      }
      for (one in .ctSynMplusRange(tok$text)) {
        items[[length(items) + 1L]] <- .ctSynMplusItem(one, src)
        pending <- c(pending, length(items))
      }
    }
  }
  items
}

# y1-y3 is y1, y2, y3, with anything after the range applied to each.
.ctSynMplusRange <- function(tokens) {
  unlist(lapply(tokens, function(tok) {
    m <- regmatches(tok, regexec(
      '^([A-Za-z_.][A-Za-z0-9_.]*?)([0-9]+)-([A-Za-z_.][A-Za-z0-9_.]*?)([0-9]+)([@*&$].*)?$',
      tok, perl = TRUE))[[1L]]
    if (!length(m) || m[2L] != m[4L]) return(tok)
    paste0(m[2L], seq(as.integer(m[3L]), as.integer(m[5L])), m[6L])
  }))
}

.ctSynMplusItem <- function(tok, src) {
  if (grepl('^s\\(', tok)) return(list(name = NA_character_, smooth = tok,
    mod = NULL))
  m <- regmatches(tok, regexec('^([^@*]+)([@*])?(.*)$', tok))[[1L]]
  name <- m[2L]
  mod <- NULL
  if (m[3L] == '@') {
    value <- suppressWarnings(as.numeric(m[4L]))
    if (is.na(value)) stop("In '", src, "', ", tok, ' needs a number after @.',
      call. = FALSE)
    mod <- list(type = 'fixed', value = value)
  } else if (m[3L] == '*') {
    mod <- list(type = if (nzchar(m[4L])) 'start' else 'free')
  }
  thresh <- NA_integer_
  if (grepl('$', name, fixed = TRUE)) {
    thresh <- suppressWarnings(as.integer(sub('^.*\\$', '', name)))
    name <- sub('\\$.*$', '', name)
  }
  if (grepl('&', name, fixed = TRUE)) {
    lag <- sub('^.*&', '', name)
    if (!identical(lag, '1')) stop("In '", src, "', ", tok, ': only a lag of ',
      '1 is part of the model syntax. Add higher order dynamics with ',
      'ctModelHigherOrder().', call. = FALSE)
    name <- sub('&.*$', '', name)
  }
  if (!.ctSynIsName(name)) stop("In '", src, "', ", tok,
    ' is not a variable name.', call. = FALSE)
  list(name = name, mod = mod, thresh = thresh)
}

.ctSynMplus <- function(x) {
  if (!is.character(x) || !length(x)) stop('syntax must be a character string.',
    call. = FALSE)
  x <- gsub('\r', '', paste(x, collapse = '\n'))
  x <- paste(vapply(strsplit(x, '\n', fixed = TRUE)[[1L]], .ctSynUncomment, '',
    mark = '!', USE.NAMES = FALSE), collapse = '\n')
  # A whole input file: only MODEL: is read.
  section <- '(?im)^[ \t]*(TITLE|DATA[ A-Z]*|VARIABLE|DEFINE|ANALYSIS|MODEL[ A-Z]*|OUTPUT|SAVEDATA|PLOT|MONTECARLO)[ \t]*:'
  found <- gregexpr(section, x, perl = TRUE)[[1L]]
  if (found[1L] > 0L) {
    heads <- toupper(trimws(sub(':$', '', trimws(regmatches(x, list(found))[[1L]]))))
    heads <- gsub('[[:space:]]+', ' ', heads)
    starts <- as.integer(found) + attr(found, 'match.length')
    ends <- c(as.integer(found)[-1L] - 1L, nchar(x))
    if (!'MODEL' %in% heads) stop('The Mplus syntax has sections but no MODEL: ',
      'section.', call. = FALSE)
    other <- setdiff(heads, 'MODEL')
    if (length(other)) message('Only the MODEL: section is read. Ignored: ',
      paste(other, collapse = ', '), '. Give the id and time columns, the ',
      'indicator types and the estimation settings as arguments to ',
      'ctModelMplus() and ctFit().')
    x <- paste(substring(x, starts[heads == 'MODEL'], ends[heads == 'MODEL']),
      collapse = ';')
  }
  x <- gsub('(%[^%\n]*%)', ';\\1;', x)
  stmts <- trimws(strsplit(x, ';', fixed = TRUE)[[1L]])
  stmts <- stmts[nzchar(stmts)]
  st <- list()
  add <- function(...) st[[length(st) + 1L]] <<- list(...)
  level <- NULL
  names_of <- function(items, s) {
    if (any(vapply(items, function(i) !is.null(i$mod) || is.na(i$name), NA))) {
      stop("'", s, "': the variables before the keyword take no modifiers.",
        call. = FALSE)
    }
    vapply(items, `[[`, '', 'name')
  }
  for (s in stmts) {
    flat <- gsub('[[:space:]]+', ' ', s)
    if (grepl('^%.*%$', s)) {
      block <- trimws(gsub('%', '', flat, fixed = TRUE))
      if (toupper(block) %in% c('WITHIN', 'OVERALL')) level <- NULL
      else if (toupper(block) == 'BETWEEN') level <- 'between'
      else if (grepl('^BETWEEN ', block, ignore.case = TRUE)) {
        level <- trimws(substring(block, 9L))
      } else stop(s, ' is not a block of a ctsem model. Blocks are %WITHIN%, ',
        '%BETWEEN% and %BETWEEN name%.', call. = FALSE)
      next
    }
    slope <- NULL
    bar <- .ctSynFindOp(s, '|')
    if (!is.null(bar)) {
      slope <- trimws(substr(s, 1L, bar$pos - 1L))
      s2 <- trimws(substring(s, bar$pos + 1L))
      if (!.ctSynIsName(slope)) stop("Could not read '", flat, "'.", call. = FALSE)
    } else s2 <- s
    if (grepl('^\\[', s2)) {
      inner <- sub('\\]', ' ', sub('^\\[', ' ', s2))
      items <- .ctSynMplusItems(inner, src = flat)
      steps <- list()
      for (it in items) {
        if (!is.na(it$thresh)) {
          steps[[it$name]] <- max(steps[[it$name]], it$thresh)
        } else add(kind = 'int', lhs = it$name, mod = it$mod, level = level,
          src = flat)
      }
      for (v in names(steps)) add(kind = 'thresh', manifest = v,
        n = steps[[v]] + 1L, src = flat)
      next
    }
    key <- regmatches(s2, regexec(
      '(?is)^(.*?)\\s+(BY|ON|WITH|PON|PWITH|XWITH)\\s+(.*)$', s2,
      perl = TRUE))[[1L]]
    if (!length(key)) {
      if (!is.null(slope)) stop("Could not read '", flat, "'.", call. = FALSE)
      for (it in .ctSynMplusItems(s2, src = flat)) add(kind = 'cov',
        lhs = it$name, rhs = it$name, mod = it$mod, level = level, src = flat)
      next
    }
    word <- toupper(key[3L])
    if (!word %in% c('BY', 'ON', 'WITH')) stop("'", flat, "': ", word,
      ' is not available. ', if (word == 'XWITH') paste('Interactions of',
        'processes are written as state dependent cells, see ?ctModel.') else
        'Write the effects one by one with ON or WITH.', call. = FALSE)
    left <- names_of(.ctSynMplusItems(key[2L], src = flat), flat)
    right <- .ctSynMplusItems(key[4L], smooth = word == 'ON', src = flat)
    if (!is.null(slope)) {
      if (word != 'ON' || length(left) != 1L || length(right) != 1L) stop("'",
        flat, "': a random slope is written as  s | y ON x;", call. = FALSE)
      add(kind = 'reg', lhs = left, rhs = right[[1L]]$name,
        mod = list(type = 'label', value = slope), level = level, src = flat)
      add(kind = 'vary', target = slope, level = NULL, on = TRUE, src = flat)
      next
    }
    if (word == 'BY') {
      if (length(left) != 1L) stop("'", flat, "': one process before BY.",
        call. = FALSE)
      for (it in right) add(kind = 'load', latent = left, manifest = it$name,
        mod = it$mod, src = flat)
    } else for (l in left) for (it in right) {
      if (word == 'ON') add(kind = 'reg', lhs = l, rhs = it$name,
        smooth = it$smooth, mod = it$mod, level = level, src = flat)
      else add(kind = 'cov', lhs = l, rhs = it$name, mod = it$mod,
        level = level, src = flat)
    }
  }
  st
}

# formula style --------------------------------------------------------------

.ctSynFormulaTerms <- function(rhs, src) {
  terms <- list()
  bars <- list()
  fn <- function(e) if (is.call(e) && is.name(e[[1L]])) as.character(e[[1L]]) else ''
  intercept <- function(mod) terms[[length(terms) + 1L]] <<-
    list(type = 'int', mod = mod)
  bar <- function(e) {
    inside <- e[[2L]]
    if (fn(inside) %in% c('|', '||')) inside <- inside[[2L]]
    if (fn(e) == '||') message('In ', src, ' random effects are read as ',
      'correlated, as all random effects in ctsem are.')
    vars <- character()
    int <- NA
    walk <- function(a, neg = FALSE) {
      if (fn(a) %in% c('+', '-') && length(a) == 3L) {
        walk(a[[2L]], neg)
        walk(a[[3L]], if (fn(a) == '-') !neg else neg)
      } else if (is.numeric(a)) int <<- (a != 0) && !neg
      else if (is.name(a)) vars <<- c(vars, as.character(a))
      else stop('In ', src, ' a random effect term lists 1 and variable names.',
        call. = FALSE)
    }
    walk(inside)
    levels <- function(g) {
      if (is.name(g)) return(as.character(g))
      if (fn(g) == '/') return(c(levels(g[[3L]]), levels(g[[2L]])))
      stop('In ', src, ' the grouping of a random effect term is an id name.',
        call. = FALSE)
    }
    bars[[length(bars) + 1L]] <<- list(vars = vars, int = int,
      levels = levels(e[[3L]]))
  }
  walk <- function(e, neg = FALSE) {
    f <- fn(e)
    if (f %in% c('+', '-') && length(e) == 3L) {
      walk(e[[2L]], neg)
      walk(e[[3L]], if (f == '-') !neg else neg)
    } else if (f == '+' && length(e) == 2L) walk(e[[2L]], neg)
    else if (f == '-' && length(e) == 2L && is.numeric(e[[2L]])) {
      intercept(list(type = 'fixed', value = if (neg) e[[2L]] else -e[[2L]]))
    } else if (f == '(') {
      if (fn(e[[2L]]) %in% c('|', '||')) bar(e[[2L]]) else walk(e[[2L]], neg)
    } else if (f == '*') {
      if (is.name(e[[2L]])) stop('In ', src, ', ', deparse(e), ': write a ',
        'label in quotes, "', as.character(e[[2L]]), '"*',
        deparse(e[[3L]]), '. Interactions are not part of the syntax. A time ',
        'independent predictor moderates an effect through a formula for ',
        'that parameter.', call. = FALSE)
      mod <- .ctSynMod(paste(deparse(e[[2L]], width.cutoff = 500L), collapse = ''))
      if (neg) {
        if (!identical(mod$type, 'fixed')) stop('In ', src, ' only a number ',
          'can follow a minus sign.', call. = FALSE)
        mod$value <- -mod$value
      }
      if (is.numeric(e[[3L]]) && e[[3L]] == 1) intercept(mod)
      else if (is.name(e[[3L]])) terms[[length(terms) + 1L]] <<-
        list(type = 'var', var = as.character(e[[3L]]), mod = mod)
      else stop('In ', src, ', ', deparse(e), ' is not a term.', call. = FALSE)
    } else if (f == 's') {
      terms[[length(terms) + 1L]] <<- list(type = 'smooth',
        text = paste(deparse(e, width.cutoff = 500L), collapse = ''))
    } else if (f == ':') {
      stop('In ', src, ', ', deparse(e), ': interactions are not part of the ',
        'syntax. A time independent predictor moderates an effect through a ',
        'formula for that parameter, such as drift_eta1_eta2 ~ age.',
        call. = FALSE)
    } else if (is.numeric(e)) {
      # 1 estimates the intercept, 0 or "- 1" removes it, a number fixes it
      if (e == 1) intercept(if (neg) list(type = 'fixed', value = 0) else
        list(type = 'free'))
      else intercept(list(type = 'fixed', value = if (neg) -e else e))
    } else if (is.name(e)) {
      if (neg) stop('In ', src, ' a variable cannot be subtracted.', call. = FALSE)
      terms[[length(terms) + 1L]] <<- list(type = 'var', var = as.character(e),
        mod = NULL)
    } else stop('In ', src, ', ', paste(deparse(e), collapse = ''),
      ' is not a term of the model syntax.', call. = FALSE)
  }
  walk(rhs)
  list(terms = terms, bars = bars)
}

.ctSynFormula <- function(formulas) {
  if (!length(formulas)) stop('ctModelFormula() needs at least one formula.',
    call. = FALSE)
  st <- list()
  add <- function(...) st[[length(st) + 1L]] <<- list(...)
  for (f in formulas) {
    src <- paste(deparse(f, width.cutoff = 500L), collapse = ' ')
    if (length(f) != 3L) stop(src, ' has no left side.', call. = FALSE)
    lhs <- f[[2L]]
    rhs <- .ctSynFormulaTerms(f[[3L]], src)
    kinds <- vapply(rhs$terms, `[[`, '', 'type')
    if (is.call(lhs) && identical(lhs[[1L]], as.name('d'))) {
      if (length(lhs) != 2L || !is.name(lhs[[2L]])) stop(src, ': d() takes ',
        'the name of one process.', call. = FALSE)
      p <- as.character(lhs[[2L]])
      add(kind = 'proc', name = p, src = src)
      if (any(kinds == 'smooth')) stop(src, ': a smooth is an effect on a ',
        'parameter. Write it in a formula for that parameter, such as cint_',
        p, ' ~ s(age).', call. = FALSE)
      for (t in rhs$terms) {
        if (t$type == 'int') add(kind = 'int', lhs = p, mod = t$mod,
          level = NULL, src = src)
        else add(kind = 'reg', lhs = p, rhs = t$var, mod = t$mod, level = NULL,
          dyn = TRUE, src = src)
      }
      for (b in rhs$bars) for (lev in b$levels) {
        if (isTRUE(b$int)) add(kind = 'vary', target = list(ref = 'int',
          name = p), level = lev, on = TRUE, src = src)
        for (v in b$vars) add(kind = 'vary', target = list(ref = 'effect',
          lhs = p, rhs = v), level = lev, on = TRUE, src = src)
      }
      next
    }
    lhsnames <- function(e) {
      if (is.name(e)) return(as.character(e))
      if (is.call(e) && is.name(e[[1L]]) &&
          as.character(e[[1L]]) %in% c('+', 'mvbind', 'cbind')) {
        return(unlist(lapply(as.list(e)[-1L], lhsnames)))
      }
      if (is.call(e) && identical(e[[1L]], as.name('|'))) stop(src, ': give ',
        'the type of an indicator with manifesttype or family, and censoring ',
        'limits with censormin and censormax.', call. = FALSE)
      stop(src, ': the left side is a variable, d(process) or a parameter ',
        'name.', call. = FALSE)
    }
    add(kind = 'plain', lhs = lhsnames(lhs), terms = rhs$terms, bars = rhs$bars,
      src = src)
  }
  st
}

# --- from statements to a model ---------------------------------------------

# Every cell the syntax can write, with ctModel()'s name for it and what it is
# when nothing is said about it.
.ctSynCells <- function(latents, manifests, tds) {
  nl <- length(latents)
  nm <- length(manifests)
  ntd <- length(tds)
  lab <- function(m) ctLabel(matrixname = m, n.latent = nl, n.manifest = nm,
    n.TDpred = ntd, n.TIpred = 0, Tpoints = 2, manifestNames = manifests,
    latentNames = latents, TDpredNames = tds, TIpredNames = character())
  mats <- list(
    LAMBDA = matrix(paste0('lambda_', rep(manifests, times = nl), '_',
      rep(latents, each = nm)), nm, nl),
    DRIFT = lab('DRIFT'), DIFFUSION = lab('DIFFUSION'),
    MANIFESTVAR = lab('MANIFESTVAR'), MANIFESTMEANS = lab('MANIFESTMEANS'),
    CINT = lab('CINT'), T0MEANS = lab('T0MEANS'), T0VAR = lab('T0VAR'))
  if (ntd) mats$TDPREDEFFECT <- lab('TDPREDEFFECT')
  cells <- do.call(rbind, lapply(names(mats), function(m) {
    x <- mats[[m]]
    data.frame(matrix = m, row = as.vector(row(x)), col = as.vector(col(x)),
      auto = as.character(as.vector(x)), stringsAsFactors = FALSE)
  }))
  upper <- cells$matrix %in% c('DIFFUSION', 'MANIFESTVAR', 'T0VAR') &
    cells$col > cells$row
  cells$auto[upper] <- NA_character_
  diagonal <- cells$row == cells$col
  cells$free <- !upper & (
    (cells$matrix %in% 'DRIFT' & diagonal) |
      cells$matrix %in% c('DIFFUSION', 'T0VAR', 'MANIFESTMEANS', 'T0MEANS') |
      (cells$matrix %in% 'MANIFESTVAR' & diagonal))
  n <- nrow(cells)
  cells$value <- ifelse(cells$free, NA_real_, 0)
  cells$name <- ifelse(cells$free, cells$auto, NA_character_)
  cells$transform <- NA_character_
  cells$sdscale <- NA_real_
  cells$mean <- NA_real_
  cells$expr <- FALSE
  cells$raw <- NA_character_
  cells$tipset <- FALSE
  out <- as.list(cells)
  out$vary <- vector('list', n)
  out$tip <- vector('list', n)
  out
}

.ctSynLabels <- function(st) {
  one <- function(mod) {
    if (is.null(mod)) return(NULL)
    if (identical(mod$type, 'label')) return(mod$value)
    if (identical(mod$type, 'cell') && !is.na(mod$fields$param) &&
        .ctSynIsName(mod$fields$param)) return(mod$fields$param)
    NULL
  }
  unique(unlist(lapply(st, function(s) c(one(s$mod),
    if (identical(s$kind, 'plain')) unlist(lapply(s$terms, function(t) one(t$mod))),
    if (identical(s$kind, 'equal')) NULL))))
}

# Which names are processes, indicators and time dependent predictors. A
# left side that is a parameter name makes a statement about that parameter,
# and the default names depend on the variables, so this is found by
# iteration.
.ctSynClassify <- function(st, order = list()) {
  labels <- .ctSynLabels(st)
  autos <- character()
  for (iteration in 1:20) {
    handles <- c(labels, autos)
    latents <- manifests <- tds <- outcomes <- character()
    for (s in st) {
      if (s$kind == 'load') {
        latents <- c(latents, s$latent)
        manifests <- c(manifests, s$manifest)
      } else if (s$kind == 'proc') latents <- c(latents, s$name)
      else if (s$kind == 'plain' && !any(s$lhs %in% handles)) {
        manifests <- c(manifests, s$lhs)
        latents <- c(latents, unlist(lapply(s$terms, function(t)
          if (t$type == 'var') t$var)))
      } else if (s$kind == 'reg' && is.null(s$level) && is.null(s$smooth) &&
          !s$lhs %in% handles) outcomes <- c(outcomes, s$lhs)
    }
    manifests <- unique(manifests)
    latents <- unique(c(latents, setdiff(outcomes, manifests)))
    for (s in st) {
      if (s$kind == 'reg' && is.null(s$level) && is.null(s$smooth) &&
          s$lhs %in% latents && !s$rhs %in% latents) tds <- c(tds, s$rhs)
    }
    tds <- unique(tds)
    new <- .ctSynCells(latents, manifests, tds)$auto
    new <- new[!is.na(new)]
    if (identical(new, autos)) break
    autos <- new
  }
  reorder <- function(found, wanted, what) {
    if (is.null(wanted)) return(found)
    if (!setequal(found, wanted)) stop(what, ' must name the variables the ',
      'syntax uses: ', paste(found, collapse = ', '), '.', call. = FALSE)
    as.character(wanted)
  }
  latents <- reorder(latents, order$latentNames, 'latentNames')
  manifests <- reorder(manifests, order$manifestNames, 'manifestNames')
  tds <- reorder(tds, order$TDpredNames, 'TDpredNames')
  both <- intersect(latents, manifests)
  if (length(both)) stop(paste(both, collapse = ', '), ' is used both as a ',
    'process and as an indicator. Give the process another name, for example ',
    'eta_', both[1L], ' =~ ', both[1L], '.', call. = FALSE)
  clash <- intersect(labels, c(latents, manifests, tds))
  if (length(clash)) stop('The label ', clash[1L], ' is also the name of a ',
    'variable. Use another label.', call. = FALSE)
  if (!length(latents)) stop('The model has no process. Define one with a ',
    'measurement statement.', call. = FALSE)
  list(latents = latents, manifests = manifests, tds = tds, labels = labels)
}

.ctSynLevel <- function(token, ids) {
  if (is.null(token)) return(NULL)
  t <- tolower(trimws(token))
  if (t %in% c('1', 'within', 'w')) return(NULL)
  if (t %in% c('between', 'b')) return(ids[1L])
  n <- suppressWarnings(as.integer(t))
  if (!is.na(n) && n >= 2L && n - 1L <= length(ids)) return(ids[n - 1L])
  if (trimws(token) %in% ids) return(trimws(token))
  stop("'", token, "' is not a level of this model. Its levels are: ",
    paste(ids, collapse = ', '), ' (the id argument).', call. = FALSE)
}

# Set a cell from a modifier. `default` is what a term without one means.
.ctSynApply <- function(cells, rows, mod, default = 'free', variance = FALSE,
  src = '') {
  if (is.null(mod)) mod <- list(type = default)
  for (k in rows) {
    if (is.na(cells$auto[k])) stop("'", src, "' is about a cell above the ",
      'diagonal of a covariance matrix.', call. = FALSE)
    free <- function(name) {
      cells$free[k] <<- TRUE
      cells$name[k] <<- name
      cells$value[k] <<- NA_real_
      cells$expr[k] <<- FALSE
    }
    # A fixed cell keeps a label it was given: it still names the cell, and
    # names the parameter again if the cell is made to vary around the value.
    fix <- function(value) {
      cells$free[k] <<- FALSE
      cells$value[k] <<- value
      cells$expr[k] <<- FALSE
    }
    if (mod$type %in% c('free', 'start')) {
      if (!cells$free[k] || cells$expr[k]) free(cells$auto[k])
    } else if (mod$type == 'one') fix(1)
    else if (mod$type == 'fixed') {
      value <- mod$value
      if (variance) {
        if (cells$row[k] != cells$col[k] && value != 0) stop("'", src, "': a ",
          'covariance can be fixed to zero. For other values give the whole ',
          'matrix to ctModel() with ctCov().', call. = FALSE)
        if (value < 0) stop("'", src, "': a variance cannot be negative.",
          call. = FALSE)
        value <- sqrt(value)
      }
      fix(value)
    } else if (mod$type == 'label') free(mod$value)
    else if (mod$type == 'cell') {
      f <- mod$fields
      p <- f$param
      if (is.na(p) || !nzchar(p)) stop("'", src, "': the cell has no name.",
        call. = FALSE)
      if (!is.na(suppressWarnings(as.numeric(p)))) fix(as.numeric(p))
      else if (.ctSynIsName(p) && !p %in% cells$vars) free(p)
      else {
        free(cells$auto[k])
        cells$expr[k] <- TRUE
        cells$raw[k] <- p
      }
      if (!is.na(f$transform)) cells$transform[k] <- f$transform
      if (!is.na(f$sdscale)) cells$sdscale[k] <- f$sdscale
      if (!is.na(f$mean)) cells$mean[k] <- f$mean
      iv <- f$indvarying
      if (is.logical(iv) && length(iv) == 1L && !is.na(iv)) {
        cells <- .ctSynVary(cells, k, NULL, iv)
      } else if (is.character(iv) && length(iv)) {
        for (lev in iv) cells <- .ctSynVary(cells, k, lev, TRUE)
      }
      if (!is.null(f$tipreds)) {
        cells$tipset[k] <- TRUE
        for (entry in f$tipreds) {
          parts <- trimws(strsplit(entry, '=', fixed = TRUE)[[1L]])
          cells$tip[[k]][parts[1L]] <- if (length(parts) > 1L) parts[2L] else 'TRUE'
        }
      }
    }
  }
  cells
}

.ctSynVary <- function(cells, rows, level, on) {
  add <- data.frame(level = if (is.null(level)) NA_character_ else level,
    on = on, stringsAsFactors = FALSE)
  for (k in rows) cells$vary[[k]] <- rbind(cells$vary[[k]], add)
  cells
}

# The indvarying field: empty when nothing was said, otherwise exactly the
# levels asked for.
.ctSynLevelField <- function(v, ids) {
  if (is.null(v)) return('')
  on <- unique(unlist(lapply(seq_len(nrow(v)), function(r) {
    if (!v$on[r]) NULL else if (is.na(v$level[r])) ids else v$level[r]
  })))
  if (!length(on)) 'FALSE' else if (all(ids %in% on)) 'TRUE' else
    paste(ids[ids %in% on], collapse = ',')
}

.ctSynNumber <- function(x) format(x, digits = 15)

.ctSynCellString <- function(cells, k, ids) {
  lev <- .ctSynLevelField(cells$vary[[k]], ids)
  varying <- !lev %in% c('', 'FALSE')
  tip <- cells$tip[[k]]
  where <- paste0(cells$matrix[k], '[', cells$row[k], ',', cells$col[k], ']')
  if (cells$expr[k]) {
    if (varying || length(tip) || !is.na(cells$transform[k])) stop(where,
      ' is an expression (', cells$raw[k], '), which takes no transform, ',
      'random effect or predictor effect. Give those to the parameters it ',
      'uses, in PARS.', call. = FALSE)
    return(cells$raw[k])
  }
  if (!cells$free[k]) {
    if (!varying) {
      if (length(tip)) stop(cells$auto[k], ' is fixed, so it cannot have a ',
        'predictor effect. Free it first.', call. = FALSE)
      return(.ctSynNumber(cells$value[k]))
    }
    # fixed and varying: random effects around a stated mean
    name <- paste0(if (is.na(cells$name[k])) cells$auto[k] else cells$name[k],
      '=', .ctSynNumber(cells$value[k]))
  } else {
    name <- cells$name[k]
    if (!is.na(cells$mean[k])) name <- paste0(name, '=',
      .ctSynNumber(cells$mean[k]))
  }
  blank <- function(x) if (is.na(x)) '' else as.character(x)
  parts <- c(name, blank(cells$transform[k]), lev,
    if (is.na(cells$sdscale[k])) '' else .ctSynNumber(cells$sdscale[k]))
  if (cells$tipset[k]) parts <- c(parts, paste(ifelse(tip == 'TRUE', names(tip),
    paste0(names(tip), '=', tip)), collapse = ','))
  keep <- if (cells$tipset[k]) 5L else if (grepl('=', name, fixed = TRUE)) 2L else 1L
  while (length(parts) > keep && !nzchar(parts[length(parts)])) {
    parts <- parts[-length(parts)]
  }
  paste(parts, collapse = '|')
}

.ctSynTypeCode <- function(x) {
  if (is.function(x)) x <- x()
  if (is.list(x) && !is.null(x$family)) x <- x$family
  if (is.numeric(x) && length(x) == 1L && x %in% 0:4) return(as.integer(x))
  word <- tolower(as.character(x)[1L])
  types <- list(`0` = c('continuous', 'gaussian', 'normal', '0'),
    `1` = c('binary', 'bernoulli', 'binomial', '1'),
    `2` = c('ordinal', 'ordered', 'cumulative', '2'),
    `3` = c('count', 'poisson', '3'), `4` = c('censored', '4'))
  hit <- names(types)[vapply(types, function(t) word %in% t, NA)]
  if (!length(hit)) stop("'", word, "' is not an indicator type. Types are ",
    'continuous, binary, ordinal, count and censored.', call. = FALSE)
  as.integer(hit)
}

# Per-indicator settings given by name, as ctModel()'s positional vectors.
.ctSynByManifest <- function(x, manifests, what, fill) {
  if (is.null(x)) return(NULL)
  out <- rep(fill, length(manifests))
  if (is.null(names(x))) {
    if (length(x) == 1L) x <- rep(x, length(manifests))
    if (length(x) != length(manifests)) stop(what, ' needs names, or one ',
      'entry for each of ', paste(manifests, collapse = ', '), '.',
      call. = FALSE)
    names(x) <- manifests
  }
  unknown <- setdiff(names(x), manifests)
  if (length(unknown)) stop(what, ' names ', paste(unknown, collapse = ', '),
    ', which the syntax does not use as an indicator.', call. = FALSE)
  for (n in names(x)) out[match(n, manifests)] <- x[[n]]
  out
}

.ctSynModel <- function(st, type, id, time, manifesttype, args,
  family = NULL) {
  if (!length(st)) stop('The model syntax has no statements.', call. = FALSE)
  if (!type %in% c('ct', 'dt')) stop("type must be 'ct' or 'dt'.", call. = FALSE)
  written <- c('LAMBDA', 'DRIFT', 'DIFFUSION', 'T0VAR', 'T0MEANS',
    'MANIFESTMEANS', 'MANIFESTVAR', 'CINT', 'TDPREDEFFECT', 'n.latent',
    'n.manifest', 'n.TDpred', 'n.TIpred')
  if (length(intersect(names(args), written))) stop(
    paste(intersect(names(args), written), collapse = ', '), ' is written by ',
    'the syntax. To set a matrix directly, use ctModel(), or edit ',
    'model$pars afterwards.', call. = FALSE)
  if (any(vapply(st, function(s) identical(s$mod$type, 'start'), NA))) message(
    'Start values in the syntax are not used. Give starting values to ',
    'ctFit(inits = ).')
  nm <- .ctSynClassify(st, args[c('latentNames', 'manifestNames', 'TDpredNames')])
  args[c('latentNames', 'manifestNames', 'TDpredNames')] <- NULL
  latents <- nm$latents
  manifests <- nm$manifests
  tds <- nm$tds
  ids <- as.character(id)
  cells <- .ctSynCells(latents, manifests, tds)
  cells$vars <- c(latents, manifests, tds)
  idx <- function(matrix, row, col) which(cells$matrix == matrix &
      cells$row == row & cells$col == col)
  covcell <- function(matrix, i, j) idx(matrix, max(i, j), min(i, j))
  ishandle <- function(x) x %in% c(nm$labels, cells$auto)
  # The cells a name stands for: a label, a default name, or in a between
  # block an indicator or process for its intercept.
  target <- function(x, between, src) {
    if (is.list(x)) {
      rows <- switch(x$ref,
        int = if (x$name %in% manifests) idx('MANIFESTMEANS',
          match(x$name, manifests), 1L) else idx('CINT', match(x$name, latents), 1L),
        effect = if (x$rhs %in% latents) idx('DRIFT', match(x$lhs, latents),
          match(x$rhs, latents)) else if (x$rhs %in% tds) idx('TDPREDEFFECT',
            match(x$lhs, latents), match(x$rhs, tds)) else integer(),
        load = idx('LAMBDA', match(x$manifest, manifests),
          match(x$latent, latents)))
      if (!length(rows)) stop("'", src, "' names an effect that is not in ",
        'the model.', call. = FALSE)
      return(rows)
    }
    rows <- which(!cells$expr & cells$name %in% x)
    if (!length(rows)) rows <- which(cells$auto %in% x)
    if (!length(rows) && between) {
      if (x %in% manifests) rows <- idx('MANIFESTMEANS', match(x, manifests), 1L)
      else if (x %in% latents) rows <- idx('CINT', match(x, latents), 1L)
    }
    if (!length(rows)) stop("In '", src, "', ", x, ' is not a parameter of ',
      'this model. Parameters are named by their label or their default ',
      'name, such as ', paste(utils::head(cells$auto[cells$free], 3L),
        collapse = ', '), '.', call. = FALSE)
    rows
  }
  role <- function(x) if (x %in% manifests) 'manifest' else if (x %in% latents)
    'process' else if (ishandle(x)) 'handle' else 'other'

  # Statements of formula style whose left side is a variable or parameter
  # name, as the statements the other styles write.
  expanded <- list()
  put <- function(...) expanded[[length(expanded) + 1L]] <<- list(...)
  for (s in st) {
    if (s$kind != 'plain') {
      expanded[[length(expanded) + 1L]] <- s
      next
    }
    kinds <- vapply(s$terms, `[[`, '', 'type')
    if (!any(ishandle(s$lhs))) {
      if ('smooth' %in% kinds) stop(s$src, ': a smooth is an effect on a ',
        'parameter. Write it in a formula for that parameter, such as mm_',
        s$lhs[1L], ' ~ s(age).', call. = FALSE)
      for (y in s$lhs) {
        for (t in s$terms) {
          if (t$type == 'var') put(kind = 'load', latent = t$var, manifest = y,
            mod = t$mod, src = s$src)
          else if (!identical(t$mod$type, 'free')) put(kind = 'int', lhs = y,
            mod = t$mod, level = NULL, src = s$src)
        }
        # a random effect term includes the intercept unless it says 0
        intlevels <- unlist(lapply(s$bars, function(b) if (!isFALSE(b$int)) b$levels))
        for (b in s$bars) for (lev in b$levels) {
          for (v in b$vars) put(kind = 'vary', target = list(ref = 'load',
            latent = v, manifest = y), level = lev, on = TRUE, src = s$src)
          put(kind = 'vary', target = list(ref = 'int', name = y), level = lev,
            on = lev %in% intlevels, src = s$src)
        }
      }
    } else {
      if (length(s$lhs) != 1L) stop(s$src, ': a formula for a parameter has ',
        'one name on the left.', call. = FALSE)
      ints <- s$terms[kinds == 'int']
      put(kind = 'int', lhs = s$lhs, level = NULL, src = s$src,
        mod = if (length(ints)) ints[[1L]]$mod else list(type = 'free'))
      for (t in s$terms[kinds != 'int']) put(kind = 'reg', lhs = s$lhs,
        rhs = if (t$type == 'var') t$var else NA_character_,
        smooth = if (t$type == 'smooth') t$text, mod = t$mod, level = NULL,
        src = s$src)
      for (b in s$bars) {
        if (length(b$vars)) stop(s$src, ': a parameter varies with (1 | id).',
          call. = FALSE)
        for (lev in b$levels) put(kind = 'vary', target = s$lhs, level = lev,
          on = !isFALSE(b$int), src = s$src)
      }
    }
  }

  # First the statements about variables, which also give the labels, then
  # those about parameters, so their order in the text does not matter.
  later <- list()
  measured <- character()
  steps <- list()
  for (s in expanded) {
    lev <- if (s$kind %in% c('reg', 'int', 'cov', 'vary')) .ctSynLevel(s$level, ids)
    between <- !is.null(lev)
    if (s$kind == 'load') {
      first <- !s$latent %in% measured
      measured <- c(measured, s$latent)
      cells <- .ctSynApply(cells, idx('LAMBDA', match(s$manifest, manifests),
        match(s$latent, latents)), s$mod, if (first) 'one' else 'free',
        src = s$src)
    } else if (s$kind == 'thresh') {
      if (!s$manifest %in% manifests) stop("'", s$src, "': ", s$manifest,
        ' is not an indicator.', call. = FALSE)
      steps[[s$manifest]] <- s$n
    } else if (s$kind == 'reg' && !between && is.null(s$smooth) &&
        role(s$lhs) == 'process') {
      if (s$rhs %in% manifests) stop("'", s$src, "': ", s$rhs, ' is an ',
        'indicator. Effects are between processes.', call. = FALSE)
      cells <- .ctSynApply(cells, target(list(ref = 'effect', lhs = s$lhs,
        rhs = s$rhs), FALSE, s$src), s$mod, src = s$src)
    } else if (s$kind == 'reg' && !between && role(s$lhs) == 'manifest') {
      stop("'", s$src, "': ", s$lhs, ' is an indicator, and effects are ',
        'between processes. A predictor of its intercept is written for the ',
        'parameter mm_', s$lhs, '.', call. = FALSE)
    } else if (s$kind == 'int' && role(s$lhs) %in% c('manifest', 'process')) {
      cells <- .ctSynApply(cells, target(list(ref = 'int', name = s$lhs), FALSE,
        s$src), s$mod, src = s$src)
    } else if (s$kind == 'cov' && !between &&
        all(c(role(s$lhs), role(s$rhs)) == 'manifest')) {
      cells <- .ctSynApply(cells, covcell('MANIFESTVAR', match(s$lhs, manifests),
        match(s$rhs, manifests)), s$mod, variance = TRUE, src = s$src)
    } else if (s$kind == 'cov' && !between &&
        all(c(role(s$lhs), role(s$rhs)) == 'process')) {
      cells <- .ctSynApply(cells, covcell('DIFFUSION', match(s$lhs, latents),
        match(s$rhs, latents)), s$mod, variance = TRUE, src = s$src)
    } else if (s$kind == 'cov' && !between &&
        !all(c(role(s$lhs), role(s$rhs)) == 'handle')) {
      stop("'", s$src, "': a covariance is between two processes or between ",
        'two indicators.', call. = FALSE)
    } else if (s$kind != 'proc') later[[length(later) + 1L]] <- s
  }
  tis <- character()
  smooths <- character()
  for (s in later) {
    lev <- if (s$kind %in% c('reg', 'int', 'cov', 'vary')) .ctSynLevel(s$level, ids)
    between <- !is.null(lev)
    if (s$kind == 'reg') {
      rows <- target(s$lhs, between, s$src)
      if (!is.null(s$smooth)) {
        call <- tryCatch(parse(text = s$smooth)[[1L]], error = function(e) NULL)
        if (is.null(call) || length(call) < 2L || !is.name(call[[2L]])) stop(
          "'", s$src, "': a smooth is written s(covariate) or ",
          's(covariate, k = 8).', call. = FALSE)
        pred <- as.character(call[[2L]])
        smooths[pred] <- paste(deparse(call), collapse = '')
        spec <- 'TRUE'
      } else {
        pred <- s$rhs
        if (pred %in% cells$vars) stop("'", s$src, "': ", pred, ' is a ',
          'variable of the model, and a predictor of a parameter is a time ',
          'independent predictor.', call. = FALSE)
        tis <- c(tis, pred)
        # a quoted label or number is a cell with nothing but its name
        quoted <- if (identical(s$mod$type, 'cell')) s$mod$fields$param
        spec <- if (is.null(s$mod) || s$mod$type %in% c('free', 'start')) 'TRUE' else
          if (s$mod$type == 'fixed') .ctSynNumber(s$mod$value) else
            if (s$mod$type == 'label') s$mod$value else
              if (!is.null(quoted) && !is.na(quoted)) quoted else stop("'", s$src,
                "': a predictor effect takes a number or a label.", call. = FALSE)
      }
      for (k in rows) {
        cells$tip[[k]][pred] <- spec
        cells$tipset[k] <- TRUE
      }
    } else if (s$kind == 'int') {
      cells <- .ctSynApply(cells, target(s$lhs, between, s$src), s$mod,
        src = s$src)
    } else if (s$kind == 'cov') {
      a <- target(s$lhs, between, s$src)
      b <- target(s$rhs, between, s$src)
      fixed <- identical(s$mod$type, 'fixed')
      if (identical(a, b)) {
        if ((fixed && s$mod$value != 0) || identical(s$mod$type, 'label') ||
            identical(s$mod$type, 'cell')) stop("'", s$src, "': the variance ",
          'of a random effect is estimated. It can be switched off with 0. ',
          'ctPar(sdscale = ) scales its prior.', call. = FALSE)
        cells <- .ctSynVary(cells, a, lev, !fixed)
      } else if (fixed) stop("'", s$src, "': random effects are always ",
        'correlated with each other. poprank and popblocks in ctFit() give ',
        'their covariance a simpler structure.', call. = FALSE)
    } else if (s$kind == 'vary') {
      cells <- .ctSynVary(cells, target(s$target, between, s$src), lev, s$on)
    } else if (s$kind == 'equal') {
      a <- target(s$a, FALSE, s$src)
      b <- target(s$b, FALSE, s$src)
      if (!all(cells$free[c(a, b)])) stop("'", s$src, "': both parameters ",
        'must be free.', call. = FALSE)
      cells$name[b] <- cells$name[a[1L]]
    }
  }

  mats <- list()
  for (m in unique(cells$matrix)) {
    rows <- which(cells$matrix == m)
    x <- matrix('0', max(cells$row[rows]), max(cells$col[rows]))
    for (k in rows) x[cells$row[k], cells$col[k]] <- .ctSynCellString(cells, k, ids)
    mats[[m]] <- x
  }

  # indicator types
  types <- rep(0L, length(manifests))
  ncat <- .ctSynByManifest(args$ncategories, manifests, 'ncategories', NA_integer_)
  if (is.null(ncat)) ncat <- rep(NA_integer_, length(manifests))
  for (v in names(steps)) {
    types[match(v, manifests)] <- if (steps[[v]] > 2L) 2L else 1L
    if (steps[[v]] > 2L) ncat[match(v, manifests)] <- steps[[v]]
  }
  given <- list(manifesttype, family)
  for (g in given[!vapply(given, is.null, NA)]) {
    if (is.function(g) || inherits(g, c('family', 'brmsfamily')) ||
        (is.list(g) && !is.null(g$family))) g <- list(g)
    codes <- .ctSynByManifest(lapply(as.list(g), .ctSynTypeCode), manifests,
      'manifesttype', NA_integer_)
    types[!is.na(codes)] <- unlist(codes[!is.na(codes)])
  }
  args$ncategories <- if (any(types == 2L)) ncat
  args$censormin <- .ctSynByManifest(args$censormin, manifests, 'censormin', NA_real_)
  args$censormax <- .ctSynByManifest(args$censormax, manifests, 'censormax', NA_real_)

  # time independent predictors: those the syntax uses, after any given
  tinames <- unique(tis)
  given <- args$TIpredNames
  if (inherits(given, 'formula')) {
    terms <- .ctSynSplit(paste(deparse(given[[2L]], width.cutoff = 500L),
      collapse = ''), '+')
    for (term in trimws(terms)) {
      if (grepl('^s[[:space:]]*\\(', term)) {
        cov <- as.character(parse(text = term)[[1L]][[2L]])
        if (!cov %in% names(smooths)) smooths[cov] <- term
      } else tinames <- unique(c(tinames, term))
    }
  } else if (!is.null(given)) tinames <- unique(c(as.character(given), tinames))
  args$TIpredNames <- if (length(smooths)) {
    stats::as.formula(paste('~', paste(c(smooths, setdiff(tinames,
      names(smooths))), collapse = ' + ')), env = baseenv())
  } else if (length(tinames)) tinames
  if (is.null(args$tipredDefault)) args$tipredDefault <- FALSE
  # A regression outcome that nothing measures is a process of its own, which
  # is also what a mistyped name becomes.
  unmeasured <- latents[colSums(mats$LAMBDA != '0') == 0L]
  if (length(unmeasured) && !isTRUE(args$silent)) message(
    paste(unmeasured, collapse = ', '), if (length(unmeasured) > 1L)
      ' are processes' else ' is a process', ' without indicators.')
  if (is.null(args$silent)) args$silent <- TRUE

  call <- c(list(type = type, id = id, time = time, manifestNames = manifests,
    latentNames = latents, manifesttype = as.numeric(types)), mats, args)
  if (length(tds)) call$TDpredNames <- tds
  call <- call[!vapply(call, is.null, NA)]
  if (isTRUE(args$silent)) {
    withCallingHandlers(do.call(ctModel, call),
      message = function(m) invokeRestart('muffleMessage'))
  } else do.call(ctModel, call)
}

# --- from a model to statements ---------------------------------------------

# The model's cells beside what the syntax assumes of them, which is what
# decides whether a cell needs a statement.
.ctSynDescribe <- function(model) {
  ids <- c(model$subjectIDname, model$groupIDnames)
  resolved <- .ctResolveLevels(model)
  pars <- resolved$pars
  latents <- model$latentNames
  manifests <- model$manifestNames
  tds <- model$TDpredNames
  smooths <- model$TIsmooths
  cells <- as.data.frame(.ctSynCells(latents, manifests, tds)[c('matrix', 'row',
    'col', 'auto', 'free')], stringsAsFactors = FALSE)
  names(cells)[5L] <- 'deffree'
  at <- match(paste(cells$matrix, cells$row, cells$col),
    paste(pars$matrix, pars$row, pars$col))
  if (anyNA(at)) stop('The model has no ', cells$matrix[is.na(at)][1L],
    ' matrix.', call. = FALSE)
  holds <- .ctSpecHoldsParameter(pars, latents, tds)[at]
  cells$kind <- ifelse(holds, 'free', ifelse(!is.na(pars$param[at]), 'expr',
    'fixed'))
  cells$name <- pars$param[at]
  cells$value <- pars$value[at]
  cells$transform <- pars$transform[at]
  cells$sdscale <- suppressWarnings(as.numeric(pars$sdscale[at]))
  cells$tfexplicit <- FALSE
  means <- model$fixedpopmeans
  cells$mean <- if (length(means)) unname(means[cells$name]) else NA_real_
  bad <- is.na(cells$auto) & !(cells$kind == 'fixed' & cells$value %in% 0)
  if (any(bad)) stop(cells$matrix[bad][1L], ' has a cell above its diagonal, ',
    'which the syntax cannot write.', call. = FALSE)
  # where each parameter varies, against the default for its matrix
  levelcols <- c('indvarying', if (length(model$groupIDnames))
    paste0('indvarying_', model$groupIDnames))
  on <- vapply(levelcols, function(cl) pars[[cl]][at] %in% TRUE,
    logical(nrow(cells)))
  if (is.null(dim(on))) on <- matrix(on, nrow = nrow(cells))
  cells$levels <- vapply(seq_len(nrow(cells)), function(k)
    paste(ids[on[k, ]], collapse = ','), '')
  intercept <- cells$matrix %in% c('MANIFESTMEANS', 'T0MEANS', 'CINT')
  cells$levelsdefault <- ifelse(cells$kind == 'free' & intercept,
    paste(ids, collapse = ','), '')
  # predictor effects, a smooth's columns as one term
  tinames <- model$TIpredNames
  smoothcols <- unlist(lapply(smooths, `[[`, 'columns'))
  cells$tip <- vapply(seq_len(nrow(cells)), function(k) {
    if (cells$kind[k] != 'free' || !length(tinames)) return('')
    spec <- vapply(tinames, function(ti)
      as.character(pars[[paste0(ti, '_effect')]][at[k]]), '')
    active <- stats::setNames(.ctTipredEffectActive(spec), tinames)
    out <- character()
    for (ti in setdiff(tinames[active], smoothcols)) {
      out <- c(out, if (spec[ti] == 'TRUE') ti else paste0(ti, '=', spec[ti]))
    }
    for (s in smooths) if (all(active[s$columns])) {
      out <- c(out, paste0('s(', s$covariate, ', k = ', s$k, ')'))
    }
    paste(out, collapse = ';')
  }, '')
  # Matrices the syntax has no statements for go back as arguments.
  asargs <- c('PARS', 'ASYMPTOTES', 'RESPONSE')
  extra <- setdiff(unique(pars$matrix), c(unique(cells$matrix), asargs,
    'THRESHOLDS', 'TDPREDEFFECT'))
  if (length(extra)) warning(paste(extra, collapse = ', '), ' is not written ',
    'by the syntax. Pass it as an argument.', call. = FALSE)
  holdsall <- .ctSpecHoldsParameter(pars, latents, tds)
  matrixargs <- list()
  for (m in intersect(asargs, pars$matrix)) {
    rows <- which(pars$matrix %in% m)
    rows <- rows[order(pars$row[rows], pars$col[rows])]
    text <- vapply(rows, function(r) {
      if (!is.na(pars$value[r]) && is.na(pars$param[r])) return(
        .ctSynNumber(pars$value[r]))
      if (m != 'PARS' || !holdsall[r]) return(pars$param[r])
      lev <- pars[r, levelcols, drop = TRUE] %in% TRUE
      ctPar(pars$param[r], transform = pars$transform[r],
        indvarying = if (all(lev)) TRUE else if (!any(lev)) FALSE else ids[lev],
        sdscale = pars$sdscale[r])
    }, '')
    quoted <- paste0('c(', paste0('"', text, '"', collapse = ', '), ')')
    matrixargs[[m]] <- if (max(pars$col[rows]) > 1L) paste0('matrix(', quoted,
      ', ', max(pars$row[rows]), ', ', max(pars$col[rows]), ', byrow = TRUE)') else
        quoted
  }
  list(cells = cells, latents = latents, manifests = manifests, tds = tds,
    tinames = tinames, smooths = smooths, ids = ids, model = model,
    matrixargs = matrixargs)
}

# What a cell needs written before its term, as a style-free modifier.
.ctSynCellMod <- function(w, k) {
  cell <- w$cells[k, ]
  if (cell$kind == 'expr') return(list(type = 'expr', text = cell$name))
  if (cell$kind == 'fixed') return(list(type = 'fixed', value = cell$value))
  scaled <- !is.na(cell$sdscale) && cell$sdscale != 1
  if (cell$tfexplicit || scaled || !is.na(cell$mean)) return(list(type = 'cell',
    name = cell$name, transform = if (cell$tfexplicit) cell$transform else NA,
    sdscale = if (scaled) cell$sdscale else NA, mean = cell$mean))
  if (identical(cell$name, cell$auto)) list(type = 'auto') else
    list(type = 'label', value = cell$name)
}

.ctSynModCtPar <- function(mod) {
  paste0('ctPar("', mod$name, '"',
    if (!is.na(mod$transform)) paste0(', transform = "', mod$transform, '"'),
    if (!is.na(mod$sdscale)) paste0(', sdscale = ', .ctSynNumber(mod$sdscale)),
    if (!is.na(mod$mean)) paste0(', mean = ', .ctSynNumber(mod$mean)), ')')
}

# A term in each style. `default` is what the reader assumes of a bare term:
# 'auto' (free under its default name) or 'one' (a first loading).
.ctSynTerm <- function(style, var, mod, default = 'auto', variance = FALSE) {
  bare <- (mod$type == 'auto' && default == 'auto') ||
    (mod$type == 'fixed' && default == 'one' && mod$value == 1)
  value <- if (mod$type == 'fixed') {
    if (variance && mod$value != 0 && !is.na(var[2L])) NA else
      .ctSynNumber(if (variance) mod$value^2 else mod$value)
  }
  if (style == 'mplus') {
    if (bare) return(list(item = var[1L], label = ''))
    return(switch(mod$type,
      auto = list(item = paste0(var[1L], '*'), label = ''),
      fixed = if (is.na(value)) list(item = var[1L],
        label = paste0('("', .ctSynNumber(mod$value), '")')) else
        list(item = paste0(var[1L], '@', value), label = ''),
      label = list(item = var[1L], label = paste0('(', mod$value, ')')),
      cell = list(item = var[1L], label = paste0('(', ctPar(mod$name,
        transform = mod$transform, sdscale = mod$sdscale, mean = mod$mean), ')')),
      expr = list(item = var[1L], label = paste0('("', mod$text, '")'))))
  }
  if (bare) return(var[1L])
  quote <- if (style == 'formula') '"' else ''
  paste0(switch(mod$type,
    auto = 'NA',
    fixed = if (is.na(value)) paste0('"', .ctSynNumber(mod$value), '"') else value,
    label = paste0(quote, mod$value, quote),
    cell = .ctSynModCtPar(mod),
    expr = paste0('"', mod$text, '"')), '*', var[1L])
}

.ctSynArgs <- function(w) {
  m <- w$model
  word <- c('continuous', 'binary', 'ordinal', 'count', 'censored')
  named <- function(x, keep) paste0('c(', paste0(w$manifests[keep], ' = ',
    x[keep], collapse = ', '), ')')
  type <- m$manifesttype
  out <- c(type = if (isTRUE(m$continuoustime)) '"ct"' else '"dt"',
    id = if (length(w$ids) > 1L) paste0('c(', paste0('"', w$ids, '"',
      collapse = ', '), ')') else paste0('"', w$ids, '"'),
    time = paste0('"', m$timeName, '"'))
  if (any(type != 0)) out['manifesttype'] <- named(paste0('"',
    word[type + 1L], '"'), type != 0)
  if (any(type == 2)) out['ncategories'] <- named(m$ncategories, type == 2)
  if (any(type == 4)) {
    if (any(is.finite(m$censormin[type == 4]))) out['censormin'] <-
        named(m$censormin, type == 4 & is.finite(m$censormin))
    if (any(is.finite(m$censormax[type == 4]))) out['censormax'] <-
        named(m$censormax, type == 4 & is.finite(m$censormax))
  }
  if (length(w$tinames)) {
    smoothcols <- unlist(lapply(w$smooths, `[[`, 'columns'))
    plain <- setdiff(w$tinames, smoothcols)
    out['TIpredNames'] <- if (length(w$smooths)) paste('~', paste(c(
      vapply(w$smooths, function(s) paste0('s(', s$covariate, ', k = ', s$k,
        ')'), ''), plain), collapse = ' + ')) else
      paste0('c(', paste0('"', plain, '"', collapse = ', '), ')')
  }
  for (m in names(w$matrixargs)) out[m] <- w$matrixargs[[m]]
  out
}

.ctSynWrite <- function(w, style) {
  cells <- w$cells
  idx <- function(matrix, row, col) which(cells$matrix == matrix &
      cells$row == row & cells$col == col)
  mod <- function(k) .ctSynCellMod(w, k)
  zero <- function(k) cells$kind[k] == 'fixed' && cells$value[k] %in% 0
  term <- function(...) .ctSynTerm(style, ...)
  lines <- character()
  out <- function(...) lines <<- c(lines, paste0(...))
  # Mplus: a label applies to what stands before it on its line, so a
  # labelled item takes a line of its own within the statement.
  mplus <- function(head, terms, open = '', close = '') {
    if (nzchar(open)) {
      for (t in terms) out(open, t$item, close, if (nzchar(t$label)) ' ',
        t$label, ';')
      return(invisible())
    }
    parts <- character()
    cur <- character()
    for (t in terms) {
      if (nzchar(t$label)) {
        if (length(cur)) parts <- c(parts, paste(cur, collapse = ' '))
        cur <- character()
        parts <- c(parts, paste(t$item, t$label))
      } else cur <- c(cur, t$item)
    }
    if (length(cur)) parts <- c(parts, paste(cur, collapse = ' '))
    out(head, paste(parts, collapse = paste0('\n', strrep(' ', nchar(head) + 2L))),
      ';')
  }
  nl <- length(w$latents)
  nm <- length(w$manifests)
  # In a formula a bare 1 estimates the intercept, so a value of 1 is 1*1.
  fnum <- function(value) if (value == 1) '1*1' else .ctSynNumber(value)

  # measurement
  measured <- character()
  loadterms <- function(i, j) {
    first <- !w$latents[j] %in% measured
    measured <<- c(measured, w$latents[j])
    term(if (style == 'formula') w$latents[j] else w$manifests[i],
      mod(idx('LAMBDA', i, j)), if (first) 'one' else 'auto')
  }
  intmod <- function(matrix, i) mod(idx(matrix, i, 1L))
  if (style == 'formula') {
    for (i in seq_len(nm)) {
      rhs <- character()
      im <- intmod('MANIFESTMEANS', i)
      if (im$type != 'auto') rhs <- if (im$type == 'fixed')
        fnum(im$value) else term('1', im)
      for (j in seq_len(nl)) if (!zero(idx('LAMBDA', i, j))) rhs <- c(rhs,
        loadterms(i, j))
      if (!length(rhs)) rhs <- '1'
      out(w$manifests[i], ' ~ ', paste(rhs, collapse = ' + '))
    }
  } else for (j in seq_len(nl)) {
    terms <- lapply(which(vapply(seq_len(nm), function(i)
      !zero(idx('LAMBDA', i, j)), NA)), loadterms, j = j)
    if (!length(terms)) next
    if (style == 'lavaan') out(w$latents[j], ' =~ ',
      paste(unlist(terms), collapse = ' + '))
    else mplus(paste0(w$latents[j], ' BY '), terms)
  }

  # dynamics
  for (i in seq_len(nl)) {
    terms <- list()
    for (j in seq_len(nl)) {
      k <- idx('DRIFT', i, j)
      if (i == j || !zero(k)) terms[[length(terms) + 1L]] <-
        term(w$latents[j], mod(k))
    }
    for (j in seq_along(w$tds)) {
      k <- idx('TDPREDEFFECT', i, j)
      if (!zero(k)) terms[[length(terms) + 1L]] <- term(w$tds[j], mod(k))
    }
    cm <- intmod('CINT', i)
    cint <- !(cm$type == 'fixed' && cm$value == 0)
    if (style == 'formula') {
      if (cint) terms <- c(list(if (cm$type == 'auto') '1' else if (
        cm$type == 'fixed') fnum(cm$value) else term('1', cm)), terms)
      out('d(', w$latents[i], ') ~ ', paste(unlist(terms), collapse = ' + '))
    } else if (style == 'lavaan') {
      out(w$latents[i], ' ~ ', paste(unlist(terms), collapse = ' + '))
      if (cint) out(w$latents[i], ' ~ ', term('1', cm))
    } else {
      mplus(paste0(w$latents[i], ' ON '), terms)
      if (cint) mplus('', list(term(w$latents[i], cm)), '[', ']')
    }
  }

  # statements about a parameter by its name
  parline <- function(k, m = mod(k)) {
    h <- cells$auto[k]
    if (style == 'formula') out(h, ' ~ ', if (m$type == 'fixed')
      fnum(m$value) else term('1', m))
    else if (style == 'lavaan') out(h, ' ~ ', term('1', m))
    else mplus('', list(term(h, m)), '[', ']')
  }
  # (co)variances and intercepts that differ from the default
  for (k in which(cells$matrix %in% c('DIFFUSION', 'MANIFESTVAR', 'T0VAR',
    'MANIFESTMEANS', 'T0MEANS') & !is.na(cells$auto))) {
    m <- mod(k)
    default <- if (cells$deffree[k]) m$type == 'auto' else zero(k)
    if (default) next
    mat <- cells$matrix[k]
    if (mat == 'MANIFESTMEANS' && style == 'formula') next
    if (style == 'formula' || mat %in% c('T0VAR', 'T0MEANS')) {
      parline(k, m)
      next
    }
    if (mat == 'MANIFESTMEANS') {
      v <- w$manifests[cells$row[k]]
      if (style == 'lavaan') out(v, ' ~ ', term('1', m))
      else mplus('', list(term(v, m)), '[', ']')
      next
    }
    vars <- if (mat == 'DIFFUSION') w$latents else w$manifests
    a <- vars[cells$row[k]]
    b <- vars[cells$col[k]]
    # free where the default is fixed needs its marker; a bare term frees
    if (m$type == 'auto') m <- list(type = 'label', value = cells$auto[k])
    offdiag <- if (a == b) NA else TRUE
    if (m$type == 'label' && identical(m$value, cells$auto[k])) {
      if (style == 'lavaan') out(a, ' ~~ ', b)
      else if (a == b) out(a, ';') else out(a, ' WITH ', b, ';')
    } else if (style == 'lavaan') {
      out(a, ' ~~ ', term(c(b, offdiag), m, variance = TRUE))
    } else if (a == b) mplus('', list(term(c(a, offdiag), m, variance = TRUE)))
    else mplus(paste0(a, ' WITH '), list(term(c(b, offdiag), m, variance = TRUE)))
  }

  # where parameters vary and what predicts them
  free <- which(cells$kind == 'free')
  first <- free[!duplicated(cells$name[free])]
  blocks <- list()
  for (k in first) {
    h <- cells$name[k]
    tips <- if (nzchar(cells$tip[k])) strsplit(cells$tip[k], ';', fixed = TRUE)[[1L]]
    varies <- cells$levels[k] != cells$levelsdefault[k]
    if (!varies && !length(tips)) next
    on <- if (nzchar(cells$levels[k])) strsplit(cells$levels[k], ',')[[1L]]
    all <- length(on) == length(w$ids)
    tipterm <- function(t) {
      p <- strsplit(t, '=', fixed = TRUE)[[1L]]
      if (grepl('^s\\(', t)) return(list(item = t, label = ''))
      if (length(p) == 1L) return(term(p, list(type = 'auto')))
      value <- suppressWarnings(as.numeric(p[2L]))
      term(p[1L], if (is.na(value)) list(type = 'label', value = p[2L]) else
        list(type = 'fixed', value = value))
    }
    if (style == 'formula') {
      rhs <- vapply(tips, function(t) { x <- tipterm(t)
        if (is.list(x)) x$item else x }, '')
      if (varies) rhs <- c(rhs, if (length(on)) paste0('(1 | ', on, ')') else
        paste0('(0 | ', w$ids[1L], ')'))
      out(h, ' ~ ', paste(rhs, collapse = ' + '))
      next
    }
    if (length(tips)) {
      terms <- lapply(tips, tipterm)
      if (style == 'lavaan') out(h, ' ~ ', paste(vapply(terms, function(x)
        if (is.list(x)) x$item else x, ''), collapse = ' + '))
      else mplus(paste0(h, ' ON '), terms)
    }
    if (!varies) next
    if (all || !length(on)) {
      if (style == 'lavaan') out(h, ' ~~ ', if (all) '' else '0*', h)
      else out(h, if (all) '' else '@0', ';')
    } else for (lev in on) blocks[[lev]] <- c(blocks[[lev]], h)
  }
  for (lev in names(blocks)) {
    out(if (style == 'lavaan') paste0('level: ', lev) else
      paste0('%BETWEEN ', lev, '%'))
    for (h in blocks[[lev]]) out(if (style == 'lavaan') paste0(h, ' ~~ ', h) else
      paste0(h, ';'))
  }

  args <- .ctSynArgs(w)
  # the order of variables, where it is not the order they are first written
  rebuilt <- tryCatch(.ctSynClassify(if (style == 'formula') .ctSynFormula(
    lapply(lines, function(l) stats::as.formula(l, env = baseenv()))) else if (
      style == 'lavaan') .ctSynLavaan(lines) else .ctSynMplus(lines)),
    error = function(e) NULL)
  for (what in c('latents', 'manifests', 'tds')) {
    if (!identical(as.character(rebuilt[[what]]), as.character(w[[what]])) &&
        length(w[[what]])) {
      args[c(latents = 'latentNames', manifests = 'manifestNames',
        tds = 'TDpredNames')[[what]]] <- paste0('c(', paste0('"', w[[what]], '"',
          collapse = ', '), ')')
    }
  }
  args <- paste0(names(args), ' = ', args, collapse = ', ')
  if (style == 'formula') {
    return(paste0('ctModelFormula(\n', paste0('  ', lines, ',', collapse = '\n'),
      '\n  ', args, ')'))
  }
  text <- paste0('  ', lines, collapse = '\n')
  fn <- if (style == 'lavaan') 'ctModelLavaan' else 'ctModelMplus'
  if (grepl("'", text, fixed = TRUE)) return(paste0(fn, '(', args, ', syntax = ',
    paste(deparse(paste0('\n', text, '\n')), collapse = ''), ')'))
  paste0(fn, '(', args, ", syntax = '\n", text, "\n')")
}
