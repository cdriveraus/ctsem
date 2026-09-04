# Context dependence of model matrix cells -----------------------------------
#
# A ctsem matrix cell may be written by an expression rather than by a bare
# parameter, and that expression can reference things the parameter vector does
# not determine. Such a cell has no single value: it has a value *at* an
# evaluation context. Every summary and plot in the package reports one anyway,
# at a context none of them names, and this file is the one place that works out
# which cells those are and what they depend on.
#
# Four kinds of dependence, which need different treatment and must not be
# lumped together:
#
#   'carrier'  `state[k]` with k > nlatent -- an `intoverpop` random effect.
#              ctsem carries an individually varying parameter as a latent state
#              with no drift and no diffusion, so this reference *is* the
#              parameter. Evaluating it at the subject's last row is not an
#              arbitrary choice but the best available estimate: the smoother has
#              seen all of that subject's data by then. Nothing here should
#              report a carrier-only cell as a reporting problem, and the tests
#              guard that explicitly -- "fixing" it would break the augmented
#              approach.
#
#   'state'    `state[k]` with k <= nlatent -- a real dynamic process. The value
#              genuinely varies along the trajectory, so a single reported number
#              is conditional on wherever it was evaluated.
#
#   'tdpred'   `tdpreds[rowi, j]` -- a time dependent predictor. Exactly the same
#              problem as 'state', and treated identically. The package currently
#              evaluates these at zero without saying so.
#
#   (time)     Not reachable from the spec language today: no rewrite rule
#              produces a `t` or `dt` reference, though the Julia engine's row
#              context carries both. Included in the evaluation context
#              abstraction anyway, since the engine already plumbs it.
#
# `ctModelStatesAndPARS()` (R/ctModelWriter.R) rewrites user-written latent and
# TD predictor names into the bracket forms above, and rewrites references to a
# PARS entry into `PARS[r,c]`. That last one makes dependence transitive: a cell
# referencing a PARS cell inherits whatever that PARS cell depends on, so the
# search below runs to a fixed point rather than reading each expression once.
#
# Both backends are served from here rather than from two implementations. The
# expressions themselves are identical -- Julia consumes the same rewritten
# strings Stan does -- so only the route to the table of cells differs.

.ctContextKindLabels <- c(
  carrier = "individually varying parameter (carrier state)",
  state = "latent state",
  tdpred = "time dependent predictor")

# Kinds that make a reported value conditional on an evaluation context. A
# carrier reference does not: see the header.
.ctContextProblemKinds <- c("state", "tdpred")

# Indices referenced as `name[<integer>...]` in an expression.
#
# `tdpreds` is written `tdpreds[rowi, j]` by ctModelStatesAndPARS and
# `ctx.tdpreds[j]` after the Julia rewrite, so the row coordinate is optional.
.ctExpressionIndices <- function(expression, name) {
  if (!length(expression)) return(integer())
  pattern <- paste0("(?:ctx\\.)?\\b", name, "\\s*\\[\\s*(?:rowi\\s*,\\s*)?([0-9]+)\\s*\\]")
  matches <- regmatches(expression, gregexpr(pattern, expression, perl = TRUE))[[1L]]
  if (!length(matches)) return(integer())
  as.integer(sub(paste0("^.*?([0-9]+)\\s*\\]$"), "\\1", matches))
}

# The (row, col) coordinates of every `PARS[r,c]` an expression references.
.ctExpressionParsRefs <- function(expression) {
  if (!length(expression)) return(NULL)
  matches <- regmatches(expression,
    gregexpr("\\bPARS\\s*\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]", expression, perl = TRUE))[[1L]]
  if (!length(matches)) return(NULL)
  coords <- regmatches(matches, regexec("([0-9]+)\\s*,\\s*([0-9]+)", matches))
  data.frame(
    row = as.integer(vapply(coords, `[`, character(1L), 2L)),
    col = as.integer(vapply(coords, `[`, character(1L), 3L)),
    stringsAsFactors = FALSE)
}

# What one expression depends on directly, ignoring PARS indirection.
.ctExpressionDependence <- function(expression, nlatent) {
  if (!length(expression) || is.na(expression) || !nzchar(expression)) return(character())
  kinds <- character()
  states <- .ctExpressionIndices(expression, "state")
  if (length(states)) {
    if (any(states <= nlatent)) kinds <- c(kinds, "state")
    if (any(states > nlatent)) kinds <- c(kinds, "carrier")
  }
  if (length(.ctExpressionIndices(expression, "tdpreds"))) kinds <- c(kinds, "tdpred")
  kinds
}

# Number of *real* latent processes, i.e. the index above which a state
# reference is a random-effect carrier rather than a dynamic process.
.ctFitNlatent <- function(fit) {
  if (inherits(fit, "ctStanModel")) return(as.integer(fit$n.latent))
  if (!is.null(fit$standata$nlatent)) return(as.integer(fit$standata$nlatent))
  spec <- .ctBackendSpec(fit)
  if (!is.null(spec$nlatent)) return(as.integer(spec$nlatent))
  model <- .ctFitModelObject(fit)
  if (!is.null(model$n.latent)) return(as.integer(model$n.latent))
  stop("Cannot determine the number of latent processes for this fit.", call. = FALSE)
}

# Every cell of every model matrix, with the expression that writes it.
#
# Backend-neutral output: data.frame(matrix, row, col, expression). The two
# routes differ only in where the rewritten expression strings are kept.
.ctContextCellTable <- function(fit) {
  # An unfitted model carries the user's own strings, with latent and TD
  # predictor *names* still in them. Run the same rewrite ctFit would, so the
  # single expression parser below sees one language rather than two.
  if (inherits(fit, "ctStanModel")) {
    pars <- ctModelStatesAndPARS(fit$pars, statenames = fit$latentNames,
      tdprednames = fit$TDpredNames)
    expression <- replace(as.character(pars$param), is.na(pars$param), "")
    # ctsem's own rewrite only touches `param`, so a reference written into a
    # `name | transform` cell is left as a bare latent name. Both backends read
    # `param` too, which makes such a cell unsupported upstream rather than
    # merely undetected -- but reporting the model as linear would be the worse
    # of the two failures, so scan the transform column as well, through the
    # same rewrite.
    if (!is.null(pars$transform)) {
      shim <- pars
      shim$param <- replace(as.character(pars$transform), is.na(pars$transform), "")
      shim <- ctModelStatesAndPARS(shim, statenames = fit$latentNames,
        tdprednames = fit$TDpredNames)
      expression <- trimws(paste(expression, shim$param))
    }
    return(data.frame(matrix = as.character(pars$matrix), row = as.integer(pars$row),
      col = as.integer(pars$col), expression = expression, stringsAsFactors = FALSE))
  }

  spec <- .ctBackendSpec(fit)
  if (!is.null(spec$parameter_table)) {
    table <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
    blank <- function(x) if (is.null(x)) rep("", nrow(table)) else replace(as.character(x), is.na(x), "")
    # `param` is blanked for dynamic cells when the parameter table is built
    # (R/ctJuliaBackend.R), the expression having moved into the three transform
    # columns -- so read those, and keep `param` too for anything not rewritten.
    expression <- trimws(paste(blank(table$predicttransform), blank(table$updatetransform),
      blank(table$tdtransform), blank(table$param)))
    return(data.frame(matrix = as.character(table$matrix), row = as.integer(table$row),
      col = as.integer(table$col), expression = expression, stringsAsFactors = FALSE))
  }

  matsetup <- fit$setup$matsetup
  if (is.null(matsetup)) matsetup <- fit$ctstanmodel$modelmats$matsetup
  if (is.null(matsetup)) stop("The fit does not carry a parameter table.", call. = FALSE)
  mats <- .ctMatricesList()$all
  names <- rep(NA_character_, max(mats))
  names[mats] <- base::names(mats)
  out <- data.frame(matrix = names[matsetup$matrix], row = as.integer(matsetup$row),
    col = as.integer(matsetup$col),
    expression = replace(as.character(matsetup$parname), is.na(matsetup$parname), ""),
    stringsAsFactors = FALSE)

  # Calcs write matrix cells from outside the parameter table: `PARS[1,1] = ...`
  # and friends. The left hand side names the cell, the right hand side carries
  # the dependence. Without these a model that does its nonlinearity in `calcs`
  # rather than in a cell expression looks linear.
  calcs <- unique(unlist(fit$ctstanmodel$modelmats$calcs))
  calcs <- calcs[!is.na(calcs) & nzchar(calcs) & grepl("=", calcs, fixed = TRUE)]
  for (calc in calcs) {
    lhs <- trimws(sub("=.*", "", calc))
    rhs <- sub("^[^=]*=", "", calc)
    target <- regexec("^([A-Za-z_][A-Za-z0-9_]*)\\s*\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]$", lhs)
    parts <- regmatches(lhs, target)[[1L]]
    if (length(parts) != 4L) next
    out <- rbind(out, data.frame(matrix = parts[2L], row = as.integer(parts[3L]),
      col = as.integer(parts[4L]), expression = rhs, stringsAsFactors = FALSE))
  }
  out[!is.na(out$matrix), , drop = FALSE]
}

#' Which model matrix cells depend on something other than the parameters
#'
#' Internal. Returns a long data frame -- one row per (cell, kind of
#' dependence) -- naming every cell whose value is written by an expression
#' referencing a latent state, a random-effect carrier state, or a time
#' dependent predictor. Zero rows for an ordinary linear model.
#'
#' Dependence is transitive through PARS: a cell referencing `PARS[r,c]`
#' inherits whatever that PARS cell depends on, so the search runs to a fixed
#' point.
#'
#' @param fit A ctStanFit or ctJuliaFit.
#' @return data.frame(matrix, row, col, kind), kind one of 'carrier', 'state',
#'   'tdpred'.
#' @noRd
.ctFitContextDependentCells <- function(fit) {
  cells <- .ctContextCellTable(fit)
  nlatent <- .ctFitNlatent(fit)
  empty <- data.frame(matrix = character(), row = integer(), col = integer(),
    kind = character(), stringsAsFactors = FALSE)
  if (!nrow(cells)) return(empty)

  direct <- lapply(cells$expression, .ctExpressionDependence, nlatent = nlatent)
  refs <- lapply(cells$expression, .ctExpressionParsRefs)

  # Propagate through PARS references until nothing changes. Bounded by the
  # number of cells, so a cyclic reference terminates rather than spinning.
  key <- paste(cells$matrix, cells$row, cells$col, sep = "\r")
  for (pass in seq_len(nrow(cells))) {
    changed <- FALSE
    for (i in seq_along(refs)) {
      if (is.null(refs[[i]])) next
      wanted <- paste("PARS", refs[[i]]$row, refs[[i]]$col, sep = "\r")
      source <- unique(unlist(direct[match(wanted, key)]))
      source <- source[!is.na(source)]
      added <- setdiff(source, direct[[i]])
      if (length(added)) {
        direct[[i]] <- c(direct[[i]], added)
        changed <- TRUE
      }
    }
    if (!changed) break
  }

  counts <- lengths(direct)
  if (!sum(counts)) return(empty)
  out <- data.frame(
    matrix = rep(cells$matrix, counts),
    row = rep(cells$row, counts),
    col = rep(cells$col, counts),
    kind = unlist(direct),
    stringsAsFactors = FALSE)
  out <- unique(out)
  out <- out[order(out$matrix, out$row, out$col, out$kind), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# The subset that makes a reported number conditional on an evaluation context.
.ctFitConditionalCells <- function(fit) {
  cells <- .ctFitContextDependentCells(fit)
  cells[cells$kind %in% .ctContextProblemKinds, , drop = FALSE]
}

# Cells referencing *both* a carrier state and something dynamic.
#
# These are the one case that cannot be reported as a subject parameter at all.
# A carrier-only cell is that subject's parameter, exactly; a dynamic-only cell
# is not a parameter and is excluded from subject parameter reporting anyway.
# A cell that is both looks like a parameter and reports a number that also
# moves with the trajectory, so its last-row value silently conflates the
# individual difference with wherever that subject happened to end up.
.ctContextConflatedCells <- function(fit) {
  cells <- .ctFitContextDependentCells(fit)
  if (!nrow(cells)) return(cells[0L, , drop = FALSE])
  key <- paste(cells$matrix, cells$row, cells$col, sep = "\r")
  carrier <- unique(key[cells$kind %in% "carrier"])
  dynamic <- unique(key[cells$kind %in% .ctContextProblemKinds])
  both <- intersect(carrier, dynamic)
  out <- unique(cells[key %in% both, c("matrix", "row", "col"), drop = FALSE])
  rownames(out) <- NULL
  out
}

#' The one sentence every context-dependence message uses
#'
#' @param cells As returned by \code{.ctFitConditionalCells}.
#' @param label Human readable name of the evaluation point, e.g. 'the T0MEANS
#'   state, with time dependent predictors at zero'.
#' @param remedy Optional trailing sentence suggesting what to do instead.
#' @return A single string, or NULL when there is nothing to say.
#' @noRd
.ctContextNote <- function(cells, label, remedy = NULL) {
  if (is.null(cells) || !nrow(cells)) return(NULL)
  # Jacobian blocks are derivatives of the model matrices, not matrices anyone
  # is shown, so naming them here would send the reader looking for output that
  # does not exist. They stay in the cells table, which is the programmatic
  # answer, and out of the sentence, which is the human one.
  reportable <- setdiff(unique(cells$matrix), .ctBackendJacobianMatrices)
  if (!length(reportable)) reportable <- unique(cells$matrix)
  matrices <- paste0(reportable, collapse = ", ")
  kinds <- paste0(unique(unname(.ctContextKindLabels[unique(cells$kind)])), collapse = " and ")
  paste0("Cells of ", matrices, " depend on the ", kinds,
    ", so they have no single value. The values reported here were evaluated at ",
    label, ", and are conditional on that.",
    if (!is.null(remedy)) paste0(" ", remedy) else "")
}

# The two evaluation points the package actually uses, named once so that every
# message spells them the same way.
.ctContextPopLabel <- "the T0MEANS state, with time dependent predictors at zero"
.ctContextSubjectLabel <- "each subject's last observed row"

# What the reader can do about it, which differs by backend: only the julia
# engine can materialise the matrices at a caller-chosen point. Retrofitting
# that into Stan's generated code is a large change to a backend that is no
# longer the development line, so its message says so rather than pretending.
.ctContextRemedy <- function(fit) {
  if (inherits(fit, "ctJuliaFit")) {
    "Pass state= (or 'mean', 'asymptotic') to evaluate them elsewhere."
  } else {
    "Evaluating them at another point requires backend='julia'."
  }
}

# Emit .ctContextNote() from a user-facing entry point. Silent for a linear
# model, and silent for a model whose only context dependence is carrier
# states. Internal callers attach the attribute without messaging, so that one
# summary() does not print the same sentence five times.
.ctContextMessage <- function(fit, label, remedy = .ctContextRemedy(fit)) {
  note <- .ctContextNote(.ctFitConditionalCells(fit), label, remedy)
  if (!is.null(note)) message(note)
  invisible(note)
}

# Attach the cells to a returned object so a caller can act on them
# programmatically rather than by parsing a message.
.ctContextAttach <- function(x, fit) {
  cells <- try(.ctFitContextDependentCells(fit), silent = TRUE)
  if (inherits(cells, "try-error")) return(x)
  attr(x, "contextDependent") <- cells
  x
}

#' Does a model or fit have context-dependent (nonlinear) matrix cells?
#'
#' A ctsem matrix cell may be specified as an expression referencing a latent
#' process or a time dependent predictor -- for instance
#' \code{DRIFT[1,1] = '-log1p(exp(param)) * eta2'}. Such a cell has no single
#' value, only a value at a given state and set of predictor values, and the
#' package's summaries and plots must therefore name the point they evaluated
#' it at.
#'
#' References to an individually varying parameter's carrier state (the
#' \code{intoverpop} representation of a random effect) do not count: those are
#' parameters, not dynamics, and are reported exactly.
#'
#' @param x A \code{ctStanModel}, \code{ctStanFit} or \code{ctJuliaFit}.
#' @return Logical.
#' @examples
#' ctModelIsNonlinear(ctstantestfit)
#' @export
ctModelIsNonlinear <- function(x) {
  cells <- try(.ctFitConditionalCells(x), silent = TRUE)
  if (inherits(cells, "try-error")) return(NA)
  nrow(cells) > 0L
}

#' Context dependence of reported model matrices
#'
#' @description
#' A ctsem matrix cell may be written as an expression referencing a latent
#' process or a time dependent predictor -- for instance
#' \code{DRIFT[1,1] = '-log1p(exp(param)) * eta2'}. Such a cell has no single
#' value. It has a value \emph{at} an evaluation context: a latent state, a set
#' of time dependent predictor values, a time and an interval.
#'
#' @details
#' ctsem's summaries and plots report a number for every cell, so for a model
#' like this they must pick a point. There are two, and which one you get
#' depends on what you asked for:
#'
#' \itemize{
#'   \item \strong{Population matrices} (\code{summary}, \code{ctSummaryMatrices},
#'     \code{ctDiscretePars(subjects='popmean')}, \code{ctTIpredEffects}) are
#'     evaluated at the population \code{T0MEANS} state, with time dependent
#'     predictors at zero.
#'   \item \strong{Subject matrices} (\code{ctSubjectPars},
#'     \code{ctDiscretePars(subjects=)}) are evaluated at each subject's last
#'     observed row.
#' }
#'
#' Only \code{\link{ctKalman}} and \code{\link{ctPredict}} avoid the choice
#' entirely, because they run the filter and re-evaluate every cell at every
#' step.
#'
#' A reference to an \emph{individually varying} parameter is not affected by
#' any of this. ctsem represents such a parameter as a carrier latent state
#' with no drift and no diffusion, so a cell reading that state simply is the
#' parameter, and the estimate at the subject's last row is the one that has
#' seen all of that subject's data -- the best available, not an arbitrary
#' point. Only references to real dynamic processes, and to time dependent
#' predictors, make a reported value conditional.
#'
#' Functions affected by this attach the cells in question as
#' \code{attr(x, 'contextDependent')}, and say so once in a message.
#' \code{\link{ctModelIsNonlinear}} answers the question directly. With
#' \code{backend='julia'} the evaluation point can be chosen; see the
#' \code{state} argument of \code{\link{ctBackendParMatrices}}.
#'
#' @name ctContextDependence
#' @seealso \code{\link{ctModelIsNonlinear}}, \code{\link{ctBackendParMatrices}}
NULL


# Choosing an evaluation point ------------------------------------------------
#
# Reporting at the T0MEANS state is a default, not a finding, and for many
# nonlinear models it is a poor one -- T0 is often nowhere near where the data
# lives. These resolve the shorthands a caller can pass instead.
#
# Everything here returns the *augmented* state vector the engine indexes, with
# carrier entries left at their population values: 'mean' and 'asymptotic' are
# statements about where the dynamic processes are, and moving a carrier would
# silently change which parameter values the matrices were built from.

.ctContextStateOptions <- c("T0MEANS", "mean", "asymptotic")

# The population T0MEANS at engine (augmented) length, used both as the default
# and as the padding for a shorter supplied state.
.ctContextBaseState <- function(fit, tipreds = NULL) {
  # suppressMessages: this is a lookup on the way to choosing an evaluation
  # point, not a report of one, and the note belongs to the caller's choice.
  as.numeric(suppressMessages(
    ctBackendParMatrices(fit, tipreds = tipreds, trim = FALSE))$T0MEANS[, 1])
}

.ctContextPadState <- function(fit, state, tipreds = NULL) {
  base <- .ctContextBaseState(fit, tipreds)
  state <- as.numeric(state)
  if (length(state) == length(base)) return(state)
  nlatent <- .ctFitNlatent(fit)
  if (length(state) != nlatent) {
    stop("state must have ", nlatent, " entries (one per latent process)",
      if (length(base) != nlatent) paste0(", or ", length(base),
        " for the augmented state the filter uses") else "",
      "; got ", length(state), ".", call. = FALSE)
  }
  # Carrier entries keep their population values: see the header.
  base[seq_len(nlatent)] <- state
  base
}

# Mean smoothed latent state over every row of every subject.
#
# Where the data actually is, as opposed to where the process started. Read
# from the filter rather than derived, so it is the same quantity ctKalman
# plots.
.ctContextMeanState <- function(fit) {
  smoothed <- suppressMessages(ctKalmanArray(fit, pointest = TRUE)$etasmooth)
  .ctContextPadState(fit, apply(smoothed, 3L, mean, na.rm = TRUE))
}

# The system's own fixed point: the state at which the deterministic change is
# zero, DRIFT(x) x + CINT(x) = 0.
#
# Newton, using the engine's own JAx as the Jacobian rather than finite
# differences -- JAx is exactly d(drift)/d(state), so one engine call per
# iteration does the whole step. For a linear model this lands on asymCINT in a
# single iteration, which is the invariant the test checks.
.ctContextAsymptoticState <- function(fit, tolerance = 1e-8, maxiter = 50L,
  tipreds = NULL) {
  nlatent <- .ctFitNlatent(fit)
  x <- .ctContextBaseState(fit, tipreds)
  discrete <- !isTRUE(.ctFitModelObject(fit)$continuoustime)
  for (iteration in seq_len(maxiter)) {
    mats <- suppressMessages(ctBackendParMatrices(fit, state = x,
      tipreds = tipreds, trim = FALSE))
    drift <- mats$DRIFT[seq_len(nlatent), seq_len(nlatent), drop = FALSE]
    cint <- as.numeric(mats$CINT[seq_len(nlatent), 1])
    jacobian <- mats$JAx[seq_len(nlatent), seq_len(nlatent), drop = FALSE]
    if (all(!is.finite(jacobian)) || all(jacobian == 0)) jacobian <- drift
    # Discrete time asks where x = A x + c, i.e. (A - I) x + c = 0.
    if (discrete) {
      drift <- drift - diag(nlatent)
      jacobian <- jacobian - diag(nlatent)
    }
    residual <- as.numeric(drift %*% x[seq_len(nlatent)]) + cint
    if (max(abs(residual)) < tolerance) return(x)
    step <- try(solve(jacobian, residual), silent = TRUE)
    if (inherits(step, "try-error")) {
      stop("No asymptotic state: the system is singular at the current iterate. ",
        "Use state='mean' or supply a state.", call. = FALSE)
    }
    x[seq_len(nlatent)] <- x[seq_len(nlatent)] - step
  }
  stop("No asymptotic state found in ", maxiter, " iterations -- the system may ",
    "have no stable fixed point. Use state='mean' or supply a state.", call. = FALSE)
}

#' Resolve a state argument to an evaluation point
#'
#' @return list(state = NULL or numeric at engine length, label = character).
#'   A NULL state means "the engine's own default", which is T0MEANS.
#' @noRd
.ctResolveState <- function(fit, state = NULL, tipreds = NULL) {
  if (is.null(state)) return(list(state = NULL, label = .ctContextPopLabel))
  if (is.character(state)) {
    state <- match.arg(state, .ctContextStateOptions)
    if (identical(state, "T0MEANS")) {
      return(list(state = NULL, label = .ctContextPopLabel))
    }
    if (identical(state, "mean")) {
      return(list(state = .ctContextMeanState(fit),
        label = "the mean smoothed latent state"))
    }
    return(list(state = .ctContextAsymptoticState(fit, tipreds = tipreds),
      label = "the system's asymptotic (fixed point) state"))
  }
  list(state = .ctContextPadState(fit, state, tipreds), label = "the supplied state")
}
