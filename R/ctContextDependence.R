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
  mats <- ctStanMatricesList()$all
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
  matrices <- paste0(unique(cells$matrix), collapse = ", ")
  kinds <- paste0(unique(unname(.ctContextKindLabels[unique(cells$kind)])), collapse = " and ")
  paste0("Cells of ", matrices, " depend on the ", kinds,
    ", so they have no single value; the values reported here were evaluated at ",
    label, " and are conditional on it.",
    if (!is.null(remedy)) paste0(" ", remedy) else "")
}

# Emit .ctContextNote() once per call site. Silent for a linear model, and
# silent for a model whose only context dependence is carrier states.
.ctContextMessage <- function(fit, label, remedy = NULL) {
  note <- .ctContextNote(.ctFitConditionalCells(fit), label, remedy)
  if (!is.null(note)) message(note)
  invisible(note)
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
