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
      tdprednames = fit$TDpredNames, tiprednames = fit$TIpredNames)
    expression <- replace(as.character(pars$param), is.na(pars$param), "")
    # ctsem's own rewrite only touches `param`, so a reference written into a
    # `name | transform` cell is left as a bare latent name. Both backends read
    # `param` too, which makes such a cell unsupported upstream rather than
    # merely undetected -- but reporting the model as linear would be the worse
    # of the two failures, so scan the transform column as well, through the
    # same rewrite. The rewrite takes PARS names from the PARS rows' `param`,
    # so the transform rows go in under no matrix and the PARS cells are added
    # as written: with the transforms in the PARS rows, a PARS cell's
    # transform `param` was read as a PARS name, and every transform came back
    # referencing the first PARS cell.
    if (!is.null(pars$transform)) {
      shim <- pars
      shim$param <- replace(as.character(pars$transform), is.na(pars$transform), "")
      shim$matrix <- rep("", nrow(shim))
      named <- fit$pars[as.character(fit$pars$matrix) %in% "PARS", , drop = FALSE]
      named$matrix <- as.character(named$matrix)
      shim <- ctModelStatesAndPARS(rbind(shim, named), statenames = fit$latentNames,
        tdprednames = fit$TDpredNames, tiprednames = fit$TIpredNames)
      expression <- trimws(paste(expression, shim$param[seq_len(nrow(pars))]))
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

# The random effects the augmented route linearises through dynamics that
# depend on the state, by name; ctFit warns of them.
#
# The augmented route makes each random effect a state of the filter. Where a
# DRIFT, CINT, DIFFUSION or TDPREDEFFECT cell reads a latent state, the
# extended Kalman filter linearises that cell at the predicted mean, and an
# effect whose spread reaches the states it reads is linearised with it, over
# the effect's whole population spread until the data narrow it: a random
# initial mean puts that spread on the states before the first observation.
# On a milk-yield model whose yield drift read fat (ten cows,
# scripts/nlgap-2026-11) the augmented route's maximum read -2843 by its filter
# and about -3467 by a particle filter; the Laplace route's -2879.
#
# An effect is linearised so where such a cell reads it, directly or through
# PARS, or where it moves the mean of a state -- a T0MEANS, DRIFT, CINT or
# TDPREDEFFECT cell in that state's row -- from which the drift carries it to
# a state such a cell reads (or the column a state dependent DRIFT cell
# multiplies). A measurement effect (MANIFESTMEANS, LAMBDA, MANIFESTVAR) moves
# no state: it reaches the states through the update alone, as the data do.
# Beside a random initial mean it keeps the first observation from pinning
# that mean down, and that combination is where the milk model's gap was:
# with the initial means alone random their sds went to zero (filter and
# particle filter 10.8 apart, 10.6 without random effects), and with the
# manifest means alone the augmented maximum was 4 above Laplace's
# (dev1, 2026-11-16). A DIFFUSION effect moves no mean. Effects on the
# measurement's state dependence are not looked at here.
.ctAugmentedNonlinearEffects <- function(model) {
  pars <- model$pars
  varying <- which(.ctVaryingRows(model, "indvarying") & !is.na(pars$param))
  if (!length(varying)) return(character())
  nlatent <- as.integer(model$n.latent)
  cells <- .ctContextCellTable(model)
  key <- paste(cells$matrix, cells$row, cells$col)
  # The rows each cell reads through PARS, itself included, closed over PARS
  # cells reading others; and the latent states those rows read.
  refs <- lapply(cells$expression, function(e) {
    r <- .ctExpressionParsRefs(e)
    if (is.null(r)) integer() else
      unique(stats::na.omit(match(paste("PARS", r$row, r$col), key)))
  })
  deps <- lapply(seq_len(nrow(cells)), function(i) {
    out <- i
    repeat {
      more <- setdiff(unlist(refs[out]), out)
      if (!length(more)) return(out)
      out <- c(out, more)
    }
  })
  reads <- lapply(deps, function(d) {
    s <- unlist(lapply(cells$expression[d], .ctExpressionIndices, name = "state"))
    unique(s[s >= 1L & s <= nlatent])
  })
  statedependent <- cells$matrix %in% c("DRIFT", "CINT", "DIFFUSION",
    "TDPREDEFFECT") & lengths(reads) > 0L
  if (!any(statedependent)) return(character())
  # The states the linearisation is taken over, and every state the drift
  # carries into them: state k is carried into by state j through a nonzero
  # DRIFT[k, j] or a mean cell of row k reading state j.
  feeds <- unique(c(unlist(reads[statedependent]),
    cells$col[statedependent & cells$matrix %in% "DRIFT"]))
  meancell <- cells$matrix %in% c("DRIFT", "CINT", "TDPREDEFFECT")
  nonzero <- !is.na(pars$param) |
    (!is.na(pars$value) & suppressWarnings(as.numeric(pars$value)) != 0)
  upstream <- lapply(seq_len(nlatent), function(k) {
    rows <- which(meancell & cells$row == k)
    unique(c(cells$col[rows[cells$matrix[rows] %in% "DRIFT" & nonzero[rows]]],
      unlist(reads[rows])))
  })
  repeat {
    more <- setdiff(unlist(upstream[feeds]), feeds)
    if (!length(more)) break
    feeds <- c(feeds, more)
  }
  moves <- cells$matrix %in% c("T0MEANS", "DRIFT", "CINT", "TDPREDEFFECT")
  reached <- vapply(varying, function(e) {
    readers <- vapply(deps, function(d) e %in% d, logical(1))
    any(readers & statedependent) || any(cells$row[readers & moves] %in% feeds)
  }, logical(1))
  unique(pars$param[varying[reached]])
}

#' The one sentence every context-dependence message uses
#'
#' @param cells As returned by \code{.ctFitConditionalCells}.
#' @param label Human readable name of the evaluation point, e.g. 'the T0MEANS
#'   state, with time dependent predictors at zero'.
#' @param remedy Optional trailing sentence suggesting what to do instead.
#' @return A single string, or NULL when there is nothing to say.
#' @noRd
# The matrix names a sentence or a figure caption may use.
#
# Jacobian blocks are derivatives of the model matrices, not matrices anyone is
# shown, so naming them sends the reader looking for output that does not exist.
# They stay in the cells table, which is the programmatic answer, and out of the
# prose, which is the human one. Shared rather than inlined, because the second
# copy of this rule was written without it and put "DRIFT/JAx" in a plot
# subtitle.
.ctContextReportableMatrices <- function(cells) {
  if (is.null(cells) || !nrow(cells)) return(character())
  reportable <- setdiff(unique(cells$matrix), .ctBackendJacobianMatrices)
  if (!length(reportable)) reportable <- unique(cells$matrix)
  reportable
}

.ctContextNote <- function(cells, label, remedy = NULL) {
  if (is.null(cells) || !nrow(cells)) return(NULL)
  matrices <- paste0(.ctContextReportableMatrices(cells), collapse = ", ")
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

# The matrices the standardiser is built from.
#
# Standardising divides by the stationary standard deviations implied by
# asymDIFFUSIONcov, which solves the (continuous or discrete) Lyapunov equation
# in DRIFT and DIFFUSION -- JAx being DRIFT's jacobian, which is what the solve
# uses once a cell depends on the state.
.ctContextStandardiserMatrices <- c('DRIFT', 'JAx', 'DIFFUSION', 'DIFFUSIONcov',
  'asymDIFFUSIONcov')

# Standardising a state dependent system is a second linearisation, on top of
# the one the reported matrix already is, and it is the one a reader is least
# likely to have in mind: the divisor is not the process's spread, it is the
# spread the system would settle to if it stayed linear at this point forever.
# Two evaluation points can therefore give curves of different heights from
# identical dynamics. Silent for a linear model, and silent for a model whose
# nonlinearity is nowhere near the dynamics -- a state dependent LAMBDA does
# not move the latent stationary variance.
.ctContextStandardiseNote <- function(fit, label = NULL) {
  cells <- try(.ctFitConditionalCells(fit), silent = TRUE)
  if (inherits(cells, 'try-error') || is.null(cells) || !nrow(cells)) return(NULL)
  hit <- cells[cells$matrix %in% .ctContextStandardiserMatrices, , drop = FALSE]
  if (!nrow(hit)) return(NULL)
  paste0('standardise=TRUE divides by the stationary standard deviations from ',
    'asymDIFFUSIONcov, and cells of ',
    paste0(.ctContextReportableMatrices(hit), collapse = ', '),
    " depend on the context, so that is the linearised system's spread",
    if (!is.null(label)) paste0(' at ', label) else '',
    " rather than the nonlinear system's own. The scale moves with the ",
    "evaluation point; standardise=FALSE reports the processes' own units.")
}

# Emit it, for the entry points that standardise.
.ctContextStandardiseMessage <- function(fit, label = NULL) {
  note <- .ctContextStandardiseNote(fit, label)
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
#' A ctsem matrix cell may be an expression of a latent process or a time
#' dependent predictor, such as
#' \code{DRIFT[1,1] = '-log1p(exp(param)) * eta2'}. Such a cell has a value
#' only at a given state and set of predictor values, so summaries and plots
#' report the point at which they evaluated it.
#'
#' References to the state that carries an individually varying parameter
#' (with \code{intoverpop='augmented'}) do not count, since these are
#' parameters, reported exactly.
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
#' A ctsem matrix cell may be an expression of a latent process or a time
#' dependent predictor, such as
#' \code{DRIFT[1,1] = '-log1p(exp(param)) * eta2'}. Such a cell has a value
#' only at a given evaluation context: a latent state, a set of time dependent
#' predictor values, a time and an interval.
#'
#' @details
#' ctsem's summaries and plots report a number for every cell, so for such a
#' model they evaluate it at one of two points, depending on the function:
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
#' \code{\link{ctKalman}} and \code{\link{ctPredict}} evaluate every cell at
#' every step of the filter, so they report no single point.
#'
#' References to an \emph{individually varying} parameter are not affected,
#' since the subject's estimate of such a parameter uses all of that
#' subject's data. Only references to dynamic processes and to
#' time dependent predictors make a reported value depend on the evaluation
#' point.
#'
#' Functions affected by this attach the cells as
#' \code{attr(x, 'contextDependent')}, and give a message naming the point.
#' \code{\link{ctModelIsNonlinear}} reports whether a model has such cells.
#' With \code{backend='julia'} the evaluation point can be chosen with the
#' \code{state} argument of \code{\link{ctSummaryMatrices}},
#' \code{\link{ctDiscretePars}}, \code{\link{ctExtract}} and
#' \code{\link{ctNetwork}}: \code{'T0MEANS'}, \code{'mean'} (the mean smoothed
#' latent state), \code{'asymptotic'} (the fixed point the processes settle
#' at), or a numeric vector with one value per latent process.
#' \code{\link{ctPhasePortrait}} and \code{\link{ctStateDependencePlot}} show
#' how the dynamics change over the state.
#'
#' @examples
#' \donttest{
#' if (isTRUE(ctJuliaStatus()$available)) {
#' # The auto effect of eta1 changes with the level of eta2, so DRIFT[1, 1]
#' # has a value only at a given state of eta2.
#' model <- ctModel(LAMBDA = diag(2),
#'   manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
#'   PARS = c('a|-log1p_exp(param)', 'b'),
#'   DRIFT = matrix(c('a + b * eta2', 'drift21', 0, 'drift22'), 2, 2))
#' model$pars$indvarying <- FALSE
#' fit <- ctFit(ctstantestdat, model)
#' ctModelIsNonlinear(fit)
#'
#' # The population matrices at the default point. The message names the
#' # point, and the attribute lists the cells that depend on it.
#' m <- ctSummaryMatrices(fit)
#' m$estimate$DRIFT
#' attr(m, 'contextDependent')
#'
#' # The same matrix at the fixed point, and at a chosen state with eta1 at 0
#' # and eta2 at 2.
#' ctSummaryMatrices(fit, state = 'asymptotic')$estimate$DRIFT
#' ctSummaryMatrices(fit, state = c(0, 2))$estimate$DRIFT
#' }
#' }
#'
#' @name ctContextDependence
#' @seealso \code{\link{ctModelIsNonlinear}}, \code{\link{ctSummaryMatrices}},
#'   \code{\link{ctPhasePortrait}}, \code{\link{ctStateDependencePlot}}
NULL


# Choosing an evaluation point ------------------------------------------------
#
# Reporting at the T0MEANS state is a default, not a finding, and for many
# nonlinear models it is a poor one -- T0 is often nowhere near where the data
# lives. These resolve the shorthands a caller can pass instead.
#
# Everything here returns the dynamic processes alone -- never the carrier
# states an intoverpop model appends to them. 'mean' and 'asymptotic' are
# statements about where the processes are, and a carrier is not a process: it
# holds an individually varying parameter, so its value belongs to the
# parameter vector the matrices are being built from. The engine fills those
# entries itself, per parameter vector; see `ctsem_parameter_matrices`.
#
# Handing them on instead was a real defect rather than a hypothetical one.
# These vectors are resolved once, from the point estimate, and then reused
# across posterior draws, so every draw was materialised with the point
# estimate's individually varying parameters. Any cell reading one came back
# identical for every draw, and the uncertainty band around it was zero wide
# with nothing said.

.ctContextStateOptions <- c("T0MEANS", "mean", "asymptotic")

# The population T0MEANS over the dynamic processes, the default evaluation
# point.
.ctContextBaseState <- function(fit, tipreds = NULL) {
  # suppressMessages: this is a lookup on the way to choosing an evaluation
  # point, not a report of one, and the note belongs to the caller's choice.
  base <- as.numeric(suppressMessages(
    ctBackendParMatrices(fit, tipreds = tipreds, trim = FALSE))$T0MEANS[, 1])
  base[seq_len(.ctFitNlatent(fit))]
}

# A supplied evaluation point: one entry per latent process, and no more. An
# augmented-length vector is no longer a second accepted spelling of the same
# request -- there is nothing a caller can say about a carrier that the
# parameter vector does not already say.
.ctContextDynamicState <- function(fit, state) {
  # Kept: a resolved 'asymptotic' is passed on as its vector, and each draw's
  # own fixed point is found from it where the draws are evaluated.
  perdraw <- attr(state, "perdraw")
  state <- as.numeric(state)
  nlatent <- .ctFitNlatent(fit)
  if (length(state) != nlatent) {
    stop("state must have ", nlatent, " entries (one per latent process); got ",
      length(state), ".", call. = FALSE)
  }
  attr(state, "perdraw") <- perdraw
  state
}

# Mean smoothed latent state over every row of every subject.
#
# Where the data actually is, as opposed to where the process started. Read
# from the filter rather than derived, so it is the same quantity ctKalman
# plots.
.ctContextMeanState <- function(fit) {
  smoothed <- suppressMessages(ctKalmanArray(fit, pointest = TRUE)$etasmooth)
  # The array may carry the augmented states; take the processes off the front
  # of it rather than asking .ctContextDynamicState to accept both lengths.
  .ctContextDynamicState(fit,
    apply(smoothed, 3L, mean, na.rm = TRUE)[seq_len(.ctFitNlatent(fit))])
}

# The system's own fixed point: the state at which the deterministic change is
# zero, DRIFT(x) x + CINT(x) = 0.
#
# Newton, using the engine's own JAx as the Jacobian rather than finite
# differences -- JAx is exactly d(drift)/d(state), so one engine call per
# iteration does the whole step. For a linear model this lands on asymCINT in a
# single iteration, which is the invariant the test checks. At the estimate,
# from the population T0MEANS.
.ctContextAsymptoticState <- function(fit, tolerance = 1e-8, maxiter = 50L,
  tipreds = NULL) {
  x <- .ctContextAsymptoticStates(fit, matrix(as.numeric(fit$estimate$raw), nrow = 1L),
    start = .ctContextBaseState(fit, tipreds), tipreds = tipreds,
    tolerance = tolerance, maxiter = maxiter)
  if (attr(x, "singular") > 0) {
    stop("No asymptotic state: the system is singular at the current iterate. ",
      "Use state='mean' or supply a state.", call. = FALSE)
  }
  if (attr(x, "unresolved") > 0) {
    stop("No asymptotic state found in ", maxiter, " iterations -- the system may ",
      "have no stable fixed point. Use state='mean' or supply a state.", call. = FALSE)
  }
  as.numeric(x[1L, ])
}

# Each draw's own fixed point, by Newton from `start` (the estimate's), every
# draw's matrices in one engine call an iteration: draws x nlatent, NaN for a
# draw with none near the start, which attribute 'unresolved' counts
# ('singular' those whose Jacobian could not be solved).
.ctContextAsymptoticStates <- function(fit, samples, start, tipreds = NULL,
  time = 0, dt = 0, tolerance = 1e-8, maxiter = 50L) {
  samples <- if (is.matrix(samples)) samples else matrix(as.numeric(samples), nrow = 1L)
  nlatent <- .ctFitNlatent(fit)
  index <- seq_len(nlatent)
  discrete <- !isTRUE(.ctFitModelObject(fit)$continuoustime)
  x <- matrix(as.numeric(start)[index], nrow(samples), nlatent, byrow = TRUE)
  open <- rep(TRUE, nrow(samples))
  singular <- rep(FALSE, nrow(samples))
  for (iteration in seq_len(maxiter)) {
    active <- which(open)
    if (!length(active)) break
    mats <- .ctBackendPopArrays(fit, samples = samples[active, , drop = FALSE],
      tipreds = tipreds, state = t(x[active, , drop = FALSE]), time = time, dt = dt)
    jacobians <- if (is.null(mats$pop_JAx)) mats$pop_DRIFT else mats$pop_JAx
    for (j in seq_along(active)) {
      d <- active[j]
      step <- .ctFixedPointStep(
        drift = matrix(mats$pop_DRIFT[j, index, index], nlatent, nlatent),
        cint = as.numeric(mats$pop_CINT[j, index, 1L]),
        jacobian = matrix(jacobians[j, index, index], nlatent, nlatent),
        x = x[d, ], discrete = discrete, tolerance = tolerance)
      if (is.null(step)) {
        open[d] <- FALSE
        singular[d] <- TRUE
      } else if (!length(step)) {
        open[d] <- FALSE
      } else x[d, ] <- x[d, ] - step
    }
  }
  unresolved <- open | singular
  x[unresolved, ] <- NaN
  attr(x, "unresolved") <- sum(unresolved)
  attr(x, "singular") <- sum(singular)
  x
}

# One Newton step towards DRIFT(x) x + CINT(x) = 0 on the drift's Jacobian:
# numeric(0) where x is there already, NULL where no step can be taken. A
# process whose rate of change is identically zero at x -- a zero row of the
# Jacobian and of the residual, as a latent predictor held constant has --
# has every value as a fixed point, so it is held where it is and the rest are
# solved given it: the ctsem book's latent TI predictor made every such model
# singular here.
.ctFixedPointStep <- function(drift, cint, jacobian, x, discrete, tolerance) {
  n <- length(x)
  if (all(!is.finite(jacobian)) || all(jacobian == 0)) jacobian <- drift
  # Discrete time asks where x = A x + c, i.e. (A - I) x + c = 0.
  if (discrete) {
    drift <- drift - diag(n)
    jacobian <- jacobian - diag(n)
  }
  residual <- as.numeric(drift %*% x) + cint
  if (!all(is.finite(residual)) || !all(is.finite(jacobian))) return(NULL)
  if (max(abs(residual)) < tolerance) return(numeric(0))
  free <- rowSums(abs(jacobian)) > 0 | abs(residual) >= tolerance
  solved <- try(solve(jacobian[free, free, drop = FALSE], residual[free]),
    silent = TRUE)
  if (inherits(solved, "try-error") || !all(is.finite(solved))) return(NULL)
  step <- numeric(n)
  step[free] <- solved
  step
}

# Is `state` asking for anything other than the default evaluation point?
#
# NULL and 'T0MEANS' are the same request, and every function here defaults to
# it, so this is the test for "the caller chose".
.ctContextStateIsDefault <- function(state) {
  if (is.null(state)) return(TRUE)
  if (!is.character(state)) return(FALSE)
  identical(match.arg(state, .ctContextStateOptions), "T0MEANS")
}

# Refuse a caller-chosen evaluation point on a backend that cannot honour it.
#
# Re-materialising the matrices somewhere other than T0MEANS means running the
# model's own cell expressions at a new state, which only the julia engine can
# do; Stan's generated quantities computed theirs once, during sampling, at the
# point the filter was at. Accepting `state=` and returning the T0MEANS matrices
# under another state's name is precisely the failure this file exists to
# prevent, so the argument is refused rather than dropped.
.ctContextRequireStateSupport <- function(fit, state, what = "state=") {
  if (.ctContextStateIsDefault(state)) return(invisible(TRUE))
  if (inherits(fit, "ctJuliaFit")) return(invisible(TRUE))
  stop(call. = FALSE, paste0(what, " needs the model matrices re-materialised at ",
    "that point by the engine that fitted the model, which requires ",
    "backend='julia'. A stan fit reports them at ", .ctContextPopLabel, "."))
}

#' Resolve a state argument to an evaluation point
#'
#' @return list(state = NULL or numeric, one entry per latent process; label =
#'   character).
#'   A NULL state means "the engine's own default", which is T0MEANS: each
#'   draw's own. 'asymptotic' is the estimate's fixed point, marked
#'   (attribute 'perdraw') so that each draw is evaluated at its own, found
#'   from it. 'mean' is the estimate's mean smoothed state for every draw:
#'   a smoother pass a draw would cost far more than the rest of the report.
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
        label = "the estimate's mean smoothed latent state"))
    }
    point <- .ctContextAsymptoticState(fit, tipreds = tipreds)
    attr(point, "perdraw") <- "asymptotic"
    return(list(state = point,
      label = "the system's asymptotic (fixed point) state"))
  }
  list(state = .ctContextDynamicState(fit, state), label = "the supplied state")
}
