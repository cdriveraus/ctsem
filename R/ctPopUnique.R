# Own variances on a reduced-rank level: `ctFit(..., popunique =)`.
#
# A level's random-effect covariance under `poprank = r` is `L L'`, `L` a
# `k x r` loading matrix. `popunique` adds an own variance to the effects it
# selects, so the covariance becomes `L L' + diag(psi^2)` -- the
# factor-analytic form, shared dimensions plus what is specific to each
# effect. With `poprank = 0` and every effect selected it is a diagonal
# covariance: independent random effects.
#
# Why it exists: person-level item intercepts are mostly shared (one level per
# affect domain) and partly specific to the item. A pure low rank can hold the
# specific part only one shared dimension at a time, so each added dimension
# keeps buying likelihood, while a full covariance leaves the covariance of
# every pair of items no person answers together to the prior. One own variance
# per item says what the data support: each is informed by everyone who
# answers that item, and cross-item covariance comes through `L` alone.
#
# In the engine an own variance is one column of the level's root,
# `F = [L | Psi]`, with a single entry in its effect's row
# (`_laplace_poploading`), transformed as a full-rank level's sd. It is built
# where the loadings are, so it reaches `'laplace'` and `'none'` (which prepare
# the same engine specification) and nothing else: the augmented route would
# need a carrier state per own variance, and stan builds no loading matrix.
# A full-rank level is refused, since an own variance on top of an unrestricted
# covariance is not identified.

# TRUE when `popunique` asks for anything: NULL and FALSE (in any shape) do not.
.ctPopuniqueRequested <- function(popunique) {
  if (is.null(popunique)) return(FALSE)
  values <- if (is.list(popunique)) unlist(popunique, use.names = FALSE) else popunique
  !(is.logical(values) && all(!values %in% TRUE))
}

# `popunique` per grouping level, as `.ctPoprankByLevel()` resolves `poprank`:
# one unnamed value applies to every level, a named one to the level it names.
# Each value is TRUE (every effect varying at that level), FALSE, or the names
# of matrices or parameters to select. A named atomic vector may name a level
# more than once, `c(subject = 'MANIFESTMEANS', subject = 'DIFFUSION')`, and a
# list may hold several names per level, `list(subject = c(...))`.
.ctPopuniqueByLevel <- function(popunique, levelnames) {
  levels <- paste(levelnames, collapse = ', ')
  named <- names(popunique)
  if (is.null(named)) {
    value <- if (is.list(popunique)) unlist(popunique, use.names = FALSE) else popunique
    return(stats::setNames(rep(list(value), length(levelnames)), levelnames))
  }
  if (any(is.na(named) | !nzchar(named))) stop(
    "every entry of a per-level popunique must name its level: ", levels,
    call. = FALSE)
  unknown <- setdiff(unique(named), levelnames)
  if (length(unknown)) stop("popunique names no level called ",
    paste(unknown, collapse = ', '), ". The levels are ", levels, ".",
    call. = FALSE)
  values <- if (is.list(popunique)) popunique else as.list(popunique)
  out <- list()
  for (i in seq_along(values)) {
    level <- named[i]
    value <- values[[i]]
    if (is.null(out[[level]])) out[[level]] <- value else {
      if (is.logical(out[[level]]) || is.logical(value)) stop("popunique gives level '",
        level, "' more than one entry, and one of them is TRUE or FALSE.", call. = FALSE)
      out[[level]] <- unique(c(out[[level]], value))
    }
  }
  out
}

# The parameters one level's `popunique` value selects, from those varying at
# that level (`column` is its flag: `indvarying`, or `indvarying_<level>`).
.ctPopuniqueSelect <- function(pars, column, value, level) {
  pars <- .ctModelCleanctspec(pars)
  flag <- if (column %in% names(pars)) pars[[column]] else rep(FALSE, nrow(pars))
  rows <- !is.na(pars$param) & is.na(pars$value) & flag %in% TRUE &
    !grepl('[][]', as.character(pars$param))
  varying <- unique(data.frame(param = as.character(pars$param[rows]),
    matrix = as.character(pars$matrix[rows]), stringsAsFactors = FALSE))
  varying <- varying[!duplicated(varying$param), , drop = FALSE]
  if (is.logical(value)) {
    if (length(value) != 1L || is.na(value)) stop("popunique for level '", level,
      "' must be TRUE, FALSE, or names of matrices or parameters.", call. = FALSE)
    return(if (isTRUE(value)) varying$param else character())
  }
  value <- as.character(value)
  known <- value %in% c(varying$matrix, varying$param)
  if (any(!known)) stop("popunique for level '", level, "' names ",
    paste(value[!known], collapse = ', '), ", which ",
    if (sum(!known) == 1L) "is" else "are",
    " neither a matrix nor a parameter varying at that level.", call. = FALSE)
  varying$param[varying$matrix %in% value | varying$param %in% value]
}

# Resolve and check `popunique` against the ranks `poprank` gave
# (`laplacerank`, a named integer vector per level, or NULL). Returns a named
# list, level -> the parameters given an own variance, holding only levels
# with any.
.ctPopuniqueResolve <- function(ctm, popunique, levelnames, levelcolumns,
  laplacerank) {
  asked <- .ctPopuniqueByLevel(popunique, levelnames)
  out <- list()
  for (level in names(asked)) {
    selected <- .ctPopuniqueSelect(ctm$pars, levelcolumns[[level]],
      asked[[level]], level)
    if (!length(selected)) next
    k <- nrow(.ctPopEffectRoles(ctm$pars, column = levelcolumns[[level]]))
    rank <- if (!is.null(laplacerank) && level %in% names(laplacerank))
      laplacerank[[level]] else NA_integer_
    if (is.na(rank) || rank >= k) stop("popunique gives own variances at level '",
      level, "', which needs a reduced rank there: an own variance on top of an ",
      "unrestricted covariance is not identified. Give poprank=c(", level, "=r) ",
      "with r below its ", k, " varying parameters, or r=0 for independent ",
      "random effects.", call. = FALSE)
    if (rank == 0L && length(selected) < k) {
      lacking <- setdiff(.ctPopEffectRoles(ctm$pars,
        column = levelcolumns[[level]])$param, selected)
      stop("poprank=0 at level '", level, "' leaves only own variances, so every ",
        "parameter varying there needs one; these have none: ",
        paste(lacking, collapse = ', '), ". Select them in popunique, or give a ",
        "rank of 1 or more.", call. = FALSE)
    }
    out[[level]] <- selected
  }
  out
}
