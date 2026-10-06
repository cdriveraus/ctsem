# Reduced-rank population covariance on the augmented route, as loadings on
# standardised dimensions written into the model around its carrier states.
# `'laplace'` and `'none'` build the same structure in the engine instead, per
# level (`.ctJuliaLevelRank()`, `_laplace_poploading`).
#
# ## The problem this solves
#
# Under `intoverpop='augmented'` an individually varying parameter becomes a
# static latent state reaching the likelihood only through the model cell it
# drives. For a cell entering the observation *mean* the Kalman update moves
# that state, so its population variance is informed. For a **variance** cell
# -- DIFFUSION, MANIFESTVAR -- the filter builds the cell from the state's mean
# only, the observation mean function's Jacobian with respect to that state is
# zero, and the update can never move it. The data then learns about such an
# effect only through its correlation with states the filter *can* update, so
# its covariance with those is identified and the split of that covariance into
# a standard deviation and correlations is not.
#
# Measured on a one-latent model with a DRIFT and a DIFFUSION random effect: the
# profile likelihood in the correlation coordinate is flat to 12 significant
# figures across `r` = 0.60 to 0.998, with the population sd compensating to
# hold `sd * r` constant to four digits. See
# `review/RANDOMEFFECTS-partial-identification-2026-09-07.md` and
# `review/REDUCEDRANK-random-effects-2026-09-07.md`.
#
# ## The structure
#
# Partition the `k` varying parameters into a **basis** `A` (`r` of them) and
# the **regressed** rest `B`. Then
#
#     theta_B = mean_B + beta theta_A
#
# so the population covariance is
#
#     Sigma = [[Sigma_AA,          Sigma_AA beta'],
#              [beta Sigma_AA, beta Sigma_AA beta']]
#
# which is rank `r` by construction. The free parameters are `Sigma_AA` in
# ctsem's existing sd/correlation form -- untouched, so its transforms, priors
# and summary rows all keep working -- plus the `|B| x r` coefficients `beta`.
# The count is `r(r+1)/2 + |B|*r`, and at `r = k` that is `k(k+1)/2`, today's
# full-rank count, so this is a strict generalisation rather than a mode.
#
# ## Coordinates: loadings, not the regression this file is named for
#
# The basis/regression partition above is how the rank and the basis are
# chosen. The fitted coordinates have been loadings since `203e223f`:
# `Sigma = L L'`, `L` k-by-r and triangular in its first r rows, on carrier
# states of unit variance -- the basis rows `M` and the regressed rows are one
# loading matrix, named `poploading_<param>_dim<j>` as the engine names it on
# the laplace and 'none' routes. The sd,
# correlation and coefficient coordinates it replaced spanned the same manifold
# but reached the optimum from 1 of 12 matched starts against 8 of 12 (see that
# commit and review/REDUCEDRANK-parameterisation-reliability-2026-09-08.md).
#
# ## What it decides rather than estimates
#
# A regressed effect has no variance independent of the basis -- the residual
# variance is fixed at zero. When `r` is the number of mean-affecting effects
# that residual is exactly the unidentified quantity, so the restriction costs
# no likelihood (measured: `+0.000000` on two models). Below that it is a real
# approximation which will inflate `beta` and distort `Sigma_AA`, because the
# fit is joint and nothing holds the basis block fixed.

# Which matrices carry a cell into the observation mean.
#
# DRIFT belongs here because the predicted observation depends on the state it
# propagates, so the Jacobian of the observation with respect to a DRIFT
# carrier state is proportional to the (data-driven, nonzero) state estimate --
# not because the process mean is nonzero. THRESHOLDS shifts an ordinal
# indicator's expectation, so it counts too. T0VAR is *not* here: an indvarying
# T0VAR cell is made redundant by `T0VARredundancies()` on the augmented route.
.ctPopMeanMatrices <- function() {
  c('T0MEANS', 'LAMBDA', 'DRIFT', 'MANIFESTMEANS', 'CINT', 'TDPREDEFFECT',
    'THRESHOLDS')
}

.ctPopVarianceMatrices <- function() c('DIFFUSION', 'MANIFESTVAR', 'T0VAR')

# Labels sitting in T0MEANS. Read before `.ctModelIntOverPop()` runs, while the
# cell still carries the label rather than a state reference.
.ctPopT0meansEffects <- function(pars) {
  pars <- .ctModelCleanctspec(pars)
  rows <- pars$matrix %in% 'T0MEANS' & !is.na(pars$param)
  unique(as.character(pars$param[rows]))
}

# Word-boundary regex for a literal parameter label.
.ctPopLabelPattern <- function(label) {
  paste0('(^|[^[:alnum:]_.])', gsub('([][{}()+*^$|\\\\?.])', '\\\\\\1', label),
    '($|[^[:alnum:]_.])')
}

# Which random effects reach the observation mean, and which do not.
#
# Returns one row per individually varying free parameter with a logical
# `mean` column. A parameter counts as mean-affecting if it appears in a cell
# of a mean-affecting matrix, following `PARS` indirection transitively -- a
# DIFFUSION parameter also referenced inside a DRIFT expression *is*
# mean-affecting and identified, and that is the case a pattern match over
# `$pars$matrix` alone gets wrong (see the review note's section 3).
#
# `pars` must already have been through `ctModelStatesAndPARS()`, so that
# cross-cell label references appear as `PARS[r,c]`, and must not yet have been
# through `.ctModelIntOverPop()`, which clears `indvarying` on the cells it
# rewrites into state references.
.ctPopEffectRoles <- function(pars, column = 'indvarying') {
  pars <- .ctModelCleanctspec(pars)
  free <- !is.na(pars$param) & is.na(pars$value)
  # `column` is the level's own flag: `indvarying` for the innermost level and
  # `indvarying_<idname>` for each one above it. Absent means no effects at
  # that level, which is a level with nothing to restrict rather than an error.
  flag <- if (column %in% names(pars)) pars[[column]] else rep(FALSE, nrow(pars))
  varying <- free & !is.na(flag) & flag
  labels <- unique(as.character(pars$param[varying]))
  # A rewritten cell (`state[3]`, `PARS[1,1] * 2`) is not a label.
  labels <- labels[!grepl('[][]', labels)]
  if (!length(labels)) {
    return(data.frame(param = character(), mean = logical(),
      stringsAsFactors = FALSE))
  }
  meanmats <- .ctPopMeanMatrices()
  text <- as.character(pars$param)
  text[is.na(text)] <- ''
  matrixname <- as.character(pars$matrix)

  # Rows referencing a given piece of text, and whether any of them is in a
  # mean-affecting matrix. Follows PARS cells onwards: a label sitting in
  # PARS[r,c] is mean-affecting if anything referencing `PARS[r,c]` is.
  reaches <- function(pattern, seen) {
    if (pattern %in% seen) return(FALSE)
    seen <- c(seen, pattern)
    hits <- grep(pattern, text)
    if (!length(hits)) return(FALSE)
    if (any(matrixname[hits] %in% meanmats)) return(TRUE)
    parsrows <- hits[matrixname[hits] %in% 'PARS']
    for (ri in parsrows) {
      ref <- sprintf('PARS\\[\\s*%d\\s*,\\s*%d\\s*\\]', pars$row[ri], pars$col[ri])
      if (reaches(ref, seen)) return(TRUE)
    }
    FALSE
  }

  data.frame(param = labels,
    mean = vapply(labels, function(l) reaches(.ctPopLabelPattern(l), character()),
      logical(1L)),
    row.names = NULL, stringsAsFactors = FALSE)
}

# Every RAWPOPVAR cell a user has stated, whatever it says.
#
# RAWPOPVAR is the specification surface for the population covariance
# (R/ctModelRawPopVar.R): a label estimates a cell and a number fixes it, and a
# label differing from the default is an equality constraint, so both are
# statements.
#
# A reduced rank honours none of them, and the reason is the parameterisation
# rather than the cell. The covariance is `Sigma = L L'` on standardised
# dimensions: below the diagonal `Sigma[i,j]` is `sum_{m<=min(i,j)} L[i,m]
# L[j,m]` and not a cell, an effect's spread is a row norm and not a cell, and a
# regressed effect has neither a spread nor correlations of its own. So the
# reduction requires the covariance free and says so when it is not.
#
# Sorting statements into the expressible and the rest was tried and removed. It
# accepted a fixed sd for the first basis effect, whose row really is one
# loading, and a zero anywhere, on the grounds that a zero covariance is exact
# -- which is true only against the first dimension, and which nothing
# implemented in either case, so the zero was accepted and then silently
# dropped.
.ctPopRegressionRawPopVarStated <- function(model) {
  popcov <- model[['RAWPOPVAR']]
  if (is.null(popcov) || !length(popcov)) return(character())
  default <- .ctModelRawPopVar(model$pars)
  if (is.null(default) || !length(default)) return(character())
  # Only parameters both matrices have. A RAWPOPVAR built for a different set
  # of varying parameters -- `pars$indvarying` edited after `ctModel()`, which
  # the tests here do -- is stale rather than stated, and refusing on it would
  # refuse a model nobody constrained. Each pair is read where the stored
  # matrix keeps it (`.ctModelRawPopVarCell()`), so a reorder of `pars` since
  # is not read as a statement either. Only a number states anything: a label,
  # whatever it says, leaves the cell free on every route.
  names <- rownames(default)[rownames(default) %in% rownames(popcov) &
    rownames(default) %in% colnames(popcov)]
  out <- character()
  for (i in seq_along(names)) for (j in seq_len(i)) {
    stated <- .ctModelRawPopVarCell(popcov, names[i], names[j])
    if (!is.finite(.ctModelRawPopVarValue(stated))) next
    out <- c(out, sprintf("RAWPOPVAR['%s', '%s'] = %s", names[i], names[j], stated))
  }
  out
}

# What `poprank` asks of each level: a list named by level, each entry
# `'auto'`, a whole number, or `NA` for that level at full rank. One unnamed
# value applies to every level in `levelnames`; a named vector to the levels it
# names, each of which must be one of them.
#
# Names decide, not length. `c(study = 2)` is one element *and* names a level,
# and reading it as "2 everywhere" would reduce every level while looking as if
# it had done what was asked -- which is what the first version of the laplace
# route did. A name that is not a level is refused for the same reason: as a
# no-op it would be a request silently not honoured. And each entry is read on
# its own, so `c(study = 'auto')` resolves the study level alone and
# `c(subject = 'auto', study = 2)` keeps the 2.
.ctPoprankByLevel <- function(poprank, levelnames) {
  levels <- paste(levelnames, collapse = ', ')
  named <- names(poprank)
  if (is.null(named)) {
    if (length(poprank) != 1L) stop(
      "a poprank per level must be named, one entry per level: ", levels,
      call. = FALSE)
    named <- levelnames
    poprank <- rep(poprank, length(levelnames))
  } else {
    if (any(is.na(named) | !nzchar(named))) stop(
      "every entry of a per-level poprank must name its level: ", levels,
      call. = FALSE)
    unknown <- setdiff(named, levelnames)
    if (length(unknown)) stop("poprank names no level called ",
      paste(unknown, collapse = ', '), ". The levels are ", levels, ".",
      call. = FALSE)
    twice <- unique(named[duplicated(named)])
    if (length(twice)) stop("poprank names ", paste(twice, collapse = ', '),
      " more than once.", call. = FALSE)
  }
  out <- lapply(as.list(poprank), function(value) {
    if (length(value) == 1L && is.na(value)) return(NA_integer_)
    if (identical(as.character(value), 'auto')) return('auto')
    number <- suppressWarnings(as.numeric(value))
    if (length(number) != 1L || is.na(number) || number != round(number)) stop(
      "poprank must be NA, 'auto', or a whole number at each level.",
      call. = FALSE)
    as.integer(number)
  })
  stats::setNames(out, named)
}

# Resolve `poprank` into a basis and a set of regressed effects.
#
# `poprank` is `NA` (no restriction, today's full-rank behaviour), `'auto'`
# (the number of mean-affecting effects, which is the largest rank the
# augmented route can identify and is a no-op when every effect is
# mean-affecting), or an integer `r` (an explicit approximation).
#
# The basis is chosen mean-affecting effects first, then in the order they
# appear. Which effects form the basis changes the coordinates but not the
# model -- any `r` effects spanning the same rank-`r` space describe the same
# population covariance -- so this is a conditioning choice, not a modelling
# one. Putting the identified effects first is what makes the retained
# coordinates the identified ones when `poprank='auto'`.
.ctPopRegressionSpec <- function(pars, poprank, explicit = TRUE, model = NULL) {
  if (is.null(poprank) || (length(poprank) == 1L && is.na(poprank))) return(NULL)
  roles <- .ctPopEffectRoles(pars)
  if (!nrow(roles)) return(NULL)

  k <- nrow(roles)
  nmean <- sum(roles$mean)
  rank <- if (identical(poprank, 'auto')) nmean else {
    rank <- suppressWarnings(as.integer(poprank))
    if (is.na(rank)) stop("poprank must be NA, 'auto', or a whole number.",
      call. = FALSE)
    rank
  }
  if (rank < 0L) stop('poprank cannot be negative.', call. = FALSE)
  if (rank > k) {
    stop('poprank is ', rank, ' but the model has only ', k,
      ' individually varying parameters.', call. = FALSE)
  }
  if (rank == 0L) {
    # Nothing reaches the observation mean, so under `intoverpop='augmented'`
    # no part of the population covariance is identified and there is no basis
    # to regress on. Refuse when the rank was asked for; when it is only the
    # default, leave the model alone -- the pre-fit warning in `ctFit()` already
    # names these parameters, and turning a call that used to run into an error
    # is not this argument's job.
    if (!isTRUE(explicit)) return(NULL)
    stop('This model has no individually varying parameter that reaches the ',
      'observation mean, so nothing about the population covariance is ',
      "identified under intoverpop='augmented'. Use intoverpop='laplace', ",
      'or give the model a random effect on a mean-affecting matrix.',
      call. = FALSE)
  }
  order <- order(!roles$mean, seq_len(k))
  basis <- roles$param[order][seq_len(rank)]
  regressed <- setdiff(roles$param[order], basis)

  # An individually varying T0MEANS cannot be a loading on the standardised
  # dimensions. T0MEANS reaches the observation mean, so such an effect is
  # always in the basis; and its carrier state is the model latent itself
  # rather than an appended one, so its cell reads its own parameter and there
  # is no `state[<carrier>]` for `.ctPopRegressionRewrite()` to substitute a
  # predictor into. Giving it a dimension of its own needs a T0MEANS cell that
  # reads another state, which is the thing the augmented arrangement exists to
  # avoid, so this is a limit of the parameterisation rather than an oversight.
  # `'laplace'` and `'none'` build their loadings in the engine and have no
  # carrier states, so they are not restricted here.
  # Only where the reduction reaches a cell: with nothing regressed the rewrite
  # returns the model untouched, so a rank equal to the number of varying
  # parameters is a no-op and there is nothing to refuse.
  if (length(regressed)) {
    t0basis <- intersect(basis, .ctPopT0meansEffects(pars))
    if (length(t0basis)) {
      if (!isTRUE(explicit)) return(NULL)
      stop('poprank cannot reduce a population covariance that includes an ',
        'individually varying T0MEANS (', paste(t0basis, collapse = ', '),
        '): its carrier state is the latent itself, so it has no dimension of ',
        'its own to load on. Use poprank=NA to estimate the full covariance, ',
        "or intoverpop='laplace'.", call. = FALSE)
    }
  }

  # The population covariance has to be free. Asked for, a statement is an
  # error naming the cells; defaulted, the user's own specification is the more
  # explicit of the two and the rank is simply not applied.
  if (!is.null(model)) {
    stated <- .ctPopRegressionRawPopVarStated(model)
    if (length(stated)) {
      if (!isTRUE(explicit)) return(NULL)
      stop('poprank needs a free population covariance, and RAWPOPVAR states: ',
        paste(stated, collapse = '; '),
        '. Under a reduced rank the covariance is a factor, where a standard ',
        'deviation is a row norm over dimensions and a covariance is a sum ',
        'over them, so none of these is a cell that can be fixed. Use ',
        'poprank=NA to estimate the full covariance, or leave RAWPOPVAR ',
        'alone.', call. = FALSE)
    }
  }
  # Above `nmean` the restriction stops being free: it starts fixing residual
  # variances the data does determine. Said once, here, rather than left for a
  # user to infer from a likelihood that moved.
  approximate <- rank < nmean
  # Where each basis effect's cells are, as coordinates, recorded now because
  # they cannot be found later. `.ctModelIntOverPop()` replaces the `param`
  # token in a basis cell with `state[j]`, so after it runs the cell's text no
  # longer mentions the effect at all -- and it rebuilds the T0VAR block, so a
  # row index taken here would not survive either.
  #
  # Every row carrying the label, which is what the regressed path also takes:
  # an effect can drive more than one cell (a DRIFT entry and its JAx mirror,
  # or a label used twice), and all of them read the same carrier state.
  basiscells <- do.call(rbind, lapply(basis, function(b) {
    rows <- which(!is.na(pars$param) & pars$param %in% b)
    if (!length(rows)) return(NULL)
    data.frame(param = b, matrix = as.character(pars$matrix[rows]),
      row = as.integer(pars$row[rows]), col = as.integer(pars$col[rows]),
      row.names = NULL, stringsAsFactors = FALSE)
  }))
  list(rank = rank, basis = basis, regressed = regressed, roles = roles,
    nmean = nmean, approximate = approximate, basiscells = basiscells,
    npar = rank * (rank + 1L) / 2L + length(regressed) * rank)
}

# What to tell the user when the rank was dropped.
#
# `poprank='auto'` is the default, so a model can lose population parameters
# without the user having asked, and the message has to carry three things: how
# much was dropped, why those particular effects, and how to get the other
# behaviour back. Kept to a few lines -- the reasoning belongs in the comment at
# the top of this file, not in every fit's output.
#
# The two reasons are reported separately, because they are not the same claim.
# An effect that never reaches the observation mean is regressed because its
# own spread is *not identified* -- nothing is lost. An effect that does reach
# the mean and is regressed anyway was dropped to meet a rank the user asked
# for, and that does lose something. A first version of this message explained
# every regressed effect with the identification reason and so told a user that
# `dr2` and `dr3` "vary only where the filter cannot see their spread", which
# is untrue of a DRIFT effect and is exactly the kind of plausible wrong
# statement this whole feature exists to remove.
.ctPopRegressionMessage <- function(spec) {
  regressed <- unique(spec$coefficients$param)
  roles <- spec$roles
  variancecell <- intersect(regressed, roles$param[!roles$mean])
  demoted <- intersect(regressed, roles$param[roles$mean])
  basis <- paste(spec$basis, collapse = ', ')
  out <- paste0('poprank: population covariance reduced to rank ', spec$rank,
    ' of ', spec$rank + length(regressed), '.')
  # The augmented route only: there the filter cannot see a variance cell's
  # own spread. (`'laplace'` and `'none'` identify those coordinates; `ctFit()`
  # says their reduction is an approximation where it applies the rank.)
  if (length(variancecell)) {
    out <- paste0(out, ' ', paste(variancecell, collapse = ', '),
      if (length(variancecell) > 1) ' vary' else ' varies',
      ' only in DIFFUSION / MANIFESTVAR, where the augmented filter cannot see',
      if (length(variancecell) > 1) ' their' else ' its', ' own spread, so ',
      if (length(variancecell) > 1) 'they are' else 'it is',
      ' estimated as a regression on ', basis, '.')
  }
  if (length(demoted)) {
    out <- paste0(out, ' ', paste(demoted, collapse = ', '),
      if (length(demoted) > 1) ' are' else ' is',
      ' also regressed on ', basis, ' to meet the requested rank, which is ',
      'below the ', spec$nmean, ' this model identifies -- an approximation, ',
      'and the retained parameters absorb what it drops.')
  }
  paste0(out, " poprank=NA estimates the full covariance instead;",
    " intoverpop='laplace' identifies it.")
}

# Turn the regressed effects off before the augmentation runs.
#
# `.ctModelIntOverPop()` creates one carrier state per individually varying
# parameter, so clearing the flag on the regressed ones is all it takes to get
# carrier states for the basis alone. Nothing else about that function changes,
# which is the point: the regression form adds a rewrite either side of it
# rather than a second augmentation path through it.
.ctPopRegressionDemote <- function(m, spec) {
  m$pars$indvarying[!is.na(m$pars$param) & m$pars$param %in% spec$regressed] <- FALSE
  m
}

# Write each regressed effect as its own mean plus a regression on the basis
# effects' carrier states.
#
# Runs *after* `.ctModelIntOverPop()`, so the basis effects already have their
# states and the coefficients can reference them directly, and *before* the
# second `ctModelStatesAndPARS()` call in `ctFit()`, so the new mean and
# coefficient parameters can be written as plain labels and be turned into
# `PARS[r,c]` references by the machinery that already does that.
#
# The basis state index comes from the T0MEANS row carrying the label, which is
# uniform across the two kinds of basis effect: an indvarying T0MEANS keeps its
# own state row, and every other effect gets an appended one whose T0MEANS row
# retains the original label.
.ctPopRegressionRewrite <- function(m, spec) {
  if (!length(spec$regressed)) return(m)
  t0 <- m$pars$matrix %in% 'T0MEANS' & m$pars$col %in% 1
  stateof <- vapply(spec$basis, function(b) {
    idx <- which(t0 & !is.na(m$pars$param) & m$pars$param %in% b)
    if (!length(idx)) NA_integer_ else as.integer(m$pars$row[idx[1L]])
  }, integer(1L))
  if (anyNA(stateof)) {
    stop('Internal error: no carrier state for basis random effect(s) ',
      paste(spec$basis[is.na(stateof)], collapse = ', '), '.', call. = FALSE)
  }

  template <- m$pars[1L, , drop = FALSE]
  effectcols <- grep('_effect$', names(m$pars), value = TRUE)
  newpars <- list()
  parsrow <- suppressWarnings(max(c(0L,
    as.integer(m$pars$row[m$pars$matrix %in% 'PARS']))))
  addpar <- function(label, effects = NULL) {
    parsrow <<- parsrow + 1L
    row <- template
    row$matrix <- 'PARS'; row$row <- parsrow; row$col <- 1L
    row$param <- label; row$value <- NA; row$transform <- 'param'
    row$indvarying <- FALSE
    if (length(effectcols)) {
      row[, effectcols] <- FALSE
      if (!is.null(effects)) row[, effectcols] <- effects
    }
    newpars[[length(newpars) + 1L]] <<- row
    invisible(NULL)
  }

  coefficients <- list()
  drivencells <- list()
  nbefore <- nrow(m$pars)
  for (p in spec$regressed) {
    # TI-predictor effects follow the parameter's *mean*, which is where they
    # already acted: a TI effect shifts a subject's raw parameter value, and
    # under the augmented route that means the population mean rather than the
    # random deviation -- `.ctModelIntOverPop()` does the same thing for a basis
    # effect, copying the flags onto the carrier state's T0MEANS row and
    # clearing them on the cell it rewrites. `.ctJuliaTIEffects()` attaches an
    # effect to any free parameter whose `param` is a plain label, which the
    # mean is and the rewritten cell is not.
    effects <- NULL
    if (length(effectcols)) {
      own <- which(!is.na(m$pars$param) & m$pars$param %in% p)
      if (length(own)) effects <- vapply(effectcols,
        function(cc) any(m$pars[own, cc] %in% TRUE), logical(1L))
    }
    # The mean keeps the parameter's own name, so the summary still has a row
    # called `df11` meaning the population mean of `df11`.
    addpar(p, effects)
    # Named as the engine names the same loading (`poploading_<param>_dim<j>`),
    # since it is the same coordinate: this row's loading on basis effect j's
    # dimension, the carrier state of `spec$basis[j]`.
    betas <- paste0('poploading_', p, '_dim', seq_along(spec$basis))
    for (b in betas) addpar(b)
    coefficients[[length(coefficients) + 1L]] <- data.frame(
      param = p, basis = spec$basis, coefficient = betas,
      state = as.integer(stateof), row.names = NULL, stringsAsFactors = FALSE)
    predictor <- paste0('(', p, ' + ',
      paste0(betas, ' * state[', stateof, ']', collapse = ' + '), ')')
    cells <- which(!is.na(m$pars$param) & m$pars$param %in% p)
    # The rows just added carry the same label, so exclude anything added here.
    # `addpar()` accumulates into `newpars` and appends after the loop, so
    # nothing added here is in `m$pars` yet and `nbefore` is what says so.
    cells <- cells[cells <= nbefore]
    if (!length(cells)) {
      stop('Internal error: no cell found for regressed random effect ', p, '.',
        call. = FALSE)
    }
    for (ri in cells) {
      drivencells[[length(drivencells) + 1L]] <- data.frame(
        param = p, matrix = as.character(m$pars$matrix[ri]),
        row = as.integer(m$pars$row[ri]), col = as.integer(m$pars$col[ri]),
        row.names = NULL, stringsAsFactors = FALSE)
      transform <- as.character(m$pars$transform[ri])
      if (is.na(transform) || !nzchar(transform)) transform <- 'param'
      # Lookarounds rather than `\b`: `param` has to be substituted as a whole
      # token, and the plain `gsub('param', ...)` the augmentation uses next
      # door would also rewrite the `param` inside a longer name.
      m$pars$param[ri] <- gsub('(?<![[:alnum:]_.])param(?![[:alnum:]_.])',
        predictor, transform, perl = TRUE)
      m$pars$transform[ri] <- NA
      m$pars$indvarying[ri] <- FALSE
      # Cleared here as well as carried above: a TI-predictor flag left on a
      # cell whose `param` is now an expression is the shape that made
      # `.ctJuliaTIEffects()` mint a coefficient with nothing to attach it to.
      if (length(effectcols)) m$pars[ri, effectcols] <- FALSE
    }
  }
  # --- the basis effects, as loadings on the standardised dimensions --------
  #
  # Each basis cell reads `tf(state[j])` by now, and the substitution target is
  # that one token. Addressed by the coordinates `.ctPopRegressionSpec()`
  # recorded before the augmentation, because the cell's text no longer
  # mentions the effect and the regressed cells contain the same token and must
  # keep it.
  loadings <- list()
  for (pos in seq_along(spec$basis)) {
    p <- spec$basis[pos]
    cells <- spec$basiscells[spec$basiscells$param %in% p, , drop = FALSE]
    if (!nrow(cells)) {
      stop('Internal error: no recorded cell for basis random effect ', p, '.',
        call. = FALSE)
    }
    effects <- NULL
    if (length(effectcols)) {
      own <- which(!is.na(m$pars$param) & m$pars$param %in% p)
      if (length(own)) effects <- vapply(effectcols,
        function(cc) any(m$pars[own, cc] %in% TRUE), logical(1L))
    }
    # The mean keeps the effect's own name, as a regressed effect's does, so a
    # summary still has a row called `dr1` meaning the population mean of dr1.
    # It moves here from the carrier's T0MEANS row, which is fixed to zero
    # below.
    addpar(p, effects)
    # Triangular: effect at position `pos` loads on dimensions 1..pos. That is
    # what removes the rotation freedom a free loading matrix would have, and
    # it keeps the parameter count at r(r+1)/2 over the basis block, which is
    # what the sd-and-correlation block had.
    own_loadings <- paste0('poploading_', p, '_dim', seq_len(pos))
    # `sdscale` multiplied this effect's population sd under the previous
    # parameterisation. Its spread is a row norm of the loading matrix now, so
    # the scale goes on the row -- which scales the spread by exactly the same
    # factor. Omitting it would leave the argument accepted and ignored.
    own_scale <- suppressWarnings(as.numeric(
      m$pars$sdscale[!is.na(m$pars$param) & m$pars$param %in% p])[1L])
    if (!is.finite(own_scale) || own_scale == 0) own_scale <- 1
    for (l in own_loadings) {
      addpar(l)
      if (own_scale != 1) {
        newpars[[length(newpars)]]$transform <-
          sprintf('%.17g * param', own_scale)
      }
    }
    loadings[[length(loadings) + 1L]] <- data.frame(
      param = p, dimension = seq_len(pos), loading = own_loadings,
      state = as.integer(stateof[seq_len(pos)]),
      row.names = NULL, stringsAsFactors = FALSE)
    predictor <- paste0('(', p, ' + ', paste0(own_loadings, ' * state[',
      stateof[seq_len(pos)], ']', collapse = ' + '), ')')
    target <- paste0('state[', stateof[pos], ']')
    for (ci in seq_len(nrow(cells))) {
      ri <- which(m$pars$matrix %in% cells$matrix[ci] &
          m$pars$row == cells$row[ci] & m$pars$col == cells$col[ci])
      if (length(ri) != 1L) {
        stop('Internal error: basis cell ', cells$matrix[ci], '[',
          cells$row[ci], ',', cells$col[ci], '] for ', p,
          ' is not uniquely locatable after augmentation.', call. = FALSE)
      }
      text <- as.character(m$pars$param[ri])
      if (is.na(text) || !grepl(target, text, fixed = TRUE)) {
        stop('Internal error: basis cell ', cells$matrix[ci], '[',
          cells$row[ci], ',', cells$col[ci], '] does not read ', target,
          ' as expected; found ', if (is.na(text)) 'NA' else text, '.',
          call. = FALSE)
      }
      m$pars$param[ri] <- gsub(target, predictor, text, fixed = TRUE)
      m$pars$indvarying[ri] <- FALSE
      if (length(effectcols)) m$pars[ri, effectcols] <- FALSE
    }
    # And the carrier becomes the dimension: mean zero, so `state[i]` is the
    # standardised deviation and nothing else.
    ti <- which(m$pars$matrix %in% 'T0MEANS' & m$pars$row == stateof[pos] &
        m$pars$col == 1)
    if (length(ti) != 1L) {
      stop('Internal error: no unique carrier T0MEANS row for ', p, '.',
        call. = FALSE)
    }
    m$pars$param[ti] <- NA_character_
    m$pars$value[ti] <- 0
    m$pars$transform[ti] <- NA_character_
    if (length(effectcols)) m$pars[ti, effectcols] <- FALSE
  }

  m$pars <- rbind(m$pars, do.call(rbind, newpars))
  m$pars[] <- lapply(m$pars, utils::type.convert, as.is = TRUE)
  spec$coefficients <- do.call(rbind, coefficients)
  spec$loadings <- do.call(rbind, loadings)
  spec$cells <- do.call(rbind, drivencells)
  spec$route <- 'augmented'
  spec$state <- stateof
  # The population block is the identity under this form, which the julia
  # augmentation reads to fix it rather than estimate it.
  spec$standardised <- TRUE
  m$popregression <- spec
  m
}
