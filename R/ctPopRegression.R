# Reduced-rank population covariance, as regression on a basis of random
# effects.
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
# ## Why regression coordinates and not a loading matrix
#
# `Sigma = L L'` with `L = [[L_A], [beta L_A]]` is a `k x r` lower-trapezoidal
# loading matrix, so this *is* the loading form -- but reparameterised, and the
# reparameterisation is the point. A free loading matrix has `r(r-1)/2` rotation
# flat directions needing a triangular restriction, and its loadings can all
# approach zero together, which is locally flat. Measured on a `k = 6`, `r = 3`
# model: `ctIdentify` reports `nweak` 0 at every evaluation point for these
# coordinates against 0..15 for free loadings, and the fitted likelihood matches
# full rank to 5.7e-9 against 1.3e-5. Here the basis effects carry their own
# identified spreads, so there is no rotation freedom and nothing to normalise.
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
.ctPopEffectRoles <- function(pars) {
  pars <- .ctModelCleanctspec(pars)
  free <- !is.na(pars$param) & is.na(pars$value)
  varying <- free & !is.na(pars$indvarying) & pars$indvarying
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
.ctPopRegressionSpec <- function(pars, poprank, explicit = TRUE) {
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
  # Above `nmean` the restriction stops being free: it starts fixing residual
  # variances the data does determine. Said once, here, rather than left for a
  # user to infer from a likelihood that moved.
  approximate <- rank < nmean
  list(rank = rank, basis = basis, regressed = regressed, roles = roles,
    nmean = nmean, approximate = approximate,
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
    betas <- paste0('beta_', p, '_', spec$basis)
    for (b in betas) addpar(b)
    coefficients[[length(coefficients) + 1L]] <- data.frame(
      param = p, basis = spec$basis, coefficient = betas,
      state = as.integer(stateof), row.names = NULL, stringsAsFactors = FALSE)
    predictor <- paste0('(', p, ' + ',
      paste0(betas, ' * state[', stateof, ']', collapse = ' + '), ')')
    cells <- which(!is.na(m$pars$param) & m$pars$param %in% p &
        !(m$pars$matrix %in% 'PARS' & m$pars$param %in% p &
            m$pars$transform %in% 'param'))
    # The rows just added carry the same label, so exclude anything added here.
    cells <- cells[cells <= nrow(m$pars)]
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
  m$pars <- rbind(m$pars, do.call(rbind, newpars))
  m$pars[] <- lapply(m$pars, utils::type.convert, as.is = TRUE)
  spec$coefficients <- do.call(rbind, coefficients)
  spec$cells <- do.call(rbind, drivencells)
  spec$state <- stateof
  m$popregression <- spec
  m
}
