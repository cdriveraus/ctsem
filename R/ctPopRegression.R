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
# Which RAWPOPVAR cells a user has stated for a regressed effect.
#
# RAWPOPVAR is the specification surface for the population covariance
# (R/ctModelRawPopVar.R), and `.ctJuliaAugmentRandomEffects()` reads it per varying
# parameter: a number there fixes a cell and a label estimates it.
#
# Worth knowing what those numbers mean, because it is not what the surface
# says. The diagonal is a population sd on the parameter's natural scale, and
# that is exact for a linear transform -- measured, 0.3 in gives 0.3 out. The
# off-diagonal is documented and error-checked as a correlation but is really
# the coordinate `constraincorsqrt1()` consumes, so 0.3 in gives a population
# correlation of 0.54 and 0.5 gives 0.79. Zero is the exception and is exact:
# uncorrelated in means uncorrelated out. None of that is this feature's doing,
# but a guard here that reasoned about "the correlation the user asked for"
# would be reasoning about the wrong number. A regressed effect has neither -- its spread and its correlations follow
# from the basis -- so anything stated about it would be **silently dropped**,
# which is the one outcome this whole feature exists to avoid. Both a fixed
# value and a relabelling count: a label differing from the default is an
# equality constraint, and that is a specification too.
# A RAWPOPVAR entry about a regressed effect cannot be honoured: its spread
# follows entirely from the basis, so there is no cell to fix.
#
# The basis effects have a separate and narrower problem under the factor
# construction, checked in `.ctPopRegressionFactorConflicts()` below rather
# than here -- the two are different claims and refusing them together refused
# statements that are perfectly honourable.
.ctPopRegressionRawPopVarConflicts <- function(model, regressed) {
  popcov <- model[['RAWPOPVAR']]
  if (is.null(popcov) || !length(popcov) || !length(regressed)) return(character())
  default <- .ctModelRawPopVar(model$pars)
  names <- rownames(popcov)
  regressed <- intersect(regressed, names)
  out <- character()
  for (r in regressed) for (other in names) {
    for (cell in unique(c(paste(r, other), paste(other, r)))) {
      coords <- strsplit(cell, ' ', fixed = TRUE)[[1]]
      stated <- .ctModelRawPopVarEntry(model, coords[1L], coords[2L])
      if (is.na(stated) || !nzchar(stated)) next
      expected <- if (!is.null(default) && all(coords %in% rownames(default)))
        as.character(default[coords[1L], coords[2L]]) else NA_character_
      if (!is.na(expected) && identical(stated, expected)) next
      # An upper-triangle zero is RAWPOPVAR's own placeholder, not a statement.
      if (identical(stated, '0') && !is.na(expected) && identical(expected, '0')) next
      out <- c(out, sprintf("RAWPOPVAR['%s', '%s'] = %s", coords[1L], coords[2L],
        stated))
    }
  }
  unique(out)
}

# What a factor construction cannot honour about the basis effects.
#
# `Sigma = M M'` with M lower triangular, so `Sigma[i,j]` for `i > j` is
# `sum_{m<=j} M[i,m] M[j,m]` -- j products, not one. Nothing below the diagonal
# is a cell of the estimand: an off-diagonal is a factor entry whose implied
# correlation depends on the rest of its row, and an effect's spread is a row
# norm rather than a cell, except for the first basis effect, which loads on
# one dimension and no other, so `|M[1,1]|` is its spread exactly.
#
# A zero is no exception. Only against the first dimension does it reduce to
# one cell (`M[i,1] M[1,1] = 0`), for any later one it is a constraint across a
# row, and in neither case does the rewrite act on it -- a stated zero leaves
# every loading free. So it is refused with the rest rather than accepted and
# dropped. Fixing loadings from a declared zero pattern is a real feature and
# would need a rotation-rigidity check to go with it; until then this is what
# is true.
.ctPopRegressionFactorConflicts <- function(model, basis) {
  popcov <- model[['RAWPOPVAR']]
  if (is.null(popcov) || !length(popcov) || length(basis) < 2L) return(character())
  names <- rownames(popcov)
  # RAWPOPVAR is constructed with a default label in every cell, so "is there
  # text here" is not the question -- a default label is the absence of a
  # statement, and treating it as one refuses every model that has not been
  # touched. Compared against the default the same way
  # `.ctPopRegressionRawPopVarConflicts()` does.
  default <- .ctModelRawPopVar(model$pars)
  isdefault <- function(stated, ...) {
    coords <- c(...)
    if (is.null(default) || !all(coords %in% rownames(default))) return(FALSE)
    identical(stated, as.character(default[coords[1L], coords[length(coords)]]))
  }
  out <- character()
  for (b in intersect(basis[-1L], names)) {
    diagonal <- .ctModelRawPopVarEntry(model, b)
    if (!is.na(diagonal) && nzchar(diagonal) && !isdefault(diagonal, b, b)) {
      out <- c(out, sprintf("RAWPOPVAR['%s', '%s'] = %s (a loading, not a standard deviation)",
        b, b, diagonal))
    }
    for (other in setdiff(names, b)) {
      stated <- .ctModelRawPopVarEntry(model, b, other)
      if (is.na(stated) || !nzchar(stated)) next
      if (isdefault(stated, b, other)) next
      out <- c(out, sprintf("RAWPOPVAR['%s', '%s'] = %s (a factor entry, not a correlation)",
        b, other, stated))
    }
  }
  unique(out)
}

.ctPopRegressionSpec <- function(pars, poprank, explicit = TRUE, model = NULL,
    augmented = TRUE) {
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
  # Laplace has no carrier states and rewrites the parameter table instead, so
  # it is not restricted here.
  # Only where the reduction reaches a cell: with nothing regressed the rewrite
  # returns the model untouched, so a rank equal to the number of varying
  # parameters is a no-op and there is nothing to refuse.
  if (isTRUE(augmented) && length(regressed)) {
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

  # A RAWPOPVAR statement about a regressed effect cannot be honoured, so it is
  # refused rather than ignored. Asked for, that is an error; defaulted, the
  # user's own specification is the more explicit statement of the two and wins.
  if (!is.null(model)) {
    conflicts <- c(.ctPopRegressionRawPopVarConflicts(model, regressed),
      .ctPopRegressionFactorConflicts(model, basis))
    if (length(conflicts)) {
      if (!isTRUE(explicit)) return(NULL)
      # Opens with the phrase it always opened with. The rest is new -- the
      # factor construction gives a second, narrower reason -- but the leading
      # clause is accurate for both and is what the tests and anything else
      # reading this message match on.
      stop('poprank would drop what RAWPOPVAR states: ',
        paste(conflicts, collapse = '; '),
        '. Under a reduced rank the population covariance is a factor. A ',
        'regressed effect has no spread of its own at all, and for a basis ',
        "effect past the first the spread is a row norm rather than a cell -- ",
        'only ', basis[1L], ' keeps a standard deviation this can state, and a ',
        'zero covariance is a constraint across a row of loadings rather than ',
        'a cell. Use poprank=NA to estimate the full covariance.',
        call. = FALSE)
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
  # The "cannot see its own spread" clause is true of the augmented filter and
  # false of everything else, so it is keyed on the route rather than always
  # said. Under laplace those coordinates are identified -- that route is the
  # one that identifies them -- and telling a laplace user their spreads were
  # invisible would be exactly backwards.
  augmented <- !identical(spec$route, 'parameters')
  if (length(variancecell) && augmented) {
    out <- paste0(out, ' ', paste(variancecell, collapse = ', '),
      if (length(variancecell) > 1) ' vary' else ' varies',
      ' only in DIFFUSION / MANIFESTVAR, where the augmented filter cannot see',
      if (length(variancecell) > 1) ' their' else ' its', ' own spread, so ',
      if (length(variancecell) > 1) 'they are' else 'it is',
      ' estimated as a regression on ', basis, '.')
  } else if (length(variancecell)) {
    out <- paste0(out, ' ', paste(variancecell, collapse = ', '),
      if (length(variancecell) > 1) ' are' else ' is',
      ' estimated as a regression on ', basis, ', with no variation',
      ' independent of ', basis, '.')
  }
  if (length(demoted)) {
    out <- paste0(out, ' ', paste(demoted, collapse = ', '),
      if (length(demoted) > 1) ' are' else ' is',
      ' also regressed on ', basis, ' to meet the requested rank, which is ',
      'below the ', spec$nmean, ' this model identifies -- an approximation, ',
      'and the retained parameters absorb what it drops.')
  }
  # And the closing sentence, for the same reason. On the augmented route the
  # restriction is free and laplace is the recommendation; on laplace the
  # restriction is the approximation and there is nothing better to point at.
  if (!augmented) {
    return(paste0(out, ' On this route those coordinates are identified, so',
      ' this is an approximation for parsimony or speed rather than a repair.',
      ' poprank=NA estimates the full covariance.'))
  }
  paste0(out, " poprank=NA estimates the full covariance instead;",
    " intoverpop='laplace' identifies it.")
}

# Make a random effect referenceable from another cell, for the routes that do
# not augment the state.
#
# On the augmented route a basis effect already has a carrier state, so a
# regressed cell can say `state[j]`. Under `intoverpop='laplace'` (and `'none'`,
# which prepares the same structure without integrating) there are no carrier
# states -- the random effect *is* the parameter -- so the only way one cell can
# refer to another's parameter is ctsem's own mechanism: the parameter lives in
# PARS, and `ctModelStatesAndPARS()` turns references to its label into
# `PARS[r,c]`.
#
# So the basis effects are moved into PARS: a row holding the raw parameter with
# the identity transform, carrying the `indvarying` flag and the `sdscale` that
# sets its population prior, and the cell it came from becomes that transform
# applied to the reference. The laplace spec builds `re_index` from the
# parameter table's indvarying entries, so a PARS row is eligible exactly as the
# original cell was, and nothing in the engine needs to know.
.ctPopRegressionToPars <- function(m, labels) {
  if (!length(labels)) return(m)
  effectcols <- grep('_effect$', names(m$pars), value = TRUE)
  parsrow <- suppressWarnings(max(c(0L,
    as.integer(m$pars$row[m$pars$matrix %in% 'PARS']))))
  template <- m$pars[1L, , drop = FALSE]
  added <- list()
  for (label in labels) {
    own <- which(!is.na(m$pars$param) & m$pars$param %in% label)
    if (!length(own)) next
    # Already a PARS entry with nothing wrapped round it: referenceable as is.
    if (all(m$pars$matrix[own] %in% 'PARS')) next
    transform <- as.character(m$pars$transform[own[1L]])
    if (is.na(transform) || !nzchar(transform)) transform <- 'param'
    sdscale <- m$pars$sdscale[own[1L]]
    parsrow <- parsrow + 1L
    row <- template
    row$matrix <- 'PARS'; row$row <- parsrow; row$col <- 1L
    row$param <- label; row$value <- NA; row$transform <- 'param'
    row$indvarying <- TRUE
    row$sdscale <- sdscale
    if (length(effectcols)) {
      row[, effectcols] <- FALSE
      # TI effects follow the parameter, which is now this row.
      for (cc in effectcols) row[[cc]] <- any(m$pars[own, cc] %in% TRUE)
    }
    for (ri in own) {
      m$pars$param[ri] <- gsub('(?<![[:alnum:]_.])param(?![[:alnum:]_.])',
        label, as.character(m$pars$transform[ri]), perl = TRUE)
      m$pars$transform[ri] <- NA
      m$pars$indvarying[ri] <- FALSE
      if (length(effectcols)) m$pars[ri, effectcols] <- FALSE
    }
    added[[length(added) + 1L]] <- row
  }
  if (length(added)) {
    m$pars <- rbind(m$pars, do.call(rbind, added))
    m$pars[] <- lapply(m$pars, utils::type.convert, as.is = TRUE)
  }
  m
}

# The rewrite for a route with no carrier states. Same structure as
# `.ctPopRegressionRewrite()`, referencing the basis parameters by label -- the
# second `ctModelStatesAndPARS()` call in `ctFit()` turns those into `PARS[r,c]`,
# which is why this has to run before it, exactly as the augmented rewrite does.
.ctPopRegressionRewriteParameters <- function(m, spec) {
  if (!length(spec$regressed)) return(m)
  m <- .ctPopRegressionToPars(m, spec$basis)
  effectcols <- grep('_effect$', names(m$pars), value = TRUE)
  template <- m$pars[1L, , drop = FALSE]
  parsrow <- suppressWarnings(max(c(0L,
    as.integer(m$pars$row[m$pars$matrix %in% 'PARS']))))
  newpars <- list()
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
    effects <- NULL
    if (length(effectcols)) {
      own <- which(!is.na(m$pars$param) & m$pars$param %in% p)
      if (length(own)) effects <- vapply(effectcols,
        function(cc) any(m$pars[own, cc] %in% TRUE), logical(1L))
    }
    addpar(p, effects)
    betas <- paste0('beta_', p, '_', spec$basis)
    for (b in betas) addpar(b)
    coefficients[[length(coefficients) + 1L]] <- data.frame(
      param = p, basis = spec$basis, coefficient = betas,
      state = NA_integer_, row.names = NULL, stringsAsFactors = FALSE)
    predictor <- paste0('(', p, ' + ',
      paste0(betas, ' * ', spec$basis, collapse = ' + '), ')')
    cells <- which(!is.na(m$pars$param) & m$pars$param %in% p)
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
      m$pars$param[ri] <- gsub('(?<![[:alnum:]_.])param(?![[:alnum:]_.])',
        predictor, transform, perl = TRUE)
      m$pars$transform[ri] <- NA
      m$pars$indvarying[ri] <- FALSE
      if (length(effectcols)) m$pars[ri, effectcols] <- FALSE
    }
  }
  m$pars <- rbind(m$pars, do.call(rbind, newpars))
  m$pars[] <- lapply(m$pars, utils::type.convert, as.is = TRUE)
  spec$coefficients <- do.call(rbind, coefficients)
  spec$cells <- do.call(rbind, drivencells)
  spec$route <- 'parameters'
  m$popregression <- spec
  m
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
    own_loadings <- paste0('L_', p, '_', seq_len(pos))
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
