# Stationary elements: a T0MEANS cell, or T0VAR cells, named 'stationary'
# start at the value the engine computes for them -- the entries of asymCINT
# and asymDIFFUSION (`_ctsem_stationary_elements_mean!`, `_cov!`) -- and every
# other initial condition stays as written.
#
# A T0VAR diagonal cell marks one latent's variance, which takes the stationary
# one with that latent's T0VAR correlations kept. An off-diagonal cell joins
# two marked latents into one block, whose covariance then takes the
# stationary block whole; a block has to be marked throughout. The engine
# places both by one congruence, so the initial covariance stays one.
#
# The moments are solved over the closure of the marked latents: they and
# every latent that drives them, through DRIFT. Its dynamics are autonomous,
# so a marked latent driven by others is fine, and a random walk or a trait
# elsewhere in the model does not get in the way -- only one in the closure
# does.
#
# Here the marked cells become fixed placeholders the engine writes over --
# T0MEANS 0, a T0VAR block the identity -- and the model carries the marks for
# the specification (`.ctJuliaObjectiveInputs`).
.ctStationaryElements <- function(ctm) {
  p <- ctm$pars
  mark <- !is.na(p$param) & p$param == 'stationary'
  if (!any(mark)) return(ctm)
  where <- function(rows) paste0(p$matrix[rows], '[', p$row[rows], ',', p$col[rows], ']',
    collapse = ', ')
  bad <- mark & !(p$matrix %in% c('T0MEANS', 'T0VAR'))
  if (any(bad)) stop("Only T0MEANS and T0VAR cells can be 'stationary' (", where(bad), ").",
    call. = FALSE)
  if (!isTRUE(as.logical(ctm$continuoustime))) stop("'stationary' initial conditions ",
    "are for continuous time models.", call. = FALSE)
  nl <- ctm$n.latent
  names <- ctm$latentNames
  means <- logical(nl)
  means[p$row[mark & p$matrix %in% 'T0MEANS']] <- TRUE

  # T0VAR blocks: a diagonal mark is a latent, an off-diagonal one joins two.
  vmark <- mark & p$matrix %in% 'T0VAR'
  variance <- logical(nl)
  variance[p$row[vmark & p$row == p$col]] <- TRUE
  block <- ifelse(variance, seq_len(nl), 0L)
  off <- which(vmark & p$row != p$col)
  for (k in off) {
    r <- p$row[k]; c <- p$col[k]
    if (!variance[r] || !variance[c]) stop("T0VAR[", r, ",", c, "] is 'stationary', so ",
      "the covariance of ", names[r], " and ", names[c], " is the stationary one, which ",
      "needs both their variances 'stationary' too (T0VAR[", r, ",", r, "] and T0VAR[",
      c, ",", c, "]).", call. = FALSE)
    block[block == block[c]] <- block[r]
  }
  marked <- function(r, c) any(vmark & p$row == max(r, c) & p$col == min(r, c))
  for (id in unique(block[block > 0L])) {
    members <- which(block == id)
    if (length(members) < 2L) next
    pairs <- utils::combn(members, 2L)
    missing <- !apply(pairs, 2L, function(rc) marked(rc[1L], rc[2L]))
    if (any(missing)) stop("The 'stationary' T0VAR cells join ",
      paste(names[members], collapse = ', '), " into one block, whose covariance is ",
      "then the stationary one throughout, but T0VAR[", max(pairs[, missing][1:2]), ",",
      min(pairs[, missing][1:2]), "] is not 'stationary'. Name it too, or leave the ",
      "off-diagonal cells as parameters to keep those covariances free.", call. = FALSE)
  }
  block[block > 0L] <- match(block[block > 0L], unique(block[block > 0L]))

  # A random initial mean's variance is the population block's, and a
  # stationary one would be added to it.
  t0m <- p$matrix %in% 'T0MEANS' & !mark & .ctVaryingRows(ctm)
  varyingmean <- vapply(seq_len(nl), function(i) any(t0m & p$row == i), logical(1))
  clash <- variance & varyingmean
  if (any(clash)) stop("T0VAR is 'stationary' for ", paste(names[clash], collapse = ', '),
    ", whose T0MEANS varies over subjects; a random initial mean takes the initial ",
    "variance its random effect implies. Name the T0MEANS cell 'stationary' too, or set ",
    "its indvarying to FALSE.", call. = FALSE)

  # The closure: the marked latents and whatever drives them.
  closure <- .ctStationaryClosure(ctm, means | variance)
  static <- .ctNeverSettles(ctm, which(closure))
  if (length(static)) stop("A 'stationary' initial condition needs the latents it depends on ",
    "to settle somewhere, and DRIFT, as fixed, leaves ", paste(static, collapse = ', '),
    " without a stationary state. Estimate T0MEANS and T0VAR for the latents they ",
    "drive instead.", call. = FALSE)
  # No noise reaching a block leaves it no stationary variance.
  diffusion <- p[p$matrix %in% 'DIFFUSION' & p$row <= nl & p$col <= nl, , drop = FALSE]
  dvalue <- suppressWarnings(as.numeric(diffusion$value))
  noisy <- unique(diffusion$row[!is.na(diffusion$param) | (!is.na(dvalue) & dvalue != 0)])
  if (any(variance) && !any(closure[noisy])) stop("A 'stationary' T0VAR cell takes the ",
    "variance DIFFUSION builds up, and no DIFFUSION reaches ",
    paste(names[variance], collapse = ', '), ".", call. = FALSE)

  # The cells become placeholders the engine writes over.
  effects <- grep('_effect$', names(p), value = TRUE)
  levelflags <- grep('^indvarying_', names(p), value = TRUE)
  t0mark <- mark | (p$matrix %in% 'T0VAR' & p$row != p$col & block[pmin(p$row, nl)] > 0L &
    block[pmin(p$col, nl)] > 0L & block[pmin(p$row, nl)] == block[pmin(p$col, nl)] &
    p$row <= nl & p$col <= nl)
  p$param[t0mark] <- NA_character_
  p$value[t0mark] <- ifelse(p$matrix[t0mark] %in% 'T0VAR' & p$row[t0mark] == p$col[t0mark],
    1, 0)
  p$transform[t0mark] <- NA_character_
  p$indvarying[t0mark] <- FALSE
  if ('sdscale' %in% names(p)) p$sdscale[t0mark] <- NA_real_
  for (e in effects) p[t0mark, e] <- 'FALSE'
  for (l in levelflags) p[t0mark, l] <- FALSE
  ctm$pars <- p
  # `t0stationary` rather than a name starting `stationary`: `$` matches a
  # prefix, so that would answer `model$stationary` on a model without it.
  ctm$t0stationary <- list(means = means, block = as.integer(block),
    closure = closure)
  ctm
}

# The names of the latents among `which` that DRIFT's fixed cells keep from
# ever settling: a block of them whose DRIFT rows are fixed and read only each
# other, with an eigenvalue whose real part is not negative -- a trait held
# constant, a random walk, or a clock such as ctModelPeriodicComponent()'s,
# which has no stationary moments for anything it drives either, so the
# engine's solve for them came out NaN with nothing said. A free cell in a
# row could move the block, so only wholly fixed rows are judged.
.ctNeverSettles <- function(ctm, which = seq_len(ctm$n.latent)) {
  nl <- ctm$n.latent
  p <- ctm$pars
  drift <- p[p$matrix %in% 'DRIFT' & p$row <= nl & p$col <= nl, , drop = FALSE]
  value <- suppressWarnings(as.numeric(drift$value))
  fixed <- is.na(drift$param) & !is.na(value)
  A <- matrix(0, nl, nl)
  A[cbind(drift$row[fixed], drift$col[fixed])] <- value[fixed]
  rowfixed <- vapply(seq_len(nl), function(i) any(drift$row == i) &&
    all(fixed[drift$row == i]), logical(1))
  block <- intersect(which, which(rowfixed))
  repeat {
    inside <- vapply(block, function(i) all(A[i, setdiff(seq_len(nl), block)] == 0), logical(1))
    if (all(inside)) break
    block <- block[inside]
  }
  out <- character()
  while (length(block)) {
    # One connected part of the block at a time.
    part <- block[1L]
    repeat {
      grown <- block[colSums(A[part, block, drop = FALSE] != 0) > 0 |
        rowSums(A[block, part, drop = FALSE] != 0) > 0 | block %in% part]
      if (setequal(grown, part)) break
      part <- grown
    }
    ev <- eigen(A[part, part, drop = FALSE], only.values = TRUE)$values
    if (any(Re(ev) >= -1e-12)) out <- c(out, ctm$latentNames[sort(part)])
    block <- setdiff(block, part)
  }
  out
}

# The latents whose stationary moments a model's 'stationary' cells need, one
# logical per latent: the latents `start` names -- by default those a T0MEANS
# or T0VAR cell named 'stationary' belongs to -- and every latent that drives
# them through DRIFT, a cell that is free or fixed away from zero. Their
# DRIFT, CINT and DIFFUSION have to be fixed within each subject, so what
# extends a model's dynamics (the score tests' terms) leaves them out.
.ctStationaryClosure <- function(ctm, start = NULL) {
  p <- ctm$pars
  nl <- ctm$n.latent
  if (is.null(start)) {
    mark <- !is.na(p$param) & p$param == 'stationary' &
      p$matrix %in% c('T0MEANS', 'T0VAR') & p$row <= nl
    start <- logical(nl)
    start[unique(c(p$row[mark], p$col[mark & p$matrix %in% 'T0VAR']))] <- TRUE
  }
  drift <- p[p$matrix %in% 'DRIFT' & p$row <= nl & p$col <= nl, , drop = FALSE]
  value <- suppressWarnings(as.numeric(drift$value))
  drives <- drift$row != drift$col & (!is.na(drift$param) | (!is.na(value) & value != 0))
  closure <- start
  repeat {
    grown <- closure
    grown[drift$col[drives & drift$row %in% which(closure)]] <- TRUE
    if (identical(grown, closure)) break
    closure <- grown
  }
  closure
}
