# Stationary elements: a T0MEANS cell or a T0VAR diagonal cell named
# 'stationary' starts that latent at the mean or the variance it settles into,
# leaving the other initial conditions as written -- its covariances with the
# other states included, which T0VAR's construction keeps a valid covariance
# whatever the variance is.
#
# For a latent no other latent drives (its DRIFT row zero off the diagonal),
# both have closed forms: the mean is -CINT / DRIFT and the standard
# deviation sd / sqrt(-2 DRIFT), sd its DIFFUSION diagonal's (a DIFFUSION
# correlation with other states does not move the latent's own variance). So
# the cells are rewritten as those expressions, before anything else reads
# the model: the latent's free DRIFT, DIFFUSION and CINT parameters are
# declared in PARS, which the expressions and the cells they came from then
# read, and everything after is the calc machinery every model uses,
# derivatives included. A latent other latents drive has its stationary
# moments in the joint system's, which is stationary = TRUE's.
.ctStationaryElements <- function(ctm) {
  p <- ctm$pars
  mark <- !is.na(p$param) & p$param == 'stationary'
  if (!any(mark)) return(ctm)
  bad <- mark & !(p$matrix %in% 'T0MEANS' | (p$matrix %in% 'T0VAR' & p$row == p$col))
  if (any(bad)) stop("Only T0MEANS cells and the diagonal of T0VAR can be 'stationary' (",
    paste0(p$matrix[bad], '[', p$row[bad], ',', p$col[bad], ']', collapse = ', '),
    "); stationary = TRUE starts every latent from its stationary distribution.",
    call. = FALSE)
  if (!isTRUE(as.logical(ctm$continuoustime))) stop("'stationary' initial conditions ",
    "are for continuous time models.", call. = FALSE)
  # Under the factor construction a DIFFUSION diagonal is not the latent's own
  # standard deviation, which the closed form reads.
  if (identical(.ctCovMatCode(ctm), 1L)) stop("'stationary' initial conditions read ",
    "DIFFUSION's diagonal as a standard deviation, which covmattransform = 'cholesky' ",
    "does not give.", call. = FALSE)
  nl <- ctm$n.latent
  latents <- sort(unique(p$row[mark]))
  at <- function(m, r, c) which(p$matrix == m & p$row == r & p$col == c)
  fixedvalue <- function(k) if (!length(k)) 0 else suppressWarnings(as.numeric(p$value[k]))
  effects <- grep('_effect$', names(p), value = TRUE)
  levelflags <- grep('^indvarying_', names(p), value = TRUE)
  readers <- c('state', 'tdpreds', ctm$latentNames, ctm$TDpredNames)
  reads <- function(text) any(vapply(readers, function(r) grepl(paste0('\\b', r, '\\b'), text),
    logical(1)))
  # The natural value of a cell, as an expression other cells can read: a
  # number, an expression already, or a free parameter declared in PARS (once)
  # with the cell reading it.
  value <- function(m, r, c, what) {
    k <- at(m, r, c)
    if (!length(k) || is.na(p$param[k])) return(format(fixedvalue(k), digits = 17))
    par <- p$param[k]
    if (reads(par)) stop("A 'stationary' initial condition for ", ctm$latentNames[r],
      " needs its ", what, " fixed within each subject, and ", m, "[", r, ",", c, "] reads ",
      "the latent states or a time dependent predictor. Estimate T0MEANS and T0VAR ",
      "instead.", call. = FALSE)
    if (!grepl('^[A-Za-z.][A-Za-z0-9._]*$', par)) return(paste0('(', par, ')'))
    if (!par %in% p$param[p$matrix %in% 'PARS']) {
      rows <- p$row[p$matrix %in% 'PARS']
      rec <- p[k, , drop = FALSE]
      rec$matrix <- 'PARS'; rec$row <- if (length(rows)) max(rows) + 1L else 1L; rec$col <- 1L
      p <<- rbind(p, rec)
    }
    # The cell now reads the declared parameter, whose transform, variation and
    # predictor effects are its own.
    p$transform[k] <<- NA_character_
    p$indvarying[k] <<- FALSE
    if ('sdscale' %in% names(p)) p$sdscale[k] <<- NA_real_
    for (e in effects) p[k, e] <<- 'FALSE'
    for (l in levelflags) p[k, l] <<- FALSE
    par
  }
  expression <- function(k, text) {
    p$param[k] <<- text; p$value[k] <<- NA_real_; p$transform[k] <<- NA_character_
    p$indvarying[k] <<- FALSE
    if ('sdscale' %in% names(p)) p$sdscale[k] <<- NA_real_
    for (e in effects) p[k, e] <<- 'FALSE'
    for (l in levelflags) p[k, l] <<- FALSE
  }
  for (i in latents) {
    others <- p[p$matrix %in% 'DRIFT' & p$row == i & p$col != i & p$col <= nl, , drop = FALSE]
    driven <- !is.na(others$param) | (!is.na(suppressWarnings(as.numeric(others$value))) &
      suppressWarnings(as.numeric(others$value)) != 0)
    if (any(driven)) stop("A 'stationary' initial condition is for a latent no other ",
      "latent drives, and ", ctm$latentNames[i], " is driven through DRIFT[", i, ",",
      paste(others$col[driven], collapse = ', '), "]: its stationary moments are the joint ",
      "system's, which stationary = TRUE uses.", call. = FALSE)
    drift <- value('DRIFT', i, i, 'persistence')
    km <- at('T0MEANS', i, 1L)
    if (length(km) && isTRUE(mark[km])) {
      cint <- value('CINT', i, 1L, 'intercept')
      if (identical(suppressWarnings(as.numeric(cint)), 0)) {
        p$param[km] <- NA_character_; p$value[km] <- 0; p$transform[km] <- NA_character_
        p$indvarying[km] <- FALSE
        for (e in effects) p[km, e] <- 'FALSE'
        for (l in levelflags) p[km, l] <- FALSE
      } else expression(km, paste0('-(', cint, ') / (', drift, ')'))
    }
    kv <- at('T0VAR', i, i)
    if (length(kv) && isTRUE(mark[kv])) {
      if (length(km) && isTRUE(as.logical(p$indvarying[km]))) stop("T0VAR[", i, ",", i,
        "] is 'stationary' while T0MEANS[", i, "] varies over subjects; a random ",
        "initial mean takes the initial variance its random effect implies.", call. = FALSE)
      sd <- value('DIFFUSION', i, i, 'innovation sd')
      expression(kv, paste0('(', sd, ') / sqrt(-2 * (', drift, '))'))
    }
  }
  ctm$pars <- p
  ctm
}
