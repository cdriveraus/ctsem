#' Raise the order of a ctsem model's latent processes.
#'
#' Each latent in \code{indices} gains a derivative latent (named \code{d}
#' and the latent's name), which drives the level with a fixed unit
#' coefficient. That makes each a damped oscillator: a second order process,
#' whose autocorrelation can dip below zero and rise again, which no first
#' order process can.
#'
#' A model from \code{ctModel(type = 'ct')} is rewritten cell by cell on its
#' parameter table, with a new inertia time \eqn{\tau > 0} for each latent
#' (\code{inertia_} and the latent's name, in \code{PARS}):
#' \deqn{\tau \ddot{x} + \dot{x} = a x + q \, dW,}
#' where \eqn{a} is the level's own drift and \eqn{q} its diffusion, both
#' keeping their names and every setting -- transform, individual variation,
#' predictor effects. The derivative's row of \code{DRIFT} reads \eqn{a/\tau}
#' and \eqn{-1/\tau}, and its diffusion \eqn{q/\tau}. As \eqn{\tau} goes to
#' zero the derivative follows the level at once and the process is the
#' first order one it was, so the new model contains the old one with every
#' old parameter at its value; from \eqn{\tau = 1/(4|a|)} on, the process
#' oscillates. The derivative starts where the process alone holds it at
#' equilibrium: mean zero, uncorrelated with every state, and variance
#' \eqn{q^2/(2\tau)}, so the inertia time is the one new parameter. With the
#' noise left on the level (\code{diffusion = FALSE}), or none, its initial
#' mean, variance and covariances are free.
#'
#' A model in the older matrix-list form from \code{ctModel(type = 'omx')} is
#' written instead as the derivative pulled back by the level through an
#' estimated coefficient (negative unless \code{explosive}), with the level's
#' own autoregression as the derivative's damping. That is the same family of
#' processes, but its first order member is reached only as the damping, the
#' coefficient and the diffusion all grow without bound.
#'
#' @param ctm ctModel, \code{type = 'ct'} or \code{type = 'omx'}.
#' @param indices Vector of integers, which latents to raise the order of.
#' @param diffusion Move the process noise from the level to its derivative?
#' @param crosseffects Move the effects of the other processes on this one to
#'   its derivative (divided by \eqn{\tau}, for a \code{type = 'ct'} model)?
#'   For a \code{type = 'omx'} model, its effects on them move as well.
#' @param cint Move the continuous intercept to the derivative (divided by
#'   \eqn{\tau}, for a \code{type = 'ct'} model)?
#' @param explosive For a \code{type = 'omx'} model, allow explosive (non
#'   equilibrium returning) processes. A \code{type = 'ct'} model's level keeps
#'   its own drift's transform, which decides that.
#'
#' @return extended ctModel
#' @export
#'
#' @examples
#' om <- ctModel(LAMBDA=diag(1,2),DRIFT=0, 
#'   MANIFESTMEANS=0,type='omx',Tpoints=4)
#'   
#' om <- ctModelHigherOrder(om,1:2)
#' print(om$DRIFT)
#' 
#' print(om$pars)
ctModelHigherOrder <- function(ctm, indices,diffusion=TRUE, crosseffects=FALSE,cint=FALSE, explosive=FALSE){
  if (inherits(ctm, 'ctStanModel')) return(.ctHigherOrderPars(ctm, indices,
    diffusion = diffusion, crosseffects = crosseffects, cint = cint,
    explosive = explosive))
  ctm$latentNames <- c(ctm$latentNames,paste0('d',ctm$latentNames[indices]))
  nl <- ctm$n.latent
  
  for(i in 1:length(indices)){
    for(m in c('DRIFT','DIFFUSION','T0VAR')){
    ctm[[m]] <- rbind(cbind(ctm[[m]],0),0)
    }
    for(m in c('CINT','T0MEANS')){
      ctm[[m]] <- rbind(ctm[[m]],0)
    }
    ctm$LAMBDA <- cbind(ctm$LAMBDA,0)
    if(ctm$n.TDpred > 0) ctm$TDPREDEFFECT <- rbind(ctm$TDPREDEFFECT,0)
  }
  
  for(i in 1:length(indices)){
    # browser()
    if(!crosseffects) ctm$DRIFT[i+nl,i+nl] <- ctm$DRIFT[indices[i],indices[i]] #move ar param / value to higher order
    if(crosseffects){
      m <- 'DRIFT'
      # browser()
      ctm[[m]][i+nl,] <- ctm[[m]][indices[i],] #move crosseffect params / value to higher order
      ctm[[m]][,i+nl] <- ctm[[m]][,indices[i]] #move crosseffect params / value to higher order
      ctm[[m]][indices[i],] <- 0 #set order 1 crosseffect to 0
      ctm[[m]][,indices[i]] <- 0 #set order 1 crosseffect to 0
    }
    ctm$DRIFT[indices[i],indices[i]] <- 0 #set ar param / value to 0
    ctm$DRIFT[indices[i],i+nl] <- 1 #set effect of higher order to 1
    ctm$DRIFT[i+nl,indices[i]] <-paste0('drift_',ctm$latentNames[i+nl],'_',
      ctm$latentNames[indices[i]],
      ifelse(!explosive,'|-log1p(exp(-param*2))-1e-6','|param*2-1')) #estimate effect of 1st order on 2nd
    if(indices[i]==tail(indices,1)) rownames(ctm$DRIFT) <- ctm$latentNames
    if(indices[i]==tail(indices,1)) colnames(ctm$DRIFT) <- ctm$latentNames
    
    for(m in c('T0VAR','DIFFUSION')){
      if(diffusion && m=='DIFFUSION'){
        ctm[[m]][i+nl,] <- ctm[[m]][indices[i],] #move diffusion params / value to higher order
        ctm[[m]][,i+nl] <- ctm[[m]][,indices[i]] #move diffusion params / value to higher order
        ctm[[m]][indices[i],] <- 0 #set order 1 diffusion to 0
        ctm[[m]][,indices[i]] <- 0 #set order 1 diffusion to 0
      }
      if(indices[i]==tail(indices,1)) colnames(ctm[[m]]) <- ctm$latentNames
      if(indices[i]==tail(indices,1)) rownames(ctm[[m]]) <- ctm$latentNames
    }
    
    for(m in c('CINT','T0MEANS')){
      # browser()
      if(m %in% 'T0MEANS') ctm[[m]][i+nl,1] <- paste0('T0mean_d_',ctm$latentNames[indices[i]])
      if(!m %in% 'T0MEANS' && cint){
        ctm[[m]][i+nl,] <- ctm[[m]][indices[i],]
        ctm[[m]][indices[i],] <- 0 #set order 1 param to 0
        }
      if(indices[i]==tail(indices,1)) rownames(ctm[[m]]) <- ctm$latentNames
    }
  }
  # browser()
  for(i in c(nl+seq_along(indices))){
    for(j in c(1:nl,nl+seq_along(indices))){
      if(i >= j) ctm$T0VAR[i,j] <- paste0('T0var_',ctm$latentNames[i],'_',
        ctm$latentNames[j])
    }
  }
  ctm$n.latent <- ctm$n.latent + length(indices)
  
  return(ctm)
}


# The `type = 'ct'` form, `tau x'' + x' = a x + q dW`, written on the parameter
# table a cell at a time, so the level's drift and diffusion keep their names
# and every setting -- transform, individual variation, sdscale, predictor
# effects -- moved into PARS, where the derivative's row reads them divided by
# the inertia time. It is the damped oscillator the matrix-list form writes,
# reparameterised so that the first order model lies at tau -> 0 along one
# coordinate with every other parameter where it was: written as the
# derivative's damping and a restoring coefficient, the first order model is
# reached only as both run off to infinity, together with the diffusion, and
# the refits ctModelImprove() made from a first order fit followed them for
# hours (raw -17000 on one simulated search) without stopping. A person-varying
# autoregression copied onto the damping also made the noise vary with the
# rate, so that limit was not the first order model either: 44 log likelihood
# units below it there. New free cells take ctsem's own defaults
# (`.ctModelDefaultFreePar`) and do not vary over subjects.
.ctHigherOrderPars <- function(ctm, indices, diffusion = TRUE, crosseffects = FALSE,
  cint = FALSE, explosive = FALSE) {
  if (!isTRUE(ctm$continuoustime)) stop("Raising the order is for continuous ",
    "time models.", call. = FALSE)
  if (isTRUE(explosive)) stop("explosive is for a model from ctModel(type = 'omx'). ",
    "In a type = 'ct' model the level's own drift keeps its transform, which decides ",
    "whether the process may be explosive.", call. = FALSE)
  nl <- ctm$n.latent
  indices <- as.integer(indices)
  if (anyDuplicated(indices) || any(indices < 1L | indices > nl)) stop(
    "indices must be distinct latents, from 1 to ", nl, ".", call. = FALSE)
  names0 <- ctm$latentNames
  ed <- .ctParsLatentEditor(ctm, paste0('d', names0[indices]))
  allnames <- ed$names
  # The lower-triangle cell between latents a and b, where a covariance's
  # parameter lives.
  low <- function(a, b) c(max(a, b), min(a, b))
  # A cell's term divided by the inertia time, or zero where the cell is.
  over <- function(m, r, c, tr, tc, tau) {
    term <- .ctHigherOrderTerm(ed, m, r, c)
    if (is.null(term)) ed$fixed(m, tr, tc, 0) else
      ed$expression(m, tr, tc, sprintf('%s / %s', term, tau))
  }
  covaries <- function(a, b) {
    rec <- ed$cell('T0VAR', low(a, b)[1], low(a, b)[2])
    !is.null(rec) && (!is.na(rec$param) || isTRUE(suppressWarnings(as.numeric(rec$value)) != 0))
  }
  for (k in seq_along(indices)) {
    i <- indices[k]; n <- nl + k
    tau <- paste0('inertia_', names0[i])
    ed$free('PARS', ed$npars() + 1L, 1L, tau, transform = 'log1p_exp(2 * param)',
      tipred = FALSE)
    over('DRIFT', i, i, n, i, tau)
    ed$expression('DRIFT', n, n, sprintf('-1 / %s', tau))
    ed$fixed('DRIFT', i, i, 0)
    ed$fixed('DRIFT', i, n, 1)
    # The other processes' effects on this one, through its derivative; its
    # effects on them stay on the level, where the first order model has them.
    if (crosseffects) for (j in seq_len(nl)) if (j != i) {
      over('DRIFT', i, j, n, j, tau)
      ed$fixed('DRIFT', i, j, 0)
    }
    q <- if (diffusion) .ctHigherOrderTerm(ed, 'DIFFUSION', i, i)
    if (diffusion) {
      over('DIFFUSION', i, i, n, n, tau)
      # The correlations of the noise with the other processes' move with it.
      for (j in seq_len(nl)) if (j != i) {
        src <- low(i, j)
        ed$copy('DIFFUSION', src[1], src[2], n, j)
      }
      for (j in seq_len(nl)) { ed$fixed('DIFFUSION', i, j, 0); ed$fixed('DIFFUSION', j, i, 0) }
    }
    if (cint) { over('CINT', i, 1, n, 1, tau); ed$fixed('CINT', i, 1, 0) }
    # The derivative starts where the process alone holds it at equilibrium:
    # mean zero, uncorrelated with every state (a stationary process is with
    # its own derivative), and the variance q^2 / (2 tau) that its noise
    # q / tau builds up at rate 1 / tau. Free, its initial mean, variance and
    # covariances fitted starting values rather than dynamics -- 1.7 log
    # likelihood units on simulated first order data, for four parameters --
    # and lost their influence as tau -> 0, four flat directions beside the
    # inertia's at the first order model. Not 'stationary', which solves for
    # the moments of everything driving the level: a daily cycle never
    # settles, so has none, and the Fisher search's PA came out NaN.
    if (!is.null(q)) {
      ed$fixed('T0MEANS', n, 1, 0)
      ed$expression('T0VAR', n, n, sprintf('%s / sqrt(2 * %s)', q, tau))
    } else {
      # Free, covarying with its level, the derivatives raised beside it, and
      # whatever the level covaries with; free with anything else -- a
      # cycle's clock -- the Fisher search's PA derivative ran its covariance
      # with the clock's sine out to a raw 8.
      ed$free('T0MEANS', n, 1, paste0('T0mean_', allnames[n]))
      lev <- c(seq_len(nl), indices)
      for (j in seq_len(n)) if (j == n || j == i || covaries(i, lev[j])) ed$free('T0VAR', n,
        j, if (j == n) paste0('T0var_', allnames[n]) else paste0('T0var_', allnames[n], '_',
          allnames[j]))
    }
  }
  ed$finish()
}

# What a cell contributes to an expression: its parameter's name, declared in
# PARS so expressions can read it (the cell is then the caller's to rewrite);
# its value or expression, bracketed; NULL where it is zero.
.ctHigherOrderTerm <- function(ed, m, r, c) {
  rec <- ed$cell(m, r, c)
  if (is.null(rec)) return(NULL)
  if (is.na(rec$param) || !nzchar(rec$param)) {
    v <- suppressWarnings(as.numeric(rec$value))
    if (is.na(v) || v == 0) return(NULL)
    return(sprintf('(%s)', format(v, digits = 17)))
  }
  if (grepl('^[A-Za-z][A-Za-z0-9_.]*$', rec$param) && rec$param != 'stationary') {
    ed$declare(m, r, c)
    return(rec$param)
  }
  sprintf('(%s)', rec$param)
}

# A model's parameter table opened to new latents (`newnames`, appended after
# the existing ones): every cell of their rows and columns set to zero, then
# edited a cell at a time by the closures returned -- `fixed`; `free`, with
# ctsem's default transform unless one is given, varying over no one, and
# predictor effects as `tipred` says; `copy`, a cell with every setting it has;
# `expression`, a cell computed from the names it reads; `declare`, a free
# cell's parameter declared in PARS so expressions can read it; `cell`, a
# cell's current row -- and `finish()` returns the model. Shared by ctModelHigherOrder(),
# ctModelSlowComponent(), ctModelPeriodicComponent() and ctModelUniqueProcess(),
# so a cell written by any is written the same way.
.ctParsLatentEditor <- function(ctm, newnames) {
  pars <- ctm$pars
  nl <- ctm$n.latent
  allnames <- c(ctm$latentNames, newnames)
  N <- length(allnames)
  effects <- grep('_effect$', names(pars), value = TRUE)
  levelflags <- grep('^indvarying_', names(pars), value = TRUE)
  template <- pars[1, , drop = FALSE]
  at <- function(m, r, c) which(pars$matrix == m & pars$row == r & pars$col == c)
  put <- function(m, r, c, rec) {
    rec$matrix <- m; rec$row <- as.integer(r); rec$col <- as.integer(c)
    k <- at(m, r, c)
    if (length(k)) pars[k, ] <<- rec[, names(pars)] else pars <<- rbind(pars, rec[, names(pars)])
  }
  fixed <- function(m, r, c, v) {
    rec <- template
    rec$param <- NA_character_; rec$value <- v; rec$transform <- NA_character_
    rec$indvarying <- FALSE
    if ('sdscale' %in% names(rec)) rec$sdscale <- NA_real_
    for (e in effects) rec[[e]] <- 'FALSE'
    for (l in levelflags) rec[[l]] <- FALSE
    put(m, r, c, rec)
  }
  # A new free cell takes the model's default predictor effects, but only of
  # the user's own predictors: one a component or the search added states its
  # effects cell by cell (`statedTIpreds`) -- a cycle's clock phases set its
  # initial means and nothing else. Defaulting those too gave the Fisher
  # search's second-order processes 38 effects of the clock on their drifts
  # and initial conditions, and a refit that crawled for hours.
  stated <- paste0(ctm$statedTIpreds, '_effect')
  # A plain name: a parameter, as opposed to an expression.
  plain <- function(x) !is.na(x) & grepl('^[A-Za-z.][A-Za-z0-9._]*$', x)
  # PARS rows that hold a parameter moved out of a matrix cell sit straight
  # after that cell (.ctInsertMovedRows does this for terms the score tests
  # and the search write into CINT): random effects are ordered as their
  # parameters are declared (.ctIndvaryingIndex), so a moved parameter keeps
  # its place among them. Those already so placed, and those `declare` adds,
  # stay so through `finish`; a PARS declaration the user wrote stays where
  # it was.
  attached <- {
    k <- which(pars$matrix %in% 'PARS' & seq_len(nrow(pars)) > 1L)
    k <- k[!pars$matrix[k - 1L] %in% 'PARS' & plain(pars$param[k - 1L]) &
      pars$param[k - 1L] == pars$param[k]]
    unique(pars$param[k])
  }
  declared <- character()
  free <- function(m, r, c, name, transform = NULL, tipred = isTRUE(ctm$tipredDefault)) {
    rec <- template
    def <- .ctModelDefaultFreePar(m, r, c, TRUE)
    rec$param <- name; rec$value <- NA_real_
    rec$transform <- if (is.null(transform)) def$transform else transform
    rec$indvarying <- FALSE
    if ('sdscale' %in% names(rec)) rec$sdscale <- 1
    for (e in effects) rec[[e]] <- as.character(isTRUE(tipred) && !e %in% stated)
    for (l in levelflags) rec[[l]] <- FALSE
    put(m, r, c, rec)
  }
  # A cell's settings put in another cell. Its parameter keeps the source's
  # place in the table (`copied`), so a random effect moved to another cell
  # keeps its place among the others and a fit's population loadings mean
  # what they did: moved by column, the Fisher search's derivative took
  # PA's innovation correlation past four effects, and its refit started on
  # a fifth loading dimension the current fit had not got.
  copied <- list()
  copy <- function(m, fr, fc, tr, tc) {
    k <- at(m, fr, fc)
    if (!length(k)) return(fixed(m, tr, tc, 0))
    rec <- pars[k, , drop = FALSE]
    if (plain(rec$param) && is.null(copied[[rec$param]]))
      copied[[rec$param]] <<- list(m = m, r = fr, c = fc)
    put(m, tr, tc, rec)
  }
  # A cell computed from the names it reads, which ctModel takes as
  # references: states, predictors, PARS parameters.
  expression <- function(m, r, c, expr) {
    rec <- template
    rec$param <- expr; rec$value <- NA_real_; rec$transform <- NA_character_
    rec$indvarying <- FALSE
    if ('sdscale' %in% names(rec)) rec$sdscale <- NA_real_
    for (e in effects) rec[[e]] <- 'FALSE'
    for (l in levelflags) rec[[l]] <- FALSE
    put(m, r, c, rec)
  }
  cell <- function(m, r, c) {
    k <- at(m, r, c)
    if (length(k) == 1L) pars[k, , drop = FALSE] else NULL
  }
  # The parameter of a free cell, declared in PARS with all its settings --
  # transform, individual variation, predictor effects -- so an expression
  # elsewhere can read it; the cell itself then refers to that declaration.
  # Nothing if PARS declares it already. Where the caller then rewrites the
  # cell, the declaration keeps that cell's place (`origin`).
  origin <- list()
  declare <- function(m, r, c) {
    rec <- cell(m, r, c)
    if (is.null(rec) || is.na(rec$param)) return(invisible(FALSE))
    if (!any(pars$matrix == 'PARS' & pars$param %in% rec$param)) {
      put('PARS', sum(pars$matrix == 'PARS') + 1L, 1L, rec)
      declared <<- c(declared, rec$param)
      origin[[rec$param]] <<- list(m = m, r = r, c = c)
    }
    invisible(TRUE)
  }
  npars <- function() sum(pars$matrix == 'PARS')
  # An initial condition at its stationary value (`.ctStationaryElements`).
  stationary <- function(m, r, c) {
    rec <- template
    rec$param <- 'stationary'; rec$value <- NA_real_; rec$transform <- NA_character_
    rec$indvarying <- FALSE
    if ('sdscale' %in% names(rec)) rec$sdscale <- NA_real_
    for (e in effects) rec[[e]] <- 'FALSE'
    for (l in levelflags) rec[[l]] <- FALSE
    put(m, r, c, rec)
  }
  for (n in seq_len(N)[-seq_len(nl)]) {
    for (j in seq_len(N)) {
      for (m in c('DRIFT', 'DIFFUSION', 'T0VAR')) { fixed(m, n, j, 0); fixed(m, j, n, 0) }
    }
    fixed('CINT', n, 1, 0)
    for (j in seq_len(ctm$n.manifest)) fixed('LAMBDA', j, n, 0)
    for (j in seq_len(ctm$n.TDpred)) fixed('TDPREDEFFECT', n, j, 0)
  }
  finish <- function() {
    mo <- match(pars$matrix, unique(c(ctm$pars$matrix, 'PARS')))
    rank <- numeric(nrow(pars))
    rank[order(mo, pars$col, pars$row)] <- seq_len(nrow(pars))
    cellat <- function(o) which(pars$matrix == o$m & pars$row == o$r & pars$col == o$c)
    for (k in which(!pars$matrix %in% 'PARS' & plain(pars$param) &
      pars$param %in% names(copied))) {
      at <- cellat(copied[[pars$param[k]]])
      if (length(at) && !identical(pars$param[at], pars$param[k])) rank[k] <- rank[at] + 0.5
    }
    moved <- which(pars$matrix %in% 'PARS' & pars$param %in% c(attached, declared))
    for (k in moved) {
      at <- which(!pars$matrix %in% 'PARS' & plain(pars$param) & pars$param == pars$param[k])
      o <- origin[[pars$param[k]]]
      if (!length(at) && !is.null(o)) at <- cellat(o)
      if (length(at)) rank[k] <- min(rank[at]) + 0.5
    }
    out <- pars[order(rank), ]
    rownames(out) <- NULL
    ctm$pars <- out
    ctm$latentNames <- allnames
    ctm$n.latent <- N
    ctm
  }
  list(fixed = fixed, free = free, copy = copy, stationary = stationary, finish = finish,
    expression = expression, cell = cell, declare = declare, npars = npars,
    names = allnames)
}
