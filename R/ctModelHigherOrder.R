#' Raise the order of a ctsem model's latent processes.
#'
#' Each latent in \code{indices} gains a derivative latent (named \code{d}
#' and the latent's name): the level is driven by its derivative with a fixed
#' unit coefficient, the derivative is pulled back by the level through an
#' estimated coefficient (negative unless \code{explosive}), and the level's
#' own autoregression becomes the derivative's damping. That makes each a
#' damped oscillator: a second order process, whose autocorrelation can dip
#' below zero and rise again, which no first order process can.
#'
#' Works on a model from \code{ctModel(type = 'ct')}, rewriting its parameter
#' table cell by cell so every other setting of each cell -- transform,
#' individual variation, predictor effects -- is kept, and on the older
#' matrix-list form from \code{ctModel(type = 'omx')}.
#'
#' @param ctm ctModel, \code{type = 'ct'} or \code{type = 'omx'}.
#' @param indices Vector of integers, which latents to raise the order of.
#' @param diffusion Shift the diffusion parameters / values to the higher order?
#' @param crosseffects Shift cross coupling parameters of the DRIFT matrix to the higher order?
#' @param cint shift continuous intercepts to higher order?
#' @param explosive Allow explosive (non equilibrium returning) processes?
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


# The `type = 'ct'` form: the same surgery on the parameter table, a cell at a
# time, so a moved cell takes its transform, individual variation, sdscale and
# predictor effects with it. New free cells take ctsem's own defaults
# (`.ctModelDefaultFreePar`) and do not vary over subjects.
.ctHigherOrderPars <- function(ctm, indices, diffusion = TRUE, crosseffects = FALSE,
  cint = FALSE, explosive = FALSE) {
  if (!isTRUE(ctm$continuoustime)) stop("Raising the order is for continuous ",
    "time models.", call. = FALSE)
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
  for (k in seq_along(indices)) {
    i <- indices[k]; n <- nl + k
    if (crosseffects) {
      for (j in seq_len(nl)) if (j != i) {
        ed$copy('DRIFT', i, j, n, j); ed$copy('DRIFT', j, i, j, n)
        ed$fixed('DRIFT', i, j, 0); ed$fixed('DRIFT', j, i, 0)
      }
    }
    ed$copy('DRIFT', i, i, n, n)
    ed$fixed('DRIFT', i, i, 0)
    ed$fixed('DRIFT', i, n, 1)
    ed$free('DRIFT', n, i, paste0('drift_', allnames[n], '_', names0[i]),
      transform = if (explosive) 'param * 2 - 1' else '-log1p_exp(-param * 2) - 1e-6')
    if (diffusion) {
      for (j in seq_len(nl)) {
        src <- low(i, j)
        if (j == i) ed$copy('DIFFUSION', i, i, n, n) else ed$copy('DIFFUSION', src[1], src[2], n, j)
      }
      for (j in seq_len(nl)) { ed$fixed('DIFFUSION', i, j, 0); ed$fixed('DIFFUSION', j, i, 0) }
    }
    if (cint) { ed$copy('CINT', i, 1, n, 1); ed$fixed('CINT', i, 1, 0) }
    ed$free('T0MEANS', n, 1, paste0('T0mean_', allnames[n]))
    for (j in seq_len(n)) ed$free('T0VAR', n, j, if (j == n) paste0('T0var_', allnames[n]) else
      paste0('T0var_', allnames[n], '_', allnames[j]))
  }
  ed$finish()
}

# A model's parameter table opened to new latents (`newnames`, appended after
# the existing ones): every cell of their rows and columns set to zero, then
# edited a cell at a time by the closures returned -- `fixed`; `free`, with
# ctsem's default transform unless one is given, varying over no one, and
# predictor effects as `tipred` says; `copy`, a cell with every setting it has
# -- and `finish()` returns the model. Shared by ctModelHigherOrder(),
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
  free <- function(m, r, c, name, transform = NULL, tipred = isTRUE(ctm$tipredDefault)) {
    rec <- template
    def <- .ctModelDefaultFreePar(m, r, c, TRUE)
    rec$param <- name; rec$value <- NA_real_
    rec$transform <- if (is.null(transform)) def$transform else transform
    rec$indvarying <- FALSE
    if ('sdscale' %in% names(rec)) rec$sdscale <- 1
    for (e in effects) rec[[e]] <- as.character(isTRUE(tipred))
    for (l in levelflags) rec[[l]] <- FALSE
    put(m, r, c, rec)
  }
  copy <- function(m, fr, fc, tr, tc) {
    k <- at(m, fr, fc)
    if (length(k)) put(m, tr, tc, pars[k, , drop = FALSE]) else fixed(m, tr, tc, 0)
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
    out <- pars[order(match(pars$matrix, unique(ctm$pars$matrix)), pars$col, pars$row), ]
    rownames(out) <- NULL
    ctm$pars <- out
    ctm$latentNames <- allnames
    ctm$n.latent <- N
    ctm
  }
  list(fixed = fixed, free = free, copy = copy, finish = finish, names = allnames)
}
