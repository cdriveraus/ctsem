#' Prior and posterior distributions of a fit's parameters
#'
#' Plots, for each free population parameter of a fit, its prior and its
#' posterior on the natural scale: the population mean, and for an individually
#' varying parameter the population standard deviation, at each level it varies
#' at. Works for fits from either backend. The prior is ctsem's
#' \code{normal(0,1)} raw prior carried through each parameter's transform --
#' the prior \code{priors = TRUE} fits with -- so on a fit made without it the
#' plot shows how far the data moved each parameter from where that prior
#' would put it. The posterior is the fit's draws: genuine posterior draws for
#' a sampled fit, the approximate ones \code{\link{ctFitUncertainty}} writes
#' for an optimized one. \code{ctStanPlotPost} is maintained as a
#' backward-compatible alias.
#'
#' @param obj A fit from \code{\link{ctFit}}, on either backend. A model on its
#' own has no posterior; \code{plot(model)} plots its priors.
#' @param rows Character string 'all', or which parameters to plot, by position
#' in the fit's free population parameters, in parameter order.
#' @param npp Integer number of parameters to show per page.
#' @param priorwidth if TRUE, plots will be scaled to show bulk of both the prior
#' and posterior distributions. If FALSE, scale is based only on the posterior.
#' @param smoothness Positive numeric -- multiplier to modify smoothness of density plots, higher is smoother but
#' can cause plots to exceed natural boundaries, such as standard deviations below zero.
#' @param priorsamples number of samples from prior to use. More is slower.
#' @param wait If true, user is prompted to continue before plotting next graph.
#' If false, graphs are plotted one after another without waiting.
#' @param plot Logical, if FALSE, ggplot objects are returned in a list instead of plotting.
#' @param ... For a stan fit, passed to the preparation of the prior draws:
#' \code{cores = x} will speed things up, where x is the number of cpu cores to
#' use. A julia fit takes nothing here.
#' @examples
#' \donttest{
#' ctPlotPosterior(ctstantestfit, rows=3:4)
#' }
#' @return With \code{plot = TRUE}, \code{NULL} invisibly, having drawn the
#'   plots; with \code{plot = FALSE}, a list of ggplot objects.
#' @aliases ctStanPlotPost
#' @export

ctPlotPosterior<-function(obj, rows='all', npp=4,priorwidth=TRUE,
  smoothness=1,priorsamples=10000,
  plot=TRUE,wait=FALSE,...){

  # inherits(), not class() %in%: a fit carries both 'ctStanFit' and 'ctFit',
  # so the membership test returns two values and the condition errors.
  if(inherits(obj, 'ctStanModel')) stop('ctPlotPosterior() needs a fit, ',
    'which has a posterior; plot(model) plots the priors of a model on its own.',
    call. = FALSE)
  if(!inherits(obj, c('ctStanFit','ctJuliaFit'))) stop('obj must be a fit from ctFit().',
    call. = FALSE)
  if(.ctFitIsJulia(obj) && ...length()) stop('ctPlotPosterior() passes ... to the ',
    'stan prior preparation; a backend=\'julia\' fit takes no further arguments.',
    call. = FALSE)
  pars <- if(.ctFitIsJulia(obj)) .ctPlotPosteriorDrawsJulia(obj, priorsamples) else
    .ctPlotPosteriorDrawsStan(obj, priorsamples, ...)

  if(rows[1]=='all') rows<-seq_along(pars)
  nplots<-ceiling(length(rows) / npp)
  if(1==99) Par.Value <- type <- quantity <- Density <- NULL
  quantity <- c('Posterior','Prior')
  xlimsindex <- if(priorwidth) 'all' else 1
  plots <- list()
  for(ploti in 1:nplots){
    dat <- data.table(quantity='',Par.Value=0, Density=0,type='',param='')
    for(ri in if(length(rows) > 1) rows[as.integer(cut_number(rows,nplots))==ploti] else rows){
      for(q in pars[[ri]]$quantities){
        dens <- ctDensityList(list(q$post[is.finite(q$post)], q$prior[is.finite(q$prior)]),
          probs=c(.05,.95),plot=FALSE,xlimsindex=xlimsindex,cut=q$cut,smoothness=smoothness)
        for(i in 1:length(dens$density)){
          dat <- rbind(dat,data.table(quantity=quantity[i],Par.Value=dens$density[[i]]$x,
            Density=dens$density[[i]]$y, type=q$type,param=pars[[ri]]$name))
        }
      }
    }
    dat <- dat[-1,]

    plots<-c(plots,list(
      ggplot(dat,aes(x=Par.Value,fill=quantity,ymax=Density,y=Density) )+
        geom_line(alpha=.3) +
        geom_ribbon(alpha=.4,ymin=0) +
        scale_fill_manual(values=c('red','blue')) +
        theme_minimal()+
        theme(legend.title = element_blank(),
          panel.grid.minor = element_line(linewidth = 0.1), panel.grid.major = element_line(linewidth = .2),
          strip.text.x = element_text(margin = margin(.05, 0, .1, 0, "cm"))) +
        facet_wrap(vars(type,param),scales='free')
    ))
  }

  if(plot) {
    firstplot=TRUE
    lapply(plots,function(x){
      if(wait && !firstplot) readline("Press [return] for next plot.")
      firstplot <<- FALSE
      suppressWarnings(print(x))
    })
    return(invisible(NULL))
  } else return(plots)
}

# What ctPlotPosterior() draws, one entry per free population parameter in
# parameter order: its name, and posterior and prior draws on the natural scale
# of its population mean and, where it varies, of its population sd. The two
# backends reach these differently; the plot is the same.
.ctPlotPosteriorDrawsStan <- function(obj, priorsamples, ...) {
  ps <- cbind(obj$setup$popsetup, obj$setup$popvalues)
  # `matrix < 11` is this plot's own restriction to the base matrices.
  ps <- ps[.ctMatsetupFreeRows(ps, defining = TRUE) & ps$matrix < 11,]
  ps <- ps[!duplicated(ps$param),]
  ps <- ps[order(ps$param),]
  # The internal, not the deprecated wrapper: `parsonly` is not part of
  # ctGenerate(fromPriors=TRUE), which returns data, and a user asking for a
  # posterior plot should not be told an internal call is deprecated.
  priors <- .ctGenerateFromPriors(cts = obj,parsonly=TRUE,nsamples=priorsamples,...)
  priors <- priors$stanfit$transformedpars
  posteriors <- ctExtract(obj)
  lapply(seq_len(nrow(ps)), function(ri) {
    pari <- ps[ri,'param']
    out <- list(name = ps$parname[ri], quantities = list(list(type = 'Pop. Mean',
      cut = TRUE, post = posteriors$popmeans[,pari], prior = priors$popmeans[,pari])))
    if(ps[ri,'indvarying'] > 0) out$quantities[[2L]] <- list(type = 'Pop. SD',
      cut = FALSE, post = posteriors$popsd[,ps$indvarying[ri]],
      prior = priors$popsd[,ps$indvarying[ri]])
    out
  })
}

# The julia fit's draws go through the same transforms summary() reads them
# with, and so do draws from the full raw prior, built over this fit's own
# parameter layout (reduced rank, own variances and grouping levels included)
# by the spec the prior warm-up uses.
.ctPlotPosteriorDrawsJulia <- function(fit, priorsamples) {
  post <- .ctBackendRawSamples(fit)
  if(is.null(post) || NROW(post) < 2L) stop('This fit carries no posterior ',
    'draws to plot; ctFitUncertainty(fit) computes them.', call. = FALSE)
  npar <- length(fit$estimate$raw)
  prior <- .ctPriorRawDraws(fit$standata, npar, priorsamples,
    spec = .ctJuliaFullPriorSpec(.ctBackendSpec(fit), fit$standata, npar))
  meanpost <- .ctBackendPopMeanSamples(fit, samples = post)$values
  meanprior <- .ctBackendPopMeanSamples(fit, samples = prior)$values
  levelspost <- .ctBackendConstrained(fit, post)$randomeffectlevels
  levelsprior <- .ctBackendConstrain(fit, prior)$randomeffectlevels
  multilevel <- length(levelspost) > 1L
  lapply(colnames(meanpost), function(p) {
    out <- list(name = p, quantities = list(list(type = 'Pop. Mean', cut = TRUE,
      post = meanpost[,p], prior = meanprior[,p])))
    for(l in seq_along(levelspost)) {
      sdpost <- levelspost[[l]]$popsd
      if(is.null(sdpost) || !p %in% colnames(sdpost)) next
      out$quantities[[length(out$quantities) + 1L]] <- list(
        type = if(multilevel) paste0('Pop. SD (', levelspost[[l]]$level, ')') else 'Pop. SD',
        cut = FALSE, post = sdpost[,p], prior = levelsprior[[l]]$popsd[,p])
    }
    out
  })
}

#' @export
ctStanPlotPost <- ctPlotPosterior

