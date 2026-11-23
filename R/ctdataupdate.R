#' Regenerate ctsem's shipped test data and/or test fits
#'
#' Regenerates \code{data/ctsemtestdat.rda}, \code{data/ctstantestfit.rda}
#' and/or \code{data/ctsemtestfit.rda}, the fixtures the examples and a number
#' of tests use. The julia fit needs Julia.
#'
#' Deliberately not exported: it overwrites data shipped with the package,
#' so it stays reachable only as \code{ctsem:::ctdataupdate()}.
#'
#' Called bare (no arguments) this behaves exactly as it always has: it asks
#' for interactive \code{T}/\code{F} confirmation and regenerates all three
#' objects. \code{confirm} and \code{what} exist so automated tooling can call
#' it without hanging on \code{readline()} and without moving
#' \code{ctsemtestdat} (which many tests pin numbers against) when only
#' \code{ctstantestfit} needs a refresh -- generation behaviour has changed
#' more than once this cycle, and refreshing the fit alone is normally what
#' that calls for.
#'
#' @param what Character vector, some of \code{c("data", "fit", "juliafit")}.
#'   Which object(s) to regenerate: the data, the stan fit
#'   \code{ctstantestfit}, the julia fit \code{ctsemtestfit}. Defaults to all
#'   three. Without \code{"data"}, the fits use the currently-shipped
#'   \code{ctsemtestdat} rather than regenerating it.
#' @param confirm Ask for interactive \code{T}/\code{F} confirmation before
#'   doing anything, exactly as this always has. Defaults to \code{TRUE};
#'   pass \code{FALSE} for a non-interactive/scripted call.
#' @param forcerecompile logical. For development purposes.
#' @keywords internal
#' @noRd
ctdataupdate <- function(what = c("data", "fit", "juliafit"), confirm = TRUE,
  forcerecompile = FALSE) {
  what <- match.arg(what, c("data", "fit", "juliafit"), several.ok = TRUE)

  run <- function() {
    set.seed(1)

    if ("data" %in% what) {
      Tpoints=10
      n.manifest=2
      n.TDpred=1
      n.TIpred=3
      n.latent=2
      n.subjects=30
      tipredEffect <- matrix(c(.5,0,0,-.7,0,2),nrow=2)
      tipredVar <- matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3)
      gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
      n.TDpred=n.TDpred,
        n.TIpred=0,
        n.manifest=n.manifest,
        MANIFESTVAR=diag(0.5,2),
        TDPREDMEANS=matrix(round(exp(rnorm(n.TDpred*(Tpoints),-1.9,1)),0),
         nrow=n.TDpred*(Tpoints)),
         TDPREDEFFECT = matrix(c(1,-1),ncol=1),
        LAMBDA=diag(1,2),
        DRIFT=matrix(c(-.3,.2,0,-.2),nrow=2),
        DIFFUSION=matrix(c(2,1,0,2),2),
        CINT=matrix(c(0,0),nrow=2),
        T0MEANS=matrix(10,ncol=1,nrow=2),
        T0VAR=diag(1,2))

      for(i in seq_len(n.subjects)){
        gm_i <- gm
        tipreds_i <- tipredVar %*% rnorm(n.TIpred)
        gm_i$CINT <- gm$CINT + tipredEffect %*% tipreds_i
        dat_i <- ctGenerate(gm_i,n=1,burnin=3,logdtsd=.4,dtmean = .3)
        dat_i[,'id'] <- i
        dat_i <- cbind(dat_i, matrix(tipreds_i, nrow=nrow(dat_i), ncol=n.TIpred, byrow=TRUE))
        if(i == 1) ctsemtestdat <- dat_i else ctsemtestdat <- rbind(ctsemtestdat, dat_i)
      }
      colnames(ctsemtestdat)[(ncol(ctsemtestdat)-n.TIpred+1):ncol(ctsemtestdat)] <- paste0('TI',1:n.TIpred)

      ctsemtestdat[2,'Y1'] <- NA
      ctsemtestdat[ctsemtestdat[,'id']==2,'TI1'] <- NA
      ctsemtestdat[2,'TD1'] <- NA

      save(ctsemtestdat,file=file.path('data','ctsemtestdat.rda'))
    } else {
      # Fits requested without the data: fit the currently-shipped data rather
      # than regenerating it, so a fit-only refresh cannot move ctsemtestdat
      # out from under tests that pin numbers against it.
      e <- new.env(parent = emptyenv())
      utils::data("ctsemtestdat", package = "ctsem", envir = e)
      ctsemtestdat <- e$ctsemtestdat
    }

    if ("fit" %in% what) {
      checkm<-ctModel(
        type='ct',
        n.latent=2,n.TDpred=1,n.TIpred=1,n.manifest=2,
        MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
        MANIFESTMEANS=0,
        DRIFT=c('dr1','dr12','dr21||||TI1','dr22'),
        DIFFUSION=c('diff11',0,'diff21','diff22||||TI1'),
        CINT=matrix(c('cint1||||TI1','cint2||||TI1'),ncol=1),
        LAMBDA=diag(2),tipredDefault=FALSE)

      ctstantestfit<-ctFit(ctsemtestdat,checkm,cores=1,inits=0,backend='stan',
        optimize = TRUE,optimcontrol=list(finishsamples=20,stochastic=T,tol=1e-5),priors=TRUE)

      ctstantestfit <- ctGenerateFromFit(ctstantestfit,nsamples = 20,fullposterior = TRUE,cores=1)

      save(ctstantestfit,file=file.path('data','ctstantestfit.rda'))
    }

    if ("juliafit" %in% what) {
      # The default model, fitted as ctFit()'s help does. It keeps its draws
      # pushed through the model's transforms, at the estimate too, so the
      # examples that summarise it run where Julia is missing (CRAN).
      # test-ctsemtestfit.R says when the engine no longer agrees with it.
      set.seed(1)
      model <- ctModel(LAMBDA = diag(2), manifestNames = c('Y1', 'Y2'),
        latentNames = c('eta1', 'eta2'), TDpredNames = 'TD1')
      ctsemtestfit <- ctFit(ctsemtestdat, model, priors = TRUE, cores = 1)
      save(ctsemtestfit, file = file.path('data', 'ctsemtestfit.rda'))
    }

    written <- file.path('data', paste0(
      c(data = 'ctsemtestdat', fit = 'ctstantestfit', juliafit = 'ctsemtestfit')[what],
      '.rda'))
    paths <- sort(Sys.glob(written))
    tools::resaveRdaFiles(paths)
    invisible(NULL)
  }

  if (isTRUE(confirm)) {
    # Unchanged from before: `if(continue)` on the raw readline() string,
    # which R coerces via the same rule `if()` always has -- "T"/"TRUE"
    # proceeds, "F"/"FALSE" is silently skipped, anything else errors.
    # `confirm = FALSE` is the escape hatch for a caller that cannot answer
    # this prompt at all; it does not change what answering it does.
    message(paste0('Updating from ',(getwd()),', continue T / F?'))
    continue <- readline()
    if(continue){
      run()
    }
  } else {
    run()
  }
}
