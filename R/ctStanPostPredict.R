# Posterior predictive checking -------------------------------------------
#
# One entry point, ctPostPredPlots(), producing two families of panel:
#
#   structure   -- does the model reproduce the joint shape of the data?
#                  Marginal densities, the value/time and value/occasion
#                  clouds, and the phase plane (rate of change against current
#                  value) that reveals nonlinear dynamics a linear drift cannot
#                  produce. These came from the old ctPostPredict().
#
#   calibration -- is the predictive distribution the right width, in the right
#                  place, everywhere? Predicted interval against observed
#                  value, interval coverage, the probability integral
#                  transform, and the PIT resolved by time interval and by
#                  subject. These came from the old ctPostPredPlots().
#
# ctPostPredict() and ctStanPostPredict() remain as aliases.
#
# Colour convention, used by every panel here: the model is blue and the
# observed data is black. The two functions merged here disagreed about which
# was which -- and the bivariate panels had the two data sources swapped
# outright, so a plot labelled 'Observed' was showing model draws.

.ctPostPredCols <- c(Model = "#2166AC", Observed = "#111111")

# Panel groups. 'all' is both, in this order.
.ctPostPredGroups <- list(
  structure = c('Density', 'ValueByTime', 'ValueByOccasion',
    'ChangeByValue', 'ChangeByTime'),
  calibration = c('PredictedVsObserved', 'IntervalCoverage', 'PIT',
    'CalibrationByInterval', 'CalibrationBySubject', 'SubjectLogLik')
)

# The explanatory notes. Each panel carries one as a ggplot caption unless
# `notes=FALSE`. They say what the panel is and, more usefully, what a
# departure from the reference means -- a plot that only says 'here is a
# number' leaves the reader to guess which direction is bad.
.ctPostPredNotes <- function(interval, ndraws) {
  pc <- paste0(round(interval * 100, 1), '%')
  outpc <- round((1 - interval) * 100, 1)
  list(
    Density = paste0(
      'Marginal distribution of each variable. Blue: values generated from the ',
      'fitted model (', ndraws, ' draws). Black: the observed data. A shift in ',
      'location or width means the model does not reproduce the marginal, ',
      'before any question of dynamics arises.'),
    ValueByTime = paste0(
      'Shaded: where the model puts 50, 80 and 95% of its predictive mass, from ',
      ndraws, ' draws. Points: the observed data, orange outside the 95% region. ',
      'A systematic drift of the points off the regions means the trend is wrong; ',
      'regions much wider or narrower than the scatter mean the variance is.'),
    ValueByOccasion = paste0(
      'As the time panel, but against occasion number within subject rather ',
      'than elapsed time. Separates a misfit that follows measurement order ',
      '(fatigue, practice, panel conditioning) from one that follows time.'),
    ChangeByValue = paste0(
      'Phase plane: rate of change against the value it started from, shaded by ',
      'the 50, 80 and 95% predictive regions. A linear model implies a ',
      'band sloping down at the drift rate, the same everywhere. Observed points ',
      'that curve, or that flatten at the extremes while the regions do not, ',
      'indicate nonlinear dynamics. See ctPhasePortrait() for the model side of ',
      'this on its own.'),
    ChangeByTime = paste0(
      'Rate of change against time, plotted at the midpoint of each interval and ',
      'shaded by the 50, 80 and 95% predictive regions. A scatter that ',
      'widens or narrows over time while the regions do not points at the ',
      'diffusion rather than the drift.'),
    PredictedVsObserved = paste0(
      'Each observation against the central ', pc, ' predictive interval for ',
      'it, ordered along the x axis by predicted median so the interval reads ',
      'as a band. Orange points fall outside their own interval; about ',
      outpc, '% should.'),
    IntervalCoverage = paste0(
      'The same exceedances as a rate, in ten equal-count bins of predicted ',
      'median, so a misfit concentrated among high or low predictions is ',
      'visible. Dashed line: the nominal ', outpc, '%. Bars: 95% Wilson ',
      'interval for each bin.'),
    PIT = paste0(
      'Probability integral transform -- for each observation, the fraction of ',
      'the ', ndraws, ' generated values falling below it. Under a correctly ',
      'calibrated model this is uniform. A hump in the middle means the ',
      'predictive distribution is too wide, peaks at both ends mean it is too ',
      'narrow, and a tilt means a location bias. Grey lines: where 95% of bin ',
      'counts should fall under uniformity.'),
    CalibrationByInterval = paste0(
      'Mean PIT against the time interval since the previous observation, in ',
      'equal-count bins; 0.5 is correct. A departure that depends on the ',
      'interval implicates the dynamics -- drift or diffusion -- rather than ',
      'the measurement model, which cannot know the interval.'),
    CalibrationBySubject = paste0(
      'Mean PIT per subject, sorted; 0.5 is correct, and subjects at the ends ',
      'are the ones the model fits worst. Grey band: the range expected if a ',
      'subject\'s own observations were independent. They are not -- a random ',
      'effect or a persistent state makes them agree -- so the true band is ',
      'wider, and this one is a guide rather than a test.'),
    SubjectLogLik = paste0(
      'Distribution across subjects of each subject total log likelihood, as a ',
      'cumulative curve. Band: the same curve from ', ndraws, ' replicate ',
      'datasets the model generated, each with the same subjects and the same ',
      'observation schedule, so it is the right reference however few subjects ',
      'there are. The observed curve should sit inside it. Left of the band ',
      'means subjects are less likely under the model than it expects -- ',
      'systematic misfit spread across occasions; right means the model is ',
      'fitting them better than it should, which usually means overfitting or ',
      'a variance component absorbing the noise. The rug marks the individual ',
      'subjects; CalibrationBySubject names which are which. The observed ',
      'totals are evaluated at the point estimate, so they carry no parameter ',
      'uncertainty while the replicates do.')
  )
}

# The observed row log-likelihood is evaluated at the point estimate while the
# generated ones vary over draws, so the LogLik facet is not symmetric between
# the two sides. Worth saying on the panels that show it rather than leaving
# the reader to assume otherwise.
.ctPostPredLLNote <- paste0(
  ' LogLik is each row\'s log-likelihood, and its distribution should match ',
  'like any other quantity -- but the observed side is evaluated at the point ',
  'estimate, so it carries no parameter uncertainty while the generated side does.')

.ctPostPredCaption <- function(g, key, notes, interval, ndraws, hasll = FALSE) {
  if (!isTRUE(notes)) return(g)
  txt <- .ctPostPredNotes(interval, ndraws)[[key]]
  if (is.null(txt)) return(g)
  if (hasll) txt <- paste0(txt, .ctPostPredLLNote)
  g + labs(caption = paste(strwrap(txt, width = 105), collapse = '\n')) +
    theme(plot.caption = element_text(hjust = 0, size = rel(.75),
      colour = 'grey25', margin = margin(t = 6)))
}


#' Create a data.table to compare data generated from a ctsem fit with the original data.
#'
#' This function allows for easy comparison of data generated from a fitted ctsem model
#' with the original data used to fit the model. It provides options to include residuals
#' in the comparison.
#'
#' @param fit A fitted ctsem model, from either backend.
#' @param residuals If set to TRUE, includes standardised prior residuals in the comparison,
#' as extra rows whose variable names are suffixed \code{' std. res.'}. Needs one Kalman
#' filter pass per generated dataset, so it is much slower than the default.
#' @param nsamples Integer or NA. Number of generated datasets to use. If the fit carries
#' more than this, draws are subsampled; if it carries none, this many are generated. NA
#' uses whatever the fit already has, or the \code{\link{ctGenerateFromFit}} default.
#'
#' @return A long data table with one row per generated draw x data row x variable, carrying
#' the generated \code{value}, the corresponding \code{obsValue}, and the \code{id},
#' \code{Time} and \code{TimeInterval} of that data row. \code{LogLik} appears as an extra
#' 'variable', holding the row-wise log likelihood.
#'
#' @seealso \code{\link{ctPostPredPlots}}, which turns this into diagnostic plots.
#'
#' @examples
#' data_comparison <- ctPostPredData(ctstantestfit)
#'
#' @export
ctPostPredData <- function(fit, residuals = FALSE, nsamples = NA){
  if(is.null(fit$generated)){
    fit <- if(is.na(nsamples)) ctGenerateFromFit(fit) else
      ctGenerateFromFit(fit, nsamples = nsamples)
  }
  if(!is.na(nsamples) && dim(fit$generated$Y)[1] > nsamples){
    keep <- sort(sample.int(dim(fit$generated$Y)[1], nsamples))
    fit$generated$Y <- fit$generated$Y[keep, , , drop = FALSE]
    fit$generated$llrow <- fit$generated$llrow[keep, , drop = FALSE]
  }

  ll <- melt(
    data.table(sample=1:nrow(fit$generated$llrow),(fit$generated$llrow)),
    id.vars = c('sample'),variable.name = 'row')[order(sample),]
  ll[['row']] <- as.numeric(ll[['row']])
  ll <- cbind(ll[,1:2],variable='LogLik',ll[,3])

  dat=fit$generated$Y
  dat <- as.data.table(dat,na.rm = FALSE)
  v1name <- colnames(dat)[!colnames(dat) %in% c('sample','row','value')]
  data.table::setnames(dat,v1name,'variable')

  # `as.data.table()` on the generated array keys the result <sample,row,variable>
  # from the character dimnames, so the table arrives in LEXICOGRAPHIC row order
  # -- 1, 10, 11, 2, 3, ... `row` is retyped to numeric on the next line, but the
  # stale key survives the retype, and setorder() then sees a request it believes
  # is already satisfied and does nothing. Silently. The result was a TimeInterval
  # diffed across a scrambled sequence for every subject whose rows span a
  # digit-count boundary, surfacing as negative time intervals in the by-interval
  # panel. Values were never misplaced -- as.data.table carries each label with
  # its own value -- so TimeInterval was the whole of the damage. Drop the key first.
  setkey(dat, NULL)
  dat[['sample']] <- as.integer(dat[['sample']])
  dat[['row']] <- as.numeric(dat[['row']])
  setorder(dat,'sample','row')
  setkey(ll, NULL)
  setorder(ll,'sample','row')
  dat <- rbind(dat,ll)

  # The observed data, the row-to-subject map, the times and the fitted row
  # likelihoods, from whichever backend produced the fit.
  idmap <- .ctFitIdMap(fit)
  rowsubject <- .ctFitRowSubject(fit)
  rowtime <- .ctFitRowTime(fit)
  truedat <- .ctFitObservedY(fit)

  dat[,id:= rowsubject[row]]
  dat[,id:=idmap[[1]][match(id,idmap[[2]])]]
  dat[,Time:=rowtime[row]]
  dat[,TimeInterval:=c(NA,diff(Time)),by=.(sample,id,variable)]

  # The cheap guard for the ordering trap above: time runs forwards within a
  # subject, so a negative interval means the table was diffed out of order.
  if(any(dat$TimeInterval < 0, na.rm = TRUE)) stop(
    'ctPostPredData(): negative time intervals -- the generated data was ordered ',
    'by row label rather than row number. This is a bug, please report it.',
    call. = FALSE)

  truell <- .ctFitObservedRowLoglik(fit)
  truell[apply(truedat,1,function(x) all(is.na(x)))] <- NA #ensure missings propagate to likelihood also
  dat=merge(dat, #generated
    melt(data.table(row=1:max(dat$row),cbind(truedat, #true
      LogLik=truell)),#fitted
      id.vars='row',value.name = 'obsValue'),by=c('row','variable'),all=T)

  if(residuals){
    ft <- fit
    stderrprior <- list()
    for(i in 1:max(ll[['sample']])){
      ft <- .ctFitReplaceY(fit, matrix(fit$generated$Y[i,,],ncol=ncol(truedat)))
      stderrprior[[i]] <- data.table(sample=i,suppressMessages(meltkalman(ctKalmanArray(ft,standardisederrors = TRUE))))[Element %in% 'errstdprior',.(sample,Row,value,Obs)]
    }
    stderrprior <- rbindlist(stderrprior)
    stderrpriorObs<-data.table(suppressMessages(
      meltkalman(ctKalmanArray(fit,standardisederrors = TRUE))))[Element %in% 'errstdprior',.(Row,value,Obs)]
    setnames(stderrpriorObs,'value','obsValue')
    stderrprior<-merge(stderrprior,stderrpriorObs)
    setnames(stderrprior,c('Row','Obs'),c('variable','row'))
    stderrprior[,variable:=paste0(variable,' std. res.')]
    # The residual rows need the same id/time columns the rest of `dat` carries,
    # derived from the row index the same way. Without them the rbind below has
    # never been able to run -- this branch was broken for every backend,
    # including stan. Ordered explicitly before differencing, for the same
    # reason the generated block above is.
    stderrprior[,id:= rowsubject[row]]
    stderrprior[,id:=idmap[[1]][match(id,idmap[[2]])]]
    stderrprior[,Time:=rowtime[row]]
    setkey(stderrprior, NULL)
    setorder(stderrprior, sample, variable, row)
    stderrprior[,TimeInterval:=c(NA,diff(Time)),by=.(sample,id,variable)]
    dat <- rbind(dat,stderrprior[,colnames(dat),with=FALSE])
  }

  return(dat)
}


# The generated and observed data as aligned arrays, for the panels that need
# per-row structure rather than the long table -- the derivatives especially,
# where a row has to be differenced against the right neighbouring row.
.ctPostPredArrays <- function(fit, nsamples = NA, datarows = 'all'){
  if(is.null(fit$generated)){
    fit <- if(is.na(nsamples)) ctGenerateFromFit(fit) else
      ctGenerateFromFit(fit, nsamples = nsamples)
  }
  Ygen <- fit$generated$Y
  if(!is.na(nsamples) && dim(Ygen)[1] > nsamples){
    Ygen <- Ygen[sort(sample.int(dim(Ygen)[1], nsamples)), , , drop = FALSE]
  }
  Yobs <- .ctFitObservedY(fit)
  time <- .ctFitRowTime(fit)
  subject <- .ctFitRowSubject(fit)

  nrows <- nrow(Yobs)
  if(identical(datarows[1], 'all')) datarows <- seq_len(nrows)
  datarows <- as.integer(datarows)
  if(any(is.na(datarows) | datarows < 1 | datarows > nrows)) stop(
    'datarows must index rows of the data: 1 to ', nrows, call. = FALSE)

  # Subset ONCE. The old ctPostPredict() subset Ygen and Ydat and then indexed
  # the results by datarows a second time, so any actual subset threw
  # 'subscript out of bounds' -- the argument had never worked.
  list(
    Ygen = Ygen[, datarows, , drop = FALSE],
    Yobs = Yobs[datarows, , drop = FALSE],
    time = time[datarows],
    subject = subject[datarows],
    varnames = colnames(Yobs),
    ndraws = dim(Ygen)[1]
  )
}


# A bivariate panel: the model's predictive region as shading, the observed data
# as points on top of it.
#
# The helper this replaces took the model as group 1 and the observed as group 2
# while labelling group 1 'Observed', so every panel it drew was reversed -- the
# red cloud was the data and the blue points were the model. It also estimated
# the model density from 300 subsampled draws (the group-1 cap) while giving the
# observed data the 50000 budget, which is backwards twice over.
#
# The shading is highest-density regions, not a raw density ramp. A ramp keyed to
# the normalised density fades the model's tails to invisible, so a model that
# does reach the outlying observations looks as though it cannot -- the reader
# sees a compact blob and points well outside it, and reads misfit that is not
# there. An HDR band says something a reader can act on instead: 95% of the
# model's predictive mass lies inside the outer band, so about 5% of the
# observations should fall outside it. That proportion goes in the subtitle,
# which turns the panel from an impression into a number.
.ctPostPredHDR <- function(mx, my, hx, hy, resolution, xlim, ylim,
  probs = c(.5, .8, .95)){
  k <- try(MASS::kde2d(mx, my, h = c(hx, hy), n = resolution,
    lims = c(xlim, ylim)), silent = TRUE)
  if(inherits(k, 'try-error') || !all(is.finite(k$z))) return(NULL)
  cell <- diff(k$x[1:2]) * diff(k$y[1:2])
  z <- as.vector(k$z)
  ord <- order(z, decreasing = TRUE)
  cum <- cumsum(z[ord]) * cell
  # Density height at which the enclosed mass first reaches each probability.
  # Anything at or above that height is inside the region.
  lev <- vapply(probs, function(p){
    i <- which(cum >= p)
    if(!length(i)) min(z) else z[ord][i[1]]
  }, numeric(1))
  list(k = k, lev = lev, probs = probs)
}

.ctPostPredDensity2d <- function(obsx, obsy, modx, mody, xlab, ylab,
  title = '', resolution = 100, maxobs = 4000, maxmod = 60000,
  trim = c(.005, .995), trimx = trim, probs = c(.5, .8, .95)){

  if(1 == 99) x <- y <- z <- Outside <- NULL

  keep <- function(a, b){ i <- which(is.finite(a) & is.finite(b)); list(x = a[i], y = b[i]) }
  o <- keep(obsx, obsy); m <- keep(modx, mody)
  if(length(o$x) < 2 || length(m$x) < 10) return(NULL)
  if(length(unique(m$x)) < 2 || length(unique(m$y)) < 2) return(NULL)

  # Limits from both sources together, so neither is clipped in favour of the
  # other. The old version set the model's values to NA outside its own 0.5/99.5
  # percentiles but left the observed untouched, which narrowed the model's
  # apparent spread relative to the data it was being compared against.
  # `trimx` is the identity for panels whose x axis is a design variable -- time,
  # or occasion number. Those are known, not generated, so there is nothing
  # extreme to clip, and trimming them drops the first and last measurement
  # occasion off the edge of the plot.
  lim <- function(a, b, tr) range(quantile(c(a, b), tr, na.rm = TRUE))
  xlim <- lim(o$x, m$x, trimx); ylim <- lim(o$y, m$y, trim)
  if(diff(xlim) <= 0 || diff(ylim) <= 0) return(NULL)

  if(length(m$x) > maxmod){ i <- sample.int(length(m$x), maxmod); m$x <- m$x[i]; m$y <- m$y[i] }

  bw <- function(v){
    h <- try(MASS::bandwidth.nrd(v), silent = TRUE)
    if(inherits(h, 'try-error') || !is.finite(h) || h <= 0) diff(range(v)) / 10 else h
  }
  hx <- bw(m$x); hy <- bw(m$y)
  if(!is.finite(hx) || hx <= 0 || !is.finite(hy) || hy <= 0) return(NULL)

  hdr <- .ctPostPredHDR(m$x, m$y, hx, hy, min(resolution, 200L), xlim, ylim, probs)
  if(is.null(hdr)) return(NULL)
  k <- hdr$k; lev <- hdr$lev

  # Filled contours of the density at those heights, rather than a raster of
  # binned cells -- the region boundary is a smooth curve, not a staircase.
  # Breaks must ascend, and the density heights descend as the enclosed mass
  # grows, so they go in reversed.
  labs_band <- paste0(round(probs * 100), '%')
  brk <- unique(c(rev(lev), max(k$z) * 1.001))
  if(length(brk) < 2) return(NULL)
  grid <- data.table(
    x = rep(k$x, times = length(k$y)),
    y = rep(k$y, each = length(k$x)),
    z = as.vector(k$z))

  gi <- pmin(pmax(findInterval(o$x, k$x), 1L), length(k$x))
  gj <- pmin(pmax(findInterval(o$y, k$y), 1L), length(k$y))
  ondens <- k$z[cbind(gi, gj)]
  ongrid <- o$x >= xlim[1] & o$x <= xlim[2] & o$y >= ylim[1] & o$y <= ylim[2]
  inside <- ongrid & ondens >= lev[length(lev)]
  pcout <- round(100 * mean(!inside), 1)
  # Derived from `probs`, not written in: the outer region is whatever the last
  # entry says it is, and so is the rate that should fall outside it.
  outerpc <- round(100 * probs[length(probs)], 1)
  nominalpc <- round(100 * (1 - probs[length(probs)]), 1)

  nobs <- length(o$x)
  od <- data.table(x = o$x, y = o$y, Outside = !inside)
  if(nobs > maxobs) od <- od[sample.int(nobs, maxobs)]

  # Light to dark with increasing density, so the 50% core reads as the centre.
  # geom_contour_filled hands back bands in ascending height, i.e. widest region
  # first, so the labels run the same way.
  fills <- grDevices::colorRampPalette(
    c('#DCE7F1', .ctPostPredCols[['Model']]))(length(brk) - 1L)

  ggplot() +
    geom_contour_filled(data = grid, aes(x = x, y = y, z = z), breaks = brk) +
    scale_fill_manual(name = 'Model predictive region', values = fills,
      labels = rev(labs_band)[seq_len(length(brk) - 1L)]) +
    geom_point(data = od, aes(x = x, y = y, colour = Outside), size = .8, alpha = .7) +
    scale_colour_manual(name = '', guide = 'none',
      values = c('FALSE' = .ctPostPredCols[['Observed']], 'TRUE' = '#D95F02')) +
    # Expanded rather than flush: a point sitting exactly on the limit -- every
    # occasion-1 observation does -- is half outside a flush panel and reads as
    # absent.
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = xlab, y = ylab, title = title,
      subtitle = paste0('Points: observed', if(nobs > maxobs)
        paste0(' (', maxobs, ' of ', nobs, ' shown)') else '',
        '. ', pcout, '% fall outside the ', outerpc, '% region (orange); about ',
        nominalpc, '% should.')) +
    theme_bw() +
    theme(legend.position = 'bottom',
      plot.subtitle = element_text(size = rel(.8), colour = 'grey30'))
}


# Difference quotients for the rate-of-change panels: for lag k, every position
# j whose partner j+k belongs to the same subject.
.ctPostPredDiff <- function(arr, v, lag, jitter = 0){
  n <- length(arr$time)
  if(n <= lag) return(NULL)
  j <- seq_len(n - lag)
  dt <- arr$time[j + lag] - arr$time[j]
  ok <- arr$subject[j] == arr$subject[j + lag] & is.finite(dt) & dt > 0
  j <- j[ok]
  if(!length(j)) return(NULL)
  dt <- arr$time[j + lag] - arr$time[j]

  obsy <- (arr$Yobs[j + lag, v] - arr$Yobs[j, v]) / dt
  obsx <- arr$Yobs[j, v]

  S <- dim(arr$Ygen)[1]
  gen <- arr$Ygen[, , which(arr$varnames == v)]
  if(is.null(dim(gen))) gen <- matrix(gen, nrow = S)
  # Divide by the intervals matched to the ACTUAL number of draws. The old code
  # built this divisor with ncol=nsamples, the requested count, so whenever the
  # fit already carried a different number the matrix recycled and the rates
  # came out scaled by the wrong intervals.
  mody <- (gen[, j + lag, drop = FALSE] - gen[, j, drop = FALSE]) /
    matrix(dt, nrow = S, ncol = length(j), byrow = TRUE)
  modx <- gen[, j, drop = FALSE]

  # Jitter, when asked for, goes on BOTH sides. The old code jittered the
  # observed rates only, which widened the data relative to the model it was
  # being compared against.
  if(jitter > 0){
    s <- sd(arr$Yobs[, v], na.rm = TRUE)
    if(is.finite(s) && s > 0){
      obsx <- obsx + rnorm(length(obsx), 0, jitter * s)
      modx <- modx + rnorm(length(modx), 0, jitter * s)
    }
  }

  mid <- (arr$time[j] + arr$time[j + lag]) / 2
  list(obsx = obsx, obsy = obsy,
    modx = as.vector(modx), mody = as.vector(mody),
    midtime = mid, modmidtime = rep(mid, each = S))
}


#' Posterior predictive checks for a ctsem fit
#'
#' Compares data generated from the fitted model against the data it was fitted to, as a
#' set of diagnostic plots. Works for both the stan and julia backends.
#' \code{ctPostPredict} and \code{ctStanPostPredict} are maintained as aliases.
#'
#' @param fit A fitted ctsem model.
#' @param panels Which panels to produce. \code{'all'} (the default), \code{'structure'},
#' \code{'calibration'}, or a character vector of panel names -- see Details.
#' @param variables Character vector of variable names to restrict to, or NULL for all.
#' \code{'LogLik'} is a valid name here.
#' @param nsamples Integer or NA. Number of generated datasets to compare against. NA uses
#' whatever the fit already carries.
#' @param datarows Integer vector of data rows to use, or \code{'all'}.
#' @param diffsize Integer vector > 0. Lags, in observations, for the rate-of-change panels.
#' One pair of panels is produced per lag.
#' @param interval Numeric in (0,1). Width of the central predictive interval used by the
#' \code{PredictedVsObserved} and \code{IntervalCoverage} panels.
#' @param resolution Positive integer. Grid resolution of the shaded model density in the
#' bivariate panels.
#' @param jitter Non-negative numeric. Jitter added to the values in the rate-of-change
#' panels, as a proportion of each variable's standard deviation, applied to observed and
#' generated alike. Useful when data are rounded or ordinal. 0, the default, adds none.
#' @param residuals Logical. Include standardised prior residuals as extra variables. Slow --
#' it needs one Kalman pass per generated dataset.
#' @param notes Logical. If TRUE (the default), each panel carries a caption saying what it
#' shows and what a departure from its reference line means. Set FALSE for bare plots.
#' @param plot Logical. If TRUE, prints the panels and returns them invisibly. If FALSE (the
#' default) returns the list without printing.
#' @param wait Logical. If TRUE and \code{plot=TRUE}, waits for a keypress between panels.
#'
#' @details
#' Two families of panel are produced.
#'
#' \strong{structure} asks whether the model reproduces the shape of the data:
#' \code{Density} (marginal distributions), \code{ValueByTime} and \code{ValueByOccasion}
#' (the model-implied density of each variable against time and against occasion number,
#' with the observed data as points on top), and \code{ChangeByValue} and
#' \code{ChangeByTime} (rate of change against current value and against time).
#' \code{ChangeByValue} is the phase plane, and is the panel that shows nonlinear dynamics:
#' a linear model implies a straight band with a constant slope, so curvature in the
#' observed points relative to the model shading is evidence the drift is not linear.
#'
#' \strong{calibration} asks whether the predictive distribution is the right width and in
#' the right place: \code{PredictedVsObserved}, \code{IntervalCoverage}, \code{PIT}, and the
#' PIT resolved by \code{CalibrationByInterval} and \code{CalibrationBySubject}. A departure
#' that depends on the time interval implicates the dynamics rather than the measurement
#' model, since the measurement model cannot know the interval.
#'
#' The row-wise log likelihood appears throughout as an extra variable, \code{LogLik} --
#' its distribution should match between model and data like any other quantity. It is
#' omitted from the rate-of-change panels, where the derivative of a log likelihood with
#' respect to time is not a meaningful quantity. Note that the observed side of the LogLik
#' comparison is evaluated at the point estimate, so it carries no parameter uncertainty
#' while the generated side does.
#'
#' \code{SubjectLogLik} asks the same question one level up: each subject's total log
#' likelihood, summed over that subject's rows, compared as a cumulative distribution
#' across subjects against the same curve from each generated dataset. Aggregating within
#' a subject first is what makes a subject the model fits badly across many occasions show
#' up, rather than being spread thin across the row-wise panels. The comparison is against
#' replicate datasets of the same size rather than a smooth density, so it is the right
#' reference however few subjects there are.
#'
#' The approximation is only as good as the number of generated datasets, so a fit carrying
#' few draws gives a coarse picture -- the PIT in particular is discrete on
#' \code{nsamples + 1} values.
#'
#' @return A named list of ggplot objects, invisibly if \code{plot=TRUE}.
#'
#' @seealso \code{\link{ctPostPredData}} for the underlying table,
#' \code{\link{ctFitCheck}} for a broader dashboard, \code{\link{ctFitCheckCov}} for the
#' lagged covariance version of the same question, and \code{\link{ctPhasePortrait}} for the
#' model-implied phase plane on its own.
#'
#' @export
#'
#' @examples
#' \donttest{
#' plots <- ctPostPredPlots(ctstantestfit)
#' print(plots$PIT)
#' }
ctPostPredPlots <- function(fit, panels = 'all', variables = NULL,
  nsamples = NA, datarows = 'all', diffsize = 1, interval = .95,
  resolution = 100, jitter = 0, residuals = FALSE, notes = TRUE,
  plot = FALSE, wait = FALSE){

  if(!inherits(fit, c('ctStanFit','ctJuliaFit'))) stop('Not a ctsem fit object', call.=FALSE)
  if(!is.numeric(interval) || length(interval) != 1 || interval <= 0 || interval >= 1)
    stop('interval must be a single number between 0 and 1', call. = FALSE)

  panels <- .ctPostPredResolvePanels(panels)

  if(1 == 99){
    value <- obsValue <- PIT <- med <- lo <- hi <- Outside <- DataType <- NULL
    Bin <- Rate <- RateLo <- RateHi <- Mid <- Expected <- Count <- NULL
    BandLo <- BandHi <- Rank <- SubjMean <- TimeInterval <- Time <- NULL
    mid <- obsF <- gensub <- obssub <- NULL
    variable <- row <- id <- nout <- Nobs <- se <- NULL
  }

  dat <- ctPostPredData(fit, residuals = residuals, nsamples = nsamples)
  if(!identical(datarows[1], 'all')) dat <- dat[row %in% as.integer(datarows)]
  ndraws <- length(unique(dat$sample))

  # Manifest variables first, residuals next, LogLik last, so every facetted
  # panel orders its strips the same way and LogLik does not land in the middle.
  allvars <- unique(as.character(dat$variable))
  mans <- .ctFitModelObject(fit)$manifestNames
  lev <- c(intersect(mans, allvars),
    sort(setdiff(allvars, c(mans, 'LogLik'))),
    intersect('LogLik', allvars))
  if(!is.null(variables)){
    bad <- setdiff(variables, lev)
    if(length(bad)) stop('variables not present: ', paste(bad, collapse = ', '),
      '. Available: ', paste(lev, collapse = ', '), call. = FALSE)
    lev <- intersect(lev, variables)
    dat <- dat[as.character(variable) %in% lev]
  }
  dat[, variable := factor(as.character(variable), levels = lev)]
  hasll <- 'LogLik' %in% lev

  # Per-observation summaries. The PIT is a mid-rank fraction -- the old
  # `sum(obsValue > value)/.N` counted only strict exceedances, which is biased
  # downward whenever values tie (rounded, ordinal or binary data), and divided
  # by the group size including comparisons that were NA.
  qlo <- (1 - interval) / 2
  dat[, c('PIT','med','lo','hi') := {
    o <- obsValue[1L]
    v <- value[is.finite(value)]
    if(!length(v)) list(NA_real_, NA_real_, NA_real_, NA_real_) else
      list(if(is.finite(o)) (sum(v < o) + .5 * sum(v == o)) / length(v) else NA_real_,
        median(v), quantile(v, qlo, names = FALSE), quantile(v, 1 - qlo, names = FALSE))
  }, by = .(variable, row)]

  # One row per observation, which is the unit every calibration panel works in.
  obs <- unique(dat[, .(variable, row, id, Time, TimeInterval, obsValue, PIT, med, lo, hi)])
  obs[, Outside := is.finite(obsValue) & (obsValue > hi | obsValue < lo)]
  obs[, Rank := rank(med, ties.method = 'first'), by = variable]

  gglist <- list()
  llpanels <- c('Density','PredictedVsObserved','IntervalCoverage','PIT',
    'CalibrationByInterval','CalibrationBySubject')
  cap <- function(g, key) .ctPostPredCaption(g, key, notes, interval, ndraws,
    hasll = hasll && key %in% llpanels)

  # -- structure ----------------------------------------------------------

  if('Density' %in% panels){
    dd <- rbind(
      dat[, .(value = value, variable, DataType = 'Model')],
      obs[, .(value = obsValue, variable, DataType = 'Observed')])
    g <- ggplot(dd[is.finite(value)], aes(x = value, colour = DataType, fill = DataType)) +
      geom_density(alpha = .25, linewidth = .8) +
      scale_colour_manual(name = '', values = .ctPostPredCols) +
      scale_fill_manual(name = '', values = .ctPostPredCols) +
      facet_wrap(vars(variable), scales = 'free') +
      theme_bw() + labs(x = 'Value', y = 'Density') +
      theme(legend.position = 'bottom')
    gglist$Density <- cap(g, 'Density')
  }

  if(any(c('ValueByTime','ValueByOccasion','ChangeByValue','ChangeByTime') %in% panels)){
    arr <- .ctPostPredArrays(fit, nsamples = nsamples, datarows = datarows)
    S <- dim(arr$Ygen)[1]
    occ <- stats::ave(seq_along(arr$subject), arr$subject, FUN = seq_along)
    dvars <- intersect(arr$varnames, lev)

    for(v in dvars){
      gen <- arr$Ygen[, , which(arr$varnames == v)]
      if(is.null(dim(gen))) gen <- matrix(gen, nrow = S)

      if('ValueByTime' %in% panels){
        g <- .ctPostPredDensity2d(arr$time, arr$Yobs[, v],
          rep(arr$time, each = S), as.vector(gen),
          xlab = 'Time', ylab = v, title = v, resolution = resolution,
          trimx = c(0, 1))
        if(!is.null(g)) gglist[[paste0('ValueByTime_', v)]] <- cap(g, 'ValueByTime')
      }
      if('ValueByOccasion' %in% panels){
        g <- .ctPostPredDensity2d(occ, arr$Yobs[, v],
          rep(occ, each = S), as.vector(gen),
          xlab = 'Occasion within subject', ylab = v, title = v,
          resolution = resolution, trimx = c(0, 1))
        if(!is.null(g)) gglist[[paste0('ValueByOccasion_', v)]] <- cap(g, 'ValueByOccasion')
      }

      for(k in diffsize){
        if(!any(c('ChangeByValue','ChangeByTime') %in% panels)) break
        d <- .ctPostPredDiff(arr, v, k, jitter = jitter)
        if(is.null(d)) next
        sfx <- if(length(diffsize) > 1) paste0('_lag', k) else ''
        if('ChangeByValue' %in% panels){
          g <- .ctPostPredDensity2d(d$obsx, d$obsy, d$modx, d$mody,
            xlab = v, ylab = paste0('d(', v, ')/dt'),
            title = paste0(v, ', lag ', k), resolution = resolution)
          if(!is.null(g)) gglist[[paste0('ChangeByValue_', v, sfx)]] <- cap(g, 'ChangeByValue')
        }
        if('ChangeByTime' %in% panels){
          g <- .ctPostPredDensity2d(d$midtime, d$obsy, d$modmidtime, d$mody,
            xlab = 'Time (interval midpoint)', ylab = paste0('d(', v, ')/dt'),
            title = paste0(v, ', lag ', k), resolution = resolution, trimx = c(0, 1))
          if(!is.null(g)) gglist[[paste0('ChangeByTime_', v, sfx)]] <- cap(g, 'ChangeByTime')
        }
      }
    }
  }

  # -- calibration --------------------------------------------------------

  if('PredictedVsObserved' %in% panels){
    pv <- obs[is.finite(med)]
    pv[, DataType := ifelse(Outside, 'Observed, outside interval', 'Observed')]
    g <- ggplot(pv, aes(x = Rank)) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = 'Model interval'), alpha = .3) +
      geom_line(aes(y = med, colour = 'Model median'), linewidth = .7) +
      geom_point(aes(y = obsValue, colour = DataType), size = .8, alpha = .65) +
      scale_fill_manual(name = '', values = c('Model interval' = .ctPostPredCols[['Model']])) +
      scale_colour_manual(name = '', values = c(
        'Model median' = .ctPostPredCols[['Model']],
        'Observed' = .ctPostPredCols[['Observed']],
        'Observed, outside interval' = '#D95F02')) +
      facet_wrap(vars(variable), scales = 'free') +
      theme_bw() +
      labs(x = 'Observations, ordered by predicted median', y = 'Value') +
      theme(legend.position = 'bottom')
    gglist$PredictedVsObserved <- cap(g, 'PredictedVsObserved')
  }

  if('IntervalCoverage' %in% panels){
    cv <- obs[is.finite(obsValue) & is.finite(med)]
    if(nrow(cv)){
      cv[, Bin := .ctPostPredBin(med, 10), by = variable]
      cs <- cv[, .(Mid = mean(med), n = .N, nout = sum(Outside)), by = .(variable, Bin)]
      ci <- .ctPostPredBinomCI(cs$nout, cs$n)
      cs[, Rate := nout / n][, RateLo := ci[[1]]][, RateHi := ci[[2]]]
      g <- ggplot(cs, aes(x = Mid)) +
        geom_hline(aes(yintercept = 1 - interval, linetype = 'Nominal rate'), colour = 'grey20') +
        geom_linerange(aes(ymin = RateLo, ymax = RateHi), colour = .ctPostPredCols[['Model']]) +
        geom_point(aes(y = Rate), colour = .ctPostPredCols[['Model']], size = 1.8) +
        scale_linetype_manual(name = '', values = c('Nominal rate' = 'dashed')) +
        facet_wrap(vars(variable), scales = 'free_x') +
        coord_cartesian(ylim = c(0, 1)) +
        theme_bw() +
        labs(x = 'Predicted median (bin mean)',
          y = paste0('Proportion outside the ', round(interval * 100, 1), '% interval')) +
        theme(legend.position = 'bottom')
      gglist$IntervalCoverage <- cap(g, 'IntervalCoverage')
    }
  }

  if('PIT' %in% panels){
    # The PIT is discrete on ndraws+1 values, so bin it on exactly that grid. A
    # kernel density with adjust=.25, which is what this panel used to be, turns
    # that discreteness into a comb of spurious peaks.
    nb <- max(2L, min(as.integer(ndraws) + 1L, 40L))
    p <- obs[is.finite(PIT)]
    if(nrow(p)){
      p[, Bin := pmin(floor(PIT * nb) + 1L, nb)]
      ps <- p[, .(Count = .N), by = .(variable, Bin)]
      ps <- ps[CJ(variable = unique(p$variable), Bin = seq_len(nb), unique = TRUE),
        on = .(variable, Bin)]
      ps[is.na(Count), Count := 0L]
      ps <- merge(ps, p[, .(Nobs = .N), by = variable], by = 'variable')
      ps[, Expected := Nobs / nb]
      ps[, BandLo := qbinom(.025, Nobs, 1 / nb)]
      ps[, BandHi := qbinom(.975, Nobs, 1 / nb)]
      ps[, Mid := (Bin - .5) / nb]
      # The band as two lines over the bars, not as a filled ribbon under them.
      # Under them it is hidden by every bar tall enough to matter; over them as
      # a fill it tints the bars inside the band, which reads as highlighting
      # exactly the ones that are fine. Lines leave the bars alone and still
      # show at a glance which ones break out.
      g <- ggplot(ps, aes(x = Mid)) +
        geom_col(aes(y = Count), fill = .ctPostPredCols[['Model']], width = 1 / nb) +
        geom_line(aes(y = BandLo), colour = 'grey25', linewidth = .4) +
        geom_line(aes(y = BandHi), colour = 'grey25', linewidth = .4) +
        geom_line(aes(y = Expected), colour = 'grey10', linetype = 'dashed') +
        facet_wrap(vars(variable), scales = 'free_y') +
        theme_bw() +
        labs(x = 'PIT: fraction of generated values below the observed one',
          y = 'Observations') +
        theme(legend.position = 'none')
      gglist$PIT <- cap(g, 'PIT')
    }
  }

  if('CalibrationByInterval' %in% panels){
    di <- obs[is.finite(PIT) & is.finite(TimeInterval) & TimeInterval > 0]
    if(nrow(di)){
      di[, Bin := .ctPostPredBin(TimeInterval, 10), by = variable]
      ds <- di[, .(Mid = mean(TimeInterval), SubjMean = mean(PIT),
        se = sd(PIT) / sqrt(.N)), by = .(variable, Bin)]
      ds[!is.finite(se), se := 0]
      g <- ggplot(ds, aes(x = Mid)) +
        geom_hline(aes(yintercept = .5, linetype = 'Correct (0.5)'), colour = 'grey20') +
        geom_linerange(aes(ymin = SubjMean - 1.96 * se, ymax = SubjMean + 1.96 * se),
          colour = .ctPostPredCols[['Model']]) +
        geom_point(aes(y = SubjMean), colour = .ctPostPredCols[['Model']], size = 1.8) +
        scale_linetype_manual(name = '', values = c('Correct (0.5)' = 'dashed')) +
        facet_wrap(vars(variable), scales = 'free_x') +
        coord_cartesian(ylim = c(0, 1)) +
        theme_bw() +
        labs(x = 'Time interval since previous observation (bin mean)', y = 'Mean PIT') +
        theme(legend.position = 'bottom')
      gglist$CalibrationByInterval <- cap(g, 'CalibrationByInterval')
    }
  }

  if('CalibrationBySubject' %in% panels){
    sb <- obs[is.finite(PIT), .(SubjMean = mean(PIT), n = .N), by = .(variable, id)]
    if(nrow(sb)){
      sb[, Rank := rank(SubjMean, ties.method = 'first'), by = variable]
      # Under independence the mean of n uniforms has sd sqrt(1/(12n)). Within a
      # subject the PITs are not independent -- a random effect or a persistent
      # state makes them agree -- so this band is narrower than the truth. Drawn
      # as a guide and labelled as one in the caption.
      sb[, BandLo := .5 - 1.96 * sqrt(1 / (12 * n))]
      sb[, BandHi := .5 + 1.96 * sqrt(1 / (12 * n))]
      g <- ggplot(sb, aes(x = Rank)) +
        geom_ribbon(aes(ymin = BandLo, ymax = BandHi), fill = 'grey70', alpha = .45) +
        geom_hline(aes(yintercept = .5, linetype = 'Correct (0.5)'), colour = 'grey20') +
        geom_point(aes(y = SubjMean), colour = .ctPostPredCols[['Model']], size = 1.4) +
        scale_linetype_manual(name = '', values = c('Correct (0.5)' = 'dashed')) +
        facet_wrap(vars(variable), scales = 'free_x') +
        coord_cartesian(ylim = c(0, 1)) +
        theme_bw() +
        labs(x = 'Subjects, ordered by mean PIT', y = 'Mean PIT') +
        theme(legend.position = 'bottom')
      gglist$CalibrationBySubject <- cap(g, 'CalibrationBySubject')
    }
  }

  if('SubjectLogLik' %in% panels && hasll){
    # Each subject's total log likelihood, summed over that subject's rows. The
    # row-wise version is in every other panel; this one aggregates within a
    # subject first, so a subject the model fits badly across many occasions
    # shows up as one point in the tail rather than being spread thin.
    #
    # Compared as ECDFs against a band of replicate datasets rather than as
    # densities. Each generated draw supplies one total per subject, so its ECDF
    # is a replicate of the observed ECDF at the same number of subjects -- the
    # band is then the right reference at any n, where a kernel density over a
    # few dozen subjects is not.
    ll <- dat[as.character(variable) == 'LogLik' & is.finite(obsValue) & is.finite(value)]
    # Both sides must total the same rows, or the totals are not comparable.
    ll <- ll[, if(all(is.finite(value)) && is.finite(obsValue[1L])) .SD, by = .(row)]
    gensub <- ll[, .(ll = sum(value)), by = .(sample, id)]
    obssub <- unique(ll[, .(row, id, obsValue)])[, .(ll = sum(obsValue)), by = id]
    if(nrow(obssub) > 2 && nrow(gensub) > 0){
      grid <- seq(min(c(gensub$ll, obssub$ll)), max(c(gensub$ll, obssub$ll)),
        length.out = 200)
      ecdfs <- vapply(split(gensub$ll, gensub$sample),
        function(v) vapply(grid, function(g0) mean(v <= g0), numeric(1)),
        numeric(length(grid)))
      if(is.null(dim(ecdfs))) ecdfs <- matrix(ecdfs, nrow = length(grid))
      band <- data.table(
        x = grid,
        lo = apply(ecdfs, 1, quantile, qlo, na.rm = TRUE),
        mid = apply(ecdfs, 1, median, na.rm = TRUE),
        hi = apply(ecdfs, 1, quantile, 1 - qlo, na.rm = TRUE),
        obsF = vapply(grid, function(g0) mean(obssub$ll <= g0), numeric(1)))
      outside <- round(100 * mean(band$obsF < band$lo | band$obsF > band$hi), 1)
      g <- ggplot(band, aes(x = x)) +
        geom_ribbon(aes(ymin = lo, ymax = hi, fill = 'Model replicates'), alpha = .3) +
        geom_line(aes(y = mid, colour = 'Model median'), linewidth = .7) +
        geom_line(aes(y = obsF, colour = 'Observed'), linewidth = .9) +
        geom_rug(data = obssub, aes(x = ll), inherit.aes = FALSE,
          colour = .ctPostPredCols[['Observed']], alpha = .6) +
        scale_fill_manual(name = '', values = c('Model replicates' = .ctPostPredCols[['Model']])) +
        scale_colour_manual(name = '', values = c(
          'Model median' = .ctPostPredCols[['Model']],
          'Observed' = .ctPostPredCols[['Observed']])) +
        coord_cartesian(ylim = c(0, 1)) +
        theme_bw() +
        labs(x = 'Total log likelihood per subject', y = 'Cumulative proportion of subjects',
          subtitle = paste0(nrow(obssub), ' subjects, ', ndraws,
            ' replicate datasets. Observed curve outside the band over ',
            outside, '% of the range.')) +
        theme(legend.position = 'bottom',
          plot.subtitle = element_text(size = rel(.8), colour = 'grey30'))
      gglist$SubjectLogLik <- cap(g, 'SubjectLogLik')
    }
  }

  if(plot){
    first <- TRUE
    for(nm in names(gglist)){
      if(wait && !first) readline('Press [return] for next plot.')
      first <- FALSE
      suppressWarnings(print(gglist[[nm]]))
    }
    return(invisible(gglist))
  }
  gglist
}

.ctPostPredResolvePanels <- function(panels){
  known <- unlist(.ctPostPredGroups, use.names = FALSE)
  if(identical(panels, 'all')) return(known)
  out <- unlist(lapply(panels, function(p)
    if(p %in% names(.ctPostPredGroups)) .ctPostPredGroups[[p]] else p), use.names = FALSE)
  bad <- setdiff(out, known)
  if(length(bad)) stop('unknown panels: ', paste(bad, collapse = ', '),
    '. Use "all", "structure", "calibration", or any of: ',
    paste(known, collapse = ', '), call. = FALSE)
  out
}

# Equal-count bins, falling back to a single bin when there are not enough
# distinct values to cut.
.ctPostPredBin <- function(x, nbins){
  br <- unique(quantile(x, seq(0, 1, length.out = nbins + 1), na.rm = TRUE))
  if(length(br) < 3) return(rep(1L, length(x)))
  as.integer(cut(x, breaks = br, include.lowest = TRUE))
}

# Wilson interval -- behaves at 0 and 1, where the normal approximation behind
# the geom_smooth ribbon this replaces ran outside [0,1].
.ctPostPredBinomCI <- function(k, n, level = .95){
  z <- qnorm(1 - (1 - level) / 2)
  p <- k / n
  d <- 1 + z^2 / n
  c1 <- (p + z^2 / (2 * n)) / d
  c2 <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(pmax(0, c1 - c2), pmin(1, c1 + c2))
}


#' Posterior predictive checks for a ctsem fit
#'
#' \code{ctPostPredict} and \code{ctStanPostPredict} are aliases for
#' \code{\link{ctPostPredPlots}}, which the two functions were merged into. They differ from
#' it only in printing the panels by default rather than returning them.
#'
#' @param fit A fitted ctsem model.
#' @param ... Passed to \code{\link{ctPostPredPlots}}.
#' @return A named list of ggplot objects, invisibly.
#' @aliases ctStanPostPredict
#' @export
#' @examples
#' \donttest{
#' ctPostPredict(ctstantestfit, wait = FALSE)
#' }
ctPostPredict <- function(fit, ...){
  args <- list(...)
  # Arguments the old ctPostPredict() took that the merged function cannot mean
  # the same thing by. Named rather than silently ignored -- an argument accepted
  # on one path and dropped on another is the trap this file has already had once.
  if('probs' %in% names(args)){
    warning('ctPostPredict(): `probs` is ignored. It never affected a plot -- it fed a ',
      'quantile array whose only use was to supply variable names. Use `interval` to set ',
      'the width of the predictive interval instead.', call. = FALSE)
    args$probs <- NULL
  }
  unknown <- setdiff(names(args), names(formals(ctPostPredPlots)))
  if(length(unknown)) stop('unused arguments: ', paste(unknown, collapse = ', '),
    '. See ?ctPostPredPlots for what this function now accepts.', call. = FALSE)
  if(is.null(args$plot)) args$plot <- TRUE
  do.call(ctPostPredPlots, c(list(fit = fit), args))
}

#' @export
ctStanPostPredict <- ctPostPredict
