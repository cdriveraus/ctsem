# Keeping a plot readable when a model quantity is enormously wider than the
# data it is drawn against.
# ------------------------------------------------------------------------------
#
# A fitted model can produce quantities with no useful upper bound. A count
# whose linear predictor carries a latent sd of two has a posterior predictive
# spanning four orders of magnitude; a weakly identified parameter has an
# uncertainty band far wider than anything observed. One value in such a tail
# sets the axis for the whole panel, and the data -- with the part of the model
# that can actually be compared against it -- ends up in a sliver at the bottom.
#
# That is not a broken model and the plot is not wrong. It is a model whose
# predictive really is that wide, which is a finding, and the panel has to stay
# readable enough to report it. So the view is cut back and what fell outside is
# said in words, which turns an unreadable axis into a number.
#
# Three rules, and they are what make this safe to apply widely:
#
#   1. The reference values -- the data, or whatever the model is being compared
#      against -- are never clipped. The comparison is the point, so data
#      outside the view would misrepresent it rather than compress it.
#   2. It engages only when it changes the picture by a lot. `slack` is how many
#      times wider the untrimmed range must be before the trimmed one is used,
#      so an ordinary fit, whose model tails sit a little outside its data, is
#      left exactly as it was. This is a guard for a few extreme cases, not a
#      new default.
#   3. It reports itself. A caller that cuts a view without saying so has
#      quietly changed what the reader is looking at, so the caption names the
#      variable, how much is off the axis and how far it reaches. That is the
#      whole of what a reader needs: the tail becomes a number rather than an
#      unreadable axis, and nothing is hidden by cutting the view.
#
# Callers cut the view per facet, by filtering or clamping their own data,
# rather than by setting a scale transformation: one transformation applies to
# a whole plot, while `scales = 'free'` lets each facet keep its own range, and
# a variable needing compression usually sits beside one that does not.
#
# Internal throughout. There is no switch for it and it is not part of the API:
# it engages only where the untrimmed panel was unreadable, and it says in the
# caption what it did, so there is nothing for a caller to have to manage.


#' @noRd
#' View limits for one variable, and how much of `values` falls outside them.
#'
#' `values` are the model quantities that may be extreme; `keep` are the
#' reference values that must stay visible, and may be empty. Returns NULL
#' limits when the guard is off, when there is nothing to cut, or when cutting
#' would not materially change the range -- callers treat that as "plot as
#' before" and do nothing.
.ctPlotView <- function(values, keep = numeric(0), trim = c(.005, .995),
  slack = 5){
  none <- list(limits = NULL, nout = 0L, max = NA_real_, n = 0L)
  v <- values[is.finite(values)]
  k <- keep[is.finite(keep)]
  if(length(v) < 2) return(none)
  full <- range(c(v, k))
  cut <- range(c(stats::quantile(v, trim, names = FALSE, na.rm = TRUE), k))
  if(!all(is.finite(cut)) || !all(is.finite(full))) return(none)
  if(diff(cut) <= 0 || diff(full) <= slack * diff(cut)) return(none)
  out <- v < cut[1] | v > cut[2]
  # The furthest value *that was cut*, not the largest in absolute value: with
  # a view cut only at one end, the biggest number can be one still on screen,
  # and reporting that as the reach would describe the wrong thing.
  list(limits = cut, nout = sum(out), n = length(v),
    max = if(any(out)) v[out][which.max(abs(v[out]))] else NA_real_)
}


#' @noRd
#' `.ctPlotView` over each group of a long table, returning one row per group
#' the guard engaged on, and NULL when it engaged on none.
#'
#' `valuecols` may name more than one column when the extreme quantity is a
#' band rather than a single series -- they are pooled, so the view is one
#' range for the group rather than a different one per edge. `refcol` names
#' what has to stay visible: the observed data where there is any, and
#' otherwise the central estimate, since a band far wider than its own median
#' is exactly the case where the median is what the reader needs to see.
.ctPlotViews <- function(dt, valuecols = 'value', refcol = 'obsValue',
  by = 'variable'){
  if(!nrow(dt) || !by %in% names(dt)) return(NULL)
  valuecols <- intersect(valuecols, names(dt))
  if(!length(valuecols)) return(NULL)
  g <- as.character(dt[[by]])
  rows <- lapply(unique(g), function(v){
    d <- dt[g == v]
    vals <- unlist(lapply(valuecols, function(nm) d[[nm]]), use.names = FALSE)
    ref <- if(!is.null(refcol) && refcol %in% names(d)) d[[refcol]] else numeric(0)
    vw <- .ctPlotView(vals, ref)
    if(is.null(vw$limits)) return(NULL)
    data.table::data.table(variable = v, lo = vw$limits[1], hi = vw$limits[2],
      nout = vw$nout, maxout = vw$max, n = vw$n)
  })
  rows <- data.table::rbindlist(rows)
  if(!nrow(rows)) NULL else rows
}


#' @noRd
#' The sentence a panel adds to its caption when its view was cut: which
#' variable, how much of the model is off the axis, how far it reaches, and how
#' to turn the cut off.
.ctPlotViewNote <- function(views){
  if(is.null(views) || !nrow(views)) return(NULL)
  num <- function(x) format(signif(x, 3), big.mark = ',', scientific = FALSE,
    trim = TRUE)
  bits <- sprintf('%s reaches %s (%s of model values are off the axis)',
    views$variable, num(views$maxout),
    paste0(signif(100 * views$nout / views$n, 2), '%'))
  # Deliberately says "what is being compared against" rather than "the data":
  # the same note is used where there is no data and the reference is the
  # central estimate, and a caption that named data there would be wrong.
  paste0(' Axis cut back to the 0.5-99.5% range of the model, and to whatever ',
    'it is compared against, because the rest is far wider than that: ',
    paste(bits, collapse = '; '), '.')
}


#' @noRd
#' Clamp the named columns of `dt` into each variable's view. For a band or a
#' line, where dropping the row would take a reference value with it and the
#' honest picture is a band drawn running off the edge.
.ctPlotViewClamp <- function(dt, views, cols, by = 'variable'){
  if(is.null(views) || !nrow(views) || !by %in% names(dt)) return(dt)
  cols <- intersect(cols, names(dt))
  if(!length(cols)) return(dt)
  out <- data.table::copy(dt)
  for(i in seq_len(nrow(views))){
    j <- which(as.character(out[[by]]) == views$variable[i])
    if(!length(j)) next
    for(nm in cols) data.table::set(out, i = j, j = nm,
      value = pmin(pmax(out[[nm]][j], views$lo[i]), views$hi[i]))
  }
  out
}


#' @noRd
#' Drop rows of `dt` whose `valuecol` lies outside its variable's view. For a
#' density, which has to be estimated on the values it is drawn over -- clamping
#' would pile the tail onto the boundary as a spike the model does not have.
.ctPlotViewFilter <- function(dt, views, valuecol = 'value', by = 'variable'){
  if(is.null(views) || !nrow(views) || !valuecol %in% names(dt) ||
      !by %in% names(dt)) return(dt)
  keep <- rep(TRUE, nrow(dt))
  v <- dt[[valuecol]]
  for(i in seq_len(nrow(views))){
    j <- which(as.character(dt[[by]]) == views$variable[i] & is.finite(v))
    if(!length(j)) next
    keep[j] <- v[j] >= views$lo[i] & v[j] <= views$hi[i]
  }
  dt[keep]
}
