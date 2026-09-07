# What a time-independent predictor effect can say.
#
# `<TI>_effect` was a logical: an effect on this parameter exists, or it does
# not. That is all a fit needs, because an effect that exists is estimated. It
# is not all a *specification* needs: generating data requires a size, and
# constraining two parameters to share an effect requires a name, and neither
# fits in a logical.
#
# So the column holds a character, written in the spec's own separator block:
#
#   'mypar||||TI1=4.3, TI2=myti2effect, TI3'
#
# and reading left to right that is: a fixed effect of 4.3 from TI1, a free
# effect from TI2 called `myti2effect` -- nameable, so another parameter can be
# constrained to the same one -- and a free effect from TI3 named automatically,
# which is what a bare name has always meant.
#
# Four states, then, and every consumer wants a different one of them, which is
# why they go through these three functions rather than testing the column
# directly:
#
#   'FALSE' / ''   no effect
#   'TRUE'         free, named automatically
#   '4.3'          fixed at that value
#   'myeffect'     free, named
#
# Character rather than numeric because "free" and "fixed at 1" are different
# things and both have to be sayable. As a number they are the same number.

#' @keywords internal
.ctTipredEffectSpec <- function(x) {
  x <- as.character(unlist(x))
  x[is.na(x)] <- "FALSE"
  trimws(x)
}

# Is there an effect at all? The question `<TI>_effect` used to answer, and the
# one nearly every consumer is actually asking.
#' @keywords internal
.ctTipredEffectActive <- function(x) {
  spec <- .ctTipredEffectSpec(x)
  active <- !spec %in% c("FALSE", "F", "", "NA", "0")
  # A numeric zero is no effect, however it was written.
  numeric <- suppressWarnings(as.numeric(spec))
  active & !(is.finite(numeric) & numeric == 0)
}

# Is it fixed, and to what? NA for anything free.
#' @keywords internal
.ctTipredEffectValue <- function(x) {
  spec <- .ctTipredEffectSpec(x)
  value <- suppressWarnings(as.numeric(spec))
  # 'TRUE' parses to NA as a number, which is what we want; a label does too.
  value[!.ctTipredEffectActive(x)] <- NA_real_
  value
}

# The name a free effect was given, or NA where it was not named or is not free.
#' @keywords internal
.ctTipredEffectLabel <- function(x) {
  spec <- .ctTipredEffectSpec(x)
  label <- rep(NA_character_, length(spec))
  named <- .ctTipredEffectActive(x) & is.na(.ctTipredEffectValue(x)) &
    !spec %in% c("TRUE", "T")
  label[named] <- spec[named]
  label
}

# Whether an effect is free, i.e. estimated rather than fixed.
#' @keywords internal
.ctTipredEffectFree <- function(x) {
  .ctTipredEffectActive(x) & is.na(.ctTipredEffectValue(x))
}

# Which effects in a model are fixed to a value, named for a message.
#' @keywords internal
.ctTipredFixedEffects <- function(model) {
  names <- model$TIpredNames
  pars <- model$pars
  if (!length(names) || is.null(pars)) return(character())
  out <- character()
  for (predictor in names) {
    column <- paste0(predictor, "_effect")
    if (is.null(pars[[column]])) next
    value <- .ctTipredEffectValue(pars[[column]])
    fixed <- which(is.finite(value))
    for (i in fixed) {
      label <- if (is.na(pars$param[i])) paste0(pars$matrix[i], "[",
        pars$row[i], ",", pars$col[i], "]") else as.character(pars$param[i])
      out <- c(out, paste0(label, "=", predictor, "=", format(value[i])))
    }
  }
  out
}
