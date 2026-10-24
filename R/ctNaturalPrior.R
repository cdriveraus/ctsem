# A flat prior on each parameter's natural scale, for the densities the
# uncertainty methods integrate over.
#
# A fit without priors (`priors = FALSE`) leaves its raw parameters flat. A
# raw parameter is the unbounded coordinate a transform maps to the model's
# value -- a standard deviation is roughly exp(raw) -- so flat on the raw scale
# is a prior proportional to 1 / sd on the sd, which piles up at zero: where
# the likelihood levels off as an sd goes to zero (a noise source the data can
# do without, a population with no individual differences) the posterior over
# the raw parameters has infinite mass and does not exist, and a sampler
# reports wherever its draws happened to reach. Flat on the natural scale
# instead is the change of variable `log |f'(x)|`, which falls off wherever `f`
# saturates towards a bound and is constant where it does not, so the bounded
# limits become proper and the unbounded ones stay as flat as they were. See
# review/UNCERTAINTY-checked-improper-2026-10-24.md.
#
# The terms are written here, symbolically, from each parameter's own
# transform, and evaluated by the engine as expressions beside its other prior
# terms (`_ctsem_natural_prior`): a user's transform -- a drift that may be
# positive, say -- gets its own Jacobian, constant if the transform is linear.

# `log |e|` for a derivative `e`, split over products, quotients and powers,
# with `log(inv_logit(z))` written as `-log1p_exp(-z)` and `log(exp(z))` as `z`,
# so that it stays finite far into a transform's saturation. `NULL` where a
# factor is not one of these and so might change sign.
.ctNaturalLogAbs <- function(e) {
  # Parentheses and a sign change nothing in `log |e|`.
  while (is.call(e) && length(e) == 2L &&
      (identical(e[[1]], as.name("(")) || identical(e[[1]], as.name("-")))) e <- e[[2]]
  if (is.numeric(e)) return(if (e == 0) NULL else log(abs(e)))
  if (!is.call(e)) return(NULL)
  f <- as.character(e[[1]])
  two <- function(op) {
    a <- .ctNaturalLogAbs(e[[2]]); b <- .ctNaturalLogAbs(e[[3]])
    if (is.null(a) || is.null(b)) NULL else call(op, a, b)
  }
  if (f == "*" && length(e) == 3L) return(two("+"))
  if (f == "/" && length(e) == 3L) return(two("-"))
  if (f == "exp" && length(e) == 2L) return(e[[2]])
  if (f == "inv_logit" && length(e) == 2L)
    return(call("-", call("log1p_exp", call("-", e[[2]]))))
  if (f == "^" && length(e) == 3L && is.numeric(e[[3]])) {
    a <- .ctNaturalLogAbs(e[[2]])
    return(if (is.null(a)) NULL else call("*", e[[3]], a))
  }
  NULL
}

# One parameter's terms from its transform `text` (reading `param[k]` and no
# other parameter): `log |f'|`, its derivative and its second derivative, as
# text reading `param[k]`; `NULL` where the transform is linear in it (a
# constant Jacobian changes nothing) or its derivative is not a product of
# factors whose sign cannot change.
.ctNaturalTerms <- function(text, k) {
  x <- as.name("x__natural")
  f <- str2lang(gsub(paste0("param\\[\\s*", k, "\\s*\\]"), "x__natural", text))
  slope <- function(e) .ctWithStableDerivRules(Deriv::Deriv(e, "x__natural",
    cache.exp = FALSE))
  d1 <- slope(f)
  if (!("x__natural" %in% all.vars(d1))) return(NULL)
  value <- .ctNaturalLogAbs(d1)
  if (is.null(value)) return(NULL)
  value <- Deriv::Simplify(value)
  d <- slope(value)
  dd <- slope(d)
  back <- function(e) gsub("x__natural", paste0("param[", k, "]"), deparse1(e,
    width.cutoff = 500L), fixed = TRUE)
  list(value = back(value), slope = back(d), curvature = back(dd))
}

# The natural-scale prior for a specification: one row per free raw parameter
# without a prior of its own whose transform is nonlinear in it alone --
# `index`, and the text of `value` (log |f'|), `slope` and `curvature` (its
# first two derivatives). A cell's transform where the parameter has a cell
# (on the augmented route a population sd is one, RAWPOPVAR); on the Laplace
# route a population sd or own variance is the level's, `log1p_exp(2x - 1)`
# times its scale, which the engine applies itself (`_level_popcov`).
# Parameters with no cell -- predictor effects, loadings, correlations, which
# the correlation bound already holds -- keep the raw scale.
.ctJuliaNaturalPrior <- function(spec) {
  npar <- .ctBackendNpar(spec)
  withprior <- if (is.null(spec$priors)) integer() else as.integer(spec$priors$index)
  table <- as.data.frame(spec$parameter_table, stringsAsFactors = FALSE)
  texts <- character(npar)
  for (k in setdiff(seq_len(npar), withprior)) {
    rows <- which(!is.na(table$parnumber) & as.integer(table$parnumber) == k)
    for (r in rows) {
      tf <- as.character(table$transform[r])
      if (is.na(tf) || !nzchar(tf)) next
      reads <- unique(as.integer(regmatches(tf, gregexpr("(?<=param\\[)\\s*[0-9]+",
        tf, perl = TRUE))[[1]]))
      if (identical(reads, k)) { texts[k] <- tf; break }
    }
  }
  for (lv in spec$laplace$levels) {
    for (field in c("sd_index", "unique_index")) {
      idx <- as.integer(lv[[field]])
      for (j in seq_along(idx)) {
        k <- idx[j]
        if (is.na(k) || k < 1L || k > npar || k %in% withprior) next
        texts[k] <- sprintf("1e-10 + %.17g * log1p_exp(2 * param[%d] - 1)",
          as.numeric(lv$sd_scale[j]), k)
      }
    }
  }
  out <- lapply(which(nzchar(texts)), function(k) {
    terms <- tryCatch(.ctNaturalTerms(texts[k], k), error = function(e) NULL)
    if (is.null(terms)) NULL else data.frame(index = k, value = terms$value,
      slope = terms$slope, curvature = terms$curvature, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) data.frame(index = integer(), value = character(),
    slope = character(), curvature = character()) else out
}

# The natural-scale prior's curvature at `x`, one value per row of `natural`,
# for the uncertainty methods' Hessian: the terms are separable, so this is
# the diagonal they add.
.ctNaturalCurvature <- function(natural, x) {
  env <- list(param = as.numeric(x), inv_logit = stats::plogis,
    log1p_exp = function(z) ifelse(z > 30, z, log1p(exp(z))))
  vapply(natural$curvature, function(text) as.numeric(eval(str2lang(text), env)),
    numeric(1), USE.NAMES = FALSE)
}
