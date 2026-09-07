# What a fit can honestly say about its own credibility.
#
# An optimiser can stop somewhere the data does not determine, and the fit then
# reports estimates and standard errors that look exactly like any other fit's.
# Observed on a six-subject model: population standard deviations driven to
# 1e-8, process parameters at plus or minus a hundred, `converged = TRUE`.
#
# The temptation is to call such a fit absurd and refuse it. That is the wrong
# move, because "absurd" is a judgement about the model and the data that this
# code is in no position to make -- a population standard deviation of zero is a
# legitimate finding (no detectable individual differences), an extreme raw
# value can be a perfectly ordinary value on a transformed scale, and a
# likelihood that is flat in some direction may still be exactly what the
# researcher wants to report.
#
# What can be said without assuming anything is narrower and more useful: which
# directions the data does not identify, and therefore which reported intervals
# do not mean what a reader will take them to mean. That is a property of the
# information matrix, checkable, and it does not require deciding whether the
# answer is sensible.

#' @keywords internal
.ctBackendIdentifiability <- function(hessian, parnames = NULL, rtol = 1e-8,
  loading = 0.25, fit = NULL, at = NULL, metric = NULL, vectors = FALSE) {
  empty <- list(nweak = 0L, condition = NA_real_, directions = list(),
    parameters = character())
  if (is.null(hessian)) return(empty)
  hessian <- as.matrix(hessian)
  if (!all(is.finite(hessian)) || nrow(hessian) != ncol(hessian)) return(empty)
  n <- nrow(hessian)
  if (is.null(parnames) || length(parnames) != n) parnames <- paste0("par", seq_len(n))

  # The information matrix, symmetrised. Eigenvalues rather than a determinant
  # or a condition number alone: which direction is flat is the actionable part,
  # and only the eigenvectors say that.
  information <- -(hessian + t(hessian)) / 2
  decomposition <- try(eigen(information, symmetric = TRUE), silent = TRUE)
  if (inherits(decomposition, "try-error")) return(empty)
  values <- decomposition$values
  scale <- max(abs(values))
  if (!is.finite(scale) || scale <= 0) return(empty)

  # A direction counts as unidentified when its curvature is negligible against
  # the sharpest direction, or when it is negative -- a negative eigenvalue of
  # the information matrix means the optimiser stopped somewhere that is not a
  # maximum in that direction at all.
  weak <- which(values <= rtol * scale)
  directions <- lapply(weak, function(k) {
    loadings <- decomposition$vectors[, k]
    involved <- order(-abs(loadings))
    involved <- involved[abs(loadings[involved]) >= loading]
    if (!length(involved)) involved <- which.max(abs(loadings))
    # The whole eigenvector, not only the loadings above the threshold. Two
    # things downstream need it: the partial-identification check below, which
    # asks whether a functional of the parameters changes along this direction,
    # and the aggregation in `ctIdentify()`, which asks how much of a
    # coordinate lies in the flat *subspace* rather than on one of its axes.
    list(eigenvalue = values[k], relative = values[k] / scale,
      parameters = parnames[involved], loadings = loadings[involved],
      vector = loadings)
  })
  # Which of those directions are a random-effect block trading its scale off
  # against its correlations -- partially rather than completely unidentified.
  # Only attempted when a caller supplies the model and the point, because it
  # takes engine calls; without them every direction reads as complete
  # non-identification, which is what this said before.
  directions <- .ctIdentifyClassify(fit, at, directions, metric = metric)
  # The eigenvectors go out only when a caller asked for them. A fit stores
  # this report, and one length-`npar` vector per flat direction is
  # `nweak * npar` doubles for something nothing downstream of a fit reads --
  # 12 MB on a 1490-parameter model with a thousand flat directions. The
  # classification above and `ctIdentify()`'s subspace aggregation are the two
  # readers, and only the second outlives this call.
  if (!isTRUE(vectors)) {
    directions <- lapply(directions, function(d) { d$vector <- NULL; d })
  }
  list(
    nweak = length(weak),
    condition = scale / max(min(values[values > 0], na.rm = TRUE), .Machine$double.xmin),
    # Negative *against the scale of the matrix*, not against zero. An
    # eigenvalue of -3e-16 where the largest is 9e5 is what a symmetric
    # eigendecomposition does to a direction whose true curvature is zero; it
    # is the flat direction `nweak` already counts, not a saddle. Counting it
    # as negative fired "this is not a maximum" on ordinary fits sitting at
    # their optimum with a gradient of 5e-10, which is how a warning worth
    # reading gets ignored.
    negative = sum(values < -rtol * scale),
    directions = directions,
    parameters = unique(unlist(lapply(directions, `[[`, "parameters"))))
}

# Partial against complete non-identification, and why the difference is worth
# the engine calls it costs.
#
# A random effect on a *variance* cell -- DIFFUSION, MANIFESTVAR -- under
# `intoverpop='augmented'` is a latent state whose value enters only the
# predicted covariance. The observation mean function's Jacobian with respect
# to it is zero, so the Kalman update can never move it, and it is learned
# about only through its correlation with states the filter can update. The
# consequence, measured: the population sd and the correlations trade off along
# a ridge whose *cross-covariances* are pinned to five significant figures
# while the sd runs over a factor of 16 and the log likelihood moves by 1.3e-06.
#
# So the data does determine something here, and the advice that fits complete
# non-identification -- fix one of the set, or remove it -- throws that
# something away. Telling the two cases apart is computable rather than
# guessable: if every coordinate a flat direction loads on is a population
# scale or correlation of one random-effect block, and the block's
# cross-covariance functionals do not change along the direction, then the
# covariances are what the data determines and their decomposition into scales
# and correlations is what it does not.
#
# The one functional that legitimately *does* change along such a direction is
# the implicated effect's own variance -- a variance is not a covariance
# between two effects, and nothing in this mechanism identifies it. Any other
# entry moving is a different problem, and then this reports nothing and the
# complete-non-identification advice stands.

# Which raw coordinates make up each level's population covariance block, and
# how to materialise that covariance.
#
# Both routes are described the same way because the identification question is
# the same on both: the augmented route holds the scales and correlations in
# `spec$random_effects` (they are T0VAR cells of the augmented model), the
# Laplace route in `spec$laplace$levels`. A block is one level's scales and
# correlations together, because a scale is identified or not *jointly with the
# correlations it multiplies*.
#' @keywords internal
.ctIdentifyBlocks <- function(spec) {
  if (is.null(spec)) return(list())
  blocks <- list()
  effects <- spec$random_effects
  if (!is.null(effects) && length(effects) && nrow(effects)) {
    sds <- effects[effects$type %in% "sd", , drop = FALSE]
    cors <- effects[effects$type %in% "correlation", , drop = FALSE]
    if (nrow(sds)) blocks[[length(blocks) + 1L]] <- list(
      route = "augmented", level = 1L, name = "subject",
      sd_index = as.integer(sds$parameter), cor_index = as.integer(cors$parameter),
      param = as.character(sds$param),
      # The carrier state each scale belongs to, which is where it sits in the
      # augmented T0VAR and therefore in the population covariance.
      state = as.integer(sds$row))
  }
  laplace <- spec$laplace
  if (!is.null(laplace) && length(laplace$levels)) {
    for (l in seq_along(laplace$levels)) {
      level <- laplace$levels[[l]]
      if (!length(level$sd_index)) next
      blocks[[length(blocks) + 1L]] <- list(
        route = "laplace", level = as.integer(l),
        name = as.character(.ctJuliaOr(level$name, l))[1L],
        sd_index = as.integer(level$sd_index),
        cor_index = as.integer(level$cor_index),
        param = as.character(level$param),
        state = seq_along(level$sd_index))
    }
  }
  blocks
}

# The population covariance of one block at one raw vector.
#' @keywords internal
.ctIdentifyPopcov <- function(fit, block, values) {
  spec <- .ctBackendSpec(fit)
  if (identical(block$route, "laplace")) {
    module <- .ctJuliaModule(spec$project)
    objective <- .ctJuliaObjective(fit)
    out <- lapply(seq_len(ncol(values)), function(column) {
      value <- try(.ctBackendJuliaValue(module$ctsem_laplace_popcov(objective,
        .ctJuliaNumericVector(as.numeric(values[, column])),
        as.integer(block$level))), silent = TRUE)
      if (inherits(value, "try-error")) return(NULL)
      as.numeric(as.matrix(value))
    })
    if (any(vapply(out, is.null, logical(1)))) return(NULL)
    return(matrix(unlist(out), ncol = ncol(values)))
  }
  # The augmented route's population covariance is the carrier block of the
  # augmented model's T0cov, and `rows` keeps everything else off the bridge.
  layout <- try(.ctBackendSummaryLayout(fit), silent = TRUE)
  if (inherits(layout, "try-error")) return(NULL)
  slot <- which(layout$matrix %in% "T0cov")
  if (!length(slot)) return(NULL)
  n <- layout$nrow[slot[1L]]
  states <- block$state
  if (!length(states) || any(!is.finite(states)) || any(states > n)) return(NULL)
  grid <- expand.grid(i = seq_along(states), j = seq_along(states))
  rows <- layout$offset[slot[1L]] +
    (states[grid$j] - 1L) * n + states[grid$i]
  flat <- try(.ctBackendParMatricesFlat(fit, values, rows = rows), silent = TRUE)
  if (inherits(flat, "try-error")) return(NULL)
  flat
}

# Central differences of every entry of the block's population covariance, over
# the block's own coordinates only.
#
# Computed once per block and point rather than once per direction, because it
# does not depend on the direction and a rank-limited model can have dozens.
#' @keywords internal
.ctIdentifyPopcovGradient <- function(fit, at, block, h = 1e-5) {
  index <- c(block$sd_index, block$cor_index)
  index <- index[is.finite(index) & index >= 1L & index <= length(at)]
  if (!length(index)) return(NULL)
  values <- matrix(as.numeric(at), nrow = length(at),
    ncol = 2L * length(index) + 1L)
  for (position in seq_along(index)) {
    values[index[position], 2L * position] <-
      values[index[position], 2L * position] + h
    values[index[position], 2L * position + 1L] <-
      values[index[position], 2L * position + 1L] - h
  }
  covariance <- .ctIdentifyPopcov(fit, block, values)
  if (is.null(covariance) || !all(is.finite(covariance))) return(NULL)
  k <- length(block$sd_index)
  if (nrow(covariance) != k * k) return(NULL)
  gradient <- matrix(0, nrow = k * k, ncol = length(index))
  for (position in seq_along(index)) {
    gradient[, position] <- (covariance[, 2L * position] -
      covariance[, 2L * position + 1L]) / (2 * h)
  }
  list(index = index, gradient = gradient, size = sqrt(rowSums(gradient^2)),
    entries = expand.grid(i = seq_len(k), j = seq_len(k)),
    value = covariance[, 1L])
}

# How far each entry of the block's population covariance is from constant
# along `direction`: zero means the direction preserves it.
#
# The direction has been checked to lie in the block already, so its components
# elsewhere are zero and cannot contribute to the derivative along it.
# Restricting the gradient to the block can therefore only make the reported
# cosine larger, which is the safe side for a check that has to be convinced
# before it will say the covariance is determined.
#' @keywords internal
.ctIdentifyPopcovCosines <- function(measured, direction) {
  step <- direction[measured$index]
  norm <- sqrt(sum(step^2))
  if (!is.finite(norm) || norm <= 0) return(NULL)
  along <- abs(as.numeric(measured$gradient %*% (step / norm)))
  ifelse(measured$size > 0, along / measured$size, NA_real_)
}

# Tag each flat direction that is one random-effect block's scale/correlation
# ridge rather than a parameter the data says nothing about.
#' @keywords internal
.ctIdentifyClassify <- function(fit, at, directions, metric = NULL,
  mass = 0.99, cosine = 1e-4) {
  if (is.null(fit) || is.null(at) || !length(directions)) return(directions)
  blocks <- try(.ctIdentifyBlocks(.ctBackendSpec(fit)), silent = TRUE)
  if (inherits(blocks, "try-error") || !length(blocks)) return(directions)
  at <- as.numeric(at)
  # One entry per block, filled the first time a direction needs it and NA
  # once it is known not to be obtainable, so a failed probe is not retried
  # for every remaining direction.
  probes <- vector("list", length(blocks))
  for (k in seq_along(directions)) {
    v <- directions[[k]]$vector
    if (length(v) != length(at)) next
    total <- sum(v^2)
    if (!is.finite(total) || total <= 0) next
    for (b in seq_along(blocks)) {
      block <- blocks[[b]]
      index <- c(block$sd_index, block$cor_index)
      index <- index[is.finite(index) & index >= 1L & index <= length(v)]
      if (!length(index) || length(block$sd_index) < 2L) next
      # Everything the direction loads on has to be in this block. A single
      # random effect (`sd_index` shorter than two) is skipped above: with no
      # partner there is no covariance for the data to determine, so its scale
      # being flat is complete non-identification and reads as such.
      if (sum(v[index]^2) < mass * total) next
      if (is.null(probes[[b]])) {
        probe <- .ctIdentifyPopcovGradient(fit, at, block)
        probes[[b]] <- if (is.null(probe)) NA else probe
      }
      if (!is.list(probes[[b]])) next
      measured <- probes[[b]]
      # `metric` is the scaling the information matrix was put in before its
      # eigenvectors were taken (`.ctIdentifyInformation`); the covariance
      # functionals are differentiated in raw coordinates, so the direction has
      # to come back to them before the two are compared.
      raw <- if (is.null(metric)) v else v / as.numeric(metric)
      measured$cosine <- .ctIdentifyPopcovCosines(measured, raw)
      if (is.null(measured$cosine)) next
      i <- measured$entries$i
      j <- measured$entries$j
      # Which sds are implicated is read off what actually moves rather than
      # off the eigenvector's loadings. At one fitted estimate the ridge
      # direction was 0.99 correlation and 0.16 scale in raw coordinates, so a
      # loading threshold of 0.25 said the scale was not involved while the
      # scale was exactly what was running -- the covariance it preserved says
      # so unambiguously and a threshold on the basis does not.
      present <- measured$size > 0 & is.finite(measured$cosine)
      moving <- present & measured$cosine >= cosine
      # A cross-covariance that moves along the direction is not determined
      # either, and then this is complete non-identification of the block, not
      # a scale trading against its correlations.
      if (any(moving & i != j)) next
      implicated <- sort(unique(i[moving]))
      if (!length(implicated)) next
      # Something has to be determined, or there is nothing to preserve and
      # nothing to say: a block with no cross-covariance to hold constant
      # tells us nothing about this direction.
      if (!any(present & i != j)) next
      directions[[k]]$partial <- list(route = block$route, level = block$level,
        block = block$name, parameters = block$param[implicated],
        partners = block$param[setdiff(seq_along(block$param), implicated)],
        cosine = max(measured$cosine[present & i != j]))
      break
    }
  }
  directions
}

# Split a set of flat directions into the parameters whose covariances the data
# still determines and the parameters it says nothing about.
#
# Two vocabularies meet here and they are not the same. `$partial` and
# `$partners` are *model* parameter names -- `diff_eta1` -- because that is
# what the sentence about population sds is about. `$covered` and
# `$structural` are *raw coordinate* names -- `popsd_diff_eta1` --  because
# that is what the eigenvectors are indexed by and what "parameters involved"
# has always listed. Comparing one against the other silently produced an
# empty intersection and the wrong message.
#' @keywords internal
.ctIdentifyPartition <- function(directions) {
  empty <- list(partial = character(), covered = character(),
    structural = character(), partners = character())
  if (!length(directions)) return(empty)
  partial <- character(); structural <- character()
  partners <- character(); covered <- character()
  for (direction in directions) {
    if (!is.null(direction$partial)) {
      partial <- c(partial, direction$partial$parameters)
      partners <- c(partners, direction$partial$partners)
      covered <- c(covered, direction$parameters)
    } else {
      structural <- c(structural, direction$parameters)
    }
  }
  covered <- unique(covered)
  list(partial = unique(partial), covered = covered,
    # A coordinate on a partially identified ridge and on a genuinely flat
    # direction elsewhere is reported under the stronger of the two.
    structural = setdiff(unique(structural), covered),
    # One direction's partner can be another's implicated scale; naming it in
    # both halves of the same sentence would be worse than dropping it here.
    partners = setdiff(unique(partners), unique(partial)))
}

# The wording, in one place, because it is said twice: before a fit by
# `print.ctIdentify` and after one by `.ctBackendIdentifyWarn`. One paragraph
# per element, so a caller can wrap it or run it together.
#' @keywords internal
.ctIdentifyAdvice <- function(partition) {
  lines <- character()
  if (length(partition$partial)) {
    many <- length(partition$partial) > 1L
    lines <- c(lines, paste0(
      "The population ", if (many) "sds of " else "sd of ",
      paste(partition$partial, collapse = ", "),
      if (many) " are" else " is",
      " not separately identified from ", if (many) "their" else "its",
      " correlations with ",
      if (length(partition$partners)) paste0(
        "the other individually varying parameters (",
        paste(partition$partners, collapse = ", "), ")") else "each other",
      ": only the covariances they generate are, and those the data ",
      "determines. Reported sds and correlations for these parameters trade ",
      "off along a ridge and will not repeat between runs. ",
      "intoverpop='laplace' identifies the ", if (many) "sds" else "sd",
      " separately, because there each subject's own random effect enters ",
      "that subject's likelihood."))
  }
  if (length(partition$structural)) {
    lines <- c(lines, paste0("Parameters involved: ",
      paste(partition$structural, collapse = ", "),
      ". These are not estimable from this data as the model stands. Fix one ",
      "of each set to a value, or remove it."))
  }
  lines
}

# Is each reported interval as wide as the curvature at the estimate supports?
#
# Two standard errors can be computed for a parameter from the same Hessian.
# The *conditional* one, `1 / sqrt(information[i, i])`, is what the curvature in
# that one coordinate supports with every other parameter held. The *marginal*
# one, `sqrt(cov[i, i])`, is what gets reported, and it is the conditional one
# divided by `sqrt(1 - R^2)`, where `R^2` is how well the other parameters
# reproduce this one in the information metric. So their ratio is a pure number
# saying how much of the reported width comes from the data and how much from
# the parameter not being separable from the rest: a ratio of 10 is `R^2` of
# .99, and 100 is .9999.
#
# Why it earns its place. A benchmark fit reached the right optimum -- the same
# log likelihood to eight decimal places as its twin, and Stan's -- and reported
# a drift interval a hundred times too wide, with a point estimate that had
# wandered with it. Nothing else about the fit looked wrong; it was caught only
# because two runs on identical data could be compared, and a user gets one run.
# This ratio separates the two cleanly: 1.0 to 1.5 across every parameter of the
# healthy fits measured here, against 1e4 on the parameter that had gone.
#
# Cheap: one diagonal and one square root, no extra evaluation of anything.
# Computed for julia fits, where the exact Hessian is already on the fit;
# nothing prevents the stan path from using it, and `ctReport()` is where that
# would show.
#
# The other half of the question: an interval far *narrower* than the truth.
#
# The ratio above only ever grows, so it cannot see the opposite failure, and
# that one is worse because it reads as a result rather than as a problem.
# `ctOptimCovFromHessian()` projects the no-curvature directions out of the
# information matrix before inverting it (see `.ctOptimIdentifiedInverse()`),
# which is the right thing to do with a direction whose variance is infinite --
# but a coordinate lying mostly *along* such a direction then inherits almost
# none of the variance that is left, and is reported with a tight interval and
# a large z instead of as undetermined. Measured on a one-latent model with
# `indvarying` on DRIFT and DIFFUSION under `intoverpop='augmented'`: the raw
# correlation between the two random effects has a dead flat profile -- the log
# likelihood is bit-identical at r = 0.597, 0.750, 0.958 and 0.998, the
# population sd compensating to hold their product fixed -- and `summary()`
# reported `mean 0.597, sd 0.009, z 65.3`. The spurious precision *grows* with
# the sample: z was 21.2 at 50 occasions and 65.3 at 200, so more data buys
# more confidence in a number the likelihood does not distinguish at all.
# Before the projection went in the same coordinate got a fabricated variance
# of 1e8; the projection fixed the blow-up and left this in its place.
#
# What says so is the coordinate's share of the projected-out subspace,
# `sum over flat k of V[i, k]^2` -- the squared length of `e_i`'s projection
# onto the null space, between 0 and 1, and unlike any single eigenvector's
# loading it does not depend on the arbitrary basis the eigendecomposition
# returns inside that subspace. Any of it that is not rounding means the
# asymptotic variance of that coordinate is infinite, whatever the projected
# covariance reports. `rtol` is therefore `.ctOptimIdentifiedInverse()`'s 1e-12
# rather than `.ctBackendIdentifiability()`'s 1e-8, so the two branches
# partition rather than overlap: a direction flatter than 1e-12 was dropped and
# lands here, one between the two tolerances was inverted into an enormous
# variance and lands in the ratio above. `nullmass` sits well above the
# 1e-16-ish leakage a well separated eigenvalue produces (measured: below
# 1e-10 on every identified coordinate of the model above) and well below the
# share a coordinate genuinely on the ridge carries, which was 0.36 for the
# smaller half of a two-coordinate ridge and 1.0 where the flat direction was
# a coordinate axis.
#' @keywords internal
.ctBackendIntervalCheck <- function(hessian, se, parnames = NULL,
  threshold = 100, rtol = 1e-12, nullmass = 1e-3) {
  empty <- list(threshold = threshold, nflagged = 0L, parameters = character(),
    nullmass = nullmass, nunidentified = 0L, unidentified = character(),
    table = data.frame(param = character(), se = numeric(),
      curvature_se = numeric(), ratio = numeric(), nullmass = numeric(),
      stringsAsFactors = FALSE))
  if (is.null(hessian) || is.null(se)) return(empty)
  hessian <- as.matrix(hessian)
  se <- as.numeric(se)
  if (nrow(hessian) != ncol(hessian) || nrow(hessian) != length(se)) return(empty)
  if (!all(is.finite(hessian))) return(empty)
  n <- length(se)
  if (is.null(parnames) || length(parnames) != n) parnames <- paste0("par", seq_len(n))
  information <- -(hessian + t(hessian)) / 2
  diagonal <- diag(information)
  # A non-positive diagonal is not a wider interval, it is no curvature at all;
  # `.ctBackendIdentifiability()` is what reports that, so it is left NA here
  # rather than counted as a ratio of infinity and reported twice.
  curvature <- ifelse(diagonal > 0, 1 / sqrt(diagonal), NA_real_)
  ratio <- se / curvature
  mass <- .ctBackendNullMass(information, rtol = rtol)
  table <- data.frame(param = as.character(parnames), se = se,
    curvature_se = curvature, ratio = ratio, nullmass = mass,
    stringsAsFactors = FALSE)
  flagged <- which(is.finite(ratio) & ratio > threshold)
  undetermined <- which(is.finite(mass) & mass >= nullmass)
  list(threshold = threshold, nflagged = length(flagged),
    parameters = as.character(parnames[flagged]),
    # The threshold goes out with the column, so a caller reading `$table` does
    # not have to know the default to read it.
    nullmass = nullmass, nunidentified = length(undetermined),
    unidentified = as.character(parnames[undetermined]),
    # Undetermined coordinates first, then the widest ratios. Ordering the whole
    # table by ratio put them at the bottom, which is where a reader stops
    # looking, and their ratio is small precisely because their interval
    # collapsed.
    table = table[order(-as.integer(is.finite(mass) & mass >= nullmass),
      -ifelse(is.finite(mass), mass, -Inf),
      -ifelse(is.finite(ratio), ratio, -Inf)), , drop = FALSE])
}

# How much of each coordinate lies in the null space of an information matrix.
#
# Zero for every coordinate when nothing is flat, and NA when the
# decomposition cannot be had -- never silently zero, because "no flat
# directions" and "could not tell" are opposite findings and one of them is a
# clean bill of health.
#
# The eigenvalues are taken first and the eigenvectors only if one of them is
# flat, because the usual answer is "nothing is flat" and a values-only
# decomposition is several times cheaper than a full one -- this runs on every
# fit, and on a model with a thousand-odd parameters the difference is seconds
# rather than milliseconds.
#' @keywords internal
.ctBackendNullMass <- function(information, rtol = 1e-12) {
  n <- nrow(information)
  values <- try(eigen(information, symmetric = TRUE,
    only.values = TRUE)$values, silent = TRUE)
  if (inherits(values, "try-error")) return(rep(NA_real_, n))
  scale <- max(values)
  if (!is.finite(scale) || scale <= 0) return(rep(NA_real_, n))
  if (!any(values <= rtol * scale)) return(rep(0, n))
  decomposition <- try(eigen(information, symmetric = TRUE), silent = TRUE)
  if (inherits(decomposition, "try-error")) return(rep(NA_real_, n))
  flat <- decomposition$values <= rtol * max(decomposition$values)
  if (!any(flat)) return(rep(0, n))
  rowSums(decomposition$vectors[, flat, drop = FALSE]^2)
}

# Population standard deviations that have collapsed to the floor of their
# transform.
#
# Reported as a finding rather than a fault. A zero population standard
# deviation says the data show no individual differences in that parameter,
# which is a result; what it also says is that the parameter sits at the edge of
# its own transform, where the curvature is zero and the reported interval is
# therefore not a confidence statement about anything.
#' @keywords internal
.ctBackendCollapsedScales <- function(fit, tolerance = 1e-6) {
  spec <- fit$model_spec
  if (is.null(spec$laplace)) return(data.frame())
  module <- .ctJuliaModule(spec$project)
  objective <- .ctJuliaObjective(fit)
  estimate <- .ctJuliaNumericVector(as.numeric(fit$estimate$raw))
  levels <- max(1L, length(spec$laplace$levels))
  rows <- list()
  for (l in seq_len(levels)) {
    covariance <- try(.ctBackendJuliaValue(module$ctsem_laplace_popcov(
      objective, estimate, as.integer(l))), silent = TRUE)
    if (inherits(covariance, "try-error")) next
    covariance <- as.matrix(covariance)
    sds <- sqrt(abs(diag(covariance)))
    name <- if (!is.null(spec$laplace$levels) &&
        !is.null(spec$laplace$levels[[l]]$name)) {
      as.character(spec$laplace$levels[[l]]$name)
    } else as.character(l)
    collapsed <- which(sds <= tolerance)
    if (length(collapsed)) {
      rows[[length(rows) + 1L]] <- data.frame(level = name,
        effect = collapsed, sd = sds[collapsed], stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# The coordinates whose reported spread the projection removed, named.
#
# Said as its own sentence because it is the one thing the surrounding warning
# does not imply: "the standard errors along those directions are arbitrary"
# prepares a reader for a number that is too big, and what they will actually
# see is a number that is too small and a z of 65. See
# `.ctBackendIntervalCheck()` for how the share is measured and why.
#' @keywords internal
.ctBackendNoWidthAdvice <- function(intervals) {
  if (is.null(intervals) || !isTRUE(intervals$nunidentified > 0L)) return("")
  named <- utils::head(intervals$unidentified, 6)
  many <- length(intervals$unidentified) > 1L
  paste0(" The reported spread for ", paste(named, collapse = ", "),
    if (length(intervals$unidentified) > 6) ", ..." else "",
    " is absent rather than small: ",
    if (many) "those coordinates lie" else "that coordinate lies",
    " in a direction with no curvature, which is left out of the inversion, ",
    "so the sd, interval and z printed for ", if (many) "them" else "it",
    " are artefacts of that projection and say nothing about the data. ",
    "summary() reports them as NA.")
}

# Say it once, at the end of a fit, in the terms a reader needs.
#' @keywords internal
.ctBackendIdentifyWarn <- function(identify, collapsed, intervals = NULL) {
  nowidth <- .ctBackendNoWidthAdvice(intervals)
  if (!is.null(identify) && identify$nweak > 0L) {
    # Which parameters are on a random-effect scale/correlation ridge and which
    # the data says nothing about, in the same words `print.ctIdentify()` uses
    # -- the advice differs between the two cases and a fit is where it is
    # most expensive to get wrong.
    partition <- .ctIdentifyPartition(identify$directions)
    warning("The data do not identify ", identify$nweak, " direction",
      if (identify$nweak > 1L) "s" else "", " of this model. The estimates ",
      "are still whatever the optimiser found, but the standard errors along ",
      "those directions are arbitrary rather than small or large, and any ",
      "interval built from them will be too. ",
      paste(.ctIdentifyAdvice(partition), collapse = " "),
      if (!length(partition$partial) && !length(partition$structural))
        paste0("Parameters involved: ",
          paste(utils::head(identify$parameters, 6), collapse = ", "),
          if (length(identify$parameters) > 6) ", ..." else "", ".") else "",
      nowidth,
      " See fit$identifiability, and ctIdentify(data, model) to check this ",
      "before spending a fit next time.", call. = FALSE)
    # Carried by the warning above rather than repeated under it: the two are
    # about one finding and a reader who has to be told twice stops reading.
    nowidth <- ""
  }
  # A null direction at 1e-12 is a flat direction at 1e-8 too, so this normally
  # travels with the warning above. It stands alone only when the two were
  # computed from different Hessians -- `ctOptimUncertainty()` re-run with a
  # different method is the case -- and then it is the more specific of the two
  # and worth saying by itself.
  if (nzchar(nowidth)) {
    warning(trimws(nowidth), " See fit$uncertainty$intervalcheck.",
      call. = FALSE)
  }
  if (!is.null(identify) && isTRUE(identify$negative > 0L)) {
    warning(identify$negative, " direction",
      if (identify$negative > 1L) "s have" else " has",
      " negative curvature at the estimate, so it is not a maximum there. ",
      "Treat the estimate as a stopping point rather than a solution.",
      call. = FALSE)
  }
  if (is.data.frame(collapsed) && nrow(collapsed)) {
    message(nrow(collapsed), " population standard deviation",
      if (nrow(collapsed) > 1L) "s were" else " was",
      " estimated at zero (level ",
      paste(unique(collapsed$level), collapse = ", "),
      "). That is a finding -- no detectable individual differences -- but the ",
      "parameter sits at the edge of its transform, where the curvature is ",
      "zero, so its reported interval is not a confidence statement. See ",
      "fit$collapsedScales.")
  }
  # Said even when no direction is flat enough to count as unidentified, which
  # is the case this exists for: a parameter can be separable in principle and
  # still have almost all of its reported width come from its entanglement with
  # the others, and then the interval moves by orders of magnitude between two
  # runs that reached the same optimum.
  if (!is.null(intervals) && isTRUE(intervals$nflagged > 0L)) {
    involved <- paste(utils::head(intervals$parameters, 6), collapse = ", ")
    if (length(intervals$parameters) > 6) involved <- paste0(involved, ", ...")
    widest <- max(intervals$table$ratio[is.finite(intervals$table$ratio)])
    warning(intervals$nflagged, " reported interval",
      if (intervals$nflagged > 1L) "s are" else " is",
      " far wider than the curvature at the estimate supports -- up to ",
      signif(widest, 3), " times the width that parameter's own curvature ",
      "gives. That width comes from the parameter not being separable from ",
      "the others rather than from the data, and it is not stable: it can ",
      "move by orders of magnitude between two fits that reach the same ",
      "optimum. Parameters involved: ", involved,
      ". See fit$uncertainty$intervalcheck.", call. = FALSE)
  }
  invisible(NULL)
}
