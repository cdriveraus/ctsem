# Score tests of relaxations of a fitted model.
#
# Every relaxation tested here is an extra coordinate whose value zero is the
# fitted model: a random effect on a parameter that does not yet vary, a time
# independent covariate's effect on a parameter, a term making the drift depend
# on the latent state. The test needs only the score in that coordinate and the
# curvature blocks that involve it, evaluated at the fit -- no refit -- and
# projects the fitted parameters out of it (Neyman's C(alpha) form), so the
# test holds at an estimate that is near, not exactly at, the maximum (priors,
# a finish that stopped on gain). The design, the literature and the
# simulations behind the choices below are in
# CT-SEM/review/SCORETESTS-misfit-design-2026-10-09.md.
#
# Two sources of the pieces:
#
#   * Per-subject: each subject's score s_i (the engine's own per-subject
#     gradient, `ctsem_subject_gradients`, which is the marginal score on every
#     route) and its Hessian H_i, by central differences of those exact
#     gradients. A covariate's effect adds to a parameter's raw value, so its
#     score is z_i s_ij exactly; a random effect's variance has score
#     (s_ij^2 + H_i,jj) / 2 exactly (Chesher, 1984). Any number of covariates
#     is then a few matrix products in R.
#   * An extended model: the fitted model with the new terms written into it at
#     zero, prepared the way the fit was (`.ctFitReplayArgs`), its gradient and
#     its Hessian columns along the new coordinates (`ctsem_hessian_columns`).
#     State dependence goes this way, because its terms enter the filter's
#     linearisation and no per-subject identity carries them.
#
# The Hessian projects and sizes; it does not set the variance. A screen meets
# models that are wrong, and when one is wrong in any way the information
# equality fails, so a variance built from the Hessian mis-scales every other
# test: an unrelated covariate's joint test rejected 27 to 44 per cent of the
# time under a missing nonlinearity. The variance is the spread of each
# subject's projected contribution (Boos, 1992), and the p-value a studentised
# bootstrap of those contributions, since a variance score is a sum of skewed
# terms; the same unrelated test then rejected 2 per cent.
#
# State dependence is written into CINT rather than into DRIFT cells. The
# engine's Jacobian is derived from DRIFT * state + CINT (`ctJacobian`), so a
# term d * eta_k * eta_m in CINT[i] is the same model as DRIFT[i,k] + d * eta_m,
# and leaves every fitted parameter where and as it was. Using the latent state
# inside the filter is what keeps the test valid: the filter evaluates the term
# at the state predicted from earlier rows, and a weight built from earlier
# rows keeps each row's score mean zero, where a smoothed or updated state does
# not (scripts/scoretests-2026-10/predictable_weights.R: a true model rejected
# at 100 per cent with those weights).

#' Score tests for relaxing a fitted model
#'
#' Tests, from one fitted model and without refitting it, whether the model
#' would fit better if it were relaxed: parameters that differ between
#' subjects (a missing random effect), parameters that depend on a time
#' independent covariate the model does not include, or dynamics that depend
#' on the latent state (a missing nonlinearity). Each test is a score
#' (Lagrange multiplier) test: the slope of the log likelihood in the
#' direction of the relaxation, at the fitted model, with the fitted
#' parameters projected out. Each also gives a one-step estimate of the size
#' of the relaxation. These estimates are for screening, so refit with the
#' relaxations that matter before interpreting their sizes.
#'
#' @section Heterogeneity: Tests letting each parameter vary over subjects,
#'   as \code{ctFit} would fit it, for each free parameter that does not vary
#'   already and is not an initial condition (\code{T0MEANS} and
#'   \code{T0VAR}). The test covers the parameter's variance over subjects
#'   together with its covariances with what already varies: each subject's
#'   initial state (the free \code{T0MEANS} parameters) and the parameters
#'   with random effects, listed in \code{with}. Even a model without random
#'   effects has initial states, so these covariances are tested there too.
#'   They matter because a varying \code{CINT} intercept that is unrelated to
#'   the initial means makes subjects drift apart from a common start, while
#'   the same intercept correlated with each subject's start is a random
#'   intercept.
#'
#'   The covariances are tested as ordinary score tests, and the variance,
#'   given them, as a one sided test (Stram & Lee, 1994), so \code{LM} is
#'   referred to an equal mixture of chi-square distributions on
#'   \code{df - 1} and \code{df} degrees of freedom. \code{z} is the one sided
#'   test of the variance alone, with its covariances held at zero, which has
#'   more power when the new effect is independent of what already varies.
#'   \code{p} is the smaller of the two tests' p-values, calibrated by
#'   resampling subjects. A negative \code{z} beside a large \code{LM} means
#'   the parameter varies only together with what already varies, as an
#'   intercept does with each subject's start. No assumption is made about how
#'   the initial states relate to the dynamics.
#'
#'   The statistic is the information matrix test, which is the score test
#'   for neglected heterogeneity (Chesher, 1984). Other misspecification can
#'   also move it. A drift that varies over subjects also moves the
#'   diffusion's test, for example, and a missing nonlinearity moves several.
#'   \code{sd} is the one-step estimate of the subjects' standard deviation on
#'   the raw scale, with the covariances estimated alongside, and
#'   \code{sd_natural} approximately the same on the parameter's own scale.
#'   They are \code{NA} where the test comes from the covariances and the
#'   variance alone is estimated as zero, as often happens for an intercept,
#'   and then only a refit can give the size. To act on a result, set
#'   \code{indvarying = TRUE} for the parameter and the \code{T0MEANS}
#'   parameters, refit, and test again. \code{\link{ctModelImprove}} does this.
#'
#'   Each \code{MANIFESTMEANS} cell the model fixes is also tested, as a
#'   random intercept on that indicator around its fixed value (rows named
#'   \code{MANIFESTMEANS[...]}). With one indicator per latent this overlaps
#'   with the covariances of the \code{CINT} rows, with more indicators it is
#'   the indicator's own intercept, and for a latent with no free intercept it
#'   is the only test of its level. To act on it, give the cell a parameter
#'   with the fixed value as its stated mean (\code{ctPar(mean=)}) and
#'   \code{indvarying = TRUE}.
#'
#'   On a reduced rank level (\code{poprank}), a parameter that varies only
#'   through the shared dimensions, without a variance of its own, is tested
#'   for variation outside those dimensions. Each such test is of the variance
#'   alone. A parameter whose test is already determined by the fitted model,
#'   as for an intercept that loads on the shared dimensions, is not tested
#'   (\code{NA}) and is left out of \code{dimension}.
#'
#'   \code{dimension} tests whether subjects vary along some combination of
#'   the tested parameters, a new dimension of variation, with its
#'   covariances with what already varies held at zero. The statistic is the
#'   largest eigenvalue of the standardised score matrix of the variances,
#'   tested against a bootstrap over subjects, and its eigenvector gives the
#'   direction, as relative loadings on the raw scale. It can detect variation
#'   spread over several parameters that no single test picks up, and is less
#'   sensitive than the best single test when one parameter carries it. It is
#'   not computed over a reduced rank level's parameters, whose blocks are
#'   tested instead.
#'
#'   On a reduced rank level, \code{blocks} tests what its blocks of shared
#'   dimensions (\code{popblocks}, or without them the whole level as one
#'   block) leave out, beside the own variances. \code{blocks$dimension}
#'   tests, within each block with more effects than dimensions, whether the
#'   block needs another shared dimension, with the block's fitted loadings
#'   projected out. The rows for single parameters, which test own variances,
#'   cannot show this. \code{blocks$covariance} tests, for each pair of
#'   blocks, the covariances between their dimensions, which the model holds
#'   at zero, with a joint two sided score test. An example is whether a
#'   response style shared by an item's thresholds covaries with the dynamic
#'   level. A pair in which one block loads on every effect of the other
#'   already has these covariances in the model, and is listed with a note but
#'   not tested. Both are calibrated by resampling subjects and adjusted
#'   together (\code{p_adj}). An effect whose population mean is stated is
#'   tested through the model with that mean free and set to its stated value,
#'   and only the fit's own parameters are projected out.
#'   \code{\link{ctModelImprove}} acts on both, adding a dimension to a block
#'   or letting blocks share dimensions.
#'
#' @section Covariates: For each covariate, tests the effects it would have on
#'   every tested parameter if it were added as a time independent predictor
#'   (on the raw scale), including effects on the population standard
#'   deviations of the parameters that vary. The \code{joint} test of all of a
#'   covariate's effects is the screen. For each parameter, \code{estimate} is
#'   the one-step effect per unit of the covariate, estimated jointly over the
#'   parameters with the other fitted parameters held at their estimates, so
#'   that a covariate that shifts one parameter is not attributed to another
#'   that is correlated with it. Holding the other parameters keeps the
#'   estimate close to what a refit would give. \code{estimate_natural} gives
#'   it approximately on the parameter's own scale, and \code{p} tests it. A
#'   covariate enters linearly (centred), as treatment indicators for a
#'   factor, or as a natural spline with \code{df} degrees of freedom, for
#'   which only \code{p} is given. Subjects with a missing value contribute at
#'   the mean. To act on a result, add the covariate to the data and to
#'   \code{TIpredNames}, and refit.
#'
#'   The model's own time independent predictors are tested as well, each on
#'   the parameters it does not yet act on, using the values the fit read.
#'   Their results have \code{joint$inmodel} set to \code{TRUE}. A predictor
#'   added with an effect on one parameter may act on another too, and that
#'   can show only once other changes are made. To act on a result, set the
#'   predictor's effect on that parameter (\code{<predictor>_effect} in
#'   \code{model$pars}) and refit. \code{\link{ctModelImprove}} does that.
#'
#'   \code{form = 'threshold'} tests for a step at an unknown cut instead. The
#'   joint test is computed at each cut from the 15th to the 85th percentile,
#'   and the largest (Andrews, 1993), which is the split a score based tree
#'   makes, is compared with the same scan with the covariate permuted over
#'   subjects (Hothorn & Lausen, 2003). This accounts for the search and needs
#'   no refitting. The effects reported are those of being above the best
#'   cut, given as \code{cut} in the joint row.
#'
#'   On a reduced rank level, a covariate's effect on the spread of a varying
#'   parameter is tested as a scaling of that parameter's loadings and own sd
#'   together, reported as \code{scale(<parameter>)}: the effect on its log
#'   standard deviation per unit of the covariate, with the correlations held
#'   fixed. With about 100 subjects these rows reject a correct model somewhat
#'   too often (6.5 to 8 per cent at .05), so read them only once the
#'   covariate's joint test has picked it out.
#'
#' @section State dependence: For each latent equation and each pair of
#'   latent states \code{k <= m}, tests a term \code{d * eta_k * eta_m} in that
#'   equation's drift, that is, a drift coefficient that varies linearly with
#'   a state (Luukkonen, Saikkonen & Terasvirta, 1988). A term
#'   \code{eta_k * eta_m} in the equation of \code{eta_i} is
#'   \code{DRIFT[i,k]} varying with \code{eta_m}, or when \code{k != m}
#'   equally \code{DRIFT[i,m]} varying with \code{eta_k}. The data cannot tell
#'   these apart, and \code{reads} gives both. \code{estimate} is the one-step
#'   change in that drift coefficient per unit of the state. It is \code{NA}
#'   when the fitted model cannot size the term, which happens when the
#'   departure is strong, and the test itself is still valid.
#'   \code{order = 3} adds the terms of degree three. The terms use the state
#'   predicted from earlier observations, as the filter does. To act on a
#'   result, write the dependence into the drift cell, for example
#'   \code{DRIFT[1,1] = 'dr11 + d * eta2'} with
#'   \code{PARS = c('dr11|-log1p_exp(-param)', 'd')}, and refit.
#'
#' @section Volatility: \code{tests = 'volatility'} tests whether the size of
#'   each latent's innovations depends on a state, by scaling
#'   \code{DIFFUSION[i, i]}, the innovation sd, by \code{exp(v * eta_m)} and
#'   testing \code{v = 0}, for every latent with innovations and every state.
#'   Variance that rises with a process's own level, as when negative affect
#'   is more variable when high, is the diagonal term. As for state
#'   dependence, the predicted state is used. To act on a result, write the
#'   scaling into the cell, for example
#'   \code{DIFFUSION[2,2] = 'sd2 * exp(v * eta2)'} with
#'   \code{PARS = c('sd2|log1p_exp(param)', 'v')}, and refit.
#'   \code{\link{ctModelImprove}} does this. These terms make the model
#'   nonlinear, so the test and a refit use the extended Kalman filter.
#'
#' @section Level: \code{tests = 'level'} tests whether a parameter of the
#'   dynamics or of the error variances depends on a subject's level, where
#'   the level is a random intercept. For each free \code{DRIFT},
#'   \code{DIFFUSION} and \code{MANIFESTVAR} parameter that does not vary over
#'   subjects, including a \code{PARS} parameter that those matrices read in
#'   an expression, and each \code{MANIFESTMEANS} or \code{CINT} intercept
#'   that does, it tests a term \code{g * intercept} added to the parameter,
#'   at \code{g = 0}, wherever the model reads the parameter. When subjects'
#'   levels are in the intercepts, the latent states are deviations from them,
#'   so the state dependence and volatility tests only see dependence on the
#'   deviation. A parameter that depends on a subject's usual level, such as
#'   persistence or innovation size that differs between people who are
#'   usually high and those usually low, otherwise shows only as
#'   heterogeneity, with no direction. These terms test that direction. The
#'   score is taken from the extended model's gradient, averaged over each
#'   subject's posterior for the intercept. \code{estimate} is the one-step
#'   change in the parameter, on its own scale, per unit of the intercept. To
#'   act on a result, let the parameter vary over subjects, so its correlation
#'   with the intercept is estimated. \code{\link{ctModelImprove}} does this.
#'   A model with no random intercept has nothing to test.
#'
#' @section Fixed cells: For each cell the model fixes at zero that could
#'   sensibly be freed (a \code{DRIFT} cross effect, a \code{DIFFUSION} or
#'   \code{MANIFESTVAR} covariance between two that have a variance, a
#'   \code{LAMBDA} cross loading, a \code{TDPREDEFFECT}), tests freeing it.
#'   This is the modification index (Sorbom, 1989), and \code{estimate} is
#'   the expected parameter change. \code{CINT} and \code{MANIFESTMEANS} are
#'   not offered, since with the other one free, either is identified only
#'   through the initial transient. Under \code{covmattransform = 'z'}, a
#'   covariance fixed at zero between variables that other covariances
#'   connect is held at zero by solving its cell. It is freed from that
#'   solved value (\code{at} in \code{cells}), and \code{estimate} is the
#'   change from it. If the solved value differs between subjects, because
#'   a covariance it depends on has random or predictor effects, the cell is
#'   not tested and is listed in \code{untested}. State dependence,
#'   volatility and fixed cells are tested on one extended model, so asking
#'   for several costs one preparation. If that model cannot be built, each
#'   family is tested separately.
#'
#' @section Stability: \code{tests = 'instability'}, run only when asked for,
#'   tests whether each parameter is stable across the data, with the data
#'   ordered by time, and by each latent's state as predicted from earlier
#'   observations. It uses each data row's score, whose partial sums behave as
#'   a Brownian bridge under the model (Nyblom, 1989). A parameter that
#'   changes over the study, or differs at high and low levels of a state,
#'   bends them. Ordered by a state, this is a test of state dependence of any
#'   form. The p-values come from a wild bootstrap over rows. It can take a
#'   long time on large data.
#'
#' @section Alike relaxations: Two relaxations whose contributions from each
#'   subject are nearly collinear (correlation of 0.9 or more) would explain
#'   the same signal, and the data cannot tell them apart (MacCallum et al.,
#'   1993). They share a letter in \code{class}, and \code{NA} means alike to
#'   none. Refit with the most plausible one, then test again.
#'
#' @section Variance and p-values: By default (\code{variance = 'robust'}) the
#'   variance of each test is the spread over subjects of their projected
#'   contributions, and its p-value comes from a bootstrap over subjects of
#'   the studentised statistic (\code{nboot} resamples, from a fixed random
#'   stream, so a call is reproducible and leaves the session's random
#'   numbers alone). This keeps each test valid when the model is wrong in
#'   other ways. \code{variance = 'model'} uses the fitted model's curvature
#'   for the variance and chi-square and normal reference distributions. It
#'   is always used with fewer than 20 subjects.
#'
#' @section Reading the results: The tests are screens. A badly misspecified
#'   model shows many significant tests at once, and the one-step sizes are
#'   local approximations, especially for drift parameters. Relax what stands
#'   out most, refit, and test again. Within each family \code{p_adj} is the
#'   Benjamini-Hochberg adjustment over the tests shown.
#'
#'   The tests use the fit's own objective, so they are for the
#'   \code{intoverpop} it was fitted with. They need an optimised fit, and
#'   refuse a sampled one. Heterogeneity and covariate tests need subjects as
#'   the independent units, so a model with a level above the subject can
#'   only be tested for state dependence. If a family cannot be tested on a
#'   model, the reason is given as a warning and in \code{info$errors}, and
#'   the other families' results stand. With \code{priors = TRUE}, the new
#'   terms are tested on the likelihood, without a prior.
#'
#'   On a Laplace fit with random dynamics beside random intercepts, the
#'   state dependence and volatility families use the one-step Laplace
#'   approximation for both the fit and its extension (see
#'   \code{optimcontrol$laplace_onestep} in \code{\link{ctFit}}). This is a
#'   fast screen whose approximation error stays in the scores, and
#'   \code{info$onestep} names the families concerned.
#'   \code{intoverpop = 'quadrature'} integrates every subject's random
#'   effects for these families instead (\code{info$quadrature}).
#'
#' @param fit An optimised fit from \code{\link{ctFit}} with
#'   \code{backend = 'julia'}.
#' @param covariates Optional data frame of subject level covariates to test,
#'   with a column named as the fit's subject id (or \code{id}) and one column
#'   per covariate. Long data with repeated rows per subject can be used if
#'   each covariate is constant within subject. Without an id column, give one
#'   row per subject in the fit's subject order.
#' @param tests Which families to run: any of \code{'heterogeneity'},
#'   \code{'covariates'}, \code{'statedependence'}, \code{'volatility'},
#'   \code{'fixed'}, \code{'level'} and \code{'instability'}. By default all
#'   but \code{'instability'}. Covariates are tested only when
#'   \code{covariates} is given or \code{'covariates'} is asked for (then the
#'   model's own predictors alone), and levels only where the model has a
#'   random intercept.
#' @param form How a covariate enters: \code{'auto'} (indicators for a factor,
#'   character or logical covariate, linear otherwise), \code{'linear'},
#'   \code{'spline'}, \code{'factor'} or \code{'threshold'} (a step at an
#'   unknown cut). One value, or one per covariate.
#' @param df Degrees of freedom of the spline basis.
#' @param parameters Optional parameter names to restrict the heterogeneity and
#'   covariate tests to.
#' @param order Highest degree of the state dependence terms. 2 tests each
#'   drift coefficient varying linearly with each state. 3 adds variation with
#'   products and squares of states, which is needed to see a symmetric
#'   dependence, such as regulation that strengthens at both extremes.
#' @param variance \code{'robust'} (the default) or \code{'model'}; see the
#'   section on variance.
#' @param nboot Bootstrap resamples for the robust p-values.
#' @param exact For \code{statedependence}, \code{volatility} and
#'   \code{fixed} with the robust variance: how many of each family's terms
#'   are corrected exactly for the estimation of the fitted parameters, which
#'   takes two evaluations of the extended model per term. Every term is first
#'   tested without that correction, which is conservative, and the strongest
#'   \code{exact} terms are then corrected (\code{projected} in the tables).
#'   \code{Inf} corrects every term.
#' @param intoverpop How the extended families (\code{statedependence},
#'   \code{volatility}, \code{fixed}, \code{level}) integrate the random
#'   effects on a Laplace fit. \code{'laplace'} (the default) uses the
#'   Laplace approximation, a fast screen. \code{'quadrature'} integrates
#'   every subject's random effects by quadrature, as
#'   \code{ctFit(intoverpop = 'quadrature')} does, which removes the Laplace
#'   error from these families and takes longer. The other families use the
#'   fit's own objective either way.
#' @param cores Engine cores for the evaluations.
#' @param verbose Report progress.
#' @return An object of class \code{ctScoreTest}, a list with
#'   \code{heterogeneity} (a data frame, one row per parameter or fixed
#'   \code{MANIFESTMEANS} cell), \code{dimension} (the largest eigenvalue, its
#'   p-value and the direction, when two or more parameters are tested),
#'   \code{blocks} (on a reduced rank level, \code{dimension}, a row per block
#'   with the one-step \code{sd} along its direction, \code{covariance}, a row
#'   per pair of blocks, and \code{directions}, the direction found for each
#'   tested block),
#'   \code{covariates} (one element per covariate, each a list of its
#'   \code{joint} test and a data frame \code{parameters}, including the
#'   model's own predictors, marked by \code{joint$inmodel}),
#'   \code{statedependence} (a list of \code{terms}, a data frame, and
#'   \code{joint}, the tests per equation and overall), \code{volatility}
#'   (the same, for the innovation sds), \code{fixed} (a list of \code{cells},
#'   their \code{joint} test, and \code{untested}, cells that could not be
#'   tested and why), \code{instability} (a list of
#'   \code{parameters}, one row per ordering and parameter, and \code{joint},
#'   per ordering), and \code{info} (how the fit integrated its random
#'   effects, the number of subjects and parameters, the variance used, and
#'   \code{errors}, why any family was not tested). Tables have a
#'   \code{class} column where relaxations are alike, and the extended
#'   families a \code{projected} column, whether a term was corrected exactly
#'   (see \code{exact}).
#' @references Boos, D. D. (1992). On generalized score tests. The American
#'   Statistician, 46, 327-333.
#'
#'   Chesher, A. (1984). Testing for neglected heterogeneity. Econometrica, 52,
#'   865-872.
#'
#'   Stram, D. O., & Lee, J. W. (1994). Variance components testing in the
#'   longitudinal mixed effects model. Biometrics, 50, 1171-1177.
#'
#'   Self, S. G., & Liang, K.-Y. (1987). Asymptotic properties of maximum
#'   likelihood estimators and likelihood ratio tests under nonstandard
#'   conditions. Journal of the American Statistical Association, 82, 605-610.
#'
#'   Hothorn, T., & Lausen, B. (2003). On the exact distribution of maximally
#'   selected rank statistics. Computational Statistics & Data Analysis, 43,
#'   121-137.
#'
#'   Andrews, D. W. K. (1993). Tests for parameter instability and structural
#'   change with unknown change point. Econometrica, 61, 821-856.
#'
#'   Nyblom, J. (1989). Testing for the constancy of parameters over time.
#'   Journal of the American Statistical Association, 84, 223-230.
#'
#'   Hall, D. B., & Praestgaard, J. T. (2001). Order-restricted score tests for
#'   homogeneity in generalised linear and nonlinear mixed models.
#'   Biometrika, 88, 739-751.
#'
#'   MacCallum, R. C., Wegener, D. T., Uchino, B. N., & Fabrigar, L. R. (1993).
#'   The problem of equivalent models in applications of covariance structure
#'   analysis. Psychological Bulletin, 114, 185-199.
#'
#'   Luukkonen, R., Saikkonen, P., & Terasvirta, T. (1988). Testing linearity
#'   against smooth transition autoregressive models. Biometrika, 75, 491-499.
#'
#'   Sorbom, D. (1989). Model modification. Psychometrika, 54, 371-384.
#' @seealso \code{vignette('score-tests')} for three worked examples.
#' @examples
#' \donttest{
#' if (isTRUE(ctJuliaStatus()$available)) {
#' # A model in which the innovations of the two processes are uncorrelated,
#' # with the DIFFUSION cell between them fixed at zero.
#' model <- ctModel(LAMBDA = diag(2),
#'   manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
#'   TDpredNames = 'TD1',
#'   DIFFUSION = matrix(c('diff_eta1', 0, 0, 'diff_eta2'), 2, 2))
#' fit <- ctFit(ctsemtestdat, model, priors = TRUE)
#'
#' # Every family of relaxations, tested from this one fit. Under "Fixed
#' # cells, freed", freeing DIFFUSION[eta2,eta1] has a small p-value, so the
#' # data are not consistent with a zero covariance.
#' st <- ctScoreTest(fit)
#' print(st)
#'
#' # One family's table: each fixed cell with its test, and an estimate of
#' # the parameter if it were freed.
#' print(st$fixed$cells)
#'
#' # Whether subject level covariates, here the time independent predictors
#' # in the data, predict any parameter. The joint test is over all
#' # parameters.
#' stc <- ctScoreTest(fit,
#'   covariates = ctsemtestdat[, c('id', 'TI1', 'TI2', 'TI3')])
#' print(stc$covariates$TI1$joint)
#' }
#' }
#' @export
ctScoreTest <- function(fit, covariates = NULL,
  tests = c('heterogeneity', 'covariates', 'statedependence', 'volatility', 'fixed',
    'level', 'instability'),
  form = c('auto', 'linear', 'spline', 'factor', 'threshold'), df = 3, parameters = NULL,
  order = 2L, variance = c('robust', 'model'), nboot = 1999L, exact = 5L,
  intoverpop = c('laplace', 'quadrature'), cores = getOption('mc.cores', 2L),
  verbose = TRUE) {
  .ctWeightsRefuse(fit, 'ctScoreTest()')
  intoverpop <- match.arg(intoverpop)
  if (!.ctFitIsJulia(fit)) stop("ctScoreTest() needs a fit from ",
    "ctFit(backend = 'julia'). Refit on that backend to use it.", call. = FALSE)
  if (!is.null(fit$sample)) stop("Score tests are made at an optimum, and this ",
    "fit was sampled. Fit with optimize = TRUE to test it.", call. = FALSE)
  if (isFALSE(fit$args$resolved$intoverstates)) stop("intoverstates = FALSE ",
    "targets the joint density of parameters and states, whose maximum is not ",
    "an estimate of anything (see ?ctFit); there is no fitted model to test.",
    call. = FALSE)
  # Before anything reaches Julia: a session this call starts is as wide as
  # `cores`, as ctFit's is. Without it a saved fit tested in a fresh R session
  # (a cluster job) ran every pass at getOption('mc.cores', 2) threads.
  .ctBackendResolveThreads(cores)
  testsgiven <- !missing(tests)
  tests <- match.arg(tests, several.ok = TRUE)
  # Stability needs a dual pass of the filter per eight parameters (on the
  # Laplace route one per truncated record), so it is run when asked for.
  if (!testsgiven) tests <- setdiff(tests, 'instability')
  # Asked for without covariates, the family tests the model's own predictors
  # on the parameters they do not act on (`.ctScoreOwnPredictors()`).
  if ('covariates' %in% tests && is.null(covariates)) {
    if (testsgiven && !length(.ctScoreOwnNames(fit))) stop("tests includes ",
      "'covariates' but no covariates were given, and the model has no time ",
      "independent predictors of its own.", call. = FALSE)
    if (!testsgiven) tests <- setdiff(tests, 'covariates')
  }
  if (!is.null(covariates) && !'covariates' %in% tests && testsgiven)
    message("covariates given but not tested: tests does not include 'covariates'.")
  # One form for every covariate unless several are given; match.arg() with
  # several.ok would return the whole default vector for a missing argument.
  form <- if (missing(form)) 'auto' else match.arg(form, several.ok = TRUE)
  say <- function(...) if (isTRUE(verbose)) message(...)
  variance <- match.arg(variance)
  if (!is.numeric(exact) || length(exact) != 1L || is.na(exact) || exact < 0)
    stop("exact is a number of terms, 0 or more (Inf for every term).", call. = FALSE)
  if (variance == 'robust' && .ctFitNsubjects(fit) < 20L) {
    say("Fewer than 20 subjects: using the model-based variance.")
    variance <- 'model'
  }
  # The bootstrap draws from its own stream, so a call gives the same answer
  # every time and leaves the caller's random numbers where they were.
  restore <- .ctStreamKeep(20261009L)
  on.exit(restore(), add = TRUE)

  spec <- fit$model_spec
  est <- .ctFitRawEstimate(fit)
  # The fit's objective is the one its run ended on, its modes held
  # (R/ctJuliaFitModes.R): a score differenced either side of the estimate
  # must not see a unit land on another maximum.
  .ctJuliaHoldModes(spec, .ctJuliaFitModes(fit))
  ptab <- .ctScoreParameterTable(fit)
  out <- list(heterogeneity = NULL, covariates = NULL, statedependence = NULL,
    volatility = NULL, fixed = NULL, level = NULL, info = list(route = .ctBackendIntOverPop(spec),
      nsubjects = .ctFitNsubjects(fit), npar = length(est)))
  out$info$variance <- variance
  out$info$nboot <- nboot
  # A family that cannot be tested on this model says why, in the result and
  # as a warning, and the others stand: the pilot's state-dependence refusal
  # once discarded a heterogeneity family that had taken twenty minutes.
  attempt <- function(label, expr) {
    ok <- tryCatch({ force(expr); TRUE }, error = function(e) {
      out$info$errors[[label]] <<- conditionMessage(e)
      warning(label, ": ", conditionMessage(e), call. = FALSE)
      FALSE
    })
    invisible(ok)
  }

  persubject <- intersect(tests, c('heterogeneity', 'covariates'))
  pieces <- NULL
  # A reduced-rank level (poprank) carries its variances as loadings and own
  # sds, so a parameter already varying can still lack a dimension of
  # variation, and a covariate's effect on its spread has no popsd coordinate
  # to act on: it acts on the scale of the parameter's row of the loadings.
  reduced <- any(grepl('^poploading_', ptab$name))
  if (length(persubject)) {
    het <- .ctScoreSelect(ptab, parameters, 'heterogeneity', reduced)
    cov <- .ctScoreSelect(ptab, parameters, 'covariates', reduced)
    if (!'heterogeneity' %in% tests) het <- integer()
    if (!'covariates' %in% tests) cov <- integer()
    # The covariates and their bases are checked before any pass is paid for:
    # a wrong id column or a non-numeric threshold is the caller's to fix, and
    # stops the call, where a family the model cannot support does not.
    bases <- NULL
    if (length(cov) && !is.null(covariates)) {
      values <- .ctScoreCovariateValues(fit, covariates)
      forms <- rep_len(form, length(values))
      bases <- Map(function(z, f, name) {
        if (identical(f, 'threshold')) {
          if (is.factor(z) || is.character(z) || is.logical(z)) stop("Covariate '", name,
            "' is not numeric; a threshold needs an ordered covariate.", call. = FALSE)
          list(threshold = z)
        } else .ctScoreBasis(z, f, df, name)
      }, values, forms, names(values))
    }
    # The model's own predictors, each on the model parameters among these it
    # does not act on; without covariates given, those parameters are all the
    # family tests.
    own <- if (length(cov)) .ctScoreOwnPredictors(fit, ptab, cov) else NULL
    if (is.null(covariates)) cov <- sort(unique(unlist(lapply(own, `[[`, 'targets'))))
    dirs <- if (reduced && length(cov) && !is.null(covariates))
      .ctScoreScaleDirections(est, ptab, parameters) else NULL
    # Each heterogeneity test takes in its covariances with the initial means
    # and the effects already varying (.ctScoreExisting): their rows, and a
    # pass along each tested parameter and effect together for the mixed
    # third derivatives the covariances' projections need.
    R <- if (length(het)) .ctScoreExisting(ptab, reduced) else integer()
    pairs <- if (length(R)) {
      g <- as.matrix(expand.grid(het, R))
      unname(g[g[, 1] != g[, 2], , drop = FALSE])
    } else NULL
    # On a reduced-rank level, the blocks of its shared dimensions
    # (`.ctScoreBlockPlan()`): each block's effects, for the test of a further
    # dimension within it, and its loadings as directions, for the test of
    # covariances between blocks, which needs each direction's pass and one
    # along every pair of directions from two blocks.
    plan <- if (reduced && 'heterogeneity' %in% tests && is.null(parameters))
      .ctScoreBlockPlan(fit, est) else NULL
    # An effect whose population mean is stated has no coordinate in the fit,
    # so where a block holds one its tests are taken on the model with those
    # means estimated, at their stated values (`.ctScoreStatedSpec()`).
    stated <- unique(unlist(lapply(plan, `[[`, 'stated')))
    bext <- if (length(stated)) tryCatch(.ctScoreStatedSpec(fit, est, stated, cores),
      error = function(e) {
        out$info$errors[['heterogeneity, blocks: stated means']] <<- conditionMessage(e)
        warning("heterogeneity, blocks: ", conditionMessage(e), "; effects with a ",
          "stated mean are left out of the block tests.", call. = FALSE)
        NULL
      })
    onext <- !is.null(bext)
    if (onext) plan <- .ctScoreBlockPlan(fit, est, bext)
    within <- unique(unlist(lapply(plan, function(b)
      if (length(b$index) > b$rank) b$index)))
    vdirs <- do.call(cbind, lapply(plan, `[[`, 'dirs'))
    vblock <- unlist(lapply(seq_along(plan), function(i) rep(i, NCOL(plan[[i]]$dirs))))
    vpairs <- .ctScoreBlockPairs(plan, vdirs, vblock)
    nested <- attr(vpairs, 'nested')
    if (!NROW(vpairs)) vpairs <- NULL
    rows <- union(union(union(het, cov), R), if (!onext) within)
    say("Score contributions: ", 2L * (length(rows) + NCOL(dirs) + NROW(pairs) +
      NCOL(vdirs) + NROW(vpairs)) + 1L, " per-subject gradient passes")
    pieces <- if (length(rows) || !onext) .ctScorePieces(fit, est, rows = rows,
      second = union(union(het, R), if (!onext) within), dirs = dirs, pairs = pairs,
      vdirs = if (!onext) vdirs, vpairs = if (!onext) vpairs, cores = cores) else
      .ctBackendWithMaxChunks(cores, list(S = .ctBackendScoreMatrix(fit, est),
        Hsum = local({ H <- .ctBackendHessian(fit, est); (H + t(H)) / 2 })))
    if (nrow(pieces$S) != .ctFitNsubjects(fit)) stop("This model has a level ",
      "above the subject, so its independent units are not subjects. ",
      "Heterogeneity and covariate tests need subjects; use ",
      "tests = 'statedependence'.", call. = FALSE)
    if (length(het)) attempt('heterogeneity', {
      items <- lapply(het, function(j) .ctScoreHetItem(pieces, ptab, j, R))
      # Random intercepts on the MANIFESTMEANS the model fixes, through the
      # model with those cells free; if that cannot be prepared, the rest
      # stand.
      if (is.null(parameters)) items <- c(items, tryCatch(
        .ctScoreFixedMeans(fit, est, pieces, R, cores = cores), error = function(e) {
          out$info$errors[['heterogeneity, fixed MANIFESTMEANS']] <<- conditionMessage(e)
          warning("heterogeneity, fixed MANIFESTMEANS: ", conditionMessage(e), call. = FALSE)
          NULL
        }))
      out$heterogeneity <- .ctScoreHeterogeneity(pieces, ptab, items, variance, nboot)
      # Not over a reduced-rank level's parameters: there it found a missing
      # dimension 10 per cent of the time where the per-parameter test of
      # variation outside the shared dimensions found it 58 (simulation, rank
      # 1 fitted to two independent effects). Its blocks are tested below.
      if (length(het) > 1L && !reduced) out$dimension <- .ctScoreDimension(pieces, ptab,
        het, nboot)
    })
    if (length(plan)) attempt('heterogeneity, blocks', {
      out$blocks <- if (onext) .ctScoreBlocks(.ctScoreStatedPieces(bext, pieces, within,
        vdirs, vpairs, cores), list(name = names(bext$x2)), plan, vblock, variance, nboot,
        nested = nested) else .ctScoreBlocks(pieces, ptab, plan, vblock, variance, nboot,
        nested = nested)
    })
    if (length(cov)) attempt('covariates', {
      items <- .ctScoreItems(pieces, ptab, cov)
      if (length(bases)) out$covariates <- Map(function(G, name)
        if (is.list(G) && !is.null(G$threshold))
        .ctScoreThreshold(pieces, items, G$threshold, name, variance, nboot) else
        .ctScoreCovariate(pieces, items, G, variance, nboot), bases, names(bases))
      itemnames <- vapply(items, `[[`, '', 'name')
      for (name in names(own)) {
        G <- tryCatch(.ctScoreBasis(own[[name]]$values, 'linear', df, name),
          error = function(e) NULL)
        if (is.null(G)) next
        res <- .ctScoreCovariate(pieces, items[match(ptab$name[own[[name]]$targets],
          itemnames)], G, variance, nboot)
        res$joint$inmodel <- TRUE
        out$covariates[[name]] <- res
      }
    })
  }
  extended <- intersect(tests, c('statedependence', 'volatility', 'fixed'))
  information <- if (!is.null(pieces)) -pieces$Hsum else NULL
  # One extended model for these families, one compile of its shape. Where
  # that model is refused, each family is tested on a model of its own, so
  # that one family's refusal leaves the others tested.
  extend <- function(families) {
    say("Preparing the model extended at zero (", paste(families, collapse = ', '), ")")
    ext <- .ctScoreExtended(fit, est, families, information, order = order,
      variance = variance, nboot = nboot, cores = cores, exact = exact,
      quadrature = intoverpop == 'quadrature')
    for (f in families) out[[f]] <<- ext[[f]]
    if (isTRUE(ext$onestep)) out$info$onestep <<- c(out$info$onestep, families)
    if (isTRUE(ext$quadrature)) out$info$quadrature <<- c(out$info$quadrature, families)
  }
  if (length(extended)) attempt(paste(extended, collapse = ' and '),
    tryCatch(extend(extended), ctScoreRefused = function(e) {
      if (length(extended) == 1L) stop(e)
      for (f in extended) attempt(f, extend(f))
    }))
  # A model of its own: its terms rewrite the cells the other families'
  # terms rewrite too.
  if ('level' %in% tests) attempt('level', {
    say("Preparing the model extended at zero (level)")
    lev <- .ctScoreExtended(fit, est, 'level', information, order = order,
      variance = variance, nboot = nboot, cores = cores, exact = exact,
      quadrature = intoverpop == 'quadrature')
    out$level <- lev$level
    if (isTRUE(lev$onestep)) out$info$onestep <- c(out$info$onestep, 'level')
    if (isTRUE(lev$quadrature)) out$info$quadrature <- c(out$info$quadrature, 'level')
  })
  if ('instability' %in% tests) attempt('instability', {
    say("Row scores for the stability tests")
    out$instability <- .ctScoreInstability(fit, est, ptab, parameters, nboot, cores)
  })
  structure(out, class = 'ctScoreTest')
}

# One row per raw parameter: its name, what kind it is, the cell it lives in,
# whether it already varies over subjects, and the slope of its transform at
# the estimate (to put one-step sizes on the parameter's own scale).
.ctScoreParameterTable <- function(fit) {
  est <- .ctFitRawEstimate(fit)
  name <- as.character(.ctFitRawParNames(fit))
  kind <- ifelse(grepl('^popsd_', name), 'popsd',
    ifelse(grepl('^(rawcor_|poploading_|popunique_|rawtipredeffect_|tipredvalue_|raw\\[)', name),
      'other', 'model'))
  base <- .ctFitBaseModel(fit)
  if (is.null(base)) base <- .ctFitModelObject(fit)
  pars <- base$pars
  at <- .ctParamRows(base, name)
  ns <- asNamespace('ctsem')
  slope <- vapply(seq_along(name), function(k) {
    if (kind[k] != 'model' || is.na(at[k])) return(NA_real_)
    tf <- pars$transform[at[k]]
    if (is.null(tf) || is.na(tf) || !nzchar(tf)) return(1)
    f <- function(x) eval(parse(text = tf), list(param = x), ns)
    h <- 1e-6 * max(1, abs(est[k]))
    out <- try((f(est[k] + h) - f(est[k] - h)) / (2 * h), silent = TRUE)
    if (inherits(out, 'try-error') || !is.finite(out)) NA_real_ else out
  }, numeric(1))
  data.frame(index = seq_along(name), name = name, kind = kind,
    matrix = ifelse(is.na(at), NA_character_, as.character(pars$matrix[at])),
    row = ifelse(is.na(at), NA_integer_, as.integer(pars$row[at])),
    col = ifelse(is.na(at), NA_integer_, as.integer(pars$col[at])),
    indvarying = !is.na(at) & as.logical(pars$indvarying[at]) %in% TRUE,
    slope = slope, stringsAsFactors = FALSE)
}

# The raw positions a family tests: model parameters for both, minus the
# initial conditions and those already varying for heterogeneity (on a
# reduced-rank level, those with a variance of their own); plus the population
# scales for covariates. `parameters` narrows either.
.ctScoreSelect <- function(ptab, parameters, family, reduced = FALSE) {
  # On a reduced-rank level a varying parameter with no variance of its own
  # (no popunique_) has none outside the shared dimensions: one of its own is
  # a boundary, testable as for a parameter that does not vary. A parameter
  # with its own variance already has that direction free, so a test there
  # asks nothing the fit did not estimate.
  own <- paste0('popunique_', ptab$name) %in% ptab$name
  keep <- if (family == 'heterogeneity') {
    ptab$kind == 'model' & (!ptab$indvarying | (reduced & !own)) &
      !ptab$matrix %in% c('T0MEANS', 'T0VAR')
  } else ptab$kind %in% c('model', 'popsd')
  if (!is.null(parameters)) {
    unknown <- setdiff(parameters, ptab$name)
    if (length(unknown)) stop("Unknown parameters: ", paste(unknown, collapse = ', '),
      ". The fit's are: ", paste(ptab$name, collapse = ', '), call. = FALSE)
    keep <- keep & ptab$name %in% parameters
  }
  ptab$index[keep]
}

# Each subject's score at the estimate, and its Hessian rows for the tested
# parameters, by central differences of the engine's exact per-subject
# gradients along the tested parameters only (2 * length(rows) + 1 passes): a
# subject's Hessian is symmetric, so its column along a tested parameter is the
# row the tests read. The same passes give, for the parameters in `second`,
# the second derivative of every score component along that parameter, which
# is the third derivative the heterogeneity test's projection needs:
# d H_i[j,j] / d theta_k = d^2 s_ik / d theta_j^2. A step of 1e-3 on the raw
# scale balances truncation (h^2) against the rounding of the second
# difference (eps / h^2). The full information is the fit's own Hessian, which
# the optimiser's finish formed already; differencing along every parameter
# for it cost 2p + 1 passes, 705 on a 352-parameter fit.
#
# `dirs` are directions a covariate acts along (H_i w only). `pairs` (j, r),
# both in `second`, are what a covariance's projection needs: along e_j + e_r
# the second difference is d2_jj + 2 d2_jr + d2_rr, so with j's and r's own
# it gives d2 s_ik / d theta_j d theta_r, summed over subjects (`Tmix`, a
# column a pair).
#
# `vdirs` are directions whose second derivatives a covariance between them
# needs (a reduced-rank level's loadings, `.ctScoreBlocks()`): along each, H_i v
# (`Hv`) and the summed second difference of the scores (`Dv`); along each pair
# of them in `vpairs`, the summed mixed second derivative v_a' d2 s v_b
# (`Vmix`), from the difference along their sum.
.ctScorePieces <- function(fit, est, rows, second, dirs = NULL, pairs = NULL,
  vdirs = NULL, vpairs = NULL, step = 1e-3, cores = 2L, score = NULL, hessian = NULL) {
  p <- length(est)
  if (!length(rows) && !NCOL(vdirs)) stop("No parameters to test.", call. = FALSE)
  Hdir <- Tmix <- Hv <- Dv <- Vmix <- NULL
  along <- function(w) step * max(1, max(abs(est[w != 0]))) / max(abs(w))
  # Each subject's scores, and the summed Hessian: the fit's, or another
  # objective's (`.ctScoreStatedPieces()`).
  scores <- if (is.null(score)) function(x) .ctBackendScoreMatrix(fit, x) else score
  .ctBackendWithMaxChunks(cores, {
    S <- scores(est)
    n <- nrow(S)
    H <- array(0, c(n, length(rows), p))
    D2 <- array(0, c(n, p, length(second)))
    for (a in seq_along(rows)) {
      k <- rows[a]
      h <- step * max(1, abs(est[k]))
      e <- est; e[k] <- est[k] + h
      plus <- scores(e)
      e[k] <- est[k] - h
      minus <- scores(e)
      H[, a, ] <- (plus - minus) / (2 * h)
      j <- match(k, second)
      if (!is.na(j)) D2[, , j] <- (plus - 2 * S + minus) / h^2
    }
    # Along each direction, the same difference gives H_i w.
    if (!is.null(dirs)) {
      Hdir <- array(0, c(n, ncol(dirs), p))
      for (d in seq_len(ncol(dirs))) {
        w <- dirs[, d]
        h <- step * max(1, max(abs(est[w != 0]))) / max(abs(w))
        plus <- scores(est + h * w)
        minus <- scores(est - h * w)
        Hdir[, d, ] <- (plus - minus) / (2 * h)
      }
    }
    if (NROW(pairs)) {
      own <- function(j) colSums(matrix(D2[, , match(j, second)], n))
      Tmix <- matrix(0, p, nrow(pairs))
      for (k in seq_len(nrow(pairs))) {
        jr <- pairs[k, ]
        w <- replace(numeric(p), jr, 1)
        h <- step * max(1, abs(est[jr]))
        d2 <- colSums(scores(est + h * w) - 2 * S +
          scores(est - h * w)) / h^2
        Tmix[, k] <- (d2 - own(jr[1]) - own(jr[2])) / 2
      }
    }
    if (NCOL(vdirs)) {
      Hv <- array(0, c(n, ncol(vdirs), p))
      Dv <- matrix(0, p, ncol(vdirs))
      for (d in seq_len(ncol(vdirs))) {
        h <- along(vdirs[, d])
        plus <- scores(est + h * vdirs[, d])
        minus <- scores(est - h * vdirs[, d])
        Hv[, d, ] <- (plus - minus) / (2 * h)
        Dv[, d] <- colSums(plus - 2 * S + minus) / h^2
      }
      if (NROW(vpairs)) {
        Vmix <- matrix(0, p, nrow(vpairs))
        for (k in seq_len(nrow(vpairs))) {
          w <- vdirs[, vpairs[k, 1]] + vdirs[, vpairs[k, 2]]
          h <- along(w)
          d2 <- colSums(scores(est + h * w) - 2 * S +
            scores(est - h * w)) / h^2
          Vmix[, k] <- (d2 - Dv[, vpairs[k, 1]] - Dv[, vpairs[k, 2]]) / 2
        }
      }
    }
    Hsum <- if (is.null(hessian)) .ctBackendHessian(fit, est) else hessian
  })
  Hsum <- (Hsum + t(Hsum)) / 2
  list(S = S, H = H, rows = rows, D2 = D2, second = second, Hsum = Hsum,
    dirs = dirs, Hdir = Hdir, pairs = pairs, Tmix = Tmix, vdirs = vdirs, Hv = Hv,
    Dv = Dv, vpairs = vpairs, Vmix = Vmix)
}

# What a new random effect's covariances are tested with: each initial mean (a
# free T0MEANS parameter, whether or not it varies yet) and the population
# mean of every parameter that varies already. Every subject has an initial
# state, so there is always something for an effect to covary with; T0VAR, or
# random T0MEANS, give it its spread. Not on a reduced-rank level, where the
# shared dimensions carry the covariances and what is tested is variation
# outside them.
.ctScoreExisting <- function(ptab, reduced = FALSE) {
  if (reduced) return(integer())
  ptab$index[ptab$kind == 'model' & (ptab$matrix %in% 'T0MEANS' | ptab$indvarying)]
}

# A tested coordinate's pieces, from the fit's own passes: its score, its
# Hessian row, its own second derivative, the summed third derivative its
# variance's projection needs, and for each effect in R the mixed one its
# covariance's projection needs.
.ctScoreHetItem <- function(pieces, ptab, j, R) {
  n <- nrow(pieces$S)
  R <- setdiff(R, j)
  H <- .ctScoreHrow(pieces, j)
  at <- match(paste(j, R), paste(pieces$pairs[, 1], pieces$pairs[, 2]))
  list(name = ptab$name[j], matrix = ptab$matrix[j], row = ptab$row[j],
    col = ptab$col[j], slope = ptab$slope[j], s = pieces$S[, j], H = H, hjj = H[, j],
    D2 = colSums(matrix(pieces$D2[, , match(j, pieces$second)], n)), R = R,
    Tmix = if (length(R)) pieces$Tmix[, at, drop = FALSE] else NULL)
}

# The statistic for a variance tested with its covariances (Stram & Lee,
# 1994; Self & Liang, 1987): U and V over (variance, covariances). The
# covariances are an ordinary score test (W); the variance, given them, is
# one sided (z). LM = W + max(0, z)^2, against an equal mixture of chi-square
# on q and q + 1, q the covariances; with none it is the one-sided z alone.
.ctScoreMixture <- function(U, V) {
  q <- length(U) - 1L
  if (q == 0L) {
    z <- if (V[1, 1] > 0) U[1] / sqrt(V[1, 1]) else NA_real_
    return(list(LM = max(0, z)^2, z = z, q = 0L))
  }
  pk <- .ctScorePinv(V[-1, -1, drop = FALSE])
  b <- as.numeric(pk$inv %*% V[-1, 1])
  W <- as.numeric(t(U[-1]) %*% pk$inv %*% U[-1])
  v <- V[1, 1] - sum(V[1, -1] * b)
  z <- if (v > 0) (U[1] - sum(b * U[-1])) / sqrt(v) else NA_real_
  list(LM = W + max(0, z, na.rm = TRUE)^2, z = z, q = pk$rank)
}

.ctScoreMixtureP <- function(LM, q) {
  # chi-square on 0 is a point mass at zero
  p0 <- ifelse(q > 0, stats::pchisq(LM, pmax(q, 1), lower.tail = FALSE), as.numeric(LM <= 0))
  0.5 * p0 + 0.5 * stats::pchisq(LM, q + 1, lower.tail = FALSE)
}

# The raw-coordinate direction that scales one varying parameter's spread on a
# reduced-rank level: its row of the loadings (`poploading_<par>_dim<q>`, raw
# times the parameter's sd scale) and its own sd if it has one
# (`popunique_<par>`, log1p_exp(2 raw - 1) times the scale; see
# `_laplace_poploading`). Multiplying that row by exp(g) moves each loading's
# raw value by g times itself and the own sd's by g f / f'. A covariate's
# effect on g is its effect on the parameter's log sd, the correlations held.
.ctScoreScaleDirections <- function(est, ptab, parameters = NULL) {
  nm <- ptab$name
  load <- grepl('^poploading_.*_dim[0-9]+$', nm)
  own <- grepl('^popunique_', nm) & !grepl('\\.', nm)
  pars <- unique(c(sub('_dim[0-9]+$', '', sub('^poploading_', '', nm[load])),
    sub('^popunique_', '', nm[own])))
  if (!is.null(parameters)) pars <- intersect(pars, parameters)
  W <- matrix(0, length(est), 0)
  for (par in pars) {
    w <- numeric(length(est))
    l <- which(load & startsWith(nm, paste0('poploading_', par, '_dim')) &
      grepl('^[0-9]+$', substring(nm, nchar(paste0('poploading_', par, '_dim')) + 1L)))
    w[l] <- est[l]
    u <- which(own & nm == paste0('popunique_', par))
    if (length(u)) {
      x <- 2 * est[u] - 1
      w[u] <- log1p(exp(x)) / (2 * stats::plogis(x))
    }
    if (any(w != 0)) W <- cbind(W, matrix(w, ncol = 1, dimnames = list(NULL, paste0('scale(', par, ')'))))
  }
  if (!ncol(W)) NULL else W
}

# What a covariate's effects act on, each with its score, its Hessian row and
# its direction in raw coordinates: the tested parameters, and on a
# reduced-rank level the scale of each varying parameter's spread.
.ctScoreItems <- function(pieces, ptab, cov) {
  p <- ncol(pieces$S)
  items <- lapply(cov, function(j) list(name = ptab$name[j], matrix = ptab$matrix[j],
    slope = ptab$slope[j], s = pieces$S[, j], H = .ctScoreHrow(pieces, j),
    dir = replace(numeric(p), j, 1)))
  if (!is.null(pieces$dirs)) for (d in seq_len(ncol(pieces$dirs))) {
    H <- pieces$Hdir[, d, , drop = TRUE]
    if (is.null(dim(H))) H <- matrix(H, nrow = nrow(pieces$S))
    items[[length(items) + 1L]] <- list(name = colnames(pieces$dirs)[d],
      matrix = 'POPCOV', slope = 1, s = as.numeric(pieces$S %*% pieces$dirs[, d]),
      H = H, dir = pieces$dirs[, d])
  }
  items
}

# A symmetric pseudo-inverse and its rank: a direction can carry no
# information here (two parameters whose scores coincide) without that being
# an error.
.ctScorePinv <- function(M, rtol = 1e-10) {
  M <- (M + t(M)) / 2
  e <- eigen(M, symmetric = TRUE)
  keep <- e$values > rtol * max(abs(e$values), .Machine$double.eps)
  inv <- e$vectors[, keep, drop = FALSE] %*%
    (t(e$vectors[, keep, drop = FALSE]) / e$values[keep])
  list(inv = inv, rank = sum(keep))
}

# A joint test of the score `U` with covariance `V`: its statistic, the rank
# it is on, and its p-value, from a studentised bootstrap of the units'
# contributions `C` when the variance is theirs, else from the chi-square.
# A robust joint test needs fewer terms than the units' contributions span:
# with n units their spread has rank n - 1 at most, and at that rank the
# statistic carries no information. On the Fisher search's final model, 42
# state dependence terms over 40 people gave LM 902 on 39 df and a bootstrap
# p of .0035, where no term came near and neither equation's joint test did
# (p .98 and .84). NA there.
.ctScoreJoint <- function(U, V, C = NULL, nboot = 1999L) {
  j <- .ctScorePinv(V)
  if (!j$rank) return(list(LM = NA_real_, df = 0L, p = NA_real_))
  if (!is.null(C) && j$rank >= nrow(as.matrix(C)) - 1L)
    return(list(LM = NA_real_, df = j$rank, p = NA_real_))
  stat <- as.numeric(t(U) %*% j$inv %*% U)
  p <- if (!is.null(C)) .ctScoreBoot(C, nboot) else
    stats::pchisq(stat, j$rank, lower.tail = FALSE)
  list(LM = stat, df = j$rank, p = p)
}

.ctScoreBH <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  out[ok] <- stats::p.adjust(p[ok], 'BH')
  out
}

# Each tested coordinate's variance over subjects, with its covariances with
# the effects in R: the score test of letting it vary as ctsem fits a new
# random effect, correlated with those already there. Per subject, adding a
# covariance Sigma over raw coordinates changes the log likelihood at zero by
# tr((s_i s_i' + H_i) Sigma) / 2 (Chesher, 1984), so the variance's
# contribution is (s_ij^2 + H_i,jj) / 2 and a covariance's s_ij s_ir + H_i,jr.
# Each is projected on the fitted parameters (C(alpha)): plus D' I^-1 s_i,
# with D its derivative in theta. Their variance, `robust`, is the spread of
# the projected contributions, which stays right when the model is wrong in
# some other way -- the case a screen meets most -- where the information
# equality the model-based variance leans on fails. `model` is a Gaussian
# score's, cov(M_ab, M_cd) = H_ac H_bd + H_ad H_bc summed, less D' I^-1 D;
# it also sets the one-step size.
#
# Why the covariances: every subject already has a random effect, their
# initial state. A CINT that varies independently of it, the initial means
# shared, is subjects drifting apart from a common start, and data whose
# subjects differ from the first occasion reject that: the book's FI-CLPM read
# z -5 on the variance alone, against a random intercept model 92 -2LL units
# better. The covariance with the initial state is a first-order parameter --
# a regression of the subject's intercept on where they start -- so it is
# tested as one, with the variance one-sided given it, and nothing is assumed
# about either: not stationarity, not a common start.
.ctScoreHeterogeneity <- function(pieces, ptab, items, variance = 'robust',
  nboot = 1999L) {
  S <- pieces$S
  n <- nrow(S)
  Iinv <- .ctScorePinv(-pieces$Hsum)$inv
  rows <- lapply(items, function(it) {
    R <- it$R
    Hr <- lapply(R, function(r) .ctScoreHrow(pieces, r))
    # D = dT/dtheta. Variance: d(s_ij^2/2) = s_ij H_i[j, ] and
    # d(H_i,jj/2)/dtheta_k = d^2 s_ik / dtheta_j^2 / 2. Covariance:
    # d(s_ij s_ir) = H_i[j, ] s_ir + s_ij H_i[r, ], and d H_i,jr / dtheta_k =
    # d^2 s_ik / dtheta_j dtheta_r. T(theta) + D' I^-1 S(theta) is T at the
    # true theta to first order, wherever it is evaluated.
    raw <- matrix(0.5 * (it$s^2 + it$hjj), n)
    D <- matrix(colSums(it$s * it$H) + 0.5 * it$D2, ncol = 1)
    for (k in seq_along(R)) {
      raw <- cbind(raw, it$s * S[, R[k]] + it$H[, R[k]])
      D <- cbind(D, colSums(it$H * S[, R[k]]) + colSums(it$s * Hr[[k]]) + it$Tmix[, k])
    }
    X <- raw + S %*% (Iinv %*% D)
    m <- ncol(X)
    info <- matrix(0, m, m)
    info[1, 1] <- sum(0.5 * it$hjj^2)
    for (k in seq_along(R)) {
      info[1, 1 + k] <- info[1 + k, 1] <- sum(it$hjj * it$H[, R[k]])
      for (l in seq_len(k)) info[1 + k, 1 + l] <- info[1 + l, 1 + k] <-
        sum(it$hjj * Hr[[k]][, R[l]] + it$H[, R[k]] * it$H[, R[l]])
    }
    info <- info - t(D) %*% Iinv %*% D
    info <- (info + t(info)) / 2
    # A covariance the fitted scores already determine says nothing and is
    # dropped; a variance that is determined leaves nothing to test.
    determined <- .ctScoreDetermined(X[, 1], raw[, 1])
    keep <- c(TRUE, vapply(seq_along(R), function(k)
      !.ctScoreDetermined(X[, 1 + k], raw[, 1 + k]), logical(1)))
    X <- X[, keep, drop = FALSE]
    info <- info[keep, keep, drop = FALSE]
    R <- R[keep[-1]]
    U <- colSums(X)
    V <- if (variance == 'robust') .ctScoreSpread(X) else info
    mix <- .ctScoreMixture(U, V)
    q <- ncol(X) - 1L
    # Where the new effect is independent of what varies, the covariances only
    # widen the test's reference: on the drift fixture of test-ctScoreTest.R,
    # p .043 against the variance alone's under .01, and a default model's
    # rows carry four. So the test is the stronger of the variance alone
    # (its covariances held at zero, z) and the joint test (LM), by their
    # smaller p, calibrated by resampling subjects as both are recomputed --
    # with the model variance, twice the smaller p.
    zalone <- function(U, V) if (V[1, 1] > 0) U[1] / sqrt(V[1, 1]) else NA_real_
    evidence <- function(U, V) {
      m <- .ctScoreMixture(U, V)
      max(-stats::pnorm(zalone(U, V), lower.tail = FALSE, log.p = TRUE),
        -log(.ctScoreMixtureP(m$LM, m$q)), na.rm = TRUE)
    }
    z <- zalone(U, V)
    p <- if (determined) NA_real_ else if (q == 0L) {
      if (variance == 'robust') .ctScoreBoot(X, nboot, onesided = TRUE) else
        stats::pnorm(z, lower.tail = FALSE)
    } else if (variance == 'robust') .ctScoreBoot(X, nboot, statistic = function(Y)
      evidence(colSums(Y), .ctScoreSpread(Y))) else
      min(1, 2 * min(stats::pnorm(z, lower.tail = FALSE), .ctScoreMixtureP(mix$LM, mix$q)))
    # The one-step variance, the covariances estimated beside it. Where the
    # projection leaves the information indefinite (its third derivatives are
    # differenced twice and can overshoot), the variance's own, unprojected.
    # Where the covariances carry the statistic the variance given them can
    # read zero or below, and that is no size: NA. Bounding it by what the
    # one-step covariances imply, k' P^-1 k, was tried and read 11 for an
    # intercept whose refit sd is under 1 (the old book data).
    ev <- eigen(info, symmetric = TRUE, only.values = TRUE)$values
    tau <- if (all(ev > 0)) solve(info, U)[1] else {
      i1 <- sum(0.5 * it$hjj^2)
      if (i1 > 0) U[1] / i1 else NA_real_
    }
    sdraw <- if (!is.finite(tau)) NA_real_ else if (tau > 0) sqrt(tau) else
      if (q > 0L) NA_real_ else 0
    if (determined) mix$LM <- z <- NA_real_
    row <- data.frame(parameter = it$name, matrix = it$matrix, row = it$row,
      col = it$col, LM = mix$LM, df = q + 1L, z = z, p = p, sd = sdraw,
      sd_natural = abs(it$slope) * sdraw,
      with = if (length(R)) paste(ptab$name[R], collapse = ', ') else NA_character_,
      stringsAsFactors = FALSE)
    # For likeness between tests, each subject's share of the statistic's
    # direction.
    attr(row, 'contribution') <- if (q == 0L) X[, 1] else
      as.numeric(X %*% (.ctScorePinv(V)$inv %*% U))
    row
  })
  contributions <- do.call(cbind, lapply(rows, attr, 'contribution'))
  out <- do.call(rbind, rows)
  out$p_adj <- .ctScoreBH(out$p)
  out$class <- .ctScoreClasses(contributions)
  out[order(out$p, -out$LM), , drop = FALSE]
}

# Whether the subjects vary along a direction no single parameter's test
# isolates -- a new dimension of person variation over several parameters, or
# on a reduced-rank level one more than its rank. Adding variance tau^2 along
# a raw direction u changes subject i's log likelihood, at tau^2 = 0, by
# (u' s_i)^2 / 2 + u' H_i u / 2: the expectation over the subject's posterior
# of the conditional score's square and Hessian is the marginal's square plus
# its Hessian. So T = sum_i (s_i s_i' + H_i) / 2 over the tested parameters is
# the score matrix of every such direction; its diagonal is the
# per-parameter test. The off-diagonal projections omit the third
# derivatives the diagonal has, which can only widen the spread and so keeps
# the test conservative. The statistic is the largest eigenvalue of T scaled
# by each diagonal's spread (invariant to each parameter's units), against a
# multiplier bootstrap over subjects; its eigenvector is the direction, read
# as loadings (Hall & Praestgaard, 2001; Zhu & Zhang, 2006). `sd` is the
# one-step sd along the direction, scaled to its largest loading of one: the
# score sum_i u'C_i u over its expected information, sum_i (u'H_i u)^2 / 2,
# as a new random effect's one-step variance is; NA where that is not
# positive.
.ctScoreDimension <- function(pieces, ptab, het, nboot = 1999L) {
  S <- pieces$S
  n <- nrow(S)
  Iinv <- .ctScoreIinv(pieces)
  Hr <- lapply(het, function(j) .ctScoreHrow(pieces, j))
  # Candidates whose own contributions the fit determines say nothing here
  # (.ctScoreDetermined), and their near-zero spread would scale the rest.
  live <- vapply(seq_along(het), function(a) {
    ja <- het[a]
    raw <- 0.5 * (S[, ja]^2 + Hr[[a]][, ja])
    D <- colSums(S[, ja] * Hr[[a]]) +
      0.5 * colSums(pieces$D2[, , match(ja, pieces$second), drop = TRUE])
    !.ctScoreDetermined(raw + as.numeric(S %*% (Iinv %*% D)), raw)
  }, logical(1))
  het <- het[live]; Hr <- Hr[live]
  K <- length(het)
  if (K < 2L) return(NULL)
  C <- array(0, c(n, K, K))
  for (a in seq_len(K)) for (b in a:K) {
    ja <- het[a]; jb <- het[b]
    base <- 0.5 * (S[, ja] * S[, jb] + Hr[[a]][, jb])
    D <- 0.5 * (colSums(S[, ja] * Hr[[b]]) + colSums(S[, jb] * Hr[[a]]))
    if (a == b) D <- colSums(S[, ja] * Hr[[a]]) +
      0.5 * colSums(pieces$D2[, , match(ja, pieces$second), drop = TRUE])
    C[, a, b] <- C[, b, a] <- base + as.numeric(S %*% (Iinv %*% D))
  }
  T <- apply(C, c(2, 3), sum)
  dsd <- vapply(seq_len(K), function(a) sqrt(.ctScoreSpread(matrix(C[, a, a]))[1, 1]), 0)
  if (any(!is.finite(dsd) | dsd <= 0)) return(NULL)
  e <- eigen(T / sqrt(outer(dsd, dsd)), symmetric = TRUE)
  lambda <- e$values[1]
  # The same statistic from any rows of the flattened contributions: each
  # diagonal's spread recomputed, so the reference is studentised as lambda is.
  diagcols <- (seq_len(K) - 1L) * K + seq_len(K)
  largest <- function(X) {
    d <- sqrt(pmax(diag(.ctScoreSpread(X[, diagcols, drop = FALSE])), 0))
    if (any(!is.finite(d) | d <= 0)) return(NA_real_)
    eigen(matrix(colSums(X), K) / sqrt(outer(d, d)), symmetric = TRUE,
      only.values = TRUE)$values[1]
  }
  u <- e$vectors[, 1] / sqrt(dsd)
  u <- u / u[which.max(abs(u))]
  score <- sum(vapply(seq_len(K), function(a) sum(u[a] * u * T[a, ]), 0))
  hu <- Reduce(`+`, lapply(seq_len(K), function(a)
    u[a] * as.numeric(Hr[[a]][, het, drop = FALSE] %*% u)))
  tau2 <- score / sum(hu^2 / 2)
  list(lambda = lambda, p = .ctScoreBoot(matrix(C, n, K * K), nboot, statistic = largest),
    direction = stats::setNames(u, ptab$name[het])[order(-abs(u))],
    sd = if (is.finite(tau2) && tau2 > 0) sqrt(tau2) else NA_real_)
}


# The subject level's blocks on a reduced-rank level (`popblocks`; without
# them the level is one block of every effect, its rank the level's): each
# block's name, rank, its effects and their raw population-mean coordinates
# (`index`), and its fitted loadings as directions in those coordinates
# (`dirs`, one column a dimension, each loading times the effect's sd scale
# as `_laplace_poploading` builds it). An effect with a stated mean has no
# coordinate in the fit (`stated` names them): given `ext`, the model with
# those means estimated (`.ctScoreStatedSpec()`), the coordinates are that
# model's and it has one. A dimension whose loadings are all zero gives no
# direction. NULL where the subject level is not of reduced rank.
.ctScoreBlockPlan <- function(fit, est, ext = NULL) {
  levels <- fit$model_spec$laplace$levels
  if (!length(levels)) return(NULL)
  lv <- levels[[1]]
  if (!length(lv$load_index) || !isTRUE(lv$rank > 0)) return(NULL)
  npar <- length(est)
  k <- length(lv$param)
  raw <- as.integer(lv$re_index)
  raw[raw > npar] <- NA_integer_
  if (!is.null(ext)) {
    free <- !is.na(raw)
    raw[free] <- ext$theta[raw[free]]
    raw[!free] <- ext$delta[match(as.character(lv$param)[!free],
      names(ext$x2)[ext$delta])]
    npar <- length(ext$x2)
  }
  entries <- .ctJuliaLoadingEntries(lv)
  L <- matrix(0, k, lv$rank)
  L[cbind(entries$row, entries$col)] <- est[lv$load_index] *
    as.numeric(lv$sd_scale)[entries$row]
  blocks <- if (length(lv$blocks)) lv$blocks else
    list(list(name = 'level', rank = lv$rank, param = lv$param))
  start <- 0L
  lapply(blocks, function(b) {
    cols <- start + seq_len(b$rank)
    start <<- start + b$rank
    members <- match(b$param, lv$param)
    free <- !is.na(raw[members])
    dirs <- matrix(0, npar, b$rank,
      dimnames = list(NULL, paste0(b$name, '.dim', seq_len(b$rank))))
    dirs[raw[members[free]], ] <- L[members[free], cols, drop = FALSE]
    dirs <- dirs[, colSums(abs(dirs)) > 0, drop = FALSE]
    list(name = b$name, rank = as.integer(b$rank), param = b$param,
      index = raw[members[free]], dirs = dirs,
      stated = as.character(lv$param)[members[!free]])
  })
}

# The inverse information a projection reads: the pieces' own where they carry
# one (`.ctScoreStatedPieces()`), else that of their summed Hessian.
.ctScoreIinv <- function(pieces) if (!is.null(pieces$Iinv)) pieces$Iinv else
  .ctScorePinv(-pieces$Hsum)$inv

# The fitted model with the stated population means of `stated` estimated, at
# their stated values (`spec$anchors$raw`), which is the fit: prepared as the
# fit was (`.ctScoreExtendedSpec()`), for the block tests of effects whose mean
# the fit states. The modes it takes are held for as long as `hold`, the
# caller's frame, runs.
.ctScoreStatedSpec <- function(fit, est, stated, cores = 2L, hold = parent.frame()) {
  base <- .ctFitBaseModel(fit)
  anchors <- fit$model_spec$anchors
  if (is.null(base) || !all(stated %in% anchors$param)) stop("the stated means of ",
    paste(stated, collapse = ', '), " are not in this fit's model.", call. = FALSE)
  model <- base
  keep <- setdiff(names(base$fixedpopmeans), stated)
  model$fixedpopmeans <- if (length(keep)) base$fixedpopmeans[keep] else NULL
  model$fixedpopmeansraw <- NULL
  .ctScoreExtendedSpec(fit, est, model, stated, 'effects with a stated mean', cores,
    hold = hold, at = stats::setNames(as.numeric(anchors$raw[match(stated,
      anchors$param)]), stated))
}

# The block tests' pieces on the model with stated means estimated (`ext`):
# each subject's scores there, along the block directions and the effects'
# coordinates (`.ctScorePieces()`), projected on the fit's own parameters
# alone -- the estimated means are held at their stated values under the
# null, as the fit holds them -- through the fit's information (`pieces`), or
# where the extension is read one-step or on the layer, its own along the
# fit's parameters.
.ctScoreStatedPieces <- function(ext, pieces, within, vdirs, vpairs, cores = 2L) {
  n2 <- length(ext$x2)
  info <- if (isTRUE(ext$onestep) || isTRUE(ext$layer))
    .ctScoreExtendedInformation(ext, cores) else -pieces$Hsum
  H <- matrix(0, n2, n2)
  H[ext$theta, ext$theta] <- -info
  score <- function(x) as.matrix(.ctJuliaGet(ext$module$ctsem_subject_gradients(
    ext$objective, .ctJuliaVector(as.numeric(x)), prior = TRUE))$scores)
  bp <- .ctScorePieces(NULL, as.numeric(ext$x2), rows = within, second = within,
    vdirs = vdirs, vpairs = vpairs, cores = cores, score = score, hessian = H)
  if (nrow(bp$S) != nrow(pieces$S)) stop("the extended model's units are not the fit's.",
    call. = FALSE)
  bp$Iinv <- matrix(0, n2, n2)
  bp$Iinv[ext$theta, ext$theta] <- .ctScorePinv(info)$inv
  bp
}

# The pairs of directions, one from each of two blocks, whose covariance the
# blocks hold at zero (`vdirs`' columns; `vblock` their blocks). A covariance
# between dimension `u` of block A and `w` of block B moves the covariance by
# u w' + w u', which is a move of B's loadings along `u` (with a turn of B's
# dimensions where its triangle fixes an entry) wherever `u` loads only on
# B's effects, and likewise the other way. Such a covariance is in the model
# already: its score at the estimate is the fitted loadings' score, zero, up
# to how far Laplace's term departs from the Gaussian identity the test reads
# it by. Tested, a style block over the thresholds and the dynamic effects
# gave LM 66 on 4 df against a dynamics block over those effects, from
# Laplace's error alone. So a pair of blocks one of whose effects contain the
# other's is not tested; the pairs left out are the attribute `nested`, as
# block indices.
.ctScoreBlockPairs <- function(plan, vdirs, vblock) {
  if (length(unique(vblock)) < 2L) return(NULL)
  g <- as.matrix(expand.grid(a = seq_along(vblock), b = seq_along(vblock)))
  g <- g[vblock[g[, 1]] < vblock[g[, 2]], , drop = FALSE]
  inside <- function(d, blk) all(which(vdirs[, d] != 0) %in% plan[[blk]]$index)
  inmodel <- vapply(seq_len(nrow(g)), function(r)
    inside(g[r, 1], vblock[g[r, 2]]) || inside(g[r, 2], vblock[g[r, 1]]), logical(1))
  out <- unname(g[!inmodel, , drop = FALSE])
  nested <- unique(cbind(vblock[g[inmodel, 1]], vblock[g[inmodel, 2]]))
  attr(out, 'nested') <- if (nrow(nested)) unname(nested) else NULL
  out
}

# What a reduced-rank level's blocks leave out, beside the per-parameter own
# variances: a further shared dimension within a block, and covariances between
# blocks.
#
# Within a block, the dimension test (`.ctScoreDimension()`) over its effects:
# whether subjects vary along a direction among them that the block's
# dimensions do not span, the fitted loadings projected out with the other
# parameters. Only where the block has more effects than dimensions; an effect
# whose own variance the fit estimates determines its diagonal and drops out.
#
# Between blocks A and B, the covariances of their dimensions: Sigma_AB = L_A
# C L_B', zero in the fit. Along `u` (a column of L_A) and `w` (one of L_B),
# a subject's log likelihood moves at C = 0 by u'(s_i s_i' + H_i) w, as for
# the per-parameter covariances (Chesher, 1984), projected on the fitted
# parameters with the mixed third derivative along `u` and `w`
# (`.ctScorePieces()`'s `Vmix`); the loadings' own dependence on the
# parameters enters with expectation zero and is left out. The r_A r_B
# covariances are tested jointly, two sided, by a studentised bootstrap of
# the subjects' contributions, as the covariates' joint tests are; with the
# model variance, against chi-square.
#
# Both tables' `p_adj` are adjusted together (Benjamini-Hochberg), as one
# family. `directions` keeps each tested block's whole direction, named by
# its effects' raw coordinates, for a refit that adds the dimension
# (ctModelImprove); the table shows its five largest loadings, and `sd` the
# one-step sd along it.
.ctScoreBlocks <- function(pieces, ptab, plan, vblock, variance = 'robust',
  nboot = 1999L, nested = NULL) {
  directions <- list()
  dimension <- do.call(rbind, lapply(plan, function(b) {
    if (length(b$index) <= b$rank) return(NULL)
    d <- .ctScoreDimension(pieces, ptab, b$index, nboot)
    if (!is.null(d)) directions[[b$name]] <<- d$direction
    top <- if (is.null(d)) character() else utils::head(names(d$direction), 5L)
    data.frame(block = b$name, rank = b$rank, effects = length(b$index),
      lambda = if (is.null(d)) NA_real_ else d$lambda,
      p = if (is.null(d)) NA_real_ else d$p,
      sd = if (is.null(d)) NA_real_ else d$sd,
      direction = if (is.null(d)) NA_character_ else
        paste0(top, ' ', signif(d$direction[top], 2), collapse = ', '),
      stringsAsFactors = FALSE)
  }))
  covariance <- NULL
  if (NROW(pieces$vpairs)) {
    S <- pieces$S
    n <- nrow(S)
    Iinv <- .ctScoreIinv(pieces)
    V <- pieces$vdirs
    Hd <- function(a) {
      M <- pieces$Hv[, a, , drop = TRUE]
      if (is.null(dim(M))) M <- matrix(M, n)
      M
    }
    sv <- S %*% V
    pairs <- pieces$vpairs
    key <- paste(vblock[pairs[, 1]], vblock[pairs[, 2]])
    covariance <- do.call(rbind, lapply(unique(key), function(kk) {
      idx <- which(key == kk)
      ab <- pairs[idx, , drop = FALSE]
      blocks <- paste(plan[[vblock[ab[1, 1]]]]$name, plan[[vblock[ab[1, 2]]]]$name,
        sep = ' x ')
      Ha <- lapply(ab[, 1], Hd)
      Hb <- lapply(ab[, 2], Hd)
      raw <- matrix(0, n, length(idx))
      D <- matrix(0, ncol(S), length(idx))
      for (k in seq_along(idx)) {
        a <- ab[k, 1]; b <- ab[k, 2]
        raw[, k] <- sv[, a] * sv[, b] + as.numeric(Ha[[k]] %*% V[, b])
        D[, k] <- colSums(Ha[[k]] * sv[, b]) + colSums(Hb[[k]] * sv[, a]) +
          pieces$Vmix[, idx[k]]
      }
      X <- raw + S %*% (Iinv %*% D)
      keep <- vapply(seq_along(idx), function(k) !.ctScoreDetermined(X[, k], raw[, k]),
        logical(1))
      if (!any(keep)) return(data.frame(blocks = blocks, LM = NA_real_, df = 0L,
        p = NA_real_, stringsAsFactors = FALSE))
      kept <- which(keep)
      X <- X[, kept, drop = FALSE]
      j <- if (variance == 'robust') .ctScoreJoint(colSums(X), .ctScoreSpread(X), X, nboot) else {
        # cov(M_ab, M_cd) = H_ac H_bd + H_ad H_bc for a Gaussian score, along
        # the directions, summed over subjects, less what estimation took.
        info <- matrix(0, length(kept), length(kept))
        for (k in seq_along(kept)) for (l in seq_len(k)) {
          a2 <- ab[kept[l], 1]; b2 <- ab[kept[l], 2]
          info[k, l] <- info[l, k] <- sum(
            as.numeric(Ha[[kept[k]]] %*% V[, a2]) * as.numeric(Hb[[kept[k]]] %*% V[, b2]) +
            as.numeric(Ha[[kept[k]]] %*% V[, b2]) * as.numeric(Hb[[kept[k]]] %*% V[, a2]))
        }
        Dk <- D[, kept, drop = FALSE]
        info <- info - t(Dk) %*% Iinv %*% Dk
        .ctScoreJoint(colSums(X), (info + t(info)) / 2, NULL, nboot)
      }
      data.frame(blocks = blocks, LM = j$LM, df = j$df, p = j$p, stringsAsFactors = FALSE)
    }))
  }
  # Pairs of blocks whose covariances are in the model already
  # (`.ctScoreBlockPairs()`): said, not tested.
  if (NROW(nested)) covariance <- rbind(if (!is.null(covariance))
    cbind(covariance, note = NA_character_, stringsAsFactors = FALSE),
    do.call(rbind, lapply(seq_len(nrow(nested)), function(r) {
      a <- plan[[nested[r, 1]]]$name; b <- plan[[nested[r, 2]]]$name
      data.frame(blocks = paste(a, b, sep = ' x '), LM = NA_real_, df = 0L,
        p = NA_real_, note = 'in the model: one block loads on every effect of the other',
        stringsAsFactors = FALSE)
    })))
  adj <- .ctScoreBH(c(dimension$p, covariance$p))
  if (NROW(dimension)) dimension$p_adj <- adj[seq_len(nrow(dimension))]
  if (NROW(covariance)) covariance$p_adj <- adj[NROW(dimension) + seq_len(nrow(covariance))]
  list(dimension = dimension, covariance = covariance, directions = directions)
}

# Whether a projected contribution is one the fitted scores already determine:
# under one per cent of its spread left once estimation is taken out (the
# guard the adequacy families use). On a reduced-rank level the variance of a
# parameter that loads on the shared dimensions is mostly estimated already,
# and an intercept's heterogeneity contributions kept 0.6 per cent; divided by
# what was left, its noise dominated the dimension test (rank 1 against 2
# found 4 per cent of the time where its parameters' own tests found 67).
.ctScoreDetermined <- function(projected, raw) {
  v <- stats::var(as.numeric(projected)); v0 <- stats::var(as.numeric(raw))
  !is.finite(v) || !is.finite(v0) || v0 <= 0 || v < 0.01 * v0
}

# Subject i's Hessian row for raw parameter j, n x p.
.ctScoreHrow <- function(pieces, j) {
  Hj <- pieces$H[, match(j, pieces$rows), , drop = TRUE]
  if (is.null(dim(Hj))) Hj <- matrix(Hj, nrow = nrow(pieces$S))
  Hj
}

# The variance of a sum of independent units' contributions (rows of C),
# estimated from their spread: n / (n - 1) times the centred cross product.
.ctScoreSpread <- function(C) {
  C <- as.matrix(C)
  n <- nrow(C)
  centred <- sweep(C, 2, colMeans(C))
  crossprod(centred) * n / (n - 1)
}

# The p-value of U' V^-1 U, U = colSums(C) and V its robust variance, by a
# studentised bootstrap over the units: rows of C centred (the null) and
# resampled, the statistic recomputed with its own variance each time. The
# chi-square (or normal) reference assumes U is Gaussian, and a variance
# score is a sum of chi-square-like terms: with 100 subjects the robust
# one-sided z rejected a true model at 2 per cent where 5 was asked, and the
# joint covariate test at 7.5 to 10.5. Studentising and resampling corrects
# the skewness to second order, and costs milliseconds: nothing is refitted.
# `onesided` takes a single column and rejects for large U. `statistic`
# replaces U' V^-1 U with another function of the contributions, recomputed the
# same way. `C` may be a list of contribution matrices over the same units, for
# the largest of their statistics -- a search, as over a threshold's cuts --
# against the largest in each resample.
#
# The draws must be studentised as the observed statistic is. A Gaussian
# multiplier with the observed variance held fixed is not: for the dimension
# test, whose diagonal sums are right-skewed and so left-skewed once
# studentised, it made the reference too wide (0 of 100 correct models
# rejected, and 2 per cent power for a missing dimension). A search over a
# covariate's cuts is the exception, for which see .ctScoreThreshold.
.ctScoreBoot <- function(C, nboot = 1999L, onesided = FALSE, statistic = NULL) {
  Cs <- if (is.list(C) && !is.data.frame(C)) lapply(C, as.matrix) else list(as.matrix(C))
  n <- nrow(Cs[[1]])
  stat <- if (!is.null(statistic)) statistic else function(X) {
    u <- colSums(X)
    v <- .ctScoreSpread(X)
    if (onesided) return(u / sqrt(v[1, 1]))
    pv <- .ctScorePinv(v)
    as.numeric(t(u) %*% pv$inv %*% u)
  }
  sup <- function(Xs) max(vapply(Xs, stat, numeric(1)))
  observed <- sup(Cs)
  if (!is.finite(observed)) return(NA_real_)
  nulls <- lapply(Cs, function(X) sweep(X, 2, colMeans(X)))
  draws <- vapply(seq_len(nboot), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    sup(lapply(nulls, function(X) X[idx, , drop = FALSE]))
  }, numeric(1))
  (1 + sum(draws >= observed, na.rm = TRUE)) / (1 + sum(is.finite(draws)))
}

# Covariate values per subject, in the fit's subject order, as a named list.
.ctScoreCovariateValues <- function(fit, covariates, own = FALSE) {
  covariates <- as.data.frame(covariates, stringsAsFactors = FALSE)
  idname <- .ctFitModelObject(fit)$subjectIDname
  ids <- .ctFitIdMap(fit)$original
  idcol <- intersect(c(idname, 'id'), names(covariates))[1]
  vars <- setdiff(names(covariates), idcol)
  if (!length(vars)) stop("covariates has no columns besides the id.", call. = FALSE)
  inmodel <- intersect(vars, .ctFitModelObject(fit)$TIpredNames)
  if (length(inmodel) && !own) stop("Already in the model as time independent ",
    "predictors: ", paste(inmodel, collapse = ', '), ". Those are tested on ",
    "the parameters they do not act on without being given; give covariates ",
    "the model does not contain.", call. = FALSE)
  if (is.na(idcol)) {
    if (nrow(covariates) != length(ids)) stop("covariates needs a '", idname,
      "' column, or exactly one row per subject in the fit's subject order (",
      length(ids), ").", call. = FALSE)
    return(as.list(covariates[vars]))
  }
  key <- as.character(covariates[[idcol]])
  out <- lapply(vars, function(v) {
    x <- covariates[[v]]
    first <- match(as.character(ids), key)
    if (anyDuplicated(key)) {
      varies <- tapply(seq_along(x), key, function(r) length(unique(x[r][!is.na(x[r])])) > 1L)
      if (any(varies)) stop("Covariate '", v, "' varies within subject. Only ",
        "time independent covariates are tested here.", call. = FALSE)
      firstvalue <- tapply(seq_along(x), key, function(r) {
        r <- r[!is.na(x[r])]
        if (length(r)) r[1] else NA_integer_
      })
      first <- as.integer(firstvalue[as.character(ids)])
    }
    x[first]
  })
  names(out) <- vars
  missingids <- sum(!as.character(ids) %in% key)
  if (missingids && !own) message(missingids, " subjects have no row in ",
    "covariates and contribute at the covariates' means.")
  out
}

# The model's own time independent predictors that the covariates family
# tests: every one but the columns of a smooth term (`s()`), whose shape is
# penalised and whose effects come as a set.
.ctScoreOwnNames <- function(fit) {
  model <- .ctFitModelObject(fit)
  setdiff(model$TIpredNames, unlist(lapply(model$TIsmooths, `[[`, 'columns')))
}

# Each of the model's own predictors (`.ctScoreOwnNames()`) with the
# parameters among `cov` it does not act on (`targets`, raw positions of model
# parameters; the specification's `ti_effects` says what it acts on) and its
# value per subject as the fit's data give it. Effects on a population sd
# are left out: a predictor acts on a parameter's value. Only predictors with
# a target are returned.
.ctScoreOwnPredictors <- function(fit, ptab, cov) {
  names <- .ctScoreOwnNames(fit)
  if (!length(names)) return(NULL)
  model <- .ctFitModelObject(fit)
  te <- fit$model_spec$ti_effects
  candidates <- cov[ptab$kind[cov] == 'model']
  targets <- lapply(names, function(z) {
    acts <- if (NROW(te)) as.integer(te$parameter[as.integer(te$predictor) ==
      match(z, model$TIpredNames)]) else integer()
    setdiff(candidates, acts)
  })
  keep <- lengths(targets) > 0L
  if (!any(keep)) return(NULL)
  names <- names[keep]
  data <- as.data.frame(.ctFitLongData(fit))
  values <- .ctScoreCovariateValues(fit, data[, c(model$subjectIDname, names),
    drop = FALSE], own = TRUE)
  stats::setNames(Map(function(v, t) list(values = v, targets = t), values[names],
    targets[keep]), names)
}

# The basis a covariate enters with, centred, missing values at the mean
# (zero after centring) -- any fixed function of the covariates is a valid
# weight, so a subject without a value contributes nothing and biases nothing.
.ctScoreBasis <- function(z, form, df, name) {
  isfactor <- is.factor(z) || is.character(z) || is.logical(z)
  if (form == 'auto') form <- if (isfactor) 'factor' else 'linear'
  if (form == 'factor') {
    f <- droplevels(as.factor(z))
    if (nlevels(f) < 2L) stop("Covariate '", name, "' has fewer than two levels.",
      call. = FALSE)
    X <- matrix(0, length(z), nlevels(f) - 1L)
    for (l in seq_len(ncol(X))) X[, l] <- as.numeric(f == levels(f)[l + 1L])
    colnames(X) <- levels(f)[-1L]
  } else {
    if (isfactor) stop("Covariate '", name, "' is not numeric; use form = ",
      "'factor'.", call. = FALSE)
    z <- as.numeric(z)
    ok <- !is.na(z)
    if (length(unique(z[ok])) < 2L) stop("Covariate '", name, "' does not vary.",
      call. = FALSE)
    if (form == 'linear') {
      X <- matrix(z, ncol = 1L, dimnames = list(NULL, 'linear'))
    } else {
      X <- matrix(NA_real_, length(z), df)
      X[ok, ] <- splines::ns(z[ok], df = df)
      colnames(X) <- paste0('ns', seq_len(df))
    }
  }
  ok <- stats::complete.cases(X)
  X[ok, ] <- sweep(X[ok, , drop = FALSE], 2, colMeans(X[ok, , drop = FALSE]))
  X[!ok, ] <- 0
  attr(X, 'name') <- name
  attr(X, 'form') <- form
  attr(X, 'nmissing') <- sum(!ok)
  X
}

# The C(alpha) test of a covariate's effects on the parameters `cov`, with the
# basis G: each subject's contribution g_i s_ij, projected on the fitted
# parameters with the per-subject Hessians weighted by G. The variance is the
# contributions' spread (`robust`) or the projected information (`model`).
#
# The effects are estimated with the fitted parameters held (Idd^-1 U*), not
# free (V^-1 U*, the classical expected parameter change). A covariate the
# model lacks is absorbed by the parameters it has -- a random effect's sd
# grows to take in a covariate that shifts its mean -- and the free one-step
# extrapolates the likelihood's curvature across that whole change. On
# simulated data where a covariate shifted a random manifest mean by 0.8:
# held 0.81, refit 0.81, free 3.0, and Newton from the fit diverged; over 100
# such datasets the free one also attributed effects to parameters the
# covariate did not touch (T0MEANS 50 per cent, a population sd 58). Holding
# them is bigIRT's post-hoc adjustment (its items held, the ability side
# estimated), here joint over all of the covariate's effects.
.ctScoreCovariate <- function(pieces, items, G, variance = 'robust',
  nboot = 1999L) {
  S <- pieces$S
  q <- ncol(G)
  k <- length(items)
  Iinv <- .ctScorePinv(-pieces$Hsum)$inv
  block <- function(a) (a - 1L) * q + seq_len(q)
  # Contributions g_ia s_ij, ordered by item then basis column.
  C <- do.call(cbind, lapply(items, function(it) G * it$s))
  U <- colSums(C)
  Idt <- matrix(0, q * k, ncol(S))
  Idd <- matrix(0, q * k, q * k)
  for (a in seq_len(k)) {
    Ha <- items[[a]]$H
    Idt[block(a), ] <- -crossprod(G, Ha)
    for (b in seq_len(k)) Idd[block(a), block(b)] <-
      -crossprod(G, G * as.numeric(Ha %*% items[[b]]$dir))
  }
  Idd <- (Idd + t(Idd)) / 2
  A <- Idt %*% Iinv
  Cstar <- C - S %*% t(A)
  Ustar <- colSums(Cstar)
  Vmodel <- Idd - A %*% t(Idt)
  Vmodel <- (Vmodel + t(Vmodel)) / 2
  V <- if (variance == 'robust') .ctScoreSpread(Cstar) else Vmodel
  joint <- .ctScoreJoint(Ustar, V, if (variance == 'robust') Cstar, nboot)
  heldinv <- .ctScorePinv(Idd)$inv
  held <- as.numeric(heldinv %*% Ustar)
  heldcov <- heldinv %*% V %*% heldinv
  # Each subject's contribution to the held estimates, for their bootstrap:
  # summed they are the estimates, and their spread is `heldcov`.
  heldC <- Cstar %*% t(heldinv)
  wald <- function(b, est, covariance) {
    if (variance == 'robust') return(.ctScoreBoot(heldC[, b, drop = FALSE], nboot))
    w <- .ctScorePinv(covariance[b, b, drop = FALSE])
    if (!w$rank) return(NA_real_)
    stats::pchisq(as.numeric(t(est[b]) %*% w$inv %*% est[b]), w$rank, lower.tail = FALSE)
  }
  rows <- lapply(seq_len(k), function(a) {
    b <- block(a)
    it <- items[[a]]
    one <- function(x) if (q == 1L) x else NA_real_
    data.frame(parameter = it$name, matrix = it$matrix,
      estimate = one(held[b]), se = one(sqrt(heldcov[b, b])),
      estimate_natural = it$slope * one(held[b]),
      p = wald(b, held, heldcov), stringsAsFactors = FALSE)
  })
  parameters <- do.call(rbind, rows)
  parameters$p_adj <- .ctScoreBH(parameters$p)
  if (q == 1L) parameters$class <- .ctScoreClasses(heldC)
  list(joint = data.frame(covariate = attr(G, 'name'), form = attr(G, 'form'),
    basis = q, LM = joint$LM, df = joint$df, p = joint$p,
    nmissing = attr(G, 'nmissing'), stringsAsFactors = FALSE),
    parameters = parameters[order(parameters$p, -abs(parameters$estimate /
      parameters$se)), , drop = FALSE])
}

# A covariate that acts above or below a threshold, at a threshold not known
# in advance: the joint test of a step at each candidate cut (the 15th to the
# 85th percentile), and the largest of them (Andrews, 1993; the split a
# score-guided tree makes, Arnold et al., 2021, 2024), against the same scan
# with the covariate permuted over subjects (maximally selected statistics,
# Hothorn & Lausen, 2003). Under the null the covariate is unrelated to the
# subject's data, so the permutation null is exact whatever the statistic's
# distribution, and the fit does not involve the covariate, so nothing is
# refitted. A bootstrap failed both ways: with each cut's variance held fixed
# a correct model was rejected 21 per cent of the time at .05, because at an
# extreme cut the variance comes from the few subjects on one side; resampled
# with that variance recomputed, the reference was so wide it found nothing.
# The effects reported are those of being above the best cut, held as for a
# linear covariate.
.ctScoreThreshold <- function(pieces, items, z, name, variance = 'robust',
  nboot = 1999L, ncuts = 15L) {
  if (is.factor(z) || is.character(z) || is.logical(z)) stop("Covariate '", name,
    "' is not numeric; a threshold needs an ordered covariate.", call. = FALSE)
  z <- as.numeric(z)
  ok <- !is.na(z)
  cuts <- unique(stats::quantile(z[ok], seq(0.15, 0.85, length.out = ncuts),
    names = FALSE, type = 1))
  cuts <- cuts[vapply(cuts, function(c) {
    a <- sum(z[ok] > c); a >= 5 && sum(ok) - a >= 5 }, logical(1))]
  if (!length(cuts)) stop("Covariate '", name, "' has too few distinct values ",
    "for a threshold.", call. = FALSE)
  # The step at each cut, centred over the subjects who have the covariate;
  # those missing it contribute at the mean.
  steps <- function(zz) lapply(cuts, function(c) {
    G <- matrix(0, length(zz), 1L, dimnames = list(NULL, 'step'))
    G[ok, 1] <- as.numeric(zz[ok] > c)
    G[ok, 1] <- G[ok, 1] - mean(G[ok, 1])
    G
  })
  Iinv <- .ctScorePinv(-pieces$Hsum)$inv
  stat <- function(C) {
    v <- .ctScorePinv(.ctScoreSpread(C))
    u <- colSums(C)
    as.numeric(t(u) %*% v$inv %*% u)
  }
  scan <- function(zz) vapply(steps(zz), function(G)
    stat(.ctScoreCovariateContributions(pieces, items, G, Iinv)), 0)
  lms <- scan(z)
  best <- which.max(lms)
  G <- steps(z)[[best]]
  attr(G, 'name') <- name; attr(G, 'form') <- 'threshold'
  attr(G, 'nmissing') <- sum(!ok)
  out <- .ctScoreCovariate(pieces, items, G, variance, nboot)
  perms <- vapply(seq_len(nboot), function(b) {
    zz <- z
    zz[ok] <- z[ok][sample.int(sum(ok))]
    max(scan(zz))
  }, 0)
  out$joint$LM <- lms[best]
  out$joint$p <- (1 + sum(perms >= lms[best])) / (1 + nboot)
  out$joint$cut <- cuts[best]
  out
}

# Each subject's projected contribution to a covariate's effects: the
# columns whose sum is the joint test's U*.
.ctScoreCovariateContributions <- function(pieces, items, G, Iinv = NULL) {
  S <- pieces$S
  q <- ncol(G); k <- length(items)
  if (is.null(Iinv)) Iinv <- .ctScorePinv(-pieces$Hsum)$inv
  C <- do.call(cbind, lapply(items, function(it) G * it$s))
  Idt <- do.call(rbind, lapply(items, function(it) -crossprod(G, it$H)))
  C - S %*% t(Idt %*% Iinv)
}

# Relaxations the data cannot tell apart at first order: those whose projected
# per-subject contributions are nearly collinear (|correlation| >= 0.9) would
# each explain the same signal, as two models that are equivalent locally
# (MacCallum et al., 1993). Connected groups get a letter; a relaxation alike
# to no other gets NA.
.ctScoreClasses <- function(C, threshold = 0.9) {
  C <- as.matrix(C)
  k <- ncol(C)
  out <- rep(NA_character_, k)
  if (k < 2L) return(out)
  sdv <- apply(C, 2, stats::sd)
  r <- suppressWarnings(stats::cor(C))
  r[!is.finite(r)] <- 0
  adj <- abs(r) >= threshold & outer(sdv > 0, sdv > 0)
  diag(adj) <- FALSE
  group <- rep(0L, k); g <- 0L
  for (a in seq_len(k)) {
    if (group[a] || !any(adj[a, ])) next
    g <- g + 1L
    todo <- a
    while (length(todo)) {
      b <- todo[1]; todo <- todo[-1]
      if (group[b]) next
      group[b] <- g
      todo <- c(todo, which(adj[b, ] & !group))
    }
  }
  out[group > 0] <- LETTERS[((group[group > 0] - 1L) %% 26L) + 1L]
  out
}

# Whether parameters are stable across the data, ordered by time and by each
# latent's predicted state: Nyblom's (1989) test from each data row's score.
# Under the model a subject's row scores are a martingale difference
# sequence, so their partial sums in any order fixed in advance -- by time, or
# by the state predicted from earlier rows -- behave as a Brownian motion; at
# the estimate they are pulled to zero at the end, a bridge. A parameter that
# drifts over the study, or differs between low and high levels of a state,
# bends its partial sums away from the bridge. L_j is the mean squared partial
# sum of parameter j's score over its total, the Cramer-von Mises form.
#
# The estimation effect is removed by projecting the partial sums with the
# information accumulated up to each row, H_t H_N^-1 S_N, with H the rows'
# outer products -- not with the share of rows, which is right only when the
# information accrues evenly along the ordering, and along a state it does
# not (a drift is informed by large states). That needs every parameter's
# row scores; above `maxall` parameters only the tested ones are used, and the
# others' estimation is ignored. Calibration: a wild bootstrap with a normal
# weight per row, which keeps the rows' variance and nothing else -- what the
# null has, since the scores are uncorrelated over rows (Hansen 1992; Zeileis
# & Hornik 2007 for the fluctuation-test family).
.ctScoreInstability <- function(fit, est, ptab, parameters = NULL, nboot = 1999L,
  cores = 2L, maxall = 48L) {
  tested <- which(ptab$kind == 'model' & !ptab$matrix %in% c('T0MEANS', 'T0VAR'))
  if (!is.null(parameters)) tested <- intersect(tested, which(ptab$name %in% parameters))
  if (!length(tested)) stop("No parameters to test for stability.", call. = FALSE)
  p <- length(est)
  all <- if (p <= maxall) seq_len(p) else tested
  R <- .ctBackendRowScores(fit, est, diag(p)[, all, drop = FALSE], cores)
  err <- .ctAdequacyErrors(fit, cores = cores)
  okrow <- apply(is.finite(R), 1, all)
  R[!okrow, ] <- 0
  Tc <- match(tested, all)
  k <- length(Tc)
  N <- sum(okrow)
  latent <- .ctFitModelObject(fit)$latentNames
  nl <- min(ncol(err$eta), length(latent))
  keys <- c(list(time = err$time), stats::setNames(lapply(seq_len(nl), function(j)
    err$eta[, j]), paste0('state ', latent[seq_len(nl)])))
  rows <- list(); joint <- list()
  for (o in names(keys)) {
    # Ties broken at random: the data's own order within a tie is the
    # subjects', which a statistic over the ordering would read as structure.
    idx <- order(keys[[o]], stats::runif(length(keys[[o]])))
    idx <- idx[okrow[idx]]
    s <- R[idx, , drop = FALSE]
    H <- crossprod(s)
    Hinv <- .ctScorePinv(H)$inv
    VT <- H[Tc, Tc, drop = FALSE]
    VTinv <- .ctScorePinv(VT)$inv
    vd <- diag(VT)
    sT <- s[, Tc, drop = FALSE]
    bridge <- function(sm) {
      B <- Hinv %*% colSums(sm)
      matrix(apply(sm[, Tc, drop = FALSE] - sT * as.numeric(s %*% B), 2, cumsum), ncol = k)
    }
    statistic <- function(St) c(colSums(St^2) / (N * vd), sum((St %*% VTinv) * St) / N)
    observed <- statistic(bridge(s))
    draws <- vapply(seq_len(nboot), function(b) statistic(bridge(s * stats::rnorm(nrow(s)))),
      observed)
    draws <- matrix(draws, nrow = length(observed))
    pv <- (1 + rowSums(draws >= observed)) / (1 + nboot)
    one <- data.frame(ordering = o, parameter = ptab$name[tested], L = observed[seq_len(k)],
      p = pv[seq_len(k)], stringsAsFactors = FALSE)
    one$p_adj <- .ctScoreBH(one$p)
    rows[[o]] <- one
    joint[[o]] <- data.frame(ordering = o, L = observed[k + 1L], df = k, p = pv[k + 1L],
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  list(parameters = out[order(out$p, -out$L), , drop = FALSE],
    joint = do.call(rbind, joint),
    note = if (p > maxall) paste0("Only the ", k, " tested parameters' row scores were ",
      "used; the other parameters' estimation is ignored, which can move the ",
      "level along a state.") else NULL)
}

# The fitted model with a term `d * eta_k * eta_m` in each latent equation's
# CINT, for k <= m (and degree three with `order = 3`), every d a new PARS
# parameter with an identity transform.
.ctScoreStateModel <- function(model, order = 2L) {
  latent <- model$latentNames
  nl <- length(latent)
  # The monomials, each as its sorted state indices: the first is the state
  # the drift coefficient multiplies, the rest what it varies with.
  monomials <- list()
  for (k in seq_len(nl)) for (m in k:nl) {
    monomials[[length(monomials) + 1L]] <- c(k, m)
    if (order >= 3L) for (l in m:nl) monomials[[length(monomials) + 1L]] <- c(k, m, l)
  }
  terms <- list()
  # Equations of latents with innovations of their own: a cycle's clock turns
  # as written, and a term in its drift kept one search's clock wandering. The
  # clock's states stay among the monomials: a drift varying with the time of
  # day is a question to ask. Nor the equations of latents a 'stationary'
  # initial condition needs fixed within each subject (`.ctStationaryClosure`):
  # a term there leaves no one stationary distribution, and the extended model
  # is refused -- which on the Fisher search took these families out of every
  # screen after the first such process was kept.
  equations <- setdiff(.ctImproveNoisyLatents(model), which(.ctStationaryClosure(model)))
  for (i in equations) for (mono in monomials) {
    terms[[length(terms) + 1L]] <- data.frame(name = paste0('ctscoreterm',
      length(terms) + 1L), equation = latent[i], i = i, k = mono[1L],
      by = paste(latent[mono[-1L]], collapse = '*'),
      term = paste(latent[mono], collapse = '*'),
      other = if (length(mono) == 2L && mono[1L] != mono[2L]) mono[2L] else NA_integer_,
      stringsAsFactors = FALSE)
  }
  terms <- do.call(rbind, terms)
  model <- .ctAddCintTerms(model, terms$i, terms$name,
    lapply(strsplit(terms$term, '*', fixed = TRUE), match, latent))
  list(model = model, terms = terms)
}

# Add terms `name * eta_k * eta_m ...` to CINT[i], each `name` a new PARS
# parameter (identity transform unless given), for the score tests at zero and
# for a refit that keeps one. A CINT cell that is a free parameter moves into
# PARS -- transform, variation and predictor effects and all -- so its raw
# coordinate is unchanged; a fixed cell keeps its value; an expression keeps its
# form. The engine's Jacobian is derived from DRIFT * state + CINT
# (`ctJacobian`), so a term here is the same model as the drift coefficient
# varying with the state.
.ctAddCintTerms <- function(model, i, names, monomials, transform = 'param') {
  pars <- model$pars
  latent <- model$latentNames
  effects <- grep('_effect$', names(pars), value = TRUE)
  levelflags <- grep('^indvarying_', names(pars), value = TRUE)
  nextpars <- max(c(0L, as.integer(pars$row[pars$matrix %in% 'PARS'])))
  template <- pars[which(pars$matrix %in% 'CINT')[1L], , drop = FALSE]
  parsrow <- function(name, from = NULL) {
    r <- if (is.null(from)) template else from
    r$matrix <- 'PARS'; r$row <- nextpars <<- nextpars + 1L; r$col <- 1L
    r$param <- name
    if (is.null(from)) {
      r$value <- NA; r$transform <- transform; r$indvarying <- FALSE
      if ('sdscale' %in% names(r)) r$sdscale <- 1
      for (e in effects) r[[e]] <- 'FALSE'
      for (l in levelflags) r[[l]] <- FALSE
    }
    r
  }
  added <- list(); moved <- list()
  for (eq in unique(i)) {
    cell <- which(pars$matrix %in% 'CINT' & pars$row == eq)
    if (length(cell) != 1L) stop("The model has no CINT cell for latent ", eq,
      ".", call. = FALSE)
    these <- which(i == eq)
    sum_i <- paste(vapply(these, function(t) paste(c(names[t], latent[monomials[[t]]]),
      collapse = ' * '), character(1)), collapse = ' + ')
    for (t in these) added[[length(added) + 1L]] <- parsrow(names[t])
    current <- pars$param[cell]
    if (is.na(current) || !nzchar(current)) {
      value <- suppressWarnings(as.numeric(pars$value[cell]))
      expr <- if (is.na(value) || value == 0) sum_i else paste0(value, ' + ', sum_i)
    } else if (grepl('^[A-Za-z.][A-Za-z0-9._]*$', current) && !current %in% latent) {
      moved[[as.character(cell)]] <- parsrow(current, from = pars[cell, , drop = FALSE])
      expr <- paste0(current, ' + ', sum_i)
    } else {
      expr <- paste0('(', current, ') + ', sum_i)
    }
    pars$param[cell] <- expr
    pars$value[cell] <- NA
    pars$transform[cell] <- NA
    pars$indvarying[cell] <- FALSE
    for (e in effects) pars[[e]][cell] <- 'FALSE'
    for (l in levelflags) pars[[l]][cell] <- FALSE
  }
  model$pars <- .ctInsertMovedRows(pars, moved, added)
  model
}

# `pars` with each row of `moved` -- a cell's own parameter taken into PARS,
# named by that cell's row -- placed straight after the cell, and the `added`
# rows at the end. Random effects are ordered as their parameters are
# declared (`.ctIndvaryingIndex`), so a moved parameter keeps its place among
# them and the extended model's correlation coordinates stay the fit's.
.ctInsertMovedRows <- function(pars, moved, added) {
  at <- as.integer(names(moved))
  out <- rbind(pars, do.call(rbind, unname(moved)), do.call(rbind, added))
  out <- out[order(c(seq_len(nrow(pars)), at + 0.5, nrow(pars) + seq_along(added))), ,
    drop = FALSE]
  rownames(out) <- NULL
  out
}

# The fitted model with each latent's innovation sd scaled by exp(sum over
# latents m of v * eta_m): a term `v` per latent with innovations and per state
# it may depend on, every v a new PARS parameter, identity transform. At zero
# the model is the fit.
.ctScoreVolatilityModel <- function(model) {
  latent <- model$latentNames
  noisy <- .ctImproveNoisyLatents(model)
  # Not the innovations of a latent a 'stationary' initial condition needs
  # fixed within each subject (`.ctScoreStateModel`); any noisy state may
  # moderate the others'.
  equations <- setdiff(noisy, which(.ctStationaryClosure(model)))
  if (!length(equations)) return(list(model = model, terms = NULL))
  grid <- expand.grid(m = noisy, i = equations)
  terms <- data.frame(name = paste0('ctscorevol', seq_len(nrow(grid))),
    equation = latent[grid$i], i = grid$i, m = grid$m, by = latent[grid$m],
    stringsAsFactors = FALSE)
  list(model = .ctAddDiffusionTerms(model, terms$i, terms$name, terms$m), terms = terms)
}

# Scale DIFFUSION[i, i], the sd of latent i's innovations, by exp(name *
# eta_m) for each term, each `name` a new PARS parameter (identity transform
# unless given). A diagonal cell that is a free parameter moves into PARS --
# transform, variation and predictor effects and all -- so its raw coordinate
# and value are unchanged; a fixed cell keeps its value; an expression keeps
# its form. Off-diagonal cells, correlations under the default covariance
# transform, are untouched, so the innovations' correlations stay as they were.
.ctAddDiffusionTerms <- function(model, i, names, m, transform = 'param') {
  pars <- model$pars
  latent <- model$latentNames
  effects <- grep('_effect$', names(pars), value = TRUE)
  levelflags <- grep('^indvarying_', names(pars), value = TRUE)
  nextpars <- max(c(0L, as.integer(pars$row[pars$matrix %in% 'PARS'])))
  template <- pars[which(pars$matrix %in% 'DIFFUSION')[1L], , drop = FALSE]
  parsrow <- function(name, from = NULL) {
    r <- if (is.null(from)) template else from
    r$matrix <- 'PARS'; r$row <- nextpars <<- nextpars + 1L; r$col <- 1L
    r$param <- name
    if (is.null(from)) {
      r$value <- NA; r$transform <- transform; r$indvarying <- FALSE
      if ('sdscale' %in% names(r)) r$sdscale <- 1
      for (e in effects) r[[e]] <- 'FALSE'
      for (l in levelflags) r[[l]] <- FALSE
    }
    r
  }
  added <- list(); moved <- list()
  for (eq in unique(i)) {
    cell <- which(pars$matrix %in% 'DIFFUSION' & pars$row == eq & pars$col == eq)
    if (length(cell) != 1L) stop("The model has no DIFFUSION cell for latent ", eq,
      ".", call. = FALSE)
    these <- which(i == eq)
    scale <- paste0('exp(', paste(vapply(these, function(t) paste(names[t], '*', latent[m[t]]),
      character(1)), collapse = ' + '), ')')
    for (t in these) added[[length(added) + 1L]] <- parsrow(names[t])
    current <- pars$param[cell]
    if (is.na(current) || !nzchar(current)) {
      expr <- paste0(as.numeric(pars$value[cell]), ' * ', scale)
    } else if (grepl('^[A-Za-z.][A-Za-z0-9._]*$', current) && !current %in% latent) {
      moved[[as.character(cell)]] <- parsrow(current, from = pars[cell, , drop = FALSE])
      expr <- paste0(current, ' * ', scale)
    } else {
      expr <- paste0('(', current, ') * ', scale)
    }
    pars$param[cell] <- expr
    pars$value[cell] <- NA
    pars$transform[cell] <- NA
    pars$indvarying[cell] <- FALSE
    for (e in effects) pars[[e]][cell] <- 'FALSE'
    for (l in levelflags) pars[[l]][cell] <- FALSE
  }
  model$pars <- .ctInsertMovedRows(pars, moved, added)
  model
}

# The fixed cells a model could free, as terms at zero: DRIFT cross effects,
# DIFFUSION and MANIFESTVAR covariances, LAMBDA cross loadings and
# TDPREDEFFECT cells fixed at zero. Each becomes a parameter with ctsem's
# default transform for its cell, kept only where that transform puts a raw
# zero at zero, so the extended model is the fit at zero. CINT and
# MANIFESTMEANS are left out: with the other free, either is identified only
# by the transient from the initial state.
#
# A covariance cell the 'z' construction solves as a zero covariance is not
# at raw zero in the fit: freed, it is the fit at its solved value (`at`,
# from `est`, the fit's raw estimate named by parameter), which with the
# other cells at the fit's values satisfies the remaining zeros, the one
# solution there is. Where that value differs by subject, or `est` is not
# given, the cell is not freed and is listed in `untested`.
.ctScoreFixedModel <- function(model, est = NULL) {
  pars <- model$pars
  latent <- model$latentNames; manifest <- model$manifestNames
  tdpred <- model$TDpredNames
  ct <- isTRUE(model$continuoustime)
  value <- suppressWarnings(as.numeric(pars$value))
  free <- is.na(pars$param) | !nzchar(pars$param)
  # A latent with no innovations of its own -- a cycle's clock, the level of a
  # higher order process -- is driven only as written: freeing a drift into it
  # would make a clock wander, and a covariance with nothing has no meaning.
  quiet <- pars$row[pars$matrix %in% 'DIFFUSION' & pars$row == pars$col & free &
    !is.na(value) & value == 0]
  # Likewise an indicator measured without error: a covariance between errors
  # with no variance is sd times sd times a correlation, zero along every
  # direction, and printed as an empty row (the book's FI-CLPM, MANIFESTVAR
  # fixed at 0).
  exact <- pars$row[pars$matrix %in% 'MANIFESTVAR' & pars$row == pars$col & free &
    !is.na(value) & value == 0]
  # Nor a drift into a latent a 'stationary' initial condition needs: it would
  # make the driving latent one too, and its dynamics, which the other
  # families extend, fixed (`.ctScoreStateModel`).
  closed <- which(.ctStationaryClosure(model))
  candidate <- which(free & !is.na(value) & value == 0 & (
    (pars$matrix %in% 'DRIFT' & pars$row != pars$col & !pars$row %in% c(quiet, closed)) |
      (pars$matrix %in% 'DIFFUSION' & pars$row > pars$col & !pars$row %in% quiet &
        !pars$col %in% quiet) |
      (pars$matrix %in% 'MANIFESTVAR' & pars$row > pars$col & !pars$row %in% exact &
        !pars$col %in% exact) |
      (pars$matrix %in% 'TDPREDEFFECT' & !pars$row %in% quiet) |
      pars$matrix %in% 'LAMBDA'))
  ns <- asNamespace('ctsem')
  dims <- function(m) switch(m, DRIFT = list(latent, latent),
    DIFFUSION = list(latent, latent), MANIFESTVAR = list(manifest, manifest),
    LAMBDA = list(manifest, latent), TDPREDEFFECT = list(latent, tdpred))
  solved <- lapply(stats::setNames(nm = c('DIFFUSION', 'MANIFESTVAR')), function(m)
    .ctScoreSolvedAt(model, est, m))
  terms <- list(); untested <- list()
  for (k in candidate) {
    m <- pars$matrix[k]; r <- pars$row[k]; c <- pars$col[k]
    tf <- .ctModelDefaultFreePar(m, r, c, ct)$transform
    atzero <- try(eval(parse(text = tf), list(param = 0), ns), silent = TRUE)
    if (inherits(atzero, 'try-error') || !isTRUE(abs(atzero) < 1e-12)) next
    d <- dims(m)
    cellname <- paste0(m, '[', d[[1]][r], ',', d[[2]][c], ']')
    at <- 0
    s <- solved[[m]]
    if (!is.null(s) && isTRUE(s$solved[r, c])) {
      # A covariance cell's default transform is the identity, so the solved
      # coordinate is the raw value.
      at <- s$value[r, c]
      if (!is.finite(at) || !isTRUE(abs(eval(parse(text = tf), list(param = at), ns) - at) <
          1e-12)) {
        untested[[length(untested) + 1L]] <- data.frame(cell = cellname,
          reason = s$reason, stringsAsFactors = FALSE)
        next
      }
    }
    name <- paste0('ctscorefixed', length(terms) + 1L)
    pars <- .ctScoreFreeCell(pars, k, name, tf)
    terms[[length(terms) + 1L]] <- data.frame(name = name, matrix = m, row = r,
      col = c, cell = cellname,
      suggested = paste0(tolower(m), '_', d[[1]][r], '_', d[[2]][c]),
      transform = tf, at = at, stringsAsFactors = FALSE)
  }
  model$pars <- pars
  list(model = model, terms = if (length(terms)) do.call(rbind, terms) else NULL,
    untested = if (length(untested)) do.call(rbind, untested) else NULL)
}

# For the covariance matrix `name` of `model`, the cells the 'z' construction
# solves as zero covariances (`.ctCovSolvedCells()`) and the coordinate each
# holds at the fit: the matrix's other cells below the diagonal at their
# fitted values, from the raw estimate `est` named by parameter, and the
# solved cells solved (`.ctCovZeroSolve()`), as the engine does. NA, with the
# `reason`, where one of those other cells varies by subject -- a random or
# predictor effect, or a cell that is no parameter of its own, such as an
# expression -- since each subject's solution then differs, and where `est`
# is not given. NULL where nothing is solved.
.ctScoreSolvedAt <- function(model, est, name) {
  solved <- .ctCovSolvedCells(model, name)
  if (!any(solved)) return(NULL)
  d <- nrow(solved)
  out <- list(solved = solved, value = matrix(NA_real_, d, d), reason = NULL)
  if (is.null(est) || is.null(names(est))) {
    out$reason <- "the fit's estimates were not given"
    return(out)
  }
  model <- .ctResolveLevels(model)
  pars <- model$pars
  varying <- .ctVaryingRows(model)
  effects <- grep('_effect$', names(pars), value = TRUE)
  A <- matrix(0, d, d)
  for (k in which(pars$matrix %in% name & pars$row > pars$col)) {
    r <- pars$row[k]; c <- pars$col[k]
    if (solved[r, c]) next
    if (is.na(pars$param[k]) || !nzchar(pars$param[k])) {
      A[r, c] <- suppressWarnings(as.numeric(pars$value[k]))
      next
    }
    # The row that holds the parameter, which for a PARS label is not this one.
    p <- pars$param[k]
    holder <- .ctParamRows(model, p)
    labels <- if (name == 'DIFFUSION') model$latentNames else model$manifestNames
    cell <- paste0(name, "[", labels[r], ",", labels[c], "]")
    if (!p %in% names(est) || is.na(holder)) {
      out$reason <- paste0("a zero covariance whose value in the fit depends on ",
        cell, ", which is not a parameter of its own")
      return(out)
    }
    if (varying[k] || varying[holder] ||
        any(.ctTipredEffectActive(unlist(pars[c(k, holder), effects])))) {
      out$reason <- paste0("a zero covariance whose value in the fit differs between ",
        "subjects, since ", cell, " does")
      return(out)
    }
    tf <- as.character(pars$transform[holder])
    if (is.na(tf) || !nzchar(tf)) tf <- 'param'
    A[r, c] <- tryCatch(as.numeric(eval(parse(text = tf), list(param = est[[p]]),
      asNamespace('ctsem'))), error = function(e) NA_real_)
  }
  A <- A + t(A)
  solution <- if (all(is.finite(A))) .ctCovZeroSolve(A, solved) else NULL
  if (is.null(solution)) {
    out$reason <- "a zero covariance whose value in the fit could not be solved"
    return(out)
  }
  out$value[solved] <- solution[solved]
  out
}

# A fixed cell written as a new parameter, shared by every extended model: no
# random effects, no predictor effects, at any level.
.ctScoreFreeCell <- function(pars, k, name, transform) {
  pars$param[k] <- name
  pars$value[k] <- NA
  pars$transform[k] <- transform
  pars$indvarying[k] <- FALSE
  for (e in grep('_effect$', names(pars), value = TRUE)) pars[[e]][k] <- 'FALSE'
  for (l in grep('^indvarying_', names(pars), value = TRUE)) pars[[l]][k] <- FALSE
  pars
}

# The fitted model with each MANIFESTMEANS cell it fixes written as a
# parameter at the fixed value (identity plus the value, so raw zero is the
# fit), for the test of a random intercept on that indicator about its stated
# mean. With LAMBDA = I such an intercept lies in what a CINT test with its
# covariances with the initial states spans; with more indicators than
# latents it is the indicator's own; where a latent has no free intercept at
# all it is the only test of its level.
.ctScoreFixedMeansModel <- function(model) {
  pars <- model$pars
  value <- suppressWarnings(as.numeric(pars$value))
  free <- is.na(pars$param) | !nzchar(pars$param)
  terms <- list()
  for (k in which(pars$matrix %in% 'MANIFESTMEANS' & free & !is.na(value))) {
    name <- paste0('ctscoremean', length(terms) + 1L)
    pars <- .ctScoreFreeCell(pars, k, name, if (value[k] == 0) 'param' else
      paste0('param + ', format(value[k], digits = 15)))
    terms[[length(terms) + 1L]] <- data.frame(name = name, row = as.integer(pars$row[k]),
      value = value[k], cell = paste0('MANIFESTMEANS[', model$manifestNames[pars$row[k]], ']'),
      stringsAsFactors = FALSE)
  }
  model$pars <- pars
  list(model = model, terms = if (length(terms)) do.call(rbind, terms) else NULL)
}

# Heterogeneity items for the MANIFESTMEANS cells the model fixes, from the
# model extended with each as a parameter at its value: per subject its score,
# Hessian row and second differences along it, and along it plus each effect
# in R for its covariances' projections. The effects' own pieces are the
# fit's, which the extended model reproduces at the fitted values; the scores
# leave out the prior, which the new coordinates would otherwise carry.
.ctScoreFixedMeans <- function(fit, est, pieces, R, step = 1e-3, cores = 2L) {
  base <- .ctFitBaseModel(fit)
  if (is.null(base)) return(NULL)
  fm <- .ctScoreFixedMeansModel(base)
  if (is.null(fm$terms)) return(NULL)
  ext <- .ctScoreExtendedSpec(fit, est, fm$model, fm$terms$name,
    'random intercepts on the fixed MANIFESTMEANS', cores)
  n <- nrow(pieces$S)
  th <- ext$theta
  own <- function(r) colSums(matrix(pieces$D2[, , match(r, pieces$second)], n))
  scores <- function(x) as.matrix(.ctJuliaGet(ext$module$ctsem_subject_gradients(
    ext$objective, .ctJuliaVector(as.numeric(x)), prior = FALSE))$scores)
  items <- NULL
  .ctBackendWithMaxChunks(cores, {
    S2 <- scores(ext$x2)
    if (nrow(S2) != n) stop("the extended model's units are not the fit's.", call. = FALSE)
    items <- lapply(seq_len(nrow(fm$terms)), function(a) {
      d <- ext$delta[a]
      e <- replace(numeric(length(ext$x2)), d, 1)
      plus <- scores(ext$x2 + step * e)
      minus <- scores(ext$x2 - step * e)
      H <- (plus - minus) / (2 * step)
      D2 <- colSums(plus - 2 * S2 + minus) / step^2
      Tmix <- matrix(vapply(R, function(r) {
        w <- e
        w[th[r]] <- 1
        h <- step * max(1, abs(est[r]))
        d2 <- colSums(scores(ext$x2 + h * w) - 2 * S2 + scores(ext$x2 - h * w)) / h^2
        (d2[th] - D2[th] - own(r)) / 2
      }, numeric(length(th))), length(th))
      list(name = fm$terms$cell[a], matrix = 'MANIFESTMEANS', row = fm$terms$row[a],
        col = 1L, slope = 1, s = S2[, d], H = H[, th, drop = FALSE], hjj = H[, d],
        D2 = D2[th], R = R, Tmix = if (length(R)) Tmix else NULL)
    })
  })
  items
}

# The fitted model rewritten as `model`, whose new parameters `names_new` are
# at zero where it is the fit: prepared as the fit was, evaluated with the
# fitted values, which must reproduce the fit, and the pieces a test reads
# from it -- its gradient there, where the fit's parameters sit in it
# (`theta`) and where the new ones do (`delta`), and its objective.
#
# At zero the extension has the fit's random effects, laid out as the fit's,
# so it takes the fit's modes: the known modes on its specification, which a
# rebuilt one lacks (without them a unit's search missed a maximum the fit's
# reached, and the extension read 0.95 nats below the fit), and the modes the
# fit's run ended with, held on both objectives for as long as `hold`, the
# caller's frame, runs (`.ctJuliaHoldModes()`), so that every difference
# reads one function, the one the estimate maximises.
.ctScoreExtendedSpec <- function(fit, est, model, names_new, label, cores = 2L,
  hold = parent.frame(), at = NULL) {
  spec <- fit$model_spec
  args <- .ctScoreNamedLevels(fit, .ctFitReplayArgs(fit))
  args$model <- model
  args$fit <- FALSE
  args$inits <- NULL
  args$datalong <- .ctFitLongData(fit)
  args$intoverpop <- .ctBackendIntOverPop(spec)
  args$verbose <- 0
  spec2 <- suppressMessages(suppressWarnings(do.call(ctFit, args)))
  if (!is.null(spec$laplace$known) && !is.null(spec2$laplace))
    spec2$laplace$known <- spec$laplace$known
  nm1 <- as.character(.ctFitRawParNames(fit))
  nm2 <- as.character(.ctBackendRawParameterNames(list(model_spec = spec2),
    .ctBackendNpar(spec2)))
  # Refused, with a class of its own, so that a caller testing several
  # families on one extended model can test them one at a time instead.
  refuse <- function(...) stop(structure(class = c('ctScoreRefused', 'error', 'condition'),
    list(message = paste0(...), call = NULL)))
  lost <- setdiff(nm1, nm2)
  if (length(lost) || anyDuplicated(nm2)) refuse("The extended model's parameters ",
    "do not line up with the fit's (", paste(lost, collapse = ', '), "), so ",
    label, " cannot be tested for this model.")
  if (!setequal(setdiff(nm2, nm1), names_new)) refuse("The extended model added ",
    "parameters other than its terms: ", paste(setdiff(setdiff(nm2, nm1), names_new),
      collapse = ', '))
  # The fit's correlation wall on the fit's parameters (`.ctCorBoundCarry()`),
  # set before any objective is built from the extension.
  spec2$corbound_blocks <- .ctCorBoundCarry(spec, nm1, nm2)
  x2 <- setNames(numeric(length(nm2)), nm2)
  x2[nm1] <- est
  # New coordinates at zero, or where `at` says: a stated mean freed sits at
  # its stated value.
  if (length(at)) x2[names(at)] <- at
  # Where the extension's Laplace term is not the fit's -- the fit's filter
  # integrates its intercept-type effects exactly given the others, and the
  # extension's terms read the state, so its filter cannot -- both are taken
  # one-step (`.ctLaplaceTerm()`): the fit's value and curvature as well as
  # the extension's, so that it reproduces the fit at zero and the projection
  # reads one objective. The fit's estimate is then not that objective's
  # maximum, which the projection allows for to first order, and the
  # one-step term's error in the new terms' slope stays in their scores: a
  # fast screen, its result marked (`onestep`). Where only one of the two has
  # effects its filter integrates, both are read on the layer, every effect a
  # Laplace coordinate (`layer`, `.ctLaplaceCoordinatesDiffer()`): on the
  # hybrid the one-step term is the layer's only where every unit has one
  # maximum, and a unit with two refused about one data set in twenty here.
  onestep <- .ctLaplaceTermsDiffer(spec, spec2)
  layer <- .ctLaplaceCoordinatesDiffer(spec, spec2)
  if (onestep) {
    fit <- .ctLaplaceOneStep(fit)
    spec2 <- .ctLaplaceOneStep(spec2)
  }
  modes <- .ctJuliaFitModes(fit)
  if (!is.null(modes) && length(modes$at) == length(nm1)) {
    .ctJuliaHoldModes(fit$model_spec, modes, layer = layer, envir = hold)
    at2 <- stats::setNames(numeric(length(nm2)), nm2)
    at2[nm1] <- modes$at
    modes$at <- as.numeric(at2)
    .ctJuliaHoldModes(spec2, modes, layer = layer, envir = hold)
  }
  value1 <- .ctJuliaEvaluate(fit, est, gradient = FALSE, cores = cores, layer = layer)$value
  ev <- .ctJuliaEvaluate(spec2, x2, gradient = TRUE, cores = cores, layer = layer)
  # With priors = TRUE every raw parameter has a normal(0, 1) prior, the new
  # terms too: at zero each adds log dnorm(0) to the objective, nothing to its
  # gradient and -1 to its own curvature. The terms are tested on the
  # likelihood, so that offset is expected and the curvature is taken back out
  # where it is read. Before this, a four-latent fit with priors refused at 40
  # terms x -0.919 = -36.76, which read as a model the rebuild had got wrong.
  priored <- isTRUE(fit$args$resolved$priors)
  offset <- if (priored) sum(stats::dnorm(x2[names_new], log = TRUE)) else 0
  if (!isTRUE(abs(ev$value - value1 - offset) <= 1e-6 * max(1, abs(value1)))) refuse(
    "The extended model does not reproduce the fit at zero (log likelihood ",
    signif(ev$value - offset, 8), " against ", signif(value1, 8), "), so its scores are ",
    "not the fit's and ", label, " cannot be tested for this model.")
  list(spec = spec2, x2 = x2, gradient = setNames(as.numeric(ev$gradient), nm2),
    theta = match(nm1, nm2), delta = match(names_new, nm2), priored = priored,
    module = .ctJuliaModule(spec$project), objective = .ctJuliaObjective(spec2, layer = layer),
    onestep = onestep, layer = layer)
}

# The fit's `popblocks` and `popunique` as the parameters they selected on the
# fitted model, level by level. They may name matrices, and an extended model
# moves a varying cell's parameter into PARS (the level terms move every
# varying MANIFESTMEANS and CINT, the drift terms CINT, the volatility terms
# the DIFFUSION diagonal), where a matrix name no longer reaches it: a block
# naming MANIFESTMEANS stopped the level family. A parameter keeps its name
# wherever it moves, so the extended model shares the fit's levels.
.ctScoreNamedLevels <- function(fit, args) {
  blocks <- .ctPopblocksRequested(args$popblocks)
  unique <- .ctPopuniqueRequested(args$popunique)
  if (!blocks && !unique) return(args)
  base <- .ctFitBaseModel(fit)
  if (is.null(base)) return(args)
  base <- .ctResolveLevels(base)
  levelnames <- c(base$subjectIDname, base$groupIDnames)
  columns <- stats::setNames(c('indvarying', if (length(base$groupIDnames))
    paste0('indvarying_', base$groupIDnames)), levelnames)
  if (blocks) args$popblocks <- lapply(.ctPopblocksResolve(base, args$popblocks,
    levelnames, columns), function(level) stats::setNames(lapply(level, function(b)
      list(rank = b$rank, effects = b$param)), vapply(level, `[[`, '', 'name')))
  if (unique) {
    asked <- .ctPopuniqueByLevel(args$popunique, levelnames)
    named <- Map(function(value, level) .ctPopuniqueSelect(base$pars, columns[[level]],
      value, level), asked, names(asked))
    named <- named[lengths(named) > 0L]
    args$popunique <- if (length(named)) named else NULL
  }
  args
}

# The fit's curvature on the extended model's objective, at the fitted values
# and the new terms at zero: its Hessian columns along the fit's own
# parameters, the information the projection reads where the extension is
# taken one-step (`.ctScoreExtendedSpec()`) and the fit's Hessian, of its own
# term, is not that objective's.
.ctScoreExtendedInformation <- function(ext, cores) {
  H <- .ctBackendWithMaxChunks(cores, as.matrix(.ctBackendJuliaValue(
    ext$module$ctsem_hessian_columns(ext$objective, .ctJuliaNumericVector(as.numeric(ext$x2)),
      .ctJuliaVector(as.integer(ext$theta))))))[ext$theta, , drop = FALSE]
  -(H + t(H)) / 2
}

# The fitted model extended at zero by every term of the families asked for --
# drift terms written into CINT, fixed cells freed -- prepared as the fit was,
# evaluated with the fitted values and the terms at zero, which must reproduce
# the fit; then its gradient, its Hessian columns along the terms, and for the
# robust variance each unit's gradient, whose term columns projected on the
# fitted parameters are the contributions. One model for both families, so
# one compile of the extended shape.
.ctScoreExtended <- function(fit, est, families, information = NULL, order = 2L,
  variance = 'robust', nboot = 1999L, cores = 2L, exact = Inf, quadrature = FALSE) {
  spec <- fit$model_spec
  base <- .ctFitBaseModel(fit)
  if (is.null(base)) stop("This fit does not carry the model it was built from, ",
    "so the extended model cannot be prepared. Refit with this version of ctsem.",
    call. = FALSE)
  model <- base
  state <- vol <- fixed <- lev <- NULL
  if ('level' %in% families) {
    if (length(families) > 1L) stop("The level terms are tested on a model of ",
      "their own.", call. = FALSE)
    lev <- .ctScoreLevelModel(model)
    model <- lev$model
  }
  if ('statedependence' %in% families) {
    state <- .ctScoreStateModel(model, order)
    model <- state$model
  }
  if ('volatility' %in% families) {
    vol <- .ctScoreVolatilityModel(model)
    model <- vol$model
  }
  if ('fixed' %in% families) {
    fixed <- .ctScoreFixedModel(model, stats::setNames(est, .ctFitRawParNames(fit)))
    model <- fixed$model
  }
  names_new <- c(state$terms$name, vol$terms$name, fixed$terms$name, lev$terms$name)
  if (!length(names_new)) return(list(statedependence = NULL, volatility = NULL,
    fixed = if (!is.null(fixed$untested)) list(untested = fixed$untested), level = NULL))
  ext <- .ctScoreExtendedSpec(fit, est, model, names_new, paste(families, collapse = ', '),
    cores, at = if (!is.null(fixed$terms)) stats::setNames(fixed$terms$at, fixed$terms$name))
  x2 <- ext$x2
  grad <- ext$gradient
  theta <- ext$theta
  delta <- ext$delta
  priored <- ext$priored
  module <- ext$module
  objective <- ext$objective
  # Taken one-step, or on the layer, the fit's own Hessian is of another
  # objective.
  if (isTRUE(ext$onestep) || isTRUE(ext$layer))
    information <- .ctScoreExtendedInformation(ext, cores)
  robust <- variance == 'robust'
  # By quadrature (`intoverpop = 'quadrature'`), every piece is the placed
  # objective's, every unit on its rule (`.ctScorePlaced()`): the gradient, the
  # units' scores, and the curvature, the fit's block of it the information.
  placed <- if (isTRUE(quadrature) && .ctBackendIntOverPop(fit$model_spec) == 'laplace')
    .ctScorePlaced(ext, cores) else NULL
  if (!is.null(placed)) {
    grad <- placed$gradient
    information <- -(placed$hessian[theta, theta, drop = FALSE] +
      t(placed$hessian[theta, theta, drop = FALSE])) / 2
  }
  units <- if (!robust) NULL else if (!is.null(placed)) placed$units else
    as.matrix(.ctBackendWithMaxChunks(cores, .ctJuliaGet(
      module$ctsem_subject_gradients(objective, .ctJuliaVector(as.numeric(x2)))))$scores)
  # Which terms are projected exactly. The projection needs the Hessian's
  # column along each term, two evaluations of the extended model apiece, and
  # those are what these families cost: 84 state-dependence terms on a
  # six-latent model took 17 minutes, where every term's score, one per unit,
  # is a single gradient pass. So with the robust variance each family's terms
  # are first ranked on their unprojected scores and only the strongest
  # `exact` projected. An unprojected statistic is the same sum (the fitted
  # parameters' scores sum to zero at the estimate) over a variance that the
  # projection, under the model, can only reduce: conservative, and it ranks
  # the terms the projection then sharpens. The model variance needs every
  # column.
  family <- rep(c('statedependence', 'volatility', 'fixed', 'level'),
    c(length(state$terms$name), length(vol$terms$name), length(fixed$terms$name),
      length(lev$terms$name)))
  chosen <- seq_along(delta)
  if (robust && length(delta) && is.finite(exact)) {
    raw <- units[, delta, drop = FALSE]
    spread <- colSums(sweep(raw, 2, colMeans(raw))^2) * nrow(raw) / (nrow(raw) - 1)
    lm0 <- ifelse(spread > 0, grad[delta]^2 / spread, 0)
    chosen <- sort(unlist(lapply(split(seq_along(delta), family), function(i)
      i[order(-lm0[i])][seq_len(min(exact, length(i)))]), use.names = FALSE))
  }
  if (is.null(information)) information <- -.ctBackendHessian(fit, est)
  Iinv <- .ctScorePinv(information)$inv
  A <- matrix(0, length(delta), length(theta))
  vm <- rep(NA_real_, length(delta))
  Vmodel <- NULL
  if (length(chosen)) {
    Hc <- if (!is.null(placed)) placed$hessian[, delta[chosen], drop = FALSE] else
      .ctBackendWithMaxChunks(cores, as.matrix(.ctBackendJuliaValue(
        module$ctsem_hessian_columns(objective, .ctJuliaNumericVector(as.numeric(x2)),
          .ctJuliaVector(as.integer(delta[chosen]))))))
    if (priored) Hc[cbind(delta[chosen], seq_along(chosen))] <-
      Hc[cbind(delta[chosen], seq_along(chosen))] + 1
    Idt <- -t(Hc[theta, , drop = FALSE])
    Idd <- -Hc[delta[chosen], , drop = FALSE]
    Idd <- (Idd + t(Idd)) / 2
    A[chosen, ] <- Idt %*% Iinv
    Vmodel <- Idd - A[chosen, , drop = FALSE] %*% t(Idt)
    Vmodel <- (Vmodel + t(Vmodel)) / 2
    vm[chosen] <- diag(Vmodel)
  }
  Ustar <- grad[delta] - as.numeric(A %*% grad[theta])
  V <- Vmodel
  Cstar <- NULL
  if (robust) {
    Cstar <- units[, delta, drop = FALSE] - units[, theta, drop = FALSE] %*% t(A)
    V <- .ctScoreSpread(Cstar)
  }
  v <- diag(V)
  per <- data.frame(LM = ifelse(v > 0, Ustar^2 / v, NA_real_),
    # One term at a time, the fitted parameters free: the expected parameter
    # change. NA where the fitted model's curvature cannot size it -- a strong
    # departure can leave the projected information indefinite -- and for a
    # term not projected.
    estimate = ifelse(!is.na(vm) & vm > 0, Ustar / vm, NA_real_),
    se = ifelse(!is.na(vm) & vm > 0, sqrt(pmax(v, 0)) / vm, NA_real_),
    projected = seq_along(delta) %in% chosen)
  per$p <- stats::pchisq(per$LM, 1, lower.tail = FALSE)
  jointtest <- function(idx, label) {
    j <- .ctScoreJoint(Ustar[idx], V[idx, idx, drop = FALSE],
      if (!is.null(Cstar)) Cstar[, idx, drop = FALSE], nboot)
    data.frame(set = label, LM = j$LM, df = j$df, p = j$p, stringsAsFactors = FALSE)
  }
  out <- list(statedependence = NULL, volatility = NULL, fixed = NULL, level = NULL)
  latent <- base$latentNames
  if (!is.null(state$terms)) {
    idx <- seq_len(nrow(state$terms))
    terms <- state$terms
    table <- data.frame(term = terms$term, equation = terms$equation,
      reads = paste0('DRIFT[', latent[terms$i], ',', latent[terms$k], '] varies with ',
        terms$by, ifelse(is.na(terms$other), '', paste0(', or DRIFT[',
          latent[terms$i], ',', latent[pmax(1L, terms$other)], '] with ',
          latent[terms$k]))),
      per[idx, c('LM', 'p', 'estimate', 'se', 'projected')], stringsAsFactors = FALSE)
    table$p_adj <- .ctScoreBH(table$p)
    if (!is.null(Cstar)) table$class <- .ctScoreClasses(Cstar[, idx, drop = FALSE])
    table$i <- terms$i; table$k <- terms$k
    table$monomial <- vapply(seq_len(nrow(terms)), function(r)
      paste(match(strsplit(terms$term[r], '*', fixed = TRUE)[[1]], latent), collapse = ','),
      character(1))
    joint <- do.call(rbind, c(lapply(unique(terms$equation), function(e)
      jointtest(idx[terms$equation == e], e)), list(jointtest(idx, 'all'))))
    names(joint)[1] <- 'equation'
    out$statedependence <- list(terms = table[order(table$p, -table$LM), , drop = FALSE],
      joint = joint)
  }
  if (!is.null(vol$terms)) {
    idx <- length(state$terms$name) + seq_len(nrow(vol$terms))
    terms <- vol$terms
    table <- data.frame(term = paste0('sd(', terms$equation, ') x exp(v * ', terms$by, ')'),
      equation = terms$equation, by = terms$by,
      per[idx, c('LM', 'p', 'estimate', 'se', 'projected')], stringsAsFactors = FALSE)
    table$p_adj <- .ctScoreBH(table$p)
    if (!is.null(Cstar)) table$class <- .ctScoreClasses(Cstar[, idx, drop = FALSE])
    table$i <- terms$i; table$m <- terms$m
    joint <- do.call(rbind, c(lapply(unique(terms$equation), function(e)
      jointtest(idx[terms$equation == e], e)), list(jointtest(idx, 'all'))))
    names(joint)[1] <- 'equation'
    out$volatility <- list(terms = table[order(table$p, -table$LM), , drop = FALSE],
      joint = joint)
  }
  if (!is.null(fixed$terms)) {
    idx <- length(state$terms$name) + length(vol$terms$name) + seq_len(nrow(fixed$terms))
    terms <- fixed$terms
    table <- data.frame(cell = terms$cell, matrix = terms$matrix, row = terms$row,
      col = terms$col, per[idx, c('LM', 'p', 'estimate', 'se', 'projected')],
      suggested = terms$suggested, transform = terms$transform, at = terms$at,
      stringsAsFactors = FALSE)
    table$p_adj <- .ctScoreBH(table$p)
    if (!is.null(Cstar)) table$class <- .ctScoreClasses(Cstar[, idx, drop = FALSE])
    out$fixed <- list(cells = table[order(table$p, -table$LM), , drop = FALSE],
      joint = jointtest(idx, 'all'))
  }
  if (!is.null(fixed$untested)) out$fixed$untested <- fixed$untested
  if (!is.null(lev$terms)) {
    idx <- length(state$terms$name) + length(vol$terms$name) + length(fixed$terms$name) +
      seq_len(nrow(lev$terms))
    terms <- lev$terms
    table <- data.frame(term = paste0(terms$parameter, ' ~ ', terms$level),
      parameter = terms$parameter, level = terms$level, matrix = terms$matrix,
      per[idx, c('LM', 'p', 'estimate', 'se', 'projected')], stringsAsFactors = FALSE)
    table$p_adj <- .ctScoreBH(table$p)
    if (!is.null(Cstar)) table$class <- .ctScoreClasses(Cstar[, idx, drop = FALSE])
    joint <- do.call(rbind, c(lapply(unique(terms$level), function(l)
      jointtest(idx[terms$level == l], l)), list(jointtest(idx, 'all'))))
    names(joint)[1] <- 'level'
    out$level <- list(terms = table[order(table$p, -table$LM), , drop = FALSE],
      joint = joint)
  }
  out$onestep <- isTRUE(ext$onestep) && is.null(placed)
  out$quadrature <- !is.null(placed)
  out
}

# The extended model's pieces with every unit's random effects integrated by
# a quadrature rule placed at its mode and curvature (the refined fit's placed
# objective, `ctsem_laplace_placed`, every unit forced onto a rule): the
# gradient at the fitted values with the new terms at zero, each unit's score
# there, and the Hessian. At zero the extension is the fitted model, so these
# carry no Laplace error; the fit's estimate is the Laplace one, which the
# projection allows for to first order.
.ctScorePlaced <- function(ext, cores) {
  module <- ext$module
  x <- .ctJuliaNumericVector(as.numeric(ext$x2))
  nunits <- .ctJuliaUnits(ext$spec)
  .ctBackendWithMaxChunks(cores, {
    cont <- do.call(module$ctsem_laplace_continuation, c(list(ext$objective, x),
      .ctLaplaceContinuationArgs(.ctLaplaceContinueDefaults, force = seq_len(nunits))))
    placed <- module$ctsem_laplace_placed(cont)
    ev <- .ctJuliaGet(module$ctsem_laplace_placed_evaluate(placed, x, gradient = TRUE))
    sg <- .ctJuliaGet(module$ctsem_subject_gradients(placed, x))
    npar <- length(ext$x2)
    H <- matrix(as.numeric(.ctBackendJuliaValue(module$ctsem_laplace_placed_hessian(placed,
      x))), npar, npar)
    list(value = as.numeric(ev$value), gradient = setNames(as.numeric(ev$gradient),
      names(ext$x2)), units = as.matrix(sg$scores), hessian = (H + t(H)) / 2)
  })
}

# The fitted model with every parameter of the dynamics and of the error
# variances that does not vary over subjects (free DRIFT, DIFFUSION and
# MANIFESTVAR parameters, and PARS parameters their expressions read) shifted,
# on its own scale, by a term times each intercept that does (MANIFESTMEANS,
# CINT): the parameter becomes `<parameter> + t1 * <intercept 1> + ...`
# wherever a cell reads it, each t a new PARS parameter at zero, identity
# transform. A cell that is the parameter is rewritten so; one whose
# expression reads it, such as a second process's drift `a_eta - d_med` held
# at a gap from the first's, has the shifted parameter put in its place,
# since a term that moved the parameter in one cell and not the other would
# test a model that no longer ties them. The parameter and the intercepts are
# declared in PARS so the cells can read them, each straight after the cell it
# came from (.ctInsertMovedRows), so the random effects keep their order. At
# zero the model is the fit. The constant part of an intercept is the
# parameter's own direction, which the test projects out, so a subject's
# value needs no centring.
.ctScoreLevelModel <- function(model) {
  pars <- model$pars
  varies <- as.logical(.ctResolveLevels(model)$pars$indvarying) %in% TRUE
  isp <- pars$matrix %in% 'PARS'
  plain <- !is.na(pars$param) & grepl('^[A-Za-z.][A-Za-z0-9._]*$', pars$param) &
    !pars$param %in% c('stationary', model$latentNames, model$TDpredNames)
  anyvaries <- function(p) any(varies[plain & pars$param == p])
  levels <- unique(pars$param[pars$matrix %in% c('MANIFESTMEANS', 'CINT') & plain & varies])
  dynamic <- pars$matrix %in% c('DRIFT', 'DIFFUSION', 'MANIFESTVAR') & !is.na(pars$param)
  declared <- unique(pars$param[isp & plain])
  tested <- unique(c(pars$param[dynamic & plain], declared[vapply(declared, function(p)
    any(.ctScoreReads(pars$param[dynamic & !plain], p)), logical(1))]))
  tested <- tested[!vapply(tested, anyvaries, logical(1))]
  if (!length(levels) || !length(tested)) return(list(model = model, terms = NULL))
  terms <- expand.grid(level = levels, parameter = tested, stringsAsFactors = FALSE)
  terms <- terms[, c('parameter', 'level')]
  terms$name <- paste0('ctscorelevel', seq_len(nrow(terms)))
  terms$matrix <- vapply(terms$parameter, function(p)
    pars$matrix[dynamic & (pars$param == p | .ctScoreReads(pars$param, p))][1], '')
  effects <- grep('_effect$', names(pars), value = TRUE)
  levelflags <- grep('^indvarying_', names(pars), value = TRUE)
  nextpars <- max(c(0L, as.integer(pars$row[isp])))
  moved <- list(); added <- list()
  for (p in c(tested, levels)) {
    if (any(isp & pars$param %in% p)) next
    cell <- which(!isp & plain & pars$param == p)[1]
    r <- pars[cell, , drop = FALSE]
    r$matrix <- 'PARS'; r$row <- nextpars <- nextpars + 1L; r$col <- 1L
    moved[[as.character(cell)]] <- r
  }
  template <- pars[which(!isp)[1], , drop = FALSE]
  for (k in seq_len(nrow(terms))) {
    r <- template
    r$matrix <- 'PARS'; r$row <- nextpars <- nextpars + 1L; r$col <- 1L
    r$param <- terms$name[k]; r$value <- NA; r$transform <- 'param'; r$indvarying <- FALSE
    if ('sdscale' %in% names(r)) r$sdscale <- 1
    for (e in effects) r[[e]] <- 'FALSE'
    for (l in levelflags) r[[l]] <- FALSE
    added[[k]] <- r
  }
  for (p in tested) {
    cells <- which(!isp & plain & pars$param == p)
    tt <- terms[terms$parameter == p, , drop = FALSE]
    shifted <- paste0(p, ' + ', paste0(tt$name, ' * ', tt$level, collapse = ' + '))
    pars$param[cells] <- shifted
    pars$value[cells] <- NA
    pars$transform[cells] <- NA
    pars$indvarying[cells] <- FALSE
    for (e in effects) pars[[e]][cells] <- 'FALSE'
    for (l in levelflags) pars[[l]][cells] <- FALSE
    # Expressions that read it, in any matrix, read the shifted parameter.
    reads <- which(!plain & !is.na(pars$param) & .ctScoreReads(pars$param, p))
    pars$param[reads] <- .ctScoreReplaceName(pars$param[reads], p,
      paste0('(', shifted, ')'))
  }
  model$pars <- .ctInsertMovedRows(pars, moved, added)
  rownames(terms) <- NULL
  list(model = model, terms = terms)
}

# Whether each expression in `x` reads the name `p`: as a whole name, not part
# of a longer one (`a_eta` is not read by `a_eta2` or `ba_eta`), and not as a
# function called by that name.
.ctScoreNamePattern <- function(p) paste0('(?<![A-Za-z0-9._])',
  gsub('.', '\\.', p, fixed = TRUE), '(?![A-Za-z0-9._]|\\s*\\()')
.ctScoreReads <- function(x, p) !is.na(x) & grepl(.ctScoreNamePattern(p), x, perl = TRUE)
.ctScoreReplaceName <- function(x, p, by) gsub(.ctScoreNamePattern(p), by, x, perl = TRUE)

#' Plot score tests
#'
#' Each family's strongest relaxations: their p-values adjusted within the
#' family (Benjamini-Hochberg), on a log scale, against .05.
#'
#' @param x A \code{ctScoreTest}.
#' @param family Families to show; by default every family tested.
#' @param top Relaxations shown per family, the strongest.
#' @param ... Unused.
#' @return A ggplot.
#' @examples
#' \donttest{
#' if (isTRUE(ctJuliaStatus()$available)) {
#' # A model in which the innovations of the two processes are uncorrelated.
#' model <- ctModel(LAMBDA = diag(2),
#'   manifestNames = c('Y1', 'Y2'), latentNames = c('eta1', 'eta2'),
#'   TDpredNames = 'TD1',
#'   DIFFUSION = matrix(c('diff_eta1', 0, 0, 'diff_eta2'), 2, 2))
#' fit <- ctFit(ctsemtestdat, model, priors = TRUE)
#' st <- ctScoreTest(fit)
#'
#' # A panel per family, with a point per relaxation at its adjusted
#' # p-value. Points left of the dashed line at .05, drawn in red, are the
#' # relaxations the data support, here freeing DIFFUSION[eta2,eta1].
#' print(plot(st))
#'
#' # Two of the families, with the three strongest relaxations of each.
#' print(plot(st, family = c('fixed', 'heterogeneity'), top = 3))
#' }
#' }
#' @export
plot.ctScoreTest <- function(x, family = NULL, top = 10, ...) {
  rows <- list()
  add <- function(fam, label, p, p_adj) if (length(label))
    rows[[length(rows) + 1L]] <<- data.frame(family = fam, label = label, p = p,
      p_adj = p_adj, stringsAsFactors = FALSE)
  h <- x$heterogeneity
  if (!is.null(h)) add('heterogeneity', h$parameter, h$p, h$p_adj)
  b <- rbind(if (NROW(x$blocks$dimension)) data.frame(label = paste0('dimension in ',
    x$blocks$dimension$block), p = x$blocks$dimension$p),
    if (NROW(x$blocks$covariance)) data.frame(label = paste0('covariance ',
      x$blocks$covariance$blocks), p = x$blocks$covariance$p))
  if (NROW(b)) add('blocks', b$label, b$p, .ctScoreBH(b$p))
  for (cv in x$covariates) add('covariates', paste0(cv$joint$covariate, ': ',
    cv$parameters$parameter), cv$parameters$p, cv$parameters$p_adj)
  s <- x$statedependence$terms
  if (!is.null(s)) add('statedependence', paste0(s$equation, ': ', s$term), s$p, s$p_adj)
  v <- x$volatility$terms
  if (!is.null(v)) add('volatility', v$term, v$p, v$p_adj)
  f <- x$fixed$cells
  if (!is.null(f)) add('fixed', f$cell, f$p, f$p_adj)
  l <- x$level$terms
  if (!is.null(l)) add('level', l$term, l$p, l$p_adj)
  i <- x$instability$parameters
  if (!is.null(i)) add('instability', paste0(i$parameter, ' (', i$ordering, ')'), i$p, i$p_adj)
  d <- do.call(rbind, rows)
  if (!is.null(family)) d <- d[d$family %in% family, , drop = FALSE]
  if (is.null(d) || !nrow(d)) stop("Nothing to plot: no relaxation was tested",
    if (!is.null(family)) " in those families", ".", call. = FALSE)
  d <- do.call(rbind, lapply(split(d, factor(d$family, levels = unique(d$family))),
    function(part) utils::head(part[order(part$p_adj, part$p), , drop = FALSE], top)))
  # One axis per family, strongest at the top.
  d$key <- factor(paste0(d$family, '|', d$label),
    levels = rev(unique(paste0(d$family, '|', d$label))))
  ggplot(d, aes(x = pmax(.data$p_adj, 1e-12), y = .data$key)) +
    geom_vline(xintercept = 0.05, linetype = 2, colour = 'grey50') +
    geom_point(aes(colour = .data$p_adj < 0.05), size = 2.5) +
    scale_colour_manual(values = c(`FALSE` = 'grey35', `TRUE` = 'red3'), guide = 'none') +
    scale_x_log10() + scale_y_discrete(labels = function(k) sub('^[^|]*[|]', '', k)) +
    facet_grid(family ~ ., scales = 'free_y', space = 'free_y') +
    labs(x = 'p adjusted within the family (dashed: .05)', y = NULL,
      title = paste0('Score tests (', x$info$route, ' route, ', x$info$nsubjects,
        ' subjects)')) +
    theme_bw() + theme(strip.text.y = element_text(angle = 0),
      plot.title.position = 'plot')
}

#' @export
print.ctScoreTest <- function(x, digits = 3, top = 10, ...) {
  # A bootstrap p-value cannot go below 1 / (nboot + 1); at that floor it
  # says only that, and the table is ordered by the statistic instead.
  floor <- if (identical(x$info$variance, 'robust') && !is.null(x$info$nboot))
    1 / (x$info$nboot + 1) else 0
  fmt <- function(d) {
    d <- utils::head(d, top)
    atfloor <- if ('p' %in% names(d)) !is.na(d$p) & d$p <= floor * (1 + 1e-9) else NULL
    num <- vapply(d, is.numeric, logical(1))
    d[num] <- lapply(d[num], function(v) signif(v, digits))
    if (any(atfloor)) {
      d$p <- format(d$p)
      d$p[atfloor] <- paste0('<', signif(floor, 2))
    }
    print(d, row.names = FALSE)
  }
  cat("Score tests at the fitted model (", x$info$route, " route, ",
    x$info$nsubjects, " subjects)\n", sep = "")
  for (e in names(x$info$errors)) cat("Not tested, ", e, ": ", x$info$errors[[e]], "\n", sep = "")
  if (length(x$info$quadrature)) cat("By quadrature, every subject's random effects ",
    "integrated: ", paste(x$info$quadrature, collapse = ', '), "\n", sep = "")
  if (length(x$info$onestep)) cat("On one-step Laplace, a fast screen: ",
    paste(x$info$onestep, collapse = ', '), " (the extended model's filter cannot ",
    "integrate the fit's intercept effects); check a small p with ",
    "intoverpop = 'quadrature'\n", sep = "")
  withclass <- function(d, cols) if (any(!is.na(d$class))) c(cols, 'class') else cols
  if (!is.null(x$heterogeneity)) {
    h <- x$heterogeneity
    cat("\nUnmodelled heterogeneity (LM: the variance with its covariances, against an\n",
      "equal mixture of chi-square on df - 1 and df; z: the variance alone, one sided;\n",
      "p: the stronger of the two; sd: the one-step estimate, raw scale):\n", sep = "")
    fmt(h[, withclass(h, c('parameter', 'LM', 'df', 'z', 'p', 'p_adj', 'sd', 'sd_natural'))])
    w <- unique(unlist(strsplit(stats::na.omit(h$with), ', ', fixed = TRUE)))
    if (length(w)) cat("Covariances with: ", paste(w, collapse = ', '), "\n", sep = "")
  }
  if (!is.null(x$dimension)) {
    d <- x$dimension
    cat("\nA new dimension of variation: largest eigenvalue ", signif(d$lambda, digits),
      ", p = ", signif(d$p, digits), "; direction ",
      paste0(names(d$direction)[seq_len(min(5, length(d$direction)))], ' ',
        signif(d$direction[seq_len(min(5, length(d$direction)))], 2), collapse = ', '),
      "\n", sep = "")
  }
  if (!is.null(x$blocks$dimension) && nrow(x$blocks$dimension)) {
    cat("\nA further shared dimension within each block (largest eigenvalue, its\n",
      "loadings projected out):\n", sep = "")
    fmt(x$blocks$dimension)
  }
  if (!is.null(x$blocks$covariance) && nrow(x$blocks$covariance)) {
    cat("\nCovariances between blocks' dimensions (joint, two sided):\n")
    fmt(x$blocks$covariance)
  }
  for (c in x$covariates) {
    j <- c$joint
    cat("\nCovariate ", j$covariate, " (", if (isTRUE(j$inmodel))
      "in the model, on the parameters it does not act on" else j$form,
      if (!is.null(j$cut)) paste0(", above ", signif(j$cut, digits)), "): joint LM = ",
      signif(j$LM, digits), ", df = ", j$df, ", p = ", signif(j$p, digits), "\n", sep = "")
    fmt(c$parameters[, withclass(c$parameters, c('parameter', 'estimate', 'se',
      'estimate_natural', 'p', 'p_adj'))])
  }
  if (!is.null(x$statedependence)) {
    cat("\nState dependence of the drift:\n")
    fmt(x$statedependence$joint)
    cat("\n")
    fmt(x$statedependence$terms[, withclass(x$statedependence$terms, c('term',
      'equation', 'LM', 'p', 'p_adj', 'estimate', 'se', 'reads'))])
  }
  if (!is.null(x$volatility)) {
    cat("\nState dependence of the innovation sds:\n")
    fmt(x$volatility$joint)
    cat("\n")
    fmt(x$volatility$terms[, withclass(x$volatility$terms, c('term', 'LM', 'p', 'p_adj',
      'estimate', 'se'))])
  }
  if (!is.null(x$instability)) {
    cat("\nStability across the data (Nyblom), by ordering:\n")
    fmt(x$instability$joint)
    cat("\n")
    fmt(x$instability$parameters[, c('ordering', 'parameter', 'L', 'p', 'p_adj')])
    if (!is.null(x$instability$note)) cat(x$instability$note, "\n")
  }
  if (!is.null(x$fixed$cells)) {
    j <- x$fixed$joint
    cat("\nFixed cells, freed: joint LM = ", signif(j$LM, digits), ", df = ", j$df,
      ", p = ", signif(j$p, digits), "\n", sep = "")
    fmt(x$fixed$cells[, withclass(x$fixed$cells, c('cell', 'LM', 'p', 'p_adj',
      'estimate', 'se'))])
  }
  if (!is.null(x$fixed$untested)) {
    u <- x$fixed$untested
    cat("\nFixed cells not tested:\n")
    for (r in seq_len(nrow(u))) cat("  ", u$cell[r], ": ", u$reason[r], "\n", sep = "")
  }
  if (!is.null(x$level)) {
    cat("\nParameters depending on a subject's level (a random intercept):\n")
    fmt(x$level$joint)
    cat("\n")
    fmt(x$level$terms[, withclass(x$level$terms, c('term', 'LM', 'p', 'p_adj',
      'estimate', 'se'))])
  }
  invisible(x)
}
