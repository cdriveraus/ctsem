# The model a fit runs, from the model as written: everything ctFit() does to
# the model that does not read the data, in the order it does it. One function
# so that whatever prepares a model -- ctFit(), the julia generator
# (`.ctGeneratePrepare()`), and the julia preparation of a model no fit has
# set up (`.ctJuliaCanonicalModel()`, `.ctJuliaAugmentRandomEffects()`) --
# runs the same code. Each of the last three carried a partial copy of these
# steps, and a step added to ctFit reached none of them.
#
# `nsubjects` is all it needs of the data: one subject can carry no individual
# variation. `checks = FALSE` leaves out what is said about the fit rather than
# the model -- the estimation messages, and the refusals and warnings about the
# estimator -- for a caller that fits nothing. Returns the prepared model and
# what ctFit() reads of the route and rank it resolved.
#' @keywords internal
.ctFitModelSetup <- function(model, backend = 'julia', intoverpop = 'auto',
  intoverstates = TRUE, optimize = TRUE, priors = FALSE,
  priorscope = if (isTRUE(priors)) 'all' else 'none', stationary = FALSE,
  binomial = FALSE, poprank = 'auto', poprankexplicit = FALSE, popblocks = NULL,
  popunique = NULL, nlcontrol = .ctNlcontrolDefaults(list(), backend),
  saveComplexPars = FALSE, nsubjects = 2L, checks = TRUE, fit = TRUE,
  estonly = FALSE) {
  ctstanmodel <- model

  # Stated population means (ctPar(mean=); R/ctFixedPopMeans.R), converted to
  # raw values from the cells' own transforms before anything here rewrites a
  # cell: the augmented route moves a varying intercept into a carrier state.
  # Done in ctFit() before it called this, the conversion never reached the
  # generator, which drew persons about ten times a MANIFESTMEANS mean.
  if(length(ctstanmodel[['fixedpopmeans']])){
    if(!backend %in% 'julia') stop("Stated population means, ctPar(mean=), ",
      "need backend='julia'.", call.=FALSE)
    ctstanmodel[['fixedpopmeansraw']] <- .ctFixedPopMeans(ctstanmodel)
  }

  # The working copy carries a logical column per level; the model the fit
  # keeps as its base (`modelbase`, `ctstanmodelbase`) stays as the user wrote
  # it. See .ctResolveLevels().
  ctm <- .ctResolveLevels(ctstanmodel)

  # A model saved before the correlation squash moved into the covariance
  # construction still carries it in its own parameter table, so it would be
  # applied twice. Both backends come through here.
  .ctCheckLegacyCovTransformModel(ctm$pars)
  # Covariances given as ctCov(), rewritten for the construction in force now.
  ctm <- .ctCovRefresh(ctm)

  if(!is.null(ctm$TIpredAuto) && ctm$TIpredAuto %in% c(1L,TRUE)){ #if auto tipred, set all effects to true
    for(tip in ctm$TIpredNames){
      ctm$pars[[paste0(tip,'_effect')]] <- 'TRUE'
    }
  }

  # What is said about the fit, not the model: not for a caller that fits
  # nothing.
  if(checks){
    # The prior has something to act on only where a level has two or more
    # varying parameters; naming it on a model with none read as a prior that
    # was not there.
    hascorrelations <- any(vapply(.ctVaryingColumns(ctm), function(cl)
      length(.ctVaryingParams(ctm, cl)) > 1L, logical(1)))
    if(optimize && priorscope %in% 'randomCorr') message(
      if(hascorrelations) paste0("Maximum likelihood estimation requested, ",
        "with priors on the random-effect correlations. priors=FALSE for none, ",
        "TRUE for all.") else "Maximum likelihood estimation requested")
    if(optimize && !priors) message("Maximum likelihood estimation requested")
    # `optimcontrol$is` is refused above, so it is NULL by the time we get here
    # and the importance-sampling wording this used to choose is unreachable.
    # Importance sampling is now `optimcontrol$uncertainty='is'`, which runs
    # after optimization rather than instead of it, so the estimation this
    # message describes is a posteriori either way.
    if(optimize && priorscope %in% 'all') message(
      "Maximum a posteriori estimation requested")
    # Naming stan here was wrong for half the fits it described: with
    # backend='julia' the engine runs its own sampler -- NUTS through a learnt
    # map on the joint posterior by default, NUTS or tempered SMC on a marginal
    # one -- and a user reading
    # "Stan's NUTS sampler" on a julia fit has no reason to believe the julia
    # sampler ran at all. Which kernel is said when sampling starts.
    if(!optimize) message("Bayesian estimation via ",
      if(identical(backend,'julia')) "the julia engine's sampler" else
        "Stan's NUTS sampler", " requested")
  }


  ###stationarity
  # The latent processes start from the distribution DRIFT, CINT and DIFFUSION
  # settle into, so T0MEANS and T0VAR are no longer parameters: their cells
  # are fixed here and the engine writes the stationary moments over them
  # (`_ctsem_stationary!`). The fixed values are placeholders nothing reads.
  if(!isTRUE(stationary) && !isFALSE(stationary)) stop(
    'stationary must be TRUE or FALSE.', call. = FALSE)
  # Single elements: T0MEANS or T0VAR cells named 'stationary', which the
  # engine starts at their stationary values (`.ctStationaryElements`). Under
  # stationary = TRUE every T0 cell is written over below, so there are none.
  ctm$t0stationary <- NULL
  if(any(ctm$pars$param %in% 'stationary') && !stationary) {
    if(!identical(backend, 'julia')) stop(
      "A 'stationary' T0MEANS or T0VAR cell needs backend = 'julia'.", call. = FALSE)
    ctm <- .ctStationaryElements(ctm)
  }
  if(stationary) {
    if(!identical(backend, 'julia')) stop(
      "stationary = TRUE needs backend = 'julia'.", call. = FALSE)
    if(!isTRUE(as.logical(ctm$continuoustime))) stop(
      'stationary = TRUE is implemented for continuous time models only.',
      call. = FALSE)
    .ctStationaryCheckDrift(ctm)
    t0 <- ctm$pars$matrix %in% c('T0MEANS', 'T0VAR')
    ctm$pars$param[t0] <- NA
    ctm$pars$value[t0] <- ifelse(ctm$pars$matrix[t0] %in% 'T0VAR' &
        ctm$pars$row[t0] == ctm$pars$col[t0], 1, 0)
    ctm$pars$transform[t0] <- NA
    ctm$pars$indvarying[t0] <- FALSE
    # 'FALSE' as a string, for the reason T0VARredundancies() gives.
    for(effect in intersect(paste0(ctm$TIpredNames, '_effect'),
      names(ctm$pars))) ctm$pars[t0, effect] <- 'FALSE'
  }
  ctm$stationary <- as.integer(stationary)


  if(nsubjects == 1 && any(ctm$pars$indvarying[is.na(ctm$pars$value)]==TRUE)){
    # is.null(ctm$fixedrawpopmeans) && is.null(ctm$fixedsubpars) & is.null(ctm$forcemultisubject)) {
    ctm$pars$indvarying <- FALSE
    message('Individual variation not possible as only 1 subject! indvarying set to FALSE on all parameters')
  }

  if(nsubjects == 1 & any(is.na(ctm$pars$value[ctm$pars$matrix %in% 'T0VAR']))){
    # is.null(ctm$fixedrawpopmeans) & is.null(ctm$fixedsubpars) & is.null(ctm$forcemultisubject)) {
    for(ri in 1:nrow(ctm$pars)){
      if(is.na(ctm$pars$value[ri]) && ctm$pars$matrix[ri] %in% 'T0VAR'){
        ctm$pars$value[ri] <- ifelse(ctm$pars$row[ri] == ctm$pars$col[ri], 1e-3, 0)
      }
    }
    message('Free T0VAR parameters fixed to diagonal matrix of 0.001 as only 1 subject - consider appropriateness!')
  }

  if(binomial){
    # A warning, as the other two deprecations are. A message is easy to miss,
    # and this one silently changes `intoverstates` and every indicator's type.
    warning('binomial argument is deprecated -- set manifesttype in the model object to 1 for binary indicators instead. It has set manifesttype=1 for every indicator.', call.=FALSE)
    # It used to set `intoverstates <- FALSE` as well, which is a leftover from
    # when binary data meant sampling the latent states rather than integrating
    # them. The very next check refuses `intoverstates=FALSE` with
    # `optimize=TRUE` -- so the documented shortcut would now put a user
    # straight into an error, under the default `optimize=TRUE`. Setting `manifesttype` directly never did that, and the
    # linearised measurement handles binary indicators with the filter intact:
    # on a three-indicator model it recovers a generating drift of -0.3 as
    # -0.279 and a diffusion of 0.8 as 0.681, both intervals containing the
    # truth.
    ctm$manifesttype[] <- 1
  }

  if(checks){
    if(!optimize && priorscope %in% c('none', 'randomCorr')) warning(
      if(priorscope %in% 'none') 'Sampling with priors=FALSE' else
        "Sampling with priors='randomCorr', which leaves the population sds without a prior",
      ": where the data leave a parameter flat its posterior is improper. ",
      "priors=TRUE puts a prior on every parameter.", call. = FALSE)
    # `intoverstates=FALSE` with `optimize=TRUE` maximises the joint density of
    # the parameters and the innovations that build the states, and that mode is
    # degenerate rather than merely biased: with an innovation per observation
    # the states re-optimise to absorb almost any change in the parameters, so
    # the profile is nearly flat -- its largest eigenvalue measured 0.05 on 15
    # subjects x 6 rows, against 22.8 for a well-determined count parameter --
    # and a joint fit reports a better log likelihood and a better conditioned
    # Hessian for changes that mean nothing. Sampling the same density is sound,
    # which is what the route is for, so this is refused and points there.
    #
    # `estonly` is the way through for someone who wants the mode regardless;
    # the warnings below still say what it does. Not for `fit=FALSE`, which
    # optimises nothing and returns the prepared model, whose joint density is a
    # fair thing to evaluate. A property of the estimator rather than of a
    # backend, so it stands for both.
    if(isTRUE(optimize) && !isTRUE(intoverstates) && isTRUE(fit) &&
        !isTRUE(estonly)) stop(
      'optimize=TRUE with intoverstates=FALSE maximises over the latent states, ',
      'and that joint mode is degenerate rather than an estimate. Use ',
      'optimize=FALSE to sample the states, or intoverstates=TRUE to integrate ',
      'them out. optimcontrol$estonly=TRUE returns the joint mode anyway, ',
      'without standard errors.', call.=FALSE)
    # Maximising over the states rather than integrating them out biases the
    # variance parameters downward -- a variance whose own realisations are
    # being chosen at the same time can always be made to look smaller -- so
    # the joint mode is not the maximum likelihood estimate. That is a property
    # of the estimator and not of a backend, so the warning stands for both.
    # Sampling the same density has no such problem, which is why it points
    # there.
    if(optimize && !intoverstates){
      # Two strengths, because the damage is not uniform. Maximising over the
      # states biases every variance parameter downward, which is the general
      # case. But a parameter that governs how tightly the transition prior
      # constrains the freely chosen states is not merely biased: the joint
      # mode buys cheaper innovations by weakening mean reversion, so such a
      # parameter runs to the flat end of its transform and stays there.
      # Measured on the engine's own fixture: a free DRIFT landed at a raw
      # value of -17 to -19, against a saturation guard at 20, at every sample
      # size from 20 to 120 observations and in six of seven seeds. The
      # engine's count-model test asserts the same runaway independently, on
      # data generated from the model. Fix those parameters and the joint mode
      # is well behaved -- the same fixture with DRIFT fixed converges to a
      # gradient norm of 1e-9 or better -- which is the case this route is for.
      jointrunaway <- ctm$pars$matrix %in% c('DRIFT','DIFFUSION') &
        is.na(ctm$pars$value)
      if(any(jointrunaway)){
        warning(
          'intoverstates=FALSE with optimize=TRUE maximises over the latent ',
          'states rather than integrating them out, and this model leaves ',
          paste0(unique(ctm$pars$matrix[jointrunaway]), collapse=' and '),
          ' free. Those parameters set how tightly the transition prior holds ',
          'the states, so the joint mode weakens them without limit to buy ',
          'cheaper innovations: they run to the end of their transform and ',
          'more data does not help. Treat their estimates as unusable. Fix ',
          'them and free only the means, or use optimize=FALSE to sample the ',
          'states, or intoverstates=TRUE to integrate them out.',
          call.=FALSE)
      } else {
        warning(
          'intoverstates=FALSE maximises over the latent states rather than ',
          'integrating them out, which biases variance parameters downward: ',
          'the joint mode is not the maximum likelihood estimate. Use ',
          'intoverstates=TRUE, or optimize=FALSE to sample the states instead.',
          call.=FALSE)
      }
      # And for a Gaussian indicator with free measurement error it is worse
      # than biased: the joint density is *unbounded*. Send the measurement
      # variance to zero and let the trajectory interpolate the data exactly,
      # and the density diverges -- there is no maximum to find, and an
      # optimiser correctly runs off toward the boundary. Measured on a
      # one-indicator model: MANIFESTVAR reached a raw value of -10.6 with a
      # gradient of 3e9, reported as not converged.
      #
      # Said here rather than left to that non-convergence, which describes
      # the symptom and not the cause. A categorical indicator has no such
      # parameter and is unaffected; so is a Gaussian one whose MANIFESTVAR
      # is fixed. Confirmed degenerate rather than merely biased, so this one
      # is a hard error rather than a warning: there is no fit to return.
      freevar <- ctm$pars$matrix %in% 'MANIFESTVAR' &
        ctm$pars$row == ctm$pars$col & is.na(ctm$pars$value)
      if(any(freevar)){
        gaussian <- ctm$manifesttype[ctm$pars$row[freevar]] %in% 0
        if(any(gaussian)) stop(
          'With intoverstates=FALSE the joint density is unbounded for a ',
          'Gaussian indicator whose MANIFESTVAR is free: the measurement ',
          'variance goes to zero and the latent trajectory interpolates the ',
          'data exactly. The optimiser will run to that boundary and report ',
          'not converged. Fix MANIFESTVAR for ',
          paste(ctm$manifestNames[unique(ctm$pars$row[freevar][gaussian])],
            collapse=', '), ', or use optimize=FALSE.', call.=FALSE)
      }
    }
  }

  # `intoverpop` selects how declared individual differences are handled.
  # TRUE, FALSE and 'auto' keep their existing meanings exactly. The two
  # character forms name the two methods: 'augmented' is the existing
  # state-augmentation, which appends a static latent state per varying
  # parameter and lets the ordinary filter integrate them out; 'laplace'
  # integrates them out per subject instead, so each subject keeps the
  # single-subject state space and the cost of a random effect stops being
  # cubic in the augmented dimension.
  #
  # `intoverpopmethod` carries that choice onwards. `intoverpop` itself stays
  # the logical that drives augmentation everywhere downstream, so the Laplace
  # route reaches the backend with an unaugmented model and every existing
  # `if(intoverpop)` keeps meaning what it meant.
  intoverpopmethod <- 'none'
  # Why 'auto' went the way it did, kept for `fit$args$resolved`; NA when the
  # route was named.
  intoverpopreason <- NA_character_
  # TRUE asks for the random effects to be integrated out, as it always has:
  # the rule 'auto' applies when optimizing, whatever `optimize` is. Before
  # 4.0.0 it meant 'augmented', the only way there was.
  if(isTRUE(intoverpop)){
    auto <- .ctIntOverPopAuto(ctm, backend = backend, optimize = TRUE,
      intoverstates = intoverstates)
    intoverpopmethod <- auto$route
    intoverpopreason <- auto$reason
    intoverpop <- identical(auto$route, 'augmented')
  } else if(is.character(intoverpop)){
    intoverpop <- match.arg(intoverpop[1], c('auto','augmented','laplace','none'))
    if(intoverpop %in% 'auto'){
      # See `.ctIntOverPopAuto()` for the rule and the measurements behind it.
      auto <- .ctIntOverPopAuto(ctm, backend = backend, optimize = optimize,
        intoverstates = intoverstates)
      intoverpopmethod <- auto$route
      intoverpopreason <- auto$reason
      intoverpop <- identical(auto$route, 'augmented')
    } else {
      intoverpopmethod <- intoverpop
      intoverpop <- identical(intoverpopmethod,'augmented')
    }
  }
  intoverpop <- isTRUE(intoverpop)
  if(intoverpop) intoverpopmethod <- 'augmented'

  if(identical(intoverpopmethod,'laplace')){
    # Reached through 'auto' too (a grouping level above the subject), where
    # naming 'laplace' alone blamed a choice the caller never made.
    if(!backend %in% 'julia') stop(
      if(!is.na(intoverpopreason)) paste0("This model needs ",
        "intoverpop='laplace' (", intoverpopreason, "), which ") else
        "intoverpop='laplace' ",
      "requires backend='julia'; the generated Stan model ",
      "does not provide the higher-order derivatives it needs.", call.=FALSE)
    if(!.ctAnyVarying(ctm)) stop(
      "intoverpop='laplace' was requested but no free parameters are marked ",
      "indvarying at any level, so there is nothing to integrate over.",
      call.=FALSE)
  }

  # An outer level on a route that cannot carry it.
  #
  # `.ctModelIntOverPop()` reads the subject level's column and nothing else,
  # so an effect at a study level (`indvarying_study` once resolved, by
  # .ctResolveLevels()) would be dropped without trace and
  # the fit would report a single-level model as though that were what was
  # asked for. Refuse by name instead. There is no reason for a lower level to
  # vary before an upper one does -- a study effect with exchangeable subjects
  # inside it is an ordinary model -- so this is about which route can
  # represent the request, not about which requests are meaningful.
  # A fixed value asked to vary between individuals, which it cannot.
  #
  # `.ctVaryingRows()` treats a cell with a value as fixed whatever its
  # `indvarying` flag says, so neither preparation route augments it and the
  # request simply evaporated -- and RAWPOPVAR used to go on offering a
  # population spread for the parameter anyway. Setting `indvarying` directly
  # on `pars` is how nearly every multilevel model here is written, so this is
  # a reachable mistake rather than a hypothetical one.
  #
  # Named cells only. A number written in a matrix has no parameter name, and
  # the common idiom `model$pars$indvarying <- TRUE` sets the flag on every row
  # including those; naming each fixed LAMBDA and T0VAR cell back at someone
  # who wrote that would be noise. A *named* parameter that also carries a
  # value got there by a deliberate assignment, which is the case worth
  # reporting.
  #
  # A warning rather than an error, matching what the model spec parser does
  # with the same contradiction written in one cell: the value is kept, the
  # individual differences are dropped, and the drop is said out loud.
  fixedvarying <- .ctVaryingColumns(ctm)
  fixedvarying <- fixedvarying[fixedvarying %in% names(ctm$pars)]
  if(length(fixedvarying)){
    flagged <- rep(FALSE, nrow(ctm$pars))
    for(cl in fixedvarying) flagged <- flagged | ctm$pars[[cl]] %in% TRUE
    inert <- flagged & !is.na(ctm$pars$value) & !is.na(ctm$pars$param)
    if(any(inert)) warning(
      'Individual differences were requested for ', 
      paste(unique(ctm$pars$param[inert]), collapse=', '),
      ', which ', if(sum(inert) > 1) 'are' else 'is',
      ' fixed to a value and so cannot vary between individuals. Fitted as ',
      'fixed. Clear the value to estimate ', 
      if(sum(inert) > 1) 'them' else 'it', ' with a population spread.',
      call.=FALSE)
  }

  # A TI effect fixed to a value ('TI1=4.3') is a known offset on its
  # parameter. The julia engine adds it; the stan program has no slot for it
  # and would fit the model without it, so it is refused there by name.
  fixedti <- .ctTipredFixedEffects(ctm)
  if(length(fixedti) && !backend %in% 'julia') stop(
    'Time independent predictor effects fixed to a value (',
    paste(utils::head(fixedti, 4), collapse=', '),
    ') are fitted by backend="julia" only.', call.=FALSE)

  outervarying <- .ctVaryingParams(ctm, .ctOuterVaryingColumns(ctm))
  if(length(outervarying)){
    named <- paste(outervarying, collapse=', ')
    if(!backend %in% 'julia') stop(
      "Random effects above the subject level (", named, ") are represented by ",
      "backend='julia' only; the generated Stan model has one grouping level.",
      call.=FALSE)
    if(intoverpop || (isTRUE(optimize) && !identical(intoverpopmethod,'laplace'))) stop(
      "Random effects above the subject level (", named, ") are integrated out ",
      "by intoverpop='laplace' only -- the augmented route gives carrier states ",
      "to subject level effects and would drop these.", call.=FALSE)
  }

  # Optimizing without integrating over the population distribution is not a
  # combination that means anything. `intoverpop=FALSE` leaves each subject's
  # random effects as free parameters of the objective, so maximizing it
  # maximizes over the effects as well as over the population parameters, and
  # the population variance it lands on is the one that makes those particular
  # effects most likely -- which is zero, or as near as the data allow. Sampling
  # is what handles that model, which is why `intoverpop='auto'` resolves to
  # FALSE exactly when `optimize` is FALSE.
  #
  # It is refused rather than quietly corrected because both readings of the
  # request are plausible -- integrate them, or sample instead -- and guessing
  # would silently answer a different question. Before this it died inside a
  # transform rendering with `invalid format '%.17g'`, which named neither.
  if(isTRUE(optimize) && !intoverpop && identical(intoverpopmethod,'none') &&
      .ctAnyVarying(ctm)) stop(
    "intoverpop=FALSE with optimize=TRUE leaves each subject's random effects ",
    "as free parameters to be maximized over, which drives the population ",
    "variance to zero rather than estimating it. Use intoverpop='augmented' ",
    "or 'laplace' to integrate them out, or optimize=FALSE to sample.",
    call.=FALSE)

  # if(optimize && !intoverpop && any(ctm$pars$indvarying[is.na(ctm$pars$value)]) &&
  #     is.null(ctm$fixedrawpopchol) && is.null(ctm$fixedsubpars)){
  #   intoverpop <- TRUE
  #   message('Setting intoverpop=TRUE to enable optimization of random effects...')
  # }

  # if(intoverpop==TRUE && !any(ctm$pars$indvarying[is.na(ctm$pars$value)])) {
  #   # message('No individual variation -- disabling intoverpop switch');
  #   intoverpop <- FALSE
  # }

  # Individual variation on a variance cell is only partially identified under
  # the augmented route, and it is not backend specific -- the filter is the
  # same on both, so this sits ahead of the stan / julia dispatch deliberately.
  #
  # The augmented route makes a random effect a static latent state with LAMBDA
  # zero on it. A mean-affecting cell (MANIFESTMEANS, CINT, T0MEANS, DRIFT)
  # reaches the observation mean, so the Kalman update moves that state and its
  # population variance is informed. A DIFFUSION or MANIFESTVAR cell is built
  # from the state's mean only: the observation mean function's Jacobian with
  # respect to it is zero, so the update can never move it, and the data learns
  # about the effect solely through its correlation with states the filter can
  # update. The covariance sd_i x sd_j x corr_ij is then identified and its
  # split into an sd and correlations is a ridge -- bounded below, unbounded
  # above. See review/RANDOMEFFECTS-partial-identification-2026-09-07.md.
  #
  # This is deliberately coarse and will occasionally fire where the parameter
  # IS identified, because a parameter also referenced inside a mean-affecting
  # expression -- a DRIFT string, say -- is identified and has no DRIFT row in
  # `$pars` to reveal it. That false positive is accepted rather than chased
  # with a substring scan of the character cells. The accurate distinction is
  # computable, not pattern-matched: `.ctBackendIdentifiability()` returns each
  # flat direction's loadings, so one helper checking that the covariance
  # functionals are orthogonal to the direction would serve both this pre-fit
  # site and the post-fit `.ctBackendIdentifyWarn()`. That is where an exact
  # version belongs.
  #
  # It must run before `.ctModelIntOverPop()` below, which clears `indvarying`
  # on the cells it rewrites into state references, and MANIFESTVAR rows for
  # non-Gaussian indicators are excluded because they are fixed further down
  # (the `errfix` block) and warning about them would contradict that message.
  if(checks && intoverpop){
    revarpars <- ctm$pars$matrix %in% c('DIFFUSION','MANIFESTVAR') &
      ctm$pars$indvarying & is.na(ctm$pars$value)
    if(any(revarpars)) revarpars <- revarpars &
      !(ctm$pars$matrix %in% 'MANIFESTVAR' &
          ctm$pars$row %in% which(ctm$manifesttype > 0 & ctm$manifesttype != 4))
    if(any(revarpars)) warning(
      "Individual variation on DIFFUSION or MANIFESTVAR is only partially ",
      "identified with intoverpop='augmented': the random effect enters only ",
      "the predicted covariance, so the filter never updates it, and the data ",
      "determines its covariance with the other random effects but not the ",
      "split of that covariance into a standard deviation and correlations. ",
      "Affected: ", paste(unique(ctm$pars$param[revarpars]), collapse=', '),
      ". intoverpop='laplace' identifies these separately.", call.=FALSE)
  }

  ctm <- ctModel0DRIFT(ctm, ctm$continuoustime) #offset 0 drift
  ctm$pars <- ctModelStatesAndPARS(ctm$pars,statenames = ctm$latentNames,tdprednames=ctm$TDpredNames,tiprednames=ctm$TIpredNames) #replace latent states and PARS with state and PAR[] refs, need this early because we rely on [] detection
  if(any(grepl('tipreds[', ctm$pars$param, fixed=TRUE)) && !backend %in% 'julia') stop(
    'A cell reads a time independent predictor (',
    paste(utils::head(unique(ctm$pars$param[grepl('tipreds[', ctm$pars$param, fixed=TRUE)]), 3),
      collapse=', '), '), which backend="julia" evaluates and the stan backend does not.',
    call.=FALSE)

  # A reduced-rank population covariance, written as a regression of the
  # remaining random effects on a basis of them. Three steps around the
  # augmentation rather than a second path through it: resolve the basis while
  # `indvarying` still says which parameters vary, clear the flag on the
  # regressed ones so `.ctModelIntOverPop()` gives carrier states to the basis
  # alone, then write the regression expressions once those states exist. The
  # rewrite deliberately lands *before* the second `ctModelStatesAndPARS()`
  # call below, so the new mean and coefficient parameters can be introduced as
  # plain labels and turned into `PARS[r,c]` references by the machinery that
  # already does exactly that. See R/ctPopRegression.R.
  # `poprank` defaults to 'auto', so the places it does not apply have to be
  # inapplicable rather than errors: only a rank the user asked for is refused.
  #
  # One structure, a loading matrix on standardised dimensions, reached two
  # ways because the routes build their population covariance in two places.
  # Under `'laplace'` and `'none'` -- which prepare the same engine
  # specification, one integrating the effects and the other sampling them --
  # it is built in the engine, per level, by one piece of code (`laplacerank`
  # below, read by `.ctJuliaLevelRank()`). Under `'augmented'` the effects are
  # carrier states in the filter, so the rank has to be written into the model
  # around `.ctModelIntOverPop()` (R/ctPopRegression.R); that route has only the
  # subject level, so the rewrite's limit to the innermost level costs it
  # nothing. `'none'` used to take a separate rewrite into regression
  # coordinates -- the form `203e223f` measured as reaching the optimum from 1
  # of 12 matched starts against 8 of 12 for loadings, and which tied levels
  # together through one coefficient.
  #
  # What the restriction *means* differs between them, and only the message says
  # so: on the augmented route it removes coordinates the filter cannot see and
  # costs no likelihood, while under laplace those coordinates are identified
  # and removing them is an approximation. The same restriction, different
  # claims.
  # The rank may be stated on the model instead, as `model$poprank`, which is
  # where it belongs for anyone who thinks of it as part of the specification --
  # `indvarying` is set that way in nearly every multilevel model in the tests,
  # so the idiom is already the house one. Not a `ctModel()` argument: that
  # list goes through `ctModelConvertOMX()` and an unknown field there is a
  # risk for no gain, where a plain assignment onto the returned model works and
  # survives.
  #
  # A rank passed to `ctFit()` wins, because an argument at the call site is the
  # more specific statement of the two; the model's value is used only when the
  # argument was left at its default.
  #
  # A rank stated on the model was asked for just as much as one passed here, so
  # from this point it counts as explicit: applied under laplace and 'none',
  # refused by name where it cannot apply. Before this the flag still read "left
  # at its default", and a laplace fit took the model's rank and then ignored
  # it, running at full rank without a word. The `args` capture above already
  # holds the flag as it was, so `ctFitUpdate()` replays the argument only when
  # one was passed, and the model it refits still carries its own rank.
  if(!poprankexplicit && !is.null(ctm[['poprank']])){
    poprank <- ctm[['poprank']]
    poprankexplicit <- TRUE
  }

  # Under `laplace` a rank restricts each level's population covariance where
  # that covariance is actually built -- in the engine, as a loading matrix --
  # rather than by rewriting the model into a basis and a set of regressions.
  # Two reasons it has to be done that way here and not the augmented way.
  #
  # It reaches every level. The rewrite works from `pars$indvarying`, which is
  # the innermost level alone, so on a burst/subject/study model it can restrict
  # the burst covariance and cannot touch the study one -- and the study level,
  # with the fewest groups, is the one that needs it.
  #
  # And it does not entangle the levels. A regressed effect is written into its
  # cell as `(p + beta * b)`, so if the basis effect `b` also varies at another
  # level then `p` inherits *that* level's deviation of `b` through the same
  # `beta`. One coefficient tying two levels together is not the model anyone
  # asked for, and nothing about the resulting fit would look wrong.
  #
  # A loading matrix per level has neither problem: each level's deviation is
  # `L_level * u_level` with its own `u`, and the levels stay independent.
  #
  # Laplace takes this route alone, whether or not it ends up restricting
  # anything: a request that resolves to full rank at every level is answered
  # by full rank, not by handing the rank on to the augmented machinery below,
  # which reads it as a single number for the innermost level.
  #
  # Not by default, and this is the whole reason a request is distinguished
  # from the default. On the augmented route the coordinates `'auto'` removes
  # cannot be identified, so removing them costs nothing and is a good default.
  # Under laplace they *are* identified, and measured on a 250 x 50 design the
  # same restriction costs 48 log likelihood units and takes the fit to the
  # boundary -- basis sd to zero with the coefficient to -612. A default that
  # does that to a user who chose the route precisely because it identifies
  # these things would be indefensible, so here it has to be asked for. The
  # same holds under 'none', below.
  levelnames <- c(ctm$subjectIDname, ctm$groupIDnames)
  laplacerank <- NULL
  # Blocks of shared dimensions (`popblocks`, R/ctPopBlocks.R) are columns of
  # the engine's loading matrix, so they ask for this route as a rank does.
  blocksasked <- .ctPopblocksRequested(popblocks)
  if(blocksasked && (!identical(backend,'julia') ||
      !intoverpopmethod %in% c('laplace','none'))) stop("popblocks requires ",
    "backend='julia' with intoverpop='laplace' or sampling (intoverpop='none'): ",
    "a block is a set of columns of the engine's loading matrix, which the ",
    "augmented route and stan do not build.", call.=FALSE)
  rankasked <- poprankexplicit && !(length(poprank)==1 && is.na(poprank))
  enginerank <- intoverpopmethod %in% c('laplace','none') && (rankasked || blocksasked)
  if(enginerank && !identical(backend,'julia')) stop(
    "poprank requires backend='julia'.", call.=FALSE)
  if(enginerank && !any(ctm$pars$indvarying[is.na(ctm$pars$value)]) &&
      !.ctAnyVarying(ctm, .ctOuterVaryingColumns(ctm))) stop(
    "poprank restricts the population covariance, so it needs a model with ",
    "individually varying parameters.", call.=FALSE)
  if(enginerank){
    levelcolumns <- stats::setNames(c('indvarying', if(length(ctm$groupIDnames))
      paste0('indvarying_', ctm$groupIDnames)), levelnames)
    asked <- if(rankasked) .ctPoprankByLevel(poprank, levelnames) else list()
    laplacerank <- vapply(names(asked), function(level){
      value <- asked[[level]]
      if(!identical(value, 'auto')) return(value)
      # `'auto'` keeps the meaning it has on the augmented route -- how many of
      # this level's effects reach the observation mean -- resolved per level
      # rather than once for the model. A level with none keeps its full
      # covariance rather than a rank of zero, which would be no variation.
      roles <- .ctPopEffectRoles(ctm$pars, column=levelcolumns[[level]])
      counted <- if(nrow(roles)) as.integer(sum(roles$mean)) else NA_integer_
      if(is.na(counted) || counted < 1L) NA_integer_ else counted
    }, integer(1L))
    laplacerank <- laplacerank[!is.na(laplacerank)]
    if(!length(laplacerank)) laplacerank <- NULL
    # Blocks (`popblocks`): a level's rank is the sum of its blocks'. A number
    # asked for the same level must agree with it; `'auto'`, the default
    # written out, defers to the blocks.
    laplaceblocks <- if(blocksasked) .ctPopblocksResolve(ctm, popblocks,
      levelnames, levelcolumns) else list()
    for(level in names(laplaceblocks)){
      total <- sum(vapply(laplaceblocks[[level]], function(b) b$rank, integer(1L)))
      if(level %in% names(laplacerank) && laplacerank[[level]] != total &&
          !identical(asked[[level]], 'auto')) stop(
        "poprank for level '", level, "' is ", laplacerank[[level]], " but its ",
        "popblocks span ", total, " dimensions. Drop poprank for that level: ",
        "its blocks give its rank.", call.=FALSE)
      laplacerank[level] <- as.integer(total)
    }
    if(length(laplaceblocks)) ctm$laplaceblocks <- laplaceblocks
    # Own variances (`popunique`, R/ctPopUnique.R), resolved here because the
    # T0MEANS floor below counts them: an effect with an own variance has a
    # dimension of its own whatever the rank.
    laplaceunique <- if(.ctPopuniqueRequested(popunique)) .ctPopuniqueResolve(ctm,
      popunique, levelnames, levelcolumns, laplacerank,
      patterned = names(laplaceblocks)) else list()
    if(length(laplaceunique)) ctm$laplaceunique <- laplaceunique
    # An effect in no block shares no dimension, so without an own variance
    # it would not vary.
    for(level in names(laplaceblocks)){
      inblock <- unique(unlist(lapply(laplaceblocks[[level]], function(b) b$param)))
      roles <- .ctPopEffectRoles(ctm$pars, column=levelcolumns[[level]])
      alone <- setdiff(setdiff(roles$param, inblock), laplaceunique[[level]])
      if(length(alone)) stop("popblocks at level '", level, "' leaves ",
        paste(alone, collapse=', '), " in no block and without an own variance,",
        " so ", if(length(alone) > 1) "they" else "it", " would not vary. Name ",
        if(length(alone) > 1) "them" else "it", " in a block, or give ",
        if(length(alone) > 1) "them own variances" else "it an own variance",
        " with popunique.", call.=FALSE)
    }
    # An innermost-level varying T0MEANS whose latent has a free T0VAR takes
    # its initial variance from this covariance instead (`T0VARredundancies()`,
    # below, disables those T0VAR cells), so a rank below their number fixes
    # the initial state along a direction: on the ?ctFit example poprank 1 cost
    # 4447 log likelihood units and reported convergence. A fixed T0VAR keeps
    # its own variance, and there a lower rank is an approximation like any
    # other. The augmented route refuses the same case.
    subjectrank <- laplacerank[ctm$subjectIDname]
    t0cells <- ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$indvarying %in% TRUE &
      is.na(ctm$pars$value)
    freevar <- ctm$pars$matrix %in% 'T0VAR' & is.na(ctm$pars$value) &
      ctm$pars$row == ctm$pars$col
    t0vary <- unique(as.character(ctm$pars$param[t0cells &
      ctm$pars$row %in% ctm$pars$row[freevar]]))
    t0shared <- setdiff(t0vary, laplaceunique[[ctm$subjectIDname]])
    if(length(subjectrank) && isTRUE(subjectrank < length(t0shared))) stop(
      "poprank ", subjectrank, " is below the ", length(t0shared),
      " individually varying T0MEANS without an own variance (",
      paste(t0shared, collapse=', '), "): ",
      "their latents take their initial variance from the population ",
      "covariance, so a lower rank fixes the initial state along a direction. ",
      "Use poprank=", length(t0shared), " or more, poprank=NA, or give them ",
      "own variances with popunique.", call.=FALSE)
    ctm$laplacerank <- laplacerank
    # Said once, for the levels it reduces: on these routes the coordinates a
    # rank drops are identified, so the restriction is an approximation.
    reduced <- vapply(names(laplacerank), function(level){
      k <- nrow(.ctPopEffectRoles(ctm$pars, column=levelcolumns[[level]]))
      m <- length(laplaceunique[[level]])
      if(k > laplacerank[[level]]) paste0(level, ' ', laplacerank[[level]],
        ' of ', k, if(m) paste0(' plus ', m, ' own variance', if(m > 1) 's'))
        else NA_character_
    }, character(1L))
    reduced <- reduced[!is.na(reduced)]
    if(length(reduced)) message('poprank: population covariance reduced to rank ',
      paste(reduced, collapse=', '), " -- an approximation under intoverpop='",
      intoverpopmethod, "'; poprank=NA estimates it in full.")
  }

  # `popunique` is built where the engine builds the loadings, so it needs that
  # route and a rank asked for there; refused by name anywhere else rather
  # than accepted and ignored.
  if(.ctPopuniqueRequested(popunique) && !enginerank){
    if(!identical(backend,'julia') || !intoverpopmethod %in% c('laplace','none'))
      stop("popunique requires backend='julia' with intoverpop='laplace' or ",
        "sampling (intoverpop='none'): an own variance is a column of the ",
        "engine's loading matrix, which the augmented route and stan do not ",
        "build.", call.=FALSE)
    stop("popunique needs a reduced rank at the levels it selects: give ",
      "poprank, e.g. poprank=c(", ctm$subjectIDname, "=2), or 0 for ",
      "independent random effects.", call.=FALSE)
  }

  popregression <- NULL
  if(intoverpop && !(length(poprank)==1 && is.na(poprank))){
    if(!identical(backend,'julia')){
      if(poprankexplicit) stop("poprank requires backend='julia'.", call.=FALSE)
    } else {
      # The augmented route has the subject level only: the rewrite works from
      # `pars$indvarying`. So a per-level rank may name that level and no
      # other, and an unnamed one is a single number for it.
      asked <- .ctPoprankByLevel(poprank,
        if(is.null(names(poprank))) ctm$subjectIDname else levelnames)
      outer <- setdiff(names(asked), ctm$subjectIDname)
      if(length(outer)) stop("poprank names level ",
        paste(outer, collapse=', '), ", but under intoverpop='augmented' ",
        "a rank restricts the subject level ('", ctm$subjectIDname,
        "') alone. A rank per level needs intoverpop='laplace', or ",
        "sampling (optimize=FALSE).", call.=FALSE)
      poprank <- asked[[ctm$subjectIDname]]
      popregression <- .ctPopRegressionSpec(ctm$pars, poprank,
        explicit=poprankexplicit, model=ctm)
      if(!is.null(popregression)) ctm <- .ctPopRegressionDemote(ctm, popregression)
    }
  }
  if(intoverpop)   ctm <- .ctModelIntOverPop(ctm) #extend system matrices for individual differences
  if(!is.null(popregression)){
    ctm <- .ctPopRegressionRewrite(ctm, popregression)
    popregression <- ctm$popregression
    if(!is.null(popregression)) message(.ctPopRegressionMessage(popregression))
  }

#   #check this *after* replacing PARS references as needed
#   if(any(duplicated(ctm$pars$param[ctm$pars$matrix %in% 'T0MEANS' &
#       !grepl('[',ctm$pars$param,fixed=TRUE) &
#       !is.na(ctm$pars$param)]))) stop(paste0(
#         'Unfortunately, duplicate T0MEANS parameters must be specified via inclusion of additional PARS matrix in ctModel: e.g.,
# ctModel(... #regular model code
# PARS=c("t0mPar||TRUE"), #specify an additional parameter called t0mPar, with random effects
# T0MEANS=c("t0mPar","t0mPar"), #insert this parameter into the T0MEANS matrix as many times as needed.
# ... #regular model code)
# '))

  #jacobian addition
  ctm$jacobian <- try(ctJacobian(ctm))
  if('try-error' %in% class(ctm$jacobian)) ctm$jacobian <- ctJacobian(ctm,simplify=FALSE)
  # ctm$jacobian <- unfoldmats( #replaces matrix references with base parameter
  #   c(listOfMatrices(ctm$pars),ctm$jacobian))
  ctm$jacobian <- ctm$jacobian[names(.ctMatricesList()$jacobian)]
  jl <- ctModelUnlist(ctm$jacobian,names(ctm$jacobian))
  jl <- jl[apply(jl,1,function(x) any(!is.na(x))),] #clean up messy leftovers of NA's
  jl2 <- as.data.frame(rbind(data.table(ctm$pars[1,]),data.table(jl),fill=TRUE))[-1,]


  for(i in 1:nrow(jl2)){ #copy base parameter transforms etc to jacobian when needed
    if(!is.na(jl2$param[i]) && jl2$param[i] %in% ctm$pars$param){
      jl2[i, !colnames(jl2) %in% list('matrix','row','col')] <-
        ctm$pars[which(ctm$pars$param %in% jl2$param[i])[1], !colnames(ctm$pars) %in% list('matrix','row','col')]
    }
  }

  ctm$pars <- rbind(ctm$pars,jl2)

  ctm$pars <- ctModelStatesAndPARS(ctm$pars,statenames = ctm$latentNames,tdprednames=ctm$TDpredNames,tiprednames=ctm$TIpredNames) #replace any new state and par refs with square bracket refs

  # The transform text as the model states it, keyed by cell, kept before
  # `ctModelTransformsToNum` replaces it with four numbers recovered from it by
  # a grid search. That search scores candidates by squared residual, so it
  # cannot see a constant small enough not to move the residual: every variance
  # diagonal carries a `1e-10` floor -- `1e-10 + 5 * log1p_exp(2 * param)` --
  # and the `round(x, 6)` that follows finishes it off. DRIFT's floor is `1e-06`
  # and survives, which is why only the variances lost theirs.
  #
  # A floorless variance reaches *exactly* zero, since `log1p_exp` is exactly
  # zero once `1 + exp(x)` rounds to one, and a zero variance is then divided
  # by and differentiated through. The julia backend reads this rather than
  # re-rendering the numbers; see `.ctJuliaParameterTable`.
  #
  # Keyed by cell and kept *off* `ctm$pars`, because parts of the Stan pipeline
  # read that frame's columns by position -- carrying it there as a character
  # column turned Stan's `pop_CINT` into NaN.
  if(is.null(ctm$transformtext) && is.character(ctm$pars$transform)){
    # A cell whose transform column already parses as a bare number (the
    # output of an earlier ctModelTransformsToNum() call, which
    # ctEBadjustModel() makes before handing its adjusted model to ctFit --
    # see R/ctEmpiricalBayesFit.R) is not descriptive source text: it is a
    # transform *code*, with the real expression already discarded, and this
    # column is character-typed regardless because it mixes such codes with
    # genuine expression strings. Recording a bare code as `text` fed
    # .ctJuliaParameterTable() a rendered transform of literally "1", with no
    # param[] reference left for the engine's adjoint to find -- see the
    # comment there. NA leaves that cell to the numeric multiplier/meanscale
    # reconstruction, which is what a bare code needs regardless of whether it
    # arrived that way from the user or from an earlier numeric reduction.
    bareCode <- !is.na(suppressWarnings(as.numeric(ctm$pars$transform)))
    text <- as.character(ctm$pars$transform)
    text[bareCode] <- NA_character_
    ctm$transformtext <- data.frame(
      matrix = as.character(ctm$pars$matrix),
      row = as.integer(ctm$pars$row), col = as.integer(ctm$pars$col),
      text = text, stringsAsFactors = FALSE)
  }
  ctm <- ctModelTransformsToNum(ctm)

  ctm$pars <- .ctModelCleanctspec(ctm$pars)

  ctm <- T0VARredundancies(ctm)

  if(!all(ctm$pars$transform[!is.na(suppressWarnings(as.integer(ctm$pars$transform)))] %in% c(0,1,2,3,4))) stop('Unknown transform specified -- integers should be 0 to 4')

  #fix binary manifestvariance

  if(any(ctm$manifesttype %in% 2) && !identical(backend, 'julia')){
    stop('Ordinal manifest variables (manifesttype 2) need backend="julia". ',
      'The stan model has no ordinal measurement, so its thresholds would ',
      'never be read; the julia filter integrates the observation over the ',
      'latent instead. Ordinal variable(s): ',
      paste(ctm$manifestNames[ctm$manifesttype %in% 2], collapse=', '), '.',
      call.=FALSE)
  }
  # Same reason as ordinal, and stated separately because the reason a stan fit
  # would be wrong is different: stan has no Poisson measurement at all, so it
  # would treat a count as a Gaussian observation of the log rate and return
  # numbers that look entirely reasonable.
  if(any(ctm$manifesttype %in% 3) && !identical(backend, 'julia')){
    stop('Count manifest variables (manifesttype 3) need backend="julia". ',
      'The stan model has no count measurement, so the observation would be ',
      'treated as Gaussian; the julia filter integrates it against the ',
      'predicted state under a Poisson log link instead. Count variable(s): ',
      paste(ctm$manifestNames[ctm$manifesttype %in% 3], collapse=', '), '.',
      call.=FALSE)
  }
  # Censored is refused on stan for the same reason as the others: stan has no
  # censored measurement, so a value pinned at a limit would be treated as an
  # ordinary observation of that number and the censoring simply ignored.
  if(any(ctm$manifesttype %in% 4) && !identical(backend, 'julia')){
    stop('Censored manifest variables (manifesttype 4) need backend="julia". ',
      'The stan model has no censored measurement, so a value at its limit ',
      'would be read as an ordinary observation; the julia filter integrates ',
      'the censored likelihood against the predicted state instead. Censored ',
      'variable(s): ',
      paste(ctm$manifestNames[ctm$manifesttype %in% 4], collapse=', '), '.',
      call.=FALSE)
  }
  # The stan programs have no asymptotes in their binary measurement, so a stan
  # fit would quietly be the two parameter model. Refused on the matrix being
  # there at all, rather than on its cells, so the message does not depend on
  # what they hold.
  if(any(ctm$pars$matrix %in% 'ASYMPTOTES') && !identical(backend, 'julia')){
    stop('ASYMPTOTES need backend="julia": the stan model has no asymptotes ',
      'in its binary measurement, so it would fit the two parameter logistic.',
      call.=FALSE)
  }
  .ctCheckAsymptotePars(ctm$pars, ctm$manifesttype)
  .ctCheckResponsePars(ctm$pars, ctm$manifesttype)

  if(any(ctm$manifesttype > 0)){ #if any non continuous variables, (with free parameters)...
    # Binary and ordinal only. For those two a normal term on the linear
    # predictor is not identified rather than merely unwanted: with a probit
    # link E_e Phi(lambda'x + mu + e) is exactly Phi((lambda'x + mu)/sqrt(1+v)),
    # so it is absorbed into the loadings and thresholds, and the logistic link
    # differs only in that the absorption is approximate. Fixing it is the only
    # coherent thing to do.
    #
    # Censored is excluded because that entry is the scale of the Gaussian
    # inside the limits. Count is excluded because the Poisson's variance is
    # locked to its mean, so the term is not absorbable anywhere: it shifts the
    # mean by v/2, which MANIFESTMEANS takes, and multiplies the variance by a
    # factor nothing else in the model can produce. It is the only parameter
    # that moves a count's variance-to-mean ratio, and leaving it out does not
    # make the data equidispersed -- it makes DIFFUSION absorb the difference,
    # which is the parameter the model is usually for.
    deterministic <- which(ctm$manifesttype %in% c(1L, 2L))
    errfix <- which(ctm$pars$matrix %in% 'MANIFESTVAR' &
        (ctm$pars$row %in% deterministic |
            ctm$pars$col %in% deterministic) &
        is.na(suppressWarnings(as.numeric(
          ctm$pars$value))))

    if(length(errfix) > 0){
      message('MANIFESTVAR fixed to 0 for binary / ordinal indicators: the ',
        'measurement link already carries the randomness.')
      ctm$pars$value[errfix] <- 1e-5
      ctm$pars[errfix,c('param','transform','multiplier','offset','meanscale','inneroffset','sdscale')] <- NA
      ctm$pars$indvarying[errfix] <- FALSE
    }

    # A count's dispersion is scalar: each count row is updated on its own,
    # sequentially and conditionally independent of the others given the state,
    # so there is no off-diagonal for it to be correlated through and the
    # engine reads only the diagonal. A free off-diagonal would therefore be
    # accepted here and ignored there, which is the failure this fixes rather
    # than reports -- the cell has one right value and it is zero.
    countrows <- which(ctm$manifesttype %in% 3L)
    offfix <- which(ctm$pars$matrix %in% 'MANIFESTVAR' &
        ctm$pars$row != ctm$pars$col &
        (ctm$pars$row %in% countrows | ctm$pars$col %in% countrows) &
        is.na(suppressWarnings(as.numeric(ctm$pars$value))))
    if(length(offfix) > 0){
      message('Fixing free off-diagonal MANIFESTVAR parameters for count ',
        'indicators to zero -- a count is scored one row at a time, so it has ',
        'no correlated measurement error.')
      ctm$pars$value[offfix] <- 0
      ctm$pars[offfix,c('param','transform','multiplier','offset','meanscale','inneroffset','sdscale')] <- NA
      ctm$pars$indvarying[offfix] <- FALSE
    }

    # A *fixed* non-zero variance on a binary indicator is left alone by the
    # block above, and it is almost certainly a mistake: the Bernoulli link
    # already supplies the randomness, so an extra measurement variance adds
    # noise on top of noise. Free ones are fixed silently because there is a
    # right answer; a deliberate value is the user's, so it is questioned
    # rather than overwritten.
    # What the approximation costs, said once, because it is invisible
    # otherwise and it is large.
    #
    # A binary observation is handled by a moment-matched Gaussian update: the
    # predicted probability from `inv_logit`, a finite-difference Jacobian, and
    # `ycov = Jy etacov Jy' + p(1-p)`. The covariance update that follows,
    # `etacov -= K Jy etacov`, is the standard EKF one and ignores the link's
    # curvature, so posterior uncertainty is understated and the understatement
    # compounds over time steps. The filter ends up believing the latent is
    # pinned down, its innovations shrink, and less process noise is needed to
    # explain them.
    #
    # Measured on a one-latent process with true DIFFUSIONcov 0.16, 50
    # subjects, 50 timepoints, observed only through binary indicators:
    #
    #   5 indicators   0.114   (0.71 of truth)
    #  10 indicators   0.107   (0.67)
    #  30 indicators   0.073   (0.46)
    #
    # DRIFT and CINT survive as long as the latent has few binary indicators;
    # the bias gets *worse* with more of them, which is the signature of a
    # systematic error rather than sampling noise -- more data makes it more
    # confidently wrong. A latent that also has a continuous indicator is
    # unaffected in the same fit.
    #
    # The sharpest evidence is `test-ctRaschExampleTest.R`, which fits the same
    # data twice -- linearised, and sampling the states exactly by HMC -- and
    # compares. The item parameters agree to within 0.012. The dynamics do not:
    #
    #   diff  (process noise)   0.254  linearised   0.373  exact
    #   drift_eta1             -0.095  linearised  -0.935  exact
    #
    # A drift of -0.095 against -0.935 is not a small bias: the linearised fit
    # describes a near random walk where the exact one finds strong mean
    # reversion. That test has been failing, and it is this.
    # Backend specific, because they no longer do the same thing. The julia
    # engine integrates the Bernoulli observation against the predicted state
    # (a mode-centred Gauss-Hermite rule on the scalar linear predictor); the
    # stan model moment-matches it to a Gaussian and linearises. On 40
    # replications with 50 subjects and 12 timepoints, RMSE for a true
    # DIFFUSION of 0.8:
    #
    #   indicators   julia   stan
    #            3   0.134   0.449
    #           10   0.075   0.297
    #           30   0.045   0.229
    #
    # Paired on the same data, every cell favours julia (p <= 0.014). Note the
    # shape of it: stan's *bias* is modest, its RMSE is three to five times
    # worse -- the linearised estimate is unstable rather than uniformly low,
    # which is why single-dataset comparisons looked so different from each
    # other.
    # Nothing is said for the julia backend. That it integrates rather than
    # linearises is a property of the engine the user cannot act on, and the
    # threshold parameterisation is now reported readably by summary()'s
    # Ordinal thresholds section rather than explained at every fit. The stan
    # warning below stays because it is a caution about estimates being
    # unreliable, which is actionable.
    if(!backend %in% 'julia'){
      message('Binary indicators use a linearised (moment-matched Gaussian) ',
        'measurement update on the stan backend, which makes DRIFT and ',
        'especially DIFFUSION unreliable for a latent seen only through them. ',
        'Over 40 replications with a true diffusion of 0.8, RMSE was 0.45 with ',
        '3 indicators and 0.23 with 30, against 0.13 and 0.04 for ',
        "backend='julia', which integrates the observation instead. Prefer ",
        'the julia backend for binary data, or treat process noise as ',
        'indicative.')
    }

    # `deterministic` rather than every non-Gaussian type, so that censored is
    # excluded here too: its MANIFESTVAR entry is the standard deviation of the
    # Gaussian inside the limits, so a fixed non-zero value there is the model
    # working as specified rather than a mistake to warn about.
    linkrows <- which(ctm$pars$matrix %in% 'MANIFESTVAR' &
        ctm$pars$row %in% deterministic &
        ctm$pars$row == ctm$pars$col)
    stated <- linkrows[!is.na(ctm$pars$value[linkrows]) &
        abs(ctm$pars$value[linkrows]) > 1e-4]
    if(length(stated)){
      warning('MANIFESTVAR is fixed to a non-zero value for indicator',
        if(length(stated) > 1) 's ' else ' ',
        paste(ctm$manifestNames[ctm$pars$row[stated]], collapse=', '),
        '. A binary or ordinal indicator gets its randomness from its ',
        'measurement link -- the Bernoulli link for binary, the cumulative ',
        'logit for ordinal -- so this adds measurement noise on top of it, and ',
        'is not separately identified from the loadings and thresholds. Set it ',
        'to 0 unless that is meant.', call.=FALSE)
    }}

  ctm$modelmats <- .ctModelMatSetup(ctm) #slow!
  ctm <- .ctCalcsList(ctm,save=saveComplexPars) #get extra calculations and adjust model spec as needed???

  #store values in ctm
  ctm$intoverpop <- as.integer(intoverpop)
  ctm$nlatentpop <- as.integer(ifelse(ctm$intoverpop ==1, max(ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS']),  ctm$n.latent))
  ctm$intoverstates <- as.integer(intoverstates)
  ctm$priors <- as.integer(priors)
  ctm$stationary <- as.integer(stationary)
  ctm$nlcontrol <- nlcontrol

  list(ctm = ctm, intoverpop = intoverpop, intoverpopmethod = intoverpopmethod,
    intoverpopreason = intoverpopreason, laplacerank = laplacerank,
    popregression = popregression)
}


# nlcontrol with its defaults filled in, as every fit has it.
#' @keywords internal
.ctNlcontrolDefaults <- function(nlcontrol, backend) {
  if(is.null(nlcontrol$maxtimestep)) nlcontrol$maxtimestep = 999999
  if(is.null(nlcontrol$Jstep)) nlcontrol$Jstep = 1e-6
  # nsubsteps = 'auto' lets the julia engine choose the number of prediction
  # substeps per observation interval from how nonlinear that interval turns
  # out to be (ctsem_auto_substeps in the engine); maxtimestep stays a ceiling
  # on the step. substeptol is the largest acceptable linearisation error, as
  # a fraction of the predicted state standard deviation.
  if(!is.null(nlcontrol$nsubsteps)) {
    if(!identical(nlcontrol$nsubsteps, 'auto')) stop("nlcontrol$nsubsteps must be NULL or 'auto'")
    if(!backend %in% 'julia') stop("nlcontrol$nsubsteps = 'auto' requires backend = 'julia'")
  }
  if(is.null(nlcontrol$substeptol)) nlcontrol$substeptol = 0.01
  if(is.null(nlcontrol$maxsubsteps)) nlcontrol$maxsubsteps = 64
  # transition = 'euler' replaces the exponential (locally linearised) step of
  # the state-explicit path (intoverstates = FALSE) with plain Euler-Maruyama,
  # a reference that shares no approximation with the filter. It needs a fine
  # maxtimestep. No effect on the filter itself.
  if(is.null(nlcontrol$transition)) nlcontrol$transition = 'exponential'
  if(!nlcontrol$transition %in% c('exponential', 'euler')) stop("nlcontrol$transition must be 'exponential' or 'euler'")
  if(nlcontrol$transition == 'euler' && !backend %in% 'julia') stop("nlcontrol$transition = 'euler' requires backend = 'julia'")
  nlcontrol
}
