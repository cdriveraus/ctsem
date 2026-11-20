# Package index

## Specify a model

- [`ctModel()`](https://cdriveraus.github.io/ctsem/reference/ctModel.md)
  : Define a ctsem model
- [`ctModelLatex()`](https://cdriveraus.github.io/ctsem/reference/ctModelLatex.md)
  : Generate and optionally compile latex equation of subject level
  ctsem model.
- [`ctModelMatrices()`](https://cdriveraus.github.io/ctsem/reference/ctModelMatrices.md)
  [`` `ctModelMatrices<-`() ``](https://cdriveraus.github.io/ctsem/reference/ctModelMatrices.md)
  : Matrix view for ctStanModel objects
- [`ctPar()`](https://cdriveraus.github.io/ctsem/reference/ctPar.md) :
  Write a parameter specification
- [`ctCov()`](https://cdriveraus.github.io/ctsem/reference/ctCov.md) :
  Write a covariance matrix as the cells a ctsem model reads
- [`ctModelHigherOrder()`](https://cdriveraus.github.io/ctsem/reference/ctModelHigherOrder.md)
  : Raise the order of a ctsem model's latent processes.
- [`ctModelSlowComponent()`](https://cdriveraus.github.io/ctsem/reference/ctModelSlowComponent.md)
  : Add slow processes that drive a model's latent processes
- [`ctModelPeriodicComponent()`](https://cdriveraus.github.io/ctsem/reference/ctModelPeriodicComponent.md)
  : Add a cycle, such as a daily rhythm, that drives a model's latent
  processes
- [`ctPeriodicPhases()`](https://cdriveraus.github.io/ctsem/reference/ctPeriodicPhases.md)
  : Add the phase predictors used by ctModelPeriodicComponent
- [`ctModelUniqueProcess()`](https://cdriveraus.github.io/ctsem/reference/ctModelUniqueProcess.md)
  : Give indicators processes of their own
- [`ctModelLatentLevel()`](https://cdriveraus.github.io/ctsem/reference/ctModelLatentLevel.md)
  : Put subjects' levels in the latent processes instead of the
  indicators
- [`ctModelConvertOMX()`](https://cdriveraus.github.io/ctsem/reference/ctModelConvertOMX.md)
  : Convert an old OpenMx-style ctsem model to the modern ctsem model
  format.
- [`ctModelIsNonlinear()`](https://cdriveraus.github.io/ctsem/reference/ctModelIsNonlinear.md)
  : Does a model or fit have context-dependent (nonlinear) matrix cells?
- [`ctIdentify()`](https://cdriveraus.github.io/ctsem/reference/ctIdentify.md)
  : Check which parameters a dataset can inform, before fitting
- [`print(`*`<ctStanModel>`*`)`](https://cdriveraus.github.io/ctsem/reference/print.ctStanModel.md)
  : Print a ctsem model
- [`plot(`*`<ctStanModel>`*`)`](https://cdriveraus.github.io/ctsem/reference/plot.ctStanModel.md)
  : Prior plotting

## Prepare data

- [`ctWideToLong()`](https://cdriveraus.github.io/ctsem/reference/ctWideToLong.md)
  : ctWideToLong Convert ctsem wide to long format
- [`ctLongToWide()`](https://cdriveraus.github.io/ctsem/reference/ctLongToWide.md)
  : ctLongToWide Restructures time series / panel data from long format
  to wide format for ctsem analysis
- [`ctWideNames()`](https://cdriveraus.github.io/ctsem/reference/ctWideNames.md)
  : ctWideNames sets default column names for wide ctsem datasets.
  Primarily intended for internal ctsem usage.
- [`ctIntervalise()`](https://cdriveraus.github.io/ctsem/reference/ctIntervalise.md)
  : Converts absolute times to intervals for wide format ctsem panel
  data
- [`ctDeintervalise()`](https://cdriveraus.github.io/ctsem/reference/ctDeintervalise.md)
  : ctDeintervalise
- [`ctDiscretiseData()`](https://cdriveraus.github.io/ctsem/reference/ctDiscretiseData.md)
  : Discretise long format continuous time (ctsem) data to specific
  timestep.

## Fit

- [`ctFit()`](https://cdriveraus.github.io/ctsem/reference/ctFit.md) :
  Fit a ctsem model
- [`ctOptimControl`](https://cdriveraus.github.io/ctsem/reference/ctOptimControl.md)
  : Optimiser settings for ctFit
- [`ctFitUpdate()`](https://cdriveraus.github.io/ctsem/reference/ctFitUpdate.md)
  : Update a fit to new data, or to the current version of ctsem
- [`ctFitUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
  [`ctOptimUncertainty()`](https://cdriveraus.github.io/ctsem/reference/ctFitUncertainty.md)
  : Compute or sample a fit's uncertainty
- [`ctEmpiricalBayesFit()`](https://cdriveraus.github.io/ctsem/reference/ctEmpiricalBayesFit.md)
  : Empirical Bayes subject-wise ctsem fits
- [`summary(`*`<ctEmpiricalBayesFit>`*`)`](https://cdriveraus.github.io/ctsem/reference/summary.ctEmpiricalBayesFit.md)
  : Summarise empirical Bayes subject-wise ctsem fits
- [`ctsem-options`](https://cdriveraus.github.io/ctsem/reference/ctsem-options.md)
  : Options ctsem reads

## Inspect a fit

- [`summary(`*`<ctJuliaFit>`*`)`](https://cdriveraus.github.io/ctsem/reference/summary.ctStanFit.md)
  [`summary(`*`<ctStanFit>`*`)`](https://cdriveraus.github.io/ctsem/reference/summary.ctStanFit.md)
  : Summarise a ctsem fit
- [`print(`*`<summary.ctStanFit>`*`)`](https://cdriveraus.github.io/ctsem/reference/print.summary.ctStanFit.md)
  : Print ctStanFit summaries
- [`plot(`*`<ctJuliaFit>`*`)`](https://cdriveraus.github.io/ctsem/reference/plot.ctJuliaFit.md)
  : Plots for ctJuliaFit objects
- [`plot(`*`<ctStanFit>`*`)`](https://cdriveraus.github.io/ctsem/reference/plot.ctStanFit.md)
  : plot.ctStanFit
- [`ctSummaryMatrices()`](https://cdriveraus.github.io/ctsem/reference/ctSummaryMatrices.md)
  : ctSummaryMatrices
- [`ctSubjectPars()`](https://cdriveraus.github.io/ctsem/reference/ctSubjectPars.md)
  : Extract an array of subject specific parameters from a fit
- [`ctExtract()`](https://cdriveraus.github.io/ctsem/reference/ctExtract.md)
  : Extract samples from a ctsem fit
- [`ctRawParnames()`](https://cdriveraus.github.io/ctsem/reference/ctRawParnames.md)
  : ctRawParnames
- [`ctPlotPosterior()`](https://cdriveraus.github.io/ctsem/reference/ctPlotPosterior.md)
  : Prior and posterior distributions of a fit's parameters
- [`ctTIpredEffects()`](https://cdriveraus.github.io/ctsem/reference/ctTIpredEffects.md)
  : Get time independent predictor effect estimates
- [`ctPredictTIP()`](https://cdriveraus.github.io/ctsem/reference/ctPredictTIP.md)
  : ctPredictTIP
- [`ctDiscretePars()`](https://cdriveraus.github.io/ctsem/reference/ctDiscretePars.md)
  : ctDiscretePars
- [`ctDiscreteParsPlot()`](https://cdriveraus.github.io/ctsem/reference/ctDiscreteParsPlot.md)
  : ctDiscreteParsPlot
- [`summary(`*`<ctDiscretePars>`*`)`](https://cdriveraus.github.io/ctsem/reference/summary.ctDiscretePars.md)
  : Summarise discrete time effect curves
- [`ctNetwork()`](https://cdriveraus.github.io/ctsem/reference/ctNetwork.md)
  : Temporal and contemporaneous networks of a continuous time model
- [`ctNetworkPlot()`](https://cdriveraus.github.io/ctsem/reference/ctNetworkPlot.md)
  : Draw the temporal and contemporaneous networks of a continuous time
  model
- [`ctPhasePortrait()`](https://cdriveraus.github.io/ctsem/reference/ctPhasePortrait.md)
  : Phase portrait of a fitted continuous time model
- [`ctStateDependencePlot()`](https://cdriveraus.github.io/ctsem/reference/ctStateDependencePlot.md)
  : How model matrix cells vary with a latent state
- [`ctContextDependence`](https://cdriveraus.github.io/ctsem/reference/ctContextDependence.md)
  : Context dependence of reported model matrices
- [`ctVarianceDecomposition()`](https://cdriveraus.github.io/ctsem/reference/ctVarianceDecomposition.md)
  : Decompose model implied variance into between person, process and
  measurement
- [`ctTracePlot()`](https://cdriveraus.github.io/ctsem/reference/ctTracePlot.md)
  : Plot how a fit reached its estimate

## Predict

- [`ctPredict()`](https://cdriveraus.github.io/ctsem/reference/ctPredict.md)
  : ctPredict
- [`plot(`*`<ctKalmanDF>`*`)`](https://cdriveraus.github.io/ctsem/reference/plot.ctKalmanDF.md)
  : Plots prediction output from ctPredict.
- [`ctKalmanArray()`](https://cdriveraus.github.io/ctsem/reference/ctKalmanArray.md)
  : Get Kalman filter estimates from a ctStanFit object
- [`ctResiduals()`](https://cdriveraus.github.io/ctsem/reference/ctResiduals.md)
  : Extract Standardized Residuals from a ctsem Fit

## Check, compare and improve

- [`ctReport()`](https://cdriveraus.github.io/ctsem/reference/ctReport.md)
  : Write a folder of summary, diagnostic and interpretation output for
  a fit
- [`ctScoreTest()`](https://cdriveraus.github.io/ctsem/reference/ctScoreTest.md)
  : Score tests for relaxing a fitted model
- [`plot(`*`<ctScoreTest>`*`)`](https://cdriveraus.github.io/ctsem/reference/plot.ctScoreTest.md)
  : Plot score tests
- [`ctFitAdequacy()`](https://cdriveraus.github.io/ctsem/reference/ctFitAdequacy.md)
  : Adequacy of a fitted model, as tests
- [`plot(`*`<ctFitAdequacy>`*`)`](https://cdriveraus.github.io/ctsem/reference/plot.ctFitAdequacy.md)
  : Plot the adequacy checks
- [`ctModelImprove()`](https://cdriveraus.github.io/ctsem/reference/ctModelImprove.md)
  : Improve a model step by step, from score tests and predictive checks
- [`plot(`*`<ctModelImprove>`*`)`](https://cdriveraus.github.io/ctsem/reference/plot.ctModelImprove.md)
  : Plot a model improvement search
- [`ctModelProposals()`](https://cdriveraus.github.io/ctsem/reference/ctModelProposals.md)
  : Changes the score tests and adequacy checks propose
- [`ctPostPredPlots()`](https://cdriveraus.github.io/ctsem/reference/ctPostPredPlots.md)
  : Posterior predictive checks for a ctsem fit
- [`ctPostPredict()`](https://cdriveraus.github.io/ctsem/reference/ctPostPredict.md)
  : Posterior predictive checks for a ctsem fit
- [`ctPostPredData()`](https://cdriveraus.github.io/ctsem/reference/ctPostPredData.md)
  : Create a data.table to compare data generated from a ctsem fit with
  the original data.
- [`ctFitCheckCov()`](https://cdriveraus.github.io/ctsem/reference/ctFitCheckCov.md)
  : Visual lagged covariance or correlation diagnostics for ctsem fits.
- [`ctFitCovCheckPlot()`](https://cdriveraus.github.io/ctsem/reference/ctFitCovCheckPlot.md)
  : ctFitCovCheckPlot
- [`ctACF()`](https://cdriveraus.github.io/ctsem/reference/ctACF.md) :
  Continuous Time Autocorrelation Function (ctACF)
- [`plotctACF()`](https://cdriveraus.github.io/ctsem/reference/plotctACF.md)
  : Plot an approximate continuous-time ACF object from ctACF
- [`ctACFresiduals()`](https://cdriveraus.github.io/ctsem/reference/ctACFresiduals.md)
  : Calculate Continuous Time Autocorrelation Function (ACF) for
  Standardized Residuals of ctsem fit.
- [`anova(`*`<ctJuliaFit>`*`)`](https://cdriveraus.github.io/ctsem/reference/anova.ctJuliaFit.md)
  [`anova(`*`<ctStanFit>`*`)`](https://cdriveraus.github.io/ctsem/reference/anova.ctJuliaFit.md)
  : Likelihood-ratio tests between nested ctsem fits
- [`ctChisqTest()`](https://cdriveraus.github.io/ctsem/reference/ctChisqTest.md)
  : Chi Square test wrapper for ctsem fit objects.
- [`ctLOO()`](https://cdriveraus.github.io/ctsem/reference/ctLOO.md) :
  Cross validation for ctsem fits: K fold by refitting, or PSIS for
  Laplace fits
- [`ctFitProfile()`](https://cdriveraus.github.io/ctsem/reference/ctFitProfile.md)
  : Profile likelihood for a ctsem julia fit
- [`plot(`*`<ctFitProfile>`*`)`](https://cdriveraus.github.io/ctsem/reference/plot.ctFitProfile.md)
  : Plot a profile likelihood
- [`print(`*`<ctFitProfile>`*`)`](https://cdriveraus.github.io/ctsem/reference/print.ctFitProfile.md)
  : Print a profile likelihood
- [`ctLaplaceCheck()`](https://cdriveraus.github.io/ctsem/reference/ctLaplaceCheck.md)
  : Measure, and optionally correct, the Laplace approximation error in
  a fit
- [`ctParticleLik()`](https://cdriveraus.github.io/ctsem/reference/ctParticleLik.md)
  : Check a julia-backend fit against a particle-filter likelihood
- [`ctCoverageCheck()`](https://cdriveraus.github.io/ctsem/reference/ctCoverageCheck.md)
  : Check interval coverage by simulation

## Generate data

- [`ctGenerate()`](https://cdriveraus.github.io/ctsem/reference/ctGenerate.md)
  : ctGenerate
- [`ctGenerateFromFit()`](https://cdriveraus.github.io/ctsem/reference/ctGenerateFromFit.md)
  : Generate data from a fitted model

## Julia engine

- [`ctJuliaInstall()`](https://cdriveraus.github.io/ctsem/reference/ctJuliaInstall.md)
  : Set up the Julia backend
- [`ctJuliaStatus()`](https://cdriveraus.github.io/ctsem/reference/ctJuliaStatus.md)
  : Report Julia backend availability
- [`ctJuliaSetup()`](https://cdriveraus.github.io/ctsem/reference/ctJuliaSetup.md)
  : Configure the Julia engine used by ctsem
- [`ctJuliaProcesses()`](https://cdriveraus.github.io/ctsem/reference/ctJuliaProcesses.md)
  : Julia processes started for ctsem
- [`ctJuliaKill()`](https://cdriveraus.github.io/ctsem/reference/ctJuliaKill.md)
  : Stop Julia processes started for ctsem
- [`ctJuliaEvaluate()`](https://cdriveraus.github.io/ctsem/reference/ctJuliaEvaluate.md)
  : Evaluate a julia model's likelihood and gradient at given parameters

## Datasets

- [`ctExample1`](https://cdriveraus.github.io/ctsem/reference/ctExample1.md)
  : ctExample1
- [`ctExample1TIpred`](https://cdriveraus.github.io/ctsem/reference/ctExample1TIpred.md)
  : ctExample1TIpred
- [`ctExample2`](https://cdriveraus.github.io/ctsem/reference/ctExample2.md)
  : ctExample2
- [`ctExample2level`](https://cdriveraus.github.io/ctsem/reference/ctExample2level.md)
  : ctExample2level
- [`ctExample3`](https://cdriveraus.github.io/ctsem/reference/ctExample3.md)
  : ctExample3
- [`ctExample4`](https://cdriveraus.github.io/ctsem/reference/ctExample4.md)
  : ctExample4
- [`longexample`](https://cdriveraus.github.io/ctsem/reference/longexample.md)
  : longexample
- [`Oscillating`](https://cdriveraus.github.io/ctsem/reference/Oscillating.md)
  : Oscillating
- [`AnomAuth`](https://cdriveraus.github.io/ctsem/reference/AnomAuth.md)
  : AnomAuth
- [`datastructure`](https://cdriveraus.github.io/ctsem/reference/datastructure.md)
  : datastructure
- [`ctstantestdat`](https://cdriveraus.github.io/ctsem/reference/ctstantestdat.md)
  : ctstantestdat
- [`ctstantestfit`](https://cdriveraus.github.io/ctsem/reference/ctstantestfit.md)
  : ctstantestfit

## Utilities

- [`ctsem`](https://cdriveraus.github.io/ctsem/reference/ctsem-package.md)
  [`ctsem-package`](https://cdriveraus.github.io/ctsem/reference/ctsem-package.md)
  : ctsem
- [`ctDocs()`](https://cdriveraus.github.io/ctsem/reference/ctDocs.md) :
  Get documentation pdf for ctsem
- [`ctCollapse()`](https://cdriveraus.github.io/ctsem/reference/ctCollapse.md)
  : ctCollapse Easily collapse an array margin using a specified
  function.
- [`ctPlotArray()`](https://cdriveraus.github.io/ctsem/reference/ctPlotArray.md)
  : Plots three dimensional y values for quantile plots
- [`ctPoly()`](https://cdriveraus.github.io/ctsem/reference/ctPoly.md) :
  Plots uncertainty bands with shading
- [`inv_logit()`](https://cdriveraus.github.io/ctsem/reference/inv_logit.md)
  : Inverse logit
- [`log1p_exp()`](https://cdriveraus.github.io/ctsem/reference/log1p_exp.md)
  : log1p_exp
- [`sdpcor2cov()`](https://cdriveraus.github.io/ctsem/reference/sdpcor2cov.md)
  : sdcor2cov
- [`test_isclose()`](https://cdriveraus.github.io/ctsem/reference/test_isclose.md)
  : Tests if 2 values are close to each other

## Stan internals

Lower level functions of the Stan backend.

- [`stanoptimis()`](https://cdriveraus.github.io/ctsem/reference/stanoptimis.md)
  : Optimize / importance sample a stan or ctStan model.
- [`stan_reinitsf()`](https://cdriveraus.github.io/ctsem/reference/stan_reinitsf.md)
  : Quickly initialise stanfit object from model and data
- [`stan_unconstrainsamples()`](https://cdriveraus.github.io/ctsem/reference/stan_unconstrainsamples.md)
  : Convert samples from a stanfit object to the unconstrained scale
- [`standatact_specificsubjects()`](https://cdriveraus.github.io/ctsem/reference/standatact_specificsubjects.md)
  : Adjust standata from ctsem to only use specific subjects
- [`stanWplot()`](https://cdriveraus.github.io/ctsem/reference/stanWplot.md)
  : Runs stan, and plots sampling information while sampling.

## Deprecated

- [`ctFitAddSamples()`](https://cdriveraus.github.io/ctsem/reference/ctFitAddSamples.md)
  : Append pseudo-posterior samples to an optimized ctStanFit object
  (deprecated)
- [`ctGenerateFromPriors()`](https://cdriveraus.github.io/ctsem/reference/ctGenerateFromPriors.md)
  : Generate data from a model's priors (deprecated)
