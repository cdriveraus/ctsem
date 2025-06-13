# ctsem News

## 13/6/2025
### 3.10.3
- Include ctFitCovCheck function for checking empirical vs model implied covariance over time. 
- Reduce parallel compute init overhead
- Minor bug fixes to ctKalman
- refactor stanoptimis, changes to carefulfit logic, remove DEoptim option.

## 9/1/2025
### 3.10.2
- Stochastic optimizer improvements
- Bootstrap uncertainty improvements -- when fitting, use argument optimcontrol=list(bootstrapUncertainty=TRUE)

## 12/8/2024
### 3.10.1
- Fix bug introduced in 3.10.0 where certain combinations of Gaussian and binary variables cause convergence difficulties and invalid results. 
- Revert unconstrained correlation change introduced in 3.10.0, it was more difficult to fit in some cases.
- Detect duplicated T0MEANS parameters and propose alternative approach.
- Add ctPredictTIP function for examining and plotting differences due to time independent predictors.


## 10/05/2024
### 3.10.0
- Fix bug when computing Jacobian of certain nonlinear models.
- Modify unconstrained correlation approach for better optimization / uncertainty quantification and clearer interpretation.
- Add ctACF function for plotting approximate continuous time auto and cross correlations.
- Bug fixes to ctKalman plots, were occasionally confused re subject ID.
- Include experimental / imperfect approach to ordinal data. 

## 30/10/2023
### 3.9.1
- Fix bug in ctKalman - was dropping certain subjects resulting in no plots / output.
- Fix fatal (i.e, erroring out) bug in certain nonlinear parameter specifications.
- Allow direct references to time dependent predictors in nonlinear specifications - now measurement model can easily depend on time varying covariates.
- Update array syntax internally to rstan 2.26+ approach -- completely this time...

## 14/9/2023
### 3.9.0
- Add progress reports for stochastic optimizer and Hessian.
- Add small noise to improve sampling performance when `inits='optimize'`.
- Include ctACF and ctACFresiduafunction for approximate continuous time auto-correlations.
- Update array syntax internally to rstan 2.26+ approach.
- Improved stochastic optimizer.

## 20/8/2023
### 3.8.1
- Correct bug in nonlinear formulations when the same state is referenced for multiple nonlinearities.
- Correct unnecessary memory usage when computing Hessian with multiple cores.
- Improve stochastic subsampling first pass optimizer.
- Simplify discrete time model computations internally.
- Performance gains and reduced memory usage via usage of matrix exponential subsets and automatic computation of dynamic error indices.

## 20/6/2023
### 3.7.6
- Deprecate `nopriors` argument to `ctStanFit`, allow `priors` argument
- Allow `inits='optimize'` argument to `ctStanFit`, to speed up sampling approach
- Allow integer values for `removeObs` argument to `ctKalman`, for N step ahead predictions
- Allow `sameInitialTimes` argument to `ctStanFit`, to generate empty observations at earliest observation time, ensuring comparability of times at T0MEANS

## 24/3/2023
### 3.7.6
- Fix compile error for some higher dimensional non-linear models
- Fix NaN gradient error for certain non-linear measurement error models
- Use Rstantools for compile specification to ensure future compatibility
- Allow standardized error output from `ctKalman`

## 1/7/2022
### 3.7.0
- Fix some plotting features, speed up random effects a little

## 9/3/2022
### 3.6.0
- Fix: Random effects standard deviations were mis-estimated for time-dependent predictor effects and diffusion parameters when fitting with optimization

## 6/12/2021
### 3.5.5
- Some edge case optimizer problems resolved
- Bug fix in discrete time plots when `observational=TRUE`. Correlations were unnecessarily squared previously

## 22/7/2021
### 3.5.4
- Improved automatic imputation of time-independent predictors, fixed bug where too many imputed values were set to zero. Care still recommended if relying on automatic imputation though!

## 16/6/2021
### 3.5.3
- Fixed bug in output of `$rawpopcorr` introduced in 3.4.3. Correlations displayed incorrectly, other parameters unaffected
- Modified correlation approach to ensure monotonicity in high dims

## 31/5/2021
### 3.5.0
- Added `ctFitMultiModel` function to simplify processing of multiple models
- Added `ctChisqTest` function for simplified likelihood ratio tests between models
- Added subsampling optimization for first pass, faster for larger models/data

## 21/4/2021
### 3.4.3
- Changed optimization scheme, first BFGS, then stochastic gradient descent
- Fixed `ctStanTIpredEffects` function, much faster
- Altered correlation matrix approach - better optimization/sampling behavior, priors can differ by index though

## 10/2/2021
### 3.4.2
- `ctStanDiscretePars` temporal dependence plots work for discrete time also
- Fix: In certain circumstances with covariate effects on duplicated parameters, the effects may not have been completely applied in the last few releases
- Other small bug fixes and efficiency improvements

## 4/12/2020
### 3.4
- `ctLOO` function for leave one out/k-fold cross-validation
- `ctCheckFit` function dramatically improved for visual model diagnostics
- Fixes to a range of edge cases when specifying more complex nonlinearities
- `ctStanDiscretePars` has improved plotting options
- `ctStanFitUpdate` function can be used to attempt to update a saved `ctStanFit` object to the current version of ctsem
- `ctSaturatedCov` function for estimating a form of saturated model as a reference. Still a bit developmental.
- General robustness / efficiency improvements.


## 10/7/2020
### 3.3.8
- Stationarity,subject specific parameter output, and linear system estimation removed (uses nonlinear in all cases, a bit slower) to try satisfy CRAN compile time checks.

## 20/6/2020
### 3.3.2
- ctLOO function to compute leave k out entropy estimates for model comparison / validation.
- Optimizer parallelisation for single subject models
- ctStanFitUpdate function to use saved fit objects created in earlier versions of ctsem.
- Multiple core memory usage reductions.

## 26/4/2020
###3.2.1
- Minor updates to suit rstan 2.2.3
- Higher dim Hessian really really fixed...

## 18/4/2020
### 3.2.0
- Change to defaults of ctStanFit -- optimization without priors is new default.
- Optimization with binary variables much improved -- specify using: 
  mymodel$manifesttype <- c(1,0,1) 
  for 3 manifest variables, 1st and last binary and second continuous.
- Improved latex output for models and fits using ctModelLatex.
- Hessian estimation in higher dimensions fixed, again.
- Experimental automatic covariate (tipred) detection -- set{ mymodel$TIpredAuto <- 1L } to try it.
- Minor plotting / other fixes

## 10/2/2020
### 3.1.1
- Reverted to non sequential measurement update introduced in 3.1.0
- Further hessian estimation fixes
- Stochastic optimizer further improved


## 21/01/2020
### 3.1.0
- Fixed Hessian estimation sensitivity introduced last release.
- ctStanKalman has option to return subject specific estimates.
- Various performance improvements -- stochastic optimizer very effective with many parameters.
- Fix rounding, off by one issue when using timestep argument for nonlinear dynamics.



## 10/12/2019
### 3.0.9
- ctStanGenerate function -- generate from a ctstanmodel and prior distribution.
- Parallel optimization improvements -- memory usage halved, more cores nearly always useful.
- Missing time independent predictors single imputed when optimizing.
- ctModelHigherOrder function to easily add higher order structure to specified model. E.g. slow changing trends / oscillations.
- various minor output / efficiency / estimation robustness improvements.

## 30/10/2019
### 3.0.8
- ctStanFit bug fix: MANIFESTVAR wrongly reported as the sqrt
- ctFit bug fix: Std errors of covariances were too wide when transformedParams=TRUE
- TI predictor effects can be specified as the 5th element of parameter string in ctModel -- "drift11 | -exp(param) | FALSE | 1 | age, gender"
- Nonlinear models: Parameter names in the additional PARS matrix, and latent variables, can be referenced directly in parameters -- "-exp(cognition * drift11)" instead of "-exp(state[2] * PARS[1,1]". PARS matrix must still be specified.
- Switched plotting to ggplot2, many changes / improvements.

## 11/9/2019
### 3.0.4
- Removed some spurious warnings generated by last release.
- Improved ctKalman function for ctStanFit results (timestep now works).
- Parallelise optimization over subjects with ctStanFit (cores=xx).
- Bug fix: time varying diffusion generated errors.
- Use stochastic optimizer to check for improvements by default
- Simplify ctModel specification -- vectors interpreted as rowwise matrices, automatic dimension detection.

## 20/8/2019
### 3.0.1
- fixed a few minor / error generating bugs introduced in previous release related to random effects optimization.

## 28/7/2019
### 3.0.0
- parameter transformations can be specified naturally without recompilation
- analytic Jacobian's used for extended Kalman filter where possible
- improved optimizer performance
- generally improved performance, particularly for nonlinear systems.

## 14/5/2019
### 2.9.5
- corrected optimization of random effects using ctStanFit
- ctModelLatex function to display within subject model equation
- custom calculations allowed in ctModel for linear and nonlinear approaches
- Nonlinearity possible for discrete time now also
- Use stochastic optimizer by default for ctStanFit -- more robust

## 12/4/2019
### 2.9.0
- Improved optimizer using stochastic gradient descent.
- fixed bug introduced re finding start values when binary variables are used.
- custom calculations can be specified and estimated using both linear and nonlinear dynamics approach.
- ctStanFit is no longer available for win32 systems.
- ctStanKalman function extracts system state over time.

## 6/11/2018
### 2.7.3
- Updated for rstan 2.18.1 compatibility
- Non-linear dynamics now handled using mixture of extended and unscented filters for improved speed.
- Priors for hierarchical variance modified so prior for total variance has consistent shape regardless of dimension.
- Optimization / importance sampling works well for many cases, see arguments using ?ctStanFit.
- ctStanPostPredict produces a range of posterior predictive plots.
- Various small non-critical bug fixes / updates. (see github for details)

## 25/6/2018
### 2.6.5
- Fixed bug in ctStanFit introduced in previous release leading to errors in handling of missing data.
- ctStanTIpredMarginal function for plotting marginal relationships between predictors and parameters.

## 1/6/2018
### 2.6.0
- Removed need for compilation of standard models.
- Unscented Kalman filter for ctStanFit:
  - Non-linear / time-varying / state dependent specifications now possible. 
  - Optimization followed by importance sampling can be used instead of sampling via Stan.
  - Most plotting functions still not working correctly for such models.
- ctCheckFit function for plotting covariance of data generated from posterior against original.
- Time independent predictors can now be used independent of random effects.
- Fix bug in summary preventing display when binary variables were used in fit.
- Allow data sets to contain both binary and continuous variables.
- More robust data import, character string id's and jumbled order of rows now manageable.
- stanWplot no longer requires shiny to be explicitly loaded.
- Fix bug in ctKalman plotting function preventing interpolation.
- Fix bug in additional summary matrices introduced in 2.5.0 -- some were transposed.
- Summary no longer returns errors when partial stationarity is set.

## 27/9/2017
### 2.5.0
Fixes: 
- stanWplot function for trace plots while sampling with stan was not working on non windows platforms.
- ctStan summary reports population standard deviations more accurately -- improved delta approach.
- various minor plotting improvements

Additions / Changes:
- ctStanFit now handles correlation matrices differently -- little substantive impact.
- ctKalman can now be used to plot individual trajectories from ctFit objects and ctStanFit objects (ctStanKalman function no longer exists).
- ctStanFit now handles missing data on covariate effects -- time dependent predictors are set to zero, time independent predictors are imputed with a normal(0,10) prior (can adjust via the $tipredsimputedprior subobject of the ctStanModel).
- ctStanFit default population standard deviation prior now changed to a regularised independence Jeffreys -- previous truncated normal approach still possible because...
- ctStanFit now accepts custom specifications for the population standard deviation -- see the $rawhypersd , $rawhypersdlowerbound, and $hypersdtransform subobjects of the ctStanModel object. 
- ctStan: Plotting covariate effects via ctStanTIpredeffects function now easier to use and more versatile -- can plot effects on discrete time matrices, for instance.
- ctFit and ctMultigroupFit data argument changed to 'dat' instead of 'datawide', and now dataform="long" argument can be used to use long format data (as per ctStan) directly. 
- additional parameter matrices shown for summary of ctStanFit objects.
- ctStanParMatrices function to compute continuous time matrices for a given model and vector of free population means.


## 16/5/2017
### 2.4.0
Fixes:
- Time dependent predictors generated errors with the frequentist Kalman filter form since 2.2.0
- With stationary set to NULL (not the default but offered in help file) for ctFit, 
t0 matrices were mistakenly set to stationary.
- Duplicated parameter names now allowed in a ctStanFit model.

Features:
- ctGenerateFromFit generates data based on a model fitted with ctFit.
- ctPostPredict generates distributions from data based on a model fitted with ctFit
and plots this against the original data.


## 6/4/2017
### 2.3.1
Fixes:
- summary: Standard errors were not reported in some cases
- ctStanFit: 2.3.0 hierarchical correlation changes were applied too broadly
- ctFit: discreteTime switch no longer gives errors when traits included
- ctFit: transformedParams=FALSE argument no longer throwing errors.
- ctStanKalman: correct handling of missing data for plotting.

## 3/3/2017
### 2.3.0
Fixes:
- TRAITVAR in frequentist ctsem was incorrectly accounting for differing time 
  intervals since v2.0.0. TRAITVAR is now (again) reported as total between subjects
  variance.
- Default quantiles on ctStanDiscretePars adjusted to 95%.
- Hierarchical correlation probabilities adjusted in ctStanFit for more consistent
  behaviour with high dimensional processes.

Changes:
- Default to unstandardised cross effects plots.

## 1/2/2017
### 2.2.0
Changes:
- Time dependent predictors now have instantaneous effect in both frequentist and 
  Bayesian approaches, and the documentation is updated to reflect this.
  Previously, no TDpreds affecting first time point in frequentist.
  Accordingly, wide data structure is changed, with an extra column 
  per predictor and predictors now sorted by time point as for indicators. 
  See vignette for example. 
- Default to 0 covariance between time dependent predictors and initial (T0) 
  latents / traits / time independent predictors. Specify matrix as 'free' 
  in ctModel to estimate instead.
- Default carefulFit = TRUE for multiple groups frequentist models (ctMultigroupFit)
- Improve optimization approach for ctStanFit - but still not reliable for random effects.

Fixes:
- Multiple time dependent predictors with multiple processes resulted in inaccurate
  estimates for TDPREDEFFECT in frequentist approach of previous versions.
- Prevent ctGenerate from auto-filling matrices to 0 variance.
- Correct oscillating example for change in tolerance in OpenMx.


## 6/1/2017
### 2.1.1
Improvements:
- improved fitting of frequentist models with ctFit and ctRefineTo, due to
  changes to carefulFit penalisation and refining approach.

Changes:
- Removed package 'PSM' from suggests field and vignette as requested by CRAN

Fixes:
- rstan 2.14 caused problems with data import for ctStanFit
- eliminated spurious warnings for ctStanFit


## 20/12/2016
### 2.1.0
Features:
- Empirical Bayes, experimental but can now optimize with hierarchical model 
  (when using the Kalman filter, as per defaults)
- Easy extraction and plotting of time independent predictor (covariate) effects,
  see ctStanTIpredEffects for example.
- Added stationary argument to ctStanFit - much more efficient than setting 
  priors on stationarity.

Bugs fixed:
- incorrect number of cores spawned for parallel sessions.
- optimize and variational bayes switches for ctStanFit did not work.
- ctKalman would break if only 1 row of data passed in.


## 18/11/2016
### 2.0.0
Features:
- Hierarchical Bayesian modeling using Stan, see ctStanFit function and 
  the vignette at https://cran.r-project.org/package=ctsem/vignettes/hierarchical.pdf

Changes:
- Defaults change: Fix CINT to 0 and free MANIFESTMEANS
- Reintroduce variable effect of TRAITVAR at T0 (more flexible but more 
    fitting problems - try MANIFESTTRAITVAR instead if problematic, or 
    use step-wise fitting approach, automated with ctRefineTo)


### ctsem 1.1.6
Features added:
- now with a change log!
- ctCompareExpectation plots expected means and covariances against model implied.
- remove log transform of drift matrix diagonal, positive drift diagonals again possible.
- ctRefineTo allows easy step wise fitting from simple to complex - faster and more robust fitting in many cases.
- ctPlot is a new function that allows more customization of plots.
- ctModel now allows time varying means to be specified.

Bugs fixed:
- corrected handling of Cholesky inputs for ctGenerate
