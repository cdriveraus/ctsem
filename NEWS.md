# ctsem News

## 27/8/2026
### 3.12.0 (development)
- Count manifest variables, on `backend='julia'`. Set `manifesttype = 3`; the observation is Poisson with a log link, so the latent process gives the log rate and no `ncategories` or thresholds are involved. The same adaptive Gauss-Hermite rule the binary and ordinal paths use integrates it against the predicted state, so the mode solve, the exact lift-back and the reverse-mode gradient all apply unchanged. Verified against closed forms (log likelihood, score and information exact) and against dense numerical integration (1e-13 for a tight prior, 4.7e-05 relative in the posterior variance at a predicted sd of 4). The rate is clamped at `exp(200)`: an infinite rate makes Newton's step `-Inf/Inf` and the mode solve return NaN rather than walk back, and on real data a rate of `exp(20)` is already half a billion events.
- Censored manifest variables, on `backend='julia'`. Set `manifesttype = 4` with `censormin` and/or `censormax`. Inside its limits the observation is the ordinary Gaussian density; at a limit it contributes the probability of being there. This is the one non-Gaussian type with a measurement error of its own, so unlike the others `ctFit` leaves its `MANIFESTVAR` free -- the others have their randomness from the link and would be adding noise on top of noise. `SpecialFunctions` becomes a declared dependency for `logerfc`, which nothing downloads that was not already present; the engine's own normal CDF is accurate to 7.5e-08 *absolute*, which is no relative accuracy once the tail probability is 1e-20, and an observation pinned at a floor while the model predicts well above it is exactly that.
- Fixed, and worth knowing about generally: ForwardDiff breaks comparison ties lexicographically on the partials, so `5.0 >= 5.0` is true while `5.0 >= Dual(5.0, [0,1,0])` is false. A censored observation sits exactly at its limit, and the limits are seeded alongside the standard deviation they share a vector with -- so during differentiation an at-limit observation took the interior branch and contributed a density where the forward pass had contributed a tail probability. MANIFESTVAR's gradient was 1% wrong with two censored observations and 263% wrong with 289, while every other parameter stayed exact to 1e-11. Censoring is now decided on primal values.
- The prior-warmed retry after a failed fit is removed. `carefulfit` warms every fit from the priors already and converged 240 of 240 in each condition measured, so the retry could only act where that had already succeeded -- while it did move answers: a fit deliberately capped at one iteration from supplied starting values came back from raw 24 at raw 12.8, reported without comment.
- Fixed, `intoverpop='laplace'` with a categorical measurement: an observation whose predicted variance had collapsed was skipped entirely by the reverse pass, not just by the forward one. The forward leaves the state alone there and adds `log P(y | eta)`, and that term depends on the linear predictor -- so on the state, LAMBDA, MANIFESTMEANS and the thresholds. Skipping it dropped all of that from the gradient, exactly where a variance is collapsing and the optimiser most needs the gradient to point somewhere. The reverse now contributes the degenerate term, and agrees with forward mode to 1e-9 on a model whose T0VAR is fixed at zero.
- `_laplace_finite` tests a dual's partials, not only its value. A seeded sweep exists for its derivatives, so accepting one whose derivatives are NaN because its value is finite defeats the check entirely -- and that is how a NaN reached the assembled gradient, where it was indistinguishable from a bad trial point. Each stage of the seeded assembly now checks what it produced, so the fallback to the nested route is taken at the point of failure rather than after the sum.
- The categorical update shares one floor on the predicted variance between its forward and reverse passes, set where `s²` cubed stops being representable rather than anywhere a real model reaches. Forward and reverse have to agree on where that boundary is, or the adjoint differentiates a function the forward never computed. A floor chosen for *accuracy* instead -- 1e-10, below which the exact shrink loses its digits -- was tried and is wrong: it fires in ordinary fits and costs 1.8 log units.
- Fixed, categorical measurement: the quadrature carries the observation's *log* likelihood rather than its probability. A category probability underflows to zero once the linear predictor is a few hundred away from the threshold bounding it, which an optimiser reaches while its parameters are still poor; the whole node, and usually the whole observation, then became impossible, the subject's likelihood `-Inf`, and the trial point invalid. Carried as a logarithm it is a large negative number and the optimiser can walk back out. On one 25-subject Laplace fit that difference was 38 of 78 trial points whose inner mode solve had nothing to work with, and a fit that gave up after one iteration 142 log units short. The node sums are accumulated relative to the largest exponent for the same reason they were before.
- Fixed, `intoverpop='laplace'`: a seeded gradient sweep that returned success could still have accumulated a non-finite contribution -- it reports whether its factorizations worked, not whether the numbers coming out of them are usable. One NaN there is the whole gradient, and the trial point was then rejected with a perfectly good objective value attached to it. The fallback to the nested route already existed for failed factorizations and now covers this too; `verbose` reports how often it was used, since a fit spending its time there is slow for a reason worth knowing.
- With those two, the 25-subject ordinal model that motivated all of this converges on **ten starting draws in ten**, with no trial point rejected for any reason -- and still ten in ten with the variance transform's `1e-10` floor deliberately removed, which is the point: the engine no longer needs the specification to keep it out of trouble. It was two in ten, with one `InexactError` and two fits reporting success 146 log units below the answer.
- Fixed, `backend='julia'` with `intoverpop='laplace'`: L-BFGS took an unscaled first step. It has no curvature history on its first iteration, so it goes downhill with whatever the initial step guess gives, and Optim's default is an alpha of one -- a first step as long as the gradient, which on this objective is routinely in the hundreds. `ctsem_optimize` was fixed for exactly this a while back; the Laplace route was left out of it. The symptom was extreme sensitivity to the starting draw, which is only `rnorm(npar, 0, 0.01)` and cannot itself explain anything.
- Fixed, both julia optimisers: convergence was reported on `Optim`'s own verdict, which is the disjunction of its x, f and g criteria -- and the first two are satisfied by a line search that stops making progress, since the step went to zero so x did not move and f did not change. That is the signature of giving up, not of arriving. One starting draw in ten of a 25-subject ordinal fit stopped after three iterations with a gradient of 681, a log likelihood 146 units below the optimum, and a verdict of success; another reported convergence at a NaN gradient. Convergence now needs a small *and finite* gradient.
- Fixed, both julia optimisers: the per-evaluation objective closure named its result `result`, and so did the enclosing function for its `Optim` result. Julia binds an assignment inside a closure to the enclosing scope's local of the same name, so every objective call was overwriting the outer one. Harmless while the outer result was only read after `Optim.optimize` returned, and not harmless the moment anything read it in between.
- Fixed: the matrix exponential threw `InexactError: Int64(NaN)` from six frames inside the reverse pass when handed a matrix with a non-finite entry, because its squaring count is `ceil(Int, log2(norm))`. A trial point that produces a non-finite DRIFT is an ordinary event -- it is what a line search is for -- and every caller already treats a non-finite objective as an invalid point. It now returns NaN and lets them.
- Fixed, categorical measurement: three separate places where a quantity was formed by subtracting two nearly equal numbers. The interval probability of an ordinal category is now a product of logistic tails rather than a difference of CDFs, and its score and information have closed forms with no division at all, so a threshold gap can close to 1e-15 with the derivatives still exact (at 1e-6 the previous form returned NaN). The posterior mean comes back as an *offset* from the prior mean and the variance is accumulated about the mode, because the filter never wants either quantity itself -- it wants `(mean - a)/s²` and `(1 - v/s²)/s²`, and computing those from numbers of order one destroys them: at a predicted variance of 1e-8 the covariance shrink came out as 1.024 where the answer is 0.439, and at 1e-10 as **-2215**, which makes the filter *widen* the covariance on an observation. And the quadrature accumulates relative to its largest weight, so `M1/Z` cannot divide by a denormal -- that produced a finite value with NaN partials at a predicted variance of 23, where nothing looks extreme.
- Fixed, categorical measurement: an observation whose predicted variance was exactly zero contributed a log likelihood of *zero* -- as if it were certain -- rather than `log P(y | eta)`. That made a vanishing variance pay: with a free T0VAR the optimiser was rewarded for driving it to zero, buying about seventy log units of nothing on twenty-five subjects with two indicators at the first occasion, and the objective was discontinuous at the point it was being driven towards. The correct value is also the continuous limit, so nothing has to know where the boundary is.
- Fixed, `intoverpop='laplace'`: the inner mode solve used an absolute tolerance, which asks for a number of digits that depends on how large the unit's objective happens to be -- the same mistake the outer `g_tol` makes. One unit stalled at `1.009e-10` against `1e-10`, missing by one percent, and because a unit that misses invalidates the *whole* trial point the outer line search lost 50 of its 89 evaluations to it. The tolerance now scales, with the absolute one kept underneath.
- Fixed, `backend='julia'`: variance transforms reached the engine without their `1e-10` floor. `ctModelTransformsToNum` recovers the four transform numbers from the model's transform string by a grid search scored on squared residual, so it cannot see a constant too small to move the residual, and the `round(x, 6)` that follows finishes it; DRIFT's `1e-06` floor survives, which is why only the variances lost theirs. The julia backend now reads the model's own transform text, recorded before that reduction. (The text is kept off `ctm$pars`, because parts of the Stan pipeline read that frame's columns by position -- carrying it there as a character column turned Stan's `pop_CINT` into NaN.)
- Measured together on a 25-subject, two-indicator ordinal model with a random CINT and a TI predictor, over ten starting draws: two draws converged before, nine after, eight of them to the same optimum the state-augmented fit finds. Before, one draw threw `InexactError`, and two reported convergence 146 log units below the answer. With the variance transform's floor deliberately removed, all ten still reach within 0.25 log units and the ones that stop short say so.
- Fixed, `backend='julia'` with binary indicators: the reverse-mode gradient was wrong whenever the model had more than one latent state. The binary measurement update conditions on `c = P lambda`, and the reverse of that adds a rank-one `cbar lambda'` to the covariance cotangent, which is not symmetric -- but the step before it treated the cotangent as symmetric and doubled one half instead of adding both. With a single latent state the cotangent is a scalar and the two agree, which is why every gradient test the binary path had passed. With two states, measured against forward mode: DRIFT and DIFFUSION about 1% out, and the second state's variance 12%. That is small enough to pass for quadrature error and quite large enough to stop an optimiser early, and it applied to any multi-latent binary model as well as to any binary model with a random effect, since `intoverpop='augmented'` expands the state by one per varying parameter. The adjoint now agrees with forward mode to 5e-15 there, and the regression is pinned with a two-state test in both the binary and ordinal suites.
- Ordinal manifest variables, on `backend='julia'`. Set `manifesttype = 2` and give `ctModel()` the `ncategories` of each ordinal variable; the observation is integrated against the predicted state under a cumulative logit, `P(y <= k | eta) = plogis(tau_k - eta)`, by the same adaptive Gauss-Hermite rule the binary path uses. Nothing about that path was Bernoulli-specific except the likelihood inside the node loop, and binary is the two-category case of the new one -- they agree bit for bit, which the tests pin. Measured against dense numerical integration over four categories, the 21-node rule is accurate to 5.3e-06 relative in the posterior variance at a predicted sd of 4, better than the binary case because an interior category is a bump rather than a step. The reverse-mode gradient covers the thresholds too, and agrees with finite differences to 4e-09.
- Ordinal thresholds live in a new `THRESHOLDS` model matrix, one row per manifest variable, and hold the first threshold followed by the *gap* to each subsequent one under a positive transform. Thresholds must increase, and `K-1` unconstrained cells get crossed by an optimizer -- at which point the category between the crossed pair has probability zero, the likelihood is `-Inf`, and no gradient points back out. Storing gaps makes the ordering hold by construction with no constraint to enforce. Variables with different category counts share the matrix; the unused cells are fixed and never read.
- `ctFit()` refuses an ordinal model on the stan backend rather than silently ignoring its thresholds, and checks ordinal data against the model: values must be consecutive integers from 1, must not exceed `ncategories`, and an empty category is warned about because the thresholds bounding it are not identified. `ctGenerate()`, `ctModelLatex()`, `print()`, `ctKalman()`, `ctPredict()`, `ctLOO()` and `summary()` all handle ordinal fits.
- Julia backend: manifest data crosses to the engine as `Matrix{Float64}` with `NaN` for missing, rather than whatever element type the data frame held. `as.matrix()` on integer columns with NAs gave `Matrix{Union{Missing,Int64}}`, on doubles `Matrix{Union{Missing,Float64}}`, and plain `Matrix{Float64}` only when nothing was missing -- so the element type of the user's data decided whether any precompiled code applied, and mostly it did not. First evaluation of a captured model shape through R: 12.8 s to 3.5 s. (The same shape in a pure Julia session is 0.017 s; the rest is a separate cause, recorded in the code.)

- Fixed, `backend='julia'` with `intoverpop='laplace'`: a fit could stop at its starting values and report convergence. The inner random-effect modes are kept between evaluations as warm starts, and an outer line search that visited a distant point could strand them somewhere Newton could not recover from -- after which the *same* parameter vector evaluated to a different number, the optimizer saw an objective that was not a function of its argument, and its first line search failed. Optim's overall verdict is a disjunction that a failed first line search satisfies trivially, so this was reported as success. Two fits in thirteen did it on a 40-subject simulation. The mode solve now retries from the origin whenever a warm start fails, a trial point whose inner solve did not converge is treated as invalid rather than accepted, and a fit that ends where it started with a non-zero gradient is reported as not converged and warns.
- `ctLaplaceCheck()` measures how much of an `intoverpop='laplace'` fit is the approximation, by recomputing the same random-effect integral with adaptive Gauss-Hermite quadrature. It reports the gap in log units and, from one Newton step against the fit's own Hessian, how far each estimate would move if the integral were taken exactly -- in standard errors, so it can be read directly. On a simulated model with a random DRIFT the population standard deviation's correction is about 0.9 standard errors, which is the whole of that parameter's under-coverage.
- Julia backend: `cores` is now a ceiling that the engine tunes within rather than a chunk count it obeys. The subject loop is not monotone in the chunk count -- the reverse sweep allocates thousands of small arrays per subject, and past a few million allocations a second the allocator is what the threads queue for -- so on a one-latent, one-indicator model 23 threads ran 3.7x *slower* than one, while the same code on twenty latents ran 2.7x faster. A handful of evaluations at fit start now measures the crossover on the model in front of it.
- Julia backend: parameter transform expressions are compiled once per distinct expression rather than once per model. Each `eval` of a transform string minted a new function type, which is a type parameter of the model object, so every model recompiled the whole filter -- 34 s on the first evaluation of a one-latent model, and another 29 s for the next model built from the same strings. Refitting the same model to new data, as cross-validation, bootstrapping and simulation studies do, is now free of that entirely.
- `ctLaplaceCheck()` handles nested groupings. A group's integral does not factor over its members, but it does factor conditionally -- given the group effect the members are independent -- so the rule recurses over the same block tree the fit already builds, re-solving each subject's conditional mode at every outer node and taking each outer block's scale from its curvature after the blocks beneath it have been eliminated. Cost is linear in the number of groups at each level and exponential only in that level's own number of random effects. With one level it is the flat rule it started as, which is what `nodes = 1` reproducing the Laplace value at one, two and three levels pins.
- Julia backend: the reverse pass no longer allocates its temporaries. The tape reuses its records between subjects rather than rebuilding them, and `_reverse_update!` and `_reverse_predict!` take every working matrix from a preallocated scratch. A 24-row subject of a one-latent model allocated 6,770 arrays per sweep and now allocates 1,673; a 184-subject evaluation went from 275 MB to 90 MB. This is a single-core saving and was *not* what stopped the threading -- cutting allocations 75% moved the parallel curve by nothing, which is how the LAPACK lock above was found instead.
- Julia backend: small dense linear algebra no longer goes to LAPACK or BLAS. Two separate traps, both process-global locks that `BLAS.set_num_threads(1)` does not touch. **Factorizations:** `potrf` on a 4x4 is almost entirely the lock OpenBLAS takes to acquire its scratch buffer -- measured at 23 threads, LAPACK runs at 0.11-0.30x of serial where an unblocked loop runs at 5.9-10.9x, and the loop is 1.7-11x faster on one core besides. **Products with a transposed view:** `mul!(C, transpose(view(A,...)), B)` cannot reach `gemm`, and runs at 0.09-0.11x threaded as well as slower serially; the buffered reverse pass was made almost entirely of that shape. The engine now has its own Cholesky, triangular solves and matrix products (`small_linalg.jl`), routes small linear systems through the LU it already had for dual numbers, and takes the Frechet exponential through `my_exp!`. Both thresholds hand the work back to LAPACK/BLAS above a cube of side sixteen, so wide models are unaffected -- measured neutral at twenty latents and 1.8x faster at one -- and `ctsem_set_small_linalg!` moves them.
- Julia backend: the subject loop now threads. A 184-subject, one-latent evaluation went from 0.86x on 23 threads to **6.67x**, and from 0.2608 s to 0.1268 s serial: **13.7x more throughput** than before this work. Wider models do better still (14.5x at twenty latents). `cores` is respected exactly as a ceiling, and `ctsem_tune_chunks!` picks the count within it.
- The engine test suite has a parallel runner, `test/runtests_parallel.jl`: one process per file, longest first, 22 minutes down to 9.5. It also removed a real fragility -- `test_quadrature.jl` reached fixtures defined in `test_laplace.jl` and worked only because `runtests.jl` includes files in sorted order into one module. Those fixtures now live in `laplace_fixtures.jl` and both suites include them.
- Julia backend: units are assigned to threads heaviest-first rather than in contiguous blocks, which matters once subjects are nested in groups of unequal size, and fits report `f_calls`, `g_calls` and `stalled` alongside `iterations`.

## 25/8/2026
### 3.12.0

- New `backend='julia'` option for `ctFit()`, using a Julia extended-Kalman-filter engine with a hand-written reverse-mode adjoint for the gradient. `ctJuliaInstall()` supplies whatever is missing -- the bridge package, Julia itself, and the engine's dependencies -- asking before it downloads anything. The engine ships inside ctsem, so no repository access is needed.
- The julia backend covers continuous and discrete time, individually varying parameters, time-independent and time-dependent predictors, state-dependent (nonlinear) model matrices, and priors. Threading splits the subject loop, giving 2-4x on six threads.
- `ctFit(backend='julia')` estimates uncertainty as part of fitting, as the stan backend does, and reads the same `optimcontrol` settings for it (`uncertainty`, `uncertaintyDraws`, `finishsamples`, `uncertaintyControl`, and `estonly` to skip). `summary()` therefore reports standard errors and intervals rather than point estimates alone.
- Uncertainty for julia fits uses an exact Hessian, obtained by differentiating the engine's reverse-mode gradient in forward mode rather than by finite differences. It agrees with the analytic result to machine precision, needs no step size, and costs `npar/chunksize` sweeps instead of `2*npar`.
- `summary()` for julia fits reports the same sections as an optimized stan fit, in the same order: standardised residual covariance, random-effects correlations, time-independent predictor effects (linearised onto the transformed parameters), system matrices, random-effects standard deviations, fixed effects, and log posterior.
- Julia fits carry their constrained draws on the fit object, as stan fits carry `transformedpars`, so summarising is a collapse over something already computed rather than a fresh pass through the transforms.
- `ctLOO()` works with julia fits, including `subjectwise`, `leaveOutN`, `keepfirstobs`, `refit=FALSE` and `casewiseApproximation`.
- `plot()`, `ctSubjectPars()`, `ctPredict()`, `ctKalman()`, `ctPredictTIP()`, `ctResiduals()`, `ctACFresiduals()`, `ctDiscretePars()`, `ctSummaryMatrices()`, `ctExtract()`, `ctGenerateFromFit()`, `ctPostPredPlots()` and `ctFitCovCheck()` all accept julia fits.
- The julia engine improves on three reported quantities: the measurement model is re-evaluated at the updated state before the filtered observation is recorded, the interval transition composes its substeps rather than recomputing an exponential from the last one, and the transition includes the time-dependent predictor impulse Jacobian.
- Fixed: population values and standard deviations for individually varying parameters implemented as state expansions were reported from the raw carrier state, without the parameter's own transform applied.

## 29/6/2026
### 3.11.0

- Function names have been simplified to remove Stan-specific wording from the main user-facing API. Use `ctFit()` rather than `ctStanFit()`, and the new non-`Stan` helper names such as `ctGenerateFromPriors()`, `ctGenerateFromFit()`, `ctPlotPosterior()`, `ctPostPredict()`, `ctSubjectPars()`, `ctTIpredEffects()`, and `ctSummaryMatrices()`. The old names remain as compatibility aliases for existing code. See the ctsem GitHub repository for further details: <https://github.com/cdriveraus/ctsem>.
- `ctFit()` now uses `model` as the model argument name. The old `ctstanmodel` argument is deprecated but still accepted for existing scripts.
- Modern fitting workflows should generally use `ctModel(type='ct')` or `ctModel(type='dt')` followed by `ctFit()`. `ctModel(type='omx')` objects are retained primarily for data generation and legacy workflows; convert them with `ctModelConvertOMX()` before fitting in the modern ctsem format.
- System matrices can be directly edited in the model object `$matrices` subobject, this will automatically update the data.frame containing full specification under `$pars` and vice-versa. 
- `ctEmpiricalBayesFit` is a new function that fits single subjects under the specified model, computes an empirical prior from those fits, then re-fits the subjects with that prior. 
- Optimized-fit uncertainty can be refreshed with `ctOptimUncertainty()`, including Hessian, local surrogate, one-step bootstrap, full refit bootstrap, sandwich, and outer product of gradient approximations.

## 27/6/2025
### 3.10.4

- Fixed bug in initial optimizer stages for nonlinear models requiring compilation.

## 18/6/2025
### 3.10.3
- Major improvements to importance sampling. 
- Include ctFitCovCheck function for checking empirical vs model implied covariance over time. 
- Reduce parallel compute init overhead
- Minor bug fixes to ctKalman
- refactor stanoptimis, changes to carefulfit logic, remove DEoptim option.

## 9/1/2025
### 3.10.2
- Stochastic optimizer improvements
- Bootstrap uncertainty improvements -- when fitting, use argument optimcontrol=list(uncertainty='bootstrap')

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
