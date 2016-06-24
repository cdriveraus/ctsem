<!-- README.md is generated from README.Rmd. Please edit that file -->
### Overview

ctsem allows for easy specification and fitting of a range of continuous and discrete time dynamic models, including multiple indicators (dynamic factor analysis), multiple, potentially higher order processes, and time dependent (varying) and independent (not varying) covariates. Classic longitudinal models like latent growth curves and latent change score models are also possible. Version 1 of ctsem provided SEM based functionality by linking to the OpenMx software, allowing mixed effects models (random means but fixed regression and variance parameters) for multiple subjects. For version 2 of the package , we include a Bayesian specification and fitting routine that uses the probabilistic programming language, via the package in R. This allows for all parameters of the dynamic model to individually vary, using an estimated population mean and variance, and any time independent covariate effects, as a prior. ctsem version 1 is documented in a forthcoming JSS publication (Driver, Voelkle, Oud, in press), and in R vignette form at <https://cran.r-project.org/web/packages/ctsem/vignettes/ctsem.pdf> , here we provide the basics for getting started with the new Bayesian approach.

### Subject Level Latent Dynamic model

This section describes the fundamental subject level model, and where appropriate, the name of the ctModel argument used to specify specific matrices.

The subject level dynamics are described by the following stochastic differential equation:
\begin{equation}
\label{eq:process1}
\mathrm{d}\eta(t) =
\bigg( 
A\eta(t) +
b +
M \chi(t)  
\bigg) \mathrm{d}t +
G \mathrm{d} W(t)  
\end{equation}
Vector \(\eta(t)\in\mathbb{R}^{v}\) represents the state of the latent processes at time \(t\). The matrix \(A \in \mathbb{R}^{v \times v}\) represents the DRIFT matrix, with auto effects on the diagonal and cross effects on the off-diagonals characterizing the temporal relationships of the processes.

The long term level of processes \(\eta(t)\) is determined by the continuous time intercept (CINT) vector \(b \in\mathbb{R}^{v}\), which (in combination with \(A\)) determines the long-term level at which the processes fluctuate around.

Time dependent predictors \(\chi(t)\) represent inputs to the system that vary over time and are independent of fluctuations in the system. The above equation shows a generalised form for time dependent predictors, that could be treated a variety of ways dependent on the assumed time course (or shape) of time dependent predictors. We use a simple impulse form, in which the predictors are treated as impacting the processes only at the instant of an observation occasion \(u\). When necessary, the evolution over time can be modeled by extending the state matrices.

\begin{equation}
\label{eq:spike}
\chi (t) = \sum_{ u \in U}  x_{u} \delta (t-t)     
\end{equation}
Here, time dependent predictors \(x_u \in \mathbb{R}^{l}\) are observed at times $ u U $. The Dirac delta function \(\delta(t-t)\) is a generalized function that is \(\infty\) at 0 and 0 elsewhere, yet has an integral of 1 (when 0 is in the range of integration). It is useful to model an impulse to a system, and here is scaled by the vector of time dependent predictors \(x_u\). The effect of these impulses on processes \(\eta(t)\) is then \(M\in \mathbb{R}^{v \times l}\) (TDPREDEFFECT).

\(W(t) \in \mathbb{R}^{v}\) represents independent Wiener processes, with a Wiener process being a random-walk in continuous time. \(dW(t)\) is meaningful in the context of stochastic differential equations, and represents the stochastic error term, an infinitesimally small increment of the Wiener process. Lower triangular matrix \(G \in \mathbb{R}^{v \times v}\) represents the effect of this noise on the change in \(\eta(t)\). \(Q\), where \(Q = GG^\top\), represents the variance-covariance matrix of the diffusion process in continuous time. The DIFFUSION matrix in ctModel is essentially G, except off-diagonal elements of DIFFUSION specify Cholesky decomposed correlation values, rather than covariance.

### Subject level measurement model

The latent process vector \(\eta(t)\) has measurement model:

\begin{equation}
\label{eq:measurement}
y(t) = \Lambda \eta(t) + h + \zeta(t)  
\quad \text{where } \zeta(t) \sim  \mathrm{N} (0, \Theta)
\end{equation}
\(y (t)\in\mathbb{R}^{c}\) is the manifest variables, \(\Lambda \in \mathbb{R}^{c \times v}\) is factor loadings (LAMBDA), \(h \in\mathbb{R}^{c}\) is the manifest intercepts (MANIFESTMEANS), and residual vector \(\zeta \in \mathbb{R}^{c}\) has covariance matrix \(\Theta \in\mathbb{R}^{c \times c}\). To specify \(\Theta\) with ctModel, the lower-triangular MANIFESTVAR matrix is used, with standard deviations on the diagonal and Cholesky decomposed correlations in the off-diagonals.

### Install software and prepare data

Install ctsem software from github repository <https://github.com/cdriveraus/ctsem> .

``` install
require('devtools')
install_github("ctsem",username='cdriveraus')
```

Prepare data in long format, each row containing one time point of data for one subject. We need a subject id column (named by default "id"), columns for manifest variables (the names of which must be given in the next step using ctModel), columns for time dependent predictors (these vary over time but have no model estimatedand are assumed to impact latent processes instantly - generally intervention or event dummy variables), and columns for time independent predictors (the value will be stable for each measurement of a particular subject). Relationships are estimated between time independent predictors and individually varying subject level parameters.

### Model specification

Specify model using . "stanct" specifies a continuous time model in Stan format, "standt" specifies discrete time, while "omx" is the classic behaviour and prepares an model. Other arguments to ctModel proceed as normal, although many matrices are not relevant for the Stan formats, either because the between subject matrices have been removed, or because time dependent and independent predictors (covariates that either change over time or don't) are now treated as fixed regressors and only require effect (or design) matrices.

``` model
checkm<-ctModel(type='stanct',
  n.latent=2, latentNames=c('eta1','eta2'),
  n.manifest=2, manifestNames=c('Y1','Y2'),
  n.TDpred=1, TDpredNames='TD1', 
  n.TIpred=3, TIpredNames=c('TI1','TI2','TI3'),
  LAMBDA=diag(2))
```

This generates a simple first order bivariate latent process model, with each process measured by a potentially noisy manifest variable. Additional complexity or restrictions may be added, the table below shows the basic arguments one may consider and their link to the dynamic model parameters. For more details see the ctsem help files or papers. Note that for the Stan implementation, ctModel requires variance covariance matrices (DIFFUSION, T0VAR, MANIFESTVAR) to be specified with standard deviations on the diagonal, and Cholesky decomposed correlations on the off diagonal. This is for computational reasons, and hopefully poses little concern for the user since in our experience these are most often set to be either free, or fixed to 0, which translates directly.

<table style="width:108%;">
<colgroup>
<col width="5%" />
<col width="5%" />
<col width="15%" />
<col width="81%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Argument</th>
<th align="left">Sign</th>
<th align="left">Default</th>
<th align="left">Meaning</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">n.manifest</td>
<td align="left"></td>
<td align="left"></td>
<td align="left">Number of manifest indicators per individual at each measurement occasion.</td>
</tr>
<tr class="even">
<td align="left">n.latent</td>
<td align="left"></td>
<td align="left"></td>
<td align="left">Number of latent processes.</td>
</tr>
<tr class="odd">
<td align="left">LAMBDA</td>
<td align="left"><span class="math inline"><strong>Λ</strong></span></td>
<td align="left"></td>
<td align="left">n.manifest <span class="math inline">×</span> n.latent loading matrix relating latent to manifest variables.</td>
</tr>
<tr class="even">
<td align="left">manifestNames</td>
<td align="left"></td>
<td align="left">Y1, Y2, etc</td>
<td align="left">n.manifest length character vector of manifest names.</td>
</tr>
<tr class="odd">
<td align="left">latentNames</td>
<td align="left"></td>
<td align="left">eta1, eta2, etc</td>
<td align="left">n.latent length character vector of latent names.</td>
</tr>
<tr class="even">
<td align="left">T0VAR</td>
<td align="left"></td>
<td align="left">free</td>
<td align="left">lower tri n.latent <span class="math inline">×</span> n.latent matrix of latent process initial cov.</td>
</tr>
<tr class="odd">
<td align="left">T0MEANS</td>
<td align="left"></td>
<td align="left">free</td>
<td align="left">n.latent <span class="math inline">×</span> 1 matrix of latent process means at first time point, T0.</td>
</tr>
<tr class="even">
<td align="left">MANIFESTMEANS</td>
<td align="left"><span class="math inline"><strong>τ</strong></span></td>
<td align="left">free</td>
<td align="left">n.manifest <span class="math inline">×</span> 1 matrix of manifest means.</td>
</tr>
<tr class="odd">
<td align="left">MANIFESTVAR</td>
<td align="left"><span class="math inline"><strong>Θ</strong></span></td>
<td align="left">free diag</td>
<td align="left">lower triangular matrix of var / cov between manifests</td>
</tr>
<tr class="even">
<td align="left">DRIFT</td>
<td align="left"></td>
<td align="left">free</td>
<td align="left">n.latent <span class="math inline">×</span> n.latent matrix of continuous auto and cross effects.</td>
</tr>
<tr class="odd">
<td align="left">CINT</td>
<td align="left"><span class="math inline"><strong>κ</strong></span></td>
<td align="left">0</td>
<td align="left">n.latent <span class="math inline">×</span> 1 matrix of continuous intercepts.</td>
</tr>
<tr class="even">
<td align="left">DIFFUSION</td>
<td align="left"><span class="math inline"><strong>Q</strong></span></td>
<td align="left">free</td>
<td align="left">lower triangular n.latent <span class="math inline">×</span> n.latent matrix of diffusion variance / covariance.</td>
</tr>
<tr class="odd">
<td align="left">n.TDpred</td>
<td align="left"></td>
<td align="left">0</td>
<td align="left">Number of time dependent predictors in the dataset.</td>
</tr>
<tr class="even">
<td align="left">TDpredNames</td>
<td align="left"></td>
<td align="left">TD1, TD2, etc</td>
<td align="left">n.TDpred length character vector of time dependent predictor names.</td>
</tr>
<tr class="odd">
<td align="left">TDPREDEFFECT</td>
<td align="left"><span class="math inline"><strong>M</strong></span></td>
<td align="left">free</td>
<td align="left">n.latent <span class="math inline">×</span> n.TDpred matrix of effects from time dependent predictors to latent processes.</td>
</tr>
<tr class="even">
<td align="left">n.TIpred</td>
<td align="left"></td>
<td align="left">0</td>
<td align="left">Number of time independent predictors.</td>
</tr>
<tr class="odd">
<td align="left">TIpredNames</td>
<td align="left"></td>
<td align="left">TI1, TI2, etc</td>
<td align="left">n.TIpred length character vector of time independent predictor names.</td>
</tr>
<tr class="even">
<td align="left">TIPREDEFFECT</td>
<td align="left"><span class="math inline"><strong>B</strong></span></td>
<td align="left">free</td>
<td align="left">n.latent <span class="math inline">×</span> n.TIpred effect matrix of time independent predictors on latent processes.</td>
</tr>
</tbody>
</table>

These matrices may all be specified using a combination of character strings to name free parameters, or numeric values to represent fixed parameters.

The parameters subobject of the created model object shows the parameter specification that will go into Stan, including both fixed and free parameters, whether the parameters vary across individuals, how the parameter is transformed from a standard normal distribution (thus setting both priors and bounds), and whether that parameter is regressed on the time independent predictors.

``` model
head(checkm$parameters,8)
```

One may modify the output model to either restrict between subject differences (set some parameters to fixed over subjects), alter the transformation used to determine the prior / bounds, or restrict which effects of time independent predictors to estimate. Plotting the original prior, making a change, and plotting the resulting prior, are shown here -- in this case we believe the latent process error for our first latent process, captured by row 1 and column 1 of the DIFFUSION matrix, to be very small, so restrict our prior accordingly to both speed and improve sampling.

``` transform
par(mfrow=c(1,2))
ctStanPlotPriors(checkm,rows=11)
checkm$parameters$transform[11]<- 'log(exp((param)*1.5)+1)*2'
ctStanPlotPriors(checkm,rows=11)
```

The plots show the prior distribution for the population mean of DIFFUSION\[1,1\] in black, as well as two possible priors for the subject level parameters. The blue prior results from assuming the population mean is two standard deviations lower than the mean for our prior, and assuming that the population standard deviation is 1, which given our prior on population standard deviations is a truncated normal(0, 0.5) distribution, is also two sd's from the base of 0. To understand better, the pre-transformation population sd prior for all subject level parameters looks like:

``` sdprior
sd<-rnorm(5000000,0,.5)
sd<-sd[sd>0]
plot(density(sd,bw=.01,n=50000),lwd=2)
```

Restrict between subject effects as desired. Unnecessary between subject effects will slow sampling, but be aware of the many parameter dependencies in these models -- restricting one parameter may sometimes lead to variation from it showing up elsewhere.

``` restrictbetween
checkm$parameters[25:28,]
checkm$parameters[25:28,]$indvarying<-FALSE
```

Also restrict time independent predictor effects in a similar way, for similar reasons. In this case, the only adverse effects of restriction are that the relationship between the predictor and variables will not be estimated, but the subject level parameters themselves should not be very different, as they are still freely estimated. Note that such effects are only estimated for individually varying parameters anyway -- so after the above change there is no need to set the tipredeffect to FALSE for the T0VAR variables, it is assumed. Instead, we restrict the tipredeffects on all parameters, and free them only for the auto effect of the first latent process.

``` restricttipred
checkm$parameters[,c('TI1_effect','TI2_effect','TI3_effect')]<-FALSE
checkm$parameters[7,c('TI1_effect','TI2_effect','TI3_effect')]<-TRUE
```

### Model fitting

Once model specification is complete, the model is fit to the data using the ctStanFit function as follows -- depending on the data, model, and number of iterations requested, this can take anywhere from a few minutes to days. Current experience suggests 500 iterations is enough to get an idea of what is going on, but more are necessary for robust inference.

``` fitting
fit<-ctStanFit(long,checkm,iter=500,chains=2,fit=T,plot=T,
  densehyper=F)
```

The plot argument allows for plotting of sampling chains in real time, which is useful for slow models to ensure that sampling is proceeding in a functional manner. The densehyper argument may be set to TRUE to estimate the priors for the correlation between parameters, which may allow somewhat better priors for subject level parameters to be estimated, but also tends to slow down sampling substantially.

### Output

After fitting, the standard rstan output functions such as summary and extract are available, and the shinystan package provides an excellent browser based interface. The parameters which are likely to be of most interest in the output all begin with an "output" prefix, followed by either "hmean" for hyper (population) mean, or "hsd" for hyper standard deviation. Subject specific parameters are denoted by the matrix they are from, then the first index represents the subject id, followed by standard matrix notation. For example, the 2nd row and 1st column of the DRIFT matrix for subject 8 is "DRIFT\[8,2,1\]". Parameters are all returned in the form used for internal calculations -- that is, variance covariance matrices are returned as such, rather than the lower-triangular standard deviation and cholesky correlation matrices required for input. The exception to this are the time independent predictor effects, prefixed with "output\_tip\_", for which a linear effect of a change of 1 on the predictor is approximated. So although "output\_tip\_TI1" is only truly linear with respect to internal parameterisations, we approximate the linear effect by averaging the effect of a score of +1 or -1 on the predictor, on the population mean. For any subject that substantially differs from the mean, or simply when precise absolute values of the effects are required (as opposed to general directions), they will need to be calculated manually.

``` output
library("shinystan")
launch_shinystan(fit)
```

### Convert from wide data

Data can be converted from the wide format data used for the OpenMx based ctsem approach as follows:

``` convertdata
#specify some ctModel called mymodel, including a Tpoints argument
long<-ctWideToLong(mydata,mymodel$Tpoints,
n.manifest=mymodel$n.manifest, manifestNames = mymodel$manifestNames, 
  n.TDpred=mymodel$n.TDpred, TDpredNames = mymodel$TDpredNames,
  n.TIpred=mymodel$n.TIpred, TIpredNames = mymodel$TIpredNames)
```
