#' Define a ctsem model
#' 
#' This function is used to specify a continuous time structural equation model, 
#' which can then be fit to data with function \code{\link{ctFit}}.
#' 
#' @param type character string. Use 'ct' (continuous time) or 'dt'
#' (discrete time) for fitting with \code{\link{ctFit}}.
#' These return a modern ctsem model object whose \code{pars} data frame can
#' also be edited through the pars-backed \code{matrices} view.
#' For legacy workflows, 'omx' returns a matrix-list model object that can
#' be edited directly and then converted via \code{\link{ctModelConvertOMX}}.
#' 
#' @param n.manifest Number of manifest indicators per individual at each measurement occasion / time point.  
#' Manifest variables are included as the first element of the wide data matrix, with all the 1:n.manifest manifest variables 
#' at time 1 followed by those of time 2, and so on.
#' 
#' @param n.latent Number of latent processes.
#' 
#' @param LAMBDA n.manifest*n.latent loading matrix relating latent to manifest variables, 
#' with latent processes 1:n.latent along the columns, and manifest variables
#' 1:n.manifest in the rows.
#' 
#' @param manifestNames n.manifest length vector of manifest variable names as they appear in the data structure, 
#' without any _Tx time point suffix that may be present in wide data.  Defaults to Y1, Y2, etc.
#' 
#' @param manifesttype n.manifest length vector of manifest variable types: 0
#' (the default) for continuous Gaussian variables, 1 for binary (coded 0/1), 2
#' for ordinal (consecutive integers from 1), 3 for counts (non-negative
#' integers, Poisson with a log link, so the latent process gives the log rate)
#' and 4 for censored (Gaussian, observed between \code{censormin} and
#' \code{censormax}). For a count the \code{MANIFESTVAR} diagonal stays free
#' and is a log-scale dispersion, giving a Poisson-lognormal observation with
#' variance \eqn{\mu + \mu^2(e^{\sigma^2}-1)}. Fix it to zero for a plain
#' Poisson. A censored value at a limit contributes the probability of being
#' at or beyond that limit. Types other than 0 need the julia backend, except
#' binary, which stan also fits.
#' @param censormin,censormax n.manifest length numeric vectors giving the lower
#' and upper censoring limits of each \code{manifesttype = 4} variable, ignored
#' for other types. The defaults are \code{-Inf} and \code{Inf}, so a variable
#' censored on one side keeps the other limit infinite. The limits are fixed,
#' known values.
#' @param ncategories n.manifest length integer vector giving the number of
#' categories of each ordinal manifest variable, ignored for other types.
#' Required when any \code{manifesttype} is 2, by \code{ctGenerate} as well as
#' \code{ctFit}.
#' 
#' @param latentNames n.latent length vector of latent variable names 
#' (used for naming parameters, defaults to eta1, eta2, etc).
#' 
#' @param id character string denoting column name containing subject identification variables. 
#' id data may be of any form, though will be coerced internally to an integer sequence rising from 1.
#' Or a character vector of column names for nested grouping levels, innermost
#' first: \code{id = c('subject','study')} places each subject in one study.
#' The \code{indvarying} field of a cell (see \code{\link{ctPar}}) sets the
#' levels its parameter varies at: \code{TRUE} for every level, \code{FALSE}
#' for none, or the names of levels. For example \code{'cint||study'} gives an
#' intercept that varies between studies, and \code{'cint||subject,study'} (or
#' \code{TRUE}) one that varies at both levels. \code{model$pars$indvarying}
#' holds these values, one per cell, and can be edited directly. A cell's
#' \code{sdscale} applies at every level it varies at. Intercepts
#' (\code{T0MEANS}, \code{MANIFESTMEANS}, \code{CINT}) that are left
#' unspecified vary at every level. Levels above the subject require
#' \code{backend='julia'} with \code{intoverpop='laplace'}, and the levels must
#' be strictly nested.
#' @param time character string denoting column name containing timing data. Timing data must be numeric.
#' 
#' @param T0VAR lower triangular n.latent*n.latent cholesky matrix of latent process initial variance / covariance. 
#' "auto" freely estimates all parameters.
#' 
#' @param T0MEANS n.latent*1 matrix of latent process means at first time point, T0. 
#' "auto" freely estimates all parameters.
#' 
#' @param MANIFESTMEANS n.manifest*1 matrix of manifest intercept parameters.  
#' "auto" frees all parameters.
#' 
#' @param MANIFESTVAR lower triangular n.manifest*n.manifest cholesky matrix of variance / covariance 
#' between manifests at each measurement occasion (i.e. measurement error / residual).  
#' "auto" freely estimates variance parameters, 
#' and fixes covariances between manifests to 0. "free" frees all values, including covariances.
#'
#' @param ASYMPTOTES Matrix with one row per binary indicator, in manifest
#' order, and two columns giving the lower and upper asymptotes of its
#' response probability. The probability is then lower + (upper - lower)
#' times the logistic of the linear predictor. \code{NA} means no asymptote
#' on that side, and the default of all \code{NA} gives the two parameter
#' logistic model. A parameter name in the lower column estimates a guessing
#' probability (three parameter model), in the upper column a ceiling below
#' one, and in both the four parameter model. A number between 0 and 1 fixes
#' the asymptote, and a name repeated across rows estimates one asymptote for
#' those items. A length 2 vector applies to every binary indicator, so
#' \code{c('guess', NA)} estimates one guessing parameter for the whole test.
#' Julia backend only.
#'
#' @param RESPONSE Matrix with one row per manifest variable and two columns,
#' setting the response scale of censored indicators (manifesttype 4). Rows
#' for other indicators keep their defaults. Column 1, \code{halfwidth}, is
#' half the width of the interval a recorded value stands for. A value
#' \code{y} means the latent response lies between the transformed
#' \code{y - halfwidth} and \code{y + halfwidth}, a value at \code{censormin}
#' that it lies below \code{censormin + halfwidth}, and a value at
#' \code{censormax} that it lies above \code{censormax - halfwidth}. The
#' default, 0, treats values as exact. For a rating scale use half the step
#' between categories. A step that differs between studies can be written as
#' a time independent predictor in the cell. Column 2, \code{power}, sets the
#' transform from the recorded scale to the latent response,
#' \eqn{G(r) = L + ((1 + r - L)^p - 1)/p}, with \eqn{L} the censormin (0 if
#' there is none). The default, 1, is linear. Values below 1 make the latent
#' response change less per recorded unit towards the top of the scale, and 0
#' gives the log transform. A parameter name estimates the power, and a name
#' shared across items or studies estimates one power from all of them.
#' \code{NA} keeps the default, and a length 2 vector fills every censored
#' row, so \code{c(0.5, 'power')} suits integer data with one estimated
#' transform. Julia backend only.
#'
#' @param DRIFT n.latent*n.latent DRIFT matrix of continuous auto and cross effects, 
#' relating the processes over time. 
#' "auto" freely estimates all parameters.
#' 
#' @param CINT n.latent * 1 matrix of latent process intercepts, allowing for non 0 
#' asymptotic levels of the latent processes. Generally only necessary for additional trends and more complex dynamics.
#' 
#' @param DIFFUSION lower triangular n.latent*n.latent cholesky matrix of diffusion process 
#' variance and covariance (latent error / dynamic innovation).
#' "auto" freely estimates all parameters.
#' 
#' @param n.TIpred Number of time independent predictors. 
#' Each TIpredictor is inserted at the right of the data matrix, after the time intervals.
#' 
#' @param TIpredNames n.TIpred length vector of time independent predictor variable names,
#' as they appear in the data structure.  Default names are TI1, TI2, etc.
#' Or a formula of such names and smooth terms, such as
#' \code{~ s(age, k = 8) + female}. A smooth term is a penalised spline of a
#' subject-level covariate, with its smoothness estimated from the data
#' (\code{backend = 'julia'}). It enters the model as \code{k - 1} predictor
#' columns, \code{age_lin} (the linear trend) and \code{age_s1}, ..., built
#' from the data when the model is fitted. To give a parameter the smooth
#' effect, name the covariate in its cell's \code{tipreds} field
#' (\code{'cint||||age'}). Under \code{tipredDefault = TRUE}, cells without a
#' \code{tipreds} field get the smooth effect if their parameter is free and
#' varies between subjects.
#'
#' @param tipredDefault Logical. TRUE sets any parameters with unspecified time independent
#' predictor effects to have effects estimated, FALSE fixes the effect to zero unless individually specified.
#' A cell's effects are unspecified only when its \code{tipreds} field is left
#' out. An empty field (\code{'mypar||||'}) means no predictor effects, which
#' excludes a parameter under \code{tipredDefault = TRUE}. The setting is
#' stored on the model as \code{$tipredDefault}.
#' 
#' @param n.TDpred Number of time dependent predictor variables in the dataset.  
#' 
#' @param TDpredNames n.TDpred length vector of time dependent predictor variable names, 
#' as they appear in the data structure, without any _Tx time point suffix that may appear in wide data.  
#' Default names are TD1, TD2, etc.
#' 
#' @param TDPREDEFFECT n.latent*n.TDpred matrix of effects from time dependent predictors to latent processes.
#' Effects from 1:n.TDpred columns TDpredictors go to 1:n.latent rows of latent processes.
#' "auto" freely estimates all parameters.
#' @param TDPREDMEANS Legacy argument for time-dependent predictor means, mainly to
#' support data generation workflows. Expected dimensions are
#' (n.TDpred * Tpoints) rows and 1 column. "auto" creates a free-label matrix.
#' @param TDPREDVAR Legacy argument for time-dependent predictor covariance
#' structure, mainly to support data generation workflows. Expected dimensions are
#' (n.TDpred * Tpoints) square. "auto" creates a free-label matrix.
#' 
#' @param Tpoints Number of time points, or measurement occasions, in the data.  This will generally be the maximum 
#' number of time points for a single individual, but may be one extra if sample relative time intervals are used, 
#' see \code{\link{ctIntervalise}}. 
#' 
#' @param PARS for types 'ct' and 'dt' only. May be of any structure, only needed to contain extra parameters for certain non-linear models.
#' @param silent Suppress all output to console.
#' @param driftScaled for type 'ct' only: which of "DIFFUSION" and "CINT" to
#' parameterise by what they add over one time unit, given the process's own
#' drift diagonal \code{a}. A diffusion sd \code{s} becomes
#' \code{s * sqrt(2 a / expm1(2 a))}, so that \code{s^2} is the variance the
#' process gains from its own noise over one time unit. An intercept \code{u}
#' becomes \code{u * a / expm1(a)}. The model and its likelihood are
#' unchanged. The rescaled parameters depend less on the drift diagonal, which
#' helps uncertainty estimation and sampling when fast mean reversion and
#' larger noise are hard to tell apart. Only the diagonal of DIFFUSION is
#' rescaled, so its correlations keep their usual meaning. Each rescaled
#' parameter, with its name, transform, random effects and predictor effects,
#' moves to a new cell of \code{PARS}, and its original cell becomes an
#' expression of that cell and the drift diagonal. Not available with
#' \code{covmattransform = 'cholesky'}.
#'
#' @details A matrix cell can be one of four things:
#'
#' \itemize{
#'   \item a \strong{number}, which fixes the cell;
#'   \item a \strong{parameter name}, which frees it. The same name in two
#'     cells is one parameter, which is how equality constraints are written;
#'   \item a name \strong{with fields}, giving the transform, the levels of
#'     \code{id} at which the parameter varies, the prior scale of that
#'     variation, and the time independent predictors that affect it;
#'   \item an \strong{expression}, computed from the latent states, the time
#'     dependent predictors and the \code{PARS} cells.
#' }
#'
#' Fields can be given in a fixed order separated by \code{|}
#' (\code{'mm||TRUE|0.5'}), or by name
#' (\code{'mm, indvarying=TRUE, sdscale=0.5'}). With nested ids,
#' \code{indvarying} may name levels, as in \code{'mm||study'} or
#' \code{'mm, indvarying=c(subject, study)'}. A cell containing \code{|} is
#' read in the ordered form. Otherwise a cell with a field of the form
#' \code{name=value} is read in the named form, and an unknown field name is
#' an error. \code{\link{ctPar}} writes these strings from named arguments.
#'
#' Every name in an expression must be a latent state, a time dependent
#' predictor, or a parameter declared in \code{PARS}. For example
#' \code{DRIFT = matrix('-exp(dr)*eta1')} needs \code{PARS = c('dr')}. To
#' transform a single parameter, the transform field is shorter:
#' \code{DRIFT = matrix('dr|-exp(param)')}. With the julia backend a cell may
#' also use a time independent predictor by name, which gives the subject's
#' value of it, as in \code{MANIFESTVAR = 'rsd + mvar'} with
#' \code{PARS = c('mvar|1e-10 + 5 * log1p_exp(2 * param)')}. Missing values of
#' such a predictor are imputed or sampled like other missing predictor values
#' (\code{optimcontrol$tipredMissing}). Population summaries evaluate the cell
#' at a predictor value of zero. A cell cannot refer to a manifest variable or
#' to the time interval. Predictor \emph{effects} on a parameter are set with
#' the \code{tipreds} field.
#'
#' A latent state or time dependent predictor name in a cell refers to that
#' state or predictor, and any number of matrices may refer to the same one.
#' A name declared in \code{PARS} likewise refers to the \code{PARS} cell
#' wherever else it appears, and that cell holds the parameter and its
#' transform.
#'
#' Any other name repeated across cells is one parameter, with one value and
#' one random effect, so it needs one transform. Cells that would transform it
#' differently are refused, including cells whose \emph{default} transforms
#' differ. The \code{DRIFT} diagonal is bounded negative and the off-diagonals
#' are not, so \code{DRIFT = matrix('a', 2, 2)} is refused. To use one
#' parameter in cells that need different transforms, declare it in
#' \code{PARS} (\code{PARS = c('a')}) and write the transform in each cell, as
#' in \code{DRIFT[1,1] = '-exp(a)'} and \code{DIFFUSION[1,1] = 'exp(a)'}. If
#' one transform suits every cell, give it in \code{PARS}
#' (\code{PARS = c('a|-exp(param)')}) and write only the name in each cell.
#' Two \code{PARS} cells with the same name are one parameter and must have
#' the same transform.
#'
#' A predictor effect belongs to the parameter, so cells holding one parameter
#' must have the same \code{tipreds}. Naming an effect
#' (\code{'mypar||||age=shared'}) gives every parameter with that effect name,
#' for that predictor, one shared coefficient. Because the coefficient acts on
#' the raw scale, parameters with different transforms cannot share it. The
#' same effect name under two predictors gives two coefficients. A number
#' instead of a name (\code{'mypar||||age=0.5'}) fixes the effect, shifting
#' each subject's raw parameter value by 0.5 times their \code{age} before the
#' transform. Fixed effects are fitted by \code{backend = 'julia'}
#' and used by \code{ctGenerate}. The stan backend refuses them.
#'
#' A cell that is exactly a latent state or time dependent predictor name is
#' read as a reference to it, with a warning, since a parameter of that name
#' would give a different model. Write \code{'state[1]'} or
#' \code{'tdpreds[rowi, 1]'} to avoid the warning. In \code{PARS} such a cell
#' is an error, because a bare name there declares a parameter. For a
#' \code{PARS} cell that depends on a state, write the state in an
#' expression, such as \code{'log1p(exp(eta1))'}.
#'
#' @seealso \code{\link{ctPar}} for the fields a matrix cell can carry,
#' \code{\link{ctFit}} for fitting the result.
#'
#' @examples
#'  ### Frequentist example:
#'  ### impulse and level change time dependent predictor 
#'  ### example from Driver, Oud, Voelkle (2015)
#'  data('ctExample2')
#'  tdpredmodel <- ctModel(type='omx', n.manifest = 2, n.latent = 3, n.TDpred = 1, 
#'  Tpoints = 8, manifestNames = c('LeisureTime', 'Happiness'), 
#'  TDpredNames = 'MoneyInt', 
#'  latentNames = c('LeisureTime', 'Happiness', 'MoneyIntLatent'),
#'  LAMBDA = matrix(c(1,0, 0,1, 0,0), ncol = 3))
#'
#'  tdpredmodel$DIFFUSION[, 3] <- 0
#'  tdpredmodel$DIFFUSION[3, ] <- 0
#'  tdpredmodel$T0VAR[3, ] <- 0
#'  tdpredmodel$T0VAR[, 3] <- 0
#'  tdpredmodel$CINT[3] <- 0
#'  tdpredmodel$T0MEANS[3] <- 0
#'  tdpredmodel$DRIFT[3, ] <- 0
#'  
#'  
#' ###Bayesian example:
#' model<-ctModel(type='ct',
#' n.latent=2, latentNames=c('eta1','eta2'),
#' n.manifest=2, manifestNames=c('Y1','Y2'),
#' n.TDpred=1, TDpredNames='TD1', 
#' n.TIpred=3, TIpredNames=c('TI1','TI2','TI3'),
#' LAMBDA=diag(2))
#'
#' 
#' @return For \code{type = 'ct'} or \code{'dt'}, a model of class
#'   \code{ctStanModel}, with one row per matrix cell in \code{$pars} and the
#'   same cells as matrices in \code{$matrices}. Either can be edited. For
#'   \code{type = 'omx'}, a list of matrices of class \code{ctsemInit}, for
#'   \code{\link{ctGenerate}} or \code{\link{ctModelConvertOMX}}.
#' @export

ctModel<-function(LAMBDA, type='ct',n.manifest = 'auto', n.latent='auto', Tpoints=NULL, 
  manifestNames='auto', manifesttype=rep(0,nrow(LAMBDA)), ncategories=NULL,
  censormin=NULL, censormax=NULL,
  latentNames='auto', id='id',time='time', silent=FALSE,
  T0VAR="auto", T0MEANS="auto", MANIFESTMEANS="auto", MANIFESTVAR="diag",
  ASYMPTOTES=NA, RESPONSE=NA,
  DRIFT="auto", CINT=0, DIFFUSION="auto",
  n.TDpred='auto', TDpredNames='auto', TDPREDEFFECT="auto", TDPREDMEANS="auto", TDPREDVAR="auto",
  n.TIpred='auto', TIpredNames='auto', tipredDefault=TRUE,
  PARS=NULL, driftScaled=character()){
  
  # keep track of whether Tpoints was explicitly provided by user
  TpointsProvided <- !is.null(Tpoints) && !is.na(Tpoints[1])

  # A covariance given as ctCov(): the cells for the construction this model
  # will use -- 'z' for ct and dt models, a Cholesky factor for the old omx
  # form -- with the covariance remembered so the cells can be rewritten if
  # covmattransform changes (.ctCovRefresh()).
  covinput <- list()
  for(nm in .ctCovSystemMatrices){
    value <- get(nm)
    if(inherits(value, 'ctCov')){
      cov <- attr(value, 'covariance')
      tf <- if(type %in% 'omx') 'cholesky' else 'z'
      cells <- .ctCovCells(cov, tf, what = nm)
      covinput[[nm]] <- list(cov = cov, cells = unname(cells), covmattransform = tf)
      assign(nm, matrix(as.numeric(cells), nrow(cells), ncol(cells)))
    }
  }
  
  #get dimensions
  if(is.null(n.manifest) || is.null(n.latent) || all(n.manifest %in% 'auto') || all(n.latent %in% 'auto')){
    if(is.matrix(LAMBDA)){
      if(!silent) message('System dimensions inferred from LAMBDA')
      n.manifest <- nrow(LAMBDA)
      n.latent <- ncol(LAMBDA)
    } else {
      if(!all(manifestNames %in% 'auto') && !all(latentNames %in% 'auto')){
        if(!silent) message('System dimensions inferred from manifestNames and latentNames')
        n.manifest <- length(manifestNames)
        n.latent <- length(latentNames)
      } else stop('LAMBDA must either a matrix, n.manifest and n.latent must be specified, or manifestNames and latentNames must be specified!')
    }
  }
  
  if(is.null(n.TDpred) || all(n.TDpred %in% 'auto')){
    if(is.matrix(TDPREDEFFECT)){
      if(!silent) message('n.TDpred inferred from TDPREDEFFECT')
      n.TDpred <- ncol(TDPREDEFFECT)
    } else if(!all(TDpredNames %in% 'auto')){
      if(!silent) message('n.TDpred inferred  inferred from TDpredNames')
      n.TDpred <- length(TDpredNames)
    } else n.TDpred <- 0
  }
  
  # Smooth terms (R/ctTIpredSmooth.R): their basis columns are the predictors.
  TIsmooths <- NULL
  if(inherits(TIpredNames, 'formula')){
    if(type %in% 'omx') stop("Smooth TI predictor terms need type = 'ct' or 'dt'.",
      call. = FALSE)
    parsed <- .ctTIpredParseFormula(TIpredNames)
    TIpredNames <- parsed$names
    TIsmooths <- parsed$smooths
  }

  if(is.null(n.TIpred) || all(n.TIpred %in% 'auto')){
    if(!all(TIpredNames %in% 'auto')){
      if(!silent) message('n.TIpred inferred  inferred from TIpredNames')
      n.TIpred <- length(TIpredNames)
    } else n.TIpred <- 0
  } 
  
  
  
  #####FUNCTIONS
  
  
  checkSymmetry<-function(x){  #this checks the symmetry of matrix x, if not symmetric it stops and outputs an error
    if(isSymmetric(x)==F){
      stop(paste0(substitute(x)," must be symmetric"))
    }
  }
  
  
  if(!type %in% c('ct','dt','stanct','standt','omx')) stop('type must be either ct, dt, or omx!')
  
  msg_matrix <- function(mat, label = "") {
    lines <- capture.output(print(mat))
    header <- if (nzchar(label)) paste0(label, ":\n") else ""
    message(header, paste(lines, collapse = "\n"))
  }
  
  
  
  ###### RUN SEQUENCE
  
  if(is.null(Tpoints)) Tpoints<-3
  # if(type=='omx' && is.null(Tpoints)) stop('Tpoints must be specified for type="omx"')

  
  #names
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  if(all(latentNames=='auto')) latentNames=paste0('eta',1:n.latent)
  
  if(length(manifestNames) != n.manifest) stop(paste0(
    "Length of manifestNames (",length(manifestNames),") does not equal n.manifest (",n.manifest,")!")) 
  if(length(latentNames) != n.latent) stop(paste0(
    "Length of latentNames (",length(latentNames),") does not equal n.latent(",n.latent,")!"))
  
  if(n.TDpred > 0){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
    if(length(TDpredNames) != n.TDpred) stop(paste0(
      "Length of TDpredNames (",length(TDpredNames),") does not equal n.TDpred(",n.TDpred,")!")) 
  } else {
    TDpredNames=c()
  }
  
  if(n.TIpred > 0){
    if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
    if(length(TIpredNames) != n.TIpred) stop(paste0(
      "Length of TIpredNames (",length(TIpredNames),') does not equal n.TIpred(',n.TIpred,')!')) 
  } else {
    TIpredNames=c()
  }
  
  mats <- .ctMatricesList()
  for(m in names(mats$base)){
    if(!exists(m, inherits = FALSE)) next
    if(!is.null(get(m, inherits = FALSE))){ #if the matrix is specified
      val <- get(m, inherits = FALSE)
      if(!all(val %in% c('auto','diag')) && !is.matrix(val)){ #and not auto or a matrix
        if(m %in% 'PARS') mat <- matrix(PARS,ncol=1) else {
          mat<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames, #set the basic structure
            manifestNames=manifestNames,latentNames=latentNames,matrixname=m,n.latent=n.latent,
            n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
        }
        counter <- 0 #then start filling with vector
        if(length(val) != length(mat) && length(val) !=1) stop(paste0(m,' needs ',length(mat),' values for ',nrow(mat),' * ',ncol(mat),' matrix, but has ',length(val)))
        if(length(val)==1 && length(mat) > 1) {
          if(!silent) message(m,' specified via single value -- filling ',nrow(mat),' * ',ncol(mat),' matrix:')
          val <- rep(val,length(mat))
        } else if(!silent) message(paste0(m,' vector spec input rowwise into ',nrow(mat),' * ',ncol(mat),' matrix:'))
        for(ri in 1:nrow(mat)){
          for(ci in 1:ncol(mat)){
            counter <- counter + 1
            mat[ri,ci] <- val[counter]
          }
        }
        
        if(!silent) msg_matrix(mat)#,right=TRUE)
        assign(m, mat)
      }
    }
  }
  
  
  # An NA cell says nothing, and the four `[1]=="auto"` tests below compare
  # against it, so a single NA in the first cell of one of those matrices came
  # out as "missing value where TRUE/FALSE needed" -- naming neither the matrix
  # nor the cell. Refuse it here, where both are still known. A cell is a
  # number, a name, or an expression; there is no reading under which NA is one
  # of those, and `0` is what an absent cell usually means.
  for(m in c('T0MEANS','T0VAR','LAMBDA','DRIFT','DIFFUSION','MANIFESTVAR',
    'MANIFESTMEANS','CINT','TDPREDEFFECT','PARS')){
    mat <- get0(m, ifnotfound=NULL)
    if(is.null(mat) || (length(mat)==1 && !is.na(mat[1]) && mat[1]=='auto')) next
    # The text "NA" -- what a character conversion or a text field gives --
    # is the same non-cell, and was read as neither a number nor a name.
    bad <- which(is.na(mat) | trimws(as.character(mat)) %in% 'NA')
    if(length(bad)){
      idx <- if(is.matrix(mat)) paste0('[',(bad[1]-1) %% nrow(mat) + 1,',',
        (bad[1]-1) %/% nrow(mat) + 1,']') else paste0('[',bad[1],']')
      stop(m,idx,' is NA. A matrix cell must be a number, a parameter name, or ',
        'an expression -- write 0 for a cell that is fixed to zero.',
        call.=FALSE)
    }
  }

  #matrices
  if(T0MEANS[1]=="auto") T0MEANS<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
    manifestNames=manifestNames,latentNames=latentNames,matrixname="T0MEANS",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(ncol(T0MEANS)>1) stop("Specified T0MEANS matrix has more than one column")
  if(nrow(T0MEANS)!=n.latent) stop ("Specified T0MEANS matrix rows not equal to n.latent")
  rownames(T0MEANS)=latentNames
  
  
  
  
  
  if(T0VAR[1]=="auto") T0VAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
    manifestNames=manifestNames,latentNames=latentNames,matrixname="T0VAR",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  
  if(nrow(T0VAR)!=n.latent | ncol(T0VAR)!=n.latent) stop("Dimensions of T0VAR matrix are not n.latent * n.latent")
  dimnames(T0VAR)=list(latentNames,latentNames)
  
  
  
  
  
  if(MANIFESTMEANS[1]=="auto") MANIFESTMEANS<- ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
    manifestNames=manifestNames,latentNames=latentNames,matrixname="MANIFESTMEANS",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(ncol(MANIFESTMEANS)>1) stop("Specified MANIFESTMEANS matrix has more than one column")
  if(nrow(MANIFESTMEANS)!=n.manifest) stop ("Specified MANIFESTMEANS matrix rows not equal to n.manifest")
  rownames(MANIFESTMEANS)=manifestNames
  
  
  

    if(all(MANIFESTVAR %in% "auto")) {
      MANIFESTVAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
        manifestNames=manifestNames,latentNames=latentNames,matrixname="MANIFESTVAR",n.latent=n.latent,
        n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    }
    
    if(all(MANIFESTVAR %in% 'diag')){
      MANIFESTVAR<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,
        manifestNames=manifestNames,latentNames=latentNames,matrixname="MANIFESTVAR",n.latent=n.latent,
        n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
      MANIFESTVAR[row(MANIFESTVAR) != col(MANIFESTVAR)] <- 0
    }
    
    if(nrow(MANIFESTVAR)!=n.manifest | ncol(MANIFESTVAR)!=n.manifest) stop("Dimensions of MANIFESTVAR matrix are not n.manifest * n.manifest")
  
  dimnames(MANIFESTVAR)=list(manifestNames,manifestNames)
  
  
  
  if(all(DRIFT %in% "auto")) DRIFT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DRIFT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(all(DRIFT %in% 'diag')) {
    DRIFT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DRIFT",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    DRIFT[row(DRIFT) != col(DRIFT)] <- 0
  }
  if(nrow(DRIFT)!= n.latent | ncol(DRIFT)!= n.latent) stop("Dimensions of DRIFT matrix are not n.latent * n.latent")
  dimnames(DRIFT)=list(latentNames,latentNames)
  
  
  
  
  if(CINT[1]=="auto") CINT<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="CINT",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(ncol(CINT)!=1 | nrow(CINT)!=n.latent) stop("Dimensions of CINT are not n.latent * 1")
  rownames(CINT)=latentNames
  
  
  
  
  if(all(DIFFUSION %in% "auto")) DIFFUSION<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DIFFUSION",n.latent=n.latent,
    n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
  if(all(DIFFUSION %in% 'diag')) {
    DIFFUSION<-ctLabel(TDpredNames=TDpredNames,TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="DIFFUSION",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    DIFFUSION[row(DIFFUSION) != col(DIFFUSION)] <- 0
  }
  if(nrow(DIFFUSION)!=n.latent | ncol(DIFFUSION) != n.latent) stop("Dimensions of DIFFUSION are not n.latent * n.latent")
  dimnames(DIFFUSION)=list(latentNames,latentNames)
  
  if(n.TDpred > 0){
    if(all(TDPREDEFFECT=="auto"))  {
      TDPREDEFFECT<-ctLabel(TDpredNames=TDpredNames,
        TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDPREDEFFECT",n.latent=n.latent,
        n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    }
    if(nrow(TDPREDEFFECT) != n.latent || ncol(TDPREDEFFECT) != n.TDpred) stop("Dimensions of TDPREDEFFECT are not n.latent * n.TDpred")
    dimnames(TDPREDEFFECT)=list(latentNames,TDpredNames)
    
    if(all(TDPREDMEANS=="auto")) TDPREDMEANS <- ctLabel(TDpredNames=TDpredNames,
      TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDPREDMEANS",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    if(all(TDPREDVAR=="auto")) TDPREDVAR <- ctLabel(TDpredNames=TDpredNames,
      TIpredNames=TIpredNames,manifestNames=manifestNames,latentNames=latentNames,matrixname="TDPREDVAR",n.latent=n.latent,
      n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,Tpoints=Tpoints)
    
    if(nrow(TDPREDMEANS) != n.TDpred * Tpoints || ncol(TDPREDMEANS) != 1) stop("Dimensions of TDPREDMEANS are not (n.TDpred*Tpoints) * 1")
    if(nrow(TDPREDVAR) != n.TDpred * Tpoints || ncol(TDPREDVAR) != n.TDpred * Tpoints) stop("Dimensions of TDPREDVAR are not (n.TDpred*Tpoints) * (n.TDpred*Tpoints)")
  } else {
    TDPREDEFFECT <- NULL
    TDPREDMEANS <- NULL
    TDPREDVAR <- NULL
  }
  
  
  
  
  
  ###model checks
  
  #lower triangular check
  for(tempmatname in c('T0VAR','MANIFESTVAR','DIFFUSION')){
    assign('tempmat',get(tempmatname))
    if(is.null(tempmat)) next
    if(all(suppressWarnings(as.numeric(tempmat[upper.tri(tempmat)])) %in% 0)) next
    warning(paste0(tempmatname,' is not lower triangular! Covariance type matrices should usually be specified in the appropriate lower-triangular form.'))
  }
  
  # mvaroffdiag=MANIFESTVAR[!diag(1,nrow(MANIFESTVAR))]
  # if(nrow(MANIFESTVAR) > 1 && !all(mvaroffdiag %in% 0)) stop('MANIFESTVAR should be diagonal!')
  
  
  if(any(dim(LAMBDA)!=c(n.manifest,n.latent))) stop("Incorrect LAMBDA structure specified - check number or rows and columns")
  dimnames(LAMBDA)=list(manifestNames,latentNames)
  
  sapply(c(manifestNames,latentNames,TDpredNames,TIpredNames, time,id),function(x){
    if(grepl('\\W',x)) stop(paste0(x,' contains symbols, variable names must be alphanumerics only please!'))
  })
  
  for(names in c('manifestNames','TIpredNames','TDpredNames','latentNames')){
    if(any(duplicated(get(names)))) stop(paste0('Duplicate names in ',names))
  }
  
  # No advice to move the location into CINT. It was needed when the link was
  # applied to the latent alone, so a manifest intercept could not shift a
  # probability; it is not now -- `adjoint_binary.jl` forms the linear predictor
  # as `mu + lambda'x` and passes that to the link, so MANIFESTMEANS does the
  # job directly for binary, count and censored indicators. The one case where
  # the advice was true is ordinal, and `.ctThresholdMatrix()` now settles that
  # by construction rather than by asking the user to.

  if(any(!manifesttype %in% 0:4)) stop('manifesttype must be 0 (continuous), 1 (binary), 2 (ordinal), 3 (count) or 4 (censored)')
  censorlimits <- .ctCheckCensorLimits(censormin, censormax, manifesttype,
    manifestNames)
  censormin <- censorlimits$min
  censormax <- censorlimits$max
  THRESHOLDS <- NULL
  ncategories <- .ctCheckNcategories(ncategories, manifesttype, manifestNames)
  if(any(manifesttype %in% 2)){
    # The first threshold is fixed at zero and the location lives in
    # MANIFESTMEANS -- see `.ctThresholdMatrix()`. Nothing to warn about
    # afterwards, because the pair can no longer both be free.
    THRESHOLDS <- .ctThresholdMatrix(ncategories, manifesttype, manifestNames)
  }
  ASYMPTOTES <- .ctCheckAsymptoteMatrix(ASYMPTOTES, manifesttype,
    manifestNames)
  RESPONSE <- .ctCheckResponseMatrix(RESPONSE, manifesttype, manifestNames)
  
  
  
  
  #rewrite last 2 pieces of code into a single element:
  completemodel <- list(`n.manifest`=n.manifest, `n.latent`=n.latent, `n.TDpred`=n.TDpred, `n.TIpred`=n.TIpred, 
    `Tpoints`=Tpoints, `LAMBDA`=LAMBDA, `manifestNames`=manifestNames, `latentNames`=latentNames, 
    `TDpredNames`=TDpredNames, `TIpredNames`=TIpredNames, `T0VAR`=T0VAR, `T0MEANS`=T0MEANS, 
    `MANIFESTMEANS`=MANIFESTMEANS, `MANIFESTVAR`=MANIFESTVAR, `DRIFT`=DRIFT, `CINT`=CINT, 
    `DIFFUSION`=DIFFUSION, `TDPREDEFFECT`=TDPREDEFFECT, `TDPREDMEANS`=TDPREDMEANS, `TDPREDVAR`=TDPREDVAR, `PARS`=PARS, 
     `id`=id, `time`=time, `manifesttype`=manifesttype,
     `censormin`=censormin, `censormax`=censormax,
     `ncategories`=ncategories,
     `THRESHOLDS`=THRESHOLDS, `ASYMPTOTES`=ASYMPTOTES, `RESPONSE`=RESPONSE)
  if(length(TIsmooths)) completemodel$TIsmooths <- TIsmooths
  
  
  
  
  if(type %in% c('ct','dt','stanct','standt')) completemodel<-ctModelConvertOMX(completemodel,type=type,tipredDefault= tipredDefault)
  if(length(driftScaled)) {
    if(type %in% 'omx') stop("driftScaled needs type = 'ct'.", call. = FALSE)
    completemodel <- .ctModelDriftScaled(completemodel, driftScaled)
  }
  if(type %in% c('omx')) class(completemodel)<-'ctsemInit'
  if(TpointsProvided) completemodel$Tpoints <- Tpoints
  if(length(covinput) && !type %in% 'omx') completemodel$covinput <- covinput
  
  return(completemodel)
}
