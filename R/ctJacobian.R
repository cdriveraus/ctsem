# ctm=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','state[3]',0, 'd12','d22',0,0,0,0),3,3),type='ct')
# ctm2=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','PARS[1,1]',0, 'd12','d22',0,0,0,0),3,3),PARS=matrix('state[3]'),type='ct')

###unfold self referential list of matrices
unfoldmats <- function(ml){
  for(mati in names(ml)){
    if(prod(dim(ml[[mati]])) > 0){ #if not a 0 matrix
      for(ri in 1:nrow(ml[[mati]])){
        for(ci in 1:ncol(ml[[mati]])){
          counter=0
          while(counter < 10 && grepl(paste0('\\b(',paste0(names(ml),collapse='|'),')\\b\\['), ml[[mati]][ri,ci])){ #if system matrix referenced, unfold
            counter = counter +1 #prevent infinite loops
            
            items = regmatches(ml[[mati]][ri,ci], #extract one or more references
              gregexpr(
                paste0('\\b(',paste0(names(ml),collapse='|'),')\\b(?=\\[).*?(?<=\\])'),
                ml[[mati]][ri,ci], perl=TRUE)
            )[[1]]
            
            for(itemi in 1:length(items)){ #replace one or more references
              ml[[mati]][ri,ci] <- gsub(pattern = items[itemi], replacement = eval(parse(text=paste0('ml$',items[itemi]))),x = ml[[mati]][ri,ci], fixed=TRUE)
            }
          }
        }
      }
    }
  }
  return(ml)
}

# Logistics in the form whose symbolic derivative stays finite.
#
# Deriv's quotient rule turns `x / (1 + exp(a))` -- and `inv_logit(a)`, whose
# body it differentiates when it has no rule for it -- into terms shaped
# `exp(a) / (1 + exp(a))^2`. Its value is right, 0 once the square overflows
# at a > 354, but any derivative on `a` then makes the partial Inf / Inf, so the
# exact Hessian of a model with a steep step in a state (a questionnaire's
# form change in age) failed, and the fit fell back to differences. So a
# Jacobian's function has every such quotient written with `inv_logit` before
# it is differentiated, and `.ctWithStableDerivRules()` differentiates
# `inv_logit` to a product of two and `log1p_exp` to one: finite at any `a`,
# given an `inv_logit` that is (the engine's, and Stan's, are). Text with
# nothing to rewrite is returned as it came, so numbers are not re-deparsed.
.ctStableLogistic <- function(text) {
  stripped <- function(e) {
    while (is.call(e) && identical(e[[1]], as.name("("))) e <- e[[2]]
    e
  }
  # `a` when `d` is `1 + exp(a)` or `exp(a) + 1`, else NULL.
  exponent <- function(d) {
    d <- stripped(d)
    if (!is.call(d) || !identical(d[[1]], as.name("+")) || length(d) != 3L) return(NULL)
    for (k in 2:3) {
      one <- d[[5L - k]]
      e <- stripped(d[[k]])
      if (is.numeric(one) && identical(as.numeric(one), 1) && is.call(e) &&
          identical(e[[1]], as.name("exp")) && length(e) == 2L) return(e[[2]])
    }
    NULL
  }
  changed <- FALSE
  walk <- function(e) {
    if (!is.call(e)) return(e)
    for (i in seq_along(e)[-1L]) e[[i]] <- walk(e[[i]])
    if (!identical(e[[1]], as.name("/")) || length(e) != 3L) return(e)
    a <- exponent(e[[3]])
    if (is.null(a)) return(e)
    changed <<- TRUE
    x <- stripped(e[[2]])
    if (is.call(x) && identical(x[[1]], as.name("exp")) && length(x) == 2L &&
        identical(x[[2]], a)) return(call("inv_logit", a))
    step <- call("inv_logit", call("-", call("(", a)))
    if (is.numeric(x) && identical(as.numeric(x), 1)) return(step)
    call("*", call("(", e[[2]]), step)
  }
  parsed <- walk(str2lang(text))
  if (!changed) return(text)
  # Fifteen digits unless a number needs more to come back exactly.
  written <- function(control) gsub(" ", "", paste(deparse(parsed,
    width.cutoff = 500L, control = control), collapse = ""), fixed = TRUE)
  out <- written(c("keepNA", "keepInteger", "niceNames", "showAttributes"))
  if (identical(str2lang(out), parsed)) out else written("digits17")
}

# `expr` with Deriv differentiating `inv_logit` and `log1p_exp` by rules that
# do not overflow (see `.ctStableLogistic()`), rather than through their bodies.
# Deriv's rule table is its own and other code may use it, so the rules are
# added for this call and the table put back as it was.
.ctWithStableDerivRules <- function(expr) {
  rules <- Deriv::drule
  added <- c("inv_logit", "log1p_exp")
  previous <- mget(added, envir = rules, ifnotfound = list(NULL), inherits = FALSE)
  on.exit(for (name in added) {
    if (is.null(previous[[name]])) {
      if (exists(name, envir = rules, inherits = FALSE)) rm(list = name, envir = rules)
    } else assign(name, previous[[name]], envir = rules)
  }, add = TRUE)
  assign("inv_logit", alist(x = inv_logit(x) * inv_logit(-x)), envir = rules)
  assign("log1p_exp", alist(x = inv_logit(x)), envir = rules)
  expr
}

ctJacobian <- function(m,types=c('J0','JAx','Jtd','Jy'),simplify=TRUE ){
  
  # get system dimension
  ndim = max(m$pars$row[m$pars$matrix%in% 'T0MEANS'])
  # 2): generate vector valued function fn = drift * state
  
  # initialize fn and state
  fn     = c()
  state   = paste0("state__", 1:ndim,'__')
  
  
  # #replace basic pars with in place system matrix references
  # for(ri in 1:nrow(m$pars)){
  #   if(grepl('[',m$pars$param[ri],fixed=TRUE) && !is.na(m$pars$transform[ri])){
  #     # m$pars$param[ri] <- gsub('param',m$pars$param[ri],m$pars$transform[ri])
  #   } else if(!grepl('[',m$pars$param[ri],fixed=TRUE) && !is.na(m$pars$param[ri])) m$pars$param[ri] <- paste0(m$pars$matrix[ri],'[',m$pars$row[ri],',',m$pars$col[ri],']')
  # }
  
  mats <- listOfMatrices(m$pars)
  matnames <- names(.ctMatricesList(unsafe=TRUE)$base)
  
  mats <- unfoldmats(mats)
  
  Jout <- list()
  for(typei in types){
    if(typei=='JAx'){
      Jrows = nrow(mats$T0MEANS)
      
      if(nrow(mats$DRIFT)!=ndim){ #append extra rows and columns to drift in case of intoverpop
        mats$DRIFT=rbind(
          cbind(mats$DRIFT,
            matrix(0,nrow(mats$DRIFT),ndim-nrow(mats$DRIFT))),
          matrix(0,ndim-nrow(mats$DRIFT),ndim))
        
        if(!m$continuoustime) diag(mats$DRIFT)[((nrow(mats$CINT)+1):ndim)] <- 1 #discrete time stability is 1 not 0
      }
      
      
      
      for (row in 1:ndim) {
        for (col in 1:(ndim)) {
          fn[row] = paste0(
            ifelse(col > 1, paste0(fn[row],' + '),''), 
            "(", mats$DRIFT[row, col], ") * state[", as.character(col), "]")
        }
        # browser()
        # if(!m$continuoustime && is.na(mats$DRIFT[row,row])) fn[row] <- paste0("state[", as.character(row), "] +",fn[row])
        if(!is.na(mats$CINT[row])) fn[row] = paste0(fn[row],' + ',mats$CINT[row]) #checking for NA because CINT is not always as large as DRIFT
      }
      # browser()
    }
    
    
    if(typei=='Jtd'){
      if(m$n.TDpred ==0) {
        Jout[[typei]] <- matrix(NA,0,0)
        next
      }
      Jrows = nrow(mats$T0MEANS)
      tdrows=nrow(mats$TDPREDEFFECT)
      fn = paste0('state[',1:ndim,']')
      fn[1:tdrows] = paste0(fn[1:tdrows],' + ', prodSymb(mats$TDPREDEFFECT, cbind(rep(1,m$n.TDpred))))
    }
    
    if(typei=='J0') {
      Jrows = nrow(mats$T0MEANS)
      t0func <- mats$T0MEANS[,1]
      t0func <- sapply(1:length(t0func), function(xi){
        out <-  paste0('state[',xi,'] + ',t0func[xi])
        return(out)
      })
      fn = prodSymb(diag(nrow(mats$T0MEANS)), matrix(t0func,ncol=1))
    }
    if(typei=='Jy') {
      Jrows = nrow(mats$MANIFESTMEANS)
      Jybase <- mats$LAMBDA
      Jybase <- cbind(Jybase,matrix(0,nrow(Jybase),ncol=nrow(mats$T0MEANS)-ncol(Jybase)))
      fn = paste0(mats$MANIFESTMEANS,' + ',prodSymb(Jybase,cbind(paste0('state[',1:ndim,']'))))
    }
    
    
    fn = gsub(" ", "", fn, fixed = TRUE) #remove spaces
    # replace state[~] by state~ for cOde Jacobian and make fn and state a named list
    names(fn) = paste0("fn", 1:length(fn))
    for(statei in 1:ndim){
      fn=gsub(paste0('\\b(state)\\[(',statei,')?\\]'),paste0('state__',statei,'__'),fn,perl = TRUE)
    }
    
    
    
    fn=gsub(' ','',fn,fixed=TRUE)
    #probably redundant now but maybe useful at some point?
    # replace remaining commas and square brackets for cOde Jacobian
    fn = gsub(",", "___comma___", fn, fixed = TRUE)
    fn = gsub("[", "___leftsquarebracket___", fn, fixed = TRUE)
    fn = gsub("]", "___rightsquarebracket___", fn, fixed = TRUE)
    
    # 3): calculate Jacobian of fn symbolically, with logistics in the form
    # whose derivative stays finite (`.ctStableLogistic()`).
    fn = vapply(fn, .ctStableLogistic, character(1))
    J = .ctWithStableDerivRules(sapply(state,function(statei){ #use Deriv package to calculate Jacobian to handle new functions
        sapply(fn, function(fni){
        Deriv::Deriv(fni, statei,cache.exp = FALSE)
      })}))
    # browser()
    # 4): create Jacobian list in STAN format
    # J = sapply(J,Simplify)
    
    J = gsub(" ", "", J, fixed = TRUE) #remove spaces
    if(simplify) J=sapply(J,Simplify) else { #remove wrapping brackets only
      J=sapply(J,function(x){
        x <- gsub('^\\((\\d+)\\)$','\\1',x)
        return(x)
      })
    }
    
    for(statei in 1:ndim){
      J=gsub(paste0('state__',statei,'__'),paste0('state[',statei,']'),J,fixed = TRUE)
    }
    
    # restore commas and square brackets for Jacobian list
    J = gsub("___rightsquarebracket___", "]", J, fixed = TRUE)
    J = gsub("___leftsquarebracket___", "[", J, fixed = TRUE)
    J = gsub("___comma___", ",", J, fixed = TRUE)
    
    
    
    Jm <- matrix(J,Jrows,ndim)
    
    #replace references to simple parameters with the matrix ref, when they are found in complex parameters
    simplepars <- unique(m$pars$param[!is.na(m$pars$param) & !grepl('\\W',m$pars$param)])
    if(length(simplepars)>0){
      for(pi in 1:length(Jm)){
        if(is.na(suppressWarnings(as.numeric(Jm[pi]))) &
            grepl('\\W', Jm[pi])
        ){
          for(spi in 1:length(simplepars)){
            modelrowmatch <- which(m$pars$param %in% simplepars[spi])[1]
            Jm[pi]<-gsub(paste0('\\<',simplepars[spi],'\\>'), 
              paste0(m$pars$matrix[modelrowmatch],'[',m$pars$row[modelrowmatch],',',m$pars$col[modelrowmatch],']'),
              Jm[pi])
          }}
      }
    }

    Jout[[typei]] <- Jm
  }#end type loop
  return(Jout)
}

