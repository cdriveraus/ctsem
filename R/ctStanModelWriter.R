#to do
# don't need statedep, use whenmat?


ctModel0DRIFT <- function(ctm,continuoustime){
  repval <- ifelse(continuoustime,-1e-6, 1-1e-6)
  checkval <- ifelse(continuoustime,0, 1)
  if(any(ctm$pars$value %in% checkval & ctm$pars$matrix %in% 'DRIFT' & ctm$pars$row == ctm$pars$col)) {
    message('Setting DRIFT diagonals of ',checkval,' to approximate ', checkval,' of ',repval)
    ctm$pars$value[ctm$pars$value %in% checkval & ctm$pars$matrix %in% 'DRIFT' & ctm$pars$row == ctm$pars$col] <- repval
  }
  return(ctm)
}

simpleStateCheck <- function(x){   #checks if system matrix elements that reference states are 'simple' -- don't contain operations
  if(grepl('\\b(state)\\b\\[\\d+\\]',x)){ #if state based 
    refmatches <- gregexpr('\\[(\\d+\\W*\\d*)\\]', x) #within sq brackets, find digits, possibly with non word chars, possibly with more digits.
    simplestate <- refmatches[[1]][1] > 0 && #if state
      length(unique(regmatches(x,refmatches)[[1]])) == 1 #reference 1 state only (possibly multiple times)
  } else simplestate <- FALSE
  return(simplestate) 
}


#' Update an already compiled and fit ctStanFit object
#' 
#' Allows one to change data and or model elements that don't require recompiling, then re fit.
#'
#' @param fit ctStanFit object
#' @param datalong data as normally passed to \code{\link{ctStanFit}}
#' @param ctstanmodel model as normally passed to \code{\link{ctStanFit}}
#' @param ... extra args for \code{\link{ctStanFit}}


ctStanUpdModel <- function(fit, datalong, ctstanmodel,...){
  
  new <-ctStanFit(datalong = datalong, ctstanmodel = ctstanmodel,fit=FALSE,...)
  
  fit$standata <- new$standata
  fit$data <- new$data
  fit$setup <- new$setup
  fit$args <- match.call()
  return(fit)
}



ctModelStatesAndPARS <- function(ctspec, statenames){ #replace latentName and par references with square bracket refs
  #detect state refs
  # 
  ln <- statenames
  for(li in c(1:length(ln))){
    for(ri in grep(paste0('\\b(',ln[li],')\\b'),ctspec$param)){
      ctspec$param[ri] <- gsub(paste0('\\b(',ln[li],')\\b'),paste0('state[',li,']'),ctspec$param[ri])
    }
  }
  #expand pars
  
  ln <- ctspec$param[ctspec$matrix %in% 'PARS' & !is.na(ctspec$param)] #get extra pars
  
  for(li in seq_along(ln)){ #for every extra par
    parmatch <- which(ctspec$param %in% ln[li] & ctspec$matrix %in% 'PARS') #which row is the par itself
    # if(grepl('taux',ln[li])) browser()
    for(ri in grep(paste0('\\b',ln[li],'\\b'),ctspec$param)){ #which rows contain the par
      if(!(ctspec$param[ri] == ln[li] & ctspec$matrix[ri]=='PARS')){ #that are not the par itself in the pars matrix# #removed limitation of referencing within PARS matrices
        # print(ctspec$param[ri])
        # print(ctspec$matrix[ri])
        # print(gsub(paste0('\\b',ln[li],'\\b([^\\[])'), #replace with PARS reference...
        # paste0('PARS[',ctspec$row[parmatch],',',ctspec$col[parmatch],']\\1'),ctspec$param[ri]))
        # browser()
        ctspec$param[ri] <- 
          gsub(paste0('\\b',ln[li],'\\b([^\\[])'), #replace with PARS reference...
            paste0('PARS[',ctspec$row[parmatch],',',ctspec$col[parmatch],']\\1'),ctspec$param[ri])
      }
    }
  }
  return(ctspec)
}

ctModelTransformsToNum<-function(ctm){
  
  fit.eqs = function(e) {
    # print(e)
    # List the types of formulas we might encounter.
    types=c(tformshapes(singletext = TRUE))
    formula.types = data.frame(
      type =0:(length(types)-1),
      formula = types,
      offset = 0,
      inneroffset = -.3, #c(0, rep(0,length(types)-1)),
      multiplier = 1,
      meanscale = 1, #c(1, rep(1,length(types)-1)),
      lsfit = NA,
      stringsAsFactors = FALSE
    )
    formula.types[formula.types$type==0,c('inneroffset','meanscale')] <- NA
    # Get some param values over multiplier wide range, and compute the corresponding y
    # values.
    
    param = c(seq(-2, 2, .1),seq(-10,10,.5),c(rnorm(10)))
    y = eval(parse(text=e)) #eval(substitute(substitute(e, list(param = param)), list(e = as.quoted(e)[[1]]))))
    keep <- abs(y) < 1e5 & !is.na(param)
    keep[is.na(keep)] <- FALSE
    param <- param[keep]
    y <- y[keep]
    # plot(param,y)
    
    # Try to fit each formula to the data.
    success <- FALSE
    for(tryi in 1:2){ #if no success try with more enthusiasm
      if(success) break
      for(i in 1:nrow(formula.types)) {
        if(i > 1 && any((formula.types$lsfit < .001) %in% TRUE)) next #found the answer already
        tftype <- formula.types$type[i]
        
        ff <- function(pars){
          multiplier=pars[1];
          offset=pars[2];
          if(tftype > 0) {
            meanscale=pars[3];
            inneroffset=pars[4]
          } else{
            meanscale=1;inneroffset=0
          }
          yest<- suppressWarnings(eval(parse(text=formula.types$formula[i])))
          res <- (sum( ((y-yest)^2)/(abs(y)+.01)))
          if(is.na(res)) res <- 1e100
          return(res)
        }
        
        ffg <- function(pars){
          b=ff(pars)
          g=try(sapply(1:length(pars),function(x) {
            parx=pars
            parx[x]<-pars[x]+1e-12;
            (b-ff(parx))/1e-12
          }))
          if('try-error' %in% class(g)) g <- rnorm(pars)
          if(any(is.na(g))) g[is.na(g)] <- rnorm(sum(is.na(g)))
          return(-g)
        }
        
        if(tryi==1){
          etmp<-gsub('\\(',' ',e)
          etmp<-gsub('\\)',' ',etmp)
          etmp<- gsub('([^0-9])\\.','\\10.',etmp) #add a zero in front of decimal only numbers
          escn <- unlist(regmatches(etmp, gregexpr('[0-9,.]+e+-?[0-9\\b[0-9]+', etmp))) #replace scientific notation
          escn<- as.numeric(gsub('\\(', '-', gsub(',', '', escn)))
          etmp <- gsub('\\b[0-9,.]+e+-?[0-9\\b[0-9]+','',etmp)
          etmp <- gsub("\\[.*\\]","",etmp)
          samplestarts <- unlist(regmatches(etmp, gregexpr('\\b[0-9,.]+\\b', etmp)))
          samplestarts <- as.numeric(gsub('\\(', '-', gsub(',', '', samplestarts)))
          samplestarts <- c(samplestarts,escn)
          samplestarts <- unique(c(samplestarts,0,1))#,rep(0,max(0,4-length(samplestarts))))
          samplestarts <- unique(c(samplestarts,-samplestarts))
          teststarts=expand.grid(rep(list(samplestarts), 4))
          if(is.na(formula.types$meanscale[i])) teststarts[,3] <- 1
          if(is.na(formula.types$inneroffset[i])) teststarts[,4] <- 0
          teststarts <- unique(teststarts)
        } else{
          
          gridbase=ifelse(tryi==1,1.5,.5)
          gs=seq(-2,2,gridbase)
          p1g <- p2g <- sort(c(-exp(gs),exp(gs)))
          if(tftype > 0){
            p3g <- p4g <- p1g
          } else {
            p3g <-1
            p4g <-0
          }
          min <- Inf
          
          teststarts <- expand.grid(p1g,p2g,p3g,p4g)
          
          teststarts <- teststarts[teststarts[,3]!=0,]
          teststarts <- teststarts[teststarts[,1]!=0,]
          teststarts <- teststarts[!duplicated(teststarts),]
        } #end try 2 setup
        # 
        res <- apply(teststarts,1,function(x) ff(unlist(x)))
        start.params <- teststarts[which(res == min(res))[1],]
        fit <- list()
        fit$par <- as.numeric(c(start.params))
        fit$f<-res[which(res == min(res))[1]]
        
        
        formula.types$offset[i] = fit$par[2] #round(coef(fit)[["offset"]])
        formula.types$multiplier[i] = fit$par[1] #round(coef(fit)[["multiplier"]])
        
        if(formula.types$type[i] >0) {
          formula.types$meanscale[i] = fit$par[3] #round(coef(fit)[["meanscale"]])
          formula.types$inneroffset[i] = fit$par[4] #round(coef(fit)[["inneroffset"]])
        } else {
          formula.types$meanscale[i] =1
          formula.types$inneroffset[i]=0
        }
        formula.types$lsfit[i] = fit$f #AIC(fit)
      }
      # 
      # if(any(formula.types$lsfit < 1e-6)) success <- TRUE else {
      #   
      #   message('Trying to determine transforms...')#
      # }
      success<-TRUE
    }
    
    return(formula.types[which(formula.types$lsfit %in% min(formula.types$lsfit,na.rm=TRUE))[1],])
  }
  # 
  rownames(ctm$pars)=1:nrow(ctm$pars)
  newrows <- which(!is.na(ctm$pars$transform))
  if(length(newrows) ==0){
    return(ctm)
  } else{
    eqs <- ctm$pars$transform[newrows]
    uniqueeqs <- unique(eqs)
    l=lapply(uniqueeqs,fit.eqs)
    eqmatch <- match(eqs, uniqueeqs)
    df=data.frame(do.call(rbind,l))
    df[1:length(eqmatch),] <- df[eqmatch,]
    df[,] <- lapply(df,function(x) if(is.numeric(x)) return(round(x,6)) else return(x))
    df$formula <- eqs
    df[df$lsfit > .1,c('offset','inneroffset','multiplier','meanscale')]<-NA
    colnames(df)[1] <- 'transform'
    df$transform[df$lsfit > .1] <- eqs[df$lsfit > .1]
    df$lsfit <- NULL
    
    #hack for NA's on unused pars
    df$inneroffset[!is.na(df$offset) & is.na(df$inneroffset)] <- 0
    df$meanscale[!is.na(df$offset) & is.na(df$meanscale)] <- 1
    
    rownames(df) <- newrows
    ctm$pars$transform <- NULL
    
    df <- data.table(Row.Names=as.integer(rownames(df)),df)
    nctspec <- data.table(Row.Names=1:nrow(ctm$pars),ctm$pars)
    nctspec <- data.frame(merge(nctspec,df,all=TRUE,no.dups = FALSE))
    
    nctspec <- nctspec[,c('matrix','row','col','param','value','transform','multiplier',
      'offset','meanscale','inneroffset','indvarying','sdscale',
      colnames(ctm$pars)[grep('_effect',colnames(ctm$pars),fixed=TRUE)]) ]
    # 
    if(any(rl(suppressWarnings(as.numeric(nctspec$transform)) >= 6))){ #if any are jacobian calcs length(tformshapes(singletext = TRUE))))) {
      nctspec$transform[rl(suppressWarnings(as.numeric(nctspec$transform)) >= 6)] <- #adjust jacobian gradients
        suppressWarnings(as.numeric(
          nctspec$transform[rl(suppressWarnings(as.numeric(nctspec$transform)) >= 6)])) + 
        (50-6)
    }
    ctm$pars <- nctspec
    
    return(ctm)
  }
  
}


ctStanModelIntOverPop <- function(m){ 
  if(sum(m$pars$indvarying) < 1) {
    message('No individual variation for ctStanModelIntOverPop to work with!')
    return(m)
  } else {
    m$pars <- ctStanModelCleanctspec(m$pars)
    t0mvaryingsimple <- m$pars$row[m$pars$indvarying & m$pars$matrix %in% 'T0MEANS'] #which t0means are indvarying
    t0mvaryingnames <- m$pars$param[m$pars$indvarying & m$pars$matrix %in% 'T0MEANS'] #names of t0means that are indvarying
    t0mnotvarying <- m$pars$row[!m$pars$indvarying & m$pars$matrix %in% 'T0MEANS']
    
    
    ivnames <- unique(m$pars$param[m$pars$indvarying & !m$pars$param %in% t0mvaryingnames]) #don't need new states for t0means
    m$latentPopNames=ivnames
    ivnamesfull <- c(t0mvaryingnames,ivnames) #for t0var naming
    nindvaryingsmall <- length(ivnames)
    
    if(length(ivnames) > 0){
      #new t0means
      t0m <- m$pars[match(ivnames, m$pars$param),,drop=FALSE]
      t0m$matrix <- 'T0MEANS'
      t0m$col <- 1
      t0m$row <- (m$n.latent+1):(m$n.latent+nindvaryingsmall)
      t0m$transform <- 'param'
      t0m$indvarying <- TRUE
      
      #new t0var 
      t0v <- m$pars[m$pars$matrix %in% 'T0VAR' & m$pars$row==1 & m$pars$col==1,,drop=FALSE]
      for(ri in 1:(m$n.latent+nindvaryingsmall)){
        for(ci in 1:(m$n.latent+nindvaryingsmall)){
          if(!(ri %in% t0mnotvarying && ci %in% t0mnotvarying)){
            t0v <- rbind(t0v,c('T0VAR',ri,ci,
              NA,  #param
              0,  #value
              NA, #transform,
              FALSE,1,rep(FALSE,m$n.TIpred)))
            m$pars <- m$pars[!(m$pars$matrix %in% 'T0VAR' & m$pars$row==ri & m$pars$col==ci),,drop=FALSE] #remove old t0var line
          }
        }}
      t0v=t0v[-1,,drop=FALSE] #remove initialisation row
      
      
      
      #reference new states
      for(ivi in ivnames){
        m$pars$indvarying[m$pars$param %in% ivi] <- FALSE
        m$pars[m$pars$param %in% ivi,paste0(m$TIpredNames,rep('_effect',m$n.TIpred))] <- FALSE
        m$pars$param[m$pars$param %in% ivi] <- sapply(which(m$pars$param %in% ivi), function(ri){
          gsub('param',paste0( 'state[',m$n.latent+match(ivi,ivnames),']'),m$pars$transform[ri]) 
        })
        m$pars$transform[m$pars$param %in% ivi] <-NA
        
      }
      
      m$pars <- rbind(m$pars, t0m,t0v)#,drift) #
      m$pars[] <- lapply(m$pars, utils::type.convert, as.is = TRUE)
      
    } #finish loop for non simple t0means indvarying
    extralatents <- seq.int(m$n.latent+1,m$n.latent+length(ivnames),length.out=length(ivnames))
    
    m$intoverpopindvaryingindex <- c(t0mvaryingsimple,extralatents)
    return(m)
  }
}




simplifystanfunction<-function(bcalc,simplify=TRUE){ #input text of list of computations, output simplified form
  if(length(bcalc)==0 || !grepl('=',bcalc)) return('')
  bcalc=gsub(' ','',bcalc)
  if(simplify && nchar(bcalc) > 200 && exists('verbose') && as.logical(verbose) ==TRUE ){
    cat(bcalc)
    message('Trying to simplify the additional calculations...')
  }
  
  bcalc=gsub('\n','',bcalc)
  bcalcs=strsplit(bcalc,';')[[1]]
  bcalcs = bcalcs[!nchar(bcalcs)==0]
  bcalcs1=gsub('=.*','',bcalcs)
  bcalcs2=gsub('.*=','',bcalcs)
  bcalcs2 = gsub(",", "___comma___", bcalcs2, fixed = TRUE)
  bcalcs2 = gsub("[", "___leftsquarebracket___", bcalcs2, fixed = TRUE)
  bcalcs2 = gsub("]", "___rightsquarebracket___", bcalcs2, fixed = TRUE)
  
  bcalcs2c=paste0('c(',paste0(bcalcs2,collapse=', '),')')
  if(simplify) scalcs2=Deriv::Deriv(bcalcs2c,nderiv=0)  else scalcs2=bcalcs2c 
  if(simplify && !grepl('\\.e',scalcs2)) scalcs2=bcalcs2c #hack to avoid weird leftovers from non reduced simplify
  
  if(scalcs2 == 'NULL') return('') else{
    
    names(scalcs2)=NULL
    scalcs2=gsub(' ','',scalcs2)
    
    scalcs2 = gsub(",", "__eqbreak__=", scalcs2, fixed = TRUE)
    scalcs2 = gsub("___rightsquarebracket___", "]", scalcs2, fixed = TRUE)
    scalcs2 = gsub("___leftsquarebracket___", "[", scalcs2, fixed = TRUE)
    scalcs2 = gsub("___comma___", ",", scalcs2, fixed = TRUE)
    
    scalcs2=gsub('\\.e(\\d+)<-','real e\\1 =',scalcs2)
    scalcs2=gsub('\\.e(\\d+)','e\\1',scalcs2)
    scalcs2=gsub('\\{','',scalcs2) #remove leading curly brace
    
    ec=gsub('\\b(c)\\b\\(.*','',scalcs2) #retain just the simplifications (both lhs and rhs)
    scalcs2=gsub('.*\\b(c)\\b\\(','=',scalcs2) #retain just the rhs equations for the pars, insert first =
    scalcs2=gsub(';','',scalcs2)
    
    scalcs2=gsub('\\}$','',scalcs2) #remove trailing curly ;
    scalcs2=gsub('\\)$',';',scalcs2) #replace trailing bracket close with ;
    # scalcs2=
    #   gsub('\\)*\\){1}$',';',scalcs2) #replace trailing bracket close with ;
    
    scalcs2 = strsplit(scalcs2,'__eqbreak__')[[1]]
    scalcs = paste0(bcalcs1,scalcs2,';\n')
    
    out = paste0('    {\n  ',paste0(ec,collapse='\n  '),'\n',paste0(scalcs,collapse=' '),'  } \n  ',collapse=' ')
    # cat(out)
    return(out)
  }
}



ctStanModelCleanctspec <-  function(ctspec){ #clean ctspec structure, non numeric transform style
  tieffects <- colnames(ctspec)[grep('_effect',colnames(ctspec),fixed=TRUE)]
  found=FALSE
  ctspec$indvarying=as.logical(ctspec$indvarying)
  if(any(ctspec$indvarying[ctspec$matrix %in% 'T0VAR'])){
    ctspec$indvarying[ctspec$matrix %in% 'T0VAR'] <- FALSE
    message('Individual variation in T0VAR parameters not possible, removed')
  }
  
  ctspec$row <- as.integer(ctspec$row)
  ctspec$col <- as.integer(ctspec$col)
  ctspec$value=as.numeric(ctspec$value)
  # ctspec$matrix=as.integer(ctspec$matrix)
  ctspec$indvarying = as.logical(ctspec$indvarying)
  ctspec$transform=as.character(ctspec$transform) #ensure consistent type in case custom tforms used
  ctspec$param=gsub(' ','',as.character(ctspec$param)) #remove spaces
  
  
  #values imply nothing happening elsewhere
  fixed <- !is.na(ctspec$value) 
  calc <-  grepl('\\+|\\*|/|-|\\[|\\]', ctspec$param)
  ctspec$indvarying[ fixed | calc ] <- FALSE #remove indvarying for calcs also
  ctspec$param[ fixed ] <- NA
  ctspec$transform[ fixed] <- NA
  if(length(tieffects) > 0)  ctspec[fixed | calc,tieffects] <- FALSE
  if(any(apply(ctspec[,c('value','param')],1,function(x) all(is.na(x))))) stop('Parameters specified as NA ! Needs a value or character label.')
  
  
  #was slow!
  # for(rowi in 1:nrow(ctspec)){ 
  #   if( !is.na(ctspec$value[rowi])) {#fixed value replacement
  #     if(any(c(!is.na(ctspec[rowi,'param']),!is.na(ctspec[rowi,'transform']),ctspec[rowi,'indvarying']))){
  #       if(ctspec[rowi,'value']!=99999) found<-TRUE
  #       ctspec[rowi,c('param','transform','indvarying')]=replacement
  #     }
  #   }
  #   if(grepl('[', ctspec$param[rowi],fixed=TRUE) || !is.na(ctspec$value[rowi])){
  #     if(ctspec$indvarying[rowi]){
  #       # message('Individual variation requested on deterministic parameter ', ctspec$param[rowi],' , setting to FALSE')
  #       ctspec$indvarying[rowi] <- FALSE
  #     }
  #     if((length(tieffects) > 0 && any(as.logical(ctspec[rowi,tieffects]) %in% TRUE)) || !is.na(ctspec$value[rowi])){
  #       # message('TI predictor effects requested on deterministic parameter ', ctspec$param[rowi],' , setting to FALSE')
  #       ctspec[rowi,tieffects] <- FALSE
  #     }
  #   } #end deterministic relations check
  #   if(all(is.na(c(ctspec[rowi,c('value','param')])))) stop('Parameters specified as NA ! Needs a value or character label.')
  # } #end row loop
  # if(found) message('Minor inconsistencies in model found - removing param name, transform and indvarying from any parameters with a value specified')
  return(ctspec)
}

ctStanMatricesList <- function(unsafe=FALSE){
  base <- c(PARS=10, T0MEANS=1,LAMBDA=2,DRIFT=3,DIFFUSION=4,MANIFESTVAR=5,MANIFESTMEANS=6, CINT=7,
    T0VAR=8,TDPREDEFFECT=9)
  jacobian = c(JAx=52,Jtd=53,Jy=54) #J0=51,
  asymptotic = c(asymCINT=21,asymDIFFUSIONcov=22)
  extra <- c(DIFFUSIONcov=31,MANIFESTcov=32,T0cov=33)
  all <- c(base,jacobian,asymptotic, extra)
  mn <- list(base=base, jacobian=jacobian, asymptotic=asymptotic, extra=extra,all=all)
  mn$driftcint <- all[names(all) %in% c('DRIFT','CINT')]
  mn$diffusion <- all[names(all) %in% c('DIFFUSION','JAx')]
  mn$tdpred <- all[names(all) %in% c('TDPREDEFFECT','Jtd')]
  mn$measurement <- all[names(all) %in% c('LAMBDA','MANIFESTMEANS','MANIFESTVAR','Jy')]
  mn$t0 <- all[names(all) %in% c('T0MEANS','T0VAR')]
  return(mn)
}


ctStanModelMatrices <-function(ctm){
  mats <- ctStanMatricesList(unsafe=TRUE)
  ctspec <- ctm$pars
  n.TIpred <- ctm$n.TIpred
  matsetup <-list()
  matvalues <-list()
  freepar <- 0
  freeparcounter <- 0
  
  indvaryingcounter <- 0
  TIPREDEFFECTsetup <- matrix(0,0,n.TIpred)
  tipredcounter <- 1
  indvar <- 0
  extratformcounter <- 0
  extratforms <- c()
  calcs<-ctm$calcs
  matsetup<-NULL
  mval<-NULL
  
  
  matlist <- listOfMatrices(ctspec) #put back into matrix form
  # matlist <- unfoldmats(matlist) #unfold references to basic state
  ctspecp <- ctModelUnlist(matlist,matnames = names(c(mats$base,mats$jacobian))) #convert back to list
  ctspecp <- ctModelStatesAndPARS(ctspecp,statenames=ctm$latentNames) #ensure PARS[,] refs are not replaced by names
  ctspec[order(ctspec$matrix, ctspec$row,ctspec$col),colnames(ctspecp)] <-  
    ctspecp[order(ctspecp$matrix, ctspecp$row,ctspecp$col),] #replace updated columns, free of direct references
  
  
  for(i in 1:nrow(ctspec)){ 
    nm = ctspec$matrix[i]
    m = mats$all[names(mats$all) %in% nm]
    when = 0 #default, static parameter
    parameter = 0 #for fixed values
    dynpar=FALSE
    copymatrix=0
    copyrow=0
    copycol=0
    simplestate <- FALSE
    tipred <- 0L
    if(!is.na(ctspec$param[i])){ #if non fixed parameter,
      
      if(grepl('\\b(state)\\b\\[\\d+\\]',ctspec$param[i])){ #if state based, 
        simplestate <- simpleStateCheck(ctspec$param[i]) #check if only 1 state we can easily handle
        
        if(simplestate){
          
          cpx <- ctspec[i,]
          cpx$transform <- gsub('\\b(state)\\b\\[\\d+\\]','param',cpx$param)
          cpx$multiplier <- cpx$offset <- cpx$meanscale <- cpx$inneroffset <- NULL
          ctspec[i,] <- ctModelTransformsToNum(list(pars=cpx))$pars
          if(is.na(suppressWarnings(as.integer(ctspec$transform[i])))) simplestate <- FALSE  #if couldn't convert to integer tform
        } #if any problems set simplestate FALSE
        
        if(simplestate){
          indvar <- 0
          parameter <- gsub('\\[|\\]','',unlist(regmatches(ctspec$param[i], #extract one or more references
            gregexpr(
              paste0('(?=\\[).*?(?<=\\])'),
              ctspec$param[i], perl=TRUE)
          )))[1]
          if(nm %in% c('JAx',names(c(mats$driftcint,mats$diffusion)),'PARS')) when = 2 #included PARS here because they can't change due to dynamics
          if(nm %in% c('PARS')) when = 100 #because needs to be expanded to all whens when computing whenvecs -- find a better way of handling this
          if(nm %in% c('Jtd',names(mats$tdpred))) when = 3
          if(nm %in% c('Jy',names(mats$measurement))) when = 4
          # if(nm %in% c('J0',names(mats$t0))) when = 1
          dynpar=TRUE
        }
        
      }
      
      
      if(!simplestate && grepl('[',ctspec$param[i],fixed=TRUE)){ #if a complex calc par add to list of calcs to be added to stan program
        calcs <- c(calcs, paste0(ctspec$matrix[i],'[',ctspec$row[i], ', ', ctspec$col[i],'] = ',
          ctspec$param[i]))
        when= -999 #never use this row of ctspec during calculations (calc is extracted and placed elsewhere)
        indvar=0
        dynpar=TRUE
      }#end calculation / copy parameters
      
      else if(!grepl('[',ctspec$param[i],fixed=TRUE)){ #if non calculation parameter (consider removing duplication / copyrow checks)
        
        if(i > 1 && any(ctspec$param[1:(i-1)] %in% ctspec$param[i])){ #and after row 1, check for duplication
          
          rowmatch <- match(ctspec$param[i], matsetup$parname)#find which freepar corresponds to duplicate
          parameter <- matsetup$param[rowmatch]
          indvar <- matsetup$indvarying[rowmatch] #ifelse(ctspec$indvarying[i],  matsetup[,'indvarying'][ match(ctspec$param[i], matsetup$parname) ],0)#and which indvar corresponds to duplicate
          tipred <- matsetup$tipred[rowmatch]
          
        } else { #if not duplicated
          TIPREDEFFECTsetup <- rbind(TIPREDEFFECTsetup,rep(0,ncol(TIPREDEFFECTsetup))) #add an extra row...
          
          if(ctspec$indvarying[i]) {
            indvaryingcounter <- indvaryingcounter + 1
            indvar <- indvaryingcounter
          } else  indvar <- 0
          
          freeparcounter <- freeparcounter + 1
          freepar <- freeparcounter
          parameter <- freepar
          
          if(is.na(suppressWarnings(as.integer(ctspec$transform[i])))) { #extra tform needed
            # stop('extra transforms presently disabled -- use PARS matrix if necessary')
            extratformcounter <- extratformcounter + 1
            extratforms <- paste0(extratforms,'if(transform==',-10-extratformcounter,') param = ',
              ctspec$transform[i],';')
            ctspec$transform[i] <- -10-extratformcounter
          }
          
          if(n.TIpred > 0) {
            TIPREDEFFECTsetup[freepar,][ ctspec[i,paste0(ctm$TIpredNames,'_effect')]==TRUE ] <- 
              tipredcounter: (tipredcounter + sum(as.integer(suppressWarnings(ctspec[i,paste0(ctm$TIpredNames,'_effect')]))) -1)
            tipredcounter<- tipredcounter + sum(as.integer(suppressWarnings(ctspec[i,paste0(ctm$TIpredNames,'_effect')])))
            tipred <- as.integer( any(TIPREDEFFECTsetup[freepar,] > 0))
          }
        }#end not duplicated loop
      } #end non calculation parameters
    } #end non fixed value loop
    # if(simplestate) 
    
    mdatnew <- data.frame(
      parname=ctspec$param[i],
      row=ctspec$row[i],
      col=ctspec$col[i],
      param=parameter, # ifelse(!is.na(ctspec$param[i]) && !grepl('[',ctspec$param[i],fixed=TRUE),freepar, 0), #freepar reference
      transform=ifelse(is.na(suppressWarnings(as.integer(ctspec$transform[i]))), -1, ctspec$transform[i]), #transform
      indvarying=ifelse(!is.na(ctspec$param[i]),indvar,0), #indvarying
      tipred=tipred, #ifelse(any(TIPREDEFFECTsetup[freepar,] > 0) && !dynpar, 1, 0), #tipredvarying
      matrix=m, #matrix reference
      when=when,#when to compute
      copymatrix=copymatrix,copyrow=copyrow,copycol=copycol,
      stringsAsFactors = FALSE)
    
    if(is.null(matsetup)) matsetup <- mdatnew else matsetup<-rbind(matsetup,mdatnew)
    mvalnew<-ctspec[i,c('value','multiplier','meanscale','offset','sdscale','inneroffset'),drop=FALSE]
    if(is.null(mval)) mval <- mvalnew else mval<-rbind(mval,mvalnew)
  }
  if(!is.null(mval)) mval[is.na(mval)] <- 99999 else mval <- array(0L,dim=c(0,6))
  
  
  
  matrixdims <- array(0L,dim=c(max(mats$all),2)) #get dims for each matrix
  sapply(mats$all, function(m) {
    matrixdims[m,]<<-as.integer(c(max(c(0,matsetup[matsetup[,'matrix'] %in% m,'row'])),
      max(c(0,matsetup[matsetup[,'matrix'] %in% m,'col']))))
  })
  
  #special matrix adjustments
  matrixdims[mats$all[names(mats$all ) == 'asymDIFFUSIONcov'],] <- 
    matrixdims[mats$all[names(mats$all ) == 'DIFFUSIONcov'],] <- 
    matrixdims[mats$all[names(mats$all ) == 'DIFFUSION'],]
  matrixdims[mats$all[names(mats$all ) == 'asymCINT'],] <- 
    matrixdims[mats$all[names(mats$all ) == 'CINT'],]
  matrixdims[mats$all[names(mats$all ) == 'MANIFESTcov'],] <- 
    matrixdims[mats$all[names(mats$all ) == 'MANIFESTVAR'],]
  matrixdims[mats$all[names(mats$all ) == 'T0cov'],] <- 
    matrixdims[mats$all[names(mats$all ) == 'T0VAR'],]
  
  #set integer types
  matsetup[,!colnames(matsetup) %in% 'parname'] <- lapply(matsetup[,!colnames(matsetup) %in% 'parname'],as.integer)
  #set numeric types
  matvalues <- data.frame(apply(mval,2,as.numeric,.drop=FALSE),stringsAsFactors = FALSE  )
  
  #copy handling -- find which matsetup row contains the reference
  for(i in 1:nrow(matsetup)){
    if(matsetup$copymatrix[i] > 0){
      matsetup$copyrow[i] <- which(matsetup$row == matsetup$copyrow[i] &
          matsetup$col == matsetup$copycol[i] &
          matsetup$matrix == matsetup$copymatrix[i])[1] #select first element in case of multiple matches
    }
  }
  
  for(i in 1:nrow(matsetup)){#check for nested copying and unnest
    if(matsetup$copyrow[i] > 0 && matsetup$copyrow[matsetup$copyrow[i]] >0){
      matsetup$copyrow[i] <- matsetup$copyrow[matsetup$copyrow[i]]
    }
  }
  
  for(i in 1:nrow(matsetup)){ #copy elements of reference row to copyrow
    if(matsetup$copyrow[i] > 0){
      matsetup[i,c('param','transform','indvarying','tipred','when')] <- 
        matsetup[matsetup$copyrow[i],c('param','transform','indvarying','tipred','when')]
    }
  }
  
  matsetup$copyrow[unique(matsetup$copyrow[matsetup$copyrow > 0])] <-
    -unique(matsetup$copyrow[matsetup$copyrow > 0]) #set rows to be copied from to neg of their own row num
  
  matsetup$copymatrix <-  matsetup$copycol <- NULL #not needed columns after processing
  
  #sort order so PARS are computed first (in case of state dependency in PARS, further dependencies on PARS must happen after)
  matvalues = rbind(matvalues[matsetup$matrix==10,],matvalues[matsetup$matrix!=10,])
  matsetup = rbind(matsetup[matsetup$matrix==10,],matsetup[matsetup$matrix!=10,])
  
  
  return(list(
    matsetup=matsetup,   matvalues=matvalues, 
    extratforms=extratforms,
    TIPREDEFFECTsetup=TIPREDEFFECTsetup,
    matrixdims=matrixdims,
    calcs=calcs))
}



ctStanCalcsList <- function(ctm, save=FALSE){  #extract any calcs from model into specific lists
  temp <- ctm$modelmats$calcs
  names(temp) <- NULL
  mats<-ctStanMatricesList()
  
  calcs <- lapply(c(list(PARS=c(PARS=10)),mats[!names(mats) %in% c('base','all','jacobian')]), function(mlist) { #add custom calculations to correct list of matrices
    out <- temp[unlist(sapply(temp, function(y) any(sapply(names(mlist), function(mli) 
      grepl(mli,gsub('=.*','',y))
    ))))]
    return(out)
  })
  
  if(save){ #if outputting nonlinear computations
    uniqueCounter <- 1
    uniqueNames <- c()
    calcs <- lapply(calcs, function(x){ #get unique equations and save these to output vector if save requested
      whicheqs <- which(!duplicated(gsub('^.* = ','',x)))
      if(length(x)){
        rhs <- gsub('^.* = ','',x)[whicheqs]
        lhs <- gsub(' = .*','',x)[whicheqs]
        neweqs <- paste0('if(si >0) calcs[rowi, ', uniqueCounter:(uniqueCounter+(length(lhs)-1)),'] = ',rhs)
        x <- c(x,neweqs)
        uniqueCounter <<- uniqueCounter+length(lhs)
        uniqueNames <<- c(uniqueNames,lhs)
      }
      return(x)
    })
    
    calcs$calcNames <- uniqueNames
  }
  
  
  ctm$modelmats$calcs <- calcs
  return(ctm)
}

ctStanModelWriter <- function(ctm, gendata, extratforms,matsetup,savemodel=TRUE, simplify=TRUE){
  #if arguments change make sure to change ctStanFit !
  
  mats <- ctStanMatricesList()
  # #check when / if PARS needs to be computed
  # for(mlist in names(mats[-1])){
  #   if(any(unlist(lapply(ctm$calcs[[mlist]], function(m) grepl('PARS',m))))) mats[[mlist]]=c(mats[[mlist]],'PARS')
  # }
  
  
  #adjust diffusion calcs to diffusionsqrt
  # ctm$calcs$diffusion <- gsub('DIFFUSION','DIFFUSIONsqrt',ctm$calcs$diffusion)
  # intoverpopdynamiccalcs <- gsub('DIFFUSION','DIFFUSIONsqrt',intoverpopdynamiccalcs)
  
  #save calcs without intoverpop for param intitialisation
  
  

  finiteJ<-function(){
    paste0('
    {
    int zeroint[1];
    row_vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(JAxfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0)  state[statei] += Jstep;
          
        //initialise PARS first, and simple PARS before complex PARS
        if(statedep[10] || whenmat[10,2]) PARS=mcalc(PARS,indparams, state,{2}, 10, matsetup, matvalues, si); 
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n\n ',collapse=' '),simplify),'
        if(statedep[3] || whenmat[3,2]) DRIFT=mcalc(DRIFT,indparams, state,{2}, 3, matsetup, matvalues, si); 
        if(statedep[7] || whenmat[7,2]) CINT=mcalc(CINT,indparams, state,{2}, 7, matsetup, matvalues, si); 
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$driftcint,';\n\n ',collapse=' '),simplify),'
      
      if(statei > 0) {
        JAx[1:nlatent,statei] =  DRIFT * state[1:nlatent]\' + CINT[,1]; //compute new change
         if(verbose>1) print("JAx ",JAx);
      }
      if(statei== 0 && size(JAxfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
        base = DRIFT * state[1:nlatent]\' + CINT[,1];
        if(verbose>1) print("base = ",base,"    sjaxinit= ",JAx);
        for(fi in JAxfinite){
          JAx[1:nlatent,fi] -= base;
          JAx[1:nlatent,fi] /= Jstep; //new - baseline change divided by stepsize
        }
      }
    }
    if(verbose>1) print("JAx ",JAx);
    }
    ')
  }
  
  
  
  
  
  
  
  ukfilterfunc<-function(ppchecking){
    out<-paste0('
  int prevrow=0;
  real prevdt=0;
  real dt=1; //initialise to make sure drift is computed on first pass
  real dtsmall;
  int dtchange=1;
  real prevtime=0;
  int T0check=0;
  matrix[nlatentpop, nlatentpop] etacov; //covariance of latent states

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] syprior;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypriorcov_sqrt = rep_matrix(0,nmanifest,nmanifest); 
  matrix[nmanifest, nmanifest] ycov; 
  
  matrix[nlatentpop,nlatentpop] eJAx = diag_matrix(rep_vector(1,nlatentpop)); //time evolved jacobian
  matrix[nlatentpop,nlatentpop] eJAxs[dosmoother ? ndatapoints : 1]; //time evolved jacobian, saved for smoother

  row_vector[nlatentpop] state = rep_row_vector(-999,nlatentpop); 
  matrix[nlatentpop,nlatentpop] JAx; //Jacobian for drift
  //matrix[nlatentpop,nlatentpop] J0; //Jacobian for t0
  matrix[nlatentpop,nlatentpop] Jtd;//diag_matrix(rep_vector(1),nlatentpop); //Jacobian for nltdpredeffect
  matrix[ nmanifest,nlatentpop] Jy;//Jacobian for measurement 
  matrix[ nmanifest,nlatentpop] Jys[dosmoother ? ndatapoints : 0];//saved Jacobian for measurement smoother
  
  ',if(!is.null(ctm$taylorheun) && ctm$taylorheun==1) paste0('
  matrix[taylorheun ? nlatentpop : 0, taylorheun ? nlatentpop : 0] Kth = rep_matrix(0.0,taylorheun ? nlatentpop : 0, taylorheun ? nlatentpop : 0); 
  matrix[taylorheun ? nlatentpop : 0, taylorheun ? nlatentpop : 0] Mth = Kth;'),'

  //linear continuous time calcs
  matrix[nlatent,nlatent] discreteDRIFT;
  vector[nlatent] discreteCINT;
  matrix[nlatent,nlatent] discreteDIFFUSION = rep_matrix(0.0,nlatent,nlatent);

  
  vector[nparams] rawindparams = rawpopmeans;
  vector[nparams] indparams;
  
  matrix[nlatentpop,nlatentpop] etacovb[3,dosmoother ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ycovb[3,dosmoother ? ndatapoints : 0];
  vector[nlatentpop] etab[3,dosmoother ? ndatapoints : 0];
  vector[nmanifest] yb[3,dosmoother ? ndatapoints : 0];

  //dynamic system matrices
  ',subjectparaminit(pop=FALSE,smats=TRUE),'
  
  asymDIFFUSIONcov = rep_matrix(0,nlatent,nlatent); //in case of derrindices need to init
  DIFFUSIONcov = rep_matrix(0,nlatent,nlatent);

  for(rowx in 0:(dokalman ? ndatapoints : 0)){
    int rowi = rowx ? rowx : 1;
    if( rowx==0 ||
      (dokalmanrows[rowi] && 
        subject[rowi] >= (firstsub - .1) &&  subject[rowi] <= (lastsub + .1))){ //if doing this row for this subject
    
    int si = rowx ? subject[rowi] : 0;
    int full = (dosmoother==1 || si ==0);
    int o[full ? nmanifest : nobs_y[rowi]]; //which obs are not missing in this row
    int o1[full ? size(whichequals(manifesttype,1,1)) : nbinary_y[rowi] ];
    int o0[full ? size(whichequals(manifesttype,1,0)) : ncont_y[rowi] ];
    
    int od[nobs_y[rowi]] = whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
    int o1d[nbinary_y[rowi] ]= whichbinary_y[rowi,1:nbinary_y[rowi]];
    int o0d[ncont_y[rowi] ]= whichcont_y[rowi,1:ncont_y[rowi]];
    
    if(!full){
      o= whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
      o1= whichbinary_y[rowi,1:nbinary_y[rowi]];
      o0= whichcont_y[rowi,1:ncont_y[rowi]];
    }
    if(full){ //needed to calculate yprior and yupd ysmooth
      for(mi in 1:nmanifest) o[mi] = mi;
      o1= whichequals(manifesttype,1,1);
      o0= whichequals(manifesttype,1,0);
    }
    
    if(prevrow != 0 && rowi != 1) T0check = (si==subject[prevrow]) ? (T0check+1) : 0; //if same subject, add one, else zero
    if(T0check > 0){
      dt = time[rowi] - time[prevrow];
      dtchange = continuoustime ? dt!=prevdt : 0; 
      prevdt = dt; //update previous dt store after checking for change
      //prevtime = time[rowi];
    }

    //if(dosmoother && prevrow!=0) eJAx[rowi,,] = eJAx[prevrow,,];
    
    if(T0check == 0) { // calculate initial matrices if this is first row for si
  
  rawindparams=rawpopmeans;
  
  if(si > 0 && nindvarying > 0 && intoverpop==0)  rawindparams[indvaryingindex] += rawpopcovchol * baseindparams[si];

  if(si > 0 &&  ntieffects > 0){
  if(nmissingtipreds > 0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipreds[si]\';
    
    if(nmissingtipreds==0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipredsdata[si]\';
  }

  // compute individual parameters that are not state dependent, either all (if si=0) or just update indvarying ones.
  indparams[whichequals(whenvecp[si ? 2 : 1], 0, 0)]= 
    parvectform(whichequals(whenvecp[si ? 2 : 1], 0, 0),rawindparams\', 
    0, matsetup, matvalues, si)\';
     
  if(whenmat[1, 5] >= (si ? 1 : 0)) T0MEANS = 
    mcalc(T0MEANS, indparams, state, {0}, 1, matsetup, matvalues, si); // base t0means to init
      
 // for(li in 1:nlatentpop) if(!is_nan(T0MEANS[li,1])) state[li] = T0MEANS[li,1]; //in case of t0 dependencies, may have missingness
  
  state=T0MEANS[,1]\';
    
  ',matcalcs('si',when=0:1, c(PARS=10),basemats=TRUE),' //initialise simple PARS then do complex PARS
  ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n\n ',collapse=' '),simplify),'
    
  ',matcalcs('si',when=0:1, mats$t0,basemats=TRUE),'
      
  ',simplifystanfunction(paste0(ctm$modelmats$calcs$t0,';\n\n ',collapse=' '),simplify),'
    for(li in 1:nlatentpop) if(is_nan(state[li])) state[li] = T0MEANS[li,1]; //finish updating state
    
    //init other system matrices (already done PARS, redo t0means in case of PARS dependencies...)
   ',matcalcs('si',when=0,
     matrices = c(mats$base[!mats$base %in% c(1,10)],mats$jacobian),basemats=TRUE),'
    
    if(verbose==2) print("DRIFT = ",DRIFT);
    if(verbose==2) print("indparams = ", indparams);
    
    
 if(si==0 || (sum(whenmat[8,]) + statedep[8]) > 0 ) { // this causes problems but shouldnt -- is t0var being adjusted each iteration when it shouldnt?
   if(intoverpop && nindvarying > 0) T0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopcovbase;
    T0cov = sdcovsqrt2cov(T0VAR,choleskymats); 

    if(intoverpop && nindvarying > 0){ //adjust cov matrix for transforms
    ',if(!gendata) paste0('if(si==0) rawpopcovchol = cholesky_decompose(makesym(T0cov[intoverpopindvaryingindex, intoverpopindvaryingindex],verbose,1));'),'
      for(ri in 1:size(matsetup)){
        if(matsetup[ri,7]==1){ //if t0means
          if(matsetup[ri,5]) { //and indvarying
            T0cov[matsetup[ri,1], ] *= matvalues[ri,2] * matvalues[ri,3]; //multiplier meanscale
            T0cov[, matsetup[ri,1] ] *=  matvalues[ri,2] * matvalues[ri,3]; //multiplier meanscale
          }
        }
      }
    }
 }
  
// if(nt0varstationary > 0) {
//   if(si==0 || (sum(whenmat[8,]) + statedep[8] + sum(whenmat[3,]) + statedep[3] + sum(whenmat[4,]) + statedep[4]) > 0 ){
//     matrix[nlatent,nlatent] stationarycov;
//     stationarycov = ksolve(DRIFT[derrind,derrind], sdcovsqrt2cov(DIFFUSION[derrind,derrind],choleskymats),verbose);
//   for(ri in 1:nt0varstationary){ 
//     T0cov[t0varstationary[ri,1],t0varstationary[ri,2] ] =  stationarycov[t0varstationary[ri,1],t0varstationary[ri,2] ];
//   }
//   }
// }
// 
// if(nt0meansstationary > 0){
//   if(si==0 || //on either pop pars only
//     ( (sum(whenmat[3,])+sum(whenmat[7,])+statedep[3]+statedep[7]) > 0) && savesubjectmatrices ){ // or for each subject
//     
//     if(continuoustime) asymCINT[,1] =  -DRIFT[1:nlatent,1:nlatent] \\ CINT[ ,1 ];
//     if(!continuoustime) asymCINT[,1] =  add_diag(-DRIFT[1:nlatent,1:nlatent],1) \\ CINT[,1 ];
// 
//     if(si==0 || (sum(whenmat[1,]) + statedep[1]) > 0) {
//       for(ri in 1:nt0meansstationary){
//         T0MEANS[t0meansstationary[ri,1]] = 
//           asymCINT[t0meansstationary[ri,1] ];
//       }
//     }
//   }
// }

    etacov=T0cov;
    } //end T0 matrices
    
if(verbose > 1) print ("below t0 row ", rowi);

    if(si==0 || (T0check>0)){ //for init or subsequent time steps when observations exist
      vector[nlatent] base;
      real intstepi = 0;
      
      dtsmall = dt / ceil(dt / maxtimestep);
      
      while(intstepi < (dt-1e-10)){
        intstepi = intstepi + dtsmall;
        ',
    finiteJ(),'
    
    if(statedep[4] || whenmat[4,2]) DIFFUSION=mcalc(DIFFUSION,indparams, state,{2}, 4, matsetup, matvalues, si); 
    if(statedep[52] || whenmat[52,2]) JAx=mcalc(JAx,indparams, state,{2}, 52, matsetup, matvalues, si); 
    ',simplifystanfunction(paste0(ctm$modelmats$calcs$diffusion,';\n\n ',collapse=' '),simplify),'
    
    if(si==0 ||statedep[4] || whenmat[4,2] || ( T0check ==1 && whenmat[4,5])){
      DIFFUSIONcov[derrind,derrind] = sdcovsqrt2cov(DIFFUSION[derrind,derrind],choleskymats);
      if(!continuoustime) discreteDIFFUSION=DIFFUSIONcov;
    }
    
    if(continuoustime && (si==0 || dtchange==1 || statedep[3]|| statedep[52] || whenmat[3,2] || //if first sub or changing every state
      (T0check == 1 && whenmat[3,5]))){ //or first time step of new sub with ind difs
      
      //discreteDRIFT = expm2(append_row(append_col(DRIFT[1:nlatent, 1:nlatent],CINT),nlplusonezerovec\') * dtsmall);
      discreteDRIFT = expmSubsets(DRIFT * dtsmall,DRIFTsubsets);
      if(!JAxDRIFTequiv){ 
        eJAx =  expmSubsets(JAx * dtsmall,JAxsubsets);
      } else eJAx[1:nlatent, 1:nlatent] = discreteDRIFT;
                             
      if(si==0 || statedep[3] || statedep[4]||statedep[52]||  //if first pass or state dependent
        whenmat[4,2] || whenmat[3,2] ||
        (T0check == 1 && (whenmat[3,5]  || whenmat[4,5]))){ //or first time step of new sub with ind difs
        asymDIFFUSIONcov[derrind,derrind] = ksolve(JAx[derrind,derrind], DIFFUSIONcov[derrind,derrind],verbose);
      }
      
      discreteDIFFUSION[derrind,derrind] =  asymDIFFUSIONcov[derrind,derrind] - 
        quad_form_sym( asymDIFFUSIONcov[derrind,derrind], eJAx[derrind,derrind]\' );
        
      for(li in 1:nlatent) if(is_nan(state[li]) || is_nan(sum(discreteDRIFT[li,]))) {
        print("Possible time step problem? Intervals too large? Try reduce maxtimestep");
      }
        
    } //end discrete drift / diffusion coef calcs based on ct
          
          
    if(continuoustime) state[1:nlatent] *= discreteDRIFT\'; 
    if(!continuoustime) state[1:nlatent] *= DRIFT\'; 
    
    if(intoverstates==1 || dosmoother==1){
      if(continuoustime){
        etacov = quad_form_sym(makesym(etacov,verbose,1), eJAx\');
        etacov[derrind,derrind] += discreteDIFFUSION[derrind,derrind]; 
      }
      if(!continuoustime){
        etacov = quad_form_sym(makesym(etacov,verbose,1), JAx\');
        etacov[ derrind, derrind ] += DIFFUSIONcov[ derrind, derrind ]; 
      }
    }
      
    if(continuoustime && dosmoother && intstepi >= (dt-1e-10)) eJAxs[rowi,,] = expmSubsets(JAx * dt,JAxsubsets); //save approximate exponentiated jacobian for smoothing
    if(!continuoustime && dosmoother) eJAxs[rowi,,] = JAx;
    
    if(size(CINTnonzero)){
      if(continuoustime){
        if(si==0 || dtchange==1 || statedep[3]|| statedep[7] || whenmat[3,2] || whenmat[7,2] || // state depenency
          (T0check == 1 && (whenmat[7,5] || whenmat[3,5]))){ //or ind difs
          discreteCINT = (DRIFT \\ (discreteDRIFT-IIlatentpop[1:nlatent,1:nlatent])) * CINT[,1];
        }
        state[1:nlatent] += discreteCINT\';
      }
    if(!continuoustime) state[CINTnonzero]+= CINT[CINTnonzero,1]\';
    } // end cint section

    } // end time step loop
  } // end non linear time update
    
    if(ntdpred > 0) {
      int nonzerotdpred = 0;
      for(tdi in 1:ntdpred) if(tdpreds[rowi,tdi] != 0.0) nonzerotdpred = 1;
      if(si==0 ||nonzerotdpred){
          
        //initialise PARS first, and simple PARS before complex PARS
        if(statedep[10] || whenmat[10,3]) PARS=mcalc(PARS,indparams, state,{3}, 10, matsetup, matvalues, si); 
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n\n ',collapse=' '),simplify),'
        
      
        if(statedep[9] || whenmat[9,3]) TDPREDEFFECT=mcalc(TDPREDEFFECT,indparams, state,{3}, 9, matsetup, matvalues, si); 
        if(statedep[53] || whenmat[53,3]) Jtd=mcalc(Jtd,indparams, state,{3}, 53, matsetup, matvalues, si); 

        ',simplifystanfunction(paste0(ctm$modelmats$calcs$tdpred,';\n\n ',collapse=' '),simplify),'

        state[1:nlatent] +=   (TDPREDEFFECT * tdpreds[rowi])\'; //tdpred effect only influences at observed time point','
        etacov = quad_form_sym(makesym(etacov,verbose,1),Jtd\');  //could speed up by detecting if non diagonal Jtd
      }
    }//end nonlinear tdpred

  if(si > 0 && intoverstates==0){ //unused states if intoverpop is specified, consider fixing...
    if(T0check==0) state += (cholesky_decompose(etacov) * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)])\';
    if(T0check>0) state[derrind] +=  (cholesky_decompose(makesym(discreteDIFFUSION[derrind,derrind],verbose,1)) * 
     (etaupdbasestates[(1+(rowi-1)*nlatentpop):(nlatent+(rowi-1)*nlatentpop)])[derrind])\';
     
   // if(T0check==0) llrow[rowi]+= multi_normal_cholesky_lpdf(
  //     etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)] | rep_vector(0,nlatent), etacov);
  //  if(T0check>0) llrow[rowi]+= multi_normal_lpdf(
  //     etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)] | rep_vector(0,nlatent), discreteDIFFUSION);
  //  state+=etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
  }

if(verbose > 1){
  print("etaprior = ", state);
  print("etapriorcov = ", etacov);
}

',if(!ppchecking) 'if(dosmoother){
  etacovb[1,rowi] = etacov; 
  etab[1,rowi] = state\';
}','

//measurement update

  if(si == 0 || nobs_y[rowi] > 0 || dosmoother){ //measurement init
    int zeroint[1];
    row_vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(Jyfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0 && (dosmoother + intoverstates) > 0)  state[statei] += Jstep;

        //initialise PARS first, and simple PARS before complex PARS
        if(statedep[10] || whenmat[10,4]) PARS=mcalc(PARS,indparams, state,{4}, 10, matsetup, matvalues, si); 
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n\n ',collapse=' '),simplify),'
        
      
        if(statedep[2] || whenmat[2,4]) LAMBDA=mcalc(LAMBDA,indparams, state,{4}, 2, matsetup, matvalues, si); 
        if(statedep[5] || whenmat[5,4]) MANIFESTVAR=mcalc(MANIFESTVAR,indparams, state,{4}, 5, matsetup, matvalues, si); 
        if(statedep[6] || whenmat[6,4]) MANIFESTMEANS=mcalc(MANIFESTMEANS,indparams, state,{4}, 6, matsetup, matvalues, si); 
        if(statedep[54] || whenmat[54,4]) Jy=mcalc(Jy,indparams, state,{4}, 54, matsetup, matvalues, si); 

        ',simplifystanfunction(paste0(ctm$modelmats$calcs$measurement,';\n\n ',collapse=' '),simplify),'
        
      if(statei > 0 && (intoverstates) > 0) {
        Jy[o,statei] =  LAMBDA[o] * state[1:nlatent]\' + MANIFESTMEANS[o,1]; //compute new change
        Jy[o1,statei] = to_vector(inv_logit(to_array_1d(Jy[o1,statei])));
        if(verbose>1) print("Jy ",Jy);
      }
      if(statei==0){
        syprior[o] = LAMBDA[o] * state[1:nlatent]\' + MANIFESTMEANS[o,1];
        syprior[o1] = to_vector(inv_logit(to_array_1d( syprior[o1] )));
        if(size(Jyfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
          if(verbose>1) print("syprior = ",syprior,"    Jyinit= ",Jy);
          for(fi in Jyfinite){
            Jy[o,fi] -= syprior[o];
            Jy[o,fi] /= Jstep; //new - baseline change divided by stepsize
          }
        }
      }
    }
    if(verbose>1) print("Jy ",Jy);
  } //end measurement init
      
  if(si==0 || whenmat[5,5] || whenmat[5,4] || statedep[5]) MANIFESTcov = sdcovsqrt2cov(MANIFESTVAR,choleskymats);

 
  if(si > 0 && (nobs_y[rowi] > 0 || dosmoother)){   //if not just inits...

      if(intoverstates==1 || dosmoother==1) { //classic kalman
        ycov[o,o] = quad_form_sym(makesym(etacov,verbose,1), Jy[o,]\') + MANIFESTcov[o,o]; // previously shifted measurement error down, but reverted
        for(wi in 1:nmanifest){ 
          // if(Y[rowi,wi] != 99999 || dosmoother==1) ycov[wi,wi] += square(MANIFESTVAR[wi,wi]);
          if(manifesttype[wi]==1 && (Y[rowi,wi] != 99999  || dosmoother==1)) ycov[wi,wi] += fabs((syprior[wi] - 1) .* (syprior[wi]));
          if(manifesttype[wi]==2 && (Y[rowi,wi] != 99999  || dosmoother==1)) ycov[wi,wi] += square(fabs((syprior[wi] - round(syprior[wi])))); 
        }
      }
        
      if(intoverstates==0 && ncont_y[rowi] > 0) ypriorcov_sqrt[o,o] = 
        cholesky_decompose(makesym(MANIFESTcov[o,o],verbose,1));

        
     
'    
  ,if(ppchecking) paste0('
  {
  int skipupd = 0;
        for(vi in 1:nobs_y[rowi]){
            if(fabs(syprior[od[vi]]) > 1e10 || is_nan(syprior[od[vi]]) || is_inf(syprior[od[vi]])) {
              skipupd = 1; 
              syprior[od[vi]] =99999;
  if(verbose > 1) print("pp syprior problem! row ", rowi);
            }
          }
        if(skipupd==0){ 
          if(ncont_y[rowi] > 0){
            ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(makesym(ycov[o0d, o0d],verbose,1)); 
            Ygen[ rowi, o0d] = syprior[o0d] + ypriorcov_sqrt[o0d,o0d] * Ygenbase[rowi,o0d];
          }
          if(nbinary_y[rowi] > 0) for(obsi in 1:size(o1d)) Ygen[rowi, o1d[obsi]] = (syprior[o1d[obsi]] > Ygenbase[rowi,o1d[obsi]]) ? 1 : 0; 
          for(vi in 1:nobs_y[rowi]) if(is_nan(Ygen[rowi,od[vi]])) {
            Ygen[rowi,od[vi]] = 99999;
            print("pp ygen problem! row ", rowi);
          }
        err[od] = Ygen[rowi,od] - syprior[od]; // prediction error
        }
}
      '), 
  
  if(!ppchecking) 'err[od] = Y[rowi,od] - syprior[od]; // prediction error','
    
      if(intoverstates==1 && size(od) > 0) {
        
         if(verbose > 1) print("before K rowi =",rowi, "  si =", si, "  state =",state, "  etacov ",etacov,
          " indparams = ", indparams,
            "  syprior[o] =",syprior[o],"  ycov[o,o] ",ycov[o,o], 
            "  PARS = ", PARS, 
            "  DRIFT =", DRIFT, " DIFFUSION =", DIFFUSION, 
            " CINT =", CINT, "  discreteCINT = ", discreteCINT, "  MANIFESTcov ", (MANIFESTcov), "  MANIFESTMEANS ", MANIFESTMEANS, 
            "  T0cov", T0cov,  " T0MEANS ", T0MEANS, "LAMBDA = ", LAMBDA, "  Jy = ",Jy,
            " discreteDRIFT = ", discreteDRIFT, "  discreteDIFFUSION ", discreteDIFFUSION, "  asymDIFFUSIONcov ", asymDIFFUSIONcov, 
            " DIFFUSIONcov = ", DIFFUSIONcov,
            " eJAx = ", eJAx,
            "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
         
        K[,od] = mdivide_right_spd(etacov * Jy[od,]\', makesym(ycov[od,od],verbose,1)); // * multiply_lower_tri_self_transpose(ycovi\');// ycov[od,od]; 
        etacov += -K[,od] * Jy[od,] * etacov; //cov update
        state +=  (K[,od] * err[od])\'; //state update
      }
      
      ',if(savemodel) 'if(dosmoother==1) {
        yb[1,rowi] = syprior[o];
        etab[2,rowi] = state\';
        ycovb[1,rowi] = ycov;
        etacovb[2,rowi] = etacov;
        ycovb[2,rowi] = quad_form_sym(makesym(etacov,verbose,1), Jy\') + MANIFESTcov;
        yb[2,rowi] = MANIFESTMEANS[o,1] + LAMBDA[o,] * state[1:nlatent]\';
        Jys[rowi,,] = Jy;
      }','
      
      
      if(verbose > 1) print(" After K rowi =",rowi, "  si =", si, "  state =",state,"  etacov ",etacov,"  K[,o] ",K[,o]);
        
  //likelihood stuff
      if(nbinary_y[rowi] > 0) ',ifelse(savemodel,'llrow[rowi]','ll'),' += sum(log(1e-10+',ifelse(gendata,'Ygen','Y'),'[rowi,o1d] .* (syprior[o1d]) + (1-',ifelse(gendata,'Ygen','Y'),'[rowi,o1d]) .* (1-syprior[o1d]))); 

      if(size(o0d) > 0 && (llsinglerow==0 || llsinglerow == rowi)){
        if(intoverstates==1) ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(ycov[o0d,o0d]); //removed makesym
         ',ifelse(savemodel,'llrow[rowi]','ll'),' +=  multi_normal_cholesky_lpdf(',ifelse(gendata,'Ygen','Y'),'[rowi,o0d] | syprior[o0d], ypriorcov_sqrt[o0d,o0d]);
         //errtrans[counter:(counter + ncont_y[rowi]-1)] = 
           //mdivide_left_tri_low(ypriorcov_sqrt[o0d,o0d], err[o0d]); //transform pred errors to standard normal dist and collect
         //ll+= -sum(log(diagonal(ypriorcov_sqrt[o0d,o0d]))); //account for transformation of scale in loglik
         //counter += ncont_y[rowi];
      }
      if(verbose > 1) print(llrow[rowi]);
      
    }//end si > 0 nobs > 0 section
    
  ',if(savemodel) '     // store system matrices
       
    if(si==0 || //on either pop pars only
      (  (sum(whenmat[3,])+sum(whenmat[7,])+statedep[3]+statedep[7]) > 0 && savesubjectmatrices) ){ // or for each subject
      if(continuoustime==1) asymCINT[,1] =  -DRIFT[1:nlatent,1:nlatent] \\ CINT[ ,1 ];
      if(continuoustime==0) asymCINT[,1] =  add_diag(-DRIFT[1:nlatent,1:nlatent],1) \\ CINT[,1 ];
    }

  
  if(!continuoustime){
    if(si==0 || //on either pop pars only
    (  (sum(whenmat[3,])+sum(whenmat[4,])+statedep[3]+statedep[4]) > 0 && savesubjectmatrices) ){ // or for each subject
  
      asymDIFFUSIONcov[ derrind, derrind ] = 
        to_matrix( (add_diag( -sqkron_prod(JAx[ derrind, derrind ], JAx[ derrind, derrind ]),1)) \\  
          to_vector(DIFFUSIONcov[ derrind, derrind ]), ndiffusion, ndiffusion);
    }
  }', '
      
    
  if(si == 0){
',if(savemodel) paste0(collectsubmats(popmats=TRUE),collapse=' '),'
  }

  
  ',if(!ppchecking && savemodel) paste0('if(si > 0 && dosmoother && (rowi==ndatapoints || subject[rowi+1] != subject[rowi])){ //at subjects last datapoint, smooth
    int sri = rowi;
    while(sri>0 && subject[sri]==si){
      if(sri==rowi) {
        etab[3,sri]=etab[2,sri];
        yb[3,sri]=yb[2,sri];
        etacovb[3,sri]=etacovb[2,sri];
        ycovb[3,sri]=ycovb[2,sri];
      } else{
        matrix[nlatentpop,nlatentpop] smoother;
        smoother = etacovb[2,sri] * eJAxs[sri+1,,]\' / makesym(etacovb[1,sri+1],verbose,1);
        etab[3,sri,]= etab[2,sri,] + smoother * (etab[3,sri+1,] - etab[1,sri+1,]);
        etacovb[3,sri]= etacovb[2,sri] + smoother * ( etacovb[3,sri+1] - etacovb[1,sri+1]) * smoother\';
        yb[3,sri,] = yb[2,sri,] + Jys[sri,,] * (etab[3,sri,] - etab[2,sri,]);
        ycovb[3,sri] = ycovb[2,sri] + Jys[sri,,] * (etacovb[3,sri] - etacovb[2,sri]) * Jys[sri,,]\';
      }
      state=etab[3,sri,]\'; //update for t0means saving
      sri += -1;
      while(sri > 0 && dokalmanrows[sri]==0) sri+= -1; //skip rows if requested

      if(savesubjectmatrices && //if getting subj matrices and 
        (sri ? subject[sri] != subject[sri+1] : 1)){ //no more rows, or change of subject
     ',paste0(collectsubmats(popmats=FALSE),collapse=' '),'
     if(sum(whenmat[1,1:5]) > 0 || statedep[1]) subj_T0MEANS[si,,1] = state\'; //t0means updated, other pars as per final time point
        }
      
    }
  } //end smoother

  '),'
 } // end si loop (includes sub 0)
  
  prevrow = rowi; //update previous row marker only after doing necessary calcs
}//end active rowi

',if(savemodel) 'if(savescores){
  ya=yb;
  ycova=ycovb;
  etaa=etab;
  etacova=etacovb;
}
ll+=sum(llrow);','

')
    
    return(out)
  }

matcalcs <- function(subjectid,when, matrices, basemats){
  paste0(sapply(matrices, function(x){
    mn=names(matrices[matrices == x])
    whena=paste0('{',paste0(when,collapse=','),'}')
    whenax=when
    whenax[whenax==0] <- 5
    whenax=paste0('{',paste0(whenax,collapse=','),'}') #remove 0 
    out = paste0(
      'if(',
      ifelse(basemats,'si==0 || ',''),
      paste0('sum(whenmat[',x,',',whenax,']) > 0 )'),
      paste0(mn,'=mcalc(',mn,',indparams, state,',whena,', ',x,', matsetup, matvalues, ',subjectid,'); \n') )
  }),collapse='')
}


subjectparaminit<- function(popmats=FALSE,smats=TRUE,matrices=c(mats$base,31, 32, 33, 21,22)){
  if(smats && popmats) stop('smats and popmats cannot both be TRUE!')
  ma <- ctStanMatricesList()$all
  out<-''
  for(mn in matrices){ #removed if(smats) 's',
    m=names(ma)[ma %in% mn]
    out <- paste0(out, '
      matrix[matrixdims[',mn,', 1], matrixdims[',mn,', 2] ] ',
      if(popmats) 'pop_',
      if(!smats && !popmats) paste0('subj_'),
      m,
      if(!smats && !popmats) paste0('[ (savesubjectmatrices && (sum(whenmat[',mn,',1:5]) || statedep[',mn,'])) ? nsubjects : 0]'),
      ';')
  }
  
  return(out)
}

collectsubmats <- function(popmats=FALSE,matrices=c(mats$base,31, 32,33,21,22)){ #'DIFFUSIONcov','MANIFESTcov','asymDIFFUSIONcov','asymCINT'
  ma <- ctStanMatricesList()$all
  out<-''
  for(mn in matrices){
    m=names(ma)[ma %in% mn]
    if(!popmats) out <- paste0(out, '
    if(sum(whenmat[',mn,',1:5]) > 0 || statedep[',mn,']) subj_',m,'[si] = ',m,';')
    
    if(popmats) out <- paste0(out, 'pop_',m,' = ',m,'; ')
  }
  return(out)
}



writemodel<-function(){
  paste0('
functions{

 int[] vecequals(int[] a, int test, int comparison){ //do indices of a match test condition?
    int check[size(a)];
    for(i in 1:size(check)) check[i] = comparison ? (test==a[i]) : (test!=a[i]);
    return(check);
  }

int[] whichequals(int[] b, int test, int comparison){  //return array of indices of b matching test condition
    int check[size(b)] = vecequals(b,test,comparison);
    int which[sum(check)];
    int counter = 1;
    if(size(b) > 0){
      for(i in 1:size(b)){
        if(check[i] == 1){
          which[counter] = i;
          counter += 1;
        }
      }
    }
    return(which);
  }

 
  matrix constraincorsqrt(matrix mat){ //converts from unconstrained lower tri matrix to cor
    int d=rows(mat);
    matrix[d,d] o;
    vector[d] ss = rep_vector(0,d);
    vector[d] s = rep_vector(0,d);
    real r;
    real r3;
    real r4;
    real r1;
    
    for(i in 1:d){
      for(j in 1:d){
        if(j > i) {
          ss[i] +=square(mat[j,i]);
          s[i] +=mat[j,i];
        }
        if(j < i){
          ss[i] += square(mat[i,j]);
          s[i] += mat[i,j];
        }
      }
      s[i] += 1e-5;
      ss[i] += 1e-5;
    }

    
    for(i in 1:d){
      o[i,i]=0;
      r1=sqrt(ss[i]);
      r3=(fabs(s[i]))/(r1)-1;
      r4=sqrt(log1p_exp(2*(fabs(s[i])-s[i]-1)-4));
      r=(r4*((r3))+1)*r4+1;
      r=(sqrt(ss[i]+r));
      for(j in 1:d){
        if(j > i)  o[i,j]=mat[j,i]/r;
        if(j < i) o[i,j] = mat[i,j] /r;
      }
      o[i,i]=sqrt(1-sum(square(o[i,]))+1e-5);
    }

    return o;
}  

  matrix sdcovsqrt2cov(matrix mat, int choleskymats){ //covariance from cholesky or unconstrained cor sq root
    if(choleskymats< 1) {
      //if(choleskymats== -1){
        return(tcrossprod(diag_pre_multiply(diagonal(mat),constraincorsqrt(mat))));
      //} else {
      //  return(quad_form_diag(constraincorsqrt(mat),diagonal(mat)));
      //}
      } else return(tcrossprod(mat));
  }

  matrix sqkron_prod(matrix mata, matrix matb){
    int d=rows(mata);
    matrix[d*d,d*d] out;
    for (k in 1:d){
      for (l in 1:d){
        for (i in 1:d){
          for (j in 1:d){
            out[ d*(i-1)+k, d*(j-1)+l ] = mata[i, j] * matb[k, l];
          }
        }
      }
    }
    return out;
  }

 
 matrix ksolve(matrix A, matrix Q, int verbose){
  int d= rows(A);
  int d2= (d*d-d)/2;
  matrix[d+d2,d+d2] O;
  vector[d+d2] triQ;
  matrix[d,d] AQ;
  int z=0; //z is row of output
  for(j in 1:d){//for column reference of solution vector
    for(i in 1:j){ //and row reference...
      if(j >= i){ //if i and j denote a covariance parameter (from upper tri)
        int y=0; //start new output row
        z+=1; //shift current output row down
        
        for(ci in 1:d){//for columns and
          for(ri in 1:d){ //rows of solution
            if(ci >= ri){ //when in upper tri (inc diag)
              y+=1; //move to next column of output
              
              if(i==j){ //if output row is for a diagonal element
                if(ri==i) O[z,y] = 2*A[ri,ci];
                if(ci==i) O[z,y] = 2*A[ci,ri];
              }
              
              if(i!=j){ //if output row is not for a diagonal element
                if(y==z) O[z,y] = A[ri,ri] + A[ci,ci]; //if column of output matches row of output, sum both A diags
                if(y!=z){ //otherwise...
                  // if solution element we refer to is related to output row...
                  if(ci==ri){ //if solution element is a variance
                    if(ci==i) O[z,y] = A[j,ci]; //if variance of solution corresponds to row of our output
                    if(ci==j) O[z,y] = A[i,ci]; //if variance of solution corresponds to col of our output
                  }
                  if(ci!=ri && (ri==i||ri==j||ci==i||ci==j)){//if solution element is a related covariance
                    //for row 1,2 / 2,1 of output, if solution row ri 1 (match) and column ci 3, we need A[2,3]
                    if(ri==i) O[z,y] = A[j,ci];
                    if(ri==j) O[z,y] = A[i,ci];
                    if(ci==i) O[z,y] = A[j,ri];
                    if(ci==j) O[z,y] = A[i,ri];
                  }
                }
              }
              if(is_nan(O[z,y])) O[z,y]=0;
            }
          }
        }
      }
    }
  }
  
  z=0; //get upper tri of Q
  for(j in 1:d){
    for(i in 1:j){
    z+=1;
    triQ[z] = Q[i,j];
    }
  }
  triQ=-O \\ triQ; //get upper tri of asymQ
  
    z=0; // put upper tri of asymQ into matrix
  for(j in 1:d){
    for(i in 1:j){
    z+=1;
    AQ[i,j] = triQ[z];
    if(i!=j) AQ[j,i] = triQ[z];
    }
  }
  
  if(verbose>1) print("AQ = ", AQ, "   triQ = ", triQ, "   O = ", O);
  
  return AQ;
}

  matrix makesym(matrix mat, int verbose, int pd){
    matrix[rows(mat),cols(mat)] out;
    for(coli in 1:cols(mat)){
    //  if(pd ==1 && mat[coli,coli] < 1e-5){
     //   out[coli,coli] = 1e-5;// 
     // } else 
      out[coli,coli] = mat[coli,coli] + 1e-10; 
      for(rowi in 1:rows(mat)){
        if(rowi > coli) {
          out[rowi,coli] = mat[rowi,coli];
          out[coli,rowi] = mat[rowi,coli];
        }
        ',if(gendata) paste0('if(is_nan(out[rowi,coli])){
          if(verbose > 0) print("nan during makesym row ", rowi, " col ", coli);
          if(rowi==coli) out[rowi,coli] = 99999;
          if(rowi!=coli) {
            out[rowi,coli] = 0;
            out[coli,rowi] = 0;
          }
        }'),'
      }
    }
    return out;
  }

  real tform(real parin, int transform, data real multiplier, data real meanscale, data real offset, data real inneroffset){
    real param=parin;
    ',tformshapes(stan=TRUE),
if(length(extratforms) > 0) paste0(extratforms,collapse=" \n"),'
    return param;
  }
  
  // improve PARS when = 100 thing here too
   row_vector parvectform(int[] which, row_vector rawpar, int when, int[,] ms, data real[,] mval, int subi){
    row_vector[size(which)] parout;
    if(size(which)){
      for(whichout in 1:size(which)){
        int done=0; //only want to tform once, may be copies
      for(ri in 1:size(ms)){ //for each row of matrix setup
        if(!done){
        if((ms[ri,8]==when || ms[ri,8]==100)  && ms[ri,3] == which[whichout]){ //if correct when and free parameter //,not a copyrow,&& ms[ri,9] < 1
          if(subi ==0 ||  //if population parameter
            (ms[ri,3] > 0 && (ms[ri,5] > 0 || ms[ri,6] > 0 || ms[ri,8] > 0)) //or there is individual variation
            ){ //otherwise repeated values
            
            parout[whichout] = tform(rawpar[ms[ri,3] ], // was: whichequals(outwhen, ms[ri,3],1)[1], which outwhen refers to correct par
              ms[ri,4], mval[ri,2], mval[ri,3], mval[ri,4], mval[ri,6] ); 
            }
           done=1;
        }
      }
      }
      }
    }
  return parout;
  }
  
  
  matrix mcalc(matrix matin, vector tfpars, row_vector states, int[] when, int m, int[,] ms, data real[,] mval, int subi){
    matrix[rows(matin),cols(matin)] matout;

    for(ri in 1:size(ms)){ //for each row of matrix setup
      int whenyes = 0;
      for(wi in 1:size(when)) if(when[wi]==ms[ri,8] || ms[ri,8]==100) whenyes = 1; //improve PARS when = 100 thing
      if(m==ms[ri,7] && whenyes){ // if correct matrix and when

        if(subi ==0 ||  //if population parameter
          (ms[ri,3] > 0 && (ms[ri,5] > 0 || ms[ri,6] > 0 || ms[ri,8] > 0)) //or there is individual variation
          ){ //otherwise repeated values (maybe this check not needed now?

          if(ms[ri,3] > 0 && ms[ri,8]==0)  matout[ms[ri,1], ms[ri,2] ] = tfpars[ms[ri,3]]; //should be already tformed
          if(ms[ri,3] > 0 && ms[ri,8]>0)  matout[ms[ri,1], ms[ri,2] ] =   //if references param and is state based
            tform(states[ms[ri,3] ], ms[ri,4], mval[ri,2], mval[ri,3], mval[ri,4], mval[ri,6] );
          if(ms[ri,3] < 1) matout[ms[ri,1], ms[ri,2] ] = mval[ri, 1]; //doing this once over all subjects unless covariance matrix -- speed ups possible here, check properly!
        }
      }
    }
    for(ri in 1:rows(matin)){ //fill holes with unchanged input matrix
      for(ci in 1:cols(matin)){
        if(is_nan(matout[ri,ci]) && !is_nan(matin[ri,ci])) matout[ri,ci] = matin[ri,ci];
      }
    }
  return(matout);
  }
  
  matrix expmSubsets(matrix m, int[,] subsets){
    int nr = rows(m);
    matrix[nr,nr] e;
    for(si in 1:size(subsets)){
      int n=0;
      for(j in 1:nr){
        if(subsets[si,j]!=0) n+=1;
      }
      e[subsets[si][1:n],subsets[si][1:n]] = matrix_exp(m[subsets[si][1:n],subsets[si][1:n]]);
    }
    return e;
  }
  
}
data {
  int<lower=0> ndatapoints;
  int<lower=1> nmanifest;
  int<lower=1> nlatent;
  int nlatentpop;
  int nsubjects;
  int<lower=0> ntipred; // number of time independent covariates
  int<lower=0> ntdpred; // number of time dependent covariates
  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipredsdata;
  int nmissingtipreds;
  int ntipredeffects;
  real<lower=0> tipredsimputedscale;
  real<lower=0> tipredeffectscale;

  vector[nmanifest] Y[ndatapoints];
  int priors;
  vector[ntdpred] tdpreds[ndatapoints];
  
  real maxtimestep;
  real time[ndatapoints];
  int subject[ndatapoints];
  int<lower=0> nparams;
  int continuoustime; // logical indicating whether to incorporate timing information
  int nindvarying; // number of subject level parameters that are varying across subjects
  int nindvaryingoffdiagonals; //number of off diagonal parameters needed for popcov matrix
  vector[nindvarying] sdscale;
  int indvaryingindex[nindvarying];
  int notindvaryingindex[nparams-nindvarying];
  
//  int nt0varstationary;
//  int nt0meansstationary;
//  int t0varstationary [nt0varstationary, 2];
//  int t0meansstationary [nt0meansstationary, 2];

  int nobs_y[ndatapoints];  // number of observed variables per observation
  int whichobs_y[ndatapoints, nmanifest]; // index of which variables are observed per observation
  int ndiffusion; //number of latents involved in system noise calcs
  int derrind[ndiffusion]; //index of which latent variables are involved in system noise calculations

  int manifesttype[nmanifest];
  int nbinary_y[ndatapoints];  // number of observed binary variables per observation
  int whichbinary_y[ndatapoints, nmanifest]; // index of which variables are observed and binary per observation
  int ncont_y[ndatapoints];  // number of observed continuous variables per observation
  int whichcont_y[ndatapoints, nmanifest]; // index of which variables are observed and continuous per observation
  
  int intoverpop;
  int statedep[',max(mats$all),'];
  int choleskymats;
  int intoverstates;
  int verbose; //level of printing during model fit
  int TIPREDEFFECTsetup[nparams, ntipred];
  int nrowmatsetup;
  int matsetup[nrowmatsetup,9];
  real matvalues[nrowmatsetup,6];
  int whenmat[',max(mats$all),',5];
  int whenvecp[2,nparams];
  int whenvecs[6,nlatentpop];
  int matrixdims[',max(mats$all),',2];
  int savescores;
  int savesubjectmatrices;
  int dokalman;
  int dokalmanrows[ndatapoints];
  int nsubsets;
  real Jstep;
  real priormod;
  int intoverpopindvaryingindex[intoverpop ? nindvarying : 0];
  int nJAxfinite;
  int JAxfinite[nJAxfinite];
  int nJyfinite;
  int Jyfinite[nJyfinite];
  int taylorheun;
  int popcovn;
  int llsinglerow;
  int laplaceprior[nparams];
  int laplaceprioronly;
  int laplacetipreds;
  int CINTnonzerosize;
  int CINTnonzero[CINTnonzerosize];
  int JAxDRIFTequiv;
  
  int nDRIFTsubsets;
  int nJAxsubsets;
  int DRIFTsubsets[nDRIFTsubsets,nlatent];
  int JAxsubsets[nJAxsubsets,nlatentpop];
}
      
transformed data{
  matrix[nlatent+nindvarying,nlatent+nindvarying] IIlatentpop = diag_matrix(rep_vector(1,nlatent+nindvarying));
  vector[nlatentpop-nlatent] nlpzerovec = rep_vector(0,nlatentpop-nlatent);
  vector[nlatent+1] nlplusonezerovec = rep_vector(0,nlatent+1);
  int tieffectindices[nparams]=rep_array(0,nparams);
  int ntieffects = 0;
  int dosmoother = savescores || savesubjectmatrices;
  
  if(ntipred >0){
    for(pi in 1:nparams){
      if(sum(TIPREDEFFECTsetup[pi,]) > .5){
      ntieffects+=1;
      tieffectindices[ntieffects] = pi;
      }
    }
  }
  
}
      
parameters{
  vector[nparams] rawpopmeans; // population level means \n','
  vector',if(!is.na(ctm$rawpopsdbaselowerbound)) paste0('<lower=',ctm$rawpopsdbaselowerbound[1],'>'),'[nindvarying] rawpopsdbase; //population level std dev
  vector[nindvaryingoffdiagonals] sqrtpcov; // unconstrained basis of correlation parameters
  vector[intoverpop ? 0 : nindvarying] baseindparams[intoverpop ? 0 : nsubjects]; //vector of subject level deviations, on the raw scale
  
  vector[ntipredeffects] tipredeffectparams; // effects of time independent covariates
  vector[nmissingtipreds] tipredsimputed;
  
  ',if(!gendata) 'vector[intoverstates ? 0 : nlatentpop*ndatapoints] etaupdbasestates; //sampled latent states posterior','
  vector[(nsubsets > 1) ? 1 : 0] subsetpar;
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying, nindvarying] rawpopcovbase;
  matrix[nindvarying, nindvarying] rawpopcov;
  matrix[nindvarying, nindvarying] rawpopcovchol;
  matrix[nindvarying, nindvarying] rawpopcorr;
  real subset = (nsubsets > 1) ? subsetpar[1] : 1.0;
  real firstsub = round(nsubjects*1.0/nsubsets*(subset-1)+1);
  real lastsub = round(nsubjects*1.0/nsubsets*(subset));
  ',if(!gendata) 'real ll = 0;
',if(!gendata & savemodel) paste0('
  vector[dokalman ? ndatapoints : 1] llrow = rep_vector(0,dokalman ? ndatapoints : 1);
  matrix[nlatentpop,nlatentpop] etacova[3,savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ycova[3,savescores ? ndatapoints : 0];
  vector[nlatentpop] etaa[3,savescores ? ndatapoints : 0];
  vector[nmanifest] ya[3,savescores ? ndatapoints : 0];
  vector[',ifelse(is.null(ctm$modelmats$calcs$calcNames),0,length(ctm$modelmats$calcs$calcNames)),'] calcs[',ifelse(is.null(ctm$modelmats$calcs$calcNames),0,'ndatapoints'),'];
  ',subjectparaminit(pop=TRUE,smats=FALSE),
  subjectparaminit(pop=FALSE,smats=FALSE),
  collapse=''),'

  matrix[ntipred ? (nmissingtipreds ? nsubjects : 0) : 0, ntipred ? (nmissingtipreds ? ntipred : 0) : 0] tipreds; //tipred values to fill from data and, when needed, imputation vector
  matrix[nparams, ntipred] TIPREDEFFECT; //design matrix of individual time independent predictor effects

  if(ntipred > 0){ 
    if(nmissingtipreds > 0){
    int counter = 0;
    for(coli in 1:cols(tipreds)){ //insert missing ti predictors
      for(rowi in 1:rows(tipreds)){
        if(tipredsdata[rowi,coli]==99999) {
          counter += 1;
          tipreds[rowi,coli] = tipredsimputed[counter];
        } else tipreds[rowi,coli] = tipredsdata[rowi,coli];
      }
    }
    }
    for(ci in 1:ntipred){ //configure design matrix
      for(ri in 1:nparams){
        if(TIPREDEFFECTsetup[ri,ci] > 0) {
          TIPREDEFFECT[ri,ci] = tipredeffectparams[TIPREDEFFECTsetup[ri,ci]];
        } else {
          TIPREDEFFECT[ri,ci] = 0;
        }
      }
    }
  }

  if(nindvarying > 0){
    int counter =0;
    rawpopsd = ',ctm$rawpopsdtransform, ' + 1e-10; // sqrts of proportions of total variance
    for(j in 1:nindvarying){
      rawpopcovbase[j,j] = rawpopsd[j]; //used with intoverpop
      for(i in 1:nindvarying){
        if(i > j){
          counter += 1;
          rawpopcovbase[i,j]=inv_logit(sqrtpcov[counter])*2-1;
          rawpopcovbase[j,i]=0;// needed to avoid nan output;
        }
      }
    }
    //if(choleskymats==0) rawpopcorr = constraincorsqrt(rawpopcovbase);
    //if(choleskymats== -1) 
    rawpopcorr = tcrossprod( constraincorsqrt(rawpopcovbase));
    rawpopcov = makesym(quad_form_diag(rawpopcorr, rawpopsd +1e-8),verbose,1);
    rawpopcovchol = cholesky_decompose(rawpopcov); 
  }//end indvarying par setup

  {
',
if(!gendata) ukfilterfunc(ppchecking=FALSE),'
  }
}
      
model{
  real priormod2 = priormod / nsubsets;
  if(intoverpop==0 && nindvarying > 0) target+= multi_normal_cholesky_lpdf(baseindparams | rep_vector(0,nindvarying), IIlatentpop[1:nindvarying,1:nindvarying]);

  if(ntipred > 0){ 
    if(priors && laplacetipreds==0) target+= priormod2 * normal_lpdf(tipredeffectparams / tipredeffectscale| 0, 1);
    if(priors && laplacetipreds==1) for(i in 1:ntipredeffects) target+= priormod2 * double_exponential_lpdf(pow(fabs(tipredeffectparams[i]),1+.1/((tipredeffectparams[i]*100)^2+.1)) / tipredeffectscale| 0, 1);
    target+= normal_lpdf(tipredsimputed| 0, tipredsimputedscale); //consider better handling of this when using subset approach
  }

  if(priors){ //if split files over subjects, just compute priors once
    for(i in 1:nparams){
      if(laplaceprior[i]==1) target+= priormod2 * double_exponential_lpdf(pow(fabs(rawpopmeans[i]) ,1+.1/((rawpopmeans[i]*100)^2+.1))|0,1);
    }
  }

  if(priors && !laplaceprioronly){ //if split files over subjects, just compute priors once
  for(i in 1:nparams){
    if(laplaceprior[i]==0) target+= priormod2 * normal_lpdf(rawpopmeans[i]|0,1);
  }
  
    if(nindvarying > 0){
      if(nindvarying >1) target+= priormod2 * normal_lpdf(sqrtpcov | 0, 1);
      target+= priormod2 * normal_lpdf(rawpopsdbase | ',gsub('normal(','',ctm$rawpopsdbase,fixed=TRUE),';
    }
  } //end pop priors section
  
  ',if(!gendata) 'if(intoverstates==0) target+= normal_lpdf(etaupdbasestates|0,1);','
  
  ',if(!gendata) 'target+= ll; \n','
  if(verbose > 0) print("lp = ", target());
}
',if(savemodel) paste0('
  generated quantities{
  vector[nparams] popmeans;
  vector[nindvarying] popsd; // = rep_vector(0,nparams);
  matrix[nindvarying,nindvarying] popcov;
  matrix[nparams,ntipred] linearTIPREDEFFECT;
',if(gendata) paste0('
  real ll = 0;
  vector[ndatapoints] llrow = rep_vector(0,ndatapoints);
  matrix[nlatentpop,nlatentpop] etacova[3,savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ycova[3,savescores ? ndatapoints : 0];
  vector[nlatentpop] etaa[3,savescores ? ndatapoints : 0];
  vector[nmanifest] ya[3,savescores ? ndatapoints : 0];
  vector[nmanifest] Ygen[ndatapoints];
  vector[nlatentpop*ndatapoints] etaupdbasestates = to_vector(normal_rng(rep_array(0.0,ndatapoints*nlatentpop),rep_array(1.0,ndatapoints*nlatentpop))); //sampled latent states posterior','
  ',subjectparaminit(pop=FALSE,smats=FALSE),'
  ',subjectparaminit(pop=TRUE,smats=FALSE)
  ,collapse=''),'

  {
    matrix[popcovn, nindvarying] x;
    if(nindvarying){
      for(ri in 1:rows(x)){
        x[ri,] = (rawpopcovchol * 
          to_vector(normal_rng(rawpopmeans[indvaryingindex],rep_vector(1,nindvarying))) )\';
      }
    }
    
    for(pi in 1:nparams){
      int found=0;
      int pr1;
      int pr2;
      real rawpoppar = rawpopmeans[pi];
      while(!found){ //currently seems useless, instead just references last match if multiple
        for(ri in 1:size(matsetup)){
          if(matsetup[ri,3]==pi && matsetup[ri,8]<=0) { //if a free parameter 
            pr1 = ri; 
            pr2=ri;// unless intoverpop, pop matrix row reference is simply current row
            found=1;
            if(intoverpop && matsetup[ri,5]) { //check if shifted
              for(ri2 in 1:size(matsetup)){ //check when state reference param of matsetup corresponds to row of t0means in current matsetup row
                if(matsetup[ri2,8]  && matsetup[ri2,3] == matsetup[ri,1] && 
                matsetup[ri2,3] > nlatent && matsetup[ri2,7] < 20) pr2 = ri2; //if param is dynamic and matches row (state ref) and is not in jacobian
                //print("ri = ",ri, " pr2 = ",pr2, " ri2 = ",ri2);
              }
            }
          }
        }
      }
        
      popmeans[pi] = tform(rawpoppar, matsetup[pr2,4], matvalues[pr2,2], matvalues[pr2,3], matvalues[pr2,4], matvalues[pr2,6] ); 
      if(matsetup[pr1,5]){ //if indvarying, transform random sample
        for(ri in 1:rows(x)){
          x[ri,matsetup[pr1,5]] = tform(x[ri,matsetup[pr1,5]],matsetup[pr2,4],matvalues[pr2,2],matvalues[pr2,3],matvalues[pr2,4],matvalues[pr2,6]);
        }
        x[,matsetup[pr1,5]] += rep_vector(-mean(x[,matsetup[pr1,5]]),rows(x));
      }
      if(ntipred > 0){
      for(tij in 1:ntipred){
        if(TIPREDEFFECTsetup[matsetup[pr1,3],tij] ==0){
          linearTIPREDEFFECT[matsetup[pr1,3],tij] = 0;
        } else {
        linearTIPREDEFFECT[matsetup[pr1,3],tij] = ( //tipred reference is from row pr1, tform reference from row pr2 in case of intoverpop
          tform(rawpoppar + TIPREDEFFECT[matsetup[pr1,3],tij] * .01, matsetup[pr2,4], matvalues[pr2,2], matvalues[pr2,3], matvalues[pr2,4], matvalues[pr2,6] ) -
          tform(rawpoppar - TIPREDEFFECT[matsetup[pr1,3],tij] * .01, matsetup[pr2,4], matvalues[pr2,2], matvalues[pr2,3], matvalues[pr2,4], matvalues[pr2,6] )
          ) /2 * 100;
        }
      }
    }
    } //end nparams loop
  
  if(nindvarying){
    popcov = crossprod(x) /(rows(x)-1);
    popsd = sqrt(diagonal(popcov));
  }
  }

',if(gendata) paste0('
{
  vector[nmanifest] Ygenbase[ndatapoints];
  Ygen = rep_array(rep_vector(99999,nmanifest),ndatapoints);
  for(mi in 1:nmanifest){
    if(manifesttype[mi]==0 || manifesttype[mi]==2) {
      Ygenbase[1:ndatapoints,mi] = normal_rng(rep_vector(0,ndatapoints),rep_vector(1,ndatapoints));
    }
    if(manifesttype[mi]==1){
      Ygenbase[1:ndatapoints,mi] =  uniform_rng(rep_vector(0,ndatapoints),rep_vector(1,ndatapoints));
    }
  }
{
',if(gendata) ukfilterfunc(ppchecking=TRUE),
  '

}}
',collapse=";\n"),'

}
'),collapse='')
}
  m <- writemodel()
  return(m)
}
