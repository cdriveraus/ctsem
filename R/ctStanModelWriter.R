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

simpleStateCheck <- function(x){   
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



ctModelStatesAndPARS <- function(ctm){ #replace latentName and par references with square bracket refs
  #detect state refs
  # 
  ln <- ctm$latentNames
  for(li in c(1:length(ln))){
    for(ri in grep(paste0('\\b(',ln[li],')\\b'),ctm$pars$param)){
      ctm$pars$param[ri] <- gsub(paste0('\\b(',ln[li],')\\b'),paste0('state[',li,']'),ctm$pars$param[ri])
    }
  }
  #expand pars
  # browser()
  ln <- ctm$pars$param[ctm$pars$matrix %in% 'PARS' & !is.na(ctm$pars$param)] #get extra pars
 
  for(li in seq_along(ln)){ #for every extra par
    parmatch <- which(ctm$pars$param %in% ln[li] & ctm$pars$matrix %in% 'PARS')
    for(ri in grep(paste0('\\b',ln[li],'\\b'),ctm$pars$param)){ #which rows contain the par
      if(!(ctm$pars$param[ri] == ln[li])){ #that are not the par itself #& ctm$pars$matrix[ri]=='PARS' #removed limitation of referencing within PARS matrices
      # print(ctm$pars$param[ri])
        ctm$pars$param[ri] <- gsub(paste0('\\b',ln[li],'\\b'), #replace with PARS reference...
          paste0('PARS[',ctm$pars$row[parmatch],',',ctm$pars$col[parmatch],']'),ctm$pars$param[ri])
      }
    }
  }
  return(ctm)
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
          yest<- eval(parse(text=formula.types$formula[i]))
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

    return(formula.types[which(formula.types$lsfit %in% min(formula.types$lsfit,na.rm=TRUE)),])
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
    
    #replaced with faster version
    # nctspec <- data.frame(merge(ctm$pars,df,by=0,all=TRUE,no.dups = FALSE))
    # nctspec <- nctspec[order(as.numeric(nctspec$Row.names)),]
    
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
    # m$pars$indvarying[m$pars$matrix %in% 'T0MEANS'] <- FALSE 
    
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
      # t0m$multiplier <- 1
      # t0m$meanscale <- 1
      # t0m$offset <- 0
      # t0m$inneroffset <- 0
      t0m$indvarying <- TRUE
      
      #new t0var 
      t0v <- m$pars[m$pars$matrix %in% 'T0VAR' & m$pars$row==1 & m$pars$col==1,,drop=FALSE]
      for(ri in 1:(m$n.latent+nindvaryingsmall)){
        for(ci in 1:(m$n.latent+nindvaryingsmall)){
          if(!(ri %in% t0mnotvarying && ci %in% t0mnotvarying)){
            t0v <- rbind(t0v,c('T0VAR',ri,ci,
              NA, #ifelse(ci > ri, NA, paste0('popcov_',ivnamesfull[ri],'_',ivnamesfull[ci])), #param
              0, #ifelse(ci > ri,0, NA), #value
              NA, #ifelse(ci == ri,1, 0), #transform,
              # NA, #ifelse(ri > m$n.latent || ci > m$n.latent, 1,m$pars$sdscale[(m$pars$matrix %in% 'T0VAR' & m$pars$row==ri & m$pars$col==ci)]), #multiplier (sdscale)
              # NA,#2,
              # NA,#ifelse(ci == ri,0, 0), #offset
              # NA,#0
              FALSE,1,rep(FALSE,m$n.TIpred)))
            m$pars <- m$pars[!(m$pars$matrix %in% 'T0VAR' & m$pars$row==ri & m$pars$col==ci),,drop=FALSE] #remove old t0var line
          }
        }}
      t0v=t0v[-1,,drop=FALSE] #remove initialisation row
      
      # #remove T0VAR lines where we will replace with rawpopcovsqrt
      # for(ri in 1:(m$n.latent+nindvaryingsmall)){
      #   for(ci in 1:(m$n.latent+nindvaryingsmall)){
      #     if(!(ri %in% t0mnotvarying && ci %in% t0mnotvarying)){
      #       m$pars[(m$pars$matrix %in% 'T0VAR' & m$pars$row==ri & m$pars$col==ci),]  <-
      #         c('T0VAR',ri,ci,NA,99999,0,99999,99999,99999,99999,FALSE,1,rep(FALSE,m$n.TIpred))
      #     }
      #   }}
      
      
      # #new drift -- find a way to get rid of this -- check jacobian?
      # drift <- m$pars[m$pars$matrix %in% 'DRIFT' & m$pars$row==1 & m$pars$col==1,]
      # for(ri in 1:(m$n.latent+nindvaryingsmall)){
      #   for(ci in 1:(m$n.latent+nindvaryingsmall)){
      #     newrow <- list('DRIFT',ri,ci,NA,ifelse(m$continuoustime,0,ifelse(ri==ci,1,0)),NA,FALSE,1)
      #     if(m$n.TIpred > 0) newrow<-c(newrow,rep(FALSE,m$n.TIpred))
      #     drift[nrow(drift)+1,] <- newrow
      #   }}
      # drift=drift[-1,,drop=FALSE]
      # drift=subset(drift,drift$row > m$n.latent | drift$col > m$n.latent)

      
      #reference new states
      for(ivi in ivnames){
        #jacobian pre compute
        # jcol = m$n.latent+match(ivi,ivnames)
        # jrow = m$pars$row[m$pars$param %in% ivi]
        # jtform = suppressWarnings(as.integer(m$pars$transform[m$pars$param %in% ivi]))
        # if(!is.na(jtform)){ #if using a basic transform not custom
        # newj = list(matrix='JAx',row=jrow,col=jcol,param=
        # }
        
        # m$pars$param[m$pars$param %in% ivi] <- tform( #transferring full transform state as per here not necessary, doesn't influence jacobian?
        #   param = paste0( 'state[',m$n.latent+match(ivi,ivnames),']'),
        #   transform = m$pars$transform[m$pars$param %in% ivi],
        #   multiplier = m$pars$multiplier[m$pars$param %in% ivi],
        #   meanscale = m$pars$meanscale[m$pars$param %in% ivi],
        #   offset = m$pars$offset[m$pars$param %in% ivi],
        #   inneroffset = m$pars$inneroffset[m$pars$param %in% ivi],singletext=TRUE)
        m$pars$indvarying[m$pars$param %in% ivi] <- FALSE
        m$pars[m$pars$param %in% ivi,paste0(m$TIpredNames,rep('_effect',m$n.TIpred))] <- FALSE
        m$pars$param[m$pars$param %in% ivi] <- sapply(which(m$pars$param %in% ivi), function(ri){
          gsub('param',paste0( 'state[',m$n.latent+match(ivi,ivnames),']'),m$pars$transform[ri]) 
        })
        m$pars$transform[m$pars$param %in% ivi] <-NA
        
      }
      
      # #collect transform and state reference together in param
      # for(ri in 1:nrow(m$pars)){
      #   if(grepl('[',m$pars$param[ri],fixed=TRUE)){
      #     m$pars$param[ri] <- gsub('param',m$pars$param[ri],m$pars$transform[ri])
      #   }
      # }
      
      # m$pars$indvarying  <-FALSE
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
  

  bcalcs=strsplit(bcalc,';')[[1]]
  bcalcs=gsub('\n','',bcalcs)
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
    
    scalcs2 = gsub(",", "; =", scalcs2, fixed = TRUE)
    scalcs2 = gsub("___rightsquarebracket___", "]", scalcs2, fixed = TRUE)
    scalcs2 = gsub("___leftsquarebracket___", "[", scalcs2, fixed = TRUE)
    scalcs2 = gsub("___comma___", ",", scalcs2, fixed = TRUE)
    
    scalcs2=gsub('\\.e(\\d+)<-','real e\\1 =',scalcs2)
    scalcs2=gsub('\\.e(\\d+)','e\\1',scalcs2)
    scalcs2=gsub(';',';  \n  ',scalcs2)
    scalcs2=gsub('\\{','',scalcs2) #remove leading curly brace
    
    ec=gsub('\\b(c)\\b\\(.*','',scalcs2) #retain just the simplifications (both lhs and rhs)
    scalcs2=gsub('.*\\b(c)\\b\\(','=',scalcs2) #retain just the rhs equations for the pars, insert first =

    scalcs2=gsub('\\}$','',scalcs2) #remove trailing curly ;
    scalcs2=gsub('\\)$',';',scalcs2) #replace trailing bracket close with ;

    scalcs2 = strsplit(scalcs2,'\n')[[1]]
    scalcs = paste0(bcalcs1,scalcs2,'\n')
    # cat(scalcs)
    out = paste0('    {\n  ',paste0(ec,collapse='  '),paste0(scalcs,collapse=' '),'  } \n  ',collapse=' ')
    return(out)
  }
}



ctStanModelCleanctspec <-  function(ctspec){ #clean ctspec structure
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
  ctspec$indvarying[ fixed | grepl('[', ctspec$param,fixed=TRUE) ] <- FALSE #remove indvarying for calcs also
  ctspec$param[ fixed ] <- NA
  ctspec$transform[ fixed] <- NA
  if(length(tieffects) > 0)  ctspec[fixed,tieffects] <- FALSE
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
  jacobian = c(J0=51,JAx=52,Jtd=53,Jy=54)
  asymptotic = c(asymCINT=21,asymDIFFUSION=22)
  extra <- c(DIFFUSIONcov=31)
  all <- c(base,jacobian,asymptotic, extra)
  mn <- list(base=base, jacobian=jacobian, asymptotic=asymptotic, extra=extra,all=all)
  mn$driftcint <- all[names(all) %in% c('DRIFT','CINT')]
  mn$diffusion <- all[names(all) %in% c('DIFFUSION','JAx')]
  mn$tdpred <- all[names(all) %in% c('TDPREDEFFECT','Jtd')]
  mn$measurement <- all[names(all) %in% c('LAMBDA','MANIFESTMEANS','MANIFESTVAR','Jy')]
  mn$t0 <- all[names(all) %in% c('T0MEANS','T0VAR','J0')]
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

  for(i in seq_along(c(mats$base,mats$jacobian))){
    m=c(mats$base,mats$jacobian)[i]
    nm=names(m)
    for(i in 1:nrow(ctspec)){ 
      if(ctspec$matrix[i] == nm) {
        when = 0 #default, static parameter
        parameter = 0 #for fixed values
        dynpar=FALSE
        copymatrix=0
        copyrow=0
        copycol=0
        simplestate <- FALSE
        if(!is.na(ctspec$param[i])){ #if non fixed parameter,
          # if(ctspec$row[i] %in% 6) 
          # 
          if(grepl('\\b(state)\\b\\[\\d+\\]',ctspec$param[i])){ #if state based 
            # 
            # print( ctspec[i,])
            simplestate <- simpleStateCheck(ctspec$param[i])
            if(simplestate){# && is.na(ctspec$transform[i])){ 
              cpx <- ctspec[i,]
              cpx$transform <- gsub('\\b(state)\\b\\[\\d+\\]','param',cpx$param)
              cpx$multiplier <- cpx$offset <- cpx$meanscale <- cpx$inneroffset <- NULL
              ctspec[i,] <- ctModelTransformsToNum(list(pars=cpx))$pars
              # print( ctspec[i,])
              if(is.na(suppressWarnings(as.integer(ctspec$transform[i])))) simplestate <- FALSE  #if couldn't convert to integer tform
            }
          }
          # 
          
          if(!simplestate && grepl('[',ctspec$param[i],fixed=TRUE)){ 
            
            if(grepl(paste0('^(', #if a direct matrix reference
              paste0('s',names(c(mats$base,mats$jacobian)),collapse='|'),
              ')\\[\\d+,\\s*\\d+\\]$'),ctspec$param[i])){ 
              copymatrix=sum(sapply(seq_along(c(mats$base,mats$jacobian)), function(x){
                ifelse(grepl(paste0('s',names(c(mats$base,mats$jacobian)[x])),ctspec$param[i]),x,0)
                }))
              copymatrix = c(mats$base,mats$jacobian)[copymatrix]
              index = as.numeric(strsplit(gsub('\\]','',
                gsub('\\[','',
                  regmatches(ctspec$param[i],regexpr('\\[\\d+,\\s*\\d+\\]$',ctspec$param[i]))
                )),split = ',\\s*')[[1]])
              copyrow=index[1]
              copycol=index[2]
              when = -999 #rely on original when rather than duplicates
              indvar <- 0
            } else { #if not a direct matrix reference
              
              calcs <- c(calcs, paste0(ctspec$matrix[i],'[',ctspec$row[i], ', ', ctspec$col[i],'] = ',
                ctspec$param[i]))
              when= -999 #never use this row of ctspec during calculations (calc is extracted and placed elsewhere)
              indvar=0
              dynpar=TRUE
            }
          }#end calculation / copy parameters
          
          else if(simplestate){
            # 
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
            if(nm %in% c('J0',names(mats$t0))) when = 1
            dynpar=TRUE
          }
          
          else if(!grepl('[',ctspec$param[i],fixed=TRUE)){ #if non calculation parameter
            if(i > 1 && any(ctspec$param[1:(i-1)] %in% ctspec$param[i])){ #and after row 1, check for duplication
              parameter <- matsetup[,'param'][ match(ctspec$param[i], matsetup$parname) ] #find which freepar corresponds to duplicate
              indvar <- ifelse(ctspec$indvarying[i],  matsetup[,'indvarying'][ match(ctspec$param[i], matsetup$parname) ],0)#and which indvar corresponds to duplicate
              
            } else { #if not duplicated
              freeparcounter <- freeparcounter + 1
              TIPREDEFFECTsetup <- rbind(TIPREDEFFECTsetup,rep(0,ncol(TIPREDEFFECTsetup))) #add an extra row...
              if(ctspec$indvarying[i]) {
                indvaryingcounter <- indvaryingcounter + 1
                indvar <- indvaryingcounter
              }
              if(!ctspec$indvarying[i]) indvar <- 0
              freepar <- freeparcounter
              parameter <- freepar
              if(is.na(suppressWarnings(as.integer(ctspec$transform[i])))) { #extra tform needed
                extratformcounter <- extratformcounter + 1
                # extratforms <- paste0(extratforms,'if(transform==',-10-extratformcounter,') out = ',
                #   ctspec$offset[i],' + ',ctspec$multiplier[i],' * (inneroffset + ',
                #   gsub('param', paste0('param * ',ctspec$meanscale[i]),ctspec$transform[i]),');')
                extratforms <- paste0(extratforms,'if(transform==',-10-extratformcounter,') out = ',
                  ctspec$transform[i],';')
                ctspec$transform[i] <- -10-extratformcounter
              }
              if(n.TIpred > 0) {
                TIPREDEFFECTsetup[freepar,][ ctspec[i,paste0(ctm$TIpredNames,'_effect')]==TRUE ] <- 
                  tipredcounter: (tipredcounter + sum(as.integer(suppressWarnings(ctspec[i,paste0(ctm$TIpredNames,'_effect')]))) -1)
                tipredcounter<- tipredcounter + sum(as.integer(suppressWarnings(ctspec[i,paste0(ctm$TIpredNames,'_effect')])))
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
          tipred=ifelse(any(TIPREDEFFECTsetup[freepar,] > 0) && !dynpar, 1, 0), #tipredvarying
          matrix=m, #matrix reference
          when=when,#when to compute
          copymatrix=copymatrix,copyrow=copyrow,copycol=copycol,
          stringsAsFactors = FALSE)

        if(is.null(matsetup)) matsetup <- mdatnew else matsetup<-rbind(matsetup,mdatnew)
        mvalnew<-ctspec[i,c('value','multiplier','meanscale','offset','sdscale','inneroffset'),drop=FALSE]
        # mvalnew <- lapply(mvalnew,function(x) ifelse(is.null(x),NA,x)) #not sure why this is necessary...
        if(is.null(mval)) mval <- mvalnew else mval<-rbind(mval,mvalnew)
      }
    }
    if(!is.null(mval)) mval[is.na(mval)] <- 99999 else mval<-array(0,dim=c(0,6))
    
  }

  matrixdims <- array(0L,dim=c(max(mats$all),2))
  sapply(mats$all, function(m) {
    matrixdims[m,]<<-as.integer(c(max(c(0,matsetup[matsetup[,'matrix'] %in% m,'row'])),
      max(c(0,matsetup[matsetup[,'matrix'] %in% m,'col']))))
  })
  #special matrix adjustments
  matrixdims[mats$all[names(mats$all ) == 'asymDIFFUSION'],] <- 
    matrixdims[mats$all[names(mats$all ) == 'DIFFUSIONcov'],] <- 
    matrixdims[mats$all[names(mats$all ) == 'DIFFUSION'],]
  matrixdims[mats$all[names(mats$all ) == 'asymCINT'],] <- 
    matrixdims[mats$all[names(mats$all ) == 'CINT'],]

  matsetup[,!colnames(matsetup) %in% 'parname'] <- lapply(matsetup[,!colnames(matsetup) %in% 'parname'],as.integer)
  matvalues <- data.frame(apply(mval,2,as.numeric,.drop=FALSE),stringsAsFactors = FALSE  )

  #copy handling
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
  # browser()
  
  for(i in 1:nrow(matsetup)){
    if(matsetup$copyrow[i] > 0){
      matsetup[i,c('param','transform','indvarying','tipred','when')] <- 
        matsetup[matsetup$copyrow[i],c('param','transform','indvarying','tipred','when')]
    }
  }
  
  # matsetup$copyrow[unique(matsetup$copyrow[matsetup$copyrow > 0])] <- 
  #   -unique(matsetup$copyrow[matsetup$copyrow > 0]) #set rows to be copied from to neg of their own row num
  
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



ctStanCalcsList <- function(ctm){  #extract any calcs from model and ensure spec is correct
  temp <- ctm$modelmats$calcs
  names(temp) <- NULL
  mats<-ctStanMatricesList()
  for(mati in unique(ctm$pars$matrix)){ #append s to matrices
    temp <- gsub(paste0('\\b',mati),paste0('s',mati),temp)
  }

  calcs <- lapply(c(list(PARS=c(PARS=10)),mats[!names(mats) %in% c('base','all','jacobian')]), function(mlist) { #add custom calculations to correct list of matrices
    out <- temp[unlist(sapply(temp, function(y) any(sapply(names(mlist), function(mli) 
      grepl(mli,gsub('=.*','',y))
      ))))]
    return(out)
  })
  ctm$modelmats$calcs <- calcs
  return(ctm)
}

# duplicateRefs <- function(ctm){
#   
#   ml <- ctStanMatricesList(ctm)
#   
#   upd=lapply(ml$all, function(m){
#     dups <-
#     
#     mats[!names(mats) %in% c('base','all','jacobian')], function(mlist) {
#   
#   calcs$diffusion2 <- c(calcs$driftcint,calcs$diffusion) #create combined vec for duplicate detection
# calcs <- lapply(calcs,function(mli){ #create matrix refs for duplicates
#   mli <- gsub(' ','',mli)
#   rhs <- gsub('.*=','',mli)
#   lhs <- gsub('=.*','',mli)
#   
#   dup <- duplicated(rhs)
#   for(i in seq_along(rhs)){
#     if(dup[i]) mli[i] <- gsub('=.*',paste0('=',lhs[match(x = rhs[i],rhs)]),mli[i])
#   }
#   return(mli)
# })
# calcs$diffusion <- tail(calcs$diffusion2,length(calcs$diffusion)) #replace with portion of combined vec
# calcs$diffusion2 <- NULL
# 
# }



jacobianelements <- function(J, when, ntdpred,matsetup,mats,textadd=NA, 
  remove='', #c('drift','simplestate', 'state','fixed','lambda'),
  returndriftonly=FALSE,
  returnlambdaonly=FALSE){ 
  out=c()
  # if(returndriftonly) 
  if(!(when == 3 && ntdpred ==0)){
    lambdaonly <- driftonly <-   simplestatesonly <- c()
    if( !returnlambdaonly){
      for(ci in 1:ncol(J)){
        for(ri in 1:nrow(J)){ #check for standard drift reference
          if('drift' %in% remove && J[ri,ci] == paste0('sDRIFT[',ri,',',ci,']')) {
            chk=TRUE
          } else {
            if('drift' %in% remove && J[ri,ci] == mats$DRIFT[ri,ci] && is.na(suppressWarnings(as.numeric(J[ri,ci])))){
              chk=TRUE
            } else chk <- FALSE
          }
          driftonly <- c(driftonly, chk) 
          simplestatesonly <- c(simplestatesonly, grepl('^\\b(state)\\b\\[\\d+\\]$',J[ri,ci]))
        }
      }
    }
    chk <- c()
    if('lambda' %in% remove && !returndriftonly){
      for(ci in 1:nrow(mats$T0MEANS)){
        for(ri in 1:nrow(mats$LAMBDA)){ #check for standard lambda reference
          if(ci <= ncol(mats$LAMBDA) && J[ri,ci] == paste0('sLAMBDA[',ri,',',ci,']')) {
            chk=TRUE
          } else {
            if(ci <= ncol(mats$LAMBDA) && 
                J[ri,ci] == mats$LAMBDA[ri,ci] && 
                is.na(suppressWarnings(as.numeric(J[ri,ci])))){
              chk=TRUE
            } else chk <- FALSE
          }
          lambdaonly <- c(lambdaonly, chk) 
        }
      }
    }
    
    out = J
    drop<-as.integer(c())
    if('fixed' %in% remove) drop=which(!is.na(suppressWarnings(as.numeric(out)))) #dropping single values
    if('lambda' %in% remove) drop = unique(c(drop,which(lambdaonly)))
    if('drift' %in% remove) drop = unique(c(drop,which(driftonly)))
    if('simplestate' %in% remove && length(simplestatesonly) > 0) drop = unique(c(drop,which(simplestatesonly)))
    
    if(all(!is.na(textadd))) out=paste0(textadd, J,';\n')
    if(length(drop) > 0) out=out[-drop] #dropping selected terms
    
    if('state' %in% remove) out = out[!grepl('state[',out,fixed=TRUE)] #else out = out[(grepl('state[',out,fixed=TRUE))]
    
    if(all(!is.na(textadd))) out=paste0(out,collapse='')
    if(when==3 && ntdpred==0) out <- c() #no need for extra jacobian lines if no predictors!
    if(returndriftonly) out <- matrix(as.integer(as.logical(driftonly)),nrow(J),ncol(J))
    if(returnlambdaonly) {
      out <- matrix(as.integer(as.logical(lambdaonly)),nrow(mats$LAMBDA),nrow(mats$T0MEANS))
      
    }
  }
  return(out)
}

ctStanModelWriter <- function(ctm, gendata, extratforms,matsetup,simplify=TRUE){
  #if arguments change make sure to change ctStanFit !
  
  mats <- ctStanMatricesList()
  # #check when / if PARS needs to be computed
  # for(mlist in names(mats[-1])){
  #   if(any(unlist(lapply(ctm$calcs[[mlist]], function(m) grepl('sPARS',m))))) mats[[mlist]]=c(mats[[mlist]],'PARS')
  # }
  
  
  #adjust diffusion calcs to diffusionsqrt
  # ctm$calcs$diffusion <- gsub('sDIFFUSION','sDIFFUSIONsqrt',ctm$calcs$diffusion)
  # intoverpopdynamiccalcs <- gsub('sDIFFUSION','sDIFFUSIONsqrt',intoverpopdynamiccalcs)
  
  #save calcs without intoverpop for param intitialisation
  
  
  #consider allowing copies across different 'when' states by dropping mlist from the final loop
  simplestatedependencies <- function(when, mlist,basemats=FALSE) {
    
    paste0('statetf[whichequals(whenvecs[',when[when!=0],'],0,0)] = 
         parvectform(size(whichequals(whenvecs[',when,'],0,0)),state, ',when,
    ', matsetup, matvalues, si, subindices, whenvecs[',when,']);
    
    ',matcalcs('si',when=when, mlist,basemats=basemats))
    
  }

  
  
  finiteJ<-function(){
    paste0('
    {
    int zeroint[1];
    vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(sJAxfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0)  state[statei] += Jstep;
      
        statetf[whichequals(whenvecs[2],0,0)] = 
          parvectform(size(whichequals(whenvecs[2],0,0)),state, 2, matsetup, matvalues, si, subindices, whenvecs[2]);

        ',matcalcs('si',when=2, c(PARS=10),basemats=FALSE),' //initialise PARS first, and simple PARS before complex PARS
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n\n ',collapse=' '),simplify),'
      
        ',matcalcs('si',when=2, mats$driftcint,basemats=FALSE),'
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$driftcint,';\n\n ',collapse=' '),simplify),'
      
      if(statei > 0) {
        sJAx[1:nlatent,statei] =  sDRIFT * state[1:nlatent] + sCINT[,1]; //compute new change
         if(verbose>1) print("sJAx ",sJAx);
      }
      if(statei== 0 && size(sJAxfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
        base = sDRIFT * state[1:nlatent] + sCINT[,1];
        if(verbose>1) print("base = ",base,"    sjaxinit= ",sJAx);
        for(fi in sJAxfinite){
          sJAx[1:nlatent,fi] = (sJAx[1:nlatent,fi] - base) / Jstep; //new - baseline change divided by stepsize
        }
      }
    }
    if(verbose>1) print("sJAx ",sJAx);
    }
    ')
  }
  
  finiteJy<-function(){
    paste0('
  {
    int zeroint[1];
    vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(sJyfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0 && (savescores + intoverstates) > 0)  state[statei] += Jstep;
      
            
        statetf[whichequals(whenvecs[4],0,0)] = 
          parvectform( size(whichequals(whenvecs[4],0,0)), state, 4, matsetup, matvalues, si, subindices, whenvecs[4]);
          
        ',matcalcs('si',when=4, c(PARS=10),basemats=FALSE),' //initialise PARS first, and simple PARS before complex PARS
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n\n ',collapse=' '),simplify),'
      
        ',matcalcs('si',when=4, mats$measurement,basemats=FALSE),'
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$measurement,';\n\n ',collapse=' '),simplify),'
        
      if(statei > 0 && (savescores + intoverstates) > 0) {
        sJy[o,statei] =  sLAMBDA[o] * state[1:nlatent] + sMANIFESTMEANS[o,1]; //compute new change
        sJy[o1,statei] = to_vector(inv_logit(to_array_1d(sJy[o1,statei])));
         if(verbose>1) print("sJy ",sJy);
      }
      if(statei==0){
        syprior[o] = sLAMBDA[o] * state[1:nlatent] + sMANIFESTMEANS[o,1];
        syprior[o1] = to_vector(inv_logit(to_array_1d( syprior[o1] )));
        if(size(sJyfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
          if(verbose>1) print("syprior = ",syprior,"    sJyinit= ",sJy);
          for(fi in sJyfinite){
            sJy[o,fi] = (sJy[o,fi] - syprior[o]) / Jstep; //new - baseline change divided by stepsize
          }
        }
      }
    }
    if(verbose>1) print("sJy ",sJy);
  }')
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
  
  matrix[nlatentpop,nlatentpop] Je[savescores ? ndatapoints : 1]; //time evolved jacobian, saved for smoother

  vector[nlatentpop] state = rep_vector(-999,nlatentpop); 
  vector[nlatentpop] statetf;
  matrix[nlatentpop,nlatentpop] sJAx; //Jacobian for drift
  matrix[nlatentpop,nlatentpop] sJ0; //Jacobian for t0
  matrix[nlatentpop,nlatentpop] sJtd;//diag_matrix(rep_vector(1),nlatentpop); //Jacobian for nltdpredeffect
  matrix[ nmanifest,nlatentpop] sJy;//Jacobian for measurement 
  
  ',if(!is.null(ctm$taylorheun) && ctm$taylorheun==1) paste0('
  matrix[taylorheun ? nlatentpop : 0, taylorheun ? nlatentpop : 0] Kth = rep_matrix(0.0,taylorheun ? nlatentpop : 0, taylorheun ? nlatentpop : 0); 
  matrix[taylorheun ? nlatentpop : 0, taylorheun ? nlatentpop : 0] Mth = Kth;'),'

  //linear continuous time calcs
  matrix[nlatent+1,nlatent+1] discreteDRIFT;
  matrix[nlatent,nlatent] discreteDIFFUSION = rep_matrix(0.0,nlatent,nlatent);

  
  vector[nparams] rawindparams = rawpopmeans;
  vector[nparams] indparams;

  //dynamic system matrices
  ',subjectparaminit(pop=FALSE,smats=TRUE),'
  
  sasymDIFFUSION = rep_matrix(0,nlatent,nlatent); //in case of derrindices need to init
  sDIFFUSIONcov = rep_matrix(0,nlatent,nlatent);

  for(si in 0:max(subject)){
  for(rowi in 1:ndatapoints){
    if( (rowi==1 && si==0) ||
      (dokalman && dokalmanrows[rowi] && subject[rowi]==si) ){ //if doing this row for this subject
    
    int full = (savescores==1 || si ==0);
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
    
    if(prevrow != 0) T0check = (si==subject[prevrow]) ? (T0check+1) : 0; //if same subject, add one, else zero
    if(T0check > 0){
      dt = time[rowi] - time[prevrow];
      dtchange = dt!=prevdt; 
      prevdt = dt; //update previous dt store after checking for change
      //prevtime = time[rowi];
    }

    if(savescores && prevrow!=0) Je[rowi,,] = Je[prevrow,,];
    
    if(T0check == 0) { // calculate initial matrices if this is first row for si
  
  rawindparams=rawpopmeans;
  
  if(si > 0 && nindvarying > 0 && intoverpop==0)  rawindparams[indvaryingindex] += rawpopc[2] * baseindparams[si];

  if(si > 0 &&  ntieffects > 0){
  if(nmissingtipreds > 0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipreds[si]\';
    
    if(nmissingtipreds==0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipredsdata[si]\';
  }

  indparams[whichequals(whenvecp[si ? 2 : 1], 0, 0)]= 
    parvectform(size(whichequals(whenvecp[si ? 2 : 1], 0, 0)),rawindparams, 
    0, matsetup, matvalues, si, subindices, whenvecp[si ? 2 : 1]);
     
  if(whenmat[1, 5] >= (si ? 1 : 0)) sT0MEANS = 
    mcalc(sT0MEANS, indparams, statetf, {0}, 1, matsetup, matvalues, si, subindices); // base t0means to init
      
  for(li in 1:nlatentpop) if(!is_nan(sT0MEANS[li,1])) state[li] = sT0MEANS[li,1]; //in case of t0 dependencies, may have missingness
  
  statetf[whichequals(whenvecs[1],0,0)] = parvectform(size(whichequals(whenvecs[1],0,0)),state, 1,
    matsetup, matvalues, si, subindices, whenvecs[1]);   
    
  ',matcalcs('si',when=0:1, c(PARS=10),basemats=TRUE),' //initialise simple PARS then do complex PARS
  ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n\n ',collapse=' '),simplify),'
    
  ',matcalcs('si',when=0:1, mats$t0,basemats=TRUE),'
      
  ',simplifystanfunction(paste0(ctm$modelmats$calcs$t0,';\n\n ',collapse=' '),simplify),'
    for(li in 1:nlatentpop) if(is_nan(state[li])) state[li] = sT0MEANS[li,1]; //finish updating state
    
    //init other system matrices (already done PARS, redo t0means in case of PARS dependencies...)
    ',paste0(sapply(1:length(c(mats$base[!mats$base %in% c(1,10)],mats$jacobian)),function(x) {
          y=c(mats$base[!mats$base %in% c(1,10)],mats$jacobian)[x]
          paste0(
          'if(whenmat[',y,',5] || si==0) ',matcalcs('si',when=0,y,basemats=TRUE))
          }),collapse='  '),'
    
    
  if(si <= (subindices[8] ? nsubjects : 0)) {
   if(intoverpop && nindvarying > 0) sT0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopc[1];
    sT0VAR = makesym(sdcovsqrt2cov(sT0VAR,choleskymats),verbose,1); 

    if(intoverpop && nindvarying > 0){ //adjust cov matrix for transforms
      for(ri in 1:size(matsetup)){
        if(matsetup[ri,7]==1){ //if t0means
          if(matsetup[ri,5]) { //and indvarying
            sT0VAR[matsetup[ri,1], ] = sT0VAR[matsetup[ri,1], ] * matvalues[ri,2] * matvalues[ri,3]* matvalues[ri,5]; //multiplier meanscale sdscale
            sT0VAR[, matsetup[ri,1] ] = sT0VAR[, matsetup[ri,1] ] * matvalues[ri,2] * matvalues[ri,3]* matvalues[ri,5]; //multiplier meanscale sdscale
          }
        }
      }
    }
  }
    
    if(verbose>1) print("nl t0var = ", sT0VAR, "   sJ0 = ", sJ0);
    etacov = quad_form(sT0VAR, sJ0\'); //probably unneeded, inefficient

    } //end T0 matrices
    
if(verbose > 1) print ("below t0 row ", rowi);

      if(si==0 || (T0check>0)){ //for init or subsequent time steps when observations exist
        vector[nlatent] base;
        real intstepi = 0;
        
        dtsmall = dt / ceil(dt / maxtimestep);
        
        while(intstepi < (dt-1e-10)){
          intstepi = intstepi + dtsmall;
          ',#simplestatedependencies(when=2,mlist=c(mats$driftcint,mats$diffusion,mats$jacobian[2])),' #now done inside finite diff loop
      finiteJ(),'
      
      ',matcalcs('si',when=2, mats$diffusion,basemats=FALSE),'
      ',simplifystanfunction(paste0(ctm$modelmats$calcs$diffusion,';\n\n ',collapse=' '),simplify),'
      if(si==0 || whenmat[4,2] || (T0check==1 && whenmat[4,5])) sDIFFUSIONcov[derrind,derrind] = sdcovsqrt2cov(sDIFFUSION[derrind,derrind],choleskymats);
      
        if(continuoustime){
          if(taylorheun==0){
            if(si==0 || dtchange==1 || statedep[3]||statedep[4]||statedep[7] || 
              (T0check == 1 && (subindices[3] + subindices[4] + subindices[7]) > 0)){
                
              ',if(1==99) 'if(difftype==2 && (statedep[3]||statedep[4])){
                matrix[ndiffusion*2,ndiffusion*2] ebA;
                matrix[ndiffusion*2,ndiffusion*2] bA;
            
                bA[1:ndiffusion,1:ndiffusion] = -sJAx[derrind,derrind];
                bA[1:ndiffusion,(1+ndiffusion):(ndiffusion*2)] = sDIFFUSIONcov[derrind,derrind];
                bA[(1+ndiffusion):(ndiffusion*2),(1+ndiffusion):(ndiffusion*2)] = sJAx[derrind,derrind]\';
                bA[(1+ndiffusion):(ndiffusion*2),1:ndiffusion] = rep_matrix(0,ndiffusion,ndiffusion);
                
                ebA = expm2(bA * dtsmall);
                Je[savescores ? rowi : 1,derrind,derrind] =  ebA[(1+ndiffusion):(ndiffusion*2),(1+ndiffusion):(ndiffusion*2)]\';
                discreteDIFFUSION[derrind,derrind] = Je[savescores ? rowi : 1,derrind,derrind] * ebA[1:ndiffusion,(1+ndiffusion):(ndiffusion*2)];
                discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec\') * dtsmall);
                if(ndiffusion < nlatent && savescores) Je[savescores ? rowi : 1,,] =  expm2(sJAx * dtsmall);
              }
                
              //if(difftype==1){
              //matrix[nlatent,nlatent] V = sDIFFUSIONcov-quad_form(sDIFFUSIONcov, Je[savescores ? rowi : 1,,]\');
              //discreteDIFFUSION = solvesyl(sJAx[1:nlatent,1:nlatent],-V,discreteDIFFUSION, rep_array(nlatent,1));
              //}
              ','
              
              if(difftype==0 || (statedep[3]==0 && statedep[4]==0)){
                discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec\') * dtsmall);
                Je[savescores ? rowi : 1,,] =  expm2(sJAx * dtsmall);
               
                if(si==0 || statedep[3]||statedep[4]|| (T0check==1 && (whenmat[4,5] || whenmat[3,5]))){
                sasymDIFFUSION[derrind,derrind] = to_matrix(  -sqkron_sumii(sJAx[derrind,derrind]) \\ 
                  to_vector(sDIFFUSIONcov[derrind,derrind]), ndiffusion,ndiffusion);
                }
                discreteDIFFUSION[derrind,derrind] =  sasymDIFFUSION[derrind,derrind] - 
                  quad_form( sasymDIFFUSION[derrind,derrind], Je[savescores ? rowi : 1, derrind,derrind]\' );
              }
            } 
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent]; // ???compute before new diffusion calcs
            if(intoverstates==1 || savescores==1){
              etacov = quad_form(etacov, Je[savescores ? rowi : 1,,]\');
              etacov[derrind,derrind] += discreteDIFFUSION[derrind,derrind]; 
            }
          }
            
           ',if(!is.null(ctm$taylorheun) && ctm$taylorheun==1) paste0('if(taylorheun==1){
                if(dtchange==1 || statedep[3]||statedep[4] || 
                (T0check == 1 && (subindices[3] + subindices[4]) > 0)){
                  Kth = (add_diag(-sJAx * (dtsmall /2),1) );
                  Mth = Kth \ (add_diag(sJAx * (dtsmall /2),1) );
                }
                state[1:nlatent] += Kth[1:nlatent,1:nlatent] \
                (sDRIFT[1:nlatent,1:nlatent] * state[1:nlatent] + sCINT[1:nlatent,1]) * dtsmall;
                etacov = quad_form_sym((etacov), Mth\');
      etacov[derrind,derrind] += (Kth[derrind,derrind] \\ sDIFFUSIONcov[derrind,derrind] / Kth[derrind,derrind]\') * dtsmall;
              '),'
              
            if(intstepi >= (dt-1e-10) && savescores) Je[rowi,,] = expm2(sJAx * dt); //save approximate exponentiated jacobian for smoothing
          }
  
          if(continuoustime==0){ 
            Je[savescores ? rowi : 1,,] = sJAx;
            if(intoverstates==1 || savescores==1){
              etacov = quad_form(etacov, sJAx\');
              etacov[ derrind, derrind ] += tcrossprod(sDIFFUSION[ derrind, derrind ]); 
            }
            discreteDIFFUSION=sDIFFUSIONcov;
            discreteDRIFT=append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec\');
            discreteDRIFT[nlatent+1,nlatent+1] = 1;
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent];
            
          }
        }
      } // end non linear time update
    
    if(ntdpred > 0) {
      int nonzerotdpred = 0;
      for(tdi in 1:ntdpred) if(tdpreds[rowi,tdi] != 0.0) nonzerotdpred = 1;
      if(nonzerotdpred){
      
        statetf[whichequals(whenvecs[3],0,0)] = 
          parvectform( size(whichequals(whenvecs[3],0,0)), state, 3, matsetup, matvalues, si, subindices, whenvecs[3]);
          
        ',matcalcs('si',when=3, c(PARS=10),basemats=FALSE),' //initialise PARS first, and simple PARS before complex PARS
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n\n ',collapse=' '),simplify),'
      
        ',matcalcs('si',when=3, mats$tdpred,basemats=FALSE),'
        ',simplifystanfunction(paste0(ctm$modelmats$calcs$tdpred,';\n\n ',collapse=' '),simplify),'

        state[1:nlatent] +=   (sTDPREDEFFECT * tdpreds[rowi]); //tdpred effect only influences at observed time point','
        if(statedep[9]) etacov = quad_form(etacov,sJtd\'); //could be optimized
      }
    }//end nonlinear tdpred

  if(intoverstates==0){
    if(T0check==0) state += cholesky_decompose(sT0VAR) * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)];
    if(T0check>0) state[derrind] +=  cholesky_decompose(makesym(discreteDIFFUSION[derrind,derrind],verbose,1)) * 
      (etaupdbasestates[(1+(rowi-1)*nlatentpop):(nlatent+(rowi-1)*nlatentpop)])[derrind];
  }

if(verbose > 1){
  print("etaprior = ", state);
  print("etapriorcov = ", etacov);
}

',if(!ppchecking) 'if(savescores){
  etacova[1,rowi] = etacov; 
  etaa[1,rowi] = state;
}','

 if ((si==0 || nobs_y[rowi] > 0 || savescores)){ //do this section for 0th subject as well to init matrices
    
      ',finiteJy(),'
 
   if(si > 0 && dokalmanrows[rowi] ==1){   //if not just inits...

      if(intoverstates==1 || savescores==1) { //classic kalman
        ycov[o,o] = quad_form(etacov, sJy[o,]\'); // + sMANIFESTVAR[o,o]; shifted measurement error down
        for(wi in 1:nmanifest){ 
          if(Y[rowi,wi] != 99999 || savescores==1) ycov[wi,wi] += square(sMANIFESTVAR[wi,wi]);
          if(manifesttype[wi]==1 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += fabs((syprior[wi] - 1) .* (syprior[wi]));
          if(manifesttype[wi]==2 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += square(fabs((syprior[wi] - round(syprior[wi])))); 
        }
      }
        
      if(intoverstates==0) { //sampled states
        if(ncont_y[rowi] > 0) for(mi in o0) ypriorcov_sqrt[mi,mi] = sqrt(sMANIFESTVAR[mi,mi]);
      }
        
     
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
        
         if(verbose > 1) print("before K rowi =",rowi, "  si =", si, "  state =",state, "  statetf = ", statetf, "  etacov ",etacov,
          " indparams = ", indparams,
            "  syprior[o] =",syprior[o],"  ycov[o,o] ",ycov[o,o], 
            "  sPARS = ", sPARS, 
            "  sDRIFT =", sDRIFT, " sDIFFUSION =", sDIFFUSION, " sCINT =", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
            "  sT0VAR", sT0VAR,  " sT0MEANS ", sT0MEANS, "sLAMBDA = ", sLAMBDA, "  sJy = ",sJy,
            " discreteDRIFT = ", discreteDRIFT, "  discreteDIFFUSION ", discreteDIFFUSION, "  sasymDIFFUSION ", sasymDIFFUSION, 
            " DIFFUSIONcov = ", sDIFFUSIONcov,
            "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
         
        K[,od] = etacov * sJy[od,]\' / makesym(ycov[od,od],verbose,1);// * multiply_lower_tri_self_transpose(ycovi\');// ycov[od,od]; 
        etacov += -K[,od] * sJy[od,] * etacov;
        state +=  (K[,od] * err[od]);
      }
      
      if(savescores==1) {
        ya[1,rowi] = syprior[o];
        etaa[2,rowi] = state;
        ycova[1,rowi] = ycov;
        etacova[2,rowi] = etacov;
        ycova[2,rowi] = quad_form(etacov, sJy\');
        for(wi in 1:nmanifest) ycova[2,rowi,wi,wi] += square(sMANIFESTVAR[wi,wi]);
        ya[2,rowi] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * state[1:nlatent];
      }
      
      
      if(verbose > 1) print(" After K rowi =",rowi, "  si =", si, "  state =",state,"  etacov ",etacov,"  K[,o] ",K[,o]);
        
  //likelihood stuff
      if(nbinary_y[rowi] > 0) llrow[rowi] += sum(log(Y[rowi,o1d] .* (syprior[o1d]) + (1-Y[rowi,o1d]) .* (1-syprior[o1d]))); 

      if(size(o0d) > 0 && (llsinglerow==0 || llsinglerow == rowi)){
        if(intoverstates==1) ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(makesym(ycov[o0d,o0d],verbose,1));
         llrow[rowi] +=  multi_normal_cholesky_lpdf(Y[rowi,o0d] | syprior[o0d], ypriorcov_sqrt[o0d,o0d]);
         //errtrans[counter:(counter + ncont_y[rowi]-1)] = 
           //mdivide_left_tri_low(ypriorcov_sqrt[o0d,o0d], err[o0d]); //transform pred errors to standard normal dist and collect
         //ll+= -sum(log(diagonal(ypriorcov_sqrt[o0d,o0d]))); //account for transformation of scale in loglik
         //counter += ncont_y[rowi];
      }
      
    }//end si > 0 nobs > 0 section
  } // end measurement init loop and dokalmanrows section here to collect matrices
    
       // store system matrices
       
  if(si <= ((subindices[3] + subindices[7])  ? nsubjects : 0)){
    if(continuoustime==1) sasymCINT[,1] =  -sDRIFT[1:nlatent,1:nlatent] \\ sCINT[ ,1 ];
    if(continuoustime==0) sasymCINT[,1] =  add_diag(-sDRIFT[1:nlatent,1:nlatent],1) \\ sCINT[,1 ];
  }
  
  if(si <= ((subindices[3] + subindices[7])  ? nsubjects : 0) && !continuoustime) sasymDIFFUSION[ derrind, derrind ] = 
    to_matrix( (add_diag( 
      -sqkron_prod(sDRIFT[ derrind, derrind ], sDRIFT[ derrind, derrind ]),1)) \\  
      to_vector(sDIFFUSIONcov[ derrind, derrind ]), ndiffusion, ndiffusion);
    
  if(savesubjectmatrices && si > 0 && (rowi==ndatapoints || subject[rowi+1] != subject[rowi])){',
  paste0(collectsubmats(popmats=FALSE),collapse=' '),'
  }
  if(si == 0){
',paste0(collectsubmats(popmats=TRUE),collapse=' '),'
  }

  
  ',if(!ppchecking) paste0('if(si > 0 && savescores && (rowi==ndatapoints || subject[rowi+1] != subject[rowi])){ //at subjects last datapoint, smooth
    int sri = rowi;
    while(sri>0 && subject[sri]==si){
      if(sri==rowi) {
        etaa[3,sri]=etaa[2,sri];
        etacova[3,sri]=etacova[2,sri];
      } else{
        matrix[nlatentpop,nlatentpop] smoother;
        smoother = etacova[2,sri] * Je[sri+1,,]\' / makesym(etacova[1,sri+1],verbose,1);
        etaa[3,sri]= etaa[2,sri] + smoother * (etaa[3,sri+1] - etaa[1,sri+1]);
        etacova[3,sri]= etacova[2,sri] + smoother * ( etacova[3,sri+1] - etacova[1,sri+1]) * smoother\';

      }
      state=etaa[3,sri];
      
      ',finiteJy(),'
      
      ya[3,sri] = syprior;
      ycova[3,sri] = quad_form(etacova[3,sri], sJy\'); 
      for(wi in 1:nmanifest){
        ycova[3,sri,wi,wi] += square(sMANIFESTVAR[wi,wi]);
        if(manifesttype[wi]==1) ycova[3,sri,wi,wi] += fabs((ya[3,sri,wi] - 1) * (ya[3,sri,wi]));
        if(manifesttype[wi]==2) ycova[3,sri,wi,wi] += square(fabs((ya[3,sri,wi] - round(ya[3,sri,wi])))); 
      }
      sri += -1;
      while(sri > 0 && dokalmanrows[sri]==0) sri+= -1; //skip rows if requested
    }
  } //end smoother

  '),'
 } // end si loop (includes sub 0)
  
  prevrow = rowi; //update previous row marker only after doing necessary calcs
}//end active rowi
} //end passive rowi
ll+=sum(llrow);
')
    if(!is.null(ctm$w32)) out <- ''
    return(out)}

matcalcs <- function(subjectid,when, matrices, basemats){
  paste0(sapply(matrices, function(x){
      mn=names(matrices[matrices == x])
      whena=paste0('{',paste0(when,collapse=','),'}')
      whenax=when
      whenax[whenax==0] <- 5
      whenax=paste0('{',paste0(whenax,collapse=','),'}') #remove 0 
      out = paste0(
        ifelse(0 %in% when,'',paste0('if(',ifelse(basemats,'si==0 ||',''),'sum(whenmat[',x,',',whenax,']) > 0)')),
      's',mn,'=mcalc(s',mn,',indparams, statetf,',whena,', ',x,', matsetup, matvalues, ',subjectid,', subindices); \n')
    }),collapse='')
}


subjectparaminit<- function(popmats=FALSE,smats=TRUE,matrices=c(mats$base,31, 21,22)){
  if(smats && popmats) stop('smats and popmats cannot both be TRUE!')
  ma <- ctStanMatricesList()$all
  out<-''
  for(mn in matrices){
    m=names(ma)[ma %in% mn]
    out <- paste0(out, '
      matrix[matrixdims[',mn,', 1], matrixdims[',mn,', 2] ] ',
      if(smats) 's',
      if(popmats) 'pop_',
      m,
      if(!smats && !popmats) paste0('[ (savesubjectmatrices && sum(whenmat[',mn,',1:5])) ? nsubjects : 0]'),
      ';')
  }
      
#     
#   paste0(
#     paste0('   matrix[matrixdims[',(mats$base),', 1], matrixdims[',(mats$base),', 2] ] ',
#       ifelse(smats,'s',''),ifelse(popmats,'pop_',''),names(mats$base),if(!smats && !popmats) paste0('[subindices[',(mats$base),']  ? (savesubjectmatrices ? nsubjects : 1) : 1]'),';',collapse=' \n   '),'
# 
#   matrix[nlatent,nlatent] ',ifelse(smats,'s',''),ifelse(popmats,'pop_',''),'asymDIFFUSION',if(!smats && !popmats) '[(subindices[3] + subindices[4])  ? (savesubjectmatrices ? nsubjects : 1) : 1]','; //stationary latent process variance
#   vector[nlatent] ',ifelse(smats,'s',''),ifelse(popmats,'pop_',''),'asymCINT',if(!smats && !popmats) '[(subindices[3] + subindices[7])  ? (savesubjectmatrices ? nsubjects : 1) : 1]','; // latent process asymptotic level
# ','matrix[nlatent, nlatent] ',ifelse(smats,'s',''),ifelse(popmats,'pop_',''), 'DIFFUSIONcov',if(!smats && !popmats) '[subindices[4] ? (savesubjectmatrices ? nsubjects : 1) : 1]',';',
#     collapse='\n')
  return(out)
}

collectsubmats <- function(popmats=FALSE,matrices=c(mats$base,31, 21,22)){ #'DIFFUSIONcov','asymDIFFUSION','asymCINT'
  ma <- ctStanMatricesList()$all
  out<-''
  for(mn in matrices){
    m=names(ma)[ma %in% mn]
    if(!popmats) out <- paste0(out, '
    if(sum(whenmat[',mn,',1:5]) > 0) ',m,'[si] = s',m,';')
    
    if(popmats) out <- paste0(out, 'pop_',m,' = s',m,'; ')
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
    matrix[rows(mat),cols(mat)] o;
  
    for(i in 1:rows(o)){ //set upper tri to lower
      for(j in min(i+1,rows(mat)):rows(mat)){
        o[j,i] =  inv_logit(mat[j,i])*2-1;  // can change cor prior here
        o[i,j] = o[j,i];
      }
      o[i,i]=.999; //avoids correlations of 1
      o[i,] /= sqrt(sum(square(o[i,]))+1e-10);
    }
    return o;
  } 

  matrix sdcovsqrt2cov(matrix mat, int cholbasis){ //covariance from cholesky or unconstrained cor sq root
    if(cholbasis==0)  {
      return(tcrossprod(diag_pre_multiply(diagonal(mat),constraincorsqrt(mat))));
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

 matrix sqkron_sumii(matrix mata){
   int d=rows(mata);
   matrix[d*d,d*d] out;
     for (l in 1:d){
       for (k in 1:d){
         for (j in 1:d){
           for (i in 1:d){
             out[i+(k-1)*d,j+(l-1)*d] = 0;
             if(i==j) out[i+(k-1)*d,j+(l-1)*d] += mata[k,l];
             if(k==l) out[i+(k-1)*d,j+(l-1)*d] += mata[i,j];
           }
         }
       }
     }
   return(out);
 } 

  matrix makesym(matrix mat, int verbose, int pd){
    matrix[rows(mat),cols(mat)] out;
    for(coli in 1:cols(mat)){
    //  if(pd ==1 && mat[coli,coli] < 1e-5){
     //   out[coli,coli] = 1e-5;// 
     // } else 
      out[coli,coli] = mat[coli,coli] + 1e-6; 
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
   vector parvectform(int n, vector rawpar, int when, int[,] ms, data real[,] mval, int subi, int[] subindices, int[] whenvec){
    vector[n] parout;
    if(n>0){
      int outwhen[size(whichequals(whenvec,0,0))] = whenvec[whichequals(whenvec,0,0)]; //outwhen is nonzero elements of whenvec
      for(ri in 1:size(ms)){ //for each row of matrix setup
        if((ms[ri,8]==when || ms[ri,8]==100) && ms[ri,9] < 1 && ms[ri,3] > 0){ //if correct when,not a copyrow, and free parameter
          if(subi ==0 ||  //if population parameter
            ( ms[ri,7] == 8 && subindices[8]) || //or a covariance parameter in an individually varying matrix
            (ms[ri,3] > 0 && (ms[ri,5] > 0 || ms[ri,6] > 0 || ms[ri,8] > 0)) //or there is individual variation
            ){ //otherwise repeated values
            
            parout[whichequals(outwhen, ms[ri,3],1)[1] ] = tform(rawpar[ms[ri,3] ], //which outwhen refers to correct par
              ms[ri,4], mval[ri,2], mval[ri,3], mval[ri,4], mval[ri,6] ); 
           }
        }
      }
    }
  return parout;
  }
  
  
  matrix mcalc(matrix matin, vector tfpars, vector tfstates, int[] when, int m, int[,] ms, data real[,] mval, int subi, int[] subindices){
    matrix[rows(matin),cols(matin)] matout;

    for(ri in 1:size(ms)){ //for each row of matrix setup
      int whenyes = 0;
      for(wi in 1:size(when)) if(when[wi]==ms[ri,8] || ms[ri,8]==100) whenyes = 1; //improve PARS when = 100 thing
      if(m==ms[ri,7] && whenyes){ // if correct matrix and when

        if(subi ==0 ||  //if population parameter
          ( ms[ri,7] == 8 && subindices[8]) || //or a covariance parameter in an individually varying matrix
          (ms[ri,3] > 0 && (ms[ri,5] > 0 || ms[ri,6] > 0 || ms[ri,8] > 0)) //or there is individual variation
          ){ //otherwise repeated values (maybe this check not needed now?

          if(ms[ri,3] > 0 && ms[ri,8]==0)  matout[ms[ri,1], ms[ri,2] ] = tfpars[ms[ri,3]]; //should be already tformed
          if(ms[ri,3] > 0 && ms[ri,8]>0)  matout[ms[ri,1], ms[ri,2] ] = tfstates[ms[ri,3]]; //should be already tformed
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
  
  
  int[] checkoffdiagzero(matrix M){
    int z[rows(M)];
    for(i in 1:rows(M)){
      z[i] = 0;
      for(j in 1:rows(M)){ // check cols and rows simultaneously
        if(i!=j && (M[i,j]!=0.0 || M[j,i] != 0.0)){
          z[i] = 1;
          break;
        }
      }
    }
    return z;
  }
  
   
  matrix expm2(matrix M){
    matrix[rows(M),rows(M)] out;
    int z0[rows(out)] = checkoffdiagzero(M);
    int z1[sum(z0)]; //contains which rowcols need full expm
    int count=1;
    for(j in 1:cols(M)){
      if(z0[j]){
        z1[count]=j;
        count+=1;
      } else {
        for(i in 1:rows(M)){
          if(i!=j){
          out[i,j] =  0; 
          out[j,i] = 0;
          } else out[i,i] = exp(M[i,j]);
        }
      }
    }
    if(size(z1)) out[z1,z1] = matrix_exp(M[z1,z1]);
    return out;
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
  int nopriors;
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
  int statedep[10];
  int choleskymats;
  int intoverstates;
  int verbose; //level of printing during model fit
  int subindices[10];
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
  real Jstep;
  real dokalmanpriormodifier;
  int intoverpopindvaryingindex[intoverpop ? nindvarying : 0];
  int nsJAxfinite;
  int sJAxfinite[nsJAxfinite];
  int nJyfinite;
  int sJyfinite[nJyfinite];
  int taylorheun;
  int difftype;
  int popcovn;
  int llsinglerow;
}
      
transformed data{
  matrix[nlatent+nindvarying,nlatent+nindvarying] IIlatentpop = diag_matrix(rep_vector(1,nlatent+nindvarying));
  vector[nlatentpop-nlatent] nlpzerovec = rep_vector(0,nlatentpop-nlatent);
  vector[nlatent+1] nlplusonezerovec = rep_vector(0,nlatent+1);
  int tieffectindices[nparams]=rep_array(0,nparams);
  int ntieffects = 0;
  
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
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying, nindvarying] rawpopc[4];

',if(!gendata) paste0('
  real ll = 0;
  vector[ndatapoints] llrow = rep_vector(0,ndatapoints);
  matrix[nlatentpop,nlatentpop] etacova[3,savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ycova[3,savescores ? ndatapoints : 0];
  vector[nlatentpop] etaa[3,savescores ? ndatapoints : 0];
  vector[nmanifest] ya[3,savescores ? ndatapoints : 0];
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
      rawpopc[1,j,j] = rawpopsd[j]; //used with intoverpop
      for(i in 1:nindvarying){
        if(i > j){
          counter += 1;
          rawpopc[1,i,j]=sqrtpcov[counter];
          rawpopc[1,j,i]=0;//sqrtpcov[counter];
        }
      }
    }
    rawpopc[3] = tcrossprod( constraincorsqrt(rawpopc[1]));
    rawpopc[4] = makesym(quad_form_diag(rawpopc[3], rawpopsd +1e-8),verbose,1);
    rawpopc[2] = cholesky_decompose(rawpopc[4]); 
  }//end indvarying par setup

  {
',
if(!gendata) ukfilterfunc(ppchecking=FALSE),'
  }
}
      
model{
  if(intoverpop==0 && nindvarying > 0) target+= multi_normal_cholesky_lpdf(baseindparams | rep_vector(0,nindvarying), IIlatentpop[1:nindvarying,1:nindvarying]);

  if(ntipred > 0){ 
    if(nopriors==0) target+= dokalmanpriormodifier * normal_lpdf(tipredeffectparams| 0, tipredeffectscale);
    target+= normal_lpdf(tipredsimputed| 0, tipredsimputedscale); //consider better handling of this when using subset approach
  }

  if(nopriors==0){ //if split files over subjects, just compute priors once
   target+= dokalmanpriormodifier * normal_lpdf(rawpopmeans|0,1);
  
    if(nindvarying > 0){
      if(nindvarying >1) target+= dokalmanpriormodifier * normal_lpdf(sqrtpcov | 0, 1);
      target+= dokalmanpriormodifier * normal_lpdf(rawpopsdbase | ',gsub('normal(','',ctm$rawpopsdbase,fixed=TRUE),';
    }
  } //end pop priors section
  
  ',if(!gendata) 'if(intoverstates==0) target+= normal_lpdf(etaupdbasestates|0,1);','
  
  ',if(!gendata) 'target+= ll; \n','
  if(verbose > 0) print("lp = ", target());
}
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
        x[ri,] = (rawpopc[2] * 
          to_vector(normal_rng(rep_vector(0,nindvarying),rep_vector(1,nindvarying))) + 
          rawpopmeans[indvaryingindex])\';
      }
    }
    
    for(pi in 1:nparams){
      int found=0;
      int pr1;
      int pr2;
      real rawpoppar = rawpopmeans[pi];
      while(!found){
        for(ri in 1:size(matsetup)){
          if(matsetup[ri,9] <=0 && matsetup[ri,3]==pi && matsetup[ri,8]<=0) { //if a free parameter 
            pr1 = ri; 
            pr2=ri;// unless intoverpop, pop matrix row reference is simply current row
            found=1;
            if(intoverpop && matsetup[ri,5]) { //check if shifted
              for(ri2 in 1:size(matsetup)){ //check when state reference param of matsetup corresponds to row of t0means in current matsetup row
                if(matsetup[ri2,8]  && matsetup[ri2,3] == matsetup[ri,1] && 
                matsetup[ri2,3] > nlatent && matsetup[ri2,7] < 20 &&
                matsetup[ri,9] <=0) pr2 = ri2; //if param is dynamic and matches row (state ref) and is not in jacobian
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
')
}
  m <- writemodel()
  return(m)
}
