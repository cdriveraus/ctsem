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
  ln <- ctm$pars$param[ctm$pars$matrix %in% 'PARS' & !is.na(ctm$pars$param)] #get extra pars
 
  for(li in seq_along(ln)){ #for every extra par
    parmatch <- which(ctm$pars$param %in% ln[li] & ctm$pars$matrix %in% 'PARS')
    for(ri in grep(paste0('\\b',ln[li],'\\b'),ctm$pars$param)){ #which rows contain the par
      if(!(ctm$pars$param[ri] == ln[li] & ctm$pars$matrix[ri]=='PARS')){ #that are not the par itself
      # print(ctm$pars$param[ri])
        ctm$pars$param[ri] <- gsub(paste0('\\b',ln[li],'\\b'), #replace with PARS reference...
          paste0('sPARS[',ctm$pars$row[parmatch],',',ctm$pars$col[parmatch],']'),ctm$pars$param[ri])
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
        # if(tftype==0) start.params <- start.params[1:2]
        # 
        # if(!any(res < 1e-12)){
        # fit = try(mize(par = unlist(start.params),
        #   fg = list(fn=ff,gr=ffg),
        #   max_iter=50*ifelse(tryi==1,1,4),abs_tol=1e-3*ifelse(tryi==1,1,1e-4),
        #   rel_tol=1e-5*ifelse(tryi==1,1,1e-4),
        #   method='BFGS'))
        # 
        # # if(tftype > 0) piseq <- 1:4 else piseq <- 1:2
        # # for(pi in piseq){
        # #   if(fit$f > .1) {
        # #     start <- unlist(start.params)
        # #     start[pi] = -start[pi]
        # #     fit = try(mize(par = start,
        # #       fg = list(fn=ff,gr=ffg),
        # #       max_iter=100,abs_tol=1e-5,rel_tol=1e-8,
        # #       method='BFGS'))
        # #   }
        # # }
        # 
        # if(fit$f < .1 && fit$f > 1e-5) {
        #   # message('close, ', round(fit$f,5))
        #   fit = try(mize(par = fit$par, #if close, refine estimate
        #     fg = list(fn=ff,gr=ffg),
        #     max_iter=200,abs_tol=1e-5,rel_tol=1e-9,
        #     method='BFGS'))
        #   # message('close2, ', round(fit$f,5))
        # }
        # }#end if need to optimize
        
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
    # nctspec <- cbind(ctm$pars[newrows,,drop=FALSE],df)
    # nctspec <- 
    nctspec <- merge(ctm$pars,df,by=0,all=TRUE,no.dups = FALSE)
    nctspec <- nctspec[order(as.numeric(nctspec$Row.names)),]
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
    m$pars=ctStanModelCleanctspec(m$pars)
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
      
      
      #new drift
      drift <- m$pars[m$pars$matrix %in% 'DRIFT' & m$pars$row==1 & m$pars$col==1,]
      for(ri in 1:(m$n.latent+nindvaryingsmall)){
        for(ci in 1:(m$n.latent+nindvaryingsmall)){
          newrow <- list('DRIFT',ri,ci,NA,ifelse(m$continuoustime,0,ifelse(ri==ci,1,0)),NA,FALSE,1)
          if(m$n.TIpred > 0) newrow<-c(newrow,rep(FALSE,m$n.TIpred))
          drift[nrow(drift)+1,] <- newrow
        }}
      drift=drift[-1,,drop=FALSE]
      drift=subset(drift,drift$row > m$n.latent | drift$col > m$n.latent)
      # JAx <- drift
      # JAx$matrix <- 'JAx'
      
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
      
      m$pars <- rbind(m$pars, t0m,t0v,drift) #
      m$pars[] <- lapply(m$pars, utils::type.convert, as.is = TRUE)
      
    } #finish loop for non simple t0means indvarying
    extralatents <- seq.int(m$n.latent+1,m$n.latent+length(ivnames),length.out=length(ivnames))
    
    m$intoverpopindvaryingindex <- c(t0mvaryingsimple,extralatents)
    return(m)
  }
}




simplifystanfunction<-function(bcalc){ #input text of list of computations, output simplified form
  if(length(bcalc)==0) return('')
  # bcalcs=paste0(paste0(c(ctm$calcs$driftcint,ctm$calcs$diffusion),';\n',collapse=' '),
  #           jacobianelements(ctm$jacobian$JAx,ntdpred=ctm$n.TDpred,matsetup=matsetup,
  #              textadd=paste0('    sJAx[',rep(1:nlatentpop,nlatentpop),', ',rep(1:nlatentpop,each=nlatentpop),'] = '),
  #              when = 2,remove = c('fixed','simplestate','drift')),collapse=';\n')
  
  bcalcs=strsplit(bcalc,';')[[1]]
  bcalcs=gsub('\n','',bcalcs)
  bcalcs = bcalcs[!nchar(bcalcs)==0]
  bcalcs1=gsub('=.*','',bcalcs)
  bcalcs2=gsub('.*=','',bcalcs)
  bcalcs2 = gsub(",", "___comma___", bcalcs2, fixed = TRUE)
  bcalcs2 = gsub("[", "___leftsquarebracket___", bcalcs2, fixed = TRUE)
  bcalcs2 = gsub("]", "___rightsquarebracket___", bcalcs2, fixed = TRUE)
  
  bcalcs2c=paste0('c(',paste0(bcalcs2,collapse=', '),')')
  
  # scalcs2=gsub(' ','',Simplify(bcalcs2c))
  
  # scalcs2=Deriv(bcalcs2c,nderiv=0) #to simplify
  scalcs2=bcalcs2c #or not to simplify
  
  if(scalcs2 == 'NULL') return('') else{
    
    names(scalcs2)=NULL
    scalcs2=gsub(' ','',scalcs2)
    
    scalcs2 = gsub(",", "; =", scalcs2, fixed = TRUE)
    scalcs2 = gsub("___rightsquarebracket___", "]", scalcs2, fixed = TRUE)
    scalcs2 = gsub("___leftsquarebracket___", "[", scalcs2, fixed = TRUE)
    scalcs2 = gsub("___comma___", ",", scalcs2, fixed = TRUE)
    
    ec=gsub('\\b(c)\\b\\(.*','',scalcs2)
    ec=gsub(';',';  \n',ec)
    scalcs2=gsub('.*\\b(c)\\b\\(','',scalcs2)
    ec=gsub(';','; ',ec)
    scalcs2=gsub(';',';',scalcs2)
    
    ec=gsub('\\{.e','real e',ec)
    ec=gsub('\\{','',ec)
    ec= gsub('\n\\.','\n  real ',ec[1])
    ec=gsub('\\.','',ec)
    
    ec=gsub('<-','= ',ec)
    ec=gsub('\n','\n  ',ec)
    # cat(ec)
    scalcs2=gsub('\\}','',scalcs2)
    scalcs2=gsub('\\)$','',scalcs2)
    # scalcs2=gsub('\\.','',scalcs2)
    
    scalcs2[1] = gsub('^\\s*(\\S{1})',' = \\1',scalcs2[1])
    scalcs2 = strsplit(scalcs2,';')[[1]]
    scalcs = paste0(bcalcs1,scalcs2,';  \n  ')
    # cat(scalcs)
    out = paste0('    {\n  ',paste0(ec,collapse=''),paste0(scalcs,collapse=''),'\n  } \n  ',collapse='')
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
  ctspec$value=as.numeric(ctspec$value)
  ctspec$transform=as.character(ctspec$transform) #ensure consistent type in case custom tforms used
  ctspec$param=gsub(' ','',as.character(ctspec$param)) #remove spces
  comparison=list(NA,NA,FALSE)
  replacement=list(NA,NA,FALSE)
  # names(comparison)=c('param','transform','indvarying')
  for(rowi in 1:nrow(ctspec)){ 
    if( !is.na(ctspec$value[rowi])) {#fixed value replacement
      if(any(c(!is.na(ctspec[rowi,'param']),!is.na(ctspec[rowi,'transform']),ctspec[rowi,'indvarying']))){
        if(ctspec[rowi,'value']!=99999) found<-TRUE
        ctspec[rowi,c('param','transform','indvarying')]=replacement
      }
    }
    if(grepl('[', ctspec$param[rowi],fixed=TRUE) || !is.na(ctspec$value[rowi])){
      if(ctspec$indvarying[rowi]){
        # message('Individual variation requested on deterministic parameter ', ctspec$param[rowi],' , setting to FALSE')
        ctspec$indvarying[rowi] <- FALSE
      }
      if((length(tieffects) > 0 && any(as.logical(ctspec[rowi,tieffects]) %in% TRUE)) || !is.na(ctspec$value[rowi])){
        # message('TI predictor effects requested on deterministic parameter ', ctspec$param[rowi],' , setting to FALSE')
        ctspec[rowi,tieffects] <- FALSE
      }
    } #end deterministic relations check
    if(all(is.na(c(ctspec[rowi,c('value','param')])))) stop('Parameters specified as NA ! Needs a value or character label.')
  } #end row loop
  # if(found) message('Minor inconsistencies in model found - removing param name, transform and indvarying from any parameters with a value specified')
  return(ctspec)
}

ctStanMatricesList <- function(unsafe=FALSE){
  m <- list()
  m$base <- c("T0MEANS","LAMBDA","DRIFT","DIFFUSION","MANIFESTVAR","MANIFESTMEANS", "CINT","T0VAR","TDPREDEFFECT",'PARS')
  m$driftcint <- c('DRIFT','CINT')
  m$diffusion <- c('DIFFUSION','JAx')
  m$tdpred <- c('TDPREDEFFECT','Jtd')
  m$measurement <- c('LAMBDA','MANIFESTMEANS','MANIFESTVAR','Jy')
  m$t0 <- c('T0MEANS','T0VAR','J0')
  # if('PARS' %in% ctm$pars$matrix) {
  #   m$base <- c(m$base, 'PARS')
  #   if(!unsafe && any(sapply(ctm$pars$param[ctm$pars$matrix %in% 'PARS'], function(x) grepl('state[', x,fixed=TRUE) ))){
  #     stop('PARS matrix cannot contain further dependencies, simple parameters only!')
      # m$driftcint <- c(m$driftcint,'PARS')
      # m$measurement <- c(m$measurement,'PARS')
  #   }
  # }
  
  mn=lapply(m, function(mli){
    out=match(x = mli,table = m$base)
    names(out) = mli
    return(out)
  })
  
  mn$jacobian <- c(51,52,53,54)
  mn$t0[3] <- 51
  mn$diffusion[2] <-52
  mn$tdpred[2] <- 53
  mn$measurement[4] <- 54
  names(mn$jacobian) = c('J0','JAx','Jtd','Jy')
  mn$asymptotic <- c(21,22)
  names(mn$asymptotic) = c('asymCINT','asymDIFFUSION')
  mn$all <- c(mn$base,mn$asymptotic,mn$jacobian)
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
  # indvaryingindex <-array(0,dim=c(0))
  indvaryingcounter <- 0
  TIPREDEFFECTsetup <- matrix(0,0,n.TIpred)
  tipredcounter <- 1
  indvar <- 0
  extratformcounter <- 0
  extratforms <- c()
  calcs<-ctm$calcs
  matsetup<-NULL
  mval<-NULL
  # colnames(matsetup) <- c('row','col','param','transform', 'indvarying','tipred', 'matrix','when')
  # colnames(mval) <- c('value','multiplier','meanscale','offset','sdscale','inneroffset')
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
          
          if(!simplestate && grepl('[',ctspec$param[i],fixed=TRUE)){ #if non simple calculation parameter
            # if(grepl('^\\b(state)\\b\\[\\d+\\]$',ctspec$param[i])){ #if a simple state reference
            # if(m %in% names(c(mats$driftcint,mats$diffusion))) when = 2
            # if(m %in% names(mats$tdpred)) when = 3
            # if(m %in% names(mats$measurement)) when = 4
            # if(m %in% names(mats$t0)) when = 1
            #   parameter = gsub('^\\b(state)\\b\\[','',ctspec$param[i]) #remove state[
            #   parameter = gsub(']','',parameter,fixed=TRUE)  #and ], to leave state reference as parameter
            #   indvar <- 0 #state varying anyway
            # } else { #if a non simple calculation
            # 
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
            } else {
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
            if(nm %in% c('JAx',names(c(mats$driftcint,mats$diffusion)))) when = 2
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

  matrixdims <- t(sapply(mats$base, function(m) {
    as.integer(c(max(c(0,matsetup[matsetup[,'matrix'] %in% m,'row'])),
      max(c(0,matsetup[matsetup[,'matrix'] %in% m,'col']))))
  }))
  matsetup[,!colnames(matsetup) %in% 'parname'] <- lapply(matsetup[,!colnames(matsetup) %in% 'parname'],as.integer)
  matvalues <- data.frame(apply(mval,2,as.numeric,.drop=FALSE),stringsAsFactors = FALSE  )

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
  
  matsetup$copyrow[unique(matsetup$copyrow[matsetup$copyrow > 0])] <- 
    -unique(matsetup$copyrow[matsetup$copyrow > 0]) #set rows to be copied from to neg of their own row num
  
  matsetup$copymatrix <-  matsetup$copycol <- NULL #not needed columns after processing
  
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

ctStanModelWriter <- function(ctm, gendata, extratforms,matsetup){
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
  
  
  #consider allowing copies across different when states by dropping mlist from the final loop
  simplestatedependencies <- function(when, mlist) {
    paste0('
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == ',when,' &&(',paste0('
            matsetup[ri,7] ==',mlist,collapse='||'),' )){ //perform calcs appropriate to this section
              real newval;
              newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
              ',paste0('if(matsetup[ri, 7] == ', (mats$all)[(mats$all) %in% mlist],') s', names(mats$all)[(mats$all) %in% mlist],'[matsetup[ ri,1], matsetup[ri,2]] = newval;', collapse = ' \n      '),'
              if(matsetup[ri,9] < 0){
                for(ri2 in 1:size(matsetup)){
                  if(matsetup[ri2,9] == ri){ //if row ri2 is a copy of original row ri
                    ',paste0('if(matsetup[ri2, 7] == ', (mats$all)[(mats$all) %in% mlist],') s', names(mats$all)[(mats$all) %in% mlist],'[matsetup[ri2,1], matsetup[ri2,2]] = newval;', collapse = ' \n      '),'
                  }
                }
              }
            }
          }')
  }
  
  
  # simplifystanfunction(paste0(paste0(c(ctm$calcs$diffusion))))
  
  finiteJ<-function(){
    paste0('
    {
    int zeroint[1];
    vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(sJAxfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0)  state[statei] += Jstep;
      ',simplestatedependencies(when=2,mlist=c(PARS=10,mats$driftcint)),'
      ',simplifystanfunction(paste0(paste0(c(ctm$modelmats$calcs$PARS,ctm$modelmats$calcs$driftcint),';\n',collapse=' '))),' 
      if(statei > 0) {
        sJAx[sJAxfinite,statei] =  sDRIFT[sJAxfinite, ] * state + append_row(sCINT[,1],nlpzerovec)[sJAxfinite]; //compute new change
         if(verbose>1) print("sJAx ",sJAx);
      }
      if(statei== 0 && size(sJAxfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
        base[sJAxfinite] = sDRIFT[sJAxfinite, ] * state + append_row(sCINT[,1],nlpzerovec)[sJAxfinite];
        if(verbose>1) print("base = ",base,"    sjaxinit= ",sJAx);
        for(fi in sJAxfinite){
          sJAx[sJAxfinite,fi] = (sJAx[sJAxfinite,fi] - base[sJAxfinite]) / Jstep; //new - baseline change divided by stepsize
        }
      }
    }
    if(verbose>1) print("sJAx ",sJAx);
    }
    ')
  }
  
  finiteJy<-function(){
    paste0('
    int zeroint[1];
    vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(sJyfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0 && (savescores + intoverstates) > 0)  state[statei] += Jstep;
      ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n',collapse=' ')),
      simplestatedependencies(when=4,mlist=c(PARS=10,mats$measurement)),
      simplifystanfunction(paste0(paste0(c(ctm$modelmats$calcs$measurement),';\n',collapse=' '))),' 
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
    ')
  }
  
  
  
  
  
  ukfilterfunc<-function(ppchecking){
    out<-paste0('
  int si = 0;
  int prevrow=0;
  real prevdt=0;
  real dt;
  real dtsmall;
  int dtchange=1;
  int T0check=0;
  int subjectcount = 0;
  int counter = 1;
  matrix[nlatentpop, nlatentpop] etacov; //covariance of latent states

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] syprior;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypriorcov_sqrt; 
  matrix[nmanifest, nmanifest] ycov; 
  
  matrix[nlatentpop,nlatentpop] Je[savescores ? ndatapoints : 1]; //time evolved jacobian, saved for smoother
  matrix[nlatent*2,nlatent*2] dQi; //covariance from jacobian

  vector[nlatentpop] state = rep_vector(-1,nlatentpop); 
  matrix[nlatentpop,nlatentpop] sJAx; //Jacobian for drift
  matrix[nlatentpop,nlatentpop] sJ0; //Jacobian for t0
  matrix[nlatentpop,nlatentpop] sJtd;//diag_matrix(rep_vector(1),nlatentpop); //Jacobian for nltdpredeffect
  matrix[ nmanifest,nlatentpop] sJy;//Jacobian for measurement 
  
  matrix[nlatentpop,nlatentpop] Kth = rep_matrix(0.0,nlatentpop,nlatentpop); 
  matrix[nlatentpop,nlatentpop] Mth = Kth;

  //linear continuous time calcs
  matrix[nlatent+1,nlatent+1] discreteDRIFT;
  matrix[nlatent,nlatent] discreteDIFFUSION = rep_matrix(0.0,nlatent,nlatent);

  //dynamic system matrices
  ',subjectparaminit(pop=FALSE,smats=TRUE),'
  
  int dokalmanrows[ndatapoints] = dokalmanrowsdata;
  
  if(doonesubject==1){
    dokalmanrows=rep_array(0,ndatapoints);
    for(i in 1:ndatapoints){
      for(subi in 1:size(onesubject)){
        if(fabs(subject[i]-onesubject[subi]) < .5){
          dokalmanrows[i] = 1; 
        }
      }
    }
  }

  for(rowi in 1:(dokalman ? ndatapoints :1)){
  if(dokalmanrows[rowi] ==1) { //used for subset selection
    int o[savescores ? nmanifest : nobs_y[rowi]]; //which obs are not missing in this row
    int o1[savescores ? size(whichequals(manifesttype,1,1)) : nbinary_y[rowi] ];
    int o0[savescores ? size(whichequals(manifesttype,1,0)) : ncont_y[rowi] ];
    
    int od[nobs_y[rowi]] = whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
    int o1d[nbinary_y[rowi] ]= whichbinary_y[rowi,1:nbinary_y[rowi]];
    int o0d[ncont_y[rowi] ]= whichcont_y[rowi,1:ncont_y[rowi]];
    
    if(!savescores){
      o= whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
      o1= whichbinary_y[rowi,1:nbinary_y[rowi]];
      o0= whichcont_y[rowi,1:ncont_y[rowi]];
    }
    if(savescores){ //needed to calculate yprior and yupd ysmooth
      for(mi in 1:nmanifest) o[mi] = mi;
      o1= whichequals(manifesttype,1,1);
      o0= whichequals(manifesttype,1,0);
    }
  
    si = subject[rowi];
    if(prevrow != 0) T0check = (si==subject[prevrow]) ? (T0check+1) : 0; //if same subject, add one, else zero
    if(T0check > 0){
      dt = time[rowi] - time[prevrow];
      dtchange = dt==prevdt ? 0 : 1; 
      prevdt = dt; //update previous dt store after checking for change
    }
    if(savescores || prevrow==0) Je[savescores ? rowi : 1] = IIlatentpop; //elements updated later
    prevrow = rowi; //update previous row marker only after doing necessary calcs
       
    if(T0check == 0) { // calculate initial matrices if this is first row for si
  
    ',subjectparscalc2(popmats=ifelse(gendata,TRUE,TRUE),subjmats=ifelse(gendata,TRUE,TRUE)),'

    etacov =  sT0VAR;
    state = sT0MEANS[,1]; //init and in case of jacobian dependencies

      if(nldynamics==0){ //initialize most parts for nl later!
        if(ntdpred > 0) state[1:nlatent] += sTDPREDEFFECT * tdpreds[rowi];
      }
    } //end T0 matrices
if(verbose > 1) print ("below t0 row ", rowi);

    if(nldynamics==0 && T0check>0){ //linear kf time update
      if(verbose > 1) print ("linear update row ", rowi);
    
      if(continuoustime ==1){
        if(dtchange==1 || (T0check == 1 && (subindices[3] + subindices[7] > 0))){ //if dtchanged or if subject variability
          discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent,1:nlatent],sCINT),nlplusonezerovec\') * dt,drcintoffdiag);
          if(!savescores) Je[1, 1:nlatent, 1:nlatent] = discreteDRIFT[1:nlatent,1:nlatent];
        }
      
        if(dtchange==1 || (T0check == 1 && (subindices[4] + subindices[3] > 0))){ //if dtchanged or if subject variability
          discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
            quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]\' );
        }
      }

      if(continuoustime==0 && T0check == 1){
        if(subjectcount == 1 || subindices[4] + subindices[3] + subindices[7] > 0){ //if first subject or variability
          discreteDRIFT=append_row(append_col(sDRIFT[1:nlatent,1:nlatent],sCINT),nlplusonezerovec\');
          discreteDRIFT[nlatent+1,nlatent+1] = 1;
          if(!savescores) Je[1, 1:nlatent, 1:nlatent] = discreteDRIFT[1:nlatent,1:nlatent];
          discreteDIFFUSION=sDIFFUSIONcov;
        }
      }
      if(savescores) Je[rowi, 1:nlatent, 1:nlatent] = discreteDRIFT[1:nlatent,1:nlatent];
      state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent];
      if(ntdpred > 0) state[1:nlatent] += sTDPREDEFFECT * tdpreds[rowi];
      if(intoverstates==1 || savescores==1) {
        etacov = quad_form(etacov, Je[savescores ? rowi : 1]\');
        if(ndiffusion > 0) etacov[1:nlatent,1:nlatent] += discreteDIFFUSION;
      }
    }//end linear time update


    if(nldynamics==1){ //nldynamics time update
      if(T0check>0){
        vector[nlatentpop] base;
        real intstepi = 0;
        dtsmall = dt / ceil(dt / maxtimestep);
        
        while(intstepi < (dt-1e-10)){
          intstepi = intstepi + dtsmall;
          ',#simplestatedependencies(when=2,mlist=c(mats$driftcint,mats$diffusion,mats$jacobian[2])),' #now done inside finite diff loop
      finiteJ(),
      simplestatedependencies(when=2,mlist=c(mats$diffusion)),
      simplifystanfunction(paste0(paste0(ctm$modelmats$calcs$diffusion,';\n',collapse=' '))),' 
      
      if(statedep[4]) sDIFFUSIONcov[derrind,derrind] = sdcovsqrt2cov(sDIFFUSION[derrind,derrind],choleskymats);
      
        if(continuoustime){
          if(taylorheun==0){
            if(dtchange==1 || statedep[3]||statedep[4]||statedep[7] || 
              (T0check == 1 && (subindices[3] + subindices[4] + subindices[7]) > 0)){
                
              ',if(1==99) 'if(difftype==2 && (statedep[3]||statedep[4])){
                matrix[ndiffusion*2,ndiffusion*2] ebA;
                matrix[ndiffusion*2,ndiffusion*2] bA;
            
                bA[1:ndiffusion,1:ndiffusion] = -sJAx[derrind,derrind];
                bA[1:ndiffusion,(1+ndiffusion):(ndiffusion*2)] = sDIFFUSIONcov[derrind,derrind];
                bA[(1+ndiffusion):(ndiffusion*2),(1+ndiffusion):(ndiffusion*2)] = sJAx[derrind,derrind]\';
                bA[(1+ndiffusion):(ndiffusion*2),1:ndiffusion] = rep_matrix(0,ndiffusion,ndiffusion);
                
                ebA = matrix_exp(bA * dtsmall);
                Je[savescores ? rowi : 1,derrind,derrind] =  ebA[(1+ndiffusion):(ndiffusion*2),(1+ndiffusion):(ndiffusion*2)]\';
                discreteDIFFUSION[derrind,derrind] = Je[savescores ? rowi : 1,derrind,derrind] * ebA[1:ndiffusion,(1+ndiffusion):(ndiffusion*2)];
                discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec\') * dtsmall,drcintoffdiag);
                if(ndiffusion < nlatent && savescores) Je[savescores ? rowi : 1] =  expm2(sJAx * dtsmall, jacoffdiag);
              }
                
              //if(difftype==1){
              //matrix[nlatent,nlatent] V = sDIFFUSIONcov-quad_form(sDIFFUSIONcov, Je[savescores ? rowi : 1]\');
              //discreteDIFFUSION = solvesyl(sJAx[1:nlatent,1:nlatent],-V,discreteDIFFUSION, rep_array(nlatent,1));
              //}
              ','
              
              
              if(difftype==0 || (statedep[3]==0 && statedep[4]==0)){
                discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec\') * dtsmall,drcintoffdiag);
                Je[savescores ? rowi : 1] =  expm2(sJAx * dtsmall, jacoffdiag);
                if(statedep[3]||statedep[4]) sasymDIFFUSION[derrind,derrind] = to_matrix(  -sqkron_sumii(sJAx[derrind,derrind]) \\ 
                  to_vector(sDIFFUSIONcov[derrind,derrind]), ndiffusion,ndiffusion);
                discreteDIFFUSION[derrind,derrind] =  sasymDIFFUSION[derrind,derrind] - quad_form( sasymDIFFUSION[derrind,derrind], Je[savescores ? rowi : 1, derrind,derrind]\' );
              }
            }  else if(savescores) Je[rowi] = Je[rowi-1]; //if not updating
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent]; // ???compute before new diffusion calcs
            if(intoverstates==1 || savescores==1){
              etacov = quad_form(etacov, Je[savescores ? rowi : 1]\');
              etacov[derrind,derrind] += discreteDIFFUSION[derrind,derrind]; 
            }
          }
            
          if(taylorheun==1){
            if(dtchange==1 || statedep[3]||statedep[4] || 
              (T0check == 1 && (subindices[3] + subindices[4]) > 0)){
                Kth = (IIlatentpop - sJAx * (dtsmall /2) );
                Mth = Kth \\ (IIlatentpop + sJAx * (dtsmall /2) );
              }
            state[1:nlatent] += Kth[1:nlatent,1:nlatent] \\
              (sDRIFT[1:nlatent,1:nlatent] * state[1:nlatent] + sCINT[1:nlatent,1]) * dtsmall;
            etacov = quad_form_sym((etacov), Mth\');
            etacov[derrind,derrind] += (Kth[derrind,derrind] \\ sDIFFUSIONcov[derrind,derrind] / Kth[derrind,derrind]\') * dtsmall;
          }
            if(intstepi >= (dt-1e-10) && savescores) Je[rowi] = expm2(sJAx * dt,jacoffdiag); //save approximate exponentiated jacobian for smoothing
          }
  
          if(continuoustime==0){ 
            Je[savescores ? rowi : 1] = sJAx;
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
    } // end nonlinear time update
  
    if(T0check==0){ //nl t0
    state = sT0MEANS[,1]; //in case of t0 dependencies, may have missingness
    ',simplifystanfunction(paste0(paste0(ctm$modelmats$calcs$PARS,';\n',collapse=' '))),
    simplestatedependencies(when=1,mlist=c(mats$t0)),'
    ',paste0(ctm$modelmats$calcs$t0,';',collapse=' '),'
      state = sT0MEANS[,1];
      etacov = quad_form(sT0VAR, sJ0\');
    }//end nonlinear t0
    
    if(ntdpred > 0) {
      int nonzerotdpred = 0;
      for(tdi in 1:ntdpred) if(tdpreds[rowi,tdi] != 0.0) nonzerotdpred = 1;
      if(nonzerotdpred){
      ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n',collapse=' ')),
      simplestatedependencies(when=3,mlist=c(mats$tdpred)),'
      ',paste0(ctm$modelmats$calcs$tdpred,';',collapse=' '),'
        state[1:nlatent] +=   (sTDPREDEFFECT * tdpreds[rowi]); //tdpred effect only influences at observed time point','
        if(statedep[9]) etacov = quad_form(etacov,sJtd\'); //could be optimized
      }
    }//end nonlinear tdpred
  } // end non linear time update

  if(intoverstates==0){
    if(T0check==0) state += cholesky_decompose(sT0VAR) * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)];
    if(T0check>0) state[1:nlatent] +=  cholesky_decompose(discreteDIFFUSION) * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatent)];
  }

if(verbose > 1){
  print("etaprior = ", state);
  print("etapriorcov = ", etacov);
}

if(savescores){
  etapriorcov[rowi] = etacov; 
  etaprior[rowi] = state;
}

 if (nobs_y[rowi] > 0 || savescores) {  // if some observations create right size matrices for missingness and calculate...
    
      ',finiteJy(),'

      if(intoverstates==1 || savescores==1) { //classic kalman
        ycov[o,o] = quad_form(etacov, sJy[o,]\'); // + sMANIFESTVAR[o,o]; shifted measurement error down
        for(wi in 1:nmanifest){ 
          if(Y[rowi,wi] != 99999 || savescores==1) ycov[wi,wi] += square(sMANIFESTVAR[wi,wi]);
          if(manifesttype[wi]==1 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += fabs((syprior[wi] - 1) .* (syprior[wi]));
          if(manifesttype[wi]==2 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += square(fabs((syprior[wi] - round(syprior[wi])))); 
        }
      }
        
      if(intoverstates==0) { //sampled states
        if(ncont_y[rowi] > 0) ypriorcov_sqrt[o0,o0] = sMANIFESTVAR[o0,o0]+1e-8;
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
        K[,od] = mdivide_right(etacov * sJy[od,]\', ycov[od,od]); 
        etacov += -K[,od] * sJy[od,] * etacov;
        state +=  (K[,od] * err[od]);
      }
      
      if(savescores==1) {
        yprior[rowi] = syprior[o];
        etaupd[rowi] = state;
        ypriorcov[rowi] = ycov;
        etaupdcov[rowi] = etacov;
        yupdcov[rowi] = quad_form(etacov, sJy\');
        for(wi in 1:nmanifest) yupdcov[rowi,wi,wi] += square(sMANIFESTVAR[wi,wi]);
        yupd[rowi] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * state[1:nlatent];
      }
      
      
      if(verbose > 1) {
          print("rowi =",rowi, "  si =", si, "  state =",state,"  etacov ",etacov,
            "  syprior =",syprior,"  ycov ",ycov, "  K ",K,
            "  sDRIFT =", sDRIFT, " sDIFFUSION =", sDIFFUSION, " sCINT =", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
            "  sT0VAR", sT0VAR,  " sT0MEANS ", sT0MEANS, "sLAMBDA = ", sLAMBDA, "  sJy = ",sJy,
            " discreteDRIFT = ", discreteDRIFT, "  discreteDIFFUSION ", discreteDIFFUSION, "  sasymDIFFUSION ", sasymDIFFUSION, 
            " DIFFUSIONcov = ", sDIFFUSIONcov,
            "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        }
  
      ',if(!ppchecking){
        'if(nbinary_y[rowi] > 0) llrow[rowi] += sum(log(Y[rowi,o1d] .* (syprior[o1d]) + (1-Y[rowi,o1d]) .* (1-syprior[o1d]))); 
  
        if(size(o0d) > 0 && (llsinglerow==0 || llsinglerow == rowi)){
           if(intoverstates==1) ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(makesym(ycov[o0d,o0d],verbose,1));
           llrow[rowi] +=  multi_normal_cholesky_lpdf(Y[rowi,o0d] | syprior[o0d], ypriorcov_sqrt[o0d,o0d]);
           //errtrans[counter:(counter + ncont_y[rowi]-1)] = 
             //mdivide_left_tri_low(ypriorcov_sqrt[o0d,o0d], err[o0d]); //transform pred errors to standard normal dist and collect
           //ll+= -sum(log(diagonal(ypriorcov_sqrt[o0d,o0d]))); //account for transformation of scale in loglik
           //counter += ncont_y[rowi];
        }
      '},'
    }//end nobs > 0 section
  
  if(savescores && (rowi==ndatapoints || subject[rowi+1] != subject[rowi])){ //at subjects last datapoint, smooth
    int sri = rowi;
    while(sri>0 && subject[sri]==si){
      if(sri==rowi) {
        etasmooth[sri]=etaupd[sri];
        etasmoothcov[sri]=etaupdcov[sri];
      } else{
        matrix[nlatentpop,nlatentpop] smoother;
        smoother = etaupdcov[sri] * Je[sri+1]\' / makesym(etapriorcov[sri+1],verbose,1);
        etasmooth[sri]= etaupd[sri] + smoother * (etasmooth[sri+1] - etaprior[sri+1]);
        etasmoothcov[sri]= etaupdcov[sri] + smoother * ( etasmoothcov[sri+1] - etapriorcov[sri+1]) * smoother\';

      }
      state=etasmooth[sri];
{
      ',finiteJy(),'
}

      ysmooth[sri] = syprior;
      ysmoothcov[sri] = quad_form(etasmoothcov[sri], sJy\'); 
      for(wi in 1:nmanifest){
        ysmoothcov[sri,wi,wi] += square(sMANIFESTVAR[wi,wi]);
        if(manifesttype[wi]==1) ysmoothcov[sri,wi,wi] += fabs((ysmooth[sri,wi] - 1) .* (ysmooth[sri,wi]));
        if(manifesttype[wi]==2) ysmoothcov[sri,wi,wi] += square(fabs((ysmooth[sri,wi] - round(ysmooth[sri,wi])))); 
      }
      sri += -1;
      while(sri > 0 && dokalmanrows[sri]==0) sri+= -1; //skip rows if requested
    }
  } //end smoother
  
  } // end dokalmanrows subset selection
}//end rowi
ll+=sum(llrow);
')
    if(!is.null(ctm$w32)) out <- ''
    return(out)}


subjectparaminit<- function(popmats=FALSE,smats=TRUE){
  if(smats && popmats) stop('smats and popmats cannot both be TRUE!')
  paste0(
    paste0('   matrix[matrixdims[',(mats$base),', 1], matrixdims[',(mats$base),', 2] ] ',
      ifelse(smats,'s',''),ifelse(popmats,'pop_',''),names(mats$base),if(!smats && !popmats) paste0('[subindices[',(mats$base),']  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]'),';',collapse=' \n   '),'

  matrix[nlatent,nlatent] ',ifelse(smats,'s',''),ifelse(popmats,'pop_',''),'asymDIFFUSION',if(!smats && !popmats) '[(subindices[3] + subindices[4])  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]','; //stationary latent process variance
  vector[nlatent] ',ifelse(smats,'s',''),ifelse(popmats,'pop_',''),'asymCINT',if(!smats && !popmats) '[(subindices[3] + subindices[7])  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]','; // latent process asymptotic level
','matrix[nlatent, nlatent] ',ifelse(smats,'s',''),ifelse(popmats,'pop_',''), 'DIFFUSIONcov',if(!smats && !popmats) '[subindices[4] ? (savesubjectmatrices ? nsubjects2 : 1) : 1]',';',
    collapse='\n')
}

subjectparscalc2 <- function(popmats=FALSE,subjmats=TRUE){
  out <- paste0(
    '
 int subjectvec[subjectcount ? 1 : 2];
 vector[nparams] rawindparams = rawpopmeans;
 subjectvec[size(subjectvec)] = si;
 if(subjectcount == 0)  subjectvec[1] = 0; // only needed for subject 0 (pop pars)
 subjectcount = subjectcount + 1;
 for(subjectveci in 1:size(subjectvec)){
  int subi = subjectvec[subjectveci];

  if(subi > 0 && nindvarying > 0 && intoverpop==0) {
    if(fixedsubpars==0) rawindparams[indvaryingindex] += rawpopcovchol * baseindparams[doonesubject ? 1 : subi];
    if(fixedsubpars==1) rawindparams[indvaryingindex] += rawpopcovchol * fixedindparams[doonesubject ? 1 : subi];
  }

  if(subi > 0 &&  ntieffects > 0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipreds[subi]\';

',
    
    paste0('
    for(ri in 1:size(matsetup)){ //for each row of matrix setup
        for(statecalcs in 0:1){ //do state based calcs after initialising t0means
        if(subi ==0 ||  //if population parameter
          ( matsetup[ri,7] == 8 && subindices[8]) || //or a covariance parameter in an individually varying matrix
          (matsetup[ri,3] > 0 && (matsetup[ri,5] > 0 || matsetup[ri,6] > 0 || matsetup[ri,8] > 0)) //or there is individual variation
          ){ //otherwise repeated values
            if( (statecalcs && matsetup[ri,8]>0) || 
              (!statecalcs && matsetup[ri,8]==0) ){ //if doing statecalcs do them, if doing static calcs do them
              real newval;
              if(matsetup[ri,3] > 0)  newval = tform(matsetup[ri,8] ? state[ matsetup[ri,3] ] : rawindparams[ matsetup[ri,3] ], //tform static pars from rawindparams, dynamic from state
                matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
               if(matsetup[ri,3] < 1) newval = matvalues[ri, 1]; //doing this once over all subjects unless covariance matrix -- speed ups possible here, check properly!
              ',paste0('if(matsetup[ri, 7] == ', c(mats$base,mats$jacobian),') s',
                names(c(mats$base,mats$jacobian)),'[matsetup[ ri,1], matsetup[ri,2]] = newval;', collapse = ' \n      '),'
                if(matsetup[ri,9] < 0){ //then send copies elsewhere
                for(ri2 in 1:size(matsetup)){
                  if(matsetup[ri2,9] == ri){ 
                  ',paste0('if(matsetup[ri2, 7] == ', c(mats$base,mats$jacobian),') s', names(c(mats$base,mats$jacobian)),'[matsetup[ri2,1], matsetup[ri2,2]] = newval;', collapse = ' \n      '),'
                  }
                }
              }
            }
          }
        state=sT0MEANS[,1];
      }
    }',collapse='\n    '),'

  // perform any whole matrix transformations, nonlinear calcs based on t0 in order to fill matrices
  ',paste0(ctm$modelmats$calcs$t0,';\n',collapse=' '),'; 
  state=sT0MEANS[,1];
  ',simplifystanfunction(paste0(ctm$modelmats$calcs$PARS,';\n',collapse=' ')),
  paste0(ctm$modelmats$calcs$tdpred,';\n',collapse=' '),';
  ',paste0(ctm$modelmats$calcs$driftcint,';\n',collapse=' '),';
  ',paste0(ctm$modelmats$calcs$diffusion,';\n',collapse=' '),';
  ',paste0(ctm$modelmats$calcs$measurement,';\n',collapse=' '),';
  
  if(subi <= (subindices[4] ? nsubjects2 : 0)) {
    sDIFFUSIONcov = sdcovsqrt2cov(sDIFFUSION,choleskymats);
  }
  if(subi <= ((subindices[3] + subindices[4])  ? nsubjects2 : 0)) {
    if(ndiffusion < nlatent) sasymDIFFUSION = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

    if(continuoustime==1) sasymDIFFUSION[ derrind, derrind] = to_matrix( 
    -sqkron_sumii(sDRIFT[ derrind, derrind ]) \\  to_vector( 
         sDIFFUSIONcov[ derrind, derrind ]), ndiffusion,ndiffusion);

    if(continuoustime==0) sasymDIFFUSION[ derrind, derrind ] = to_matrix( (IIlatent2 - 
      sqkron_prod(sDRIFT[ derrind, derrind ], sDRIFT[ derrind, derrind ])) \\  to_vector(sDIFFUSIONcov[ derrind, derrind ]), ndiffusion, ndiffusion);
  } //end asymdiffusion loops

     if(subi <= (subindices[8] ? nsubjects2 : 0)) {
     if(intoverpop && nindvarying > 0) sT0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopcovsqrt;
      sT0VAR = makesym(sdcovsqrt2cov(sT0VAR,choleskymats),verbose,1); 
      if(nt0varstationary > 0) {
        for(ri in 1:nt0varstationary){ 
          sT0VAR[t0varstationary[ri,1],t0varstationary[ri,2] ] =  sasymDIFFUSION[t0varstationary[ri,1],t0varstationary[ri,2] ];
        }
      }
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
    
    if(subi <= ((subindices[3] + subindices[7])  ? nsubjects2 : 0)){
      if(continuoustime==1) sasymCINT =  -sDRIFT[1:nlatent,1:nlatent] \\ sCINT[ ,1 ];
      if(continuoustime==0) sasymCINT =  (IIlatent - sDRIFT[1:nlatent,1:nlatent]) \\ sCINT[,1 ];
    }
    
    if(nt0meansstationary > 0){
      if(subi <= (subindices[1] ? nsubjects2 : 0)) {
        for(ri in 1:nt0meansstationary){
          sT0MEANS[t0meansstationary[ri,1] , 1] = 
            sasymCINT[t0meansstationary[ri,1] ];
        }
      }
    }
    

    
  ',if(subjmats) collectsubmats(popmats=FALSE),
    if(popmats) paste0('
  if(subi == 0){
',collectsubmats(popmats=TRUE),'
  }
',collapse=''),'
    
} // end subject matrix creation
  ',collapse='\n')
  
  return(out)
}

collectsubmats <- function(popmats=FALSE,matrices=c(names(mats$base),'DIFFUSIONcov','asymDIFFUSION','asymCINT')){
  out<-ifelse(popmats,'', "if(subi == 0 || savesubjectmatrices){ \n")
  for(m in matrices){
    if(!popmats & m %in% names(mats$base)) out <-paste0(out, 'if( (subindices[',grep(m,matrices)[1],'] > 0 && (subi > 0 || savesubjectmatrices==0) ) || 
      (subindices[',which(matrices == m),'] == 0 && subi==0) ) ',m,'[(savesubjectmatrices && subindices[',which(matrices == m),']) ? subi : 1] = s',m,'; \n')
    if(popmats) out <- paste0(out, 'pop_',m,' = s',m,'; ')
  }
  if(!popmats) out <- paste0(out,'
    if( (subindices[4] > 0 && (subi > 0 || savesubjectmatrices==0) ) || 
        (subindices[4] == 0 && subi==0) ) DIFFUSIONcov[(savesubjectmatrices && subindices[4]) ? subi : 1] = sDIFFUSIONcov; 
    if( ((subindices[3] + subindices[4]) > 0 && (subi > 0 || savesubjectmatrices==0) ) || 
        ((subindices[3] + subindices[4]) == 0 && subi==0) ) asymDIFFUSION[(savesubjectmatrices && (subindices[3] + subindices[4]) ) ? subi : 1] = sasymDIFFUSION; 
    if( ((subindices[3] + subindices[7]) > 0 && (subi > 0 || savesubjectmatrices==0) ) || 
        ((subindices[3] + subindices[7]) == 0 && subi==0) ) asymCINT[(savesubjectmatrices && (subindices[3] + subindices[7]) ) ? subi : 1] = sasymCINT; 
        ')
  if(!popmats) out <- paste0(out,' \n }')
  
  return(out)
}



writemodel<-function(){
  paste0('
functions{

 int[] vecequals(int[] a, int test, int comparison){ //do indices of a match test condition?
    int check[size(a)];
    for(i in 1:size(check)){
      if(comparison) check[i] = (test==a[i]) ? 1 : 0;
      if(comparison==0) check[i] = (test==a[i]) ? 0 :1;
    }
    return(check);
  }

int[] whichequals(int[] b, int test, int comparison){  //return array of indices of b matching test condition
    int bsize = size(b);
    int check[bsize] = vecequals(b,test,comparison);
    int whichsize = sum(check);
    int which[whichsize];
    int counter = 1;
    if(bsize > 0){
    for(i in 1:bsize){
      if(check[i] == 1){
        which[counter] = i;
        counter += 1;
      }
    }
    }
    return(which);
  }

  //matrix solvesyl(matrix A, matrix C); //using form F and -V from Wahlstrom Axelsson Gustafsson 2014

   matrix expm2(matrix M,int[] z){
    matrix[rows(M),rows(M)] out;
    int z1[sum(z)];
    int z0[rows(M)-sum(z)];
    int cz1 = 1;
    int cz0 = 1;
    for(i in 1:rows(M)){
      if(z[i] == 1){
        z1[cz1] = i;
        cz1 += 1;
      }
      if(z[i] == 0){
        z0[cz0] = i;
        cz0 += 1;
      }
    }
    if(size(z1) > 0) out[z1,z1] = matrix_exp(M[z1,z1]);
    if(size(z0) > 0){
      out[z0,] = rep_matrix(0,size(z0),rows(M));
      out[,z0] = rep_matrix(0,rows(M),size(z0));
      for(i in 1:size(z0)) out[z0[i],z0[i]] = exp(M[z0[i],z0[i]]);
    }
    return out;
   }

   matrix constraincorsqrt(matrix mat){ //converts from unconstrained lower tri matrix to cor
    matrix[rows(mat),cols(mat)] o;
  
    for(i in 1:rows(o)){ //set upper tri to lower
      for(j in min(i+1,rows(mat)):rows(mat)){
        o[j,i] =  inv_logit(mat[j,i])*2-1;  // can change cor prior here
        o[i,j] = o[j,i];
      }
      o[i,i]=1; // change to adjust prior for correlations
      o[i,] = o[i,] / sqrt(sum(square(o[i,]))+1e-10);
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
    matrix[rows(mata)*rows(matb),cols(mata)*cols(matb)] out;
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
      if(pd ==1 && mat[coli,coli] < 1e-5){
        out[coli,coli] = 1e-5;// 
      } else out[coli,coli] = mat[coli,coli]; 
      for(rowi in coli:rows(mat)){
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
  int nldynamics;
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

  int nt0varstationary;
  int nt0meansstationary;
  int t0varstationary [nt0varstationary, 2];
  int t0meansstationary [nt0meansstationary, 2];

  int nobs_y[ndatapoints];  // number of observed variables per observation
  int whichobs_y[ndatapoints, nmanifest]; // index of which variables are observed per observation
  int ndiffusion; //number of latents involved in system noise calcs
  int derrind[ndiffusion]; //index of which latent variables are involved in system noise calculations
  int drcintoffdiag[nlatent+1];

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
  int matrixdims[',length(mats$base),',2];
  int savescores;
  int savesubjectmatrices;
  int fixedsubpars;
  vector[fixedsubpars ? nindvarying : 0] fixedindparams[fixedsubpars ? nsubjects : 0];
  int dokalman;
  int dokalmanrowsdata[ndatapoints];
  real Jstep;
  real dokalmanpriormodifier;
  int intoverpopindvaryingindex[intoverpop ? nindvarying : 0];
  int nsJAxfinite;
  int sJAxfinite[nsJAxfinite];
  int nJyfinite;
  int sJyfinite[nJyfinite];
  int taylorheun;
  int difftype;
  int jacoffdiag[nlatentpop];
  int njacoffdiagindex;
  int jacoffdiagindex[njacoffdiagindex];
  int popcovn;
  int llsinglerow;
  int doonesubject;
}
      
transformed data{
  matrix[nlatent,nlatent] IIlatent= diag_matrix(rep_vector(1,nlatent));
  matrix[nlatent*nlatent,nlatent*nlatent] IIlatent2 = diag_matrix(rep_vector(1,nlatent*nlatent));
  matrix[nlatentpop,nlatentpop] IIlatentpop = diag_matrix(rep_vector(1,nlatentpop));
  matrix[nindvarying,nindvarying] IIindvar = diag_matrix(rep_vector(1,nindvarying));
  vector[nlatentpop-nlatent] nlpzerovec = rep_vector(0,nlatentpop-nlatent);
  vector[nlatent+1] nlplusonezerovec = rep_vector(0,nlatent+1);
  int nsubjects2 = doonesubject ? 1 : nsubjects;
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
  vector[fixedsubpars ? 0 : (intoverpop ? 0 : nindvarying)] baseindparams[fixedsubpars ? 0 : (intoverpop ? 0 :  (doonesubject ? 1 : nsubjects2) )]; //vector of subject level deviations, on the raw scale
  
  vector[ntipredeffects] tipredeffectparams; // effects of time independent covariates
  vector[nmissingtipreds] tipredsimputed;
  //vector[ (( (ntipredeffects-1) * (1-nopriors) ) > 0) ? 1 : 0] tipredglobalscalepar;
  
  vector[intoverstates ? 0 : nlatentpop*ndatapoints] etaupdbasestates; //sampled latent states posterior
  real onesubject[doonesubject ? doonesubject : 0]; //allows multiple specific
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying, nindvarying] rawpopcovsqrt; 
  matrix[nindvarying, nindvarying] rawpopcovchol; 
  matrix[nindvarying,nindvarying] rawpopcorr;
  matrix[nindvarying,nindvarying] rawpopcov;
',if(!gendata) paste0('
  real ll = 0;
  vector[ndatapoints] llrow = rep_vector(0.0,ndatapoints);
  matrix[nlatentpop,nlatentpop] etapriorcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etaupdcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etasmoothcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ypriorcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] yupdcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ysmoothcov[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaprior[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaupd[savescores ? ndatapoints : 0];
  vector[nlatentpop] etasmooth[savescores ? ndatapoints : 0];
  vector[nmanifest] yprior[savescores ? ndatapoints : 0];
  vector[nmanifest] yupd[savescores ? ndatapoints : 0];
  vector[nmanifest] ysmooth[savescores ? ndatapoints : 0];
  ',subjectparaminit(pop=FALSE,smats=FALSE),'
  ',subjectparaminit(pop=TRUE,smats=FALSE)
  ,collapse=''),'

  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipreds; //tipred values to fill from data and, when needed, imputation vector
  matrix[nparams, ntipred] TIPREDEFFECT; //design matrix of individual time independent predictor effects
  //real tipredglobalscale = 1.0;
  
  //if( ((ntipredeffects-1) * (1-nopriors))  > 0)  tipredglobalscale = exp(tipredglobalscalepar[1]);

  if(ntipred > 0){ 
    int counter = 0;
    for(coli in 1:cols(tipreds)){ //insert missing ti predictors
      for(rowi in 1:rows(tipreds)){
        if(tipredsdata[rowi,coli]==99999) {
          counter += 1;
          tipreds[rowi,coli] = tipredsimputed[counter];
        } else tipreds[rowi,coli] = tipredsdata[rowi,coli];
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
      rawpopcovsqrt[j,j] = rawpopsd[j]; //used with intoverpop
      for(i in 1:nindvarying){
        if(i > j){
          counter += 1;
          rawpopcovsqrt[i,j]=sqrtpcov[counter];
          rawpopcovsqrt[j,i]=0;//sqrtpcov[counter];
        }
      }
    }
    rawpopcorr = tcrossprod( constraincorsqrt(rawpopcovsqrt));
    rawpopcov = makesym(quad_form_diag(rawpopcorr, rawpopsd +1e-8),verbose,1);
    rawpopcovchol = cholesky_decompose(rawpopcov); 
  }//end indvarying par setup

  {
',
if(!gendata) ukfilterfunc(ppchecking=FALSE),'
  }
}
      
model{
  if(doonesubject==0 ||onesubject[1] > .5){ 
    if(intoverpop==0 && fixedsubpars == 1 && nindvarying > 0) target+= multi_normal_cholesky_lpdf(fixedindparams | rep_vector(0,nindvarying),IIindvar);
    if(intoverpop==0 && fixedsubpars == 0 && nindvarying > 0) target+= multi_normal_cholesky_lpdf(baseindparams | rep_vector(0,nindvarying), IIindvar);
  }
  if(doonesubject==0 ||onesubject[1] < .5){ 
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
  }
  
  if(intoverstates==0) target+= normal_lpdf(etaupdbasestates|0,1);
  
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
  vector[ndatapoints] llrow = rep_vector(0.0,ndatapoints);
  matrix[nlatentpop,nlatentpop] etapriorcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etaupdcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etasmoothcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ypriorcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] yupdcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ysmoothcov[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaprior[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaupd[savescores ? ndatapoints : 0];
  vector[nlatentpop] etasmooth[savescores ? ndatapoints : 0];
  vector[nmanifest] yprior[savescores ? ndatapoints : 0];
  vector[nmanifest] yupd[savescores ? ndatapoints : 0];
  vector[nmanifest] ysmooth[savescores ? ndatapoints : 0];
  vector[nmanifest] Ygen[ndatapoints];
  ',subjectparaminit(pop=FALSE,smats=FALSE),'
  ',subjectparaminit(pop=TRUE,smats=FALSE)
  ,collapse=''),'

  {
    matrix[popcovn, nindvarying] x;
    if(nindvarying){
      for(ri in 1:rows(x)){
        x[ri,] = (rawpopcovchol * 
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
          if(matsetup[ri,9] <=0 && matsetup[ri,3]==pi && matsetup[ri,8]==0) { //if a free parameter 
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
