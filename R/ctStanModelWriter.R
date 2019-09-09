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
    # Get some param values over multiplier wide range, and compute the corresponding y
    # values.
    param = c(seq(-2, 2, .1),seq(-10,10,.5),c(rnorm(10)))
    y = eval(parse(text=e)) #eval(substitute(substitute(e, list(param = param)), list(e = as.quoted(e)[[1]]))))
    param <- param[abs(y) < 1e5]; 
    y <- y[abs(y) < 1e5]
    
    # Try to fit each formula to the data.
    for(i in 1:nrow(formula.types)) {
      start.params = list(multiplier = 1.0, offset = 0)
      if(!is.na(formula.types$meanscale[i])) {
        start.params[["meanscale"]] = 1.0
      }
      if(!is.na(formula.types$inneroffset[i])) {
        start.params[["inneroffset"]] = 0
      }
      
      ff <- function(pars){
        multiplier=pars[1];
        offset=pars[2];
        if(i > 1) {
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
        # g=try(numDeriv::grad(ff,pars,method='simple',
        #   method.args=list(eps=1e-8,inneroffset=1e-10,r=2)
        # ),silent=TRUE)
        g=try(sapply(1:length(pars),function(x) {
          parx=pars;parx[x]<-pars[x]+1e-8;
          (b-ff(parx))/1e-8
        }))
        if(class(g)=='try-error') g <- rnorm(pars)
        if(any(is.na(g))) g[is.na(g)] <- rnorm(sum(is.na(g)))
        return(-g)
      }
      # 
      fit = try(mize(par = unlist(start.params),
        fg = list(fn=ff,gr=ffg),
        max_iter=100,abs_tol=1e-3,rel_tol=1e-5,
        method='BFGS'))
      if(fit$f > .1) {
        start.params[1] = -1
        fit = try(mize(par = unlist(start.params),
          fg = list(fn=ff,gr=ffg),
          max_iter=100,abs_tol=1e-3,rel_tol=1e-5,
          method='BFGS'))
      }
      
      if(fit$f < .1 && fit$f > 1e-5) {
        # message('close, ', round(fit$f,3))
        fit = try(mize(par = unlist(start.params), #if close, refine estimate
          fg = list(fn=ff,gr=ffg),
          max_iter=500,abs_tol=1e-5,rel_tol=1e-6,
          method='BFGS'))
      }

      formula.types$offset[i] = fit$par[2] #round(coef(fit)[["offset"]])
      formula.types$multiplier[i] = fit$par[1] #round(coef(fit)[["multiplier"]])
      if(!is.na(formula.types$meanscale[i])) {
        formula.types$meanscale[i] = fit$par[3] #round(coef(fit)[["meanscale"]])
      }
      if(!is.na(formula.types$inneroffset[i])) {
        formula.types$inneroffset[i] = fit$par[4] #round(coef(fit)[["inneroffset"]])
      }
      formula.types$lsfit[i] = fit$f #AIC(fit)
    }
    # Return the values we found.
    # print(formula.types)
    # return(formula.types %>%
    #     filter(lsfit == min(lsfit)) %>%
    #     mutate(inneroffset = coalesce(inneroffset, 0),
    #       meanscale = coalesce(meanscale, 1)) %>%
    #     select(type, offset, inneroffset, multiplier, meanscale,lsfit))
    return(formula.types[which(formula.types$lsfit %in% min(formula.types$lsfit)),])
  }
  
  newrows <- which(!is.na(ctm$pars$transform))
  eqs <- ctm$pars$transform[newrows]
  uniqueeqs <- unique(eqs)
  
    
  l=lapply(uniqueeqs,fit.eqs)
  eqmatch <- match(eqs, uniqueeqs)
  df=data.frame(do.call(rbind,l))
  df[1:length(eqmatch),] <- df[eqmatch,]
  df[,] <- lapply(df,function(x) if(is.numeric(x)) return(round(x,2)) else return(x))
  df$formula <- eqs
  df[df$lsfit > .1,c('offset','inneroffset','multiplier','meanscale')]<-NA
  colnames(df)[1] <- 'transform'
  df$transform[df$lsfit > .1] <- eqs[df$lsfit > .1]
  df$lsfit <- NULL
  rownames(df) <- newrows
  
  ctm$pars$transform <- NULL
   # nctspec <- cbind(ctm$pars[newrows,,drop=FALSE],df)
   # nctspec <- 
     nctspec <- merge(ctm$pars,df,by=0,all=TRUE,no.dups = FALSE)
     nctspec <- nctspec[order(as.numeric(nctspec$Row.names)),]
     nctspec <- nctspec[,c('matrix','row','col','param','value','transform','multiplier',
     'offset','meanscale','inneroffset','indvarying','sdscale',
     colnames(ctm$pars)[grep('_effect',colnames(ctm$pars),fixed=TRUE)]) ]
   ctm$pars <- nctspec
  return(ctm)
  
}


ctStanModelIntOverPop <- function(m){ 
  if(sum(m$pars$indvarying) < 1) {
    message('No individual variation for ctStanModelIntOverPop to work with!')
    return(m)
  } else {
    m$pars=ctStanModelCleanctspec(m$pars)
    
    t0mvaryingsimple <- m$pars$row[m$pars$indvarying & m$pars$matrix %in% 'T0MEANS' & m$pars$transform==0] #which t0means are indvarying and not transformed
    t0mvaryingnames <- m$pars$param[m$pars$indvarying & m$pars$matrix %in% 'T0MEANS'& m$pars$transform==0] #which t0means are indvarying
    t0mnotvarying <- m$pars$row[!m$pars$indvarying & m$pars$matrix %in% 'T0MEANS']
    # m$pars$indvarying[m$pars$matrix %in% 'T0MEANS'] <- FALSE 
    
    ivnames <- unique(m$pars$param[m$pars$indvarying & !m$pars$param %in% t0mvaryingnames]) #don't need new states for t0means
    m$latentPopNames=ivnames
    ivnamesfull <- c(t0mvaryingnames,ivnames) #for t0var naming
    nindvaryingsmall <- length(ivnames)
    
    #new t0means
    t0m <- m$pars[m$pars$param %in% ivnames,]
    t0m$matrix <- 'T0MEANS'
    t0m$col <- 1
    t0m$row <- (m$n.latent+1):(m$n.latent+nindvaryingsmall)
    t0m$transform <- 0
    t0m$multiplier <- 1
    t0m$meanscale <- 1
    t0m$offset <- 0
    t0m$inneroffset <- 0
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
            NA, #ifelse(ri > m$n.latent || ci > m$n.latent, 1,m$pars$sdscale[(m$pars$matrix %in% 'T0VAR' & m$pars$row==ri & m$pars$col==ci)]), #multiplier (sdscale)
            NA,#2,
            NA,#ifelse(ci == ri,0, 0), #offset
            NA,#0
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
        newrow <- list('DRIFT',ri,ci,NA,ifelse(m$continuoustime,0,1),0,1,1,0,0,FALSE,1)
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
      m$pars$param[m$pars$param %in% ivi] <- paste0( 'state[',m$n.latent+match(ivi,ivnames),']')
    }
    
    # m$pars$indvarying  <-FALSE
    
    m$pars <- rbind(m$pars, t0m,t0v,drift) #
    m$pars[] <- lapply(m$pars, utils::type.convert, as.is = TRUE)
    
    m$intoverpopindvaryingindex <- c(t0mvaryingsimple,(m$n.latent+1):(m$n.latent+length(ivnames)))
    return(m)
  }
}




simplifystanfunction<-function(bcalc){ #input text of list of computations, output simplified form
  
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
    scalcs2=gsub('\\.','',scalcs2)
    
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
    if(grepl('[', ctspec$param[rowi],fixed=TRUE)){
      
      if(ctspec$indvarying[rowi]){
        message('Individual variation requested on deterministic parameter ', ctspec$param[rowi],' , setting to FALSE')
        ctspec$indvarying[rowi] <- FALSE
      }
      if(length(tieffects) > 0 && any(as.logical(ctspec[rowi,tieffects]))){
        message('TI predictor effects requested on deterministic parameter ', ctspec$param[rowi],' , setting to FALSE')
        ctspec[rowi,tieffects] <- FALSE
      }
    } #end deterministic relations check
    if(all(is.na(c(ctspec[rowi,c('value','param')])))) stop('Parameters specified as NA ! Needs a value or character label.')
  } #end row loop
  # if(found) message('Minor inconsistencies in model found - removing param name, transform and indvarying from any parameters with a value specified')
  return(ctspec)
}

ctStanMatricesList <- function(ctm,unsafe=FALSE){
  m <- list()
  m$base <- c("T0MEANS","LAMBDA","DRIFT","DIFFUSION","MANIFESTVAR","MANIFESTMEANS", "CINT","T0VAR","TDPREDEFFECT",'PARS')
  m$driftcint <- c('DRIFT','CINT')
  m$diffusion <- 'DIFFUSION'
  m$tdpred <- 'TDPREDEFFECT'
  m$measurement <- c('LAMBDA','MANIFESTMEANS','MANIFESTVAR')
  m$t0 <- c('T0MEANS','T0VAR')
  # if('PARS' %in% ctm$pars$matrix) {
  #   m$base <- c(m$base, 'PARS')
  #   if(!unsafe && any(sapply(ctm$pars$param[ctm$pars$matrix %in% 'PARS'], function(x) grepl('state[', x,fixed=TRUE) ))){
  #     stop('PARS matrix cannot contain further dependencies, simple parameters only!')
  #     # m$driftcint <- c(m$driftcint,'PARS')
  #     # m$measurement <- c(m$measurement,'PARS')
  #   }
  # }
  
  mn=lapply(m, function(mli){
    out=match(x = mli,table = m$base)
    names(out) = mli
    return(out)
  })
  
  mn$jacobian <- c(51,52,53,54)
  names(mn$jacobian) = c('J0','JAx','Jtd','Jy')
  mn$asymptotic <- c(21,22)
  names(mn$asymptotic) = c('asymCINT','asymDIFFUSION')
  mn$all <- c(mn$base,mn$asymptotic,mn$jacobian)
  return(mn)
}


ctStanModelMatrices <-function(ctm){
  mats <- ctStanMatricesList(ctm,unsafe=TRUE)
  ctspec <- ctm$pars
  n.TIpred <- ctm$n.TIpred
  matsetup <-list()
  matvalues <-list()
  freepar <- 0
  freeparcounter <- 0
  indvaryingindex <-array(0,dim=c(0))
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
  for(m in names(mats$base)){
    for(i in 1:nrow(ctspec)){ 
      if(ctspec$matrix[i] == m) {
        when = 0 #default, static parameter
        parameter = 0 #for fixed values
        if(!is.na(ctspec$param[i])){ #if non fixed parameter,
          
          if(grepl('[',ctspec$param[i],fixed=TRUE)){ #if calculation parameter
            if(grepl('^\\b(state)\\b\\[\\d+\\]$',ctspec$param[i])){ #if a simple state reference
              if(m %in% names(c(mats$driftcint,mats$diffusion))) when = 2
              if(m %in% names(mats$tdpred)) when = 3
              if(m %in% names(mats$measurement)) when = 4
              if(m %in% names(mats$t0)) when = 1
              parameter = gsub('^\\b(state)\\b\\[','',ctspec$param[i]) #remove state[
              parameter = gsub(']','',parameter,fixed=TRUE)  #and ], to leave state reference as parameter
              indvar <- 0 #state varying anyway
            } else { #if a non simple calculation
              calcs <- c(calcs, paste0(ctspec$matrix[i],'[',ctspec$row[i], ', ', ctspec$col[i],'] = ',
                ctspec$param[i]))
              when= -999 #never use this row of ctspec during calculations (calc is extracted and placed elsewhere)
              indvar=0
            }
          } #end calculation parameters
          
          if(!grepl('[',ctspec$param[i],fixed=TRUE)){ #if non calculation parameter
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
                extratforms <- paste0(extratforms,'if(transform==',-10-extratformcounter,') out = ',
                  ctspec$offset[i],' + ',ctspec$multiplier[i],' * (inneroffset + ',
                  gsub('param', paste0('param * ',ctspec$meanscale[i]),ctspec$transform[i]),');')
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
        
        mdatnew <- data.frame(
          parname=ctspec$param[i],
          row=ctspec$row[i],
          col=ctspec$col[i],
          param=parameter, # ifelse(!is.na(ctspec$param[i]) && !grepl('[',ctspec$param[i],fixed=TRUE),freepar, 0), #freepar reference
          transform=ifelse(is.na(as.integer(ctspec$transform[i])), -1, ctspec$transform[i]), #transform
          indvarying=ifelse(!is.na(ctspec$param[i]),indvar,0), #indvarying
          tipred=ifelse(any(TIPREDEFFECTsetup[freepar,] > 0), 1, 0), #tipredvarying
          matrix=which(names(mats$base)==m), #matrix reference
          when=when,#when to compute
          stringsAsFactors = FALSE)
        # rownames(mdatnew) <- 
        if(is.null(matsetup)) matsetup <- mdatnew else matsetup<-rbind(matsetup,mdatnew)
        
        mvalnew<-ctspec[i,c('value','multiplier','meanscale','offset','sdscale','inneroffset'),drop=FALSE]
        if(is.null(mval)) mval <- mvalnew else mval<-rbind(mval,mvalnew)
      }
    }
    if(!is.null(mval)) mval[is.na(mval)] <- 99999 else mval<-array(0,dim=c(0,6))
    
  }
  
  # matrixdims <- t(sapply(matsetup, function(m) as.integer(c(max(c(0,m[,1])),max(c(0,m[,2]))))))
  matrixdims <- t(sapply(mats$base, function(m) {
    as.integer(c(max(c(0,matsetup[matsetup[,'matrix'] %in% m,'row'])),
      max(c(0,matsetup[matsetup[,'matrix'] %in% m,'col']))))
  }))
  matsetup[,!colnames(matsetup) %in% 'parname'] <- lapply(matsetup[,!colnames(matsetup) %in% 'parname'],as.integer)
  matvalues <- data.frame(apply(mval,2,as.numeric,.drop=FALSE),stringsAsFactors = FALSE  )
  # matsetup <- data.frame(apply(do.call(rbind,matsetup),c(1,2),as.integer,.drop=FALSE),stringsAsFactors = FALSE )
  # matvalues <- data.frame(apply(do.call(rbind,matvalues),c(1,2),as.numeric,.drop=FALSE),stringsAsFactors = FALSE  )
  
  
  
  #add jacobian fixed values
  
  if(!is.null(ctm$jacobian)){
    
    Jvalues=data.frame(matrix(0,sum(!is.na(suppressWarnings(as.numeric(unlist(ctm$jacobian))))),
      ncol(matvalues),dimnames=list(NULL,colnames(matvalues))))
    Jsetup=data.frame(matrix(0,sum(!is.na(suppressWarnings(as.numeric(unlist(ctm$jacobian))))),
      ncol(matsetup),dimnames=list(NULL,colnames(matsetup))))
    counter=1
    for(jmati in 1:length(ctm$jacobian)){
      if(!is.null( dim(ctm$jacobian[[jmati]]))){
        for(ri in 1:nrow(ctm$jacobian[[jmati]])){
          for(ci in 1:ncol(ctm$jacobian[[jmati]])){
            if(!is.na(suppressWarnings(as.numeric(ctm$jacobian[[jmati]][ri,ci])))){ #if jacobian has a fixed value
              Jvalues[counter,'value'] <- ctm$jacobian[[jmati]][ri,ci]
              Jsetup[counter,c('row','col','param','transform','matrix')] <- list(ri,ci,0,-1,50+jmati)
              counter <- counter + 1
            }
            # if(grepl('Jtform___',ctm$jacobian[[jmati]][ri,ci])){ #if jacobian has a simple state reference form
            #   # Jvalues[counter,'value'] <- ctm$jacobian[[jmati]][ri,ci]
            #   
            #   
            #   matsetrow = which(rownames(matsetup) %in% gsub(']','.',gsub('[','.',gsub('Jtform___','',ctm$jacobian[[jmati]][ri,ci]),fixed=TRUE),fixed=TRUE))
            #   if(length(matsetrow) > 1) warning('matsetrow duplicates!')
            #   Jsetup[counter,c('param','transform')] <- matsetup[matsetrow,c('param','transform')]
            #   Jsetup[counter,c('row','col','when')] <- c(ri,ci,-1)
            #   Jvalues[counter, c('multiplier','meanscale','offset','inneroffset')] <- matvalues[matsetrow, c('multiplier','meanscale','offset','inneroffset')]
            #   counter <- counter + 1
            # }
            
          }
        }
      }
    }
    
    # for(ri in which(matsetup$when > 0 & #for nonlinear rows with a simple state reference (but not measurement related)
    #     # matsetup$transform > 0 & #even 0 transforms relevant due to multiplier
    #     matsetup$when != 4 & #no measurement jacobian yet
    #     matsetup$matrix != 4)){ #no diffusion matrix jacobian yet
    #   ms=matsetup[ri,] #ensure that the jacobian transform is also applied
    #   mv=matvalues[ri,]
    #   message('problem section, recompilation maybe needed!')
    #   newJrow <- list(row=ms$row,col=ms$param,param=ms$param,transform=ms$transform+50,
    #     indvarying=0,tipred=0,matrix=50+ms$when,when=ms$when)
    #   Jsetup[nrow(Jsetup)+1,] <- newJrow
    #   Jvalues[nrow(Jvalues)+1,] <- mv
    # }
    
    
    # Jsetup$when <- rep(-1,nrow(Jsetup))
    if(nrow(Jsetup) > 0) Jsetup$parname <- paste0('J',Jsetup$matrix,'__',Jsetup$row,'_',Jsetup$col)
    
    matsetup <- rbind(matsetup,Jsetup)
    matvalues <- rbind(matvalues,Jvalues)
    
  }#end jacobian additions
  
  # popvalues <- data.frame(matvalues[matsetup[,'param'] !=0 & matsetup[,'when']>=0,,drop=FALSE])
  # popsetup <- matsetup[matsetup[,'param'] !=0 & matsetup[,'when']>=0,,drop=FALSE]
  # popsetup <- data.frame(popsetup)
  # rownames(popsetup) <- NULL
  # rownames(popvalues) <- NULL
  # popsetup <- data.frame(parname,lapply(popsetup,as.integer),stringsAsFactors = FALSE)
  # popvalues <- data.frame(parname,lapply(popvalues,as.numeric),stringsAsFactors = FALSE)
  
  return(list(
    matsetup=matsetup,   matvalues=matvalues, 
    extratforms=extratforms,
    TIPREDEFFECTsetup=TIPREDEFFECTsetup,
    matrixdims=matrixdims,
    calcs=calcs))
}



ctStanCalcsList <- function(ctm){
  #extract any calcs from model and ensure spec is correct
  temp <- ctm$calcs# temp <- sapply(ctm$calcs,function(x) gsub(' ','',Simplify(x),fixed=TRUE))
  names(temp) <- NULL
  # calcs<-list()
  mats<-ctStanMatricesList(ctm)
  # calcindices <- grep('\\]|\\[',ctm$pars$param)
  # if(length(calcindices) > 0){
  #   for(ci in calcindices){
  #     temp <- c(temp,paste0(ctm$pars$matrix[ci],
  #       '[',ctm$pars$row[ci], ', ', ctm$pars$col[ci],'] = ',
  #       ctm$pars$param[ci]),';')
  #   }
  #   ctm$pars$value[calcindices] <- 99999
  #   ctm$pars$param[calcindices] <- NA
  #   if(ctm$n.TIpred > 0) ctm$pars[calcindices,c('indvarying',paste0(ctm$TIpredNames,'_effect'))] <- FALSE
  # }
  
  for(mati in unique(ctm$pars$matrix)){
    temp <- gsub(mati,paste0('s',mati),temp)
  }
  calcs <- lapply(mats[!names(mats) %in% c('base','all')], function(mlist) { #add custom calculations to correct list of matrices
    out <- temp[unlist(sapply(temp, function(y) any(sapply(names(mlist), function(mli) grepl(mli,y)))))]
    return(out)
  })
  
  ctm$calcs <- calcs
  
  return(ctm)
}


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
  
  mats <- ctStanMatricesList(ctm)
  nlatentpop <- max(ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS'])
  nmanifest <- ctm$n.manifest
  mx <- listOfMatrices(ctm$pars)
  # #check when / if PARS needs to be computed
  # for(mlist in names(mats[-1])){
  #   if(any(unlist(lapply(ctm$calcs[[mlist]], function(m) grepl('sPARS',m))))) mats[[mlist]]=c(mats[[mlist]],'PARS')
  # }
  
  
  #adjust diffusion calcs to diffusionsqrt
  # ctm$calcs$diffusion <- gsub('sDIFFUSION','sDIFFUSIONsqrt',ctm$calcs$diffusion)
  # intoverpopdynamiccalcs <- gsub('sDIFFUSION','sDIFFUSIONsqrt',intoverpopdynamiccalcs)
  
  #save calcs without intoverpop for param intitialisation
  
  
  
  simplestatedependencies <- function(when, mlist) {
    paste0('
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == ',when,'){ //perform calcs appropriate to this section
            real newval;
            newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
            ',paste0('if(matsetup[ri, 7] == ', (mats$all)[(mats$all) %in% mlist],') s', names(mats$all)[(mats$all) %in% mlist],'[matsetup[ ri,1], matsetup[ri,2]] = newval;', collapse = ' \n      '),'
            }
          }')
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
      ',simplestatedependencies(when=2,mlist=c(mats$driftcint,mats$diffusion,mats$jacobian[2])),'
      ',simplifystanfunction(paste0(paste0(c(ctm$calcs$driftcint,ctm$calcs$diffusion),';\n',collapse=' '))),' 
      if(statei > 0) {
        sJAx[sJAxfinite,statei] =  (sDRIFT[sJAxfinite, ] * state + append_row(sCINT[,1],rep_vector(0,nlatentpop-nlatent)))[sJAxfinite]; //compute new change
         if(verbose>1) print("sJAx ",sJAx);
      }
      if(statei== 0 && size(sJAxfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
        base[sJAxfinite] = (sDRIFT[sJAxfinite, ] * state + append_row(sCINT[,1],rep_vector(0,nlatentpop-nlatent)))[sJAxfinite];
        if(verbose>1) print("base = ",base,"    sjaxinit= ",sJAx);
        for(fi in sJAxfinite){
        //print("fi!!!!! ",fi);
          sJAx[sJAxfinite,fi] = (sJAx[sJAxfinite,fi] - base[sJAxfinite]) / Jstep; //new - baseline change divided by stepsize
        }
      }
    }
    if(verbose>1) print("sJAx ",sJAx);
    }
    ')
  }
  
  
  
  
  ukfilterfunc<-function(ppchecking){
    out<-paste0('
  int si = 0;
  int subjectcount = 0;
  int counter = 0;
  matrix[nlatentpop, nlatentpop] etacov; //covariance of latent states
  real timei = 0;
  real dt = 0;
  int dtchange;
  real integrationsteps;
  real dtsmall;
  real prevdt = 0;
  int T0check;

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] yprior;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypriorcov_sqrt; 
  matrix[nmanifest, nmanifest] ycov; 
  
  matrix[nlatentpop,nlatentpop] Je[ndatapoints]; //time evolved jacobian, saved for smoother
  matrix[nlatent*2,nlatent*2] dQi; //covariance from jacobian

  vector[nmanifest+nmanifest+ (savescores ? nmanifest*2+nlatent*2 : 0)] kout[ndatapoints];

  vector[nlatentpop] state = rep_vector(-1,nlatentpop); 
  matrix[nlatentpop,nlatentpop] sJAx; //Jacobian for drift
  matrix[nlatentpop,nlatentpop] sJ0; //Jacobian for t0
  matrix[nlatentpop,nlatentpop] sJtd;//diag_matrix(rep_vector(1),nlatentpop); //Jacobian for nltdpredeffect
  matrix[ nmanifest,nlatentpop] Jy[ndatapoints];//store Jacobian for measurement over time
  matrix[ nmanifest,nlatentpop] sJy;//Jacobian for measurement 
  matrix[nmanifest,nlatent] tLAMBDA[ndatapoints]; // store lambda time varying for smoother

  //linear continuous time calcs
  matrix[nlatent+1,nlatent+1] discreteDRIFT;
  matrix[nlatent,nlatent] discreteDIFFUSION;

  //dynamic system matrices
  ',subjectparaminit(pop=FALSE,smats=TRUE),'

  if(nldynamics==0) discreteDIFFUSION = rep_matrix(0,nlatent,nlatent); //in case some elements remain zero due to derrind

  if(savescores) kout = rep_array(rep_vector(99999,rows(kout[1])),ndatapoints);
  
  for(rowi in 1:(dokalman ? ndatapoints :1)){
  if(dokalmanrows[rowi] ==1) { //used for subset selection
    matrix[nldynamics ? nlatentpop : 0, ukffull ? 2*nlatentpop +2 : nlatentpop + 2 ] ukfstates; //sampled states relevant for dynamics
    matrix[nldynamics ? nmanifest : 0 , ukffull ? 2*nlatentpop +2 : nlatentpop + 2] ukfmeasures; // expected measures based on sampled states

    T0check = ( (si == subject[rowi]) ? (T0check + 1) : 0 ) ; //if same subject, add 1 to t0check, else set to 0
    if(T0check > 0) prevdt = dt;
    dt = T0check ? time[rowi] - timei : 0;
    timei = time[rowi]; //must come after dt!
    si=subject[rowi]; //only update subject after t0check!

    
    if(T0check == 0) { // calculate initial matrices if this is first row for si
  
    ',subjectparscalc2(popmats=ifelse(gendata,TRUE,TRUE),subjmats=ifelse(gendata,TRUE,TRUE)),'

    etacov =  sT0VAR;
    state = sT0MEANS[,1]; //init and in case of jacobian dependencies

      if(nldynamics==0){ //initialize most parts for nl later!
        if(ntdpred > 0) state[1:nlatent] += sTDPREDEFFECT * tdpreds[rowi];
      }
    } //end T0 matrices
if(verbose > 1) print ("below t0 row ", rowi);
   
    if(T0check >0)  dtchange = ( (prevdt-dt) == 0.0) ? 0 : 1;
      
    if(nldynamics==0 && T0check>0){ //linear kf time update
      if(verbose > 1) print ("linear update row ", rowi);
    
      if(continuoustime ==1){
        if(dtchange==1 || (T0check == 1 && (DRIFTsubindex + CINTsubindex > 0))){ //if dtchanged or if subject variability
          discreteDRIFT = matrix_exp(append_row(append_col(sDRIFT,sCINT),rep_matrix(0,1,nlatent+1)) * dt);
        }
      
        if(dtchange==1 || (T0check == 1 && (DIFFUSIONsubindex + DRIFTsubindex > 0))){ //if dtchanged or if subject variability
          discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
            quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]\' );
          if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(makesym(discreteDIFFUSION,verbose,1));
        }
      }

      if(continuoustime==0 && T0check == 1){
        if(subjectcount == 1 || DIFFUSIONsubindex + DRIFTsubindex + CINTsubindex > 0){ //if first subject or variability
          discreteDRIFT=append_row(append_col(sDRIFT,sCINT),rep_matrix(0,1,nlatent+1));
          discreteDRIFT[nlatent+1,nlatent+1] = 1;
          discreteDIFFUSION=sDIFFUSIONcov;
          if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(makesym(discreteDIFFUSION,verbose,1));
        }
      }

      Je[rowi] = discreteDRIFT[1:nlatent,1:nlatent];
      state[1:nlatent] = (discreteDRIFT * append_row(state,1.0))[1:nlatent];
      if(ntdpred > 0) state[1:nlatent] += sTDPREDEFFECT * tdpreds[rowi];
      if(intoverstates==1) {
        etacov = quad_form(etacov, discreteDRIFT[1:nlatent,1:nlatent]\');
        if(ndiffusion > 0) etacov += discreteDIFFUSION;
      }
    }//end linear time update


    if(nldynamics==1){ //nldynamics time update
      if(T0check>0){
        vector[nlatentpop] base;
        real intstepi = 0;
        dtsmall = dt / ceil(dt / maxtimestep);
        
        while(intstepi < dt){
          intstepi = intstepi + dtsmall;
          
          ',#simplestatedependencies(when=2,mlist=c(mats$driftcint,mats$diffusion,mats$jacobian[2])),' #now done inside finite diff loop
      finiteJ(),
      simplifystanfunction(paste0(
        paste0(ctm$calcs$diffusion,';\n',collapse=' '), 
        # paste0(c(ctm$calcs$driftcint,ctm$calcs$diffusion),';\n',collapse=' '), #removed because driftcint added during finite diff loop
        jacobianelements(ctm$jacobian$JAx,ntdpred=ctm$n.TDpred,matsetup=matsetup,mats=mx,
          textadd=paste0('    sJAx[',rep(1:nlatentpop,nlatentpop),', ',rep(1:nlatentpop,each=nlatentpop),'] = '),
          when = 2,remove = c('fixed','drift'))
        ,collapse=';')),'
             
             if(verbose>1) print("sJAx ",sJAx);
      for(ri in 1:nlatentpop){ 
        for(ci in 1:nlatentpop){
          if(sJAxdrift[ri,ci]) sJAx[ri,ci]=sDRIFT[ri,ci]; //set jacobian to drift where appropriate
        }
      }
          
          
          if(continuoustime==1){
            if(dtchange==1 || statedependence[2] || (T0check == 1 && (DRIFTsubindex + CINTsubindex > 0))){
              Je[rowi]= matrix_exp(sJAx * dtsmall);
              if(verbose > 1) print("Je = ", Je[rowi]);
              discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),rep_vector(0,nlatent+1)\') * dtsmall,drcintoffdiag);
              if(verbose > 1) print("discreteDRIFT = ", discreteDRIFT);
            } else Je[rowi] = Je[rowi-1]; //temporary hack to avoid nans
            if(dtchange==1 || statedependence[2] || (T0check == 1 && (DRIFTsubindex + DIFFUSIONsubindex + CINTsubindex) > 0)){
              sasymDIFFUSION = to_matrix(  -kronsum(sJAx[1:nlatent,1:nlatent]) \\ to_vector(tcrossprod(sDIFFUSION)), nlatent,nlatent);
              discreteDIFFUSION =  sasymDIFFUSION - quad_form( sasymDIFFUSION, Je[rowi, 1:nlatent,1:nlatent]\' );
            }
            if(verbose>1) print("sJAx ",sJAx);
            if(verbose > 1) print("rowi = ",rowi, "state = ", state);
            if(verbose > 1)  print("etacov = ",etacov," sasymDIFFUSION = ",sasymDIFFUSION," sDIFFUSION = ",sDIFFUSION);
            if(verbose > 1) print("sJAx = ",sJAx);
            etacov = quad_form(etacov, Je[rowi]\');
            etacov[1:nlatent,1:nlatent] += discreteDIFFUSION; //may need improving
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent];
          }

        if(continuoustime==0){ 
          Je[rowi] = sJAx;
          etacov = quad_form(etacov, sJAx\');
          sasymDIFFUSION[ derrind, derrind ] = to_matrix( (IIlatent2 - 
            sqkron_prod(sDRIFT[ derrind, derrind ], sDRIFT[ derrind, derrind ])) \\  to_vector(tcrossprod(sDIFFUSION[ derrind, derrind ])), ndiffusion, ndiffusion);
          etacov[1:nlatent,1:nlatent] += tcrossprod(sDIFFUSION); //may need improving re sDIFFUSION
          discreteDRIFT=append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),rep_matrix(0,1,nlatent+1));
          discreteDRIFT[nlatent+1,nlatent+1] = 1;
          state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent];
        }
      }
    } // end of non t0 time update
  
    if(T0check==0){ //nl t0
    state = sT0MEANS[,1]; //in case of t0 dependencies, may have missingness
    if(verbose > 1) print("rowi = ",rowi, "  state = ",sT0MEANS);
    if(verbose > 1) print("etacov = ",etacov);
    ',simplestatedependencies(when=1,mlist=c(mats$t0,mats$jacobian[1])),'
    ',paste0(ctm$calcs$t0,';',collapse=' '),'
    ',jacobianelements(ctm$jacobian$J0,ntdpred=ctm$n.TDpred,matsetup=matsetup,mats=mx,
      textadd=paste0('    sJ0[',rep(1:nlatentpop,nlatentpop),', ',rep(1:nlatentpop,each=nlatentpop),'] = '),
      when = 1,remove = c('fixed')),'
      state = sT0MEANS[,1];
      etacov= quad_form(sT0VAR, sJ0\');
    if(verbose > 1) print("rowi = ",rowi,"  state = ",sT0MEANS);
    if(verbose > 1) print("sJ0 = ",sJ0);
    if(verbose > 1) print("etacov = ",etacov);
    
    } 
    if(ntdpred > 0) {
    int nonzerotdpred = 0;
    for(tdi in 1:ntdpred) if(tdpreds[rowi,tdi] != 0.0) nonzerotdpred = 1;
    if(nonzerotdpred){
    ',simplestatedependencies(when=3,mlist=c(mats$tdpred,mats$jacobian[3])),'
    ',paste0(ctm$calcs$tdpred,';',collapse=' '),
      jacobianelements(ctm$jacobian$Jtd,ntdpred=ctm$n.TDpred,matsetup=matsetup,mats=mx,
        textadd=paste0('    sJtd[',rep(1:nlatentpop,nlatentpop),', ',rep(1:nlatentpop,each=nlatentpop),'] = '),
        when = 3,remove = c('fixed')),'
      state[1:nlatent] +=   (sTDPREDEFFECT * tdpreds[rowi]); //tdpred effect only influences at observed time point','
      if(verbose > 1)  print("state = ", state);
      if(verbose > 1)  print("etacov = ",etacov);
 if(verbose > 1) print("sJtd = ",sJtd);
      etacov = quad_form(etacov,sJtd\');
    }
    }
  } // end non linear time update


  if(savescores==1) {
    kout[rowi,(nmanifest*4+1):(nmanifest*4+nlatent)] = state[1:nlatent];
    etaprior[rowi] = state;
    etapriorcov[rowi]=etacov;
    if(nobs_y[rowi] == 0) etaupdcov[rowi]=etacov;
  }
  if(verbose > 1) print("etaprior = ", state, " etapriorcov = ",etacov);

    if(intoverstates==0 && nldynamics == 0) {
      if(T0check==0) state += cholesky_decompose(sT0VAR) * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)];
      if(T0check>0) state +=  discreteDIFFUSION * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)];
    }


    if(nlmeasurement==0 && T0check == 0) sJy[ ,1:nlatent] = sLAMBDA;
    if (nobs_y[rowi] > 0 || savescores) {  // if some observations create right size matrices for missingness and calculate...
    
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
      

      if(nlmeasurement==1){
      ',simplestatedependencies(when=4,mlist=c(mats$measurement,mats$jacobian[4])),'
      ',paste0(ctm$calcs$measurement,';',collapse=' '),
      jacobianelements(ctm$jacobian$Jy,matsetup=matsetup,mats=mx,ntdpred=0,
        textadd=paste0('    sJy[',rep(1:nmanifest,nlatentpop),', ',rep(1:nlatentpop,each=nmanifest),'] = '),
        when = 4,remove = c('fixed','lambda')),'
        
        for(ri in 1:nmanifest){ 
          for(ci in 1:nlatentpop){
            if(sJylambda[ri,ci]) sJy[ ri,ci]=sLAMBDA[ri,ci]; //set jacobian to lambda where appropriate
          }
        }
      }
          
        if(intoverstates==1) { //classic kalman
          yprior[o] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * state[1:nlatent];
          if(nbinary_y[rowi] > 0) yprior[o1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * state[1:nlatent])));
          if(verbose > 1) print ("sMANIFESTVAR[o,o] = ",sMANIFESTVAR[o,o])
          if(verbose > 1) print ("etacov[1:nlatent,1:nlatent] = ",etacov[1:nlatent,1:nlatent])
          if(verbose > 1) print ("sJy[o,]\' = ",sJy[o,]\');
          ycov[o,o] = quad_form(etacov, sJy[o,]\') + sMANIFESTVAR[o,o];
          for(wi in 1:nmanifest){ 
            if(manifesttype[wi]==1 && Y[rowi,wi] != 99999) ycov[wi,wi] += fabs((yprior[wi] - 1) .* (yprior[wi]));
            if(manifesttype[wi]==2 && Y[rowi,wi] != 99999) ycov[wi,wi] += square(fabs((yprior[wi] - round(yprior[wi])))); 
          }
        }
        
        if(intoverstates==0) { //sampled states
          if(ncont_y[rowi] > 0) {
            yprior[o0] = sMANIFESTMEANS[o0,1] + sJy[o0,] * state;
            ypriorcov_sqrt[o0,o0] = sqrt(sMANIFESTVAR[o0,o0]);
          }
          if(nbinary_y[rowi] > 0) yprior[o1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * state[1:nlatent])));
        }
        
     
'    
      ,if(ppchecking) paste0('
  {
  int skipupd = 0;
        for(vi in 1:nobs_y[rowi]){
            if(fabs(yprior[od[vi]]) > 1e10 || is_nan(yprior[od[vi]]) || is_inf(yprior[od[vi]])) {
              skipupd = 1; 
              yprior[od[vi]] =99999;
  if(verbose > 1) print("pp yprior problem! row ", rowi);
            }
          }
        if(skipupd==0){ 
          if(ncont_y[rowi] > 0){
            ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(makesym(ycov[o0d, o0d],verbose,1)); 
            Ygen[ rowi, o0d] = yprior[o0d] + ypriorcov_sqrt[o0d,o0d] * Ygenbase[rowi,o0d];
          }
          if(nbinary_y[rowi] > 0) for(obsi in 1:size(o1d)) Ygen[rowi, o1d[obsi]] = (yprior[o1d[obsi]] > Ygenbase[rowi,o1d[obsi]]) ? 1 : 0; 
          for(vi in 1:nobs_y[rowi]) if(is_nan(Ygen[rowi,od[vi]])) {
            Ygen[rowi,od[vi]] = 99999;
print("pp ygen problem! row ", rowi);
          }
        if(nlmeasurement==0){ //linear measurement
          if(intoverstates==1) { //classic kalman
            for(wi in 1:nmanifest){ 
              if(manifesttype[wi]> 0 && Y[rowi,wi] != 99999) Ygen[ rowi, wi] = round(Ygen[ rowi, wi]);
            }
          }
        }
        err[od] = Ygen[rowi,od] - yprior[od]; // prediction error
        }
}
      '), 
      
      if(!ppchecking) 'err[od] = Y[rowi,od] - yprior[od]; // prediction error','
    
      if(intoverstates==1 && size(od) > 0) {
        K[,od] = mdivide_right(etacov * sJy[od,]\', ycov[od,od]); 
        etacov += -K[,od] * sJy[od,] * etacov;
        state +=  (K[,od] * err[od]);
      }
      
      if(savescores==1) {
        int tmpindex[nmanifest] = o;
        for(oi in 1:nmanifest) tmpindex[oi] +=  nmanifest*2;
        kout[rowi,tmpindex] = err[o];
        for(oi in 1:nmanifest) tmpindex[oi] +=  nmanifest;
        kout[rowi,tmpindex] = yprior[o];
        etaupd[rowi] = state;
        ypriorcov[rowi] = ycov;
        etaupdcov[rowi] = etacov;
        yupdcov[rowi] = quad_form(etacov, sJy\') + sMANIFESTVAR;
        yupd[rowi] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * state[1:nlatent];
        ysmoothcov[rowi] = sMANIFESTVAR; // add the rest later
        ysmooth[rowi] = sMANIFESTMEANS[,1]; // add the rest later
        Jy[rowi] = sJy;
        tLAMBDA[rowi] = sLAMBDA;
      }
      
      
      if(verbose > 1) {
          print("rowi =",rowi, "  si =", si, "  state =",state,"  etacov ",etacov,
            "  yprior =",yprior,"  ycov ",ycov, "  K ",K,
            "  sDRIFT =", sDRIFT, " sDIFFUSION =", sDIFFUSION, " sCINT =", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
            "  sT0VAR", sT0VAR,  " sT0MEANS ", sT0MEANS, "sLAMBDA = ", sLAMBDA, "  sJy = ",sJy,
            " discreteDRIFT = ", discreteDRIFT, "  discreteDIFFUSION ", discreteDIFFUSION, "  sasymDIFFUSION ", sasymDIFFUSION, 
            "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        }
  
      ',if(!ppchecking){
        'if(nbinary_y[rowi] > 0) kout[rowi,o1d] =  Y[rowi,o1d] .* (yprior[o1d]) + (1-Y[rowi,o1d]) .* (1-yprior[o1d]); 
  
        if(size(o0d) > 0){
          int tmpindex[ncont_y[rowi]] = o0d;
          for(oi in 1:ncont_y[rowi]) tmpindex[oi] +=  nmanifest;
           if(intoverstates==1) ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(makesym(ycov[o0d,o0d],verbose,1));
           //ll+=  multi_normal_cholesky_lpdf(Y[rowi] | yprior, ypriorcov_sqrt[o0d,o0d]);
           kout[rowi,o0d] = mdivide_left_tri_low(ypriorcov_sqrt[o0d,o0d], err[o0d]); //transform pred errors to standard normal dist and collect
           kout[rowi,tmpindex] = log(diagonal(ypriorcov_sqrt[o0d,o0d])); //account for transformation of scale in loglik
        }
      '},'
    }//end nobs > 0 section
  if(savescores==1) {
    kout[rowi,(nmanifest*4+nlatent+1):(nmanifest*4+nlatent+nlatent)] = state[1:nlatent];
  }','
  
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
      ysmoothcov[sri] += quad_form(etasmoothcov[sri], Jy[sri]\'); //already added manifestvar
      ysmooth[sri] += tLAMBDA[sri] * etasmooth[sri,1:nlatent];
      sri = sri -1;
    }
  } //end smoother
  
  } // end dokalmanrows subset selection
}//end rowi
')
    if(!is.null(ctm$w32)) out <- ''
    return(out)}

kalmanll <- function(){ out <-'
  if(sum(nbinary_y) > 0) {
    vector[sum(nbinary_y)] binaryll;
    counter = 1;
    for(ri in 1:ndatapoints){
    if(dokalmanrows[ri]==1){
      int o1[nbinary_y[ri]] = whichbinary_y[ri,1:nbinary_y[ri]]; //which indicators are observed and binary
      binaryll[counter:(counter + nbinary_y[ri]-1)] = kout[ri,o1];
      counter+= nbinary_y[ri];
    }
    }
    ll+= sum(log(binaryll[1:(counter-1)]));
  }

  if(( sum(ncont_y) > 0)) {
    vector[sum(ncont_y)] errtrans[2];
    counter = 1;
    for(ri in 1:ndatapoints){
    if(dokalmanrows[ri]==1){
      int o0[ncont_y[ri]] = whichcont_y[ri,1:ncont_y[ri]]; //which indicators are observed and continuous
      errtrans[1,counter:(counter + ncont_y[ri]-1)] = kout[ri, o0];
      for(oi in 1:ncont_y[ri]) o0[oi] +=  nmanifest; //modify o0 index
      errtrans[2,counter:(counter + ncont_y[ri]-1)] = kout[ri, o0];
      counter+= ncont_y[ri];
    }
    }
    ll += normal_lpdf(errtrans[1,1:(counter-1)]|0,1) - sum(errtrans[2,1:(counter-1)]);
  }
if(savescores) kalman = kout;
'
if(!is.null(ctm$w32)) out <- ''
return(out)}

subjectparaminit<- function(popmats=FALSE,smats=TRUE){
  if(smats && popmats) stop('smats and popmats cannot both be TRUE!')
  paste0(
    paste0('   matrix[matrixdims[',(mats$base),', 1], matrixdims[',(mats$base),', 2] ] ',
      ifelse(smats,'s',''),ifelse(popmats,'pop_',''),names(mats$base),if(!smats && !popmats) paste0('[',names(mats$base),'subindex  ? nsubjects : 1]'),';',collapse=' \n   '),'

  matrix[nlatent,nlatent] ',ifelse(smats,'s',''),ifelse(popmats,'pop_',''),'asymDIFFUSION',if(!smats && !popmats) '[asymDIFFUSIONsubindex ? nsubjects : 1]','; //stationary latent process variance
  vector[nlatent] ',ifelse(smats,'s',''),ifelse(popmats,'pop_',''),'asymCINT',if(!smats && !popmats) '[asymCINTsubindex ? nsubjects : 1]','; // latent process asymptotic level
','matrix[nlatent, nlatent] ',ifelse(smats,'s',''),ifelse(popmats,'pop_',''), 'DIFFUSIONcov',if(!smats && !popmats) '[DIFFUSIONcovsubindex ? nsubjects : 1]',';',
    collapse='\n')
}

stateparscalc <- function(){
  paste0('
    for(ri in 1:size(matsetup)){ //for each row of matrix setup
        real newval;
        if(matsetup[ri,3] > 0) newval = 
          tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 

        ',paste0('if(matsetup[ri, 7] == ', (mats$base),') s',
          names(mats$base),'[matsetup[ ri,1], matsetup[ri,2]] = newval;', collapse = ' \n      '),'
    }',collapse='\n    ')
}

subjectparscalc2 <- function(popmats=FALSE,subjmats=TRUE){
  out <- paste0(
    '
 int subjectvec[subjectcount ? 1 : 2];
 subjectvec[size(subjectvec)] = si;
 if(subjectcount == 0)  subjectvec[1] = 0; // only needed for subject 0 (pop pars)
 subjectcount = subjectcount + 1;
 for(subjectveci in 1:size(subjectvec)){
  int subi = subjectvec[subjectveci];
  vector[nparams] rawindparams;
  vector[nparams] tipredaddition = rep_vector(0,nparams);
  vector[nparams] indvaraddition = rep_vector(0,nparams);

  if(subi > 0 && nindvarying > 0 && intoverpop==0) {
    if(fixedsubpars==0) indvaraddition[indvaryingindex] = rawpopcovsqrt * baseindparams[subi];
    if(fixedsubpars==1) indvaraddition[indvaryingindex] = rawpopcovsqrt * fixedindparams[subi];
  }
  
  if(subi > 0 &&  ntipred > 0) tipredaddition = TIPREDEFFECT * tipreds[subi]\';

  rawindparams = rawpopmeans + tipredaddition + indvaraddition;
',
    
    paste0('
    for(ri in 1:size(matsetup)){ //for each row of matrix setup
    for(statecalcs in 0:1){
        if(subi ==0 ||  //if population parameter
          (matsetup[ri,7]==4 && DIFFUSIONsubindex) ||( matsetup[ri,7] == 8 && T0VARsubindex) || //or a covariance parameter in an individually varying matrix
          (matsetup[ri,3] > 0 && (matsetup[ri,5] > 0 || matsetup[ri,6] > 0)) //or there is individual variation
          ){ //otherwise repeated values
            if( (statecalcs && matsetup[ri,8]>0) || (!statecalcs && matsetup[ri,8]==0) ){ //if doing statecalcs do them, if doing static calcs do them
              real newval;
              if(matsetup[ri,3] > 0)  newval = tform(matsetup[ri,8] ? state[ matsetup[ri,3] ] : rawindparams[ matsetup[ri,3] ], //tform static pars from rawindparams, dynamic from state
                matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
               if(matsetup[ri,3] < 1) newval = matvalues[ri, 1]; //doing this once over all subjects unless covariance matrix -- speed ups possible here, check properly!
              ',paste0('if(matsetup[ri, 7] == ', c(mats$base,mats$jacobian),') s',
                names(c(mats$base,mats$jacobian)),'[matsetup[ ri,1], matsetup[ri,2]] = newval;', collapse = ' \n      '),'
            }
          }
        state=sT0MEANS[,1];
      }
    }',collapse='\n    '),'

  // perform any whole matrix transformations, nonlinear calcs based on t0 in order to fill matrices
  ',paste0(ctm$calcs$t0,';\n',collapse=' '),'; 
  state=sT0MEANS[,1];
  ',paste0(ctm$calcs$tdpred,';\n',collapse=' '),';
  ',paste0(ctm$calcs$driftcint,';\n',collapse=' '),';
  ',paste0(ctm$calcs$diffusion,';\n',collapse=' '),';
  ',paste0(ctm$calcs$measurement,';\n',collapse=' '),';
  
  if(subi <= (DIFFUSIONsubindex ? nsubjects : 0)) {
    sDIFFUSIONcov = sdcovsqrt2cov(sDIFFUSION,nldynamics);
  }
  if(subi <= (asymDIFFUSIONsubindex ? nsubjects : 0)) {
      if(ndiffusion < nlatent) sasymDIFFUSION = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

      if(continuoustime==1) sasymDIFFUSION[ derrind, derrind] = to_matrix( 
      -kronsum(sDRIFT[ derrind, derrind ]) \\  to_vector( 
           sDIFFUSIONcov[ derrind, derrind ]), ndiffusion,ndiffusion);

      if(continuoustime==0) sasymDIFFUSION[ derrind, derrind ] = to_matrix( (IIlatent2 - 
        sqkron_prod(sDRIFT[ derrind, derrind ], sDRIFT[ derrind, derrind ])) \\  to_vector(sDIFFUSIONcov[ derrind, derrind ]), ndiffusion, ndiffusion);
    } //end asymdiffusion loops

      if(subi <= (MANIFESTVARsubindex ? nsubjects : 0)) {
         for(ri in 1:nmanifest) sMANIFESTVAR[ri,ri] = square(sMANIFESTVAR[ri,ri]);
      }
         
    if(subi <= (T0VARsubindex ? nsubjects : 0)) {
    if(intoverpop){
      sT0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopcovsqrt;
    }
      sT0VAR = makesym(sdcovsqrt2cov(sT0VAR,nldynamics),verbose,1);
    if(intoverpop){ //adjust cov matrix for transforms
      for(ri in 1:size(matsetup)){
        if(matsetup[ri,7]==1){ //if t0means
          if(matsetup[ri,5]) { //and indvarying
            sT0VAR[matsetup[ri,1], ] = sT0VAR[matsetup[ri,1], ] * matvalues[ri,2] * matvalues[ri,3]* matvalues[ri,5]; //multiplier meanscale sdscale
            sT0VAR[, matsetup[ri,1] ] = sT0VAR[, matsetup[ri,1] ] * matvalues[ri,2] * matvalues[ri,3]* matvalues[ri,5]; //multiplier meanscale sdscale
          }
        }
      }
    }
      if(nt0varstationary > 0) {
        for(ri in 1:nt0varstationary){ 
          sT0VAR[t0varstationary[ri,1],t0varstationary[ri,2] ] =  sasymDIFFUSION[t0varstationary[ri,1],t0varstationary[ri,2] ];
        }
      }
    }
    
      if(subi <= (asymCINTsubindex ? nsubjects : 0)){
        if(continuoustime==1) sasymCINT =  -sDRIFT[1:nlatent,1:nlatent] \\ sCINT[ ,1 ];
        if(continuoustime==0) sasymCINT =  (IIlatent - sDRIFT[1:nlatent,1:nlatent]) \\ sCINT[,1 ];
      }
    
    if(nt0meansstationary > 0){
      if(subi <= (T0MEANSsubindex ? nsubjects : 0)) {
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
  out<-'if(savesubjectmatrices){'
  for(m in matrices){
    if(!popmats) out <- paste0(out, 'if( (', m,'subindex > 0 && subi > 0) || (',m,'subindex == 0 && subi==0) ) ',m,'[',m,'subindex ? subi : 1] = s',m,'; \n')
    if(popmats) out <- paste0(out, 'pop_',m,' = s',m,'; \n }')
  }
  
  return(out)
}



writemodel<-function(){
  paste0('
functions{
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
    out[z0,] = rep_matrix(0,size(z0),rows(M));
    out[,z0] = rep_matrix(0,rows(M),size(z0));
    if(size(z0) > 0) for(i in 1:size(z0)) out[z0[i],z0[i]] = exp(M[z0[i],z0[i]]);
    return out;
  }

   matrix constraincorsqrt(matrix mat){ //converts from unconstrained lower tri matrix to cor
    matrix[rows(mat),cols(mat)] o;
    vector[rows(mat)] s;
  
    for(i in 1:rows(o)){ //set upper tri to lower
      for(j in min(i+1,rows(mat)):rows(mat)){
        o[j,i] =  inv_logit(mat[j,i])*2-1;  // can change cor prior here
        o[i,j] = o[j,i];
      }
      o[i,i]=1; // change to adjust prior for correlations
    }

    for(i in 1:rows(o)){
      s[i] = inv_sqrt(o[i,] * o[,i]);
      if(is_inf(s[i])) s[i]=0;
    }
    return diag_pre_multiply(s,o);
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

  matrix kronsum(matrix mata){
    matrix[rows(mata),rows(mata)] II = diag_matrix(rep_vector(1,rows(mata)));
    return sqkron_prod(mata, II) + sqkron_prod(II, mata );
  }

  matrix makesym(matrix mat, int verbose, int pd){
    matrix[rows(mat),cols(mat)] out;
    for(coli in 1:cols(mat)){
      if(pd ==1){ // && mat[coli,coli] < 1e-5
        //if(verbose > 0) print("diagonal too low (",mat[coli,coli],") during makesym row ", coli, " col ", coli);
        out[coli,coli] = mat[coli,coli] + 1e-5;
      } else out[coli,coli] = mat[coli,coli]; 
      for(rowi in coli:rows(mat)){
        if(rowi > coli) {
          out[rowi,coli] = mat[rowi,coli]; //(mat[coli,rowi] + ) *.5;
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

  real tform(real param, int transform, data real multiplier, data real meanscale, data real offset, data real inneroffset){
    real out;
    ',tformshapes(),tformshapes(jacobian=TRUE),
    if(length(extratforms) > 0) paste0(extratforms,collapse=" \n"),'
    return out;
  }
  
  
  real Jtform(real param, int transform, data real multiplier, data real meanscale, data real offset, data real inneroffset){
    real out;
    ',tformshapes(jacobian=TRUE),
    if(length(extratforms) > 0) paste0(extratforms,collapse=" \n"),'
    return out;
  }
  
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

}
data {
  int<lower=0> ndatapoints;
  int<lower=1> nmanifest;
  int<lower=1> nlatent;
  int nlatentpop;
  int<lower=1> nsubjects;
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
  int ndiffusion; //number of latents involved in covariance calcs
  int derrind[ndiffusion]; //index of which latent variables are involved in covariance calculations
  int drcintoffdiag[nlatent+1];

  int manifesttype[nmanifest];
  int nbinary_y[ndatapoints];  // number of observed binary variables per observation
  int whichbinary_y[ndatapoints, nmanifest]; // index of which variables are observed and binary per observation
  int ncont_y[ndatapoints];  // number of observed continuous variables per observation
  int whichcont_y[ndatapoints, nmanifest]; // index of which variables are observed and continuous per observation
  
  int intoverpop;
  int statedependence[4];
  real ukfspread;
  int ukffull;
  int nlmeasurement;
  int intoverstates;
  int verbose; //level of printing during model fit

  ',paste0(unlist(lapply(c(names(mats$base),'asymCINT','asymDIFFUSION','DIFFUSIONcov'),function(mati) paste0('int ',mati,'subindex;',collapse='\n'))),collapse='\n'),'
  int TIPREDEFFECTsetup[nparams, ntipred];
  int nrowmatsetup;
  int matsetup[nrowmatsetup,8];
  real matvalues[nrowmatsetup,6];
  int matrixdims[',length(mats$base),',2];
  int savescores;
  int savesubjectmatrices;
  int fixedsubpars;
  vector[fixedsubpars ? nindvarying : 0] fixedindparams[fixedsubpars ? nsubjects : 0];
  int dokalman;
  int dokalmanrows[ndatapoints];
  real Jstep;
  real dokalmanpriormodifier;
  int intoverpopindvaryingindex[intoverpop ? nindvarying : 0];
  int sJAxdrift[nlatentpop,nlatentpop];
  int sJylambda[nmanifest,nlatentpop];
  int nsJAxfinite;
  int sJAxfinite[nsJAxfinite];
}
      
transformed data{
  matrix[nlatent,nlatent] IIlatent= diag_matrix(rep_vector(1,nlatent));
  matrix[nlatent*nlatent,nlatent*nlatent] IIlatent2 = diag_matrix(rep_vector(1,nlatent*nlatent));
  matrix[nindvarying,nindvarying] IIindvar = diag_matrix(rep_vector(1,nindvarying));
}
      
parameters {
  vector[nparams] rawpopmeans; // population level means \n','
  vector',if(!is.na(ctm$rawpopsdbaselowerbound)) paste0('<lower=',ctm$rawpopsdbaselowerbound[1],'>'),'[nindvarying] rawpopsdbase; //population level std dev
  vector[nindvaryingoffdiagonals] sqrtpcov; // unconstrained basis of correlation parameters
  vector[fixedsubpars ? 0 : (intoverpop ? 0 : nindvarying)] baseindparams[fixedsubpars ? 0 : (intoverpop ? 0 : nsubjects)]; //vector of subject level deviations, on the raw scale
  
  vector[ntipredeffects] tipredeffectparams; // effects of time independent covariates
  vector[nmissingtipreds] tipredsimputed;
  //vector[ (( (ntipredeffects-1) * (1-nopriors) ) > 0) ? 1 : 0] tipredglobalscalepar;
  
  vector[intoverstates ? 0 : nlatentpop*ndatapoints] etaupdbasestates; //sampled latent states posterior
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying, nindvarying] rawpopcovsqrt; 
',if(!gendata) paste0('
  real ll = 0;
  vector[nmanifest+nmanifest+ (savescores ? nmanifest*2+nlatent*2 : 0)] kalman[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etapriorcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etaupdcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etasmoothcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ypriorcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] yupdcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ysmoothcov[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaprior[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaupd[savescores ? ndatapoints : 0];
  vector[nlatentpop] etasmooth[savescores ? ndatapoints : 0];
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
    rawpopsd = ',ctm$rawpopsdtransform, '; // sqrts of proportions of total variance
    for(j in 1:nindvarying){
      rawpopcovsqrt[j,j] = 1;
      for(i in 1:nindvarying){
        if(i > j){
          counter += 1;
          rawpopcovsqrt[i,j]=sqrtpcov[counter];
          rawpopcovsqrt[j,i]=sqrtpcov[counter];
        }
      }
    }
 rawpopcovsqrt = cholesky_decompose(makesym(tcrossprod(diag_pre_multiply(rawpopsd, 
      constraincorsqrt(rawpopcovsqrt))),verbose,1)); 
  }//end indvarying par setup

  {',
if(!gendata) ukfilterfunc(ppchecking=FALSE),
'if(dokalman==1){',
if(!gendata) kalmanll(),'
    }
  }
}
      
model{
  if(intoverpop==0 && fixedsubpars == 1) target+= multi_normal_cholesky_lpdf(fixedindparams | rep_vector(0,nindvarying),IIindvar);

  if(nopriors==0 && subject[1]==1){ //if split files over subjects, just compute priors once
   target+= dokalmanpriormodifier * normal_lpdf(rawpopmeans|0,1);
  
    if(ntipred > 0){ 
      target+= dokalmanpriormodifier * normal_lpdf(tipredeffectparams| 0, tipredeffectscale);
      target+= normal_lpdf(tipredsimputed| 0, tipredsimputedscale); //consider better handling of this when using subset approach
      //target+= dokalmanpriormodifier * normal_lpdf(tipredglobalscalepar | 0-log(ntipred),log(square(ntipred)));
    }
    
    if(nindvarying > 0){
      if(nindvarying >1) target+= dokalmanpriormodifier * normal_lpdf(sqrtpcov | 0, 1);
      if(intoverpop==0 && fixedsubpars == 0) target+= multi_normal_cholesky_lpdf(baseindparams | rep_vector(0,nindvarying), IIindvar);
      target+= dokalmanpriormodifier * normal_lpdf(rawpopsdbase | ',gsub('normal(','',ctm$rawpopsdbase,fixed=TRUE),';
    }
    //llp +=  log(dokalmanpriormodifier);
  } //end pop priors section
  
  if(intoverstates==0) target+= normal_lpdf(etaupdbasestates|0,1);
  
  ',if(!gendata) 'target+= ll; \n','
  if(verbose > 0) print("lp = ", target());
}
generated quantities{
  vector[nparams] popmeans;
  vector[nparams] popsd = rep_vector(0,nparams);
  matrix[nindvarying,nindvarying] rawpopcov = tcrossprod(rawpopcovsqrt);
  matrix[nindvarying,nindvarying] rawpopcorr = quad_form_diag(rawpopcov,inv_sqrt(diagonal(rawpopcov)));
  matrix[nparams,ntipred] linearTIPREDEFFECT;
',if(gendata) paste0('
  real ll = 0;
  vector[nmanifest+nmanifest+ (savescores ? nmanifest*2+nlatent*2 : 0)] kalman[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etapriorcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etaupdcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etasmoothcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ypriorcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] yupdcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ysmoothcov[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaprior[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaupd[savescores ? ndatapoints : 0];
  vector[nlatentpop] etasmooth[savescores ? ndatapoints : 0];
  vector[nmanifest] yupd[savescores ? ndatapoints : 0];
  vector[nmanifest] ysmooth[savescores ? ndatapoints : 0];
  vector[nmanifest] Ygen[ndatapoints];
  ',subjectparaminit(pop=FALSE,smats=FALSE),'
  ',subjectparaminit(pop=TRUE,smats=FALSE)
  ,collapse=''),'

{
vector[nparams] rawpopsdfull;
rawpopsdfull[indvaryingindex] = sqrt(diagonal(rawpopcov)); //base for calculations

    for(ri in 1:size(matsetup)){
      if(matsetup[ri,3] && matsetup[ri,8]==0) { //if a free parameter 
        real rawpoppar = rawpopmeans[matsetup[ri,3] ];
        int pr = ri; // unless intoverpop, pop matrix row reference is simply current row
        
        if(intoverpop && matsetup[ri,5]) { //removed ri transform of rawpop because t0means only transforms once -- if non identity state tform in future, change this!
          for(ri2 in 1:size(matsetup)){ //check when state reference param of matsetup corresponds to row of t0means in current matsetup row
            if(matsetup[ri2,8]  && matsetup[ri2,3] == matsetup[ri,1]) pr = ri2;
            //print("ri = ",ri, " pr = ",pr, " ri2 = ",ri2);
          }
        }
        
        popmeans[matsetup[ ri,3]] = tform(rawpoppar, matsetup[pr,4], matvalues[pr,2], matvalues[pr,3], matvalues[pr,4], matvalues[pr,6] ); 

        popsd[matsetup[ ri,3]] = matsetup[ ri,5] ? //if individually varying
          fabs(tform( //compute sd
            rawpoppar  + rawpopsdfull[matsetup[ ri,3]], matsetup[pr,4], matvalues[pr,2], matvalues[pr,3], matvalues[pr,4], matvalues[pr,6]) -
           tform(
            rawpoppar  - rawpopsdfull[matsetup[ ri,3]], matsetup[pr,4], matvalues[pr,2], matvalues[pr,3], matvalues[pr,4], matvalues[pr,6]) ) /2 : 
          0; //else zero

        if(ntipred > 0){
          for(tij in 1:ntipred){
            if(TIPREDEFFECTsetup[matsetup[ri,3],tij] ==0) {
              linearTIPREDEFFECT[matsetup[ri,3],tij] = 0;
            } else {
            linearTIPREDEFFECT[matsetup[ri,3],tij] = ( //tipred reference is from row ri, tform reference from row pr in case of intoverpop
              tform(rawpoppar + TIPREDEFFECT[matsetup[ri,3],tij] * .01, matsetup[pr,4], matvalues[pr,2], matvalues[pr,3], matvalues[pr,4], matvalues[pr,6] ) -
              tform(rawpoppar - TIPREDEFFECT[matsetup[ri,3],tij] * .01, matsetup[pr,4], matvalues[pr,2], matvalues[pr,3], matvalues[pr,4], matvalues[pr,6] )
              ) /2 * 100;
            }
         }
        }
      }
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
