# ctm=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','state[3]',0, 'd12','d22',0,0,0,0),3,3),type='stanct')
# ctm2=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','PARS[1,1]',0, 'd12','d22',0,0,0,0),3,3),PARS=matrix('state[3]'),type='stanct')

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

#replace inverse logit
inv_logit_gsub <- function(x){
  try=0
  while(try < 20 && any(grepl('\\<inv_logit\\((.*)\\)',x) | grepl('\\blog1p_exp\\((.*)\\)',x))){ 
    try <- try + 1
    x <- gsub('\\<inv_logit\\((.*)\\)','1/\\(1+exp\\(-\\(\\1\\)\\)\\)',x)
    x <- gsub('\\<log1p_exp\\((.*)\\)','log1p\\(exp\\(\\1\\)\\)',x)
  }
  return(x)
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
  
  m$pars$param <- inv_logit_gsub(m$pars$param) #replace inv_logit with known functions for differentiation
  
  mats <- listOfMatrices(m$pars)
  matnames <- names(ctStanMatricesList(unsafe=TRUE)$base)

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
    
    
    # browser()
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
    
    # 3): calculate Jacobian of fn symbolically
    J  = jacobianSymb(fn, state)
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
    
    #this was creating direct references but not needed here
    # for(j in 1:ncol(Jm)){
    #   for(i in 1:nrow(Jm)){
    #     if(is.na(suppressWarnings(as.numeric(Jm[i,j])))){
    #       for(mi in 1:length(mats)){
    #         if(Jm[i,j] %in% mats[[mi]]){ 
    #           arrind <- arrayInd(which(mats[[mi]] %in% Jm[i,j]),dim(mats[[mi]]))
    #           for(ari in 1:nrow(arrind)){
    #             Jm[i,j] <- paste0(names(mats)[mi],'[',arrind[ari,1],',',arrind[ari,2],']') #removed 's',
    #           }
    #         }
    #       }
    #     }
    #   }
    # }
    Jout[[typei]] <- Jm
  }#end type loop
  return(Jout)
}

