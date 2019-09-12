# ctm=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','state[3]',0, 'd12','d22',0,0,0,0),3,3),type='stanct')
# ctm2=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','PARS[1,1]',0, 'd12','d22',0,0,0,0),3,3),PARS=matrix('state[3]'),type='stanct')

ctJacobian <- function(m,types=c('J0','JAx','Jtd','Jy') ){
  
  # require(cOde)
  # require(Deriv)
  
  # get system dimension
  ndim = max(m$pars$row[m$pars$matrix%in% 'T0MEANS'])
  
  # 2): generate vector valued function fn = drift * state
  
  # initialize fn and state
  fn     = c()
  state   = paste0("state__", 1:ndim,'__')

  #replace system matrix references
  matlist <- listOfMatrices(m$pars)
  matnames <- names(ctStanMatricesList(unsafe=TRUE)$base)
 
  
  # for(mati in matnames){ 
  for(ri in 1:nrow(m$pars)){
    counter <- 0
    if(grepl('^\\b(state)\\b\\[\\d+\\]$',m$pars$param[ri]) && 
        !is.na(as.integer(m$pars$transform[ri]))){ #if a simple state reference, and integer transform, include transforms
      m$pars$param[ri] <- tform(
        param = m$pars$param[ri],
        transform = as.integer(m$pars$transform[ri]),
        multiplier = m$pars$multiplier[ri],
        meanscale = m$pars$meanscale[ri],
        offset = m$pars$offset[ri],
        inneroffset = m$pars$inneroffset[ri],singletext=TRUE)
      # 
      # m$pars$param[ri] <- paste0('Jtform___',m$pars$param[ri],' * ',m$pars$param[ri])
    }
    while(counter < 10 && grepl(paste0('\\b(',paste0(names(matlist),collapse='|'),')\\b\\['), m$pars$param[ri])){ #if system matrix referenced, unfold
      counter = counter +1 #prevent infinite loops
      
      items = regmatches(m$pars$param[ri], #extract one or more references
        gregexpr(
          paste0('\\b(',paste0(names(matlist),collapse='|'),')\\b(?=\\[).*?(?<=\\])'), 
          m$pars$param[ri], perl=TRUE)
      )[[1]]
      
      for(itemi in 1:length(items)){ #replace one or more references
        m$pars$param[ri] <- gsub(pattern = items[itemi], replacement = eval(parse(text=paste0('matlist$',items[itemi]))),x = m$pars$param[ri], fixed=TRUE)
      }
    }
  }
  
  mats<-listOfMatrices(m$pars)
  Jout <- list()
  for(typei in types){
    if(typei=='JAx'){
      Jrows = nrow(mats$T0MEANS)
      for (row in 1:ndim) {
        for (col in 1:(ndim)) {
          fn[row] = paste0(ifelse(col > 1, paste0(fn[row],' + '),''), "(", mats$DRIFT[row, col], ") * state[", as.character(col), "]")
        }
        if(!is.na(mats$CINT[row])) fn[row] = paste0(fn[row],' + ',mats$CINT[row]) #checking for NA because CINT is not always as large as DRIFT
      }
      # 
    }

    
    if(typei=='Jtd'){
      if(m$n.TDpred ==0) {
        Jout[[typei]] <- c('')
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
        if(is.na(suppressWarnings(as.numeric(t0func[xi]))) && !grepl('[',t0func[xi],fixed=TRUE)) out <- 
          paste0('state[',xi,']') else out <- t0func[xi]
          return(out)
      })
      fn = sapply(prodSymb(diag(nrow(mats$T0MEANS)), matrix(t0func,ncol=1)),Simplify)
    }
    if(typei=='Jy') {
      Jrows = nrow(mats$MANIFESTMEANS)
      Jybase <- mats$LAMBDA
      Jybase <- cbind(Jybase,matrix(0,nrow(Jybase),ncol=nrow(mats$T0MEANS)-ncol(Jybase)))
      fn = paste0(mats$MANIFESTMEANS,' + ',prodSymb(Jybase,cbind(paste0('state[',1:ndim,']'))))
      # fn = sapply(prodSymb(diag(nrow(mats$T0MEANS)), matrix(t0func,ncol=1)),Simplify)
    }
    fn = gsub(" ", "", fn, fixed = TRUE) #remove spaces
    # replace state[~] by state~ for cOde Jacobian and make fn and state a named list
    names(fn) = paste0("fn", 1:length(fn))
    for(statei in 1:ndim){
      fn=gsub(paste0('\\b(state)\\[(',statei,')?\\]'),paste0('state__',statei,'__'),fn,perl = TRUE)
    }
    
    fn=gsub(' ','',sapply(fn,Simplify),fixed=TRUE)
    #probably redundant now but maybe useful at some point?
    # replace remaining commas and square brackets for cOde Jacobian
    fn = gsub(",", "___comma___", fn, fixed = TRUE)
    fn = gsub("[", "___leftsquarebracket___", fn, fixed = TRUE)
    fn = gsub("]", "___rightsquarebracket___", fn, fixed = TRUE)
    fn=sapply(fn,Simplify)

    # 3): calculate Jacobian of fn symbolically
    J  = jacobianSymb(fn, state)
    # 4): create Jacobian list in STAN format
    # J = sapply(J,Simplify)

    J = gsub(" ", "", J, fixed = TRUE) #remove spaces
    
    for(statei in 1:ndim){
      J=gsub(paste0('state__',statei,'__'),paste0('state[',statei,']'),J,fixed = TRUE)
    }
    
    # restore commas and square brackets for Jacobian list
    J = gsub("___rightsquarebracket___", "]", J, fixed = TRUE)
    J = gsub("___leftsquarebracket___", "[", J, fixed = TRUE)
    J = gsub("___comma___", ",", J, fixed = TRUE)
    
    
    mm=ctStanModelMatrices(m)
    Js=J #replace parameter labels with matrix references
    for(x in 1:nrow(mm$matsetup)){
      Js=gsub(
        pattern = paste0("\\b",mm$matsetup$parname[x],"\\b"),
        replacement = paste0('s',matnames[mm$matsetup$matrix[x]],'[',
          mm$matsetup$row[x],',',
          mm$matsetup$col[x],']'),
        x = Js)
      
      for(mati in matnames) { #append s to system matrices
        Js=gsub(pattern = paste0("\\b",mati,"\\b"), replacement = paste0('s',mati),  x = Js)
      }
    }
    Jout[[typei]] <- matrix(Js,Jrows,ndim)
  }#end type loop
  # print(Jout)
  return(Jout)
}

