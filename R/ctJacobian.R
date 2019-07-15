# ctm=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','state[3]',0, 'd12','d22',0,0,0,0),3,3),type='stanct')
# ctm2=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','PARS[1,1]',0, 'd12','d22',0,0,0,0),3,3),PARS=matrix('state[3]'),type='stanct')
#'@import Deriv Simplify prodSymb
#'@import cOde jacobianSymb

ctJacobian <- function(m,types=c('JAx','Jtd') ){
  
  require(cOde)
  require(Deriv)
  
  
  # get system dimension
  ndim = max(m$pars$row[m$pars$matrix%in% 'T0MEANS'])
  
  # 2): generate vector valued function fn = drift * state
  
  # initialize fn and state
  fn     = c()
  state   = paste0("state__", 1:ndim,'__')
  
  #replace system matrix references
  matlist <- listOfMatrices(m$pars)
  matnames <- ctStanMatricesList(m,unsafe=TRUE)$base
  
  # for(mati in matnames){ 
  for(ri in 1:nrow(m$pars)){
    counter <- 0
    while(counter < 10 && grepl(paste0('\\b(',paste0(names(matlist),collapse='|'),')\\b\\['), m$pars$param[ri])){ #if system matrix referenced
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
  # }
  
  
  mats<-listOfMatrices(m$pars)
  Jout <- list()
  for(typei in types){
    
    if(typei=='JAx'){
      # create fn by multiplying drift matrix with state vector
      # fn = Simplify(prodSymb(mats$DRIFT, cbind(state)))
      for (row in 1:ndim) {
        for (col in 1:(ndim)) {
          fn[row] = paste0(ifelse(col > 1, paste0(fn[row],' + '),''), "(", mats$DRIFT[row, col], ") * state[", as.character(col), "]")
        }
        if(!is.na(mats$CINT[row])) fn[row] = paste0(fn[row],' + ',mats$CINT[row]) #checking for NA because CINT is not always as large as DRIFT
      }
      fn=sapply(fn,Simplify)
    }

    # if(typei=='JJAx'){
    #   fn=Jout$JAx
    # }
    
    if(typei=='Jtd'){
      if(m$n.TDpred ==0) {
        Jout[[typei]] <- ''
        next
      }
      fn = Simplify(prodSymb(mats$TDPREDEFFECT, cbind(rep(1,m$n.TDpred))))
      # # create fn by multiplying drift matrix with state vector
      # for (row in 1:nrow(mats$TDPREDEFFECT)) {
      #   for (col in 1:(ncol(mats$TDPREDEFFECT))) {
      #     fn[row] = paste0(ifelse(col > 1, paste0(fn[row],' + '),''), "(", mats$TDPREDEFFECT[row, col], ")")
      #   }
      # }
    }
    
    fn = gsub(" ", "", fn, fixed = TRUE) #remove spaces
    # replace state[~] by state~ for cOde Jacobian and make fn and state a named list
    names(fn) = paste0("fn", 1:length(fn))
    for(statei in 1:ndim){
    fn=gsub(paste0('\\b(state)\\[(',statei,')?\\]'),paste0('state__',statei,'__'),fn,perl = TRUE)
    }

    
    #probably redundant now but maybe useful at some point?
    # replace remaining commas and square brackets for cOde Jacobian
    fn = gsub(",", "___comma___", fn, fixed = TRUE)
    fn = gsub("[", "___leftsquarebracket___", fn, fixed = TRUE)
    fn = gsub("]", "___rightsquarebracket___", fn, fixed = TRUE)
    # 
    
    
    # 3): calculate Jacobian of fn symbolically
    J  = jacobianSymb(fn, state)
    
    
    # 4): create Jacobian list in STAN format
    J_STAN = sapply(J,Simplify)
    J_STAN = gsub(" ", "", J_STAN, fixed = TRUE) #remove spaces
    
    for(statei in 1:ndim){
      J_STAN=gsub(paste0('state__',statei,'__'),paste0('state[',statei,']'),J_STAN,fixed = TRUE)
    }

    # restore commas and square brackets for Jacobian list
    J_STAN = gsub("___rightsquarebracket___", "]", J_STAN, fixed = TRUE)
    J_STAN = gsub("___leftsquarebracket___", "[", J_STAN, fixed = TRUE)
    J_STAN = gsub("___comma___", ",", J_STAN, fixed = TRUE)
    
    
    mm=ctStanModelMatrices(m)
    
    Js=J_STAN #replace parameter labels with matrix references
    for(x in 1:nrow(mm$popsetup)){
      Js=gsub(
        pattern = paste0("\\b",mm$popsetup$parname[x],"\\b"),
        replacement = paste0('s',matnames[mm$popsetup$matrix[x]],'[',
          mm$popsetup$row[x],',',
          mm$popsetup$col[x],']'),
        x = Js)
      
      for(mati in matnames) { #append s to system matrices
        Js=gsub(pattern = paste0("\\b",mati,"\\b"), replacement = paste0('s',mati),  x = Js)
      }
    }
    Jout[[typei]] <- matrix(Js,ndim,ndim)
  }#end type loop
  print(Jout$JAx)
  return(Jout)
}

