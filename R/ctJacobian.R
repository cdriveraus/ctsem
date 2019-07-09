# ctm=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','state[3]',0, 'd12','d22',0,0,0,0),3,3),type='stanct')
# ctm2=ctModel(LAMBDA=cbind(diag(2),0), DRIFT=matrix(c('d11','PARS[1,1]',0, 'd12','d22',0,0,0,0),3,3),PARS=matrix('state[3]'),type='stanct')


ctJacobian <- function(m){

require(cOde)

# get system dimension
n = m$n.latent

# declare drift matrix
drift = matrix(NA, n, n)

# fill drift matrix entries
for(row in 1:nrow(m$pars)){
  
  if(m$pars$matrix[row] %in% "DRIFT") {
    
    if (is.na(m$pars$value[row])) {
      
      drift[m$pars$row[row], m$pars$col[row]] = m$pars$param[row]
      
    }
    
    else {
      
      drift[m$pars$row[row], m$pars$col[row]] = m$pars$value[row]
      
    }
    
  }
  
}



# 2): generate vector valued function f = drift * state

# initialize f and state
f     = numeric(n)
state = numeric(n)

#replace system matrix references
matlist <- listOfMatrices(m$pars)
for(mati in ctStanMatricesList(m,unsafe=TRUE)$base){
  for(ri in 1:nrow(m)){
    if(grepl(paste0(mati,'['), m$pars$param[ri],fixed=TRUE)){ #if system matrix referenced

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
}
    

# create f by multiplying drift matrix with state vector
for (row in 1:n) {
  
  for (col in 1:(n-1)) {
    
    f[row] = paste0(f[row], "(", drift[row, col], ") * state[", as.character(col), "]", " + ")
    
  }
  
  # substring removes "0" from initiation
  f[row] = substring(paste0(f[row], "(", drift[row, n], ") * state[", as.character(n), "]"), 2)
  
}



# replace state[~] by state~ for cOde Jacobian and make f and state a named list
for (i in 1:n) {
  
  state[i]    = paste0("state", as.character(i))
  
  names(f)[i] = paste0("f", as.character(i))
  
  f           = gsub(paste0("state[", as.character(i), "]"),
                     paste0("state", as.character(i)),
                     f, fixed = TRUE)
  
}

# replace remaining commas and square brackets for cOde Jacobian
f = gsub(",", "___comma___", f, fixed = TRUE)
f = gsub("[", "___leftsquarebracket___", f, fixed = TRUE)
f = gsub("]", "___rightsquarebracket___", f, fixed = TRUE)



# 3): calculate Jacobian of f symbolically
J  = jacobianSymb(f, state)
df = matrix(J, nrow = n, ncol = n)



# 4): create Jacobian list in STAN format
J_STAN = J

# replace state~ by state[~] for Jacobian list
for (i in 1:n) {
  
  J_STAN = gsub(paste0("state", as.character(i)), 
                paste0("state[", as.character(i), "]"), 
                J_STAN, fixed = TRUE)
  
}

# restore commas and square brackets for Jacobian list
J_STAN = gsub("___rightsquarebracket___", "]", J_STAN, fixed = TRUE)
J_STAN = gsub("___leftsquarebracket___", "[", J_STAN, fixed = TRUE)
J_STAN = gsub("___comma___", ",", J_STAN, fixed = TRUE)

return(matrix(J_STAN,n,n))
}
