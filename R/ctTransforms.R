fQinf <- function(A,G){
  Ahatch=A %x% diag(1,nrow(A)) + 
    diag(1,nrow(A)) %x% A
  Qinf<-matrix(-solve(Ahatch , c(G %*% t(G))), nrow=nrow(A))
  try(dimnames(Qinf)<-dimnames(G))
  return(Qinf)
}

fdtQ <- function(Qinf, dtA) Qinf - (dtA %*% Qinf %*% t(dtA ))

fdtQe <- function(A, Q,dt=1){
  d=nrow(A)
  bA <- matrix(0,d*2,d*2)
  bA[1:d,1:d] <- -(A)
  bA[1:d,(d+1):(d*2)] <- Q
  bA[(d+1):(d*2),(d+1):(d*2)] <- t(A)
  
  ebA <- expm(bA %x% dt)
  dtQ<-t(ebA[(d+1):(d*2),(d+1):(d*2)]) %*% ebA[1:d,(d+1):(d*2)]
  return(dtQ)
}

fAstd <- function(A, G){
  d=nrow(A)
  asymDIFFUSION<-fQinf(A,G)
  standardiser <- rep(sqrt(diag(asymDIFFUSION)),each=d) / rep(sqrt(diag(asymDIFFUSION)),times=d)
  Astd<-A * standardiser
  try(dimnames(Astd)<-dimnames(A))
  return(Astd)
}

fdtAstd <- function(A, G, times){
  d=nrow(A)
  Qinf<-fQinf(A,G)
  standardiser <- rep(sqrt(diag(Qinf)),each=d) / rep(sqrt(diag(Qinf)),times=d)
  dtAstd<-lapply(times, function(x) expm::expm(A*x) * standardiser)
  return(dtAstd)
}

fdtA <- function(A, times){
  d=nrow(A)
  Astd<-lapply(times, function(x) expm::expm(A*x))
  return(Astd)
}


fAstd2 <- function(A, G,Jstep=1e-3){
  d=nrow(A)
  J <- diag(Inf,d)
  Qinf <- fQinf(A,G)
  for(i in 1:d){
    for(j in 1:d){
      # if(i!=j){
      As <- A
      As[i,j] <- As[i,j] + Jstep*sign(As[i,j])
      Qinfs <- fQinf(As,G)
      # J[i,j] <- (sum(diag(t(chol(Qinfs))))-sum(diag(t(chol(Qinf)))))/Jstep
      J[i,j] <- (sum(diag(t(chol(Qinfs))))-sum(diag(t(chol(Qinf)))))/Jstep
      # }
    }
  }
  return(J)
}


# A <- matrix(c(-1,.1,0,0,-1,.1,.1,0,-1),3,3)
# G <- matrix(c(1,0,0, 0,2,0, 0,0,.1),3,3)
# fQinf(A,G)
# # fdtAstd(A,G,1)
# A
# G
# fAstd(A,G)
# fAstd2(A,G)
# 
# 
# m <- matrix(c(2,0,0, 0,2,2, 0,2,2),3,3)
# sum(diag(chol(m)))

