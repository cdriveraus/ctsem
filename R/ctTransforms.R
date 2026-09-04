fQinf <- function(A,G){
  Ahatch=A %x% diag(1,nrow(A)) +
    diag(1,nrow(A)) %x% A
  Qinf<-matrix(-solve(Ahatch , c(G %*% t(G))), nrow=nrow(A))
  try(dimnames(Qinf)<-dimnames(G))
  return(Qinf)
}

fdtQ <- function(Qinf, dtA) Qinf - (dtA %*% Qinf %*% t(dtA ))
