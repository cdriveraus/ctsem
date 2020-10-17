functions{
  
    matrix constraincorsqrt(vector rawcor, int r){ //converts from unconstrained lower tri vec to cor sqrt
    int counter = 0;
    vector[r] ones = rep_vector(1, r);
    matrix[r, r] o;
    
    for(i in 2:r){ //constrain and set upper tri to lower
      for(j in 1:i - 1){
        counter += 1;
        o[i, j] = inv_logit(rawcor[counter]) * 2 - 1;  // can change cor prior here
        o[j, i] = o[i, j];
      }
    }
    for (i in 1:r){
      o[i, i] = 1;
      o[i, ] /= sqrt(sum(square(o[i, ])) + 1e-10);
    }

    return o;
  }
}
data{
  int d;
  int n;
  matrix[n,d] dat;
  int obs[n,d];
  int nobs[n];
  int reg;
}
parameters{
  vector[d] mu;
  vector[d] logsd;
  vector[(d * d - d) / 2] rawcor;
}
transformed parameters{
  matrix[d, d] mcor = tcrossprod(constraincorsqrt(rawcor,d));
  matrix[d,d] covm = diag_matrix(exp(logsd)+1e-10) * mcor * diag_matrix(exp(logsd)+1e-10);
  //matrix[d,d] cholm = cholesky_decompose(covm);
  real corprior = normal_lpdf(rawcor| 0, 1);
  real sdprior = normal_lpdf(logsd | 0, 10);
  vector[n] llrow=rep_vector(0,n);
    for(i in 1:n){
    if(nobs[i]>0){
      llrow[i]= multi_normal_lpdf(dat[i,obs[i,1:nobs[i]]] | mu[obs[i,1:nobs[i]]], 
      covm[obs[i,1:nobs[i]], obs[i,1:nobs[i]]]);
    }
  }
}
model{
  target += sum(llrow);
  if(reg)  target+= corprior + sdprior;
}
  
  
  
  
