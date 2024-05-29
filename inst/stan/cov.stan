functions{
  
  matrix constraincorsqrt(matrix rawcor) {
    int d = rows(rawcor);
    matrix[d, d] mcholcor = rep_matrix(0, d, d);
    mcholcor[1, 1] = 1;
    if (d > 1) {
      for (coli in 1:d) {
        for (rowi in coli:d) {
          if (coli == 1 && rowi > 1) mcholcor[rowi, coli] = rawcor[rowi, coli];
          if (coli > 1) {
            if (rowi == coli) mcholcor[rowi, coli] = prod(sqrt(1 - rawcor[rowi, 1:(coli - 1)]^2));
            if (rowi > coli) mcholcor[rowi, coli] = rawcor[rowi, coli] * prod(sqrt(1 - rawcor[rowi, 1:(coli - 1)]^2));
          }
        }
      }
    }
    return mcholcor;
  }
  
  
}
data{
  int d;
  int n;
  matrix[n,d] dat;
  array[n,d] int obs;
  array[n] int nobs;
  real reg;
  int corpriortype;
  int indep;
}
parameters{
  vector[d] mu;
  vector[d] logsd;
  vector[indep ? 0 : (d * d - d) %/% 2] rawcor;
}
transformed parameters{
  matrix[d, d] mcor = diag_matrix(rep_vector(1,d));
  matrix[d, d] mcorbase = rep_matrix(0, d, d);
  matrix[d,d] covm;
  matrix[d,d] cholm;
  //matrix[d,d] cholm = cholesky_decompose(covm);
  real corprior=0;
  real sdprior = normal_lpdf(logsd | mean(logsd), 10);
  vector[n] llrow=rep_vector(0,n);
  
  if(!indep){
    int counter=0;
    for(i in 1:d){
      for(j in 1:d){
        if(i > j){
          counter+=1;
          mcorbase[i,j]=inv_logit(rawcor[counter])*2-1;
        }
      }
    }
    //print(mcorbase);
    
    mcor=tcrossprod(constraincorsqrt(mcorbase));
    //print(mcor);
    if(reg != 0){
      if(corpriortype==1)  corprior=normal_lpdf(rawcor| 0, 1); //mean(abs(rawcor))
      if(corpriortype==2) corprior= normal_lpdf(to_vector(mcor) | 0, 1);
      if(corpriortype==3) corprior= normal_lpdf(to_vector(eigenvalues_sym(mcor)) | 0, 1);
    }
  }
  
  covm = quad_form_diag(mcor, log1p_exp(logsd));
  cholm = cholesky_decompose(covm);
  
  for(i in 1:n){
    if(nobs[i]>0){
      llrow[i]= multi_normal_lpdf(dat[i,obs[i,1:nobs[i]]] | mu[obs[i,1:nobs[i]]],
      covm[obs[i,1:nobs[i]], obs[i,1:nobs[i]]]);
    }
  }
}
model{
  target += sum(llrow);
  if(reg!=0)  target+= reg*corprior + reg*sdprior;
}
