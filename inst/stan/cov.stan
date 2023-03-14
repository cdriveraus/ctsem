functions{

  matrix constraincorsqrt(vector rawcor, int d){ //converts from unconstrained lower tri vec to cor sqrt
  int counter = 0;
  matrix[d,d] o;
  vector[d] ss = rep_vector(0,d);
  vector[d] s = rep_vector(0,d);
  real r;
  real r3;
  real r4;
  real r1;
  real r2;

  for(i in 1:d){ //set upper tri to lower
  for(j in 1:d){
    if(j > i){
      counter+=1;
      o[j,i] =  rawcor[counter];//inv_logit(rawcor[counter])*2-1; //divide by i for approx whole matrix equiv priors
    }
  }
  }

  for(i in 1:d){
    for(j in 1:d){
      if(j > i) {
        ss[i] +=square(o[j,i]);
        s[i] +=o[j,i];
      }
      if(j < i){
        ss[i] += square(o[i,j]);
        s[i] += o[i,j];
      }
    }
    s[i]+=1e-5;
    ss[i]+=1e-5;
  }


  for(i in 1:d){
    o[i,i]=0;
    r1=sqrt(ss[i]);
    r2=s[i];

    r3=(fabs(r2))/(r1)-1;
    r4=sqrt(log1p_exp(2*(fabs(r2)-r2-1)-4));
    r=(r4*((r3))+1)*r4+1;
    r=(sqrt(ss[i]+r));
    for(j in 1:d){
      if(j > i)  o[i,j]=o[j,i]/r;
      if(j < i) o[i,j] = o[i,j] /r;
    }
    o[i,i]=sqrt(1-sum(square(o[i,]))+1e-5);
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
  real reg;
  int corpriortype;
  int indep;
}
parameters{
  vector[d] mu;
  vector[d] logsd;
  vector[indep ? 0 : (d * d - d) / 2] rawcor;
}
transformed parameters{
  matrix[d, d] mcor = diag_matrix(rep_vector(1,d));
  matrix[d,d] covm;
  //matrix[d,d] cholm = cholesky_decompose(covm);
  real corprior=0;
  real sdprior = normal_lpdf(logsd | mean(logsd), 10);
  vector[n] llrow=rep_vector(0,n);

  if(!indep){
    mcor=tcrossprod(constraincorsqrt(rawcor,d));
    if(corpriortype==1)  corprior=normal_lpdf(rawcor| 0, 1); //mean(fabs(rawcor))
    if(corpriortype==2) corprior= normal_lpdf(to_vector(mcor) | 0, 1);
    if(corpriortype==3) corprior= normal_lpdf(eigenvalues_sym(mcor) | 0, 1);
  }

  covm = diag_matrix(exp(logsd)+1e-5) * mcor * diag_matrix(exp(logsd)+1e-5);

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
