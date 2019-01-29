data {
  int<lower=0> ndatapoints;
}
parameters {
 real sigma;
}
model{
sigma ~ normal(0,1);
}

