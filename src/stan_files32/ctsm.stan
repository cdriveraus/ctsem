
functions{

   matrix constraincorsqrt(matrix mat){ //converts from unconstrained lower tri matrix to cor
    matrix[rows(mat),cols(mat)] o;
    vector[rows(mat)] s;
  
    for(i in 1:rows(o)){ //set upper tri to lower
      for(j in min(i+1,rows(mat)):rows(mat)){
        o[j,i] =  inv_logit(mat[j,i])*2-1;  // can change cor prior here
        o[i,j] = o[j,i];
      }
      o[i,i]=1; // change to adjust prior for correlations
    }

    for(i in 1:rows(o)){
      s[i] = inv_sqrt(o[i,] * o[,i]);
      if(is_inf(s[i])) s[i]=0;
    }
    return diag_pre_multiply(s,o);
  } 

  matrix sdcovsqrt2cov(matrix mat, int cholbasis){ //covariance from cholesky or unconstrained cor sq root
    if(cholbasis==0)  {
      return(tcrossprod(diag_pre_multiply(diagonal(mat),constraincorsqrt(mat))));
    } else return(tcrossprod(mat));
  }

  matrix sqkron_prod(matrix mata, matrix matb){
    int d=rows(mata);
    matrix[rows(mata)*rows(matb),cols(mata)*cols(matb)] out;
    for (k in 1:d){
      for (l in 1:d){
        for (i in 1:d){
          for (j in 1:d){
            out[ d*(i-1)+k, d*(j-1)+l ] = mata[i, j] * matb[k, l];
          }
        }
      }
    }
    return out;
  }

  matrix kronsum(matrix mata){
    matrix[rows(mata),rows(mata)] II = diag_matrix(rep_vector(1,rows(mata)));
    return sqkron_prod(mata, II) + sqkron_prod(II, mata );
  }

  vector colMeans(matrix mat){
    vector[cols(mat)] out;
    for(i in 1:cols(mat)){
      out[i] = mean(mat[,i]);
    }
    return out;
  }

  matrix cov_of_matrix(matrix mat){
    vector[cols(mat)] means = colMeans(mat);
    matrix[rows(mat), cols(mat)] centered;
    matrix[cols(mat), cols(mat)] covm;
    for (coli in 1:cols(mat)){
      for (ri in 1:rows(mat)){
        centered[ri,coli] = mat[ri,coli] - means[coli];
      }
    }
    covm = crossprod(centered) / (rows(mat)-1);
    return covm; 
  }

  matrix crosscov(matrix a, matrix b){
    matrix[rows(a),cols(a)] da;
    matrix[rows(b),cols(b)] db;
    matrix[cols(a),cols(b)] out;
  
    da = a - rep_matrix( (colMeans(a))',rows(a));
    db = b - rep_matrix( (colMeans(b))',rows(b));
    out = da' * db ./ (rows(a)-1.0);
    return out;
  }


  matrix makesym(matrix mat, int verbose, int pd){
    matrix[rows(mat),cols(mat)] out;
    for(coli in 1:cols(mat)){
      if(pd ==1){ // && mat[coli,coli] < 1e-5
        //if(verbose > 0) print("diagonal too low (",mat[coli,coli],") during makesym row ", coli, " col ", coli);
        out[coli,coli] = mat[coli,coli] + 1e-5;
      } else out[coli,coli] = mat[coli,coli]; 
      for(rowi in coli:rows(mat)){
        if(rowi > coli) {
          out[rowi,coli] = mat[rowi,coli]; //(mat[coli,rowi] + ) *.5;
          out[coli,rowi] = mat[rowi,coli];
        }
        
      }
    }
    return out;
  }

  real tform(real param, int transform, data real multiplier, data real meanscale, data real offset, data real inneroffset){
    real out;
  
      
  if(transform==0) out = param * meanscale * multiplier +inneroffset + offset; 

  if(transform==1) out = log(1+(exp(param * meanscale+inneroffset))) * multiplier + offset ; 

  if(transform==2) out = exp(param * meanscale+inneroffset) * multiplier + offset; 

  if(transform==3) out = inv_logit(param*meanscale+inneroffset) * multiplier + offset; 

  if(transform==4) out = ((param*meanscale+inneroffset)^3)*multiplier + offset; 



    return out;
  }

}
data {
  int<lower=0> ndatapoints;
  int<lower=1> nmanifest;
  int<lower=1> nlatent;
  int<lower=1> nsubjects;
  int<lower=0> ntipred; // number of time independent covariates
  int<lower=0> ntdpred; // number of time dependent covariates

  int T0check[ndatapoints]; // logical indicating which rows are the first for each subject
  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipredsdata;
  int nmissingtipreds;
  int ntipredeffects;
  real<lower=0> tipredsimputedscale;
  real<lower=0> tipredeffectscale;

  vector[nmanifest] Y[ndatapoints];
  int nopriors;
  int nldynamics;
  vector[ntdpred] tdpreds[ntdpred ? ndatapoints : 0];
  
  real dT[ndatapoints]; // time intervals
  real dTsmall[ndatapoints];
  int integrationsteps[ndatapoints] ; // time steps needed between time intervals for integration
  int subject[ndatapoints];
  int<lower=0> nparams;
  int continuoustime; // logical indicating whether to incorporate timing information
  int nindvarying; // number of subject level parameters that are varying across subjects
  int nindvaryingoffdiagonals; //number of off diagonal parameters needed for popcov matrix
  vector[nindvarying] sdscale;
  int indvaryingindex[nindvarying];
  int notindvaryingindex[nparams-nindvarying];

  int nt0varstationary;
  int nt0meansstationary;
  int t0varstationary [nt0varstationary, 2];
  int t0meansstationary [nt0meansstationary, 2];

  int nobs_y[ndatapoints];  // number of observed variables per observation
  int whichobs_y[ndatapoints, nmanifest]; // index of which variables are observed per observation
  int ndiffusion; //number of latents involved in covariance calcs
  int derrind[ndiffusion]; //index of which latent variables are involved in covariance calculations

  int manifesttype[nmanifest];
  int nbinary_y[ndatapoints];  // number of observed binary variables per observation
  int whichbinary_y[ndatapoints, nmanifest]; // index of which variables are observed and binary per observation
  int ncont_y[ndatapoints];  // number of observed continuous variables per observation
  int whichcont_y[ndatapoints, nmanifest]; // index of which variables are observed and continuous per observation
  
  int intoverpop;
  int ukf;
  real ukfspread;
  int ukffull;
  int nlmeasurement;
  int intoverstates;
  int verbose; //level of printing during model fit

  int T0MEANSsubindex[nsubjects];
int LAMBDAsubindex[nsubjects];
int DRIFTsubindex[nsubjects];
int DIFFUSIONsubindex[nsubjects];
int MANIFESTVARsubindex[nsubjects];
int MANIFESTMEANSsubindex[nsubjects];
int CINTsubindex[nsubjects];
int T0VARsubindex[nsubjects];
int TDPREDEFFECTsubindex[nsubjects];
int asymCINTsubindex[nsubjects];
int asymDIFFUSIONsubindex[nsubjects];
  int T0MEANSsetup_rowcount;
int LAMBDAsetup_rowcount;
int DRIFTsetup_rowcount;
int DIFFUSIONsetup_rowcount;
int MANIFESTVARsetup_rowcount;
int MANIFESTMEANSsetup_rowcount;
int CINTsetup_rowcount;
int T0VARsetup_rowcount;
int TDPREDEFFECTsetup_rowcount;
  int T0MEANSsetup[T0MEANSsetup_rowcount,6 ];
int LAMBDAsetup[LAMBDAsetup_rowcount,6 ];
int DRIFTsetup[DRIFTsetup_rowcount,6 ];
int DIFFUSIONsetup[DIFFUSIONsetup_rowcount,6 ];
int MANIFESTVARsetup[MANIFESTVARsetup_rowcount,6 ];
int MANIFESTMEANSsetup[MANIFESTMEANSsetup_rowcount,6 ];
int CINTsetup[CINTsetup_rowcount,6 ];
int T0VARsetup[T0VARsetup_rowcount,6 ];
int TDPREDEFFECTsetup[TDPREDEFFECTsetup_rowcount,6 ];
  matrix[T0MEANSsetup_rowcount, 6] T0MEANSvalues;
matrix[LAMBDAsetup_rowcount, 6] LAMBDAvalues;
matrix[DRIFTsetup_rowcount, 6] DRIFTvalues;
matrix[DIFFUSIONsetup_rowcount, 6] DIFFUSIONvalues;
matrix[MANIFESTVARsetup_rowcount, 6] MANIFESTVARvalues;
matrix[MANIFESTMEANSsetup_rowcount, 6] MANIFESTMEANSvalues;
matrix[CINTsetup_rowcount, 6] CINTvalues;
matrix[T0VARsetup_rowcount, 6] T0VARvalues;
matrix[TDPREDEFFECTsetup_rowcount, 6] TDPREDEFFECTvalues;
  int TIPREDEFFECTsetup[nparams, ntipred];
  int nmatrixslots;
  int popsetup[nmatrixslots,6];
  real popvalues[nmatrixslots,6];
  int savescores;
  int nlatentpop;
  int fixedsubpars;
  vector[fixedsubpars ? nindvarying : 0] fixedindparams[fixedsubpars ? nsubjects : 0];
  int dokalman;
}
      
transformed data{
  matrix[nlatent,nlatent] IIlatent= diag_matrix(rep_vector(1,nlatent));
  matrix[nlatent*nlatent,nlatent*nlatent] IIlatent2 = diag_matrix(rep_vector(1,nlatent*nlatent));
  matrix[nindvarying,nindvarying] IIindvar = diag_matrix(rep_vector(1,nindvarying));
  real asquared =  square(2.0/sqrt(0.0+nlatentpop) * ukfspread);
  real sqrtukfadjust = sqrt(0.0+nlatentpop +( asquared * (nlatentpop  + 0.5) - (nlatentpop) ) );
}
      
parameters {
  vector[nparams] rawpopmeans; // population level means 

  vector[nindvarying] rawpopsdbase; //population level std dev
  vector[nindvaryingoffdiagonals] sqrtpcov; // unconstrained basis of correlation parameters
  vector[fixedsubpars ? 0 : (intoverpop ? 0 : nindvarying)] baseindparams[fixedsubpars ? 0 : (intoverpop ? 0 : nsubjects)]; //vector of subject level deviations, on the raw scale
  
  vector[ntipredeffects] tipredeffectparams; // effects of time independent covariates
  vector[nmissingtipreds] tipredsimputed;
  
  vector[intoverstates ? 0 : nlatent*ndatapoints] etaupdbasestates; //sampled latent states posterior
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying,nindvarying] rawpopcovsqrt; 
  real ll;
  vector[nmanifest+nmanifest+ (savescores ? nmanifest*2+nlatentpop*2 : 0)] kalman[savescores ? ndatapoints : 0];
  matrix[ T0MEANSsetup_rowcount ? max(T0MEANSsetup[,1]) : 0, T0MEANSsetup_rowcount ? max(T0MEANSsetup[,2]) : 0 ] T0MEANS[T0MEANSsubindex[nsubjects]];
matrix[ LAMBDAsetup_rowcount ? max(LAMBDAsetup[,1]) : 0, LAMBDAsetup_rowcount ? max(LAMBDAsetup[,2]) : 0 ] LAMBDA[LAMBDAsubindex[nsubjects]];
matrix[ DRIFTsetup_rowcount ? max(DRIFTsetup[,1]) : 0, DRIFTsetup_rowcount ? max(DRIFTsetup[,2]) : 0 ] DRIFT[DRIFTsubindex[nsubjects]];
matrix[ DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,1]) : 0, DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,2]) : 0 ] DIFFUSION[DIFFUSIONsubindex[nsubjects]];
matrix[ MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,1]) : 0, MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,2]) : 0 ] MANIFESTVAR[MANIFESTVARsubindex[nsubjects]];
matrix[ MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,1]) : 0, MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,2]) : 0 ] MANIFESTMEANS[MANIFESTMEANSsubindex[nsubjects]];
matrix[ CINTsetup_rowcount ? max(CINTsetup[,1]) : 0, CINTsetup_rowcount ? max(CINTsetup[,2]) : 0 ] CINT[CINTsubindex[nsubjects]];
matrix[ T0VARsetup_rowcount ? max(T0VARsetup[,1]) : 0, T0VARsetup_rowcount ? max(T0VARsetup[,2]) : 0 ] T0VAR[T0VARsubindex[nsubjects]];
matrix[ TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,1]) : 0, TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,2]) : 0 ] TDPREDEFFECT[TDPREDEFFECTsubindex[nsubjects]];

  matrix[nlatent,nlatent] asymDIFFUSION[asymDIFFUSIONsubindex[nsubjects]]; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] asymCINT[asymCINTsubindex[nsubjects]]; // latent process asymptotic level

  
  matrix[ T0MEANSsetup_rowcount ? max(T0MEANSsetup[,1]) : 0, T0MEANSsetup_rowcount ? max(T0MEANSsetup[,2]) : 0 ] pop_T0MEANS;
matrix[ LAMBDAsetup_rowcount ? max(LAMBDAsetup[,1]) : 0, LAMBDAsetup_rowcount ? max(LAMBDAsetup[,2]) : 0 ] pop_LAMBDA;
matrix[ DRIFTsetup_rowcount ? max(DRIFTsetup[,1]) : 0, DRIFTsetup_rowcount ? max(DRIFTsetup[,2]) : 0 ] pop_DRIFT;
matrix[ DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,1]) : 0, DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,2]) : 0 ] pop_DIFFUSION;
matrix[ MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,1]) : 0, MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,2]) : 0 ] pop_MANIFESTVAR;
matrix[ MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,1]) : 0, MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,2]) : 0 ] pop_MANIFESTMEANS;
matrix[ CINTsetup_rowcount ? max(CINTsetup[,1]) : 0, CINTsetup_rowcount ? max(CINTsetup[,2]) : 0 ] pop_CINT;
matrix[ T0VARsetup_rowcount ? max(T0VARsetup[,1]) : 0, T0VARsetup_rowcount ? max(T0VARsetup[,2]) : 0 ] pop_T0VAR;
matrix[ TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,1]) : 0, TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,2]) : 0 ] pop_TDPREDEFFECT;

  matrix[nlatent,nlatent] pop_asymDIFFUSION; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] pop_asymCINT; // latent process asymptotic level

  

  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipreds; //tipred values to fill from data and, when needed, imputation vector
  matrix[nparams, ntipred] TIPREDEFFECT; //design matrix of individual time independent predictor effects

  if(ntipred > 0){ 
    int counter = 0;
    for(coli in 1:cols(tipreds)){ //insert missing ti predictors
      for(rowi in 1:rows(tipreds)){
        if(tipredsdata[rowi,coli]==99999) {
          counter += 1;
          tipreds[rowi,coli] = tipredsimputed[counter];
        } else tipreds[rowi,coli] = tipredsdata[rowi,coli];
      }
    }
    for(ci in 1:ntipred){ //configure design matrix
      for(ri in 1:nparams){
        if(TIPREDEFFECTsetup[ri,ci] > 0) {
          TIPREDEFFECT[ri,ci] = tipredeffectparams[TIPREDEFFECTsetup[ri,ci]];
        } else {
          TIPREDEFFECT[ri,ci] = 0;
        }
      }
    }
  }

  if(nindvarying > 0){
    int counter =0;
    rawpopsd = exp(2*rawpopsdbase-1) .* sdscale; // sqrts of proportions of total variance
    for(j in 1:nindvarying){
      rawpopcovsqrt[j,j] = 1;
      for(i in 1:nindvarying){
        if(i > j){
          counter += 1;
          rawpopcovsqrt[i,j]=sqrtpcov[counter];
          rawpopcovsqrt[j,i]=sqrtpcov[counter];
        }
      }
    }
 rawpopcovsqrt = cholesky_decompose(makesym(tcrossprod(diag_pre_multiply(rawpopsd, 
      constraincorsqrt(rawpopcovsqrt))),verbose,1)); 
  }//end indvarying par setup

  ll=0;
if(dokalman==1){
}

}
      
model{

  if(intoverpop==0 && fixedsubpars == 1) fixedindparams ~ multi_normal_cholesky(rep_vector(0,nindvarying),IIindvar);

  if(nopriors==0){
   target += normal_lpdf(rawpopmeans|0,1);
  
    if(ntipred > 0){ 
      tipredeffectparams ~ normal(0,tipredeffectscale);
      tipredsimputed ~ normal(0,tipredsimputedscale);
    }
    
    if(nindvarying > 0){
      if(nindvarying >1) sqrtpcov ~ normal(0,1);
      if(intoverpop==0 && fixedsubpars == 0) baseindparams ~ multi_normal_cholesky(rep_vector(0,nindvarying),IIindvar);
      rawpopsdbase ~ normal(0,1);
    }
  } //end pop priors section
  
  if(intoverstates==0) etaupdbasestates ~ normal(0,1);
  
  target += ll;
  if(verbose > 0) print("lp = ", target());
}
generated quantities{
  vector[nparams] popmeans;
  vector[nparams] popsd = rep_vector(0,nparams);
  matrix[nindvarying,nindvarying] rawpopcov = tcrossprod(rawpopcovsqrt);
  matrix[nindvarying,nindvarying] rawpopcorr = quad_form_diag(rawpopcov,inv_sqrt(diagonal(rawpopcov)));
  matrix[nparams,ntipred] linearTIPREDEFFECT;
  

{
vector[nparams] rawpopsdfull;
rawpopsdfull[indvaryingindex] = sqrt(diagonal(rawpopcov)); //base for calculations

    for(ri in 1:dims(popsetup)[1]){
      if(popsetup[ri,3] !=0) {

        popmeans[popsetup[ ri,3]] = tform(rawpopmeans[popsetup[ri,3] ], popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4], popvalues[ri,6] ); 

        popsd[popsetup[ ri,3]] = popsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[popsetup[ri,3] ]  + rawpopsdfull[popsetup[ ri,3]], popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4], popvalues[ri,6]) -
           tform(
            rawpopmeans[popsetup[ri,3] ]  - rawpopsdfull[popsetup[ ri,3]], popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4], popvalues[ri,6]) ) /2 : 
          0; 

        if(ntipred > 0){
          for(tij in 1:ntipred){
            if(TIPREDEFFECTsetup[popsetup[ri,3],tij] ==0) {
              linearTIPREDEFFECT[popsetup[ri,3],tij] = 0;
            } else {
            linearTIPREDEFFECT[popsetup[ri,3],tij] = (
              tform(rawpopmeans[popsetup[ri,3] ] + TIPREDEFFECT[popsetup[ri,3],tij] * .01, popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4], popvalues[ri,6] ) -
              tform(rawpopmeans[popsetup[ri,3] ] - TIPREDEFFECT[popsetup[ri,3],tij] * .01, popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4], popvalues[ri,6] )
              ) /2 * 100;
            }
         }
        }
      }
    }
}



}
