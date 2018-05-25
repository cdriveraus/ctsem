
functions{

matrix covsqrt2corsqrt(matrix mat, int invert){ //converts from lower partial sd matrix to cor
      matrix[rows(mat),cols(mat)] o;
      vector[rows(mat)] s;
    o=mat;

    for(i in 1:rows(o)){ //set upper tri to lower
for(j in min(i+1,rows(mat)):rows(mat)){
o[j,i] = inv_logit(o[j,i])*2-1;  // can change cor prior here
o[i,j] = o[j,i];
}
      o[i,i]=1; // change to adjust prior for correlations
    }

if(invert==1) o = inverse(o);


  for(i in 1:rows(o)){
      s[i] = inv_sqrt(o[i,] * o[,i]);
    if(is_inf(s[i])) s[i]=0;
    }
      o= diag_pre_multiply(s,o);
return o;
 }

matrix cholspd(matrix a){
    matrix[rows(a),rows(a)] l;
    for(j in 1:cols(a)){
      for(i in j:rows(a)){
        if(i != j) {
          l[i,j] = a[i,j];
          l[j,i] = a[i,j];
        }
        if(j == i)  l[j,j] = a[j,j] + 1e-10;
         // if(a[j,j] <=1e-10){
         //   l[j,j] = 1e-10; 
        //  } else 
        //    l[j,j] = a[j,j]; //print("Negative variance of ", l[j,j], " at row ", j);
        //}
      }
    }
    return cholesky_decompose(l);
}

  matrix discreteDIFFUSIONcalc(matrix DR, matrix DI, real dt){
    matrix[rows(DR)+rows(DI),rows(DR)+rows(DI)]  DRDI;
    matrix[rows(DR),rows(DR)] out;
    int d;
    
    d=rows(DR);
    DRDI[1:d,1:d] = -DR;
    DRDI[1:d,(d+1):(d*2)] = DI;
    DRDI[(d+1):(d*2), (d+1):(d*2)] = DR';
    DRDI[(d+1):(d*2), 1:d] = rep_matrix(0,d,d);
    DRDI = matrix_exp(DRDI * dt);
    out = DRDI[(d+1):(d*2), (d+1):(d*2)]' * DRDI[1:d, (d+1):(d*2)];
    return out;
  }
        

  matrix matrix_diagexp(matrix in){
    matrix[rows(in),rows(in)] out;
    for(i in 1:rows(in)){
      for(j in 1:rows(in)){
        if(i==j) out[i,i] = exp(in[i,i]);
        if(i!=j) out[i,j] = 0;
      }
    }
  return out;
  }




  matrix sdcovsqrt2cov(matrix mat, int cholesky){ //converts from lower partial sd and diag sd to cov or cholesky cov
    matrix[rows(mat),rows(mat)] out;

    for(k in 1:cols(mat)){
      for(j in 1:rows(mat)){
        if(j > k) out[j,k] = mat[j,k];
        if(k > j) out[j,k] = mat[k,j];
        if(k==j) out[j,k] = mat[j,k];
      }
    }
    if(cholesky==0) out = tcrossprod(out);
    return(out);
  }
      

  matrix kron_prod(matrix mata, matrix matb){
    int m;
    int p;
    int n;
    int q;
    matrix[rows(mata)*rows(matb),cols(mata)*cols(matb)] C;
    m=rows(mata);
    p=rows(matb);
    n=cols(mata);
    q=cols(matb);
    for (k in 1:p){
      for (l in 1:q){
        for (i in 1:m){
          for (j in 1:n){
            C[p*(i-1)+k,q*(j-1)+l] = mata[i,j]*matb[k,l];
          }
        }
      }
    }
    return C;
  }

  matrix cov_of_matrix(matrix mat){
    vector[cols(mat)] means;
    matrix[rows(mat), cols(mat)] centered;
    matrix[cols(mat), cols(mat)] covm;
    for (coli in 1:cols(mat)){
      means[coli] = mean(mat[,coli]);
      for (rowi in 1:rows(mat))  {
        centered[rowi,coli] = mat[rowi,coli] - means[coli];
      }
    }
    covm = crossprod(centered) / (rows(mat)-1);
    for(j in 1:rows(covm)){
      covm[j,j] = covm[j,j] + 1e-8;
    }
    return covm; 
  }

  vector colMeans(matrix mat){
    vector[cols(mat)] out;
    for(i in 1:cols(mat)){
      out[i] = mean(mat[,i]);
    }
    return out;
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

  matrix chol(matrix a){
    matrix[rows(a),rows(a)] l;
    for(j in 1:cols(a)){
      for(i in 1:rows(a)){
        if(j==i) {
          if(j == 1){ 
            if(a[j,j] <=0) { 
              l[j,j] = 1e-8;
              print("Negative variance ", a[j,j], " set to 1e-8 for Cholesky decomp");
            }
            if(a[j,j] > 0) l[j,j] = sqrt(a[j,j]);
          }
          if(j > 1) l[j,j] = sqrt(a[j,j] - dot_self(l[i, 1 : j-1]));
        }
        if(i > j) l[i,j] = ( a[i,j] - dot_product( l[ i, 1:(j-1) ], l[j, 1:(j-1)]) ) / l[j,j];
        if(j > i) l[i,j] = 0;
      }
    }
    return l;
  }

  real tform(real param, int transform, real multiplier, real meanscale, real offset){
    real out;
  
      if(transform==0) out = param * meanscale * multiplier + offset; 

  if(transform==1) out = log(1+(exp(param * meanscale))) * multiplier + offset ; 

  if(transform==2) out = exp(param * meanscale) * multiplier + offset; 

  if(transform==3) out = inv_logit(param*meanscale) * multiplier + offset; 

  if(transform==4) out = ((param*meanscale)^3)*multiplier + offset; 



    return out;
  }

  matrix cov2cors(matrix M){
    matrix[rows(M),cols(M)] o;
    vector[rows(M)] isd;

    isd = inv_sqrt(diagonal(M));
    o = quad_form_diag(M,isd);
    return(o);
  }


}
data {
  int<lower=0> ndatapoints;
  int<lower=1> nmanifest;
  int<lower=1> nlatent;
  int<lower=1> nsubjects;
  int<lower=0> ntipred; // number of time independent covariates
  int<lower=0> ntdpred; // number of time dependent covariates

  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipredsdata;
  int nmissingtipreds;
  int ntipredeffects;
  
  vector[nmanifest] Y[ndatapoints];
  int nopriors;
  int lineardynamics;
  vector[ntdpred] tdpreds[ntdpred ? ndatapoints : 0];
  
  real dT[ndatapoints]; // time intervals
  real dTsmall[ndatapoints];
  int driftdiagonly; //can we simplify matrix exponential to univariate?
  int binomial; //binary data only
  int integrationsteps[ndatapoints] ; // time steps needed between time intervals for integration
  int driftindex[ndatapoints]; //which discreteDRIFT matrix to use for each time point
  int diffusionindex[ndatapoints]; //which discreteDIFFUSION matrix to use for each time point
  int cintindex[ndatapoints]; //which discreteCINT matrix to use for each time point
  int subject[ndatapoints];
  int<lower=0> nparams;
  int T0check[ndatapoints]; // logical indicating which rows are the first for each subject
  int continuoustime; // logical indicating whether to incorporate timing information
  int nindvarying; // number of subject level parameters that are varying across subjects
  int nindvaryingoffdiagonals; //number of off diagonal parameters needed for popcov matrix
  int notindvaryingindex[nparams-nindvarying];
  int indvaryingindex[nindvarying];
  vector[nindvarying] sdscale;

  int nt0varstationary;
  int nt0meansstationary;
  int t0varstationary [nt0varstationary, 2];
  int t0meansstationary [nt0meansstationary, 2];

  int<lower = 0, upper = nmanifest> nobs_y[ndatapoints];  // number of observed variables per observation
  int<lower = 0, upper = nmanifest> whichobs_y[ndatapoints, nmanifest]; // index of which variables are observed per observation
  int<lower=0,upper=nlatent> ndiffusion; //number of latents involved in covariance calcs
  int<lower=0,upper=nlatent> derrind[ndiffusion]; //index of which latent variables are involved in covariance calculations

  int manifesttype[nmanifest];
  int<lower = 0, upper = nmanifest> nbinary_y[ndatapoints];  // number of observed binary variables per observation
  int<lower = 0, upper = nmanifest> whichbinary_y[ndatapoints, nmanifest]; // index of which variables are observed and binary per observation
  int<lower = 0, upper = nmanifest> ncont_y[ndatapoints];  // number of observed continuous variables per observation
  int<lower = 0, upper = nmanifest> whichcont_y[ndatapoints, nmanifest]; // index of which variables are observed and continuousper observation
  
  int ukfpop;
  int ukf;
  int intoverstates;
  int ngenerations; //number of samples of random data to generate
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
  int T0MEANSsetup[T0MEANSsetup_rowcount,5 ];
int LAMBDAsetup[LAMBDAsetup_rowcount,5 ];
int DRIFTsetup[DRIFTsetup_rowcount,5 ];
int DIFFUSIONsetup[DIFFUSIONsetup_rowcount,5 ];
int MANIFESTVARsetup[MANIFESTVARsetup_rowcount,5 ];
int MANIFESTMEANSsetup[MANIFESTMEANSsetup_rowcount,5 ];
int CINTsetup[CINTsetup_rowcount,5 ];
int T0VARsetup[T0VARsetup_rowcount,5 ];
int TDPREDEFFECTsetup[TDPREDEFFECTsetup_rowcount,5 ];
  matrix[T0MEANSsetup_rowcount, 5] T0MEANSvalues;
matrix[LAMBDAsetup_rowcount, 5] LAMBDAvalues;
matrix[DRIFTsetup_rowcount, 5] DRIFTvalues;
matrix[DIFFUSIONsetup_rowcount, 5] DIFFUSIONvalues;
matrix[MANIFESTVARsetup_rowcount, 5] MANIFESTVARvalues;
matrix[MANIFESTMEANSsetup_rowcount, 5] MANIFESTMEANSvalues;
matrix[CINTsetup_rowcount, 5] CINTvalues;
matrix[T0VARsetup_rowcount, 5] T0VARvalues;
matrix[TDPREDEFFECTsetup_rowcount, 5] TDPREDEFFECTvalues;
  int TIPREDEFFECTsetup[nparams, ntipred];
}
      
transformed data{
  matrix[nlatent,nlatent] IIlatent;
  matrix[nlatent*nlatent,nlatent*nlatent] IIlatent2;
  int nlatentpop;
  //int ncont_y[ndatapoints];
  //int whichcont_y[ndatapoints, nmanifest];

  //ncont_y = nobs_y;
  //whichcont_y = whichobs_y;

  nlatentpop = ukfpop ? nlatent + nindvarying : nlatent;
  IIlatent = diag_matrix(rep_vector(1,nlatent));
  IIlatent2 = diag_matrix(rep_vector(1,nlatent*nlatent));
}
      
parameters {
  vector[nparams] rawpopmeans; // population level means 

  vector[nindvarying] rawpopsdbase; //population level std dev
  vector[nindvaryingoffdiagonals] sqrtpcov;
  vector[ukfpop ? 0 : nindvarying*nsubjects] baseindparams; //vector of subject level deviations, on the raw scale
  
  vector[ntipredeffects] tipredeffectparams; // effects of time independent covariates
  vector[nmissingtipreds] tipredsimputed;
  
  vector[intoverstates ? 0 : nlatent*ndatapoints] etaupdbasestates; //sampled latent states posterior
  //real<lower=1e-5,upper=5> ukfscale;
  
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  //matrix[nindvarying,nindvarying] rawpopcov;
  matrix[nindvarying,nindvarying] rawpopcorrsqrt;
  matrix[nindvarying,nindvarying] rawpopcovsqrt; 

  
  matrix[ T0MEANSsetup_rowcount ? max(T0MEANSsetup[,1]) : 0, T0MEANSsetup_rowcount ? max(T0MEANSsetup[,2]) : 0 ] T0MEANS[T0MEANSsubindex[nsubjects]];
matrix[ LAMBDAsetup_rowcount ? max(LAMBDAsetup[,1]) : 0, LAMBDAsetup_rowcount ? max(LAMBDAsetup[,2]) : 0 ] LAMBDA[LAMBDAsubindex[nsubjects]];
matrix[ DRIFTsetup_rowcount ? max(DRIFTsetup[,1]) : 0, DRIFTsetup_rowcount ? max(DRIFTsetup[,2]) : 0 ] DRIFT[DRIFTsubindex[nsubjects]];
matrix[ DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,1]) : 0, DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,2]) : 0 ] DIFFUSION[DIFFUSIONsubindex[nsubjects]];
matrix[ MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,1]) : 0, MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,2]) : 0 ] MANIFESTVAR[MANIFESTVARsubindex[nsubjects]];
matrix[ MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,1]) : 0, MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,2]) : 0 ] MANIFESTMEANS[MANIFESTMEANSsubindex[nsubjects]];
matrix[ CINTsetup_rowcount ? max(CINTsetup[,1]) : 0, CINTsetup_rowcount ? max(CINTsetup[,2]) : 0 ] CINT[CINTsubindex[nsubjects]];
matrix[ T0VARsetup_rowcount ? max(T0VARsetup[,1]) : 0, T0VARsetup_rowcount ? max(T0VARsetup[,2]) : 0 ] T0VAR[T0VARsubindex[nsubjects]];
matrix[ TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,1]) : 0, TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,2]) : 0 ] TDPREDEFFECT[TDPREDEFFECTsubindex[nsubjects]];

  matrix[nlatent,nlatent] asymDIFFUSION[ lineardynamics ? asymDIFFUSIONsubindex[nsubjects] : 0]; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] asymCINT[nt0meansstationary ? asymCINTsubindex[nsubjects] : 0]; // latent process asymptotic level
  
  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipreds; //tipred values to fill from data and, when needed, imputation vector
  matrix[nparams, ntipred] TIPREDEFFECT; //design matrix of individual time independent predictor effects

  if(ntipred > 0){ 
    {
    int counter;
      counter = 0;
      for(coli in 1:cols(tipreds)){ //insert missing ti predictors
        for(rowi in 1:rows(tipreds)){
          if(tipredsdata[rowi,coli]==99999) {
            counter = counter + 1;
            tipreds[rowi,coli] = tipredsimputed[counter];
          } else tipreds[rowi,coli] = tipredsdata[rowi,coli];
        }
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
    int counter;
    rawpopsd = exp(rawpopsdbase * 2 -2) .* sdscale;
    counter=0;
    for(j in 1:nindvarying){
      rawpopcovsqrt[j,j] = 1;
      for(i in 1:nindvarying){
        if(i > j){
          counter=counter+1;
          rawpopcovsqrt[i,j]=sqrtpcov[counter];
          rawpopcovsqrt[j,i]=sqrtpcov[counter];
        }
      }
    }
  rawpopcorrsqrt = covsqrt2corsqrt(rawpopcovsqrt,0);// cov2cors(tcrossprod(rawpopcovsqrt));
  //rawpopcov=quad_form_diag(rawpopcorr,rawpopsd);
  rawpopcovsqrt = diag_pre_multiply(rawpopsd, rawpopcorrsqrt); //chol(rawpopcov);
  }//end indvarying par setup

{
  vector[nparams] rawindparams;
  rawindparams = rawpopmeans;
  for(si in 1:nsubjects){

    if(ntipred==0 && nindvarying > 0 && ukfpop ==0) rawindparams[indvaryingindex] = 
      rawpopmeans[indvaryingindex] + rawpopcovsqrt * baseindparams[(1+(si-1)*nindvarying):(si*nindvarying)];

    if(ntipred > 0  && nindvarying > 0 && ukfpop ==0) rawindparams[indvaryingindex] = 
      rawpopmeans[indvaryingindex] + rawpopcovsqrt * baseindparams[(1+(si-1)*nindvarying):(si*nindvarying)] +
      TIPREDEFFECT[indvaryingindex,] * tipreds[si]';

  if(si <= T0MEANSsubindex[nsubjects]){
    for(ri in 1:size(T0MEANSsetup)){
      T0MEANS[si, T0MEANSsetup[ ri,1], T0MEANSsetup[ri,2]] = T0MEANSsetup[ri,3] ? tform(rawindparams[ T0MEANSsetup[ri,3] ], T0MEANSsetup[ri,4], T0MEANSvalues[ri,2], T0MEANSvalues[ri,3], T0MEANSvalues[ri,4] ) : T0MEANSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= LAMBDAsubindex[nsubjects]){
    for(ri in 1:size(LAMBDAsetup)){
      LAMBDA[si, LAMBDAsetup[ ri,1], LAMBDAsetup[ri,2]] = LAMBDAsetup[ri,3] ? tform(rawindparams[ LAMBDAsetup[ri,3] ], LAMBDAsetup[ri,4], LAMBDAvalues[ri,2], LAMBDAvalues[ri,3], LAMBDAvalues[ri,4] ) : LAMBDAvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= DRIFTsubindex[nsubjects]){
    for(ri in 1:size(DRIFTsetup)){
      DRIFT[si, DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = DRIFTsetup[ri,3] ? tform(rawindparams[ DRIFTsetup[ri,3] ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ) : DRIFTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= DIFFUSIONsubindex[nsubjects]){
    for(ri in 1:size(DIFFUSIONsetup)){
      DIFFUSION[si, DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = DIFFUSIONsetup[ri,3] ? tform(rawindparams[ DIFFUSIONsetup[ri,3] ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ) : DIFFUSIONvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= MANIFESTVARsubindex[nsubjects]){
    for(ri in 1:size(MANIFESTVARsetup)){
      MANIFESTVAR[si, MANIFESTVARsetup[ ri,1], MANIFESTVARsetup[ri,2]] = MANIFESTVARsetup[ri,3] ? tform(rawindparams[ MANIFESTVARsetup[ri,3] ], MANIFESTVARsetup[ri,4], MANIFESTVARvalues[ri,2], MANIFESTVARvalues[ri,3], MANIFESTVARvalues[ri,4] ) : MANIFESTVARvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= MANIFESTMEANSsubindex[nsubjects]){
    for(ri in 1:size(MANIFESTMEANSsetup)){
      MANIFESTMEANS[si, MANIFESTMEANSsetup[ ri,1], MANIFESTMEANSsetup[ri,2]] = MANIFESTMEANSsetup[ri,3] ? tform(rawindparams[ MANIFESTMEANSsetup[ri,3] ], MANIFESTMEANSsetup[ri,4], MANIFESTMEANSvalues[ri,2], MANIFESTMEANSvalues[ri,3], MANIFESTMEANSvalues[ri,4] ) : MANIFESTMEANSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= CINTsubindex[nsubjects]){
    for(ri in 1:size(CINTsetup)){
      CINT[si, CINTsetup[ ri,1], CINTsetup[ri,2]] = CINTsetup[ri,3] ? tform(rawindparams[ CINTsetup[ri,3] ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ) : CINTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= T0VARsubindex[nsubjects]){
    for(ri in 1:size(T0VARsetup)){
      T0VAR[si, T0VARsetup[ ri,1], T0VARsetup[ri,2]] = T0VARsetup[ri,3] ? tform(rawindparams[ T0VARsetup[ri,3] ], T0VARsetup[ri,4], T0VARvalues[ri,2], T0VARvalues[ri,3], T0VARvalues[ri,4] ) : T0VARvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= TDPREDEFFECTsubindex[nsubjects]){
    for(ri in 1:size(TDPREDEFFECTsetup)){
      TDPREDEFFECT[si, TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = TDPREDEFFECTsetup[ri,3] ? tform(rawindparams[ TDPREDEFFECTsetup[ri,3] ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ) : TDPREDEFFECTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      


  // perform any whole matrix transformations 
    
  if(si <= DIFFUSIONsubindex[nsubjects]) DIFFUSION[si] = sdcovsqrt2cov(DIFFUSION[si], lineardynamics * intoverstates ? 0 : 1);

  if(lineardynamics==1 && ndiffusion > 0){
    if(si <= asymDIFFUSIONsubindex[nsubjects]) {
      if(ndiffusion < nlatent) asymDIFFUSION[si] = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

      if(continuoustime==1) asymDIFFUSION[si, derrind, derrind] = to_matrix( 
      -( kron_prod( DRIFT[DRIFTsubindex[si], derrind, derrind ], IIlatent[ derrind, derrind ]) + 
         kron_prod(IIlatent[ derrind, derrind ], DRIFT[ DRIFTsubindex[si], derrind, derrind ]) ) \ 
      to_vector( DIFFUSION[ DIFFUSIONsubindex[si], derrind, derrind ]), ndiffusion,ndiffusion);

      if(continuoustime==0) asymDIFFUSION[si, derrind, derrind] = to_matrix( (IIlatent2[ derrind, derrind ] - 
        kron_prod(DRIFT[ DRIFTsubindex[si], derrind, derrind  ], 
          DRIFT[ DRIFTsubindex[si], derrind, derrind  ])) * 
        to_vector(DIFFUSION[ DIFFUSIONsubindex[si], derrind, derrind  ]) , ndiffusion, ndiffusion);
    } //end asymdiffusion loops
  }
          
    if(nt0meansstationary > 0){
      if(si <= asymCINTsubindex[nsubjects]){
        if(continuoustime==1) asymCINT[si] =  -DRIFT[ DRIFTsubindex[si] ] \ CINT[ CINTsubindex[si], ,1 ];
        if(continuoustime==0) asymCINT[si] =  (IIlatent - DRIFT[ DRIFTsubindex[si] ]) \ CINT[ CINTsubindex[si], ,1 ];
      }
    }

          
    if(binomial==0){
      if(si <= MANIFESTVARsubindex[nsubjects]) {
         for(ri in 1:nmanifest) MANIFESTVAR[si,ri,ri] = square(MANIFESTVAR[si,ri,ri]);
      }
    }
          
          
    if(si <= T0VARsubindex[nsubjects]) {
      T0VAR[si] = sdcovsqrt2cov(T0VAR[si],lineardynamics * intoverstates ? 0 : 1);

      if(nt0varstationary > 0) for(rowi in 1:nt0varstationary){
        T0VAR[si,t0varstationary[rowi,1],t0varstationary[rowi,2] ] = 
          asymDIFFUSION[si,t0varstationary[rowi,1],t0varstationary[rowi,2] ];
        T0VAR[si,t0varstationary[rowi,2],t0varstationary[rowi,1] ] = 
          asymDIFFUSION[si,t0varstationary[rowi,2],t0varstationary[rowi,1] ];
      }
    }

    
    if(nt0meansstationary > 0){
      if(si <= T0MEANSsubindex[nsubjects]) {
        for(rowi in 1:nt0meansstationary){
          T0MEANS[si,t0meansstationary[rowi,1] , 1] = 
            asymCINT[ asymCINTsubindex[si], t0meansstationary[rowi,1] ];
        }
      }
    }
  }
} //end subject loop
  

}
      
model{
  real ll;

  if(nopriors==0){
    target += normal_lpdf(rawpopmeans|0,1);
  
    if(ntipred > 0){ 
     tipredeffectparams ~ normal(0,1); 
     tipredsimputed ~ normal(0,10);
    }
    
    if(nindvarying > 0){
      if(nindvarying >1) sqrtpcov ~ normal(0,1);
      if(ukfpop==0) baseindparams ~ normal(0,1);
      rawpopsdbase ~ normal(0,1);
    }

  } //end pop priors section
  
  if(intoverstates==0)etaupdbasestates ~ normal(0,1);
  
  ll = 0;{
  int si;
  int counter;
  vector[nlatentpop] etaprior[ndatapoints]; //prior for latent states
  vector[nlatentpop] etaupd[ndatapoints]; //updated latent states
  matrix[nlatentpop, nlatentpop] etapriorcov[ndatapoints]; //prior for covariance of latent states
  matrix[nlatentpop, nlatentpop] etaupdcov[ndatapoints]; //updated covariance of latent states

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] ypred;
  vector[ukf ? nmanifest : 0] ystate;
  matrix[nmanifest, nmanifest] ypredcov;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypredcov_sqrt; 

  
  //likelihood
  vector[intoverstates ? sum(nobs_y): sum(ncont_y)] errtrans; // collection of prediction errors transformed to standard normal
  vector[intoverstates ? sum(nobs_y): sum(ncont_y)] errscales; // collection of prediction error scaling factors
  int cobscount; // counter summing over number of non missing observations in each row
  int nobsi; 

  //ukf
  matrix[ukf ? nlatentpop :0,ukf ? nlatentpop :0] sigpoints;
  vector[ukf ? nlatent :0] state; //dynamic portion of current states
  real dynerror; //dynamic error variable
  real k;
  real asquared;
  real l;
  real sqrtukfadjust;
  int ndynerror; // number of variance elements to include
  vector[lineardynamics ? 0 : nlatent] rkstates[lineardynamics ? 0 : 5]; //runge kutta integration steps

  //linear continuous time calcs
  matrix[lineardynamics ? nlatent : 0,lineardynamics ? nlatent : 0] discreteDRIFT;
  vector[lineardynamics ? nlatent : 0] discreteCINT;
  matrix[lineardynamics ? nlatent : 0, lineardynamics ? nlatent : 0] discreteDIFFUSION;

  // create simple, modifiable forms of the system matrices for easier use in the filter
  matrix[nlatent,1] sT0MEANS;
  matrix[nlatent,nlatent] sT0VAR;
  matrix[nlatent,nlatent] sDIFFUSION; 
  matrix[nlatent,nlatent] sasymDIFFUSION; 
  matrix[nlatent,nlatent] sDRIFT; 
  matrix[nlatent,1] sCINT;
  matrix[nmanifest,nmanifest] sMANIFESTVAR; 
  matrix[nmanifest,1] sMANIFESTMEANS;
  matrix[nmanifest,nlatent] sLAMBDA;
  matrix[ntdpred ? nlatent : 0,ntdpred] sTDPREDEFFECT;

  //ukf approximation parameters
  if(ukf==1) k = 0.5;
  
  cobscount=0; //running total of observed indicators treated as continuous

  for(rowi in 1:ndatapoints){
    int o[nobs_y[rowi]]; //which indicators are observed
    int o1[nbinary_y[rowi]]; //which indicators are observed and binary
    int o0[ncont_y[rowi]]; //which indicators are observed and continuous

    matrix[ukf ? nlatentpop : 0, 2*(nlatentpop+ (T0check[rowi] ? nlatent : ndiffusion)) +2 ] ukfstates; //sampled states relevant for dynamics
    matrix[ukf ? nmanifest : 0 , 2*(nlatentpop+(T0check[rowi] ? nlatent : ndiffusion))+2] ukfmeasures; // expected measures based on sampled states

    o = whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
    si=subject[rowi];
    nobsi = nobs_y[rowi]; //number of obs this row

    o1 = whichbinary_y[rowi,1:nbinary_y[rowi]];
    o0 = whichcont_y[rowi,1:ncont_y[rowi]];
    
    if(rowi!=1 && intoverstates==1) cobscount += nobs_y[rowi-1]; // number of non missing observations, treated as gaussian, until now
    if(rowi!=1 && intoverstates==0) cobscount += ncont_y[rowi-1]; // number of non missing observations, treated as gaussian, until now

    if(ukf==1){ //ukf approximation parameters
      if(T0check[rowi] == 1) { ndynerror = nlatent; } else ndynerror = ndiffusion;
      if(T0check[rowi]==1 || ( ndiffusion < nlatent && T0check[rowi-1]==1)) {
        asquared =  2.0/sqrt(nlatentpop+ndynerror) * 1e-1;
        l = asquared * (nlatentpop + ndynerror + k) - (nlatentpop+ndynerror); 
        sqrtukfadjust = sqrt(nlatentpop + ndynerror +l);
      }
    }

    if(T0check[rowi] == 1) { // calculate initial matrices if this is first row for si

      if( si ==1 || (ukfpop==0 && nindvarying > 0) ){ //whenever we need to get new subject parameters
        sT0MEANS = T0MEANS[ T0MEANSsubindex[si]];
        sT0VAR = T0VAR[ T0VARsubindex[si] ];
        sDRIFT = DRIFT[ DRIFTsubindex[si] ];
        sDIFFUSION = DIFFUSION[ DIFFUSIONsubindex[si] ];
        sCINT = CINT[ CINTsubindex[si] ];
        sLAMBDA = LAMBDA[ LAMBDAsubindex[si] ];
        sMANIFESTMEANS = MANIFESTMEANS[ MANIFESTMEANSsubindex[si] ];
        sMANIFESTVAR = MANIFESTVAR[ MANIFESTVARsubindex[si] ];
        sTDPREDEFFECT = TDPREDEFFECT[ TDPREDEFFECTsubindex[si] ];
        
        if(1==99 && lineardynamics==1 && (rowi==1 || asymDIFFUSIONsubindex[si] != asymDIFFUSIONsubindex[si-1])){ 
          sasymDIFFUSION[derrind,derrind] = to_matrix( 
            -( kron_prod(sDRIFT[derrind,derrind], IIlatent[derrind,derrind]) + kron_prod(IIlatent[derrind,derrind], sDRIFT[derrind,derrind]) ) \ 
            to_vector(sDIFFUSION[derrind,derrind] ), 
            ndiffusion, ndiffusion);
        }
      }

      if(ukf==1){
        etaprior[rowi,] = rep_vector(0,nlatentpop); // because some values stay zero
        sigpoints = rep_matrix(0, nlatentpop,nlatentpop);
      
        if(ukfpop==1) {
          if(ntipred ==0) etaprior[rowi, (nlatent+1):(nlatentpop)] = rawpopmeans[indvaryingindex];
          if(ntipred >0) etaprior[rowi, (nlatent+1):(nlatentpop)] = rawpopmeans[indvaryingindex] + TIPREDEFFECT[indvaryingindex] * tipreds[si]';
          sigpoints[(nlatent+1):(nlatentpop), (nlatent+1):(nlatentpop)] = rawpopcovsqrt * sqrtukfadjust;
        }
      }

      if(ukf==0){
      if(ntdpred ==0)etaprior[rowi] = sT0MEANS[,1]; //prior for initial latent state
      if(ntdpred > 0) etaprior[rowi] = TDPREDEFFECT[ TDPREDEFFECTsubindex[si] ] * tdpreds[rowi] + sT0MEANS[,1];
      etapriorcov[rowi] =  T0VAR[ T0VARsubindex[si] ];
      }

    } //end T0 matrices

    if(lineardynamics==1 && ukf==0 && T0check[rowi]==0){ //linear kf time update
    
      if(continuoustime ==1){
        int dtchange = 0;
        if(si==1 && T0check[rowi -1] == 1) {
          dtchange = 1;
        } else if(T0check[rowi-1] == 1 && dT[rowi-2] != dT[rowi]){
          dtchange = 1;
        } else if(T0check[rowi-1] == 0 && dT[rowi-1] != dT[rowi]) dtchange = 1;
        
        if(dtchange==1 || (T0check[rowi-1]==1 && si <= DRIFTsubindex[si])){
          if(driftdiagonly==1) discreteDRIFT = matrix_diagexp(sDRIFT * dT[rowi]);
          if(driftdiagonly==0) discreteDRIFT = matrix_exp(sDRIFT * dT[rowi]);
        }
        if(dtchange==1 || (T0check[rowi-1]==1 && (si <= CINTsubindex[si] || si <= DRIFTsubindex[si]))){
          discreteCINT = sDRIFT \ (discreteDRIFT - IIlatent) * sCINT[,1];
        }
    
        if(dtchange==1 || (T0check[rowi-1]==1 && (si <= DIFFUSIONsubindex[si]|| si <= DRIFTsubindex[si]))){
          //discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
            //quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]' );
          discreteDIFFUSION[derrind, derrind] = discreteDIFFUSIONcalc(DRIFT[ DRIFTsubindex[si], derrind, derrind], sDIFFUSION[derrind, derrind], dT[rowi]);
        }
      }
  
      if(continuoustime==0 && T0check[rowi-1] == 1){
        discreteDRIFT=sDRIFT;
        discreteCINT=sCINT[,1];
        discreteDIFFUSION=sDIFFUSION;
      }

      etaprior[rowi] = discreteDRIFT * etaupd[rowi-1] + discreteCINT;
      if(intoverstates==1) {
        if(ndiffusion ==0) etapriorcov[rowi] = quad_form(etaupdcov[rowi-1], discreteDRIFT');
        if(ndiffusion > 0) etapriorcov[rowi,derrind,derrind] = quad_form(etaupdcov[rowi-1], discreteDRIFT') + discreteDIFFUSION[derrind,derrind];
      }
    }//end linear time update

    if(ukf==1){ //ukf time update

      if(T0check[rowi]==1) dynerror = sqrtukfadjust;
      if(T0check[rowi]==0 && lineardynamics==0) dynerror = sqrtukfadjust / sqrt(dT[rowi]); //Weiner process variance adjustment
  
      if(T0check[rowi]==0){ //compute updated sigpoints
        sigpoints = cholspd(etaupdcov[rowi-1,,]) * sqrtukfadjust;
        etaprior[rowi,] = etaupd[rowi-1,];
      }
  
      //configure ukf states
      for(statei in 1:cols(ukfstates)){
        if(statei > (2+nlatentpop) && statei <= 2+2*nlatentpop){
          ukfstates[, statei] = etaprior[rowi,] - sigpoints[,statei-(2+nlatentpop)];
        } else
        if(statei > 2 && statei <= 2+2*nlatentpop){
          ukfstates[, statei] = etaprior[rowi,] + sigpoints[,statei-2]; 
        } else
          ukfstates[, statei] = etaprior[rowi,]; 
      }
  
      for(statei in 2:cols(ukfstates) ){ //for each ukf state sample
    
        if(T0check[rowi]==1){
     
          if(statei <= 2+2*nlatentpop+1){
          
      if(ukfpop==1){
        
    for(ri in 1:size(T0MEANSsetup)){
      if(T0MEANSsetup[ ri,5] > 0){ 
        sT0MEANS[T0MEANSsetup[ ri,1], T0MEANSsetup[ri,2]] = tform(ukfstates[nlatent +T0MEANSsetup[ri,5], statei ], T0MEANSsetup[ri,4], T0MEANSvalues[ri,2], T0MEANSvalues[ri,3], T0MEANSvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(T0VARsetup)){
      if(T0VARsetup[ ri,5] > 0){ 
        sT0VAR[T0VARsetup[ ri,1], T0VARsetup[ri,2]] = tform(ukfstates[nlatent +T0VARsetup[ri,5], statei ], T0VARsetup[ri,4], T0VARvalues[ri,2], T0VARvalues[ri,3], T0VARvalues[ri,4] ); 
      }
    }};
          }
  
          state = sT0MEANS[,1];
          if(statei > (2+2*nlatentpop+ndynerror)) {
            state += -sT0VAR[ , statei - (2+2*nlatentpop+ndynerror) ] * dynerror; 
          } else
          if(statei > (2+2*nlatentpop))  state += sT0VAR[ , statei - (2+2*nlatentpop) ] * dynerror; 
  
        } 
    
        if(T0check[rowi]==0){
          state = ukfstates[1:nlatent, statei];
    
          if(lineardynamics==0){
           for(stepi in 1:integrationsteps[rowi]){ //for each euler integration step
              rkstates[5] = state; //store initial states for this integration step
    
              for(ki in 1:4){ //runge kutta integration within euler scheme
                if(ki==2 || ki==3) state=rkstates[5] + dTsmall[rowi] /2 * rkstates[ki-1];
                if(ki==4) state = rkstates[5] + dTsmall[rowi] * rkstates[3];
                  if(ukfpop==1){

    for(ri in 1:size(DRIFTsetup)){
      if(DRIFTsetup[ ri,5] > 0){ 
        sDRIFT[DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = 
          tform(ukfstates[nlatent +DRIFTsetup[ri,5], statei ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(DIFFUSIONsetup)){
      if(DIFFUSIONsetup[ ri,5] > 0){ 
        sDIFFUSION[DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = 
          tform(ukfstates[nlatent +DIFFUSIONsetup[ri,5], statei ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(CINTsetup)){
      if(CINTsetup[ ri,5] > 0){ 
        sCINT[CINTsetup[ ri,1], CINTsetup[ri,2]] = 
          tform(ukfstates[nlatent +CINTsetup[ri,5], statei ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(TDPREDEFFECTsetup)){
      if(TDPREDEFFECTsetup[ ri,5] > 0){ 
        sTDPREDEFFECT[TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = 
          tform(ukfstates[nlatent +TDPREDEFFECTsetup[ri,5], statei ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ); 
      }
    }};
  
                if(statei <= (2+2*nlatentpop) ) {
                  rkstates[ki] = sDRIFT * state + sCINT[,1];
                } else if(statei <= (2+2*nlatentpop + ndynerror) ){
                  rkstates[ki] = sDRIFT * state + sCINT[,1] + sDIFFUSION[ , derrind[statei - (2+2*nlatentpop)] ] * dynerror; 
                } else rkstates[ki] = sDRIFT * state + sCINT[,1] - sDIFFUSION[ , derrind[statei - (2+2*nlatentpop + ndynerror)] ] * dynerror;
  
              }
              state = (rkstates[5] + dTsmall[rowi]/6  *(rkstates[1]+2*rkstates[2]+2*rkstates[3]+rkstates[4])); //integrate over rk steps
            }
          } //end nonlinear time update
    
    
          if(lineardynamics==1){ //this could be much more efficient...
              if(ukfpop==1){

    for(ri in 1:size(DRIFTsetup)){
      if(DRIFTsetup[ ri,5] > 0){ 
        sDRIFT[DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = 
          tform(ukfstates[nlatent +DRIFTsetup[ri,5], statei ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(DIFFUSIONsetup)){
      if(DIFFUSIONsetup[ ri,5] > 0){ 
        sDIFFUSION[DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = 
          tform(ukfstates[nlatent +DIFFUSIONsetup[ri,5], statei ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(CINTsetup)){
      if(CINTsetup[ ri,5] > 0){ 
        sCINT[CINTsetup[ ri,1], CINTsetup[ri,2]] = 
          tform(ukfstates[nlatent +CINTsetup[ri,5], statei ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(TDPREDEFFECTsetup)){
      if(TDPREDEFFECTsetup[ ri,5] > 0){ 
        sTDPREDEFFECT[TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = 
          tform(ukfstates[nlatent +TDPREDEFFECTsetup[ri,5], statei ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ); 
      }
    }};
            if(statei <= 2+2*nlatentpop+1){ //because after this its only noise variables
              if(continuoustime ==1){ 
                if(driftdiagonly==1) discreteDRIFT = matrix_diagexp(sDRIFT * dT[rowi]);
                if(driftdiagonly==0) discreteDRIFT = matrix_exp(sDRIFT * dT[rowi]);
                discreteCINT = sDRIFT \ (discreteDRIFT - IIlatent) * sCINT[,1];
                //discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
                  //quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]' );
                discreteDIFFUSION[derrind, derrind] = discreteDIFFUSIONcalc(sDRIFT[derrind, derrind], sDIFFUSION[derrind, derrind],dT[rowi]);;
              }
              if(continuoustime==0){
                discreteDRIFT=sDRIFT;
                discreteCINT=sCINT[,1];
                discreteDIFFUSION=sDIFFUSION;
              }
            }
            state = discreteDRIFT * state + discreteCINT;
            if(statei > (2+2*nlatentpop) && statei <= (2+2*nlatentpop+ndynerror) )  state += discreteDIFFUSION[ , derrind[statei - (2+2*nlatentpop)] ] * dynerror; 
            if(statei > (2+2*nlatentpop+ndynerror)) state += -discreteDIFFUSION[ , derrind[statei - (2+2*nlatentpop+ndynerror)] ] * dynerror;
          }
        }  // end of if not t0 section (time update)
    
          if(ukfpop==1){

    for(ri in 1:size(DRIFTsetup)){
      if(DRIFTsetup[ ri,5] > 0){ 
        sDRIFT[DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = 
          tform(ukfstates[nlatent +DRIFTsetup[ri,5], statei ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(DIFFUSIONsetup)){
      if(DIFFUSIONsetup[ ri,5] > 0){ 
        sDIFFUSION[DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = 
          tform(ukfstates[nlatent +DIFFUSIONsetup[ri,5], statei ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(CINTsetup)){
      if(CINTsetup[ ri,5] > 0){ 
        sCINT[CINTsetup[ ri,1], CINTsetup[ri,2]] = 
          tform(ukfstates[nlatent +CINTsetup[ri,5], statei ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(TDPREDEFFECTsetup)){
      if(TDPREDEFFECTsetup[ ri,5] > 0){ 
        sTDPREDEFFECT[TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = 
          tform(ukfstates[nlatent +TDPREDEFFECTsetup[ri,5], statei ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ); 
      }
    }};
        if(ntdpred > 0) state +=  (sTDPREDEFFECT * tdpreds[rowi]); //tdpred effect only influences at observed time point
        ukfstates[1:nlatent, statei] = state; //now contains time updated state
        if(statei==2) ukfstates[1:nlatent, 1] = state; //mean goes in twice for weighting
      }

      etaprior[rowi] = colMeans(ukfstates');
      etapriorcov[rowi] = cov_of_matrix(ukfstates') / asquared;
    } //end ukf time update

    if(intoverstates==0 && lineardynamics == 1) {
      if(T0check[rowi]==1) etaupd[rowi] = etaprior[rowi] +  sT0VAR * etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
      if(T0check[rowi]==0) etaupd[rowi] = etaprior[rowi] +  sDIFFUSION * etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
    }
  
    if(nobsi==0 && intoverstates==1) {
      etaupdcov[rowi] = etapriorcov[rowi];
      etaupd[rowi] = etaprior[rowi];
    }

    if (nobsi > 0) {  // if some observations create right size matrices for missingness and calculate...
  
      int cindex[intoverstates ? nobsi : ncont_y[rowi]];

      if(intoverstates==0) cindex = o0;
      if(intoverstates==1) cindex = o; //treat all obs as continuous gaussian

      if(ukf==0){ //non ukf measurement
        if(intoverstates==1) { //classic kalman
          ypred[o] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * etaprior[rowi];
          ypredcov[o,o] = quad_form(etapriorcov[rowi], sLAMBDA[o,]') + sMANIFESTVAR[o,o];
          for(wi in 1:nmanifest){ 
            if(manifesttype[wi]==1 && Y[rowi,wi] != 99999) ypredcov[wi,wi] = ypredcov[wi,wi] + fabs((ypred[wi] - 1) .* (ypred[wi])); 
          }
          K[,o] = mdivide_right(etapriorcov[rowi] * sLAMBDA[o,]', ypredcov[o,o]); 
          etaupdcov[rowi] = (IIlatent - K[,o] * sLAMBDA[o,]) * etapriorcov[rowi];
        }
        if(intoverstates==0) { //sampled states
          if(ncont_y[rowi] > 0) ypred[o0] = sMANIFESTMEANS[o0,1] + sLAMBDA[o0,] * etaupd[rowi];
          if(nbinary_y[rowi] > 0) ypred[o1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * etaupd[rowi])));
          ypredcov[o,o] = sMANIFESTVAR[o,o];
        }
      }
  

      if(ukf==1){ //ukf measurement
        for(statei in 2:cols(ukfmeasures)){
          state = ukfstates[ 1:nlatent, statei];
          
      if(ukfpop==1){
        
    for(ri in 1:size(LAMBDAsetup)){
      if(LAMBDAsetup[ ri,5] > 0){ 
        sLAMBDA[LAMBDAsetup[ ri,1], LAMBDAsetup[ri,2]] = tform(ukfstates[nlatent +LAMBDAsetup[ri,5], statei ], LAMBDAsetup[ri,4], LAMBDAvalues[ri,2], LAMBDAvalues[ri,3], LAMBDAvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(MANIFESTMEANSsetup)){
      if(MANIFESTMEANSsetup[ ri,5] > 0){ 
        sMANIFESTMEANS[MANIFESTMEANSsetup[ ri,1], MANIFESTMEANSsetup[ri,2]] = tform(ukfstates[nlatent +MANIFESTMEANSsetup[ri,5], statei ], MANIFESTMEANSsetup[ri,4], MANIFESTMEANSvalues[ri,2], MANIFESTMEANSvalues[ri,3], MANIFESTMEANSvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(MANIFESTVARsetup)){
      if(MANIFESTVARsetup[ ri,5] > 0){ 
        sMANIFESTVAR[MANIFESTVARsetup[ ri,1], MANIFESTVARsetup[ri,2]] = tform(ukfstates[nlatent +MANIFESTVARsetup[ri,5], statei ], MANIFESTVARsetup[ri,4], MANIFESTVARvalues[ri,2], MANIFESTVARvalues[ri,3], MANIFESTVARvalues[ri,4] ); 
      }
    }};

  
          if(ncont_y[rowi] > 0) ukfmeasures[o0 , statei] = sMANIFESTMEANS[o0,1] + sLAMBDA[o0,] * state;
          if(nbinary_y[rowi] > 0) {
            ukfmeasures[o1 , statei] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * state)));
          }
          if(statei==2) { //temporary measure to get mean in twice
          if(ncont_y[rowi] > 0) ukfmeasures[o0 , 1] = sMANIFESTMEANS[o0,1] + sLAMBDA[o0,] * state;
          if(nbinary_y[rowi] > 0) {
            ukfmeasures[o1 , 1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * state)));
          }
          }
        } //end ukf measurement loop
    
        ypred[o] = colMeans(ukfmeasures[o,]'); 
        ypredcov[o,o] = cov_of_matrix(ukfmeasures[o,]') /asquared + sMANIFESTVAR[o,o];
        for(wi in 1:nmanifest){ 
          if(manifesttype[wi]==1 && Y[rowi,wi] != 99999) ypredcov[wi,wi] = ypredcov[wi,wi] + fabs((ypred[wi] - 1) .* (ypred[wi])); //sMANIFESTVAR[wi,wi] + (merror[wi] / cols(ukfmeasures) +1e-8);
        }
        K[,o] = mdivide_right(crosscov(ukfstates', ukfmeasures[o,]') /asquared, ypredcov[o,o]); 
        etaupdcov[rowi] = etapriorcov[rowi] - quad_form(ypredcov[o,o],  K[,o]');
      }

      
  
      err[o] = Y[rowi,o] - ypred[o]; // prediction error
      if(intoverstates==1) etaupd[rowi,] = etaprior[rowi,] + (K[,o] * err[o]);

  
      
      if(intoverstates==0 && nbinary_y[rowi] > 0) ll += sum(log( Y[rowi,o1] .* (ypred[o1]) + (1-Y[rowi,o1]) .* (1-ypred[o1])));

      if(verbose > 1) {
        print("rowi ",rowi, "  si ", si, "  etaprior[rowi] ",etaprior[rowi],"  etapriorcov[rowi] ",etapriorcov[rowi],
          "  etaupd[rowi] ",etaupd[rowi],"  etaupdcov[rowi] ",etaupdcov[rowi],"  ypred ",ypred,"  ypredcov ",ypredcov, "  K ",K,
          "  sDRIFT ", sDRIFT, " sDIFFUSION ", sDIFFUSION, " sCINT ", sCINT, "  sMANIFESTVAR ", sMANIFESTVAR, "  rawpopsd ", rawpopsd,
          "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        if(lineardynamics==1) print("discreteDRIFT ",discreteDRIFT,"  discreteCINT ", discreteCINT, "  discreteDIFFUSION ", discreteDIFFUSION)
      }
      if(verbose > 2) print("ukfstates ", ukfstates, "  ukfmeasures ", ukfmeasures);

      if(size(cindex) > 0){
         ypredcov_sqrt[cindex,cindex]=cholspd(ypredcov[cindex,cindex]);
         errtrans[(cobscount+1):(cobscount+size(cindex))] = mdivide_left_tri_low(ypredcov_sqrt[cindex,cindex], err[cindex]); //transform pred errors to standard normal dist and collect
         errscales[(cobscount+1):(cobscount+size(cindex))] = log(diagonal(ypredcov_sqrt[cindex,cindex])); //account for transformation of scale in loglik
      }
    }//end nobs > 0 section
  }//end rowi

  if(intoverstates==1 || sum(ncont_y) > 0) ll = ll + normal_lpdf(errtrans|0,1) - sum(errscales);
}
  target += ll;
  
  if(verbose > 0) print("lp = ", target());
  
}
generated quantities{
  vector[nparams] popmeans;
  vector[nparams] popsd;
  matrix[nindvarying,nindvarying] rawpopcov;
  matrix[nindvarying,nindvarying] rawpopcorr;
  matrix[ T0MEANSsetup_rowcount ? max(T0MEANSsetup[,1]) : 0, T0MEANSsetup_rowcount ? max(T0MEANSsetup[,2]) : 0 ] pop_T0MEANS[T0MEANSsubindex[1]];
matrix[ LAMBDAsetup_rowcount ? max(LAMBDAsetup[,1]) : 0, LAMBDAsetup_rowcount ? max(LAMBDAsetup[,2]) : 0 ] pop_LAMBDA[LAMBDAsubindex[1]];
matrix[ DRIFTsetup_rowcount ? max(DRIFTsetup[,1]) : 0, DRIFTsetup_rowcount ? max(DRIFTsetup[,2]) : 0 ] pop_DRIFT[DRIFTsubindex[1]];
matrix[ DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,1]) : 0, DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,2]) : 0 ] pop_DIFFUSION[DIFFUSIONsubindex[1]];
matrix[ MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,1]) : 0, MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,2]) : 0 ] pop_MANIFESTVAR[MANIFESTVARsubindex[1]];
matrix[ MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,1]) : 0, MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,2]) : 0 ] pop_MANIFESTMEANS[MANIFESTMEANSsubindex[1]];
matrix[ CINTsetup_rowcount ? max(CINTsetup[,1]) : 0, CINTsetup_rowcount ? max(CINTsetup[,2]) : 0 ] pop_CINT[CINTsubindex[1]];
matrix[ T0VARsetup_rowcount ? max(T0VARsetup[,1]) : 0, T0VARsetup_rowcount ? max(T0VARsetup[,2]) : 0 ] pop_T0VAR[T0VARsubindex[1]];
matrix[ TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,1]) : 0, TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,2]) : 0 ] pop_TDPREDEFFECT[TDPREDEFFECTsubindex[1]];

  matrix[nlatent,nlatent] asympop_DIFFUSION[ lineardynamics ? asymDIFFUSIONsubindex[1] : 0]; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] asympop_CINT[nt0meansstationary ? asymCINTsubindex[1] : 0]; // latent process asymptotic level
  

vector[nmanifest] Ygen[ngenerations, ndatapoints];
for(geni in 1:ngenerations) Ygen[geni,,] = rep_array(rep_vector(0,nmanifest), ndatapoints);
for(geni in 1:ngenerations){

  int si;
  int counter;
  vector[nlatentpop] etaprior[ndatapoints]; //prior for latent states
  vector[nlatentpop] etaupd[ndatapoints]; //updated latent states
  matrix[nlatentpop, nlatentpop] etapriorcov[ndatapoints]; //prior for covariance of latent states
  matrix[nlatentpop, nlatentpop] etaupdcov[ndatapoints]; //updated covariance of latent states

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] ypred;
  vector[ukf ? nmanifest : 0] ystate;
  matrix[nmanifest, nmanifest] ypredcov;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypredcov_sqrt; 

  
  //likelihood
  vector[intoverstates ? sum(nobs_y): sum(ncont_y)] errtrans; // collection of prediction errors transformed to standard normal
  vector[intoverstates ? sum(nobs_y): sum(ncont_y)] errscales; // collection of prediction error scaling factors
  int cobscount; // counter summing over number of non missing observations in each row
  int nobsi; 

  //ukf
  matrix[ukf ? nlatentpop :0,ukf ? nlatentpop :0] sigpoints;
  vector[ukf ? nlatent :0] state; //dynamic portion of current states
  real dynerror; //dynamic error variable
  real k;
  real asquared;
  real l;
  real sqrtukfadjust;
  int ndynerror; // number of variance elements to include
  vector[lineardynamics ? 0 : nlatent] rkstates[lineardynamics ? 0 : 5]; //runge kutta integration steps

  //linear continuous time calcs
  matrix[lineardynamics ? nlatent : 0,lineardynamics ? nlatent : 0] discreteDRIFT;
  vector[lineardynamics ? nlatent : 0] discreteCINT;
  matrix[lineardynamics ? nlatent : 0, lineardynamics ? nlatent : 0] discreteDIFFUSION;

  // create simple, modifiable forms of the system matrices for easier use in the filter
  matrix[nlatent,1] sT0MEANS;
  matrix[nlatent,nlatent] sT0VAR;
  matrix[nlatent,nlatent] sDIFFUSION; 
  matrix[nlatent,nlatent] sasymDIFFUSION; 
  matrix[nlatent,nlatent] sDRIFT; 
  matrix[nlatent,1] sCINT;
  matrix[nmanifest,nmanifest] sMANIFESTVAR; 
  matrix[nmanifest,1] sMANIFESTMEANS;
  matrix[nmanifest,nlatent] sLAMBDA;
  matrix[ntdpred ? nlatent : 0,ntdpred] sTDPREDEFFECT;

  //ukf approximation parameters
  if(ukf==1) k = 0.5;
  
  cobscount=0; //running total of observed indicators treated as continuous

  for(rowi in 1:ndatapoints){
    int o[nobs_y[rowi]]; //which indicators are observed
    int o1[nbinary_y[rowi]]; //which indicators are observed and binary
    int o0[ncont_y[rowi]]; //which indicators are observed and continuous

    matrix[ukf ? nlatentpop : 0, 2*(nlatentpop+ (T0check[rowi] ? nlatent : ndiffusion)) +2 ] ukfstates; //sampled states relevant for dynamics
    matrix[ukf ? nmanifest : 0 , 2*(nlatentpop+(T0check[rowi] ? nlatent : ndiffusion))+2] ukfmeasures; // expected measures based on sampled states

    o = whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
    si=subject[rowi];
    nobsi = nobs_y[rowi]; //number of obs this row

    o1 = whichbinary_y[rowi,1:nbinary_y[rowi]];
    o0 = whichcont_y[rowi,1:ncont_y[rowi]];
    
    if(rowi!=1 && intoverstates==1) cobscount += nobs_y[rowi-1]; // number of non missing observations, treated as gaussian, until now
    if(rowi!=1 && intoverstates==0) cobscount += ncont_y[rowi-1]; // number of non missing observations, treated as gaussian, until now

    if(ukf==1){ //ukf approximation parameters
      if(T0check[rowi] == 1) { ndynerror = nlatent; } else ndynerror = ndiffusion;
      if(T0check[rowi]==1 || ( ndiffusion < nlatent && T0check[rowi-1]==1)) {
        asquared =  2.0/sqrt(nlatentpop+ndynerror) * 1e-1;
        l = asquared * (nlatentpop + ndynerror + k) - (nlatentpop+ndynerror); 
        sqrtukfadjust = sqrt(nlatentpop + ndynerror +l);
      }
    }

    if(T0check[rowi] == 1) { // calculate initial matrices if this is first row for si

      if( si ==1 || (ukfpop==0 && nindvarying > 0) ){ //whenever we need to get new subject parameters
        sT0MEANS = T0MEANS[ T0MEANSsubindex[si]];
        sT0VAR = T0VAR[ T0VARsubindex[si] ];
        sDRIFT = DRIFT[ DRIFTsubindex[si] ];
        sDIFFUSION = DIFFUSION[ DIFFUSIONsubindex[si] ];
        sCINT = CINT[ CINTsubindex[si] ];
        sLAMBDA = LAMBDA[ LAMBDAsubindex[si] ];
        sMANIFESTMEANS = MANIFESTMEANS[ MANIFESTMEANSsubindex[si] ];
        sMANIFESTVAR = MANIFESTVAR[ MANIFESTVARsubindex[si] ];
        sTDPREDEFFECT = TDPREDEFFECT[ TDPREDEFFECTsubindex[si] ];
        
        if(1==99 && lineardynamics==1 && (rowi==1 || asymDIFFUSIONsubindex[si] != asymDIFFUSIONsubindex[si-1])){ 
          sasymDIFFUSION[derrind,derrind] = to_matrix( 
            -( kron_prod(sDRIFT[derrind,derrind], IIlatent[derrind,derrind]) + kron_prod(IIlatent[derrind,derrind], sDRIFT[derrind,derrind]) ) \ 
            to_vector(sDIFFUSION[derrind,derrind] ), 
            ndiffusion, ndiffusion);
        }
      }

      if(ukf==1){
        etaprior[rowi,] = rep_vector(0,nlatentpop); // because some values stay zero
        sigpoints = rep_matrix(0, nlatentpop,nlatentpop);
      
        if(ukfpop==1) {
          if(ntipred ==0) etaprior[rowi, (nlatent+1):(nlatentpop)] = rawpopmeans[indvaryingindex];
          if(ntipred >0) etaprior[rowi, (nlatent+1):(nlatentpop)] = rawpopmeans[indvaryingindex] + TIPREDEFFECT[indvaryingindex] * tipreds[si]';
          sigpoints[(nlatent+1):(nlatentpop), (nlatent+1):(nlatentpop)] = rawpopcovsqrt * sqrtukfadjust;
        }
      }

      if(ukf==0){
      if(ntdpred ==0)etaprior[rowi] = sT0MEANS[,1]; //prior for initial latent state
      if(ntdpred > 0) etaprior[rowi] = TDPREDEFFECT[ TDPREDEFFECTsubindex[si] ] * tdpreds[rowi] + sT0MEANS[,1];
      etapriorcov[rowi] =  T0VAR[ T0VARsubindex[si] ];
      }

    } //end T0 matrices

    if(lineardynamics==1 && ukf==0 && T0check[rowi]==0){ //linear kf time update
    
      if(continuoustime ==1){
        int dtchange = 0;
        if(si==1 && T0check[rowi -1] == 1) {
          dtchange = 1;
        } else if(T0check[rowi-1] == 1 && dT[rowi-2] != dT[rowi]){
          dtchange = 1;
        } else if(T0check[rowi-1] == 0 && dT[rowi-1] != dT[rowi]) dtchange = 1;
        
        if(dtchange==1 || (T0check[rowi-1]==1 && si <= DRIFTsubindex[si])){
          if(driftdiagonly==1) discreteDRIFT = matrix_diagexp(sDRIFT * dT[rowi]);
          if(driftdiagonly==0) discreteDRIFT = matrix_exp(sDRIFT * dT[rowi]);
        }
        if(dtchange==1 || (T0check[rowi-1]==1 && (si <= CINTsubindex[si] || si <= DRIFTsubindex[si]))){
          discreteCINT = sDRIFT \ (discreteDRIFT - IIlatent) * sCINT[,1];
        }
    
        if(dtchange==1 || (T0check[rowi-1]==1 && (si <= DIFFUSIONsubindex[si]|| si <= DRIFTsubindex[si]))){
          //discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
            //quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]' );
          discreteDIFFUSION[derrind, derrind] = discreteDIFFUSIONcalc(DRIFT[ DRIFTsubindex[si], derrind, derrind], sDIFFUSION[derrind, derrind], dT[rowi]);
        }
      }
  
      if(continuoustime==0 && T0check[rowi-1] == 1){
        discreteDRIFT=sDRIFT;
        discreteCINT=sCINT[,1];
        discreteDIFFUSION=sDIFFUSION;
      }

      etaprior[rowi] = discreteDRIFT * etaupd[rowi-1] + discreteCINT;
      if(intoverstates==1) {
        if(ndiffusion ==0) etapriorcov[rowi] = quad_form(etaupdcov[rowi-1], discreteDRIFT');
        if(ndiffusion > 0) etapriorcov[rowi,derrind,derrind] = quad_form(etaupdcov[rowi-1], discreteDRIFT') + discreteDIFFUSION[derrind,derrind];
      }
    }//end linear time update

    if(ukf==1){ //ukf time update

      if(T0check[rowi]==1) dynerror = sqrtukfadjust;
      if(T0check[rowi]==0 && lineardynamics==0) dynerror = sqrtukfadjust / sqrt(dT[rowi]); //Weiner process variance adjustment
  
      if(T0check[rowi]==0){ //compute updated sigpoints
        sigpoints = cholspd(etaupdcov[rowi-1,,]) * sqrtukfadjust;
        etaprior[rowi,] = etaupd[rowi-1,];
      }
  
      //configure ukf states
      for(statei in 1:cols(ukfstates)){
        if(statei > (2+nlatentpop) && statei <= 2+2*nlatentpop){
          ukfstates[, statei] = etaprior[rowi,] - sigpoints[,statei-(2+nlatentpop)];
        } else
        if(statei > 2 && statei <= 2+2*nlatentpop){
          ukfstates[, statei] = etaprior[rowi,] + sigpoints[,statei-2]; 
        } else
          ukfstates[, statei] = etaprior[rowi,]; 
      }
  
      for(statei in 2:cols(ukfstates) ){ //for each ukf state sample
    
        if(T0check[rowi]==1){
     
          if(statei <= 2+2*nlatentpop+1){
          
      if(ukfpop==1){
        
    for(ri in 1:size(T0MEANSsetup)){
      if(T0MEANSsetup[ ri,5] > 0){ 
        sT0MEANS[T0MEANSsetup[ ri,1], T0MEANSsetup[ri,2]] = tform(ukfstates[nlatent +T0MEANSsetup[ri,5], statei ], T0MEANSsetup[ri,4], T0MEANSvalues[ri,2], T0MEANSvalues[ri,3], T0MEANSvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(T0VARsetup)){
      if(T0VARsetup[ ri,5] > 0){ 
        sT0VAR[T0VARsetup[ ri,1], T0VARsetup[ri,2]] = tform(ukfstates[nlatent +T0VARsetup[ri,5], statei ], T0VARsetup[ri,4], T0VARvalues[ri,2], T0VARvalues[ri,3], T0VARvalues[ri,4] ); 
      }
    }};
          }
  
          state = sT0MEANS[,1];
          if(statei > (2+2*nlatentpop+ndynerror)) {
            state += -sT0VAR[ , statei - (2+2*nlatentpop+ndynerror) ] * dynerror; 
          } else
          if(statei > (2+2*nlatentpop))  state += sT0VAR[ , statei - (2+2*nlatentpop) ] * dynerror; 
  
        } 
    
        if(T0check[rowi]==0){
          state = ukfstates[1:nlatent, statei];
    
          if(lineardynamics==0){
           for(stepi in 1:integrationsteps[rowi]){ //for each euler integration step
              rkstates[5] = state; //store initial states for this integration step
    
              for(ki in 1:4){ //runge kutta integration within euler scheme
                if(ki==2 || ki==3) state=rkstates[5] + dTsmall[rowi] /2 * rkstates[ki-1];
                if(ki==4) state = rkstates[5] + dTsmall[rowi] * rkstates[3];
                  if(ukfpop==1){

    for(ri in 1:size(DRIFTsetup)){
      if(DRIFTsetup[ ri,5] > 0){ 
        sDRIFT[DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = 
          tform(ukfstates[nlatent +DRIFTsetup[ri,5], statei ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(DIFFUSIONsetup)){
      if(DIFFUSIONsetup[ ri,5] > 0){ 
        sDIFFUSION[DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = 
          tform(ukfstates[nlatent +DIFFUSIONsetup[ri,5], statei ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(CINTsetup)){
      if(CINTsetup[ ri,5] > 0){ 
        sCINT[CINTsetup[ ri,1], CINTsetup[ri,2]] = 
          tform(ukfstates[nlatent +CINTsetup[ri,5], statei ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(TDPREDEFFECTsetup)){
      if(TDPREDEFFECTsetup[ ri,5] > 0){ 
        sTDPREDEFFECT[TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = 
          tform(ukfstates[nlatent +TDPREDEFFECTsetup[ri,5], statei ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ); 
      }
    }};
  
                if(statei <= (2+2*nlatentpop) ) {
                  rkstates[ki] = sDRIFT * state + sCINT[,1];
                } else if(statei <= (2+2*nlatentpop + ndynerror) ){
                  rkstates[ki] = sDRIFT * state + sCINT[,1] + sDIFFUSION[ , derrind[statei - (2+2*nlatentpop)] ] * dynerror; 
                } else rkstates[ki] = sDRIFT * state + sCINT[,1] - sDIFFUSION[ , derrind[statei - (2+2*nlatentpop + ndynerror)] ] * dynerror;
  
              }
              state = (rkstates[5] + dTsmall[rowi]/6  *(rkstates[1]+2*rkstates[2]+2*rkstates[3]+rkstates[4])); //integrate over rk steps
            }
          } //end nonlinear time update
    
    
          if(lineardynamics==1){ //this could be much more efficient...
              if(ukfpop==1){

    for(ri in 1:size(DRIFTsetup)){
      if(DRIFTsetup[ ri,5] > 0){ 
        sDRIFT[DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = 
          tform(ukfstates[nlatent +DRIFTsetup[ri,5], statei ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(DIFFUSIONsetup)){
      if(DIFFUSIONsetup[ ri,5] > 0){ 
        sDIFFUSION[DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = 
          tform(ukfstates[nlatent +DIFFUSIONsetup[ri,5], statei ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(CINTsetup)){
      if(CINTsetup[ ri,5] > 0){ 
        sCINT[CINTsetup[ ri,1], CINTsetup[ri,2]] = 
          tform(ukfstates[nlatent +CINTsetup[ri,5], statei ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(TDPREDEFFECTsetup)){
      if(TDPREDEFFECTsetup[ ri,5] > 0){ 
        sTDPREDEFFECT[TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = 
          tform(ukfstates[nlatent +TDPREDEFFECTsetup[ri,5], statei ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ); 
      }
    }};
            if(statei <= 2+2*nlatentpop+1){ //because after this its only noise variables
              if(continuoustime ==1){ 
                if(driftdiagonly==1) discreteDRIFT = matrix_diagexp(sDRIFT * dT[rowi]);
                if(driftdiagonly==0) discreteDRIFT = matrix_exp(sDRIFT * dT[rowi]);
                discreteCINT = sDRIFT \ (discreteDRIFT - IIlatent) * sCINT[,1];
                //discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
                  //quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]' );
                discreteDIFFUSION[derrind, derrind] = discreteDIFFUSIONcalc(sDRIFT[derrind, derrind], sDIFFUSION[derrind, derrind],dT[rowi]);;
              }
              if(continuoustime==0){
                discreteDRIFT=sDRIFT;
                discreteCINT=sCINT[,1];
                discreteDIFFUSION=sDIFFUSION;
              }
            }
            state = discreteDRIFT * state + discreteCINT;
            if(statei > (2+2*nlatentpop) && statei <= (2+2*nlatentpop+ndynerror) )  state += discreteDIFFUSION[ , derrind[statei - (2+2*nlatentpop)] ] * dynerror; 
            if(statei > (2+2*nlatentpop+ndynerror)) state += -discreteDIFFUSION[ , derrind[statei - (2+2*nlatentpop+ndynerror)] ] * dynerror;
          }
        }  // end of if not t0 section (time update)
    
          if(ukfpop==1){

    for(ri in 1:size(DRIFTsetup)){
      if(DRIFTsetup[ ri,5] > 0){ 
        sDRIFT[DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = 
          tform(ukfstates[nlatent +DRIFTsetup[ri,5], statei ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(DIFFUSIONsetup)){
      if(DIFFUSIONsetup[ ri,5] > 0){ 
        sDIFFUSION[DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = 
          tform(ukfstates[nlatent +DIFFUSIONsetup[ri,5], statei ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(CINTsetup)){
      if(CINTsetup[ ri,5] > 0){ 
        sCINT[CINTsetup[ ri,1], CINTsetup[ri,2]] = 
          tform(ukfstates[nlatent +CINTsetup[ri,5], statei ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(TDPREDEFFECTsetup)){
      if(TDPREDEFFECTsetup[ ri,5] > 0){ 
        sTDPREDEFFECT[TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = 
          tform(ukfstates[nlatent +TDPREDEFFECTsetup[ri,5], statei ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ); 
      }
    }};
        if(ntdpred > 0) state +=  (sTDPREDEFFECT * tdpreds[rowi]); //tdpred effect only influences at observed time point
        ukfstates[1:nlatent, statei] = state; //now contains time updated state
        if(statei==2) ukfstates[1:nlatent, 1] = state; //mean goes in twice for weighting
      }

      etaprior[rowi] = colMeans(ukfstates');
      etapriorcov[rowi] = cov_of_matrix(ukfstates') / asquared;
    } //end ukf time update

    if(intoverstates==0 && lineardynamics == 1) {
      if(T0check[rowi]==1) etaupd[rowi] = etaprior[rowi] +  sT0VAR * etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
      if(T0check[rowi]==0) etaupd[rowi] = etaprior[rowi] +  sDIFFUSION * etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
    }
  
    if(nobsi==0 && intoverstates==1) {
      etaupdcov[rowi] = etapriorcov[rowi];
      etaupd[rowi] = etaprior[rowi];
    }

    if (nobsi > 0) {  // if some observations create right size matrices for missingness and calculate...
  
      int cindex[intoverstates ? nobsi : ncont_y[rowi]];

      if(intoverstates==0) cindex = o0;
      if(intoverstates==1) cindex = o; //treat all obs as continuous gaussian

      if(ukf==0){ //non ukf measurement
        if(intoverstates==1) { //classic kalman
          ypred[o] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * etaprior[rowi];
          ypredcov[o,o] = quad_form(etapriorcov[rowi], sLAMBDA[o,]') + sMANIFESTVAR[o,o];
          for(wi in 1:nmanifest){ 
            if(manifesttype[wi]==1 && Y[rowi,wi] != 99999) ypredcov[wi,wi] = ypredcov[wi,wi] + fabs((ypred[wi] - 1) .* (ypred[wi])); 
          }
          K[,o] = mdivide_right(etapriorcov[rowi] * sLAMBDA[o,]', ypredcov[o,o]); 
          etaupdcov[rowi] = (IIlatent - K[,o] * sLAMBDA[o,]) * etapriorcov[rowi];
        }
        if(intoverstates==0) { //sampled states
          if(ncont_y[rowi] > 0) ypred[o0] = sMANIFESTMEANS[o0,1] + sLAMBDA[o0,] * etaupd[rowi];
          if(nbinary_y[rowi] > 0) ypred[o1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * etaupd[rowi])));
          ypredcov[o,o] = sMANIFESTVAR[o,o];
        }
      }
  

      if(ukf==1){ //ukf measurement
        for(statei in 2:cols(ukfmeasures)){
          state = ukfstates[ 1:nlatent, statei];
          
      if(ukfpop==1){
        
    for(ri in 1:size(LAMBDAsetup)){
      if(LAMBDAsetup[ ri,5] > 0){ 
        sLAMBDA[LAMBDAsetup[ ri,1], LAMBDAsetup[ri,2]] = tform(ukfstates[nlatent +LAMBDAsetup[ri,5], statei ], LAMBDAsetup[ri,4], LAMBDAvalues[ri,2], LAMBDAvalues[ri,3], LAMBDAvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(MANIFESTMEANSsetup)){
      if(MANIFESTMEANSsetup[ ri,5] > 0){ 
        sMANIFESTMEANS[MANIFESTMEANSsetup[ ri,1], MANIFESTMEANSsetup[ri,2]] = tform(ukfstates[nlatent +MANIFESTMEANSsetup[ri,5], statei ], MANIFESTMEANSsetup[ri,4], MANIFESTMEANSvalues[ri,2], MANIFESTMEANSvalues[ri,3], MANIFESTMEANSvalues[ri,4] ); 
      }
    }

    for(ri in 1:size(MANIFESTVARsetup)){
      if(MANIFESTVARsetup[ ri,5] > 0){ 
        sMANIFESTVAR[MANIFESTVARsetup[ ri,1], MANIFESTVARsetup[ri,2]] = tform(ukfstates[nlatent +MANIFESTVARsetup[ri,5], statei ], MANIFESTVARsetup[ri,4], MANIFESTVARvalues[ri,2], MANIFESTVARvalues[ri,3], MANIFESTVARvalues[ri,4] ); 
      }
    }};

  
          if(ncont_y[rowi] > 0) ukfmeasures[o0 , statei] = sMANIFESTMEANS[o0,1] + sLAMBDA[o0,] * state;
          if(nbinary_y[rowi] > 0) {
            ukfmeasures[o1 , statei] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * state)));
          }
          if(statei==2) { //temporary measure to get mean in twice
          if(ncont_y[rowi] > 0) ukfmeasures[o0 , 1] = sMANIFESTMEANS[o0,1] + sLAMBDA[o0,] * state;
          if(nbinary_y[rowi] > 0) {
            ukfmeasures[o1 , 1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * state)));
          }
          }
        } //end ukf measurement loop
    
        ypred[o] = colMeans(ukfmeasures[o,]'); 
        ypredcov[o,o] = cov_of_matrix(ukfmeasures[o,]') /asquared + sMANIFESTVAR[o,o];
        for(wi in 1:nmanifest){ 
          if(manifesttype[wi]==1 && Y[rowi,wi] != 99999) ypredcov[wi,wi] = ypredcov[wi,wi] + fabs((ypred[wi] - 1) .* (ypred[wi])); //sMANIFESTVAR[wi,wi] + (merror[wi] / cols(ukfmeasures) +1e-8);
        }
        K[,o] = mdivide_right(crosscov(ukfstates', ukfmeasures[o,]') /asquared, ypredcov[o,o]); 
        etaupdcov[rowi] = etapriorcov[rowi] - quad_form(ypredcov[o,o],  K[,o]');
      }

      
if(verbose > 1) {
print("rowi ",rowi, "  si ", si, "  etaprior[rowi] ",etaprior[rowi],"  etapriorcov[rowi] ",etapriorcov[rowi],
          "  etaupd[rowi] ",etaupd[rowi],"  etaupdcov[rowi] ",etaupdcov[rowi],"  ypred ",ypred,"  ypredcov ",ypredcov, "  K ",K,
          "  sDRIFT ", sDRIFT, " sDIFFUSION ", sDIFFUSION, " sCINT ", sCINT, "  sMANIFESTVAR ", sMANIFESTVAR, "  rawpopsd ", rawpopsd,
          "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        if(lineardynamics==1) print("discreteDRIFT ",discreteDRIFT,"  discreteCINT ", discreteCINT, "  discreteDIFFUSION ", discreteDIFFUSION)
}
if(verbose > 2) print("ukfstates ", ukfstates, "  ukfmeasures ", ukfmeasures);

          ypredcov_sqrt[cindex,cindex]=cholspd(ypredcov[o0,o0]); //use o0, or cindex?
          if(ncont_y[rowi] > 0) Ygen[geni, rowi, o0] = multi_normal_cholesky_rng(ypred[o], ypredcov_sqrt[o0,o0]);
          if(nbinary_y[rowi] > 0) for(obsi in 1:size(o1)) Ygen[geni, rowi, o1[obsi]] = bernoulli_rng(ypred[o1[obsi]]);
      
  
      err[o] = Y[rowi,o] - ypred[o]; // prediction error
      if(intoverstates==1) etaupd[rowi,] = etaprior[rowi,] + (K[,o] * err[o]);

  
      
    }//end nobs > 0 section
  }//end rowi

  
}

{
  vector[nparams] rawindparams;
  rawindparams = rawpopmeans;
  for(si in 1:1){

    if(ntipred==0 && nindvarying > 0 && ukfpop ==0) rawindparams[indvaryingindex] = 
      rawpopmeans[indvaryingindex] + rawpopcovsqrt * baseindparams[(1+(si-1)*nindvarying):(si*nindvarying)];

    if(ntipred > 0  && nindvarying > 0 && ukfpop ==0) rawindparams[indvaryingindex] = 
      rawpopmeans[indvaryingindex] + rawpopcovsqrt * baseindparams[(1+(si-1)*nindvarying):(si*nindvarying)] +
      TIPREDEFFECT[indvaryingindex,] * tipreds[si]';

  if(si <= T0MEANSsubindex[1]){
    for(ri in 1:size(T0MEANSsetup)){
      pop_T0MEANS[si, T0MEANSsetup[ ri,1], T0MEANSsetup[ri,2]] = T0MEANSsetup[ri,3] ? tform(rawindparams[ T0MEANSsetup[ri,3] ], T0MEANSsetup[ri,4], T0MEANSvalues[ri,2], T0MEANSvalues[ri,3], T0MEANSvalues[ri,4] ) : T0MEANSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= LAMBDAsubindex[1]){
    for(ri in 1:size(LAMBDAsetup)){
      pop_LAMBDA[si, LAMBDAsetup[ ri,1], LAMBDAsetup[ri,2]] = LAMBDAsetup[ri,3] ? tform(rawindparams[ LAMBDAsetup[ri,3] ], LAMBDAsetup[ri,4], LAMBDAvalues[ri,2], LAMBDAvalues[ri,3], LAMBDAvalues[ri,4] ) : LAMBDAvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= DRIFTsubindex[1]){
    for(ri in 1:size(DRIFTsetup)){
      pop_DRIFT[si, DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = DRIFTsetup[ri,3] ? tform(rawindparams[ DRIFTsetup[ri,3] ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ) : DRIFTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= DIFFUSIONsubindex[1]){
    for(ri in 1:size(DIFFUSIONsetup)){
      pop_DIFFUSION[si, DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = DIFFUSIONsetup[ri,3] ? tform(rawindparams[ DIFFUSIONsetup[ri,3] ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ) : DIFFUSIONvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= MANIFESTVARsubindex[1]){
    for(ri in 1:size(MANIFESTVARsetup)){
      pop_MANIFESTVAR[si, MANIFESTVARsetup[ ri,1], MANIFESTVARsetup[ri,2]] = MANIFESTVARsetup[ri,3] ? tform(rawindparams[ MANIFESTVARsetup[ri,3] ], MANIFESTVARsetup[ri,4], MANIFESTVARvalues[ri,2], MANIFESTVARvalues[ri,3], MANIFESTVARvalues[ri,4] ) : MANIFESTVARvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= MANIFESTMEANSsubindex[1]){
    for(ri in 1:size(MANIFESTMEANSsetup)){
      pop_MANIFESTMEANS[si, MANIFESTMEANSsetup[ ri,1], MANIFESTMEANSsetup[ri,2]] = MANIFESTMEANSsetup[ri,3] ? tform(rawindparams[ MANIFESTMEANSsetup[ri,3] ], MANIFESTMEANSsetup[ri,4], MANIFESTMEANSvalues[ri,2], MANIFESTMEANSvalues[ri,3], MANIFESTMEANSvalues[ri,4] ) : MANIFESTMEANSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= CINTsubindex[1]){
    for(ri in 1:size(CINTsetup)){
      pop_CINT[si, CINTsetup[ ri,1], CINTsetup[ri,2]] = CINTsetup[ri,3] ? tform(rawindparams[ CINTsetup[ri,3] ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ) : CINTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= T0VARsubindex[1]){
    for(ri in 1:size(T0VARsetup)){
      pop_T0VAR[si, T0VARsetup[ ri,1], T0VARsetup[ri,2]] = T0VARsetup[ri,3] ? tform(rawindparams[ T0VARsetup[ri,3] ], T0VARsetup[ri,4], T0VARvalues[ri,2], T0VARvalues[ri,3], T0VARvalues[ri,4] ) : T0VARvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      

  if(si <= TDPREDEFFECTsubindex[1]){
    for(ri in 1:size(TDPREDEFFECTsetup)){
      pop_TDPREDEFFECT[si, TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = TDPREDEFFECTsetup[ri,3] ? tform(rawindparams[ TDPREDEFFECTsetup[ri,3] ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ) : TDPREDEFFECTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      


  // perform any whole matrix transformations 
    
  if(si <= DIFFUSIONsubindex[1]) pop_DIFFUSION[si] = sdcovsqrt2cov(pop_DIFFUSION[si], lineardynamics * intoverstates ? 0 : 1);

  if(lineardynamics==1 && ndiffusion > 0){
    if(si <= asymDIFFUSIONsubindex[1]) {
      if(ndiffusion < nlatent) asympop_DIFFUSION[si] = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

      if(continuoustime==1) asympop_DIFFUSION[si, derrind, derrind] = to_matrix( 
      -( kron_prod( pop_DRIFT[DRIFTsubindex[si], derrind, derrind ], IIlatent[ derrind, derrind ]) + 
         kron_prod(IIlatent[ derrind, derrind ], pop_DRIFT[ DRIFTsubindex[si], derrind, derrind ]) ) \ 
      to_vector( pop_DIFFUSION[ DIFFUSIONsubindex[si], derrind, derrind ]), ndiffusion,ndiffusion);

      if(continuoustime==0) asympop_DIFFUSION[si, derrind, derrind] = to_matrix( (IIlatent2[ derrind, derrind ] - 
        kron_prod(pop_DRIFT[ DRIFTsubindex[si], derrind, derrind  ], 
          pop_DRIFT[ DRIFTsubindex[si], derrind, derrind  ])) * 
        to_vector(pop_DIFFUSION[ DIFFUSIONsubindex[si], derrind, derrind  ]) , ndiffusion, ndiffusion);
    } //end asymdiffusion loops
  }
          
    if(nt0meansstationary > 0){
      if(si <= asymCINTsubindex[1]){
        if(continuoustime==1) asympop_CINT[si] =  -pop_DRIFT[ DRIFTsubindex[si] ] \ pop_CINT[ CINTsubindex[si], ,1 ];
        if(continuoustime==0) asympop_CINT[si] =  (IIlatent - pop_DRIFT[ DRIFTsubindex[si] ]) \ pop_CINT[ CINTsubindex[si], ,1 ];
      }
    }

          
    if(binomial==0){
      if(si <= MANIFESTVARsubindex[1]) {
         for(ri in 1:nmanifest) pop_MANIFESTVAR[si,ri,ri] = square(pop_MANIFESTVAR[si,ri,ri]);
      }
    }
          
          
    if(si <= T0VARsubindex[1]) {
      pop_T0VAR[si] = sdcovsqrt2cov(pop_T0VAR[si],lineardynamics * intoverstates ? 0 : 1);

      if(nt0varstationary > 0) for(rowi in 1:nt0varstationary){
        pop_T0VAR[si,t0varstationary[rowi,1],t0varstationary[rowi,2] ] = 
          asympop_DIFFUSION[si,t0varstationary[rowi,1],t0varstationary[rowi,2] ];
        pop_T0VAR[si,t0varstationary[rowi,2],t0varstationary[rowi,1] ] = 
          asympop_DIFFUSION[si,t0varstationary[rowi,2],t0varstationary[rowi,1] ];
      }
    }

    
    if(nt0meansstationary > 0){
      if(si <= T0MEANSsubindex[1]) {
        for(rowi in 1:nt0meansstationary){
          pop_T0MEANS[si,t0meansstationary[rowi,1] , 1] = 
            asympop_CINT[ asymCINTsubindex[si], t0meansstationary[rowi,1] ];
        }
      }
    }
  }
} //end subject loop
  

rawpopcorr = tcrossprod(rawpopcorrsqrt);
rawpopcov = tcrossprod(rawpopcovsqrt);

popsd = rep_vector(0,nparams);
popsd[indvaryingindex] = rawpopsd; //base to begin calculations
 
    for(ri in 1:size(T0MEANSsetup)){
      if(T0MEANSsetup[ri,3] !=0) {

        popmeans[T0MEANSsetup[ ri,3]] = tform(rawpopmeans[T0MEANSsetup[ri,3] ], T0MEANSsetup[ri,4], T0MEANSvalues[ri,2], T0MEANSvalues[ri,3], T0MEANSvalues[ri,4] ); 

        popsd[T0MEANSsetup[ ri,3]] = T0MEANSsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[T0MEANSsetup[ri,3] ], T0MEANSsetup[ri,4], T0MEANSvalues[ri,2], T0MEANSvalues[ri,3], T0MEANSvalues[ri,4] + popsd[T0MEANSsetup[ ri,3]]) -
           tform(
            rawpopmeans[T0MEANSsetup[ri,3] ], T0MEANSsetup[ri,4], T0MEANSvalues[ri,2], T0MEANSvalues[ri,3], T0MEANSvalues[ri,4] - popsd[T0MEANSsetup[ ri,3]]) / 2) : 0; 
      }
    }
      

    for(ri in 1:size(LAMBDAsetup)){
      if(LAMBDAsetup[ri,3] !=0) {

        popmeans[LAMBDAsetup[ ri,3]] = tform(rawpopmeans[LAMBDAsetup[ri,3] ], LAMBDAsetup[ri,4], LAMBDAvalues[ri,2], LAMBDAvalues[ri,3], LAMBDAvalues[ri,4] ); 

        popsd[LAMBDAsetup[ ri,3]] = LAMBDAsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[LAMBDAsetup[ri,3] ], LAMBDAsetup[ri,4], LAMBDAvalues[ri,2], LAMBDAvalues[ri,3], LAMBDAvalues[ri,4] + popsd[LAMBDAsetup[ ri,3]]) -
           tform(
            rawpopmeans[LAMBDAsetup[ri,3] ], LAMBDAsetup[ri,4], LAMBDAvalues[ri,2], LAMBDAvalues[ri,3], LAMBDAvalues[ri,4] - popsd[LAMBDAsetup[ ri,3]]) / 2) : 0; 
      }
    }
      

    for(ri in 1:size(DRIFTsetup)){
      if(DRIFTsetup[ri,3] !=0) {

        popmeans[DRIFTsetup[ ri,3]] = tform(rawpopmeans[DRIFTsetup[ri,3] ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ); 

        popsd[DRIFTsetup[ ri,3]] = DRIFTsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[DRIFTsetup[ri,3] ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] + popsd[DRIFTsetup[ ri,3]]) -
           tform(
            rawpopmeans[DRIFTsetup[ri,3] ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] - popsd[DRIFTsetup[ ri,3]]) / 2) : 0; 
      }
    }
      

    for(ri in 1:size(DIFFUSIONsetup)){
      if(DIFFUSIONsetup[ri,3] !=0) {

        popmeans[DIFFUSIONsetup[ ri,3]] = tform(rawpopmeans[DIFFUSIONsetup[ri,3] ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ); 

        popsd[DIFFUSIONsetup[ ri,3]] = DIFFUSIONsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[DIFFUSIONsetup[ri,3] ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] + popsd[DIFFUSIONsetup[ ri,3]]) -
           tform(
            rawpopmeans[DIFFUSIONsetup[ri,3] ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] - popsd[DIFFUSIONsetup[ ri,3]]) / 2) : 0; 
      }
    }
      

    for(ri in 1:size(MANIFESTVARsetup)){
      if(MANIFESTVARsetup[ri,3] !=0) {

        popmeans[MANIFESTVARsetup[ ri,3]] = tform(rawpopmeans[MANIFESTVARsetup[ri,3] ], MANIFESTVARsetup[ri,4], MANIFESTVARvalues[ri,2], MANIFESTVARvalues[ri,3], MANIFESTVARvalues[ri,4] ); 

        popsd[MANIFESTVARsetup[ ri,3]] = MANIFESTVARsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[MANIFESTVARsetup[ri,3] ], MANIFESTVARsetup[ri,4], MANIFESTVARvalues[ri,2], MANIFESTVARvalues[ri,3], MANIFESTVARvalues[ri,4] + popsd[MANIFESTVARsetup[ ri,3]]) -
           tform(
            rawpopmeans[MANIFESTVARsetup[ri,3] ], MANIFESTVARsetup[ri,4], MANIFESTVARvalues[ri,2], MANIFESTVARvalues[ri,3], MANIFESTVARvalues[ri,4] - popsd[MANIFESTVARsetup[ ri,3]]) / 2) : 0; 
      }
    }
      

    for(ri in 1:size(MANIFESTMEANSsetup)){
      if(MANIFESTMEANSsetup[ri,3] !=0) {

        popmeans[MANIFESTMEANSsetup[ ri,3]] = tform(rawpopmeans[MANIFESTMEANSsetup[ri,3] ], MANIFESTMEANSsetup[ri,4], MANIFESTMEANSvalues[ri,2], MANIFESTMEANSvalues[ri,3], MANIFESTMEANSvalues[ri,4] ); 

        popsd[MANIFESTMEANSsetup[ ri,3]] = MANIFESTMEANSsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[MANIFESTMEANSsetup[ri,3] ], MANIFESTMEANSsetup[ri,4], MANIFESTMEANSvalues[ri,2], MANIFESTMEANSvalues[ri,3], MANIFESTMEANSvalues[ri,4] + popsd[MANIFESTMEANSsetup[ ri,3]]) -
           tform(
            rawpopmeans[MANIFESTMEANSsetup[ri,3] ], MANIFESTMEANSsetup[ri,4], MANIFESTMEANSvalues[ri,2], MANIFESTMEANSvalues[ri,3], MANIFESTMEANSvalues[ri,4] - popsd[MANIFESTMEANSsetup[ ri,3]]) / 2) : 0; 
      }
    }
      

    for(ri in 1:size(CINTsetup)){
      if(CINTsetup[ri,3] !=0) {

        popmeans[CINTsetup[ ri,3]] = tform(rawpopmeans[CINTsetup[ri,3] ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ); 

        popsd[CINTsetup[ ri,3]] = CINTsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[CINTsetup[ri,3] ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] + popsd[CINTsetup[ ri,3]]) -
           tform(
            rawpopmeans[CINTsetup[ri,3] ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] - popsd[CINTsetup[ ri,3]]) / 2) : 0; 
      }
    }
      

    for(ri in 1:size(T0VARsetup)){
      if(T0VARsetup[ri,3] !=0) {

        popmeans[T0VARsetup[ ri,3]] = tform(rawpopmeans[T0VARsetup[ri,3] ], T0VARsetup[ri,4], T0VARvalues[ri,2], T0VARvalues[ri,3], T0VARvalues[ri,4] ); 

        popsd[T0VARsetup[ ri,3]] = T0VARsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[T0VARsetup[ri,3] ], T0VARsetup[ri,4], T0VARvalues[ri,2], T0VARvalues[ri,3], T0VARvalues[ri,4] + popsd[T0VARsetup[ ri,3]]) -
           tform(
            rawpopmeans[T0VARsetup[ri,3] ], T0VARsetup[ri,4], T0VARvalues[ri,2], T0VARvalues[ri,3], T0VARvalues[ri,4] - popsd[T0VARsetup[ ri,3]]) / 2) : 0; 
      }
    }
      

    for(ri in 1:size(TDPREDEFFECTsetup)){
      if(TDPREDEFFECTsetup[ri,3] !=0) {

        popmeans[TDPREDEFFECTsetup[ ri,3]] = tform(rawpopmeans[TDPREDEFFECTsetup[ri,3] ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ); 

        popsd[TDPREDEFFECTsetup[ ri,3]] = TDPREDEFFECTsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[TDPREDEFFECTsetup[ri,3] ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] + popsd[TDPREDEFFECTsetup[ ri,3]]) -
           tform(
            rawpopmeans[TDPREDEFFECTsetup[ri,3] ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] - popsd[TDPREDEFFECTsetup[ ri,3]]) / 2) : 0; 
      }
    }
      

}
