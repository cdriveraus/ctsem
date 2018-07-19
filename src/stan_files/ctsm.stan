
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
        l[i,j] = (a[i,j] + a[j,i])/2;
        l[j,i] = l[i,j];
      }
      if(j == i){
        l[j,j] = a[j,j] + 1e-3;
        if(l[j,j] <=1e-3) l[j,j] = 1e-3; 
        if(l[j,j] > 1e10) l[j,j] = -99999;
      }
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




  matrix sdcovsqrt2cov(matrix mat, int msqrt){ //converts from lower partial sd and diag sd to cov or cholesky cov
    matrix[rows(mat),rows(mat)] out;

    if(msqrt==1){
      for(k in 1:cols(mat)){
        for(j in 1:rows(mat)){
          if(j > k) out[j,k] = mat[j,k];
          if(k > j) out[j,k] = mat[k,j];
          if(k==j) out[j,k] = mat[j,k];
        }
      }
    }
  
    if(msqrt==0){
      out = covsqrt2corsqrt(mat, 0);
      out = diag_pre_multiply(diagonal(mat), out);
    }

    if(msqrt==0) out = tcrossprod(out);
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
      for (ri in 1:rows(mat))  {
        centered[ri,coli] = mat[ri,coli] - means[coli];
      }
    }
    covm = crossprod(centered) / (rows(mat)-1);
    for(j in 1:rows(covm)){
      covm[j,j] = covm[j,j] + 1e-6;
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

  matrix makesym(matrix mat){
    matrix[rows(mat),cols(mat)] out;

    for(coli in 1:cols(mat)){
      for(rowi in coli:rows(mat)){
        if(rowi > coli) {
          out[rowi,coli] = mat[rowi,coli]; //(mat[coli,rowi] + ) *.5;
          out[coli,rowi] = mat[rowi,coli];
        }
        if(rowi==coli) out[rowi,coli] = mat[rowi,coli] + 1e-5;
      }
    }
    return out;
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
  real<lower=0> tipredsimputedscale;
  real<lower=0> tipredeffectscale;
  
  
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
  int<lower = 0, upper = nmanifest> whichcont_y[ndatapoints, nmanifest]; // index of which variables are observed and continuous per observation
  
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
int PARSsubindex[nsubjects];
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
int PARSsetup_rowcount;
  int T0MEANSsetup[T0MEANSsetup_rowcount,6 ];
int LAMBDAsetup[LAMBDAsetup_rowcount,6 ];
int DRIFTsetup[DRIFTsetup_rowcount,6 ];
int DIFFUSIONsetup[DIFFUSIONsetup_rowcount,6 ];
int MANIFESTVARsetup[MANIFESTVARsetup_rowcount,6 ];
int MANIFESTMEANSsetup[MANIFESTMEANSsetup_rowcount,6 ];
int CINTsetup[CINTsetup_rowcount,6 ];
int T0VARsetup[T0VARsetup_rowcount,6 ];
int TDPREDEFFECTsetup[TDPREDEFFECTsetup_rowcount,6 ];
int PARSsetup[PARSsetup_rowcount,6 ];
  matrix[T0MEANSsetup_rowcount, 5] T0MEANSvalues;
matrix[LAMBDAsetup_rowcount, 5] LAMBDAvalues;
matrix[DRIFTsetup_rowcount, 5] DRIFTvalues;
matrix[DIFFUSIONsetup_rowcount, 5] DIFFUSIONvalues;
matrix[MANIFESTVARsetup_rowcount, 5] MANIFESTVARvalues;
matrix[MANIFESTMEANSsetup_rowcount, 5] MANIFESTMEANSvalues;
matrix[CINTsetup_rowcount, 5] CINTvalues;
matrix[T0VARsetup_rowcount, 5] T0VARvalues;
matrix[TDPREDEFFECTsetup_rowcount, 5] TDPREDEFFECTvalues;
matrix[PARSsetup_rowcount, 5] PARSvalues;
  int TIPREDEFFECTsetup[nparams, ntipred];
  int nmatrixslots;
  int popsetup[nmatrixslots,6];
  real popvalues[nmatrixslots,5];
  int savescores;
}
      
transformed data{
  matrix[nlatent,nlatent] IIlatent;
  matrix[nlatent*nlatent,nlatent*nlatent] IIlatent2;
  int nlatentpop;

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
  real ll;
  vector[savescores ? nlatentpop : 0] etaprior_out[savescores ? ndatapoints : 0];
  vector[savescores ? nlatentpop : 0] etaupd_out[savescores ? ndatapoints : 0];

  
  
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
    rawpopsd = log(1+exp(2*rawpopsdbase)) .* sdscale;
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
  vector[nt0meansstationary ? nlatent : 0] sasymCINT;
  matrix[nmanifest,nmanifest] sMANIFESTVAR; 
  matrix[nmanifest,1] sMANIFESTMEANS;
  matrix[nmanifest,nlatent] sLAMBDA;
  matrix[ntdpred ? nlatent : 0,ntdpred] sTDPREDEFFECT;
  matrix[PARSsetup_rowcount ? max(PARSsetup[,1]) : 0 ,PARSsetup_rowcount ? max(PARSsetup[,2]) : 0] sPARS;

  //ukf approximation parameters
  if(ukf==1) k = 0.5;

  if(lineardynamics) discreteDIFFUSION = rep_matrix(0,nlatent,nlatent); //in case some elements remain zero due to derrind
  
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
        asquared =  2.0/sqrt(0.0+nlatentpop+ndynerror) * 1e-1;
        l = asquared * (nlatentpop + ndynerror + k) - (nlatentpop+ndynerror); 
        sqrtukfadjust = sqrt(0.0+nlatentpop + ndynerror +l);
      }
    }

    if(T0check[rowi] == 1) { // calculate initial matrices if this is first row for si

    {
  vector[nparams] rawindparams;
  vector[nparams] tipredaddition;
  vector[nparams] indvaraddition;
  
  if(si==1 || (si > 1 && (nindvarying >0 || ntipred > 0))){
    tipredaddition = rep_vector(0,nparams);
    indvaraddition = rep_vector(0,nparams);

    if(nindvarying > 0 && ukfpop==0) indvaraddition[indvaryingindex] = rawpopcovsqrt * baseindparams[(1+(si-1)*nindvarying):(si*nindvarying)];
  
    if(ntipred > 0) tipredaddition = TIPREDEFFECT * tipreds[si]';
  
    rawindparams = rawpopmeans + tipredaddition + indvaraddition;
  }

  if(si <= T0MEANSsubindex[nsubjects]){
    for(ri in 1:size(T0MEANSsetup)){
      if(si==1 || T0MEANSsetup[ri,5] > 0 || T0MEANSsetup[ri,6] > 0){
        sT0MEANS[T0MEANSsetup[ ri,1], T0MEANSsetup[ri,2]] = T0MEANSsetup[ri,3] ? tform(rawindparams[ T0MEANSsetup[ri,3] ], T0MEANSsetup[ri,4], T0MEANSvalues[ri,2], T0MEANSvalues[ri,3], T0MEANSvalues[ri,4] ) : T0MEANSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= LAMBDAsubindex[nsubjects]){
    for(ri in 1:size(LAMBDAsetup)){
      if(si==1 || LAMBDAsetup[ri,5] > 0 || LAMBDAsetup[ri,6] > 0){
        sLAMBDA[LAMBDAsetup[ ri,1], LAMBDAsetup[ri,2]] = LAMBDAsetup[ri,3] ? tform(rawindparams[ LAMBDAsetup[ri,3] ], LAMBDAsetup[ri,4], LAMBDAvalues[ri,2], LAMBDAvalues[ri,3], LAMBDAvalues[ri,4] ) : LAMBDAvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= DRIFTsubindex[nsubjects]){
    for(ri in 1:size(DRIFTsetup)){
      if(si==1 || DRIFTsetup[ri,5] > 0 || DRIFTsetup[ri,6] > 0){
        sDRIFT[DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = DRIFTsetup[ri,3] ? tform(rawindparams[ DRIFTsetup[ri,3] ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ) : DRIFTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= DIFFUSIONsubindex[nsubjects]){
    for(ri in 1:size(DIFFUSIONsetup)){
      if(si==1 || DIFFUSIONsetup[ri,5] > 0 || DIFFUSIONsetup[ri,6] > 0){
        sDIFFUSION[DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = DIFFUSIONsetup[ri,3] ? tform(rawindparams[ DIFFUSIONsetup[ri,3] ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ) : DIFFUSIONvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= MANIFESTVARsubindex[nsubjects]){
    for(ri in 1:size(MANIFESTVARsetup)){
      if(si==1 || MANIFESTVARsetup[ri,5] > 0 || MANIFESTVARsetup[ri,6] > 0){
        sMANIFESTVAR[MANIFESTVARsetup[ ri,1], MANIFESTVARsetup[ri,2]] = MANIFESTVARsetup[ri,3] ? tform(rawindparams[ MANIFESTVARsetup[ri,3] ], MANIFESTVARsetup[ri,4], MANIFESTVARvalues[ri,2], MANIFESTVARvalues[ri,3], MANIFESTVARvalues[ri,4] ) : MANIFESTVARvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= MANIFESTMEANSsubindex[nsubjects]){
    for(ri in 1:size(MANIFESTMEANSsetup)){
      if(si==1 || MANIFESTMEANSsetup[ri,5] > 0 || MANIFESTMEANSsetup[ri,6] > 0){
        sMANIFESTMEANS[MANIFESTMEANSsetup[ ri,1], MANIFESTMEANSsetup[ri,2]] = MANIFESTMEANSsetup[ri,3] ? tform(rawindparams[ MANIFESTMEANSsetup[ri,3] ], MANIFESTMEANSsetup[ri,4], MANIFESTMEANSvalues[ri,2], MANIFESTMEANSvalues[ri,3], MANIFESTMEANSvalues[ri,4] ) : MANIFESTMEANSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= CINTsubindex[nsubjects]){
    for(ri in 1:size(CINTsetup)){
      if(si==1 || CINTsetup[ri,5] > 0 || CINTsetup[ri,6] > 0){
        sCINT[CINTsetup[ ri,1], CINTsetup[ri,2]] = CINTsetup[ri,3] ? tform(rawindparams[ CINTsetup[ri,3] ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ) : CINTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= T0VARsubindex[nsubjects]){
    for(ri in 1:size(T0VARsetup)){
      if(si==1 || T0VARsetup[ri,5] > 0 || T0VARsetup[ri,6] > 0){
        sT0VAR[T0VARsetup[ ri,1], T0VARsetup[ri,2]] = T0VARsetup[ri,3] ? tform(rawindparams[ T0VARsetup[ri,3] ], T0VARsetup[ri,4], T0VARvalues[ri,2], T0VARvalues[ri,3], T0VARvalues[ri,4] ) : T0VARvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= TDPREDEFFECTsubindex[nsubjects]){
    for(ri in 1:size(TDPREDEFFECTsetup)){
      if(si==1 || TDPREDEFFECTsetup[ri,5] > 0 || TDPREDEFFECTsetup[ri,6] > 0){
        sTDPREDEFFECT[TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = TDPREDEFFECTsetup[ri,3] ? tform(rawindparams[ TDPREDEFFECTsetup[ri,3] ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ) : TDPREDEFFECTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= PARSsubindex[nsubjects]){
    for(ri in 1:size(PARSsetup)){
      if(si==1 || PARSsetup[ri,5] > 0 || PARSsetup[ri,6] > 0){
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = PARSsetup[ri,3] ? tform(rawindparams[ PARSsetup[ri,3] ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ) : PARSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      


  // perform any whole matrix transformations 
    
  if(si <= DIFFUSIONsubindex[nsubjects] && lineardynamics) sDIFFUSION = sdcovsqrt2cov(sDIFFUSION, 0);

    if(si <= asymDIFFUSIONsubindex[nsubjects]) {
      if(ndiffusion < nlatent) sasymDIFFUSION = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

      if(continuoustime==1) sasymDIFFUSION[ derrind, derrind] = to_matrix( 
      -( kron_prod( sDRIFT[ derrind, derrind ], IIlatent[ derrind, derrind ]) + 
         kron_prod(IIlatent[ derrind, derrind ], sDRIFT[ derrind, derrind ]) ) \ 
      to_vector( sDIFFUSION[ derrind, derrind ] + IIlatent[ derrind, derrind ] * 1e-5), ndiffusion,ndiffusion);

      if(continuoustime==0) sasymDIFFUSION[derrind, derrind] = to_matrix( (IIlatent2[ derrind, derrind ] - 
        kron_prod(sDRIFT[derrind, derrind  ], 
          sDRIFT[derrind, derrind  ])) * 
        to_vector(sDIFFUSION[derrind, derrind  ]) , ndiffusion, ndiffusion);
    } //end asymdiffusion loops
          
    if(nt0meansstationary > 0){
      if(si <= asymCINTsubindex[nsubjects]){
        if(continuoustime==1) sasymCINT =  -sDRIFT \ sCINT[ ,1 ];
        if(continuoustime==0) sasymCINT =  (IIlatent - sDRIFT) \ sCINT[,1 ];
      }
    }

          
    if(binomial==0){
      if(si <= MANIFESTVARsubindex[nsubjects]) {
         for(ri in 1:nmanifest) sMANIFESTVAR[ri,ri] = square(sMANIFESTVAR[ri,ri]);
      }
    }
          
          
    if(si <= T0VARsubindex[nsubjects]) {
      if(lineardynamics * intoverstates !=0) sT0VAR = sdcovsqrt2cov(sT0VAR,0);
      if(nt0varstationary > 0) for(ri in 1:nt0varstationary){
        sT0VAR[t0varstationary[ri,1],t0varstationary[ri,2] ] = 
          sasymDIFFUSION[t0varstationary[ri,1],t0varstationary[ri,2] ];
        sT0VAR[t0varstationary[ri,2],t0varstationary[ri,1] ] = 
          sasymDIFFUSION[t0varstationary[ri,2],t0varstationary[ri,1] ];
      }
    }

    
    if(nt0meansstationary > 0){
      if(si <= T0MEANSsubindex[nsubjects]) {
        for(ri in 1:nt0meansstationary){
          sT0MEANS[t0meansstationary[ri,1] , 1] = 
            sasymCINT[t0meansstationary[ri,1] ];
        }
      }
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
      if(ntdpred > 0) etaprior[rowi] = sTDPREDEFFECT * tdpreds[rowi] + sT0MEANS[,1];
      etapriorcov[rowi] =  sT0VAR;
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
          discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
            quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]' );
          //discreteDIFFUSION[derrind, derrind] = discreteDIFFUSIONcalc(DRIFT[ DRIFTsubindex[si], derrind, derrind], sDIFFUSION[derrind, derrind], dT[rowi]);
          if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(discreteDIFFUSION);
        }
      }
  
      if(continuoustime==0 && T0check[rowi-1] == 1){
        discreteDRIFT=sDRIFT;
        discreteCINT=sCINT[,1];
        discreteDIFFUSION=sDIFFUSION;
        if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(discreteDIFFUSION);
      }

      if(ntdpred == 0) etaprior[rowi] = discreteDRIFT * etaupd[rowi-1] + discreteCINT;
      if(ntdpred > 0) etaprior[rowi] = discreteDRIFT * etaupd[rowi-1] + discreteCINT + sTDPREDEFFECT * tdpreds[rowi];
      if(intoverstates==1) {
        if(ndiffusion ==0) etapriorcov[rowi] = quad_form(etaupdcov[rowi-1], discreteDRIFT');
        if(ndiffusion > 0) etapriorcov[rowi] = quad_form(etaupdcov[rowi-1], discreteDRIFT') + discreteDIFFUSION;
      }
    }//end linear time update

    if(ukf==1){ //ukf time update

      if(T0check[rowi]==1) dynerror = sqrtukfadjust;
      if(T0check[rowi]==0 && lineardynamics==0) dynerror = sqrtukfadjust / sqrt(dT[rowi]); //Weiner process variance adjustment
  
      if(T0check[rowi]==0){ //compute updated sigpoints
      sigpoints = cholesky_decompose(makesym(etaupdcov[rowi-1,,] * sqrtukfadjust));
      
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
          ;
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
                ;
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
    }

    for(ri in 1:size(PARSsetup)){
      if(PARSsetup[ ri,5] > 0){ 
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = 
          tform(ukfstates[nlatent +PARSsetup[ri,5], statei ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ); 
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
            ;
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
    }

    for(ri in 1:size(PARSsetup)){
      if(PARSsetup[ ri,5] > 0){ 
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = 
          tform(ukfstates[nlatent +PARSsetup[ri,5], statei ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ); 
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
    
        ;
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
    }

    for(ri in 1:size(PARSsetup)){
      if(PARSsetup[ ri,5] > 0){ 
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = 
          tform(ukfstates[nlatent +PARSsetup[ri,5], statei ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ); 
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
      if(T0check[rowi]==0) etaupd[rowi] = etaprior[rowi] +  discreteDIFFUSION * etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
    }
  
    if(nobsi==0 && intoverstates==1 ) {
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
            if(manifesttype[wi]==2 && Y[rowi,wi] != 99999) ypredcov[wi,wi] = ypredcov[wi,wi] + square(fabs((ypred[wi] - round(ypred[wi])))); 
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
          ;
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
    }

    for(ri in 1:size(PARSsetup)){
      if(PARSsetup[ ri,5] > 0){ 
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = tform(ukfstates[nlatent +PARSsetup[ri,5], statei ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ); 
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
          if(manifesttype[wi]==2 && Y[rowi,wi] != 99999) ypredcov[wi,wi] = ypredcov[wi,wi] + square(fabs((ypred[wi] - round(ypred[wi])))); 
        
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
          "  sDRIFT ", sDRIFT, " sDIFFUSION ", sDIFFUSION, " sCINT ", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
          "  sT0VAR", sT0VAR, 
          "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        if(lineardynamics==1) print("discreteDRIFT ",discreteDRIFT,"  discreteCINT ", discreteCINT, "  discreteDIFFUSION ", discreteDIFFUSION)
      }
      if(verbose > 2) print("ukfstates ", ukfstates, "  ukfmeasures ", ukfmeasures);

      if(size(cindex) > 0){
         ypredcov_sqrt[cindex,cindex]=cholesky_decompose(makesym(ypredcov[cindex,cindex]));
         errtrans[(cobscount+1):(cobscount+size(cindex))] = mdivide_left_tri_low(ypredcov_sqrt[cindex,cindex], err[cindex]); //transform pred errors to standard normal dist and collect
         errscales[(cobscount+1):(cobscount+size(cindex))] = log(diagonal(ypredcov_sqrt[cindex,cindex])); //account for transformation of scale in loglik
      }
  

    }//end nobs > 0 section
  
  }//end rowi
if(savescores==1) {
  etaprior_out=etaprior;
  etaupd_out = etaupd;
}

  if((intoverstates==1 || sum(ncont_y) > 0)) ll = ll + normal_lpdf(errtrans|0,1) - sum(errscales);
}

}
      
model{

  if(nopriors==0){
    target += normal_lpdf(rawpopmeans|0,1);
  
    if(ntipred > 0){ 
     tipredeffectparams ~ normal(0,tipredeffectscale);
     tipredsimputed ~ normal(0,tipredsimputedscale);
    }
    
    if(nindvarying > 0){
      if(nindvarying >1) sqrtpcov ~ normal(0,1);
      if(ukfpop==0) baseindparams ~ normal(0,1);
      rawpopsdbase ~ normal(0,1);
    }

  } //end pop priors section
  
  if(intoverstates==0)etaupdbasestates ~ normal(0,1);
  
  target += ll;
  
  if(verbose > 0) print("lp = ", target());
  
}
generated quantities{
  vector[nparams] popmeans;
  vector[nparams] popsd;
  matrix[nindvarying,nindvarying] rawpopcov;
  matrix[nindvarying,nindvarying] rawpopcorr;
  matrix[nparams,ntipred] linearTIPREDEFFECT;

  matrix[ T0MEANSsetup_rowcount ? max(T0MEANSsetup[,1]) : 0, T0MEANSsetup_rowcount ? max(T0MEANSsetup[,2]) : 0 ] pop_T0MEANS[T0MEANSsubindex[1]];
matrix[ LAMBDAsetup_rowcount ? max(LAMBDAsetup[,1]) : 0, LAMBDAsetup_rowcount ? max(LAMBDAsetup[,2]) : 0 ] pop_LAMBDA[LAMBDAsubindex[1]];
matrix[ DRIFTsetup_rowcount ? max(DRIFTsetup[,1]) : 0, DRIFTsetup_rowcount ? max(DRIFTsetup[,2]) : 0 ] pop_DRIFT[DRIFTsubindex[1]];
matrix[ DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,1]) : 0, DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,2]) : 0 ] pop_DIFFUSION[DIFFUSIONsubindex[1]];
matrix[ MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,1]) : 0, MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,2]) : 0 ] pop_MANIFESTVAR[MANIFESTVARsubindex[1]];
matrix[ MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,1]) : 0, MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,2]) : 0 ] pop_MANIFESTMEANS[MANIFESTMEANSsubindex[1]];
matrix[ CINTsetup_rowcount ? max(CINTsetup[,1]) : 0, CINTsetup_rowcount ? max(CINTsetup[,2]) : 0 ] pop_CINT[CINTsubindex[1]];
matrix[ T0VARsetup_rowcount ? max(T0VARsetup[,1]) : 0, T0VARsetup_rowcount ? max(T0VARsetup[,2]) : 0 ] pop_T0VAR[T0VARsubindex[1]];
matrix[ TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,1]) : 0, TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,2]) : 0 ] pop_TDPREDEFFECT[TDPREDEFFECTsubindex[1]];
matrix[ PARSsetup_rowcount ? max(PARSsetup[,1]) : 0, PARSsetup_rowcount ? max(PARSsetup[,2]) : 0 ] pop_PARS[PARSsubindex[1]];

  matrix[nlatent,nlatent] asympop_DIFFUSION[asymDIFFUSIONsubindex[1]]; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] asympop_CINT[asymCINTsubindex[1]]; // latent process asymptotic level
  

  matrix[ T0MEANSsetup_rowcount ? max(T0MEANSsetup[,1]) : 0, T0MEANSsetup_rowcount ? max(T0MEANSsetup[,2]) : 0 ] T0MEANS[T0MEANSsubindex[nsubjects]];
matrix[ LAMBDAsetup_rowcount ? max(LAMBDAsetup[,1]) : 0, LAMBDAsetup_rowcount ? max(LAMBDAsetup[,2]) : 0 ] LAMBDA[LAMBDAsubindex[nsubjects]];
matrix[ DRIFTsetup_rowcount ? max(DRIFTsetup[,1]) : 0, DRIFTsetup_rowcount ? max(DRIFTsetup[,2]) : 0 ] DRIFT[DRIFTsubindex[nsubjects]];
matrix[ DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,1]) : 0, DIFFUSIONsetup_rowcount ? max(DIFFUSIONsetup[,2]) : 0 ] DIFFUSION[DIFFUSIONsubindex[nsubjects]];
matrix[ MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,1]) : 0, MANIFESTVARsetup_rowcount ? max(MANIFESTVARsetup[,2]) : 0 ] MANIFESTVAR[MANIFESTVARsubindex[nsubjects]];
matrix[ MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,1]) : 0, MANIFESTMEANSsetup_rowcount ? max(MANIFESTMEANSsetup[,2]) : 0 ] MANIFESTMEANS[MANIFESTMEANSsubindex[nsubjects]];
matrix[ CINTsetup_rowcount ? max(CINTsetup[,1]) : 0, CINTsetup_rowcount ? max(CINTsetup[,2]) : 0 ] CINT[CINTsubindex[nsubjects]];
matrix[ T0VARsetup_rowcount ? max(T0VARsetup[,1]) : 0, T0VARsetup_rowcount ? max(T0VARsetup[,2]) : 0 ] T0VAR[T0VARsubindex[nsubjects]];
matrix[ TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,1]) : 0, TDPREDEFFECTsetup_rowcount ? max(TDPREDEFFECTsetup[,2]) : 0 ] TDPREDEFFECT[TDPREDEFFECTsubindex[nsubjects]];
matrix[ PARSsetup_rowcount ? max(PARSsetup[,1]) : 0, PARSsetup_rowcount ? max(PARSsetup[,2]) : 0 ] PARS[PARSsubindex[nsubjects]];

  matrix[nlatent,nlatent] asymDIFFUSION[asymDIFFUSIONsubindex[nsubjects]]; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] asymCINT[asymCINTsubindex[nsubjects]]; // latent process asymptotic level
  

vector[nmanifest] Ygen[ngenerations, ndatapoints];
for(geni in 1:ngenerations) Ygen[geni,,] = rep_array(rep_vector(99999,nmanifest), ndatapoints);
for(geni in 0:ngenerations){

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
  vector[nt0meansstationary ? nlatent : 0] sasymCINT;
  matrix[nmanifest,nmanifest] sMANIFESTVAR; 
  matrix[nmanifest,1] sMANIFESTMEANS;
  matrix[nmanifest,nlatent] sLAMBDA;
  matrix[ntdpred ? nlatent : 0,ntdpred] sTDPREDEFFECT;
  matrix[PARSsetup_rowcount ? max(PARSsetup[,1]) : 0 ,PARSsetup_rowcount ? max(PARSsetup[,2]) : 0] sPARS;

  //ukf approximation parameters
  if(ukf==1) k = 0.5;

  if(lineardynamics) discreteDIFFUSION = rep_matrix(0,nlatent,nlatent); //in case some elements remain zero due to derrind
  
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
        asquared =  2.0/sqrt(0.0+nlatentpop+ndynerror) * 1e-1;
        l = asquared * (nlatentpop + ndynerror + k) - (nlatentpop+ndynerror); 
        sqrtukfadjust = sqrt(0.0+nlatentpop + ndynerror +l);
      }
    }

    if(T0check[rowi] == 1) { // calculate initial matrices if this is first row for si

    {
  vector[nparams] rawindparams;
  vector[nparams] tipredaddition;
  vector[nparams] indvaraddition;
  
  if(si==1 || (si > 1 && (nindvarying >0 || ntipred > 0))){
    tipredaddition = rep_vector(0,nparams);
    indvaraddition = rep_vector(0,nparams);

    if(nindvarying > 0 && ukfpop==0) indvaraddition[indvaryingindex] = rawpopcovsqrt * baseindparams[(1+(si-1)*nindvarying):(si*nindvarying)];
  
    if(ntipred > 0) tipredaddition = TIPREDEFFECT * tipreds[si]';
  
    rawindparams = rawpopmeans + tipredaddition + indvaraddition;
  }

  if(si <= T0MEANSsubindex[nsubjects]){
    for(ri in 1:size(T0MEANSsetup)){
      if(si==1 || T0MEANSsetup[ri,5] > 0 || T0MEANSsetup[ri,6] > 0){
        sT0MEANS[T0MEANSsetup[ ri,1], T0MEANSsetup[ri,2]] = T0MEANSsetup[ri,3] ? tform(rawindparams[ T0MEANSsetup[ri,3] ], T0MEANSsetup[ri,4], T0MEANSvalues[ri,2], T0MEANSvalues[ri,3], T0MEANSvalues[ri,4] ) : T0MEANSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= LAMBDAsubindex[nsubjects]){
    for(ri in 1:size(LAMBDAsetup)){
      if(si==1 || LAMBDAsetup[ri,5] > 0 || LAMBDAsetup[ri,6] > 0){
        sLAMBDA[LAMBDAsetup[ ri,1], LAMBDAsetup[ri,2]] = LAMBDAsetup[ri,3] ? tform(rawindparams[ LAMBDAsetup[ri,3] ], LAMBDAsetup[ri,4], LAMBDAvalues[ri,2], LAMBDAvalues[ri,3], LAMBDAvalues[ri,4] ) : LAMBDAvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= DRIFTsubindex[nsubjects]){
    for(ri in 1:size(DRIFTsetup)){
      if(si==1 || DRIFTsetup[ri,5] > 0 || DRIFTsetup[ri,6] > 0){
        sDRIFT[DRIFTsetup[ ri,1], DRIFTsetup[ri,2]] = DRIFTsetup[ri,3] ? tform(rawindparams[ DRIFTsetup[ri,3] ], DRIFTsetup[ri,4], DRIFTvalues[ri,2], DRIFTvalues[ri,3], DRIFTvalues[ri,4] ) : DRIFTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= DIFFUSIONsubindex[nsubjects]){
    for(ri in 1:size(DIFFUSIONsetup)){
      if(si==1 || DIFFUSIONsetup[ri,5] > 0 || DIFFUSIONsetup[ri,6] > 0){
        sDIFFUSION[DIFFUSIONsetup[ ri,1], DIFFUSIONsetup[ri,2]] = DIFFUSIONsetup[ri,3] ? tform(rawindparams[ DIFFUSIONsetup[ri,3] ], DIFFUSIONsetup[ri,4], DIFFUSIONvalues[ri,2], DIFFUSIONvalues[ri,3], DIFFUSIONvalues[ri,4] ) : DIFFUSIONvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= MANIFESTVARsubindex[nsubjects]){
    for(ri in 1:size(MANIFESTVARsetup)){
      if(si==1 || MANIFESTVARsetup[ri,5] > 0 || MANIFESTVARsetup[ri,6] > 0){
        sMANIFESTVAR[MANIFESTVARsetup[ ri,1], MANIFESTVARsetup[ri,2]] = MANIFESTVARsetup[ri,3] ? tform(rawindparams[ MANIFESTVARsetup[ri,3] ], MANIFESTVARsetup[ri,4], MANIFESTVARvalues[ri,2], MANIFESTVARvalues[ri,3], MANIFESTVARvalues[ri,4] ) : MANIFESTVARvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= MANIFESTMEANSsubindex[nsubjects]){
    for(ri in 1:size(MANIFESTMEANSsetup)){
      if(si==1 || MANIFESTMEANSsetup[ri,5] > 0 || MANIFESTMEANSsetup[ri,6] > 0){
        sMANIFESTMEANS[MANIFESTMEANSsetup[ ri,1], MANIFESTMEANSsetup[ri,2]] = MANIFESTMEANSsetup[ri,3] ? tform(rawindparams[ MANIFESTMEANSsetup[ri,3] ], MANIFESTMEANSsetup[ri,4], MANIFESTMEANSvalues[ri,2], MANIFESTMEANSvalues[ri,3], MANIFESTMEANSvalues[ri,4] ) : MANIFESTMEANSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= CINTsubindex[nsubjects]){
    for(ri in 1:size(CINTsetup)){
      if(si==1 || CINTsetup[ri,5] > 0 || CINTsetup[ri,6] > 0){
        sCINT[CINTsetup[ ri,1], CINTsetup[ri,2]] = CINTsetup[ri,3] ? tform(rawindparams[ CINTsetup[ri,3] ], CINTsetup[ri,4], CINTvalues[ri,2], CINTvalues[ri,3], CINTvalues[ri,4] ) : CINTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= T0VARsubindex[nsubjects]){
    for(ri in 1:size(T0VARsetup)){
      if(si==1 || T0VARsetup[ri,5] > 0 || T0VARsetup[ri,6] > 0){
        sT0VAR[T0VARsetup[ ri,1], T0VARsetup[ri,2]] = T0VARsetup[ri,3] ? tform(rawindparams[ T0VARsetup[ri,3] ], T0VARsetup[ri,4], T0VARvalues[ri,2], T0VARvalues[ri,3], T0VARvalues[ri,4] ) : T0VARvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= TDPREDEFFECTsubindex[nsubjects]){
    for(ri in 1:size(TDPREDEFFECTsetup)){
      if(si==1 || TDPREDEFFECTsetup[ri,5] > 0 || TDPREDEFFECTsetup[ri,6] > 0){
        sTDPREDEFFECT[TDPREDEFFECTsetup[ ri,1], TDPREDEFFECTsetup[ri,2]] = TDPREDEFFECTsetup[ri,3] ? tform(rawindparams[ TDPREDEFFECTsetup[ri,3] ], TDPREDEFFECTsetup[ri,4], TDPREDEFFECTvalues[ri,2], TDPREDEFFECTvalues[ri,3], TDPREDEFFECTvalues[ri,4] ) : TDPREDEFFECTvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      

  if(si <= PARSsubindex[nsubjects]){
    for(ri in 1:size(PARSsetup)){
      if(si==1 || PARSsetup[ri,5] > 0 || PARSsetup[ri,6] > 0){
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = PARSsetup[ri,3] ? tform(rawindparams[ PARSsetup[ri,3] ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ) : PARSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      


  // perform any whole matrix transformations 
    
  if(si <= DIFFUSIONsubindex[nsubjects] && lineardynamics) sDIFFUSION = sdcovsqrt2cov(sDIFFUSION, 0);

    if(si <= asymDIFFUSIONsubindex[nsubjects]) {
      if(ndiffusion < nlatent) sasymDIFFUSION = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

      if(continuoustime==1) sasymDIFFUSION[ derrind, derrind] = to_matrix( 
      -( kron_prod( sDRIFT[ derrind, derrind ], IIlatent[ derrind, derrind ]) + 
         kron_prod(IIlatent[ derrind, derrind ], sDRIFT[ derrind, derrind ]) ) \ 
      to_vector( sDIFFUSION[ derrind, derrind ] + IIlatent[ derrind, derrind ] * 1e-5), ndiffusion,ndiffusion);

      if(continuoustime==0) sasymDIFFUSION[derrind, derrind] = to_matrix( (IIlatent2[ derrind, derrind ] - 
        kron_prod(sDRIFT[derrind, derrind  ], 
          sDRIFT[derrind, derrind  ])) * 
        to_vector(sDIFFUSION[derrind, derrind  ]) , ndiffusion, ndiffusion);
    } //end asymdiffusion loops
          
    if(nt0meansstationary > 0){
      if(si <= asymCINTsubindex[nsubjects]){
        if(continuoustime==1) sasymCINT =  -sDRIFT \ sCINT[ ,1 ];
        if(continuoustime==0) sasymCINT =  (IIlatent - sDRIFT) \ sCINT[,1 ];
      }
    }

          
    if(binomial==0){
      if(si <= MANIFESTVARsubindex[nsubjects]) {
         for(ri in 1:nmanifest) sMANIFESTVAR[ri,ri] = square(sMANIFESTVAR[ri,ri]);
      }
    }
          
          
    if(si <= T0VARsubindex[nsubjects]) {
      if(lineardynamics * intoverstates !=0) sT0VAR = sdcovsqrt2cov(sT0VAR,0);
      if(nt0varstationary > 0) for(ri in 1:nt0varstationary){
        sT0VAR[t0varstationary[ri,1],t0varstationary[ri,2] ] = 
          sasymDIFFUSION[t0varstationary[ri,1],t0varstationary[ri,2] ];
        sT0VAR[t0varstationary[ri,2],t0varstationary[ri,1] ] = 
          sasymDIFFUSION[t0varstationary[ri,2],t0varstationary[ri,1] ];
      }
    }

    
    if(nt0meansstationary > 0){
      if(si <= T0MEANSsubindex[nsubjects]) {
        for(ri in 1:nt0meansstationary){
          sT0MEANS[t0meansstationary[ri,1] , 1] = 
            sasymCINT[t0meansstationary[ri,1] ];
        }
      }
    }
  }
  

    T0MEANS[T0MEANSsubindex[si]] = sT0MEANS; 
LAMBDA[LAMBDAsubindex[si]] = sLAMBDA; 
DRIFT[DRIFTsubindex[si]] = sDRIFT; 
DIFFUSION[DIFFUSIONsubindex[si]] = sDIFFUSION; 
MANIFESTVAR[MANIFESTVARsubindex[si]] = sMANIFESTVAR; 
MANIFESTMEANS[MANIFESTMEANSsubindex[si]] = sMANIFESTMEANS; 
CINT[CINTsubindex[si]] = sCINT; 
T0VAR[T0VARsubindex[si]] = sT0VAR; 
TDPREDEFFECT[TDPREDEFFECTsubindex[si]] = sTDPREDEFFECT; 
PARS[PARSsubindex[si]] = sPARS; 
asymDIFFUSION[asymDIFFUSIONsubindex[si]] = sasymDIFFUSION; 
asymCINT[asymCINTsubindex[si]] = sasymCINT; 


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
      if(ntdpred > 0) etaprior[rowi] = sTDPREDEFFECT * tdpreds[rowi] + sT0MEANS[,1];
      etapriorcov[rowi] =  sT0VAR;
      }

    } //end T0 matrices

if(geni > 0){

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
          discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
            quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]' );
          //discreteDIFFUSION[derrind, derrind] = discreteDIFFUSIONcalc(DRIFT[ DRIFTsubindex[si], derrind, derrind], sDIFFUSION[derrind, derrind], dT[rowi]);
          if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(discreteDIFFUSION);
        }
      }
  
      if(continuoustime==0 && T0check[rowi-1] == 1){
        discreteDRIFT=sDRIFT;
        discreteCINT=sCINT[,1];
        discreteDIFFUSION=sDIFFUSION;
        if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(discreteDIFFUSION);
      }

      if(ntdpred == 0) etaprior[rowi] = discreteDRIFT * etaupd[rowi-1] + discreteCINT;
      if(ntdpred > 0) etaprior[rowi] = discreteDRIFT * etaupd[rowi-1] + discreteCINT + sTDPREDEFFECT * tdpreds[rowi];
      if(intoverstates==1) {
        if(ndiffusion ==0) etapriorcov[rowi] = quad_form(etaupdcov[rowi-1], discreteDRIFT');
        if(ndiffusion > 0) etapriorcov[rowi] = quad_form(etaupdcov[rowi-1], discreteDRIFT') + discreteDIFFUSION;
      }
    }//end linear time update

    if(ukf==1){ //ukf time update

      if(T0check[rowi]==1) dynerror = sqrtukfadjust;
      if(T0check[rowi]==0 && lineardynamics==0) dynerror = sqrtukfadjust / sqrt(dT[rowi]); //Weiner process variance adjustment
  
      if(T0check[rowi]==0){ //compute updated sigpoints
      
      sigpoints = chol(makesym(etaupdcov[rowi-1,,] * sqrtukfadjust));
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
          ;
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
                ;
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
    }

    for(ri in 1:size(PARSsetup)){
      if(PARSsetup[ ri,5] > 0){ 
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = 
          tform(ukfstates[nlatent +PARSsetup[ri,5], statei ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ); 
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
            ;
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
    }

    for(ri in 1:size(PARSsetup)){
      if(PARSsetup[ ri,5] > 0){ 
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = 
          tform(ukfstates[nlatent +PARSsetup[ri,5], statei ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ); 
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
    
        ;
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
    }

    for(ri in 1:size(PARSsetup)){
      if(PARSsetup[ ri,5] > 0){ 
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = 
          tform(ukfstates[nlatent +PARSsetup[ri,5], statei ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ); 
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
      if(T0check[rowi]==0) etaupd[rowi] = etaprior[rowi] +  discreteDIFFUSION * etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
    }
  
    if(nobsi==0 && intoverstates==1 ) {
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
            if(manifesttype[wi]==2 && Y[rowi,wi] != 99999) ypredcov[wi,wi] = ypredcov[wi,wi] + square(fabs((ypred[wi] - round(ypred[wi])))); 
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
          ;
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
    }

    for(ri in 1:size(PARSsetup)){
      if(PARSsetup[ ri,5] > 0){ 
        sPARS[PARSsetup[ ri,1], PARSsetup[ri,2]] = tform(ukfstates[nlatent +PARSsetup[ri,5], statei ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ); 
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
          if(manifesttype[wi]==2 && Y[rowi,wi] != 99999) ypredcov[wi,wi] = ypredcov[wi,wi] + square(fabs((ypred[wi] - round(ypred[wi])))); 
        
        }
        K[,o] = mdivide_right(crosscov(ukfstates', ukfmeasures[o,]') /asquared, ypredcov[o,o]); 
        etaupdcov[rowi] = etapriorcov[rowi] - quad_form(ypredcov[o,o],  K[,o]');
      }

      
if(verbose > 1) {
print("rowi ",rowi, "  si ", si, "  etaprior[rowi] ",etaprior[rowi],"  etapriorcov[rowi] ",etapriorcov[rowi],
          "  etaupd[rowi] ",etaupd[rowi],"  etaupdcov[rowi] ",etaupdcov[rowi],"  ypred ",ypred,"  ypredcov ",ypredcov, "  K ",K,
          "  sDRIFT ", sDRIFT, " sDIFFUSION ", sDIFFUSION, " sCINT ", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
          "  sT0VAR", sT0VAR,
          "  rawpopsd ", rawpopsd, "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        if(lineardynamics==1) print("discreteDRIFT ",discreteDRIFT,"  discreteCINT ", discreteCINT, "  discreteDIFFUSION ", discreteDIFFUSION)
}
if(verbose > 2) print("ukfstates ", ukfstates, "  ukfmeasures ", ukfmeasures);
        
        ypredcov_sqrt[cindex,cindex]=chol(ypredcov[cindex, cindex]); //use o0, or cindex?
        for(vi in 1:nobsi){
          if(fabs(ypred[o[vi]]) > 1e10 || is_nan(ypred[o[vi]]) || is_inf(ypred[o[vi]])) {
            nobsi = 0; //set nobsi to 0 to skip update steps
            ypred[o[vi]] =99999;
          }
        }
        if(nobsi > 0){ //check nobsi again in case of problems
          if(ncont_y[rowi] > 0) Ygen[geni, rowi, cindex] = multi_normal_cholesky_rng(ypred[cindex], ypredcov_sqrt[cindex,cindex]);
          if(nbinary_y[rowi] > 0) for(obsi in 1:size(o1)) Ygen[geni, rowi, o1[obsi]] = bernoulli_rng(ypred[o1[obsi]]);
          for(vi in 1:nobsi) if(is_nan(Ygen[geni,rowi,o[vi]])) {
            Ygen[geni,rowi,o[vi]] = 99999;
            nobsi = 0;
print("pp problem2! row ", rowi);
          }
          err[o] = Ygen[geni,rowi,o] - ypred[o]; // prediction error
        }
      
  
      

      if(intoverstates==1) etaupd[rowi,] = etaprior[rowi,] + (K[,o] * err[o]);
  
      

    }//end nobs > 0 section
  } //end if geni >0 section
  }//end rowi


  
}

{
  vector[nparams] rawindparams;
  vector[nparams] tipredaddition;
  vector[nparams] indvaraddition;
  rawindparams = rawpopmeans;
  tipredaddition = rep_vector(0,nparams);
  indvaraddition = rep_vector(0,nparams);

  for(si in 1:1){

  
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
      

  if(si <= PARSsubindex[1]){
    for(ri in 1:size(PARSsetup)){
      pop_PARS[si, PARSsetup[ ri,1], PARSsetup[ri,2]] = PARSsetup[ri,3] ? tform(rawindparams[ PARSsetup[ri,3] ], PARSsetup[ri,4], PARSvalues[ri,2], PARSvalues[ri,3], PARSvalues[ri,4] ) : PARSvalues[ri,1]; //either transformed, scaled and offset free par, or fixed value
    }
  }
      


  // perform any whole matrix transformations 
    
  if(si <= DIFFUSIONsubindex[1] && lineardynamics !=0 ) pop_DIFFUSION[si] = sdcovsqrt2cov(pop_DIFFUSION[si], 0);

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
      if(lineardynamics * intoverstates !=0) pop_T0VAR[si] = sdcovsqrt2cov(pop_T0VAR[si],0);
      if(nt0varstationary > 0) for(ri in 1:nt0varstationary){
        pop_T0VAR[si,t0varstationary[ri,1],t0varstationary[ri,2] ] = 
          asympop_DIFFUSION[si,t0varstationary[ri,1],t0varstationary[ri,2] ];
        pop_T0VAR[si,t0varstationary[ri,2],t0varstationary[ri,1] ] = 
          asympop_DIFFUSION[si,t0varstationary[ri,2],t0varstationary[ri,1] ];
      }
    }

    
    if(nt0meansstationary > 0){
      if(si <= T0MEANSsubindex[1]) {
        for(ri in 1:nt0meansstationary){
          pop_T0MEANS[si,t0meansstationary[ri,1] , 1] = 
            asympop_CINT[ asymCINTsubindex[si], t0meansstationary[ri,1] ];
        }
      }
    }
  }
} //end subject loop
  

rawpopcorr = tcrossprod(rawpopcorrsqrt);
rawpopcov = tcrossprod(rawpopcovsqrt);

popsd = rep_vector(0,nparams);
{
vector[nparams] rawpopsdfull;
rawpopsdfull[indvaryingindex] = rawpopsd; //base for calculations

    for(ri in 1:dims(popsetup)[1]){
      if(popsetup[ri,3] !=0) {

        popmeans[popsetup[ ri,3]] = tform(rawpopmeans[popsetup[ri,3] ], popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4] ); 

        popsd[popsetup[ ri,3]] = popsetup[ ri,5] ? 
          fabs(tform(
            rawpopmeans[popsetup[ri,3] ]  + rawpopsdfull[popsetup[ ri,3]], popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4]) -
           tform(
            rawpopmeans[popsetup[ri,3] ]  - rawpopsdfull[popsetup[ ri,3]], popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4] - popsd[popsetup[ ri,3]]) ) /2 : 
          0; 

        if(ntipred > 0){
          for(tij in 1:ntipred){
            if(TIPREDEFFECTsetup[popsetup[ri,3],tij] ==0) {
              linearTIPREDEFFECT[popsetup[ri,3],tij] = 0;
            } else {
            linearTIPREDEFFECT[popsetup[ri,3],tij] = (
              tform(rawpopmeans[popsetup[ri,3] ] + TIPREDEFFECT[popsetup[ri,3],tij] * .01, popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4] ) -
              tform(rawpopmeans[popsetup[ri,3] ] - TIPREDEFFECT[popsetup[ri,3],tij] * .01, popsetup[ri,4], popvalues[ri,2], popvalues[ri,3], popvalues[ri,4] )
              ) /2 * 100;
            }
         }
        }
      }
    }
}
}
