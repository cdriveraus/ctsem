
functions{
   matrix expm2(matrix M,int[] z){
    matrix[rows(M),rows(M)] out;
    int z1[sum(z)];
    int z0[rows(M)-sum(z)];
    int cz1 = 1;
    int cz0 = 1;
    for(i in 1:rows(M)){
      if(z[i] == 1){
        z1[cz1] = i;
        cz1 += 1;
      }
      if(z[i] == 0){
        z0[cz0] = i;
        cz0 += 1;
      }
    }
    if(size(z1) > 0) out[z1,z1] = matrix_exp(M[z1,z1]);
    out[z0,] = rep_matrix(0,size(z0),rows(M));
    out[,z0] = rep_matrix(0,rows(M),size(z0));
    if(size(z0) > 0) for(i in 1:size(z0)) out[z0[i],z0[i]] = exp(M[z0[i],z0[i]]);
    return out;
  }

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
    if(transform==0) out = inneroffset + meanscale * multiplier * param + offset;
if(transform==1) out = multiplier * log1p(exp(inneroffset + meanscale * param)) + offset;
if(transform==2) out = multiplier * exp(inneroffset + meanscale * param) + offset;
if(transform==3) out = multiplier * exp(inneroffset + meanscale * param)/(1 + exp(param)) +offset;
if(transform==4) out = multiplier * (inneroffset + meanscale * param)^3 + offset;
if(transform==50) out = meanscale*multiplier;
if(transform==51) out = multiplier*(exp(inneroffset+meanscale*param)*meanscale/(1+exp(inneroffset+meanscale*param)));
if(transform==52) out = multiplier*(exp(inneroffset+meanscale*param)*meanscale);
if(transform==53) out = multiplier*(exp(inneroffset+meanscale*param)*meanscale)/(1+exp(param))-multiplier*exp(inneroffset+meanscale*param)*exp(param)/(1+exp(param))^2;
if(transform==54) out = multiplier*(3*(meanscale*(inneroffset+meanscale*param)^2));

    return out;
  }
  
  
  real Jtform(real param, int transform, data real multiplier, data real meanscale, data real offset, data real inneroffset){
    real out;
    if(transform==50) out = meanscale*multiplier;
if(transform==51) out = multiplier*(exp(inneroffset+meanscale*param)*meanscale/(1+exp(inneroffset+meanscale*param)));
if(transform==52) out = multiplier*(exp(inneroffset+meanscale*param)*meanscale);
if(transform==53) out = multiplier*(exp(inneroffset+meanscale*param)*meanscale)/(1+exp(param))-multiplier*exp(inneroffset+meanscale*param)*exp(param)/(1+exp(param))^2;
if(transform==54) out = multiplier*(3*(meanscale*(inneroffset+meanscale*param)^2));

    return out;
  }

}
data {
  int<lower=0> ndatapoints;
  int<lower=1> nmanifest;
  int<lower=1> nlatent;
  int nlatentpop;
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
  int nldynamics;
  vector[ntdpred] tdpreds[ndatapoints];
  
  real maxtimestep;
  real time[ndatapoints];
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
  int drcintoffdiag[nlatent+1];

  int manifesttype[nmanifest];
  int nbinary_y[ndatapoints];  // number of observed binary variables per observation
  int whichbinary_y[ndatapoints, nmanifest]; // index of which variables are observed and binary per observation
  int ncont_y[ndatapoints];  // number of observed continuous variables per observation
  int whichcont_y[ndatapoints, nmanifest]; // index of which variables are observed and continuous per observation
  
  int intoverpop;
  real ukfspread;
  int ukffull;
  int nlmeasurement;
  int intoverstates;
  int verbose; //level of printing during model fit

  int T0MEANSsubindex;
int LAMBDAsubindex;
int DRIFTsubindex;
int DIFFUSIONsubindex;
int MANIFESTVARsubindex;
int MANIFESTMEANSsubindex;
int CINTsubindex;
int T0VARsubindex;
int TDPREDEFFECTsubindex;
int PARSsubindex;
int asymCINTsubindex;
int asymDIFFUSIONsubindex;
int DIFFUSIONcovsubindex;
  int TIPREDEFFECTsetup[nparams, ntipred];
  int nrowpopsetup;
  int nrowmatsetup;
  int popsetup[nrowpopsetup,8];
  int matsetup[nrowmatsetup,8];
  real popvalues[nrowpopsetup,6];
  real matvalues[nrowmatsetup,6];
  int matrixdims[10,2];
  int savescores;
  int fixedsubpars;
  vector[fixedsubpars ? nindvarying : 0] fixedindparams[fixedsubpars ? nsubjects : 0];
  int dokalman;
  int dokalmanrows[ndatapoints];
  real Jstep;
  real dokalmanpriormodifier;
  int intoverpopindvaryingindex[nindvarying];
  int sJAxdrift[nlatentpop,nlatentpop];
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
  //vector[ (( (ntipredeffects-1) * (1-nopriors) ) > 0) ? 1 : 0] tipredglobalscalepar;
  
  vector[intoverstates ? 0 : nlatentpop*ndatapoints] etaupdbasestates; //sampled latent states posterior
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying, nindvarying] rawpopcovsqrt; 

  real ll = 0;
  vector[nmanifest+nmanifest+ (savescores ? nmanifest*2+nlatent*2 : 0)] kalman[savescores ? ndatapoints : 0];
  matrix[nlatent,nlatent] etapriorcov[savescores ? ndatapoints : 0];
  matrix[nlatent,nlatent] etaupdcov[savescores ? ndatapoints : 0];
     matrix[matrixdims[1, 1], matrixdims[1, 2] ] T0MEANS[T0MEANSsubindex  ? nsubjects : 1]; 
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] LAMBDA[LAMBDAsubindex  ? nsubjects : 1]; 
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] DRIFT[DRIFTsubindex  ? nsubjects : 1]; 
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] DIFFUSION[DIFFUSIONsubindex  ? nsubjects : 1]; 
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] MANIFESTVAR[MANIFESTVARsubindex  ? nsubjects : 1]; 
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] MANIFESTMEANS[MANIFESTMEANSsubindex  ? nsubjects : 1]; 
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] CINT[CINTsubindex  ? nsubjects : 1]; 
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] T0VAR[T0VARsubindex  ? nsubjects : 1]; 
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] TDPREDEFFECT[TDPREDEFFECTsubindex  ? nsubjects : 1]; 
      matrix[matrixdims[10, 1], matrixdims[10, 2] ] PARS[PARSsubindex  ? nsubjects : 1];

  matrix[nlatent,nlatent] asymDIFFUSION[asymDIFFUSIONsubindex ? nsubjects : 1]; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] asymCINT[asymCINTsubindex ? nsubjects : 1]; // latent process asymptotic level
matrix[nlatent, nlatent] DIFFUSIONcov[DIFFUSIONcovsubindex ? nsubjects : 1];
     matrix[matrixdims[1, 1], matrixdims[1, 2] ] pop_T0MEANS; 
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] pop_LAMBDA; 
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] pop_DRIFT; 
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] pop_DIFFUSION; 
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] pop_MANIFESTVAR; 
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] pop_MANIFESTMEANS; 
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] pop_CINT; 
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] pop_T0VAR; 
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] pop_TDPREDEFFECT; 
      matrix[matrixdims[10, 1], matrixdims[10, 2] ] pop_PARS;

  matrix[nlatent,nlatent] pop_asymDIFFUSION; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] pop_asymCINT; // latent process asymptotic level
matrix[nlatent, nlatent] pop_DIFFUSIONcov;

  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipreds; //tipred values to fill from data and, when needed, imputation vector
  matrix[nparams, ntipred] TIPREDEFFECT; //design matrix of individual time independent predictor effects
  //real tipredglobalscale = 1.0;
  
  //if( ((ntipredeffects-1) * (1-nopriors))  > 0)  tipredglobalscale = exp(tipredglobalscalepar[1]);

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

  {
  int si = 0;
  int subjectcount = 0;
  int counter = 0;
  vector[nlatentpop] eta; //latent states
  matrix[nlatentpop, nlatentpop] etacov; //covariance of latent states
  real timei = 0;
  real dt = 0;
  real integrationsteps;
  real dtsmall;
  real prevdt = 0;
  int T0check;

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] ypred;
  matrix[nmanifest, nmanifest] ypredcov;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypredcov_sqrt; 

  vector[nmanifest+nmanifest+ (savescores ? nmanifest*2+nlatent*2 : 0)] kout[ndatapoints];

  matrix[nldynamics ? nlatentpop :0,nldynamics ? nlatentpop :0] sigpoints;
  vector[nlatentpop] state = rep_vector(-1,nlatentpop); 
  matrix[nlatentpop,nlatentpop] sJAx; //Jacobian for drift
  matrix[nlatentpop,nlatentpop] sJ0; //Jacobian for t0
  matrix[nlatentpop,nlatentpop] sJtd;//diag_matrix(rep_vector(1),nlatentpop); //Jacobian for nltdpredeffect

  //linear continuous time calcs
  matrix[nlatent+1,nlatent+1] discreteDRIFT;
  matrix[nlatent,nlatent] discreteDIFFUSION;

  //dynamic system matrices
     matrix[matrixdims[1, 1], matrixdims[1, 2] ] sT0MEANS; 
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] sLAMBDA; 
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] sDRIFT; 
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] sDIFFUSION; 
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] sMANIFESTVAR; 
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] sMANIFESTMEANS; 
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] sCINT; 
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] sT0VAR; 
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] sTDPREDEFFECT; 
      matrix[matrixdims[10, 1], matrixdims[10, 2] ] sPARS;

  matrix[nlatent,nlatent] sasymDIFFUSION; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] sasymCINT; // latent process asymptotic level
matrix[nlatent, nlatent] sDIFFUSIONcov;

  if(nldynamics==0) discreteDIFFUSION = rep_matrix(0,nlatent,nlatent); //in case some elements remain zero due to derrind

  if(savescores) kout = rep_array(rep_vector(99999,rows(kout[1])),ndatapoints);
  
  for(rowi in 1:(dokalman ? ndatapoints :1)){
  if(dokalmanrows[rowi] ==1) { //used for subset selection
    matrix[nldynamics ? nlatentpop : 0, ukffull ? 2*nlatentpop +2 : nlatentpop + 2 ] ukfstates; //sampled states relevant for dynamics
    matrix[nldynamics ? nmanifest : 0 , ukffull ? 2*nlatentpop +2 : nlatentpop + 2] ukfmeasures; // expected measures based on sampled states

    T0check = ( (si == subject[rowi]) ? (T0check + 1) : 0 ) ; //if same subject, add 1 to t0check, else set to 0
    if(T0check > 0) prevdt = dt;
    dt = T0check ? time[rowi] - timei : 0;
    timei = time[rowi]; //must come after dt!
    si=subject[rowi]; //only update subject after t0check!

    
    if(T0check == 0) { // calculate initial matrices if this is first row for si
  
    
 int subjectvec[subjectcount ? 1 : 2];
 subjectvec[size(subjectvec)] = si;
 if(subjectcount == 0)  subjectvec[1] = 0; // only needed for subject 0 (pop pars)
 subjectcount = subjectcount + 1;
 for(subjectveci in 1:size(subjectvec)){
  int subi = subjectvec[subjectveci];
  vector[nparams] rawindparams;
  vector[nparams] tipredaddition = rep_vector(0,nparams);
  vector[nparams] indvaraddition = rep_vector(0,nparams);

  if(subi > 0 && nindvarying > 0 && intoverpop==0) {
    if(fixedsubpars==0) indvaraddition[indvaryingindex] = rawpopcovsqrt * baseindparams[subi];
    if(fixedsubpars==1) indvaraddition[indvaryingindex] = rawpopcovsqrt * fixedindparams[subi];
  }
  
  if(subi > 0 &&  ntipred > 0) tipredaddition = TIPREDEFFECT * tipreds[subi]';

  rawindparams = rawpopmeans + tipredaddition + indvaraddition;

    for(ri in 1:size(matsetup)){ //for each row of matrix setup
    for(statecalcs in 0:1){
        if(subi ==0 ||  //if population parameter
          (matsetup[ri,7]==8 && DIFFUSIONsubindex) ||( matsetup[ri,7] == 4 && T0VARsubindex) || //or a covariance parameter in an individually varying matrix
          (matsetup[ri,3] > 0 && (matsetup[ri,5] > 0 || matsetup[ri,6] > 0)) //or there is individual variation
          ){ //otherwise repeated values
            if( (statecalcs && matsetup[ri,8]>0) || (!statecalcs && matsetup[ri,8]==0) ){ //if doing statecalcs do them, if doing static calcs do them
              real newval;
              if(matsetup[ri,3] > 0)  newval = tform(matsetup[ri,8] ? state[ matsetup[ri,3] ] : rawindparams[ matsetup[ri,3] ], //tform static pars from rawindparams, dynamic from state
                matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
               if(matsetup[ri,3] < 1) newval = matvalues[ri, 1]; //doing this once over all subjects unless covariance matrix -- speed ups possible here, check properly!
              if(matsetup[ri, 7] == 1) sT0MEANS[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 2) sLAMBDA[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 3) sDRIFT[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 4) sDIFFUSION[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 5) sMANIFESTVAR[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 6) sMANIFESTMEANS[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 7) sCINT[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 8) sT0VAR[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 9) sTDPREDEFFECT[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 10) sPARS[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 51) sJ0[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 52) sJAx[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 53) sJtd[matsetup[ ri,1], matsetup[ri,2]] = newval;
            }
          }
        state=sT0MEANS[,1];
      }
    }

  // perform any whole matrix transformations, nonlinear calcs based on t0 in order to fill matrices
  ;
; 
  state=sT0MEANS[,1];
  ;
;
  ;
;
  ;
;
  ;
;
  
  if(subi <= (DIFFUSIONsubindex ? nsubjects : 0)) {
    sDIFFUSIONcov = sdcovsqrt2cov(sDIFFUSION,nldynamics);
  }
  if(subi <= (asymDIFFUSIONsubindex ? nsubjects : 0)) {
      if(ndiffusion < nlatent) sasymDIFFUSION = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

      if(continuoustime==1) sasymDIFFUSION[ derrind, derrind] = to_matrix( 
      -kronsum(sDRIFT[ derrind, derrind ]) \  to_vector( 
           sDIFFUSIONcov[ derrind, derrind ]), ndiffusion,ndiffusion);

      if(continuoustime==0) sasymDIFFUSION = to_matrix( (IIlatent2 - 
        sqkron_prod(sDRIFT, sDRIFT)) *  to_vector(sDIFFUSIONcov[ derrind, derrind ]), ndiffusion, ndiffusion);
    } //end asymdiffusion loops

      if(subi <= (MANIFESTVARsubindex ? nsubjects : 0)) {
         for(ri in 1:nmanifest) sMANIFESTVAR[ri,ri] = square(sMANIFESTVAR[ri,ri]);
      }
         
    if(subi <= (T0VARsubindex ? nsubjects : 0)) {
    if(intoverpop){
      sT0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopcovsqrt;
    }
      sT0VAR = makesym(sdcovsqrt2cov(sT0VAR,nldynamics),verbose,1);
      if(nt0varstationary > 0) {
        for(ri in 1:nt0varstationary){ 
          sT0VAR[t0varstationary[ri,1],t0varstationary[ri,2] ] =  sasymDIFFUSION[t0varstationary[ri,1],t0varstationary[ri,2] ];
        }
      }
    }
    
    if(nt0meansstationary > 0){
      if(subi <= (asymCINTsubindex ? nsubjects : 0)){
        if(continuoustime==1) sasymCINT =  -sDRIFT \ sCINT[ ,1 ];
        if(continuoustime==0) sasymCINT =  (IIlatent - sDRIFT) \ sCINT[,1 ];
      }
      if(subi <= (T0MEANSsubindex ? nsubjects : 0)) {
        for(ri in 1:nt0meansstationary){
          sT0MEANS[t0meansstationary[ri,1] , 1] = 
            sasymCINT[t0meansstationary[ri,1] ];
        }
      }
    }
  if( (T0MEANSsubindex > 0 && subi > 0) || (T0MEANSsubindex == 0 && subi==0) ) T0MEANS[T0MEANSsubindex ? subi : 1] = sT0MEANS; 
if( (LAMBDAsubindex > 0 && subi > 0) || (LAMBDAsubindex == 0 && subi==0) ) LAMBDA[LAMBDAsubindex ? subi : 1] = sLAMBDA; 
if( (DRIFTsubindex > 0 && subi > 0) || (DRIFTsubindex == 0 && subi==0) ) DRIFT[DRIFTsubindex ? subi : 1] = sDRIFT; 
if( (DIFFUSIONsubindex > 0 && subi > 0) || (DIFFUSIONsubindex == 0 && subi==0) ) DIFFUSION[DIFFUSIONsubindex ? subi : 1] = sDIFFUSION; 
if( (MANIFESTVARsubindex > 0 && subi > 0) || (MANIFESTVARsubindex == 0 && subi==0) ) MANIFESTVAR[MANIFESTVARsubindex ? subi : 1] = sMANIFESTVAR; 
if( (MANIFESTMEANSsubindex > 0 && subi > 0) || (MANIFESTMEANSsubindex == 0 && subi==0) ) MANIFESTMEANS[MANIFESTMEANSsubindex ? subi : 1] = sMANIFESTMEANS; 
if( (CINTsubindex > 0 && subi > 0) || (CINTsubindex == 0 && subi==0) ) CINT[CINTsubindex ? subi : 1] = sCINT; 
if( (T0VARsubindex > 0 && subi > 0) || (T0VARsubindex == 0 && subi==0) ) T0VAR[T0VARsubindex ? subi : 1] = sT0VAR; 
if( (TDPREDEFFECTsubindex > 0 && subi > 0) || (TDPREDEFFECTsubindex == 0 && subi==0) ) TDPREDEFFECT[TDPREDEFFECTsubindex ? subi : 1] = sTDPREDEFFECT; 
if( (PARSsubindex > 0 && subi > 0) || (PARSsubindex == 0 && subi==0) ) PARS[PARSsubindex ? subi : 1] = sPARS; 
if( (DIFFUSIONcovsubindex > 0 && subi > 0) || (DIFFUSIONcovsubindex == 0 && subi==0) ) DIFFUSIONcov[DIFFUSIONcovsubindex ? subi : 1] = sDIFFUSIONcov; 
if( (asymDIFFUSIONsubindex > 0 && subi > 0) || (asymDIFFUSIONsubindex == 0 && subi==0) ) asymDIFFUSION[asymDIFFUSIONsubindex ? subi : 1] = sasymDIFFUSION; 
if( (asymCINTsubindex > 0 && subi > 0) || (asymCINTsubindex == 0 && subi==0) ) asymCINT[asymCINTsubindex ? subi : 1] = sasymCINT; 

  if(subi == 0){
pop_T0MEANS = sT0MEANS; 
pop_LAMBDA = sLAMBDA; 
pop_DRIFT = sDRIFT; 
pop_DIFFUSION = sDIFFUSION; 
pop_MANIFESTVAR = sMANIFESTVAR; 
pop_MANIFESTMEANS = sMANIFESTMEANS; 
pop_CINT = sCINT; 
pop_T0VAR = sT0VAR; 
pop_TDPREDEFFECT = sTDPREDEFFECT; 
pop_PARS = sPARS; 
pop_DIFFUSIONcov = sDIFFUSIONcov; 
pop_asymDIFFUSION = sasymDIFFUSION; 
pop_asymCINT = sasymCINT; 

  }

} // end subject matrix creation
  

    etacov =  sT0VAR;
    state = sT0MEANS[,1]; //init and in case of jacobian dependencies

      if(nldynamics==0){ //initialize most parts for nl later!
        eta = sT0MEANS[,1]; //could use state instead of eta if ukf is dropped
        ;
        if(ntdpred > 0) eta[1:nlatent] += sTDPREDEFFECT * tdpreds[rowi];
      }

    } //end T0 matrices

   


    if(nldynamics==1){ //nldynamics time update
      if(T0check>0){
        vector[nlatent] base;
        real intstepi = 0;
        dtsmall = dt / ceil(dt / maxtimestep);
        state = eta;
        
        while(intstepi < dt){
          intstepi = intstepi + dtsmall;
          
          
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == 2){ //perform calcs appropriate to this section
            real newval;
            newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
            if(matsetup[ri, 7] == 3) sDRIFT[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 4) sDIFFUSION[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 7) sCINT[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 52) sJAx[matsetup[ ri,1], matsetup[ri,2]] = newval;
            }
          }
          
             
      for(ri in 1:nlatentpop){ 
        for(ci in 1:nlatentpop){
          if(sJAxdrift[ri,ci]) sJAx[ri,ci]=sDRIFT[ri,ci]; //set jacobian to drift where appropriate
        }
      }
          
          
          if(continuoustime==1){
            matrix[nlatentpop,nlatentpop] Je;
            matrix[nlatent*2,nlatent*2] dQi;
            Je= matrix_exp(sJAx * dtsmall) ;
            discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),rep_vector(0,nlatent+1)') * dtsmall,drcintoffdiag);
            sasymDIFFUSION = to_matrix(  -kronsum(sJAx[1:nlatent,1:nlatent]) \ to_vector(tcrossprod(sDIFFUSION)), nlatent,nlatent);
            discreteDIFFUSION =  sasymDIFFUSION - quad_form( sasymDIFFUSION, Je[1:nlatent,1:nlatent]' );
                  if(verbose > 1) print("rowi = ",rowi, "state = ", state);
      if(verbose > 1)  print("etacov = ",etacov," sasymDIFFUSION = ",sasymDIFFUSION," sDIFFUSION = ",sDIFFUSION);
 if(verbose > 1) print("sJAx = ",sJAx);
            etacov = quad_form(etacov, Je');
            etacov[1:nlatent,1:nlatent] += discreteDIFFUSION; //may need improving
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent];
          }

        if(continuoustime==0){ 
          for(nli in (nlatent+1):nlatentpop) sJAx[nli,nli] = 1;
            etacov = quad_form(etacov, sJAx');
            etacov[1:nlatent,1:nlatent] += sDIFFUSIONcov; //may need improving re sDIFFUSION
            discreteDRIFT=append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),rep_matrix(0,1,nlatent+1));
            discreteDRIFT[nlatent+1,nlatent+1] = 1;
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent];
          }
        }
      } // end of non t0 time update
  
    if(T0check==0){ //nl t0
    state = sT0MEANS[,1]; //in case of t0 dependencies, may have missingness
    if(verbose > 1) print("rowi = ",rowi, "  state = ",sT0MEANS);
    if(verbose > 1) print("etacov = ",etacov);
    
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == 1){ //perform calcs appropriate to this section
            real newval;
            newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
            if(matsetup[ri, 7] == 1) sT0MEANS[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 8) sT0VAR[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 51) sJ0[matsetup[ ri,1], matsetup[ri,2]] = newval;
            }
          }
    ;
    
      state = sT0MEANS[,1];
      etacov= quad_form_sym(sT0VAR, sJ0');
    if(verbose > 1) print("rowi = ",rowi,"  state = ",sT0MEANS);
    if(verbose > 1) print("sJ0 = ",sJ0);
    if(verbose > 1) print("etacov = ",etacov);
    
    } 
    if(ntdpred > 0) {
    
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == 3){ //perform calcs appropriate to this section
            real newval;
            newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
            if(matsetup[ri, 7] == 9) sTDPREDEFFECT[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 53) sJtd[matsetup[ ri,1], matsetup[ri,2]] = newval;
            }
          }
    ;
      state[1:nlatent] +=   (sTDPREDEFFECT * tdpreds[rowi]); //tdpred effect only influences at observed time point
      if(verbose > 1)  print("state = ", state);
      if(verbose > 1)  print("etacov = ",etacov);
 if(verbose > 1) print("sJtd = ",sJtd);
      etacov = quad_form_sym(etacov,sJtd');
     }
     eta=state;
  } // end non linear time update


  if(savescores==1) {
    kout[rowi,(nmanifest*4+1):(nmanifest*4+nlatent)] = eta[1:nlatent];
    etapriorcov[rowi,]=etacov[1:nlatent,1:nlatent];
    if(nobs_y[rowi] == 0) etaupdcov[rowi,]=etacov[1:nlatent,1:nlatent];
  }
  if(verbose > 1) print("etaprior = ", eta, " etapriorcov = ",etacov);

    if(intoverstates==0 && nldynamics == 0) {
      if(T0check==0) eta += cholesky_decompose(sT0VAR) * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)];
      if(T0check>0) eta +=  discreteDIFFUSION * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)];
    }

    if (nobs_y[rowi] > 0) {  // if some observations create right size matrices for missingness and calculate...
    
      int o[nobs_y[rowi]]= whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
      int o1[nbinary_y[rowi]]= whichbinary_y[rowi,1:nbinary_y[rowi]];
      int o0[ncont_y[rowi]]= whichcont_y[rowi,1:ncont_y[rowi]];

      if(nlmeasurement==0){ //linear measurement
      ;

        if(intoverstates==1) { //classic kalman
          ypred[o] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * eta[1:nlatent];
          if(nbinary_y[rowi] > 0) ypred[o1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * eta[1:nlatent])));
          ypredcov[o,o] = quad_form(etacov[1:nlatent,1:nlatent], sLAMBDA[o,]') + sMANIFESTVAR[o,o];
          for(wi in 1:nmanifest){ 
            if(manifesttype[wi]==1 && Y[rowi,wi] != 99999) ypredcov[wi,wi] += fabs((ypred[wi] - 1) .* (ypred[wi]));
            if(manifesttype[wi]==2 && Y[rowi,wi] != 99999) ypredcov[wi,wi] += square(fabs((ypred[wi] - round(ypred[wi])))); 
          }
          K[,o] = mdivide_right(etacov * append_row(sLAMBDA[o,]',rep_matrix(0,nlatentpop-nlatent,nmanifest)[,o]), ypredcov[o,o]); 
          etacov += -K[,o] * append_col(sLAMBDA[o,],rep_matrix(0,nmanifest,nlatentpop-nlatent)[o,]) * etacov;
        }
        if(intoverstates==0) { //sampled states
          if(ncont_y[rowi] > 0) {
            ypred[o0] = sMANIFESTMEANS[o0,1] + sLAMBDA[o0,] * eta[1:nlatent];
            ypredcov_sqrt[o0,o0] = sqrt(sMANIFESTVAR[o0,o0]);
          }
          if(nbinary_y[rowi] > 0) ypred[o1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * eta[1:nlatent])));
        }
      } 

      if(nlmeasurement==1){ //ukf measurement
        matrix[nmanifest,cols(ukfmeasures)] merrorstates;
        sigpoints = cholesky_decompose(makesym(etacov,verbose,1)) * sqrtukfadjust;
        
        for(statei in 2:cols(ukfmeasures) ){ //for each ukf state sample
          state = eta; 
          if(statei > (2+nlatentpop)){
            state += -sigpoints[,statei-(2+nlatentpop)];
          } else {
            if(statei > 2) state += sigpoints[,statei-2]; 
          }
        ukfstates[,statei] = state;
        
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == 4){ //perform calcs appropriate to this section
            real newval;
            newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
            if(matsetup[ri, 7] == 2) sLAMBDA[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 5) sMANIFESTVAR[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 6) sMANIFESTMEANS[matsetup[ ri,1], matsetup[ri,2]] = newval;
            }
          }

          ;

  
          ukfmeasures[, statei] = sMANIFESTMEANS[,1] + sLAMBDA * state[1:nlatent];
          if(nbinary_y[rowi] > 0) {
            ukfmeasures[o1 , statei] = to_vector(inv_logit(to_array_1d(ukfmeasures[o1 , statei])));
          }
        
          merrorstates[,statei] = sqrt(diagonal(sMANIFESTVAR));
          for(wi in 1:nmanifest){ 
            if(manifesttype[wi]==1 && Y[rowi,wi] != 99999) merrorstates[wi,statei] += fabs((ukfmeasures[wi,statei] - 1) .* (ukfmeasures[wi,statei])); 
            if(manifesttype[wi]==2 && Y[rowi,wi] != 99999) merrorstates[wi,statei] += square(fabs((ukfmeasures[wi,statei]- round(ukfmeasures[wi,statei])))); 
          }
        } //end state loop 
        if(ukffull == 1) {
          merrorstates[,1] = merrorstates[,2];
          ukfmeasures[ , 1] = ukfmeasures [,2];
          ukfstates[ , 1] = ukfstates [,2];
          ypred[o] = colMeans(ukfmeasures[o,]'); 
          ypredcov[o,o] = cov_of_matrix(ukfmeasures[o,]') /asquared + diag_matrix(square(colMeans(merrorstates[o,]'))); //
          K[,o] = mdivide_right(crosscov(ukfstates', ukfmeasures[o,]') /asquared, ypredcov[o,o]); 
        }
        if(ukffull == 0){
          ypred[o] = ukfmeasures[o,2];
          for(ci in 3:cols(ukfmeasures)) ukfmeasures[o,ci] += -ukfmeasures[o,2];
          for(ci in 3:cols(ukfstates)) ukfstates[,ci] += -ukfstates[,2];
          ypredcov[o,o] = tcrossprod(ukfmeasures[o,3:(nlatentpop+2)]) /asquared / (nlatentpop +.5) + diag_matrix(square(merrorstates[o,2]));
          K[,o] = mdivide_right(ukfstates[,3:cols(ukfstates)] * ukfmeasures[o,3:cols(ukfmeasures)]' /asquared / (nlatentpop+.5), ypredcov[o,o]); 
        }
        etacov +=  - quad_form(ypredcov[o,o],  K[,o]');
      } //end nlmeasurement


err[o] = Y[rowi,o] - ypred[o]; // prediction error
    
      if(savescores==1) {
        int tmpindex[nobs_y[rowi]] = o;
        for(oi in 1:ncont_y[rowi]) tmpindex[oi] +=  nmanifest*2;
        kout[rowi,tmpindex] = err[o];
        for(oi in 1:ncont_y[rowi]) tmpindex[oi] +=  nmanifest;
        kout[rowi,tmpindex] = ypred[o];
        etaupdcov[rowi,]=etacov[1:nlatent,1:nlatent];
      }
      if(intoverstates==1) eta +=  (K[,o] * err[o]);
  
      if(nbinary_y[rowi] > 0) kout[rowi,o1] =  Y[rowi,o1] .* (ypred[o1]) + (1-Y[rowi,o1]) .* (1-ypred[o1]); //intoverstates==0 && 
  
        if(verbose > 1) {
          print("rowi =",rowi, "  si =", si, "  eta =",eta,"  etacov ",etacov,
            "  ypred =",ypred,"  ypredcov ",ypredcov, "  K ",K,
            "  sDRIFT =", sDRIFT, " sDIFFUSION =", sDIFFUSION, " sCINT =", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
            "  sT0VAR", sT0VAR,  " sT0MEANS ", sT0MEANS, "sLAMBDA = ", sLAMBDA, 
            "discreteDRIFT ", discreteDRIFT, "  discreteDIFFUSION ", discreteDIFFUSION, "  sasymDIFFUSION ", sasymDIFFUSION, 
            "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        }
        if(verbose > 2) print("ukfstates ", ukfstates, "  ukfmeasures ", ukfmeasures);
  
        if(size(o0) > 0){
          int tmpindex[ncont_y[rowi]] = o0;
          for(oi in 1:ncont_y[rowi]) tmpindex[oi] +=  nmanifest;
           if(intoverstates==1) ypredcov_sqrt[o0,o0]=cholesky_decompose(makesym(ypredcov[o0,o0],verbose,1));
           //print("ypredcov = ", ypredcov, " err = ", err);
           //ll+=  multi_normal_cholesky_lpdf(Y[rowi] | ypred, ypredcov_sqrt);
           kout[rowi,o0] = mdivide_left_tri_low(ypredcov_sqrt[o0,o0], err[o0]); //transform pred errors to standard normal dist and collect
           kout[rowi,tmpindex] = log(diagonal(ypredcov_sqrt[o0,o0])); //account for transformation of scale in loglik
        }
      
    }//end nobs > 0 section
  if(savescores==1) kout[rowi,(nmanifest*4+nlatent+1):(nmanifest*4+nlatent+nlatent)] = eta[1:nlatent];
    } // end dokalmanrows subset selection
}//end rowi
if(dokalman==1){
  if(sum(nbinary_y) > 0) {
    vector[sum(nbinary_y)] binaryll;
    counter = 1;
    for(ri in 1:ndatapoints){
    if(dokalmanrows[ri]==1){
      int o1[nbinary_y[ri]] = whichbinary_y[ri,1:nbinary_y[ri]]; //which indicators are observed and binary
      binaryll[counter:(counter + nbinary_y[ri]-1)] = kout[ri,o1];
      counter+= nbinary_y[ri];
    }
    }
    ll+= sum(log(binaryll[1:(counter-1)]));
  }

  if(( sum(ncont_y) > 0)) {
    vector[sum(ncont_y)] errtrans[2];
    counter = 1;
    for(ri in 1:ndatapoints){
    if(dokalmanrows[ri]==1){
      int o0[ncont_y[ri]] = whichcont_y[ri,1:ncont_y[ri]]; //which indicators are observed and continuous
      errtrans[1,counter:(counter + ncont_y[ri]-1)] = kout[ri, o0];
      for(oi in 1:ncont_y[ri]) o0[oi] +=  nmanifest; //modify o0 index
      errtrans[2,counter:(counter + ncont_y[ri]-1)] = kout[ri, o0];
      counter+= ncont_y[ri];
    }
    }
    ll += normal_lpdf(errtrans[1,1:(counter-1)]|0,1) - sum(errtrans[2,1:(counter-1)]);
  }
if(savescores) kalman = kout;

    }
  }
}
      
model{
  if(intoverpop==0 && fixedsubpars == 1) target+= multi_normal_cholesky_lpdf(fixedindparams | rep_vector(0,nindvarying),IIindvar);

  if(nopriors==0){
   target+= dokalmanpriormodifier * normal_lpdf(rawpopmeans|0,1);
  
    if(ntipred > 0){ 
      target+= dokalmanpriormodifier * normal_lpdf(tipredeffectparams| 0, tipredeffectscale);
      target+= normal_lpdf(tipredsimputed| 0, tipredsimputedscale); //consider better handling of this when using subset approach
      //target+= dokalmanpriormodifier * normal_lpdf(tipredglobalscalepar | 0-log(ntipred),log(square(ntipred)));
    }
    
    if(nindvarying > 0){
      if(nindvarying >1) target+= dokalmanpriormodifier * normal_lpdf(sqrtpcov | 0, 1);
      if(intoverpop==0 && fixedsubpars == 0) target+= multi_normal_cholesky_lpdf(baseindparams | rep_vector(0,nindvarying), IIindvar);
      target+= dokalmanpriormodifier * normal_lpdf(rawpopsdbase | 0,1);
    }
    //llp +=  log(dokalmanpriormodifier);
  } //end pop priors section
  
  if(intoverstates==0) target+= normal_lpdf(etaupdbasestates|0,1);
  
  target+= ll; 

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
      if(popsetup[ri,3] && popsetup[ri,8]==0) { //if a free parameter //or state dependent parameter maybe?

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
