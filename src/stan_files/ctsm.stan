
functions{

 int[] vecequals(int[] a, int test, int comparison){ //do indices of a match test condition?
    int check[size(a)];
    for(i in 1:size(check)){
      if(comparison) check[i] = (test==a[i]) ? 1 : 0;
      if(comparison==0) check[i] = (test==a[i]) ? 0 :1;
    }
    return(check);
  }

int[] whichequals(int[] b, int test, int comparison){  //return array of indices of b matching test condition
    int bsize = size(b);
    int check[bsize] = vecequals(b,test,comparison);
    int whichsize = sum(check);
    int which[whichsize];
    int counter = 1;
    if(bsize > 0){
    for(i in 1:bsize){
      if(check[i] == 1){
        which[counter] = i;
        counter += 1;
      }
    }
    }
    return(which);
  }

  //matrix solvesyl(matrix A, matrix C); //using form F and -V from Wahlstrom Axelsson Gustafsson 2014

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
    if(size(z0) > 0){
      out[z0,] = rep_matrix(0,size(z0),rows(M));
      out[,z0] = rep_matrix(0,rows(M),size(z0));
      for(i in 1:size(z0)) out[z0[i],z0[i]] = exp(M[z0[i],z0[i]]);
    }
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
      s[i] = inv_sqrt(o[i,] * o[,i] + 1e-10);
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

  matrix kronsum(matrix mata, matrix II){
    return sqkron_prod(mata, II) + sqkron_prod(II, mata );
  }

  matrix makesym(matrix mat, int verbose, int pd){
    matrix[rows(mat),cols(mat)] out;
    for(coli in 1:cols(mat)){
      out[coli,coli] = mat[coli,coli]; 
      if(pd ==1 && out[coli,coli] < 1e-5){
        //if(verbose > 0) print("diagonal too low (",mat[coli,coli],") during makesym row ", coli, " col ", coli);
        out[coli,coli] = 1e-5 + fmin(0.0,out[coli,coli]);
      }  
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
  int statedependence[4];
  int multiplicativenoise;
  int choleskymats;
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
  int nrowmatsetup;
  int matsetup[nrowmatsetup,9];
  real matvalues[nrowmatsetup,6];
  int matrixdims[10,2];
  int savescores;
  int savesubjectmatrices;
  int fixedsubpars;
  vector[fixedsubpars ? nindvarying : 0] fixedindparams[fixedsubpars ? nsubjects : 0];
  int dokalman;
  int dokalmanrows[ndatapoints];
  real Jstep;
  real dokalmanpriormodifier;
  int intoverpopindvaryingindex[intoverpop ? nindvarying : 0];
  int nsJAxfinite;
  int sJAxfinite[nsJAxfinite];
  int taylorheun;
  int difftype;
  int jacoffdiag[nlatentpop];
  int njacoffdiagindex;
  int jacoffdiagindex[njacoffdiagindex];
  int popcovn;
}
      
transformed data{
  matrix[nlatent,nlatent] IIlatent= diag_matrix(rep_vector(1,nlatent));
  matrix[nlatent*nlatent,nlatent*nlatent] IIlatent2 = diag_matrix(rep_vector(1,nlatent*nlatent));
  matrix[nlatentpop,nlatentpop] IIlatentpop = diag_matrix(rep_vector(1,nlatentpop));
  matrix[nindvarying,nindvarying] IIindvar = diag_matrix(rep_vector(1,nindvarying));
  vector[ndatapoints] dtsmall = rep_vector(0,ndatapoints);
  vector[ndatapoints] dt = rep_vector(0,ndatapoints);
  int T0check[ndatapoints] = rep_array(-1,ndatapoints);
  int dtchange[ndatapoints] = rep_array(-1,ndatapoints);
  vector[nlatentpop-nlatent] nlpzerovec = rep_vector(0,nlatentpop-nlatent);
  vector[nlatent+1] nlplusonezerovec = rep_vector(0,nlatent+1);
  
  { //dt calcs
    int si = 0;
    real prevdt;
    int T0checktemp=0;
    real dttemp=0;
    real timei=0.0;
    for(rowi in 1:(dokalman ? ndatapoints :1)){
      if(dokalmanrows[rowi] ==1) { //used for subset selection
      T0check[rowi] = ( (si == subject[rowi]) ? (T0checktemp + 1) : 0 ) ; //if same subject, add 1 to t0check, else set to 0
      T0checktemp = T0check[rowi];
      if(T0check[rowi] > 0) prevdt = dttemp;
      dt[rowi] = T0check[rowi] ? time[rowi] - timei : 0;
      dttemp=dt[rowi];
      timei = time[rowi]; //must come after dt!
      dtsmall[rowi] = dt[rowi] / ceil(dt[rowi] / maxtimestep);
      si = subject[rowi];
      if(T0check[rowi] >0)  dtchange[rowi] = ( (prevdt-dt[rowi]) == 0.0) ? 0 : 1;
      
      }
    }
  }
  
}
      
parameters{
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
  matrix[nindvarying, nindvarying] rawpopcovchol; 
  matrix[nindvarying,nindvarying] rawpopcorr;
  matrix[nindvarying,nindvarying] rawpopcov;

  real ll = 0;
  vector[1] llrow[savescores ? ndatapoints : 0] = rep_array(rep_vector(0.0,1),savescores ? ndatapoints : 0);
  matrix[nlatentpop,nlatentpop] etapriorcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etaupdcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etasmoothcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ypriorcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] yupdcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ysmoothcov[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaprior[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaupd[savescores ? ndatapoints : 0];
  vector[nlatentpop] etasmooth[savescores ? ndatapoints : 0];
  vector[nmanifest] yprior[savescores ? ndatapoints : 0];
  vector[nmanifest] yupd[savescores ? ndatapoints : 0];
  vector[nmanifest] ysmooth[savescores ? ndatapoints : 0];
     matrix[matrixdims[1, 1], matrixdims[1, 2] ] T0MEANS[T0MEANSsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1]; 
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] LAMBDA[LAMBDAsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1]; 
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] DRIFT[DRIFTsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1]; 
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] DIFFUSION[DIFFUSIONsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1]; 
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] MANIFESTVAR[MANIFESTVARsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1]; 
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] MANIFESTMEANS[MANIFESTMEANSsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1]; 
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] CINT[CINTsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1]; 
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] T0VAR[T0VARsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1]; 
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] TDPREDEFFECT[TDPREDEFFECTsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1]; 
      matrix[matrixdims[10, 1], matrixdims[10, 2] ] PARS[PARSsubindex  ? (savesubjectmatrices ? nsubjects : 1) : 1];

  matrix[nlatent,nlatent] asymDIFFUSION[asymDIFFUSIONsubindex ? (savesubjectmatrices ? nsubjects : 1) : 1]; //stationary latent process variance
  vector[nlatent] asymCINT[asymCINTsubindex ? (savesubjectmatrices ? nsubjects : 1) : 1]; // latent process asymptotic level
matrix[nlatent, nlatent] DIFFUSIONcov[DIFFUSIONcovsubindex ? (savesubjectmatrices ? nsubjects : 1) : 1];
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
  vector[nlatent] pop_asymCINT; // latent process asymptotic level
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
    rawpopsd = log1p(exp(2*rawpopsdbase-1)) .* sdscale + 1e-10; // sqrts of proportions of total variance
    for(j in 1:nindvarying){
      rawpopcovsqrt[j,j] = rawpopsd[j]; //used with intoverpop
      for(i in 1:nindvarying){
        if(i > j){
          counter += 1;
          rawpopcovsqrt[i,j]=sqrtpcov[counter];
          rawpopcovsqrt[j,i]=0;//sqrtpcov[counter];
        }
      }
    }
    rawpopcorr = tcrossprod( constraincorsqrt(rawpopcovsqrt));
    rawpopcov = makesym(quad_form_diag(rawpopcorr, rawpopsd +1e-8),verbose,1);
    rawpopcovchol = cholesky_decompose(rawpopcov); 
  }//end indvarying par setup

  {
  int si = 0;
  int subjectcount = 0;
  int counter = 1;
  matrix[nlatentpop, nlatentpop] etacov; //covariance of latent states

  //measurement 
  vector[nmanifest] err;
  vector[sum(ncont_y)] errtrans = rep_vector(0, sum(ncont_y)); //to collect normalised errors
  vector[nmanifest] syprior;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypriorcov_sqrt; 
  matrix[nmanifest, nmanifest] ycov; 
  
  matrix[nlatentpop,nlatentpop] Je[savescores ? ndatapoints : 1]; //time evolved jacobian, saved for smoother
  matrix[nlatent*2,nlatent*2] dQi; //covariance from jacobian

  vector[nlatentpop] state = rep_vector(-1,nlatentpop); 
  matrix[nlatentpop,nlatentpop] sJAx; //Jacobian for drift
  matrix[nlatentpop,nlatentpop] sJ0; //Jacobian for t0
  matrix[nlatentpop,nlatentpop] sJtd;//diag_matrix(rep_vector(1),nlatentpop); //Jacobian for nltdpredeffect
  matrix[ nmanifest,nlatentpop] Jy[ndatapoints];//store Jacobian for measurement over time
  matrix[ nmanifest,nlatentpop] sJy;//Jacobian for measurement 
  matrix[nmanifest,nlatent] tLAMBDA[ndatapoints]; // store lambda time varying for smoother

  //linear continuous time calcs
  matrix[nlatent+1,nlatent+1] discreteDRIFT;
  matrix[nlatent,nlatent] discreteDIFFUSION = rep_matrix(0.0,nlatent,nlatent);

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
  vector[nlatent] sasymCINT; // latent process asymptotic level
matrix[nlatent, nlatent] sDIFFUSIONcov;

  for(rowi in 1:(dokalman ? ndatapoints :1)){
  if(dokalmanrows[rowi] ==1) { //used for subset selection
    si = subject[rowi];

    if(savescores || rowi==1) Je[savescores ? rowi : 1] = IIlatentpop; //elements updated later
       
    if(T0check[rowi] == 0) { // calculate initial matrices if this is first row for si
  
    
 int subjectvec[subjectcount ? 1 : 2];
 vector[nparams] rawindparams;
 vector[nparams] tipredaddition = rep_vector(0,nparams);
 vector[nparams] indvaraddition = rep_vector(0,nparams);
 subjectvec[size(subjectvec)] = si;
 if(subjectcount == 0)  subjectvec[1] = 0; // only needed for subject 0 (pop pars)
 subjectcount = subjectcount + 1;
 for(subjectveci in 1:size(subjectvec)){
  int subi = subjectvec[subjectveci];

  if(subi > 0 && nindvarying > 0 && intoverpop==0) {
    if(fixedsubpars==0) indvaraddition[indvaryingindex] = rawpopcovchol * baseindparams[subi];
    if(fixedsubpars==1) indvaraddition[indvaryingindex] = rawpopcovchol * fixedindparams[subi];
  }
  
  if(subi > 0 &&  ntipred > 0) tipredaddition = TIPREDEFFECT * tipreds[subi]';

  rawindparams = rawpopmeans + tipredaddition + indvaraddition;

    for(ri in 1:size(matsetup)){ //for each row of matrix setup
        for(statecalcs in 0:1){
        if(subi ==0 ||  //if population parameter
          ( matsetup[ri,7] == 8 && T0VARsubindex) || //or a covariance parameter in an individually varying matrix
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
      if(matsetup[ri, 7] == 54) sJy[matsetup[ ri,1], matsetup[ri,2]] = newval;
                if(matsetup[ri,9] < 0){ //then send copies elsewhere
                for(ri2 in 1:size(matsetup)){
                  if(matsetup[ri2,9] == ri){ 
                  if(matsetup[ri2, 7] == 1) sT0MEANS[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 2) sLAMBDA[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 3) sDRIFT[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 4) sDIFFUSION[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 5) sMANIFESTVAR[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 6) sMANIFESTMEANS[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 7) sCINT[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 8) sT0VAR[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 9) sTDPREDEFFECT[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 10) sPARS[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 51) sJ0[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 52) sJAx[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 53) sJtd[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 54) sJy[matsetup[ri2,1], matsetup[ri2,2]] = newval;
                  }
                }
              }
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
    sDIFFUSIONcov = sdcovsqrt2cov(sDIFFUSION,choleskymats);
  }
  if(subi <= (asymDIFFUSIONsubindex ? nsubjects : 0)) {
    if(ndiffusion < nlatent) sasymDIFFUSION = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

    if(continuoustime==1) sasymDIFFUSION[ derrind, derrind] = to_matrix( 
    -kronsum(sDRIFT[ derrind, derrind ],IIlatentpop[derrind,derrind]) \  to_vector( 
         sDIFFUSIONcov[ derrind, derrind ]), ndiffusion,ndiffusion);

    if(continuoustime==0) sasymDIFFUSION[ derrind, derrind ] = to_matrix( (IIlatent2 - 
      sqkron_prod(sDRIFT[ derrind, derrind ], sDRIFT[ derrind, derrind ])) \  to_vector(sDIFFUSIONcov[ derrind, derrind ]), ndiffusion, ndiffusion);
  } //end asymdiffusion loops

     if(subi <= (T0VARsubindex ? nsubjects : 0)) {
     if(intoverpop) sT0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopcovsqrt;
      sT0VAR = makesym(sdcovsqrt2cov(sT0VAR,choleskymats),verbose,1); 
      if(nt0varstationary > 0) {
        for(ri in 1:nt0varstationary){ 
          sT0VAR[t0varstationary[ri,1],t0varstationary[ri,2] ] =  sasymDIFFUSION[t0varstationary[ri,1],t0varstationary[ri,2] ];
        }
      }
      if(intoverpop){ //adjust cov matrix for transforms
        for(ri in 1:size(matsetup)){
          if(matsetup[ri,7]==1){ //if t0means
            if(matsetup[ri,5]) { //and indvarying
              sT0VAR[matsetup[ri,1], ] = sT0VAR[matsetup[ri,1], ] * matvalues[ri,2] * matvalues[ri,3]* matvalues[ri,5]; //multiplier meanscale sdscale
              sT0VAR[, matsetup[ri,1] ] = sT0VAR[, matsetup[ri,1] ] * matvalues[ri,2] * matvalues[ri,3]* matvalues[ri,5]; //multiplier meanscale sdscale
            }
          }
        }
      }
    }
    
    if(subi <= (asymCINTsubindex ? nsubjects : 0)){
      if(continuoustime==1) sasymCINT =  -sDRIFT[1:nlatent,1:nlatent] \ sCINT[ ,1 ];
      if(continuoustime==0) sasymCINT =  (IIlatent - sDRIFT[1:nlatent,1:nlatent]) \ sCINT[,1 ];
    }
    
    if(nt0meansstationary > 0){
      if(subi <= (T0MEANSsubindex ? nsubjects : 0)) {
        for(ri in 1:nt0meansstationary){
          sT0MEANS[t0meansstationary[ri,1] , 1] = 
            sasymCINT[t0meansstationary[ri,1] ];
        }
      }
    }
  if(subi == 0 || savesubjectmatrices){ 
if( (T0MEANSsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (T0MEANSsubindex == 0 && subi==0) ) T0MEANS[(savesubjectmatrices && T0MEANSsubindex) ? subi : 1] = sT0MEANS; 
if( (LAMBDAsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (LAMBDAsubindex == 0 && subi==0) ) LAMBDA[(savesubjectmatrices && LAMBDAsubindex) ? subi : 1] = sLAMBDA; 
if( (DRIFTsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (DRIFTsubindex == 0 && subi==0) ) DRIFT[(savesubjectmatrices && DRIFTsubindex) ? subi : 1] = sDRIFT; 
if( (DIFFUSIONsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (DIFFUSIONsubindex == 0 && subi==0) ) DIFFUSION[(savesubjectmatrices && DIFFUSIONsubindex) ? subi : 1] = sDIFFUSION; 
if( (MANIFESTVARsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (MANIFESTVARsubindex == 0 && subi==0) ) MANIFESTVAR[(savesubjectmatrices && MANIFESTVARsubindex) ? subi : 1] = sMANIFESTVAR; 
if( (MANIFESTMEANSsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (MANIFESTMEANSsubindex == 0 && subi==0) ) MANIFESTMEANS[(savesubjectmatrices && MANIFESTMEANSsubindex) ? subi : 1] = sMANIFESTMEANS; 
if( (CINTsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (CINTsubindex == 0 && subi==0) ) CINT[(savesubjectmatrices && CINTsubindex) ? subi : 1] = sCINT; 
if( (T0VARsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (T0VARsubindex == 0 && subi==0) ) T0VAR[(savesubjectmatrices && T0VARsubindex) ? subi : 1] = sT0VAR; 
if( (TDPREDEFFECTsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (TDPREDEFFECTsubindex == 0 && subi==0) ) TDPREDEFFECT[(savesubjectmatrices && TDPREDEFFECTsubindex) ? subi : 1] = sTDPREDEFFECT; 
if( (PARSsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (PARSsubindex == 0 && subi==0) ) PARS[(savesubjectmatrices && PARSsubindex) ? subi : 1] = sPARS; 
if( (DIFFUSIONcovsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (DIFFUSIONcovsubindex == 0 && subi==0) ) DIFFUSIONcov[(savesubjectmatrices && DIFFUSIONcovsubindex) ? subi : 1] = sDIFFUSIONcov; 
if( (asymDIFFUSIONsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (asymDIFFUSIONsubindex == 0 && subi==0) ) asymDIFFUSION[(savesubjectmatrices && asymDIFFUSIONsubindex) ? subi : 1] = sasymDIFFUSION; 
if( (asymCINTsubindex > 0 && (subi > 0 || savesubjectmatrices==0) ) || (asymCINTsubindex == 0 && subi==0) ) asymCINT[(savesubjectmatrices && asymCINTsubindex) ? subi : 1] = sasymCINT; 
 
 }
  if(subi == 0){
pop_T0MEANS = sT0MEANS; pop_LAMBDA = sLAMBDA; pop_DRIFT = sDRIFT; pop_DIFFUSION = sDIFFUSION; pop_MANIFESTVAR = sMANIFESTVAR; pop_MANIFESTMEANS = sMANIFESTMEANS; pop_CINT = sCINT; pop_T0VAR = sT0VAR; pop_TDPREDEFFECT = sTDPREDEFFECT; pop_PARS = sPARS; pop_DIFFUSIONcov = sDIFFUSIONcov; pop_asymDIFFUSION = sasymDIFFUSION; pop_asymCINT = sasymCINT; 
  }

} // end subject matrix creation
  

    etacov =  sT0VAR;
    state = sT0MEANS[,1]; //init and in case of jacobian dependencies

      if(nldynamics==0){ //initialize most parts for nl later!
        if(ntdpred > 0) state[1:nlatent] += sTDPREDEFFECT * tdpreds[rowi];
      }
    } //end T0 matrices
if(verbose > 1) print ("below t0 row ", rowi);

    if(nldynamics==0 && T0check[rowi]>0){ //linear kf time update
      if(verbose > 1) print ("linear update row ", rowi);
    
      if(continuoustime ==1){
        if(dtchange[rowi]==1 || (T0check[rowi] == 1 && (DRIFTsubindex + CINTsubindex > 0))){ //if dtchanged or if subject variability
          discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent,1:nlatent],sCINT),nlplusonezerovec') * dt[rowi],drcintoffdiag);
          if(!savescores) Je[1, 1:nlatent, 1:nlatent] = discreteDRIFT[1:nlatent,1:nlatent];
        }
      
        if(dtchange[rowi]==1 || (T0check[rowi] == 1 && (DIFFUSIONsubindex + DRIFTsubindex > 0))){ //if dtchanged or if subject variability
          discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
            quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]' );
          if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(makesym(discreteDIFFUSION,verbose,1));
        }
      }

      if(continuoustime==0 && T0check[rowi] == 1){
        if(subjectcount == 1 || DIFFUSIONsubindex + DRIFTsubindex + CINTsubindex > 0){ //if first subject or variability
          discreteDRIFT=append_row(append_col(sDRIFT[1:nlatent,1:nlatent],sCINT),nlplusonezerovec');
          discreteDRIFT[nlatent+1,nlatent+1] = 1;
          if(!savescores) Je[1, 1:nlatent, 1:nlatent] = discreteDRIFT[1:nlatent,1:nlatent];
          discreteDIFFUSION=sDIFFUSIONcov;
          if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(makesym(discreteDIFFUSION,verbose,1));
        }
      }
      if(savescores) Je[rowi, 1:nlatent, 1:nlatent] = discreteDRIFT[1:nlatent,1:nlatent];
      state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent];
      if(ntdpred > 0) state[1:nlatent] += sTDPREDEFFECT * tdpreds[rowi];
      if(intoverstates==1) {
        etacov = quad_form(etacov, Je[savescores ? rowi : 1]');
        if(ndiffusion > 0) etacov[1:nlatent,1:nlatent] += discreteDIFFUSION;

      }
    }//end linear time update


    if(nldynamics==1){ //nldynamics time update
      if(T0check[rowi]>0){
        vector[nlatentpop] base;
        real intstepi = 0;
        
        while(intstepi < (dt[rowi]-1e-10)){
          intstepi = intstepi + dtsmall[rowi];
          
    {
    int zeroint[1];
    vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(sJAxfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0)  state[statei] += Jstep;
      
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == 2 &&(
            matsetup[ri,7] ==3||
            matsetup[ri,7] ==7 )){ //perform calcs appropriate to this section
              real newval;
              newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
              if(matsetup[ri, 7] == 3) sDRIFT[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 7) sCINT[matsetup[ ri,1], matsetup[ri,2]] = newval;
              if(matsetup[ri,9] < 0){
                for(ri2 in 1:size(matsetup)){
                  if(matsetup[ri2,9] == ri){ //if row ri2 is a copy of original row ri
                    if(matsetup[ri2, 7] == 3) sDRIFT[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 7) sCINT[matsetup[ri2,1], matsetup[ri2,2]] = newval;
                  }
                }
              }
            }
          }
          {
  ;  
  
  } 
   
      if(statei > 0) {
        sJAx[sJAxfinite,statei] =  sDRIFT[sJAxfinite, ] * state + append_row(sCINT[,1],nlpzerovec)[sJAxfinite]; //compute new change
         if(verbose>1) print("sJAx ",sJAx);
      }
      if(statei== 0 && size(sJAxfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
        base[sJAxfinite] = sDRIFT[sJAxfinite, ] * state + append_row(sCINT[,1],nlpzerovec)[sJAxfinite];
        if(verbose>1) print("base = ",base,"    sjaxinit= ",sJAx);
        for(fi in sJAxfinite){
          sJAx[sJAxfinite,fi] = (sJAx[sJAxfinite,fi] - base[sJAxfinite]) / Jstep; //new - baseline change divided by stepsize
        }
      }
    }
    if(verbose>1) print("sJAx ",sJAx);
    }
    
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == 2 &&(
            matsetup[ri,7] ==4||
            matsetup[ri,7] ==52 )){ //perform calcs appropriate to this section
              real newval;
              newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
              if(matsetup[ri, 7] == 4) sDIFFUSION[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 52) sJAx[matsetup[ ri,1], matsetup[ri,2]] = newval;
              if(matsetup[ri,9] < 0){
                for(ri2 in 1:size(matsetup)){
                  if(matsetup[ri2,9] == ri){ //if row ri2 is a copy of original row ri
                    if(matsetup[ri2, 7] == 4) sDIFFUSION[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 52) sJAx[matsetup[ri2,1], matsetup[ri2,2]] = newval;
                  }
                }
              }
            }
          }    {
  ;  
  
  } 
   
      
      if(multiplicativenoise) sDIFFUSIONcov[derrind,derrind] = sdcovsqrt2cov(sDIFFUSION[derrind,derrind],choleskymats);
      
        if(continuoustime){
          if(taylorheun==0){
            if(dtchange[rowi]==1 || statedependence[2] || 
              (T0check[rowi] == 1 && (DRIFTsubindex + DIFFUSIONsubindex + CINTsubindex) > 0)){
                
              if(difftype==2){
                matrix[ndiffusion*2,ndiffusion*2] ebA;
                matrix[ndiffusion*2,ndiffusion*2] bA;
            
                bA[1:ndiffusion,1:ndiffusion] = -sJAx[derrind,derrind];
                bA[1:ndiffusion,(1+ndiffusion):(ndiffusion*2)] = sDIFFUSIONcov[derrind,derrind];
                bA[(1+ndiffusion):(ndiffusion*2),(1+ndiffusion):(ndiffusion*2)] = sJAx[derrind,derrind]';
                bA[(1+ndiffusion):(ndiffusion*2),1:ndiffusion] = rep_matrix(0,ndiffusion,ndiffusion);
                
                ebA = matrix_exp(bA * dtsmall[rowi]);
                Je[savescores ? rowi : 1,derrind,derrind] =  ebA[(1+ndiffusion):(ndiffusion*2),(1+ndiffusion):(ndiffusion*2)]';
                discreteDIFFUSION[derrind,derrind] = Je[savescores ? rowi : 1,derrind,derrind] * ebA[1:ndiffusion,(1+ndiffusion):(ndiffusion*2)];
                discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec') * dtsmall[rowi],drcintoffdiag);
                if(ndiffusion < nlatentpop) Je[savescores ? rowi : 1] =  expm2(sJAx * dtsmall[rowi], jacoffdiag);
              } else if(savescores) Je[rowi] = Je[rowi-1];
                
              //if(difftype==1){
              //matrix[nlatent,nlatent] V = sDIFFUSIONcov-quad_form(sDIFFUSIONcov, Je[savescores ? rowi : 1]');
              //discreteDIFFUSION = solvesyl(sJAx[1:nlatent,1:nlatent],-V,discreteDIFFUSION, rep_array(nlatent,1));
              //}
              
              //sasymDIFFUSION[derrind,derrind] = to_matrix(  -kronsum(sJAx[derrind,derrind],IIlatentpop[derrind,derrind]) \ to_vector(sDIFFUSIONcov[derrind,derrind]), ndiffusion,ndiffusion);
              //discreteDIFFUSION[derrind,derrind] =  sasymDIFFUSION[derrind,derrind] - quad_form( sasymDIFFUSION[derrind,derrind], Je[savescores ? rowi : 1, derrind,derrind]' );
            }
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent]; // ???compute before new diffusion calcs
            etacov = quad_form(etacov, Je[savescores ? rowi : 1]');
            etacov[derrind,derrind] += discreteDIFFUSION[derrind,derrind]; 
          }
            
          if(taylorheun==1){
            matrix[nlatentpop,nlatentpop] Kth = inverse(IIlatentpop - sJAx * (dtsmall[rowi] /2) );
            matrix[nlatentpop,nlatentpop] Mth = Kth * (IIlatentpop + sJAx * (dtsmall[rowi] /2) );
            state[1:nlatent] = state[1:nlatent] + Kth[1:nlatent,1:nlatent] *
              (sDRIFT[1:nlatent,1:nlatent] * state[1:nlatent] + sCINT[1:nlatent,1]) * dtsmall[rowi];
            etacov = quad_form(etacov, Mth');
            etacov[derrind,derrind] += quad_form(Kth[derrind,derrind],sDIFFUSIONcov[derrind,derrind]) * dtsmall[rowi];
          }
            if(intstepi >= (dt[rowi]-1e-10) && savescores) Je[rowi] = expm2(sJAx * dt[rowi],jacoffdiag); //save approximate exponentiated jacobian for smoothing
          }
  
          if(continuoustime==0){ 
            Je[savescores ? rowi : 1] = sJAx;
            etacov = quad_form(etacov, sJAx');
            etacov[ derrind, derrind ] += tcrossprod(sDIFFUSION[ derrind, derrind ]); //may need improving re sDIFFUSION
            discreteDRIFT=append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec');
            discreteDRIFT[nlatent+1,nlatent+1] = 1;
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent];
          }
      }
    } // end nonlinear time update
  
    if(T0check[rowi]==0){ //nl t0
    state = sT0MEANS[,1]; //in case of t0 dependencies, may have missingness
    
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == 1 &&(
            matsetup[ri,7] ==1||
            matsetup[ri,7] ==8||
            matsetup[ri,7] ==51 )){ //perform calcs appropriate to this section
              real newval;
              newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
              if(matsetup[ri, 7] == 1) sT0MEANS[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 8) sT0VAR[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 51) sJ0[matsetup[ ri,1], matsetup[ri,2]] = newval;
              if(matsetup[ri,9] < 0){
                for(ri2 in 1:size(matsetup)){
                  if(matsetup[ri2,9] == ri){ //if row ri2 is a copy of original row ri
                    if(matsetup[ri2, 7] == 1) sT0MEANS[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 8) sT0VAR[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 51) sJ0[matsetup[ri2,1], matsetup[ri2,2]] = newval;
                  }
                }
              }
            }
          }
    ;
      state = sT0MEANS[,1];
      etacov = quad_form(sT0VAR, sJ0');
    }//end nonlinear t0
    
    if(ntdpred > 0) {
      int nonzerotdpred = 0;
      for(tdi in 1:ntdpred) if(tdpreds[rowi,tdi] != 0.0) nonzerotdpred = 1;
      if(nonzerotdpred){
      
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == 3 &&(
            matsetup[ri,7] ==9||
            matsetup[ri,7] ==53 )){ //perform calcs appropriate to this section
              real newval;
              newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
              if(matsetup[ri, 7] == 9) sTDPREDEFFECT[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 53) sJtd[matsetup[ ri,1], matsetup[ri,2]] = newval;
              if(matsetup[ri,9] < 0){
                for(ri2 in 1:size(matsetup)){
                  if(matsetup[ri2,9] == ri){ //if row ri2 is a copy of original row ri
                    if(matsetup[ri2, 7] == 9) sTDPREDEFFECT[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 53) sJtd[matsetup[ri2,1], matsetup[ri2,2]] = newval;
                  }
                }
              }
            }
          }
      ;
        state[1:nlatent] +=   (sTDPREDEFFECT * tdpreds[rowi]); //tdpred effect only influences at observed time point
        etacov = quad_form(etacov,sJtd'); //could be optimized
      }
    }//end nonlinear tdpred
  } // end non linear time update

  if(intoverstates==0 && nldynamics == 0) {
    if(T0check[rowi]==0) state += cholesky_decompose(sT0VAR) * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)];
    if(T0check[rowi]>0) state +=  discreteDIFFUSION * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)];
  }

if(verbose > 1){
  print("etaprior = ", state);
  print("etapriorcov = ", etacov);
}

if(savescores){
  etapriorcov[rowi] = etacov; 
  etaprior[rowi] = state;
}

 if (nobs_y[rowi] > 0 || savescores) {  // if some observations create right size matrices for missingness and calculate...
    
      int o[savescores ? nmanifest : nobs_y[rowi]]; //which obs are not missing in this row
      int o1[savescores ? size(whichequals(manifesttype,1,1)) : nbinary_y[rowi] ];
      int o0[savescores ? size(whichequals(manifesttype,1,0)) : ncont_y[rowi] ];
      
      int od[nobs_y[rowi]] = whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
      int o1d[nbinary_y[rowi] ]= whichbinary_y[rowi,1:nbinary_y[rowi]];
      int o0d[ncont_y[rowi] ]= whichcont_y[rowi,1:ncont_y[rowi]];
      
      if(!savescores){
        o= whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
        o1= whichbinary_y[rowi,1:nbinary_y[rowi]];
        o0= whichcont_y[rowi,1:ncont_y[rowi]];
      }
      if(savescores){ //needed to calculate yprior and yupd ysmooth
        for(mi in 1:nmanifest) o[mi] = mi;
        o1= whichequals(manifesttype,1,1);
        o0= whichequals(manifesttype,1,0);
      }
      
      if(nlmeasurement==1){
      
          for(ri in 1:size(matsetup)){ //for each row of matrix setup
            if(matsetup[ri,3] > 0 && matsetup[ri,8] == 4 &&(
            matsetup[ri,7] ==2||
            matsetup[ri,7] ==6||
            matsetup[ri,7] ==5||
            matsetup[ri,7] ==54 )){ //perform calcs appropriate to this section
              real newval;
              newval = tform(state[ matsetup[ri,3] ], matsetup[ri,4], matvalues[ri,2], matvalues[ri,3], matvalues[ri,4], matvalues[ri,6] ); 
              if(matsetup[ri, 7] == 2) sLAMBDA[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 5) sMANIFESTVAR[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 6) sMANIFESTMEANS[matsetup[ ri,1], matsetup[ri,2]] = newval; 
      if(matsetup[ri, 7] == 54) sJy[matsetup[ ri,1], matsetup[ri,2]] = newval;
              if(matsetup[ri,9] < 0){
                for(ri2 in 1:size(matsetup)){
                  if(matsetup[ri2,9] == ri){ //if row ri2 is a copy of original row ri
                    if(matsetup[ri2, 7] == 2) sLAMBDA[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 5) sMANIFESTVAR[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 6) sMANIFESTMEANS[matsetup[ri2,1], matsetup[ri2,2]] = newval; 
      if(matsetup[ri2, 7] == 54) sJy[matsetup[ri2,1], matsetup[ri2,2]] = newval;
                  }
                }
              }
            }
          }
      ;
      }
      
      syprior[o] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * state[1:nlatent];
      if(intoverstates==1) { //classic kalman
        if(nbinary_y[rowi] > 0) syprior[o1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * state[1:nlatent])));
        ycov[o,o] = quad_form(etacov, sJy[o,]'); // + sMANIFESTVAR[o,o]; shifted measurement error down
        for(wi in 1:nmanifest){ 
          if(Y[rowi,wi] != 99999 || savescores==1) ycov[wi,wi] += square(sMANIFESTVAR[wi,wi]);
          if(manifesttype[wi]==1 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += fabs((syprior[wi] - 1) .* (syprior[wi]));
          if(manifesttype[wi]==2 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += square(fabs((syprior[wi] - round(syprior[wi])))); 
        }
      }
        
      if(intoverstates==0) { //sampled states
        if(ncont_y[rowi] > 0) ypriorcov_sqrt[o0,o0] = sMANIFESTVAR[o0,o0]+1e-8;
        if(nbinary_y[rowi] > 0) syprior[o1] = to_vector(inv_logit(to_array_1d(syprior[o1])));
      }
        
     
err[od] = Y[rowi,od] - syprior[od]; // prediction error
    
      if(intoverstates==1 && size(od) > 0) {
        K[,od] = mdivide_right(etacov * sJy[od,]', ycov[od,od]); 
        etacov += -K[,od] * sJy[od,] * etacov;
        state +=  (K[,od] * err[od]);
      }
      
      if(savescores==1) {
        yprior[rowi] = syprior[o];
        etaupd[rowi] = state;
        ypriorcov[rowi] = ycov;
        etaupdcov[rowi] = etacov;
        yupdcov[rowi] = quad_form(etacov, sJy') + sMANIFESTVAR;
        yupd[rowi] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * state[1:nlatent];
        ysmoothcov[rowi] = sMANIFESTVAR; // add the rest later
        ysmooth[rowi] = sMANIFESTMEANS[,1]; // add the rest later
        Jy[rowi] = sJy;
        tLAMBDA[rowi] = sLAMBDA;
      }
      
      
      if(verbose > 1) {
          print("rowi =",rowi, "  si =", si, "  state =",state,"  etacov ",etacov,
            "  syprior =",syprior,"  ycov ",ycov, "  K ",K,
            "  sDRIFT =", sDRIFT, " sDIFFUSION =", sDIFFUSION, " sCINT =", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
            "  sT0VAR", sT0VAR,  " sT0MEANS ", sT0MEANS, "sLAMBDA = ", sLAMBDA, "  sJy = ",sJy,
            " discreteDRIFT = ", discreteDRIFT, "  discreteDIFFUSION ", discreteDIFFUSION, "  sasymDIFFUSION ", sasymDIFFUSION, 
            " DIFFUSIONcov = ", sDIFFUSIONcov,
            "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        }
  
      if(nbinary_y[rowi] > 0) ll+= sum(log(Y[rowi,o1d] .* (syprior[o1d]) + (1-Y[rowi,o1d]) .* (1-syprior[o1d]))); 
  
        if(size(o0d) > 0){
           if(intoverstates==1) ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(makesym(ycov[o0d,o0d],verbose,1));
           if(savescores) llrow[rowi,1] =  multi_normal_cholesky_lpdf(Y[rowi,o0d] | syprior[o0d], ypriorcov_sqrt[o0d,o0d]);
           errtrans[counter:(counter + ncont_y[rowi]-1)] = mdivide_left_tri_low(ypriorcov_sqrt[o0d,o0d], err[o0d]); //transform pred errors to standard normal dist and collect
           ll+= -sum(log(diagonal(ypriorcov_sqrt[o0d,o0d]))); //account for transformation of scale in loglik
           counter += ncont_y[rowi];
        }
      
    }//end nobs > 0 section
  
  if(savescores && (rowi==ndatapoints || subject[rowi+1] != subject[rowi])){ //at subjects last datapoint, smooth
    int sri = rowi;
    while(sri>0 && subject[sri]==si){
      if(sri==rowi) {
        etasmooth[sri]=etaupd[sri];
        etasmoothcov[sri]=etaupdcov[sri];
      } else{
        matrix[nlatentpop,nlatentpop] smoother;
        smoother = etaupdcov[sri] * Je[sri+1]' / makesym(etapriorcov[sri+1],verbose,1);
        etasmooth[sri]= etaupd[sri] + smoother * (etasmooth[sri+1] - etaprior[sri+1]);
        etasmoothcov[sri]= etaupdcov[sri] + smoother * ( etasmoothcov[sri+1] - etapriorcov[sri+1]) * smoother';

      }
      ysmoothcov[sri] += quad_form(etasmoothcov[sri], Jy[sri]'); //already added manifestvar
      ysmooth[sri] += tLAMBDA[sri] * etasmooth[sri,1:nlatent];
      sri += -1;
      while(sri > 0 && dokalmanrows[sri]==0) sri+= -1; //skip rows if requested
    }
  } //end smoother
  
  } // end dokalmanrows subset selection
}//end rowi
if(dokalman==1){
 
  if(( sum(ncont_y) > 0)) {
  if(verbose > 1) print("errtrans = ", errtrans);
    ll += normal_lpdf(errtrans|0,1);
  }

    }
  }
}
      
model{
  if(intoverpop==0 && fixedsubpars == 1) target+= multi_normal_cholesky_lpdf(fixedindparams | rep_vector(0,nindvarying),IIindvar);

    if(ntipred > 0){ 
      if(nopriors==0) target+= dokalmanpriormodifier * normal_lpdf(tipredeffectparams| 0, tipredeffectscale);
      target+= normal_lpdf(tipredsimputed| 0, tipredsimputedscale); //consider better handling of this when using subset approach
    }

  if(nopriors==0){ //if split files over subjects, just compute priors once
   target+= dokalmanpriormodifier * normal_lpdf(rawpopmeans|0,1);
  
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
  vector[nindvarying] popsd; // = rep_vector(0,nparams);
  matrix[nindvarying,nindvarying] popcov;
  matrix[nparams,ntipred] linearTIPREDEFFECT;


  {
    matrix[popcovn, nindvarying] x;
    if(nindvarying){
      for(ri in 1:rows(x)){
        x[ri,] = (rawpopcovchol * 
          to_vector(normal_rng(rep_vector(0,nindvarying),rep_vector(1,nindvarying))) + 
          rawpopmeans[indvaryingindex])';
      }
    }
    
    for(pi in 1:nparams){
      int found=0;
      int pr1;
      int pr2;
      real rawpoppar = rawpopmeans[pi];
      while(!found){
        for(ri in 1:size(matsetup)){
          if(matsetup[ri,9] <=0 && matsetup[ri,3]==pi && matsetup[ri,8]==0) { //if a free parameter 
            pr1 = ri; 
            pr2=ri;// unless intoverpop, pop matrix row reference is simply current row
            found=1;
            if(intoverpop && matsetup[ri,5]) { //check if shifted
              for(ri2 in 1:size(matsetup)){ //check when state reference param of matsetup corresponds to row of t0means in current matsetup row
                if(matsetup[ri2,8]  && matsetup[ri2,3] == matsetup[ri,1] && 
                matsetup[ri2,3] > nlatent && matsetup[ri2,7] < 20 &&
                matsetup[ri,9] <=0) pr2 = ri2; //if param is dynamic and matches row (state ref) and is not in jacobian
                //print("ri = ",ri, " pr2 = ",pr2, " ri2 = ",ri2);
              }
            }
          }
        }
      }
        
      popmeans[pi] = tform(rawpoppar, matsetup[pr2,4], matvalues[pr2,2], matvalues[pr2,3], matvalues[pr2,4], matvalues[pr2,6] ); 
      if(matsetup[pr1,5]){ //if indvarying, transform random sample
        for(ri in 1:rows(x)){
          x[ri,matsetup[pr1,5]] = tform(x[ri,matsetup[pr1,5]],matsetup[pr2,4],matvalues[pr2,2],matvalues[pr2,3],matvalues[pr2,4],matvalues[pr2,6]);
        }
        x[,matsetup[pr1,5]] += rep_vector(-mean(x[,matsetup[pr1,5]]),rows(x));
      }
      if(ntipred > 0){
      for(tij in 1:ntipred){
        if(TIPREDEFFECTsetup[matsetup[pr1,3],tij] ==0){
          linearTIPREDEFFECT[matsetup[pr1,3],tij] = 0;
        } else {
        linearTIPREDEFFECT[matsetup[pr1,3],tij] = ( //tipred reference is from row pr1, tform reference from row pr2 in case of intoverpop
          tform(rawpoppar + TIPREDEFFECT[matsetup[pr1,3],tij] * .01, matsetup[pr2,4], matvalues[pr2,2], matvalues[pr2,3], matvalues[pr2,4], matvalues[pr2,6] ) -
          tform(rawpoppar - TIPREDEFFECT[matsetup[pr1,3],tij] * .01, matsetup[pr2,4], matvalues[pr2,2], matvalues[pr2,3], matvalues[pr2,4], matvalues[pr2,6] )
          ) /2 * 100;
        }
      }
    }
    } //end nparams loop
  
  if(nindvarying){
    popcov = crossprod(x) /(rows(x)-1);
    popsd = sqrt(diagonal(popcov));
  }
  }



}
