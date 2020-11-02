
functions{

 int[] vecequals(int[] a, int test, int comparison){ //do indices of a match test condition?
    int check[size(a)];
    for(i in 1:size(check)) check[i] = comparison ? (test==a[i]) : (test!=a[i]);
    return(check);
  }

int[] whichequals(int[] b, int test, int comparison){  //return array of indices of b matching test condition
    int check[size(b)] = vecequals(b,test,comparison);
    int which[sum(check)];
    int counter = 1;
    if(size(b) > 0){
      for(i in 1:size(b)){
        if(check[i] == 1){
          which[counter] = i;
          counter += 1;
        }
      }
    }
    return(which);
  }

 
   matrix constraincorsqrt(matrix mat){ //converts from unconstrained lower tri matrix to cor
    matrix[rows(mat),cols(mat)] o;
  
    for(i in 1:rows(o)){ //set upper tri to lower
      for(j in min(i+1,rows(mat)):rows(mat)){
        o[j,i] =  inv_logit(mat[j,i])*2-1;  // can change cor prior here
        o[i,j] = o[j,i];
      }
      o[i,i]=1; // change to adjust prior for correlations
      o[i,] = o[i,] / sqrt(sum(square(o[i,]))+1e-10);
    }
    return o;
  } 

  matrix sdcovsqrt2cov(matrix mat, int cholbasis){ //covariance from cholesky or unconstrained cor sq root
    if(cholbasis==0)  {
      return(tcrossprod(diag_pre_multiply(diagonal(mat),constraincorsqrt(mat))));
    } else return(tcrossprod(mat));
  }

  matrix sqkron_prod(matrix mata, matrix matb){
    int d=rows(mata);
    matrix[d*d,d*d] out;
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

 matrix sqkron_sumii(matrix mata){
   int d=rows(mata);
   matrix[d*d,d*d] out;
     for (l in 1:d){
       for (k in 1:d){
         for (j in 1:d){
           for (i in 1:d){
             out[i+(k-1)*d,j+(l-1)*d] = 0;
             if(i==j) out[i+(k-1)*d,j+(l-1)*d] += mata[k,l];
             if(k==l) out[i+(k-1)*d,j+(l-1)*d] += mata[i,j];
           }
         }
       }
     }
   return(out);
 } 

  matrix makesym(matrix mat, int verbose, int pd){
    matrix[rows(mat),cols(mat)] out;
    for(coli in 1:cols(mat)){
    //  if(pd ==1 && mat[coli,coli] < 1e-5){
     //   out[coli,coli] = 1e-5;// 
     // } else 
      out[coli,coli] = mat[coli,coli] + 1e-6; 
      for(rowi in 1:rows(mat)){
        if(rowi > coli) {
          out[rowi,coli] = mat[rowi,coli];
          out[coli,rowi] = mat[rowi,coli];
        }
        
      }
    }
    return out;
  }

  real tform(real parin, int transform, data real multiplier, data real meanscale, data real offset, data real inneroffset){
    real param=parin;
    if(meanscale!=1.0) param *= meanscale; 
if(inneroffset != 0.0) param += inneroffset; 
if(transform==0) param = param;
if(transform==1) param = (log1p_exp(param));
if(transform==2) param = (exp(param));
if(transform==3) param = (1/(1+exp(-param)));
if(transform==4) param = ((param)^3);
if(transform==5) param = log1p(param);
if(transform==50) param = meanscale;
if(transform==51) param = 1/(1+exp(-param));
if(transform==52) param = exp(param);
if(transform==53) param = 1/(1+exp(-param))-(exp(param)^2)/(1+exp(param))^2;
if(transform==54) param = 3*param^2;
if(transform==55) param = 1/(1+param);

if(multiplier != 1.0) param *=multiplier;
if(transform < 49 && offset != 0.0) param+=offset;
    return param;
  }
  
   vector parvectform(int n, vector rawpar, int when, int[,] ms, data real[,] mval, int subi, int[] subindices, int[] whenvec){
    vector[n] parout;
    if(n>0){
      int outwhen[size(whichequals(whenvec,0,0))] = whenvec[whichequals(whenvec,0,0)]; //outwhen is nonzero elements of whenvec
      for(ri in 1:size(ms)){ //for each row of matrix setup
        if(ms[ri,8]==when && ms[ri,9] < 1 && ms[ri,3] > 0){ //if correct when,not a copyrow, and free parameter
          if(subi ==0 ||  //if population parameter
            ( ms[ri,7] == 8 && subindices[8]) || //or a covariance parameter in an individually varying matrix
            (ms[ri,3] > 0 && (ms[ri,5] > 0 || ms[ri,6] > 0 || ms[ri,8] > 0)) //or there is individual variation
            ){ //otherwise repeated values
            
            parout[whichequals(outwhen, ms[ri,3],1)[1] ] = tform(rawpar[ms[ri,3] ], //which outwhen refers to correct par
              ms[ri,4], mval[ri,2], mval[ri,3], mval[ri,4], mval[ri,6] ); 
           }
        }
      }
    }
  return parout;
  }
  
  
  matrix mcalc(matrix matin, vector tfpars, vector tfstates, int[] when, int m, int[,] ms, data real[,] mval, int subi, int[] subindices){
    matrix[rows(matin),cols(matin)] matout;

    for(ri in 1:size(ms)){ //for each row of matrix setup
      int whenyes = 0;
      for(wi in 1:size(when)) if(when[wi]==ms[ri,8]) whenyes = 1;
      if(m==ms[ri,7] && whenyes){ // if correct matrix and when

        if(subi ==0 ||  //if population parameter
          ( ms[ri,7] == 8 && subindices[8]) || //or a covariance parameter in an individually varying matrix
          (ms[ri,3] > 0 && (ms[ri,5] > 0 || ms[ri,6] > 0 || ms[ri,8] > 0)) //or there is individual variation
          ){ //otherwise repeated values (maybe this check not needed now?

          if(ms[ri,3] > 0 && ms[ri,8]==0)  matout[ms[ri,1], ms[ri,2] ] = tfpars[ms[ri,3]]; //should be already tformed
          if(ms[ri,3] > 0 && ms[ri,8]>0)  matout[ms[ri,1], ms[ri,2] ] = tfstates[ms[ri,3]]; //should be already tformed
          if(ms[ri,3] < 1) matout[ms[ri,1], ms[ri,2] ] = mval[ri, 1]; //doing this once over all subjects unless covariance matrix -- speed ups possible here, check properly!
        }
      }
    }
    for(ri in 1:rows(matin)){ //fill holes with unchanged input matrix
      for(ci in 1:cols(matin)){
        if(is_nan(matout[ri,ci]) && !is_nan(matin[ri,ci])) matout[ri,ci] = matin[ri,ci];
      }
    }
  return(matout);
  }
  
  
  int[] checkoffdiagzero(matrix M){
    int z[rows(M)];
    for(i in 1:rows(M)){
      z[i] = 0;
      for(j in 1:rows(M)){ // check cols and rows simultaneously
        if(i!=j && (M[i,j]!=0.0 || M[j,i] != 0.0)){
          z[i] = 1;
          break;
        }
      }
    }
    return z;
  }
  
   
  matrix expm2(matrix M){
    matrix[rows(M),rows(M)] out;
    int z0[rows(out)] = checkoffdiagzero(M);
    int z1[sum(z0)]; //contains which rowcols need full expm
    int count=1;
    for(j in 1:cols(M)){
      if(z0[j]){
        z1[count]=j;
        count+=1;
      } else {
        for(i in 1:rows(M)){
          if(i!=j){
          out[i,j] =  0; 
          out[j,i] = 0;
          } else out[i,i] = exp(M[i,j]);
        }
      }
    }
    if(size(z1)) out[z1,z1] = matrix_exp(M[z1,z1]);
    return out;
  }
  
}
data {
  int<lower=0> ndatapoints;
  int<lower=1> nmanifest;
  int<lower=1> nlatent;
  int nlatentpop;
  int nsubjects;
  int<lower=0> ntipred; // number of time independent covariates
  int<lower=0> ntdpred; // number of time dependent covariates
  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipredsdata;
  int nmissingtipreds;
  int ntipredeffects;
  real<lower=0> tipredsimputedscale;
  real<lower=0> tipredeffectscale;

  vector[nmanifest] Y[ndatapoints];
  int nopriors;
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

  int nobs_y[ndatapoints];  // number of observed variables per observation
  int whichobs_y[ndatapoints, nmanifest]; // index of which variables are observed per observation
  int ndiffusion; //number of latents involved in system noise calcs
  int derrind[ndiffusion]; //index of which latent variables are involved in system noise calculations

  int manifesttype[nmanifest];
  int nbinary_y[ndatapoints];  // number of observed binary variables per observation
  int whichbinary_y[ndatapoints, nmanifest]; // index of which variables are observed and binary per observation
  int ncont_y[ndatapoints];  // number of observed continuous variables per observation
  int whichcont_y[ndatapoints, nmanifest]; // index of which variables are observed and continuous per observation
  
  int intoverpop;
  int statedep[10];
  int choleskymats;
  int intoverstates;
  int verbose; //level of printing during model fit
  int subindices[10];
  int TIPREDEFFECTsetup[nparams, ntipred];
  int nrowmatsetup;
  int matsetup[nrowmatsetup,9];
  real matvalues[nrowmatsetup,6];
  int whenmat[54,5];
  int whenvecp[3,nparams];
  int whenvecs[6,nlatentpop];
  int matrixdims[54,2];
  int savescores;
  int savesubjectmatrices;
  int dokalman;
  int dokalmanrows[ndatapoints];
  real Jstep;
  real dokalmanpriormodifier;
  int intoverpopindvaryingindex[intoverpop ? nindvarying : 0];
  int nsJAxfinite;
  int sJAxfinite[nsJAxfinite];
  int nJyfinite;
  int sJyfinite[nJyfinite];
  int taylorheun;
  int difftype;
  int popcovn;
  int llsinglerow;
}
      
transformed data{
  matrix[nlatent+nindvarying,nlatent+nindvarying] IIlatentpop = diag_matrix(rep_vector(1,nlatent+nindvarying));
  vector[nlatentpop-nlatent] nlpzerovec = rep_vector(0,nlatentpop-nlatent);
  vector[nlatent+1] nlplusonezerovec = rep_vector(0,nlatent+1);
  int tieffectindices[nparams]=rep_array(0,nparams);
  int ntieffects = 0;
  
  if(ntipred >0){
    for(pi in 1:nparams){
      if(sum(TIPREDEFFECTsetup[pi,]) > .5){
      ntieffects+=1;
      tieffectindices[ntieffects] = pi;
      }
    }
  }
  
}
      
parameters{
  vector[nparams] rawpopmeans; // population level means 

  vector[nindvarying] rawpopsdbase; //population level std dev
  vector[nindvaryingoffdiagonals] sqrtpcov; // unconstrained basis of correlation parameters
  vector[intoverpop ? 0 : nindvarying] baseindparams[intoverpop ? 0 : nsubjects]; //vector of subject level deviations, on the raw scale
  
  vector[ntipredeffects] tipredeffectparams; // effects of time independent covariates
  vector[nmissingtipreds] tipredsimputed;
  
  vector[intoverstates ? 0 : nlatentpop*ndatapoints] etaupdbasestates; //sampled latent states posterior
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying, nindvarying] rawpopc[4];


  real ll = 0;
  vector[savescores ? ndatapoints : 1] llrow = rep_vector(0.0,savescores ? ndatapoints : 1);
  matrix[nlatentpop,nlatentpop] etacova[3,savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ycova[3,savescores ? ndatapoints : 0];
  vector[nlatentpop] etaa[3,savescores ? ndatapoints : 0];
  vector[nmanifest] ya[3,savescores ? ndatapoints : 0];
  
      matrix[matrixdims[10, 1], matrixdims[10, 2] ] pop_PARS;
      matrix[matrixdims[1, 1], matrixdims[1, 2] ] pop_T0MEANS;
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] pop_LAMBDA;
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] pop_DRIFT;
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] pop_DIFFUSION;
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] pop_MANIFESTVAR;
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] pop_MANIFESTMEANS;
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] pop_CINT;
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] pop_T0VAR;
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] pop_TDPREDEFFECT;
      matrix[matrixdims[31, 1], matrixdims[31, 2] ] pop_DIFFUSIONcov;
      matrix[matrixdims[21, 1], matrixdims[21, 2] ] pop_asymCINT;
      matrix[matrixdims[22, 1], matrixdims[22, 2] ] pop_asymDIFFUSION;
      matrix[matrixdims[10, 1], matrixdims[10, 2] ] PARS[ (savesubjectmatrices && sum(whenmat[10,1:5])) ? nsubjects : 0];
      matrix[matrixdims[1, 1], matrixdims[1, 2] ] T0MEANS[ (savesubjectmatrices && sum(whenmat[1,1:5])) ? nsubjects : 0];
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] LAMBDA[ (savesubjectmatrices && sum(whenmat[2,1:5])) ? nsubjects : 0];
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] DRIFT[ (savesubjectmatrices && sum(whenmat[3,1:5])) ? nsubjects : 0];
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] DIFFUSION[ (savesubjectmatrices && sum(whenmat[4,1:5])) ? nsubjects : 0];
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] MANIFESTVAR[ (savesubjectmatrices && sum(whenmat[5,1:5])) ? nsubjects : 0];
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] MANIFESTMEANS[ (savesubjectmatrices && sum(whenmat[6,1:5])) ? nsubjects : 0];
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] CINT[ (savesubjectmatrices && sum(whenmat[7,1:5])) ? nsubjects : 0];
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] T0VAR[ (savesubjectmatrices && sum(whenmat[8,1:5])) ? nsubjects : 0];
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] TDPREDEFFECT[ (savesubjectmatrices && sum(whenmat[9,1:5])) ? nsubjects : 0];
      matrix[matrixdims[31, 1], matrixdims[31, 2] ] DIFFUSIONcov[ (savesubjectmatrices && sum(whenmat[31,1:5])) ? nsubjects : 0];
      matrix[matrixdims[21, 1], matrixdims[21, 2] ] asymCINT[ (savesubjectmatrices && sum(whenmat[21,1:5])) ? nsubjects : 0];
      matrix[matrixdims[22, 1], matrixdims[22, 2] ] asymDIFFUSION[ (savesubjectmatrices && sum(whenmat[22,1:5])) ? nsubjects : 0];

  matrix[ntipred ? (nmissingtipreds ? nsubjects : 0) : 0, ntipred ? (nmissingtipreds ? ntipred : 0) : 0] tipreds; //tipred values to fill from data and, when needed, imputation vector
  matrix[nparams, ntipred] TIPREDEFFECT; //design matrix of individual time independent predictor effects

  if(ntipred > 0){ 
    if(nmissingtipreds > 0){
    int counter = 0;
    for(coli in 1:cols(tipreds)){ //insert missing ti predictors
      for(rowi in 1:rows(tipreds)){
        if(tipredsdata[rowi,coli]==99999) {
          counter += 1;
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
    int counter =0;
    rawpopsd = log1p_exp(2*rawpopsdbase-1) .* sdscale + 1e-10; // sqrts of proportions of total variance
    for(j in 1:nindvarying){
      rawpopc[1,j,j] = rawpopsd[j]; //used with intoverpop
      for(i in 1:nindvarying){
        if(i > j){
          counter += 1;
          rawpopc[1,i,j]=sqrtpcov[counter];
          rawpopc[1,j,i]=0;//sqrtpcov[counter];
        }
      }
    }
    rawpopc[3] = tcrossprod( constraincorsqrt(rawpopc[1]));
    rawpopc[4] = makesym(quad_form_diag(rawpopc[3], rawpopsd +1e-8),verbose,1);
    rawpopc[2] = cholesky_decompose(rawpopc[4]); 
  }//end indvarying par setup

  {

  int si = 0;
  int prevrow=0;
  real prevdt=0;
  real dt=1; //initialise to make sure drift is computed on first pass
  real dtsmall;
  int dtchange=1;
  int T0check=0;
  int counter = 1;
  matrix[nlatentpop, nlatentpop] etacov; //covariance of latent states

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] syprior;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypriorcov_sqrt = rep_matrix(0,nmanifest,nmanifest); 
  matrix[nmanifest, nmanifest] ycov; 
  
  matrix[nlatentpop,nlatentpop] Je[savescores ? ndatapoints : 1]; //time evolved jacobian, saved for smoother

  vector[nlatentpop] state = rep_vector(-999,nlatentpop); 
  vector[nlatentpop] statetf;
  matrix[nlatentpop,nlatentpop] sJAx; //Jacobian for drift
  matrix[nlatentpop,nlatentpop] sJ0; //Jacobian for t0
  matrix[nlatentpop,nlatentpop] sJtd;//diag_matrix(rep_vector(1),nlatentpop); //Jacobian for nltdpredeffect
  matrix[ nmanifest,nlatentpop] sJy;//Jacobian for measurement 
  matrix[ nmanifest,nlatentpop] Jysaved[savescores ? ndatapoints : 0];
  
  

  //linear continuous time calcs
  matrix[nlatent+1,nlatent+1] discreteDRIFT;
  matrix[nlatent,nlatent] discreteDIFFUSION = rep_matrix(0.0,nlatent,nlatent);

  
  vector[nparams] rawindparams = rawpopmeans;
  vector[nparams] indparams;

  //dynamic system matrices
  
      matrix[matrixdims[10, 1], matrixdims[10, 2] ] sPARS;
      matrix[matrixdims[1, 1], matrixdims[1, 2] ] sT0MEANS;
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] sLAMBDA;
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] sDRIFT;
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] sDIFFUSION;
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] sMANIFESTVAR;
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] sMANIFESTMEANS;
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] sCINT;
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] sT0VAR;
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] sTDPREDEFFECT;
      matrix[matrixdims[31, 1], matrixdims[31, 2] ] sDIFFUSIONcov;
      matrix[matrixdims[21, 1], matrixdims[21, 2] ] sasymCINT;
      matrix[matrixdims[22, 1], matrixdims[22, 2] ] sasymDIFFUSION;
  
  sasymDIFFUSION = rep_matrix(0,nlatent,nlatent); //in case of derrindices need to init
  sDIFFUSIONcov = rep_matrix(0,nlatent,nlatent);

  for(rowi in 1:(dokalman ? ndatapoints :1)){
    int nsubs = prevrow ? 1 : 2; //first pass through needs 2 subjects for sub 0 (pop pars)
    int rowsubs[nsubs]; //container for 1 or 2 subject refs
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
  
    si = subject[rowi];
    rowsubs[size(rowsubs)] = si; //setup 2 subject index container
    if(prevrow==0) rowsubs[1] = 0;
    
    if(prevrow != 0) T0check = (si==subject[prevrow]) ? (T0check+1) : 0; //if same subject, add one, else zero
    if(T0check > 0){
      dt = time[rowi] - time[prevrow];
      dtchange = dt!=prevdt; 
      prevdt = dt; //update previous dt store after checking for change
    }
    //if(savescores || prevrow==0) Je[savescores ? rowi : 1] = IIlatentpop[1:nlatentpop,1:nlatentpop]; //elements updated later
    if(savescores && prevrow!=0) Je[rowi,,] = Je[prevrow,,];
    
    
  for(subi in rowsubs){
    
    if(T0check == 0) { // calculate initial matrices if this is first row for si
  
  rawindparams=rawpopmeans;
  
  if(subi > 0 && nindvarying > 0 && intoverpop==0)  rawindparams[indvaryingindex] += rawpopc[2] * baseindparams[subi];

  if(subi > 0 &&  ntieffects > 0){
  if(nmissingtipreds > 0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipreds[subi]';
    
    if(nmissingtipreds==0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipredsdata[subi]';
  }

  indparams[whichequals(whenvecp[subi ? 2 : 1], 0, 0)]= 
    parvectform(size(whichequals(whenvecp[subi ? 2 : 1], 0, 0)),rawindparams, 
    0, matsetup, matvalues, subi, subindices, whenvecp[subi ? 2 : 1]);
     
  if(whenmat[1, 5] >= (subi ? 1 : 0)) sT0MEANS = 
    mcalc(sT0MEANS, indparams, statetf, {0}, 1, matsetup, matvalues, subi, subindices); // base t0means to init
      
  for(li in 1:nlatentpop) if(!is_nan(sT0MEANS[li,1])) state[li] = sT0MEANS[li,1]; //in case of t0 dependencies, may have missingness
  
  statetf[whichequals(whenvecs[1],0,0)] = parvectform(size(whichequals(whenvecs[1],0,0)),state, 1,
    matsetup, matvalues, si, subindices, whenvecs[1]);   
    
  
    
  sT0MEANS=mcalc(sT0MEANS,indparams, statetf,{0,1}, 1, matsetup, matvalues, si, subindices); 
sT0VAR=mcalc(sT0VAR,indparams, statetf,{0,1}, 8, matsetup, matvalues, si, subindices); 
sJ0=mcalc(sJ0,indparams, statetf,{0,1}, 51, matsetup, matvalues, si, subindices); 

      
  
    for(li in 1:nlatentpop) if(is_nan(state[li])) state[li] = sT0MEANS[li,1]; //finish updating state
    
    //init other system matrices
    if(whenmat[10,5] || subi==0) sPARS=mcalc(sPARS,indparams, statetf,{0}, 10, matsetup, matvalues, subi, subindices); 
  if(whenmat[1,5] || subi==0) sT0MEANS=mcalc(sT0MEANS,indparams, statetf,{0}, 1, matsetup, matvalues, subi, subindices); 
  if(whenmat[2,5] || subi==0) sLAMBDA=mcalc(sLAMBDA,indparams, statetf,{0}, 2, matsetup, matvalues, subi, subindices); 
  if(whenmat[3,5] || subi==0) sDRIFT=mcalc(sDRIFT,indparams, statetf,{0}, 3, matsetup, matvalues, subi, subindices); 
  if(whenmat[4,5] || subi==0) sDIFFUSION=mcalc(sDIFFUSION,indparams, statetf,{0}, 4, matsetup, matvalues, subi, subindices); 
  if(whenmat[5,5] || subi==0) sMANIFESTVAR=mcalc(sMANIFESTVAR,indparams, statetf,{0}, 5, matsetup, matvalues, subi, subindices); 
  if(whenmat[6,5] || subi==0) sMANIFESTMEANS=mcalc(sMANIFESTMEANS,indparams, statetf,{0}, 6, matsetup, matvalues, subi, subindices); 
  if(whenmat[7,5] || subi==0) sCINT=mcalc(sCINT,indparams, statetf,{0}, 7, matsetup, matvalues, subi, subindices); 
  if(whenmat[8,5] || subi==0) sT0VAR=mcalc(sT0VAR,indparams, statetf,{0}, 8, matsetup, matvalues, subi, subindices); 
  if(whenmat[9,5] || subi==0) sTDPREDEFFECT=mcalc(sTDPREDEFFECT,indparams, statetf,{0}, 9, matsetup, matvalues, subi, subindices); 
  if(whenmat[51,5] || subi==0) sJ0=mcalc(sJ0,indparams, statetf,{0}, 51, matsetup, matvalues, subi, subindices); 
  if(whenmat[52,5] || subi==0) sJAx=mcalc(sJAx,indparams, statetf,{0}, 52, matsetup, matvalues, subi, subindices); 
  if(whenmat[53,5] || subi==0) sJtd=mcalc(sJtd,indparams, statetf,{0}, 53, matsetup, matvalues, subi, subindices); 
  if(whenmat[54,5] || subi==0) sJy=mcalc(sJy,indparams, statetf,{0}, 54, matsetup, matvalues, subi, subindices); 

    
    
  if(subi <= (subindices[8] ? nsubjects : 0)) {
   if(intoverpop && nindvarying > 0) sT0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopc[1];
    sT0VAR = makesym(sdcovsqrt2cov(sT0VAR,choleskymats),verbose,1); 

    if(intoverpop && nindvarying > 0){ //adjust cov matrix for transforms
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
    
    if(verbose>1) print("nl t0var = ", sT0VAR, "   sJ0 = ", sJ0);
    etacov = quad_form(sT0VAR, sJ0'); //probably unneeded, inefficient

    } //end T0 matrices
    
if(verbose > 1) print ("below t0 row ", rowi);

      if(subi==0 || T0check>0){
        vector[nlatent] base;
        real intstepi = 0;
        dtsmall = dt / ceil(dt / maxtimestep);
        
        while(intstepi < (dt-1e-10)){
          intstepi = intstepi + dtsmall;
          
    {
    int zeroint[1];
    vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(sJAxfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0)  state[statei] += Jstep;
      statetf[whichequals(whenvecs[2],0,0)] = 
         parvectform(size(whichequals(whenvecs[2],0,0)),state, 2, matsetup, matvalues, si, subindices, whenvecs[2]);
    
    if(sum(whenmat[10,{2}]) > 0)sPARS=mcalc(sPARS,indparams, statetf,{2}, 10, matsetup, matvalues, subi, subindices); 
if(sum(whenmat[3,{2}]) > 0)sDRIFT=mcalc(sDRIFT,indparams, statetf,{2}, 3, matsetup, matvalues, subi, subindices); 
if(sum(whenmat[7,{2}]) > 0)sCINT=mcalc(sCINT,indparams, statetf,{2}, 7, matsetup, matvalues, subi, subindices); 

       
      
      if(statei > 0) {
        sJAx[1:nlatent,statei] =  sDRIFT * state[1:nlatent] + sCINT[,1]; //compute new change
         if(verbose>1) print("sJAx ",sJAx);
      }
      if(statei== 0 && size(sJAxfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
        base = sDRIFT * state[1:nlatent] + sCINT[,1];
        if(verbose>1) print("base = ",base,"    sjaxinit= ",sJAx);
        for(fi in sJAxfinite){
          sJAx[1:nlatent,fi] = (sJAx[1:nlatent,fi] - base) / Jstep; //new - baseline change divided by stepsize
        }
      }
    }
    if(verbose>1) print("sJAx ",sJAx);
    }
    statetf[whichequals(whenvecs[2],0,0)] = 
         parvectform(size(whichequals(whenvecs[2],0,0)),state, 2, matsetup, matvalues, si, subindices, whenvecs[2]);
    
    if(sum(whenmat[4,{2}]) > 0)sDIFFUSION=mcalc(sDIFFUSION,indparams, statetf,{2}, 4, matsetup, matvalues, subi, subindices); 
if(sum(whenmat[52,{2}]) > 0)sJAx=mcalc(sJAx,indparams, statetf,{2}, 52, matsetup, matvalues, subi, subindices); 
 
      if(subi==0 || whenmat[4,2] || (T0check==1 && whenmat[4,5])) sDIFFUSIONcov[derrind,derrind] = sdcovsqrt2cov(sDIFFUSION[derrind,derrind],choleskymats);
      
        if(continuoustime){
          if(taylorheun==0){
            if(subi==0 || dtchange==1 || statedep[3]||statedep[4]||statedep[7] || 
              (T0check == 1 && (subindices[3] + subindices[4] + subindices[7]) > 0)){
                
              
              
              if(difftype==0 || (statedep[3]==0 && statedep[4]==0)){
                discreteDRIFT = expm2(append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec') * dtsmall);
                Je[savescores ? rowi : 1,,] =  expm2(sJAx * dtsmall);
               
                if(subi==0 || statedep[3]||statedep[4]|| (T0check==1 && (whenmat[4,5] || whenmat[3,5]))){
                sasymDIFFUSION[derrind,derrind] = to_matrix(  -sqkron_sumii(sJAx[derrind,derrind]) \ 
                  to_vector(sDIFFUSIONcov[derrind,derrind]), ndiffusion,ndiffusion);
                }
                discreteDIFFUSION[derrind,derrind] =  sasymDIFFUSION[derrind,derrind] - 
                  quad_form( sasymDIFFUSION[derrind,derrind], Je[savescores ? rowi : 1, derrind,derrind]' );
              }
            } 
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent]; // ???compute before new diffusion calcs
            if(intoverstates==1 || savescores==1){
              etacov = quad_form(etacov, Je[savescores ? rowi : 1,,]');
              etacov[derrind,derrind] += discreteDIFFUSION[derrind,derrind]; 
            }
          }
            
           
              
            if(intstepi >= (dt-1e-10) && savescores) Je[rowi,,] = expm2(sJAx * dt); //save approximate exponentiated jacobian for smoothing
          }
  
          if(continuoustime==0){ 
            Je[savescores ? rowi : 1,,] = sJAx;
            if(intoverstates==1 || savescores==1){
              etacov = quad_form(etacov, sJAx');
              etacov[ derrind, derrind ] += tcrossprod(sDIFFUSION[ derrind, derrind ]); 
            }
            discreteDIFFUSION=sDIFFUSIONcov;
            discreteDRIFT=append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec');
            discreteDRIFT[nlatent+1,nlatent+1] = 1;
            state[1:nlatent] = (discreteDRIFT * append_row(state[1:nlatent],1.0))[1:nlatent];
            
          }
        }
      } // end non linear time update
    
    if(ntdpred > 0) {
      int nonzerotdpred = 0;
      for(tdi in 1:ntdpred) if(tdpreds[rowi,tdi] != 0.0) nonzerotdpred = 1;
      if(nonzerotdpred){
      statetf[whichequals(whenvecs[3],0,0)] = 
         parvectform(size(whichequals(whenvecs[3],0,0)),state, 3, matsetup, matvalues, si, subindices, whenvecs[3]);
    
    if(sum(whenmat[9,{3}]) > 0)sTDPREDEFFECT=mcalc(sTDPREDEFFECT,indparams, statetf,{3}, 9, matsetup, matvalues, subi, subindices); 
if(sum(whenmat[53,{3}]) > 0)sJtd=mcalc(sJtd,indparams, statetf,{3}, 53, matsetup, matvalues, subi, subindices); 

      
        state[1:nlatent] +=   (sTDPREDEFFECT * tdpreds[rowi]); //tdpred effect only influences at observed time point
        if(statedep[9]) etacov = quad_form(etacov,sJtd'); //could be optimized
      }
    }//end nonlinear tdpred

  if(intoverstates==0){
    if(T0check==0) state += cholesky_decompose(sT0VAR) * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)];
    if(T0check>0) state[derrind] +=  cholesky_decompose(makesym(discreteDIFFUSION[derrind,derrind],verbose,1)) * 
      (etaupdbasestates[(1+(rowi-1)*nlatentpop):(nlatent+(rowi-1)*nlatentpop)])[derrind];
  }

if(verbose > 1){
  print("etaprior = ", state);
  print("etapriorcov = ", etacov);
}

if(savescores){
  etacova[1,rowi] = etacov; 
  etaa[1,rowi] = state;
}

 if ((subi==0 || nobs_y[rowi] > 0 || savescores) && dokalmanrows[rowi] ==1){ //do this section for 0th subject as well to init matrices
    
      
    int zeroint[1];
    vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(sJyfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0 && (savescores + intoverstates) > 0)  state[statei] += Jstep;
      statetf[whichequals(whenvecs[4],0,0)] = 
         parvectform(size(whichequals(whenvecs[4],0,0)),state, 4, matsetup, matvalues, si, subindices, whenvecs[4]);
    
    if(sum(whenmat[10,{4}]) > 0)sPARS=mcalc(sPARS,indparams, statetf,{4}, 10, matsetup, matvalues, subi, subindices); 
if(sum(whenmat[2,{4}]) > 0)sLAMBDA=mcalc(sLAMBDA,indparams, statetf,{4}, 2, matsetup, matvalues, subi, subindices); 
if(sum(whenmat[5,{4}]) > 0)sMANIFESTVAR=mcalc(sMANIFESTVAR,indparams, statetf,{4}, 5, matsetup, matvalues, subi, subindices); 
if(sum(whenmat[6,{4}]) > 0)sMANIFESTMEANS=mcalc(sMANIFESTMEANS,indparams, statetf,{4}, 6, matsetup, matvalues, subi, subindices); 
if(sum(whenmat[54,{4}]) > 0)sJy=mcalc(sJy,indparams, statetf,{4}, 54, matsetup, matvalues, subi, subindices); 
 
      if(statei > 0 && (savescores + intoverstates) > 0) {
        sJy[o,statei] =  sLAMBDA[o] * state[1:nlatent] + sMANIFESTMEANS[o,1]; //compute new change
        sJy[o1,statei] = to_vector(inv_logit(to_array_1d(sJy[o1,statei])));
         if(verbose>1) print("sJy ",sJy);
      }
      if(statei==0){
        syprior[o] = sLAMBDA[o] * state[1:nlatent] + sMANIFESTMEANS[o,1];
        syprior[o1] = to_vector(inv_logit(to_array_1d( syprior[o1] )));
        if(size(sJyfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
          if(verbose>1) print("syprior = ",syprior,"    sJyinit= ",sJy);
          for(fi in sJyfinite){
            sJy[o,fi] = (sJy[o,fi] - syprior[o]) / Jstep; //new - baseline change divided by stepsize
          }
        }
      }
    }
    if(verbose>1) print("sJy ",sJy);
    
      
 } // end measurement init
 
   if(subi > 0 && dokalmanrows[rowi] ==1 && (nobs_y[rowi] > 0 || savescores)){   // if some observations create right size matrices for missingness and calculate...

      if(intoverstates==1 || savescores==1) { //classic kalman
        ycov[o,o] = quad_form(etacov, sJy[o,]'); // + sMANIFESTVAR[o,o]; shifted measurement error down
        for(wi in 1:nmanifest){ 
          if(Y[rowi,wi] != 99999 || savescores==1) ycov[wi,wi] += square(sMANIFESTVAR[wi,wi]);
          if(manifesttype[wi]==1 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += fabs((syprior[wi] - 1) .* (syprior[wi]));
          if(manifesttype[wi]==2 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += square(fabs((syprior[wi] - round(syprior[wi])))); 
        }
      }
        
      if(intoverstates==0) { //sampled states
        if(ncont_y[rowi] > 0) for(mi in o0) ypriorcov_sqrt[mi,mi] = sqrt(sMANIFESTVAR[mi,mi]);
      }
        
     
err[od] = Y[rowi,od] - syprior[od]; // prediction error
    
      if(intoverstates==1 && size(od) > 0) {
       // matrix[size(od),size(od)] ycovi;
        
         if(verbose > 1) {
          print("init rowi =",rowi, "  si =", si, "  state =",state,"  etacov ",etacov,
          " indparams = ", indparams,
            "  syprior =",syprior,"  ycov ",ycov, "  K ",K,
            "  sDRIFT =", sDRIFT, " sDIFFUSION =", sDIFFUSION, " sCINT =", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
            "  sT0VAR", sT0VAR,  " sT0MEANS ", sT0MEANS, "sLAMBDA = ", sLAMBDA, "  sJy = ",sJy,
            " discreteDRIFT = ", discreteDRIFT, "  discreteDIFFUSION ", discreteDIFFUSION, "  sasymDIFFUSION ", sasymDIFFUSION, 
            " DIFFUSIONcov = ", sDIFFUSIONcov,
            "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
         }
        //ycovi= inverse(ypriorcov_sqrt[o0d,o0d]); //inverse avoids weird singularity when dividing by same matrix, ie no merror
        K[,od] = etacov * sJy[od,]' / makesym(ycov[od,od],verbose,1);// * multiply_lower_tri_self_transpose(ycovi');// ycov[od,od]; 
        etacov += -K[,od] * sJy[od,] * etacov;
        state +=  (K[,od] * err[od]);
      }
      
      if(savescores==1) {
        ya[1,rowi] = syprior[o];
        etaa[2,rowi] = state;
        ycova[1,rowi] = ycov;
        etacova[2,rowi] = etacov;
        ycova[2,rowi] = quad_form(etacov, sJy');
        for(wi in 1:nmanifest) ycova[2,rowi,wi,wi] += square(sMANIFESTVAR[wi,wi]);
        ya[2,rowi] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * state[1:nlatent];
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
  
      if(nbinary_y[rowi] > 0) llrow[savescores ? rowi : 1] += sum(log(Y[rowi,o1d] .* (syprior[o1d]) + (1-Y[rowi,o1d]) .* (1-syprior[o1d]))); 
  
        if(size(o0d) > 0 && (llsinglerow==0 || llsinglerow == rowi)){
          ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(makesym(ycov[o0d,o0d],verbose,1));
           llrow[savescores ? rowi : 1] +=  multi_normal_cholesky_lpdf(Y[rowi,o0d] | syprior[o0d], ypriorcov_sqrt[o0d,o0d]);
           //errtrans[counter:(counter + ncont_y[rowi]-1)] = 
             //mdivide_left_tri_low(ypriorcov_sqrt[o0d,o0d], err[o0d]); //transform pred errors to standard normal dist and collect
           //ll+= -sum(log(diagonal(ypriorcov_sqrt[o0d,o0d]))); //account for transformation of scale in loglik
           //counter += ncont_y[rowi];
        }
      
      if(savescores) Jysaved[rowi] = sJy;
    }//end nobs > 0 section
    
       // store system matrices
       
  if(subi <= ((subindices[3] + subindices[7])  ? nsubjects : 0)){
    if(continuoustime==1) sasymCINT[,1] =  -sDRIFT[1:nlatent,1:nlatent] \ sCINT[ ,1 ];
    if(continuoustime==0) sasymCINT[,1] =  add_diag(-sDRIFT[1:nlatent,1:nlatent],1) \ sCINT[,1 ];
  }
  
  if(subi <= ((subindices[3] + subindices[7])  ? nsubjects : 0) && !continuoustime) sasymDIFFUSION[ derrind, derrind ] = 
    to_matrix( (add_diag( 
      -sqkron_prod(sDRIFT[ derrind, derrind ], sDRIFT[ derrind, derrind ]),1)) \  
      to_vector(sDIFFUSIONcov[ derrind, derrind ]), ndiffusion, ndiffusion);


    
  if(savesubjectmatrices && subi > 0 && (rowi==ndatapoints || subject[rowi+1] != subject[rowi])){
    if(sum(whenmat[10,1:5]) > 0) PARS[subi] = sPARS;
    if(sum(whenmat[1,1:5]) > 0) T0MEANS[subi] = sT0MEANS;
    if(sum(whenmat[2,1:5]) > 0) LAMBDA[subi] = sLAMBDA;
    if(sum(whenmat[3,1:5]) > 0) DRIFT[subi] = sDRIFT;
    if(sum(whenmat[4,1:5]) > 0) DIFFUSION[subi] = sDIFFUSION;
    if(sum(whenmat[5,1:5]) > 0) MANIFESTVAR[subi] = sMANIFESTVAR;
    if(sum(whenmat[6,1:5]) > 0) MANIFESTMEANS[subi] = sMANIFESTMEANS;
    if(sum(whenmat[7,1:5]) > 0) CINT[subi] = sCINT;
    if(sum(whenmat[8,1:5]) > 0) T0VAR[subi] = sT0VAR;
    if(sum(whenmat[9,1:5]) > 0) TDPREDEFFECT[subi] = sTDPREDEFFECT;
    if(sum(whenmat[31,1:5]) > 0) DIFFUSIONcov[subi] = sDIFFUSIONcov;
    if(sum(whenmat[21,1:5]) > 0) asymCINT[subi] = sasymCINT;
    if(sum(whenmat[22,1:5]) > 0) asymDIFFUSION[subi] = sasymDIFFUSION;
  }
  if(subi == 0){
pop_PARS = sPARS; pop_T0MEANS = sT0MEANS; pop_LAMBDA = sLAMBDA; pop_DRIFT = sDRIFT; pop_DIFFUSION = sDIFFUSION; pop_MANIFESTVAR = sMANIFESTVAR; pop_MANIFESTMEANS = sMANIFESTMEANS; pop_CINT = sCINT; pop_T0VAR = sT0VAR; pop_TDPREDEFFECT = sTDPREDEFFECT; pop_DIFFUSIONcov = sDIFFUSIONcov; pop_asymCINT = sasymCINT; pop_asymDIFFUSION = sasymDIFFUSION; 
  }
    
  } // end extra subi == 0 for subject matrix creation

  
  if(savescores && (rowi==ndatapoints || subject[rowi+1] != subject[rowi])){ //at subjects last datapoint, smooth
    int sri = rowi;
    int subi = si; //copy needed
    while(sri>0 && subject[sri]==si){
      if(sri==rowi) {
        etaa[3,sri]=etaa[2,sri];
        etacova[3,sri]=etacova[2,sri];
      } else{ // could be improved by recomputing jacobian
        matrix[nlatentpop,nlatentpop] smoother;
        smoother = etacova[2,sri] * Je[sri+1,,]' / makesym(etacova[1,sri+1],verbose,1);
        etaa[3,sri]= etaa[2,sri] + smoother * (etaa[3,sri+1] - etaa[1,sri+1]);
        etacova[3,sri]= etacova[2,sri] + smoother * ( etacova[3,sri+1] - etacova[1,sri+1]) * smoother';

      }
      state=etaa[3,sri];
      ya[3,sri] = syprior;
      ycova[3,sri] = quad_form(etacova[3,sri], Jysaved[rowi]'); //could be improved by recomputing jacobian
      for(wi in 1:nmanifest){
        ycova[3,sri,wi,wi] += square(sMANIFESTVAR[wi,wi]);
        if(manifesttype[wi]==1) ycova[3,sri,wi,wi] += fabs((ya[3,sri,wi] - 1) * (ya[3,sri,wi]));
        if(manifesttype[wi]==2) ycova[3,sri,wi,wi] += square(fabs((ya[3,sri,wi] - round(ya[3,sri,wi])))); 
      }
      sri += -1;
      while(sri > 0 && dokalmanrows[sri]==0) sri+= -1; //skip rows if requested
    }
  } //end smoother

  
  
  prevrow = rowi; //update previous row marker only after doing necessary calcs
}//end rowi
ll+=sum(llrow);

  }
}
      
model{
  if(intoverpop==0 && nindvarying > 0) target+= multi_normal_cholesky_lpdf(baseindparams | rep_vector(0,nindvarying), IIlatentpop[1:nindvarying,1:nindvarying]);

  if(ntipred > 0){ 
    if(nopriors==0) target+= dokalmanpriormodifier * normal_lpdf(tipredeffectparams| 0, tipredeffectscale);
    target+= normal_lpdf(tipredsimputed| 0, tipredsimputedscale); //consider better handling of this when using subset approach
  }

  if(nopriors==0){ //if split files over subjects, just compute priors once
   target+= dokalmanpriormodifier * normal_lpdf(rawpopmeans|0,1);
  
    if(nindvarying > 0){
      if(nindvarying >1) target+= dokalmanpriormodifier * normal_lpdf(sqrtpcov | 0, 1);
      target+= dokalmanpriormodifier * normal_lpdf(rawpopsdbase | 0,1);
    }
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
        x[ri,] = (rawpopc[2] * 
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
          if(matsetup[ri,9] <=0 && matsetup[ri,3]==pi && matsetup[ri,8]<=0) { //if a free parameter 
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
