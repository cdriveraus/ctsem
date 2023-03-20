
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
    int d=rows(mat);
    matrix[d,d] o;
    vector[d] ss = rep_vector(0,d);
    vector[d] s = rep_vector(0,d);
    real r;
    real r3;
    real r4;
    real r1;
    
    for(i in 1:d){
      for(j in 1:d){
        if(j > i) {
          ss[i] +=square(mat[j,i]);
          s[i] +=mat[j,i];
        }
        if(j < i){
          ss[i] += square(mat[i,j]);
          s[i] += mat[i,j];
        }
      }
      s[i] += 1e-5;
      ss[i] += 1e-5;
    }

    
    for(i in 1:d){
      o[i,i]=0;
      r1=sqrt(ss[i]);
      r3=(fabs(s[i]))/(r1)-1;
      r4=sqrt(log1p_exp(2*(fabs(s[i])-s[i]-1)-4));
      r=(r4*((r3))+1)*r4+1;
      r=(sqrt(ss[i]+r));
      for(j in 1:d){
        if(j > i)  o[i,j]=mat[j,i]/r;
        if(j < i) o[i,j] = mat[i,j] /r;
      }
      o[i,i]=sqrt(1-sum(square(o[i,]))+1e-5);
    }

    return o;
}  

  matrix sdcovsqrt2cov(matrix mat, int choleskymats){ //covariance from cholesky or unconstrained cor sq root
    if(choleskymats< 1) {
      //if(choleskymats== -1){
        return(tcrossprod(diag_pre_multiply(diagonal(mat),constraincorsqrt(mat))));
      //} else {
      //  return(quad_form_diag(constraincorsqrt(mat),diagonal(mat)));
      //}
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

 
 matrix ksolve(matrix A, matrix Q, int verbose){
  int d= rows(A);
  int d2= (d*d-d)/2;
  matrix[d+d2,d+d2] O;
  vector[d+d2] triQ;
  matrix[d,d] AQ;
  int z=0; //z is row of output
  for(j in 1:d){//for column reference of solution vector
    for(i in 1:j){ //and row reference...
      if(j >= i){ //if i and j denote a covariance parameter (from upper tri)
        int y=0; //start new output row
        z+=1; //shift current output row down
        
        for(ci in 1:d){//for columns and
          for(ri in 1:d){ //rows of solution
            if(ci >= ri){ //when in upper tri (inc diag)
              y+=1; //move to next column of output
              
              if(i==j){ //if output row is for a diagonal element
                if(ri==i) O[z,y] = 2*A[ri,ci];
                if(ci==i) O[z,y] = 2*A[ci,ri];
              }
              
              if(i!=j){ //if output row is not for a diagonal element
                if(y==z) O[z,y] = A[ri,ri] + A[ci,ci]; //if column of output matches row of output, sum both A diags
                if(y!=z){ //otherwise...
                  // if solution element we refer to is related to output row...
                  if(ci==ri){ //if solution element is a variance
                    if(ci==i) O[z,y] = A[j,ci]; //if variance of solution corresponds to row of our output
                    if(ci==j) O[z,y] = A[i,ci]; //if variance of solution corresponds to col of our output
                  }
                  if(ci!=ri && (ri==i||ri==j||ci==i||ci==j)){//if solution element is a related covariance
                    //for row 1,2 / 2,1 of output, if solution row ri 1 (match) and column ci 3, we need A[2,3]
                    if(ri==i) O[z,y] = A[j,ci];
                    if(ri==j) O[z,y] = A[i,ci];
                    if(ci==i) O[z,y] = A[j,ri];
                    if(ci==j) O[z,y] = A[i,ri];
                  }
                }
              }
              if(is_nan(O[z,y])) O[z,y]=0;
            }
          }
        }
      }
    }
  }
  
  z=0; //get upper tri of Q
  for(j in 1:d){
    for(i in 1:j){
    z+=1;
    triQ[z] = Q[i,j];
    }
  }
  triQ=-O \ triQ; //get upper tri of asymQ
  
    z=0; // put upper tri of asymQ into matrix
  for(j in 1:d){
    for(i in 1:j){
    z+=1;
    AQ[i,j] = triQ[z];
    if(i!=j) AQ[j,i] = triQ[z];
    }
  }
  
  if(verbose>1) print("AQ = ", AQ, "   triQ = ", triQ, "   O = ", O);
  
  return AQ;
}

  matrix makesym(matrix mat, int verbose, int pd){
    matrix[rows(mat),cols(mat)] out;
    for(coli in 1:cols(mat)){
    //  if(pd ==1 && mat[coli,coli] < 1e-5){
     //   out[coli,coli] = 1e-5;// 
     // } else 
      out[coli,coli] = mat[coli,coli] + 1e-10; 
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
  
  // improve PARS when = 100 thing here too
   row_vector parvectform(int[] which, row_vector rawpar, int when, int[,] ms, data real[,] mval, int subi){
    row_vector[size(which)] parout;
    if(size(which)){
      for(whichout in 1:size(which)){
        int done=0; //only want to tform once, may be copies
      for(ri in 1:size(ms)){ //for each row of matrix setup
        if(!done){
        if((ms[ri,8]==when || ms[ri,8]==100)  && ms[ri,3] == which[whichout]){ //if correct when and free parameter //,not a copyrow,&& ms[ri,9] < 1
          if(subi ==0 ||  //if population parameter
            (ms[ri,3] > 0 && (ms[ri,5] > 0 || ms[ri,6] > 0 || ms[ri,8] > 0)) //or there is individual variation
            ){ //otherwise repeated values
            
            parout[whichout] = tform(rawpar[ms[ri,3] ], // was: whichequals(outwhen, ms[ri,3],1)[1], which outwhen refers to correct par
              ms[ri,4], mval[ri,2], mval[ri,3], mval[ri,4], mval[ri,6] ); 
            }
           done=1;
        }
      }
      }
      }
    }
  return parout;
  }
  
  
  matrix mcalc(matrix matin, vector tfpars, row_vector tfstates, int[] when, int m, int[,] ms, data real[,] mval, int subi){
    matrix[rows(matin),cols(matin)] matout;

    for(ri in 1:size(ms)){ //for each row of matrix setup
      int whenyes = 0;
      for(wi in 1:size(when)) if(when[wi]==ms[ri,8] || ms[ri,8]==100) whenyes = 1; //improve PARS when = 100 thing
      if(m==ms[ri,7] && whenyes){ // if correct matrix and when

        if(subi ==0 ||  //if population parameter
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
  
//  int nt0varstationary;
//  int nt0meansstationary;
//  int t0varstationary [nt0varstationary, 2];
//  int t0meansstationary [nt0meansstationary, 2];

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
  int statedep[54];
  int choleskymats;
  int intoverstates;
  int verbose; //level of printing during model fit
  int TIPREDEFFECTsetup[nparams, ntipred];
  int nrowmatsetup;
  int matsetup[nrowmatsetup,9];
  real matvalues[nrowmatsetup,6];
  int whenmat[54,5];
  int whenvecp[2,nparams];
  int whenvecs[6,nlatentpop];
  int matrixdims[54,2];
  int savescores;
  int savesubjectmatrices;
  int dokalman;
  int dokalmanrows[ndatapoints];
  int nsubsets;
  real Jstep;
  real priormod;
  int intoverpopindvaryingindex[intoverpop ? nindvarying : 0];
  int nJAxfinite;
  int JAxfinite[nJAxfinite];
  int nJyfinite;
  int Jyfinite[nJyfinite];
  int taylorheun;
  int popcovn;
  int llsinglerow;
  int laplaceprior[nparams];
  int laplaceprioronly;
  int laplacetipreds;
  int CINTnonzerosize;
  int CINTnonzero[CINTnonzerosize];
  int JAxDRIFTequiv;
}
      
transformed data{
  matrix[nlatent+nindvarying,nlatent+nindvarying] IIlatentpop = diag_matrix(rep_vector(1,nlatent+nindvarying));
  vector[nlatentpop-nlatent] nlpzerovec = rep_vector(0,nlatentpop-nlatent);
  vector[nlatent+1] nlplusonezerovec = rep_vector(0,nlatent+1);
  int tieffectindices[nparams]=rep_array(0,nparams);
  int ntieffects = 0;
  int dosmoother = savescores || savesubjectmatrices;
  
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
  vector[(nsubsets > 1) ? 1 : 0] subsetpar;
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying, nindvarying] rawpopcovbase;
  matrix[nindvarying, nindvarying] rawpopcov;
  matrix[nindvarying, nindvarying] rawpopcovchol;
  matrix[nindvarying, nindvarying] rawpopcorr;
  real subset = (nsubsets > 1) ? subsetpar[1] : 1.0;
  real firstsub = round(nsubjects*1.0/nsubsets*(subset-1)+1);
  real lastsub = round(nsubjects*1.0/nsubsets*(subset));
  real ll = 0;

  vector[dokalman ? ndatapoints : 1] llrow = rep_vector(0,dokalman ? ndatapoints : 1);
  matrix[nlatentpop,nlatentpop] etacova[3,savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ycova[3,savescores ? ndatapoints : 0];
  vector[nlatentpop] etaa[3,savescores ? ndatapoints : 0];
  vector[nmanifest] ya[3,savescores ? ndatapoints : 0];
  vector[0] calcs[0];
  
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
      matrix[matrixdims[32, 1], matrixdims[32, 2] ] pop_MANIFESTcov;
      matrix[matrixdims[33, 1], matrixdims[33, 2] ] pop_T0cov;
      matrix[matrixdims[21, 1], matrixdims[21, 2] ] pop_asymCINT;
      matrix[matrixdims[22, 1], matrixdims[22, 2] ] pop_asymDIFFUSIONcov;
      matrix[matrixdims[10, 1], matrixdims[10, 2] ] subj_PARS[ (savesubjectmatrices && (sum(whenmat[10,1:5]) || statedep[10])) ? nsubjects : 0];
      matrix[matrixdims[1, 1], matrixdims[1, 2] ] subj_T0MEANS[ (savesubjectmatrices && (sum(whenmat[1,1:5]) || statedep[1])) ? nsubjects : 0];
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] subj_LAMBDA[ (savesubjectmatrices && (sum(whenmat[2,1:5]) || statedep[2])) ? nsubjects : 0];
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] subj_DRIFT[ (savesubjectmatrices && (sum(whenmat[3,1:5]) || statedep[3])) ? nsubjects : 0];
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] subj_DIFFUSION[ (savesubjectmatrices && (sum(whenmat[4,1:5]) || statedep[4])) ? nsubjects : 0];
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] subj_MANIFESTVAR[ (savesubjectmatrices && (sum(whenmat[5,1:5]) || statedep[5])) ? nsubjects : 0];
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] subj_MANIFESTMEANS[ (savesubjectmatrices && (sum(whenmat[6,1:5]) || statedep[6])) ? nsubjects : 0];
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] subj_CINT[ (savesubjectmatrices && (sum(whenmat[7,1:5]) || statedep[7])) ? nsubjects : 0];
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] subj_T0VAR[ (savesubjectmatrices && (sum(whenmat[8,1:5]) || statedep[8])) ? nsubjects : 0];
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] subj_TDPREDEFFECT[ (savesubjectmatrices && (sum(whenmat[9,1:5]) || statedep[9])) ? nsubjects : 0];
      matrix[matrixdims[31, 1], matrixdims[31, 2] ] subj_DIFFUSIONcov[ (savesubjectmatrices && (sum(whenmat[31,1:5]) || statedep[31])) ? nsubjects : 0];
      matrix[matrixdims[32, 1], matrixdims[32, 2] ] subj_MANIFESTcov[ (savesubjectmatrices && (sum(whenmat[32,1:5]) || statedep[32])) ? nsubjects : 0];
      matrix[matrixdims[33, 1], matrixdims[33, 2] ] subj_T0cov[ (savesubjectmatrices && (sum(whenmat[33,1:5]) || statedep[33])) ? nsubjects : 0];
      matrix[matrixdims[21, 1], matrixdims[21, 2] ] subj_asymCINT[ (savesubjectmatrices && (sum(whenmat[21,1:5]) || statedep[21])) ? nsubjects : 0];
      matrix[matrixdims[22, 1], matrixdims[22, 2] ] subj_asymDIFFUSIONcov[ (savesubjectmatrices && (sum(whenmat[22,1:5]) || statedep[22])) ? nsubjects : 0];

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
      rawpopcovbase[j,j] = rawpopsd[j]; //used with intoverpop
      for(i in 1:nindvarying){
        if(i > j){
          counter += 1;
          rawpopcovbase[i,j]=inv_logit(sqrtpcov[counter])*2-1;
          rawpopcovbase[j,i]=0;// needed to avoid nan output;
        }
      }
    }
    //if(choleskymats==0) rawpopcorr = constraincorsqrt(rawpopcovbase);
    //if(choleskymats== -1) 
    rawpopcorr = tcrossprod( constraincorsqrt(rawpopcovbase));
    rawpopcov = makesym(quad_form_diag(rawpopcorr, rawpopsd +1e-8),verbose,1);
    rawpopcovchol = cholesky_decompose(rawpopcov); 
  }//end indvarying par setup

  {

  int prevrow=0;
  real prevdt=0;
  real dt=1; //initialise to make sure drift is computed on first pass
  real dtsmall;
  int dtchange=1;
  real prevtime=0;
  int T0check=0;
  matrix[nlatentpop, nlatentpop] etacov; //covariance of latent states

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] syprior;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypriorcov_sqrt = rep_matrix(0,nmanifest,nmanifest); 
  matrix[nmanifest, nmanifest] ycov; 
  
  matrix[nlatentpop,nlatentpop] eJAx = diag_matrix(rep_vector(1,nlatentpop)); //time evolved jacobian
  matrix[nlatentpop,nlatentpop] eJAxs[dosmoother ? ndatapoints : 1]; //time evolved jacobian, saved for smoother

  row_vector[nlatentpop] state = rep_row_vector(-999,nlatentpop); 
  row_vector[nlatentpop] statetf;
  matrix[nlatentpop,nlatentpop] JAx; //Jacobian for drift
  //matrix[nlatentpop,nlatentpop] J0; //Jacobian for t0
  matrix[nlatentpop,nlatentpop] Jtd;//diag_matrix(rep_vector(1),nlatentpop); //Jacobian for nltdpredeffect
  matrix[ nmanifest,nlatentpop] Jy;//Jacobian for measurement 
  matrix[ nmanifest,nlatentpop] Jys[dosmoother ? ndatapoints : 0];//saved Jacobian for measurement smoother
  
  

  //linear continuous time calcs
  matrix[nlatent,nlatent] discreteDRIFT;
  vector[nlatent] discreteCINT;
  matrix[nlatent,nlatent] discreteDIFFUSION = rep_matrix(0.0,nlatent,nlatent);

  
  vector[nparams] rawindparams = rawpopmeans;
  vector[nparams] indparams;
  
  matrix[nlatentpop,nlatentpop] etacovb[3,dosmoother ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ycovb[3,dosmoother ? ndatapoints : 0];
  vector[nlatentpop] etab[3,dosmoother ? ndatapoints : 0];
  vector[nmanifest] yb[3,dosmoother ? ndatapoints : 0];

  //dynamic system matrices
  
      matrix[matrixdims[10, 1], matrixdims[10, 2] ] PARS;
      matrix[matrixdims[1, 1], matrixdims[1, 2] ] T0MEANS;
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] LAMBDA;
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] DRIFT;
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] DIFFUSION;
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] MANIFESTVAR;
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] MANIFESTMEANS;
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] CINT;
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] T0VAR;
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] TDPREDEFFECT;
      matrix[matrixdims[31, 1], matrixdims[31, 2] ] DIFFUSIONcov;
      matrix[matrixdims[32, 1], matrixdims[32, 2] ] MANIFESTcov;
      matrix[matrixdims[33, 1], matrixdims[33, 2] ] T0cov;
      matrix[matrixdims[21, 1], matrixdims[21, 2] ] asymCINT;
      matrix[matrixdims[22, 1], matrixdims[22, 2] ] asymDIFFUSIONcov;
  
  asymDIFFUSIONcov = rep_matrix(0,nlatent,nlatent); //in case of derrindices need to init
  DIFFUSIONcov = rep_matrix(0,nlatent,nlatent);

  for(rowx in 0:(dokalman ? ndatapoints : 0)){
    int rowi = rowx ? rowx : 1;
    if( rowx==0 ||
      (dokalmanrows[rowi] && 
        subject[rowi] >= (firstsub - .1) &&  subject[rowi] <= (lastsub + .1))){ //if doing this row for this subject
    
    int si = rowx ? subject[rowi] : 0;
    int full = (dosmoother==1 || si ==0);
    int o[full ? nmanifest : nobs_y[rowi]]; //which obs are not missing in this row
    int o1[full ? size(whichequals(manifesttype,1,1)) : nbinary_y[rowi] ];
    int o0[full ? size(whichequals(manifesttype,1,0)) : ncont_y[rowi] ];
    
    int od[nobs_y[rowi]] = whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
    int o1d[nbinary_y[rowi] ]= whichbinary_y[rowi,1:nbinary_y[rowi]];
    int o0d[ncont_y[rowi] ]= whichcont_y[rowi,1:ncont_y[rowi]];
    
    if(!full){
      o= whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
      o1= whichbinary_y[rowi,1:nbinary_y[rowi]];
      o0= whichcont_y[rowi,1:ncont_y[rowi]];
    }
    if(full){ //needed to calculate yprior and yupd ysmooth
      for(mi in 1:nmanifest) o[mi] = mi;
      o1= whichequals(manifesttype,1,1);
      o0= whichequals(manifesttype,1,0);
    }
    
    if(prevrow != 0 && rowi != 1) T0check = (si==subject[prevrow]) ? (T0check+1) : 0; //if same subject, add one, else zero
    if(T0check > 0){
      dt = time[rowi] - time[prevrow];
      dtchange = dt!=prevdt; 
      prevdt = dt; //update previous dt store after checking for change
      //prevtime = time[rowi];
    }

    //if(dosmoother && prevrow!=0) eJAx[rowi,,] = eJAx[prevrow,,];
    
    if(T0check == 0) { // calculate initial matrices if this is first row for si
  
  rawindparams=rawpopmeans;
  
  if(si > 0 && nindvarying > 0 && intoverpop==0)  rawindparams[indvaryingindex] += rawpopcovchol * baseindparams[si];

  if(si > 0 &&  ntieffects > 0){
  if(nmissingtipreds > 0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipreds[si]';
    
    if(nmissingtipreds==0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipredsdata[si]';
  }

  // compute individual parameters that are not state dependent, either all (if si=0) or just update indvarying ones.
  indparams[whichequals(whenvecp[si ? 2 : 1], 0, 0)]= 
    parvectform(whichequals(whenvecp[si ? 2 : 1], 0, 0),rawindparams', 
    0, matsetup, matvalues, si)';
     
  if(whenmat[1, 5] >= (si ? 1 : 0)) T0MEANS = 
    mcalc(T0MEANS, indparams, statetf, {0}, 1, matsetup, matvalues, si); // base t0means to init
      
 // for(li in 1:nlatentpop) if(!is_nan(T0MEANS[li,1])) state[li] = T0MEANS[li,1]; //in case of t0 dependencies, may have missingness
  
  state=T0MEANS[,1]';
  
  statetf[whichequals(whenvecs[1],0,0)] = parvectform(whichequals(whenvecs[1],0,0),state, 1,
    matsetup, matvalues, si);   
    
  if(si==0 || sum(whenmat[10,{5,1}]) > 0 )PARS=mcalc(PARS,indparams, statetf,{0,1}, 10, matsetup, matvalues, si); 
 //initialise simple PARS then do complex PARS
  
    
  if(si==0 || sum(whenmat[1,{5,1}]) > 0 )T0MEANS=mcalc(T0MEANS,indparams, statetf,{0,1}, 1, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[8,{5,1}]) > 0 )T0VAR=mcalc(T0VAR,indparams, statetf,{0,1}, 8, matsetup, matvalues, si); 

      
  
    for(li in 1:nlatentpop) if(is_nan(state[li])) state[li] = T0MEANS[li,1]; //finish updating state
    
    //init other system matrices (already done PARS, redo t0means in case of PARS dependencies...)
   if(si==0 || sum(whenmat[2,{5}]) > 0 )LAMBDA=mcalc(LAMBDA,indparams, statetf,{0}, 2, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[3,{5}]) > 0 )DRIFT=mcalc(DRIFT,indparams, statetf,{0}, 3, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[4,{5}]) > 0 )DIFFUSION=mcalc(DIFFUSION,indparams, statetf,{0}, 4, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[5,{5}]) > 0 )MANIFESTVAR=mcalc(MANIFESTVAR,indparams, statetf,{0}, 5, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[6,{5}]) > 0 )MANIFESTMEANS=mcalc(MANIFESTMEANS,indparams, statetf,{0}, 6, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[7,{5}]) > 0 )CINT=mcalc(CINT,indparams, statetf,{0}, 7, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[8,{5}]) > 0 )T0VAR=mcalc(T0VAR,indparams, statetf,{0}, 8, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[9,{5}]) > 0 )TDPREDEFFECT=mcalc(TDPREDEFFECT,indparams, statetf,{0}, 9, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[52,{5}]) > 0 )JAx=mcalc(JAx,indparams, statetf,{0}, 52, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[53,{5}]) > 0 )Jtd=mcalc(Jtd,indparams, statetf,{0}, 53, matsetup, matvalues, si); 
if(si==0 || sum(whenmat[54,{5}]) > 0 )Jy=mcalc(Jy,indparams, statetf,{0}, 54, matsetup, matvalues, si); 

    
    if(verbose==2) print("DRIFT = ",DRIFT);
    if(verbose==2) print("indparams = ", indparams);
    
    
 if(si==0 || (sum(whenmat[8,]) + statedep[8]) > 0 ) { // this causes problems but shouldnt -- is t0var being adjusted each iteration when it shouldnt?
   if(intoverpop && nindvarying > 0) T0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopcovbase;
    T0cov = sdcovsqrt2cov(T0VAR,choleskymats); 

    if(intoverpop && nindvarying > 0){ //adjust cov matrix for transforms
    if(si==0) rawpopcovchol = cholesky_decompose(makesym(T0cov[intoverpopindvaryingindex, intoverpopindvaryingindex],verbose,1));
      for(ri in 1:size(matsetup)){
        if(matsetup[ri,7]==1){ //if t0means
          if(matsetup[ri,5]) { //and indvarying
            T0cov[matsetup[ri,1], ] *= matvalues[ri,2] * matvalues[ri,3]; //multiplier meanscale
            T0cov[, matsetup[ri,1] ] *=  matvalues[ri,2] * matvalues[ri,3]; //multiplier meanscale
          }
        }
      }
    }
 }
  
// if(nt0varstationary > 0) {
//   if(si==0 || (sum(whenmat[8,]) + statedep[8] + sum(whenmat[3,]) + statedep[3] + sum(whenmat[4,]) + statedep[4]) > 0 ){
//     matrix[nlatent,nlatent] stationarycov;
//     stationarycov = ksolve(DRIFT[derrind,derrind], sdcovsqrt2cov(DIFFUSION[derrind,derrind],choleskymats),verbose);
//   for(ri in 1:nt0varstationary){ 
//     T0cov[t0varstationary[ri,1],t0varstationary[ri,2] ] =  stationarycov[t0varstationary[ri,1],t0varstationary[ri,2] ];
//   }
//   }
// }
// 
// if(nt0meansstationary > 0){
//   if(si==0 || //on either pop pars only
//     ( (sum(whenmat[3,])+sum(whenmat[7,])+statedep[3]+statedep[7]) > 0) && savesubjectmatrices ){ // or for each subject
//     
//     if(continuoustime) asymCINT[,1] =  -DRIFT[1:nlatent,1:nlatent] \ CINT[ ,1 ];
//     if(!continuoustime) asymCINT[,1] =  add_diag(-DRIFT[1:nlatent,1:nlatent],1) \ CINT[,1 ];
// 
//     if(si==0 || (sum(whenmat[1,]) + statedep[1]) > 0) {
//       for(ri in 1:nt0meansstationary){
//         T0MEANS[t0meansstationary[ri,1]] = 
//           asymCINT[t0meansstationary[ri,1] ];
//       }
//     }
//   }
// }

    etacov=T0cov;
    } //end T0 matrices
    
if(verbose > 1) print ("below t0 row ", rowi);

      if(si==0 || (T0check>0)){ //for init or subsequent time steps when observations exist
        vector[nlatent] base;
        real intstepi = 0;
        
        dtsmall = dt / ceil(dt / maxtimestep);
        
        while(intstepi < (dt-1e-10)){
          intstepi = intstepi + dtsmall;
          
    {
    int zeroint[1];
    row_vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(JAxfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0)  state[statei] += Jstep;
      
        statetf[whichequals(whenvecs[2],0,0)] = 
          parvectform(whichequals(whenvecs[2],0,0),state, 2, matsetup, matvalues, si);
          
        //initialise PARS first, and simple PARS before complex PARS
        if(statedep[10] || whenmat[10,2]) PARS=mcalc(PARS,indparams, statetf,{2}, 10, matsetup, matvalues, si); 
        
        if(statedep[3] || whenmat[3,2]) DRIFT=mcalc(DRIFT,indparams, statetf,{2}, 3, matsetup, matvalues, si); 
        if(statedep[7] || whenmat[7,2]) CINT=mcalc(CINT,indparams, statetf,{2}, 7, matsetup, matvalues, si); 
        
      
      if(statei > 0) {
        JAx[1:nlatent,statei] =  DRIFT * state[1:nlatent]' + CINT[,1]; //compute new change
         if(verbose>1) print("JAx ",JAx);
      }
      if(statei== 0 && size(JAxfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
        base = DRIFT * state[1:nlatent]' + CINT[,1];
        if(verbose>1) print("base = ",base,"    sjaxinit= ",JAx);
        for(fi in JAxfinite){
          JAx[1:nlatent,fi] -= base;
          JAx[1:nlatent,fi] /= Jstep; //new - baseline change divided by stepsize
        }
      }
    }
    if(verbose>1) print("JAx ",JAx);
    }
    
      
      if(statedep[4] || whenmat[4,2]) DIFFUSION=mcalc(DIFFUSION,indparams, statetf,{2}, 4, matsetup, matvalues, si); 
      if(statedep[52] || whenmat[52,2]) JAx=mcalc(JAx,indparams, statetf,{2}, 52, matsetup, matvalues, si); 
      
      
      if(si==0 ||statedep[4] || whenmat[4,2] || ( T0check ==1 && whenmat[4,5])){
        DIFFUSIONcov[derrind,derrind] = sdcovsqrt2cov(DIFFUSION[derrind,derrind],choleskymats);
        if(!continuoustime) discreteDIFFUSION=DIFFUSIONcov;
      }
      
        if(continuoustime){
        
            if(si==0 || dtchange==1 || statedep[3]|| statedep[52] || whenmat[3,2] || //if first sub or changing every state
              (T0check == 1 && whenmat[3,5])){ //or first time step of new sub with ind difs
              
                //discreteDRIFT = expm2(append_row(append_col(DRIFT[1:nlatent, 1:nlatent],CINT),nlplusonezerovec') * dtsmall);
                discreteDRIFT = expm2(DRIFT * dtsmall);
                
                if(!JAxDRIFTequiv){ 
                  eJAx =  expm2(JAx * dtsmall);
                } else eJAx[1:nlatent, 1:nlatent] = discreteDRIFT;
                               
            if(si==0 || statedep[3] || statedep[4]||statedep[52]||  //if first pass or state dependent
              whenmat[4,2] || whenmat[3,2] ||
              (T0check == 1 && (whenmat[3,5]  || whenmat[4,5]))){ //or first time step of new sub with ind difs
              asymDIFFUSIONcov[derrind,derrind] = ksolve(JAx[derrind,derrind], DIFFUSIONcov[derrind,derrind],verbose);
            }
            
            discreteDIFFUSION[derrind,derrind] =  asymDIFFUSIONcov[derrind,derrind] - 
              quad_form_sym( asymDIFFUSIONcov[derrind,derrind], eJAx[derrind,derrind]' );
              
            } //end discrete coef calcs
            
            for(li in 1:nlatent) if(is_nan(state[li]) || is_nan(sum(discreteDRIFT[li,]))) {
            print("Possible time step problem? Intervals too large? Try reduce maxtimestep");
            }
            
            state[1:nlatent] *= discreteDRIFT'; // ???compute before new diffusion calcs
            
            if(size(CINTnonzero)){
            
              if(si==0 || dtchange==1 || statedep[3]|| statedep[7] || whenmat[3,2] || whenmat[7,2] || // state depenency
                (T0check == 1 && (whenmat[7,5] || whenmat[3,5]))){ //or ind difs
                discreteCINT = (DRIFT \ (discreteDRIFT-IIlatentpop[1:nlatent,1:nlatent])) * CINT[,1];
              }
              
            state[1:nlatent] += discreteCINT';
            } // end cint section
            
            if(intoverstates==1 || dosmoother==1){
              etacov = quad_form_sym(makesym(etacov,verbose,1), eJAx');
              etacov[derrind,derrind] += discreteDIFFUSION[derrind,derrind]; 
            }
              
            if(intstepi >= (dt-1e-10) && dosmoother) eJAxs[rowi,,] = expm2(JAx * dt); //save approximate exponentiated jacobian for smoothing
        } // end continuous time section
        
        if(continuoustime==0){ 
          if(dosmoother) eJAxs[rowi,,] = JAx;
          if(intoverstates==1 || dosmoother==1){
            etacov = quad_form_sym(makesym(etacov,verbose,1), JAx');
            etacov[ derrind, derrind ] += DIFFUSIONcov[ derrind, derrind ]; 
          }
          state[1:nlatent] *= DRIFT';
          state[CINTnonzero]+= CINT[CINTnonzero,1]';
        } // end non continuous time coefs
        
        } // end time step loop
      } // end non linear time update
    
    if(ntdpred > 0) {
      int nonzerotdpred = 0;
      for(tdi in 1:ntdpred) if(tdpreds[rowi,tdi] != 0.0) nonzerotdpred = 1;
      if(si==0 ||nonzerotdpred){
      
        statetf[whichequals(whenvecs[3],0,0)] = 
          parvectform( whichequals(whenvecs[3],0,0), state, 3, matsetup, matvalues, si);
          
        //initialise PARS first, and simple PARS before complex PARS
        if(statedep[10] || whenmat[10,3]) PARS=mcalc(PARS,indparams, statetf,{3}, 10, matsetup, matvalues, si); 
        
        
      
        if(statedep[9] || whenmat[9,3]) TDPREDEFFECT=mcalc(TDPREDEFFECT,indparams, statetf,{3}, 9, matsetup, matvalues, si); 
        if(statedep[53] || whenmat[53,3]) Jtd=mcalc(Jtd,indparams, statetf,{3}, 53, matsetup, matvalues, si); 

        

        state[1:nlatent] +=   (TDPREDEFFECT * tdpreds[rowi])'; //tdpred effect only influences at observed time point
        etacov = quad_form_sym(makesym(etacov,verbose,1),Jtd');  //could speed up by detecting if non diagonal Jtd
      }
    }//end nonlinear tdpred

  if(si > 0 && intoverstates==0){ //unused states if intoverpop is specified, consider fixing...
    if(T0check==0) state += (cholesky_decompose(etacov) * etaupdbasestates[(1+(rowi-1)*nlatentpop):(rowi*nlatentpop)])';
    if(T0check>0) state[derrind] +=  (cholesky_decompose(makesym(discreteDIFFUSION[derrind,derrind],verbose,1)) * 
     (etaupdbasestates[(1+(rowi-1)*nlatentpop):(nlatent+(rowi-1)*nlatentpop)])[derrind])';
     
   // if(T0check==0) llrow[rowi]+= multi_normal_cholesky_lpdf(
  //     etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)] | rep_vector(0,nlatent), etacov);
  //  if(T0check>0) llrow[rowi]+= multi_normal_lpdf(
  //     etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)] | rep_vector(0,nlatent), discreteDIFFUSION);
  //  state+=etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
  }

if(verbose > 1){
  print("etaprior = ", state);
  print("etapriorcov = ", etacov);
}

if(dosmoother){
  etacovb[1,rowi] = etacov; 
  etab[1,rowi] = state';
}

//measurement update

  if(si == 0 || nobs_y[rowi] > 0 || dosmoother){ //measurement init
    int zeroint[1];
    row_vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(Jyfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0 && (dosmoother + intoverstates) > 0)  state[statei] += Jstep;
      
            
        statetf[whichequals(whenvecs[4],0,0)] = 
          parvectform( whichequals(whenvecs[4],0,0), state, 4, matsetup, matvalues, si);
          
        //initialise PARS first, and simple PARS before complex PARS
        if(statedep[10] || whenmat[10,4]) PARS=mcalc(PARS,indparams, statetf,{4}, 10, matsetup, matvalues, si); 
        
        
      
        if(statedep[2] || whenmat[2,4]) LAMBDA=mcalc(LAMBDA,indparams, statetf,{4}, 2, matsetup, matvalues, si); 
        if(statedep[5] || whenmat[5,4]) MANIFESTVAR=mcalc(MANIFESTVAR,indparams, statetf,{4}, 5, matsetup, matvalues, si); 
        if(statedep[6] || whenmat[6,4]) MANIFESTMEANS=mcalc(MANIFESTMEANS,indparams, statetf,{4}, 6, matsetup, matvalues, si); 
        if(statedep[54] || whenmat[54,4]) Jy=mcalc(Jy,indparams, statetf,{4}, 54, matsetup, matvalues, si); 

        
        
      if(statei > 0 && (intoverstates) > 0) {
        Jy[o,statei] =  LAMBDA[o] * state[1:nlatent]' + MANIFESTMEANS[o,1]; //compute new change
        Jy[o1,statei] = to_vector(inv_logit(to_array_1d(Jy[o1,statei])));
        if(verbose>1) print("Jy ",Jy);
      }
      if(statei==0){
        syprior[o] = LAMBDA[o] * state[1:nlatent]' + MANIFESTMEANS[o,1];
        syprior[o1] = to_vector(inv_logit(to_array_1d( syprior[o1] )));
        if(size(Jyfinite) ) { //only need these calcs if there are finite differences to do -- otherwise loop just performs system calcs.
          if(verbose>1) print("syprior = ",syprior,"    Jyinit= ",Jy);
          for(fi in Jyfinite){
            Jy[o,fi] -= syprior[o];
            Jy[o,fi] /= Jstep; //new - baseline change divided by stepsize
          }
        }
      }
    }
    if(verbose>1) print("Jy ",Jy);
  } //end measurement init
      
  if(si==0 || whenmat[5,5] || whenmat[5,4] || statedep[5]) MANIFESTcov = sdcovsqrt2cov(MANIFESTVAR,choleskymats);

 
  if(si > 0 && (nobs_y[rowi] > 0 || dosmoother)){   //if not just inits...

      if(intoverstates==1 || dosmoother==1) { //classic kalman
        ycov[o,o] = quad_form_sym(makesym(etacov,verbose,1), Jy[o,]') + MANIFESTcov[o,o]; // previously shifted measurement error down, but reverted
        for(wi in 1:nmanifest){ 
          // if(Y[rowi,wi] != 99999 || dosmoother==1) ycov[wi,wi] += square(MANIFESTVAR[wi,wi]);
          if(manifesttype[wi]==1 && (Y[rowi,wi] != 99999  || dosmoother==1)) ycov[wi,wi] += fabs((syprior[wi] - 1) .* (syprior[wi]));
          if(manifesttype[wi]==2 && (Y[rowi,wi] != 99999  || dosmoother==1)) ycov[wi,wi] += square(fabs((syprior[wi] - round(syprior[wi])))); 
        }
      }
        
      if(intoverstates==0 && ncont_y[rowi] > 0) ypriorcov_sqrt[o,o] = 
        cholesky_decompose(makesym(MANIFESTcov[o,o],verbose,1));

        
     
err[od] = Y[rowi,od] - syprior[od]; // prediction error
    
      if(intoverstates==1 && size(od) > 0) {
        
         if(verbose > 1) print("before K rowi =",rowi, "  si =", si, "  state =",state, "  statetf = ", statetf, "  etacov ",etacov,
          " indparams = ", indparams,
            "  syprior[o] =",syprior[o],"  ycov[o,o] ",ycov[o,o], 
            "  PARS = ", PARS, 
            "  DRIFT =", DRIFT, " DIFFUSION =", DIFFUSION, 
            " CINT =", CINT, "  discreteCINT = ", discreteCINT, "  MANIFESTcov ", (MANIFESTcov), "  MANIFESTMEANS ", MANIFESTMEANS, 
            "  T0cov", T0cov,  " T0MEANS ", T0MEANS, "LAMBDA = ", LAMBDA, "  Jy = ",Jy,
            " discreteDRIFT = ", discreteDRIFT, "  discreteDIFFUSION ", discreteDIFFUSION, "  asymDIFFUSIONcov ", asymDIFFUSIONcov, 
            " DIFFUSIONcov = ", DIFFUSIONcov,
            " eJAx = ", eJAx,
            "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
         
        K[,od] = mdivide_right_spd(etacov * Jy[od,]', makesym(ycov[od,od],verbose,1)); // * multiply_lower_tri_self_transpose(ycovi');// ycov[od,od]; 
        etacov += -K[,od] * Jy[od,] * etacov; //cov update
        state +=  (K[,od] * err[od])'; //state update
      }
      
      if(dosmoother==1) {
        yb[1,rowi] = syprior[o];
        etab[2,rowi] = state';
        ycovb[1,rowi] = ycov;
        etacovb[2,rowi] = etacov;
        ycovb[2,rowi] = quad_form_sym(makesym(etacov,verbose,1), Jy') + MANIFESTcov;
        yb[2,rowi] = MANIFESTMEANS[o,1] + LAMBDA[o,] * state[1:nlatent]';
        Jys[rowi,,] = Jy;
      }
      
      
      if(verbose > 1) print(" After K rowi =",rowi, "  si =", si, "  state =",state,"  etacov ",etacov,"  K[,o] ",K[,o]);
        
  //likelihood stuff
      if(nbinary_y[rowi] > 0) llrow[rowi] += sum(log(1e-10+Y[rowi,o1d] .* (syprior[o1d]) + (1-Y[rowi,o1d]) .* (1-syprior[o1d]))); 

      if(size(o0d) > 0 && (llsinglerow==0 || llsinglerow == rowi)){
        if(intoverstates==1) ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(ycov[o0d,o0d]); //removed makesym
         llrow[rowi] +=  multi_normal_cholesky_lpdf(Y[rowi,o0d] | syprior[o0d], ypriorcov_sqrt[o0d,o0d]);
         //errtrans[counter:(counter + ncont_y[rowi]-1)] = 
           //mdivide_left_tri_low(ypriorcov_sqrt[o0d,o0d], err[o0d]); //transform pred errors to standard normal dist and collect
         //ll+= -sum(log(diagonal(ypriorcov_sqrt[o0d,o0d]))); //account for transformation of scale in loglik
         //counter += ncont_y[rowi];
      }
      if(verbose > 1) print(llrow[rowi]);
      
    }//end si > 0 nobs > 0 section
    
       // store system matrices
       
    if(si==0 || //on either pop pars only
      (  (sum(whenmat[3,])+sum(whenmat[7,])+statedep[3]+statedep[7]) > 0 && savesubjectmatrices) ){ // or for each subject
      if(continuoustime==1) asymCINT[,1] =  -DRIFT[1:nlatent,1:nlatent] \ CINT[ ,1 ];
      if(continuoustime==0) asymCINT[,1] =  add_diag(-DRIFT[1:nlatent,1:nlatent],1) \ CINT[,1 ];
    }

  
  if(!continuoustime){
    if(si==0 || //on either pop pars only
    (  (sum(whenmat[3,])+sum(whenmat[4,])+statedep[3]+statedep[4]) > 0 && savesubjectmatrices) ){ // or for each subject
  
      asymDIFFUSIONcov[ derrind, derrind ] = 
        to_matrix( (add_diag( -sqkron_prod(JAx[ derrind, derrind ], JAx[ derrind, derrind ]),1)) \  
          to_vector(DIFFUSIONcov[ derrind, derrind ]), ndiffusion, ndiffusion);
    }
  }
      
    
  if(si == 0){
pop_PARS = PARS; pop_T0MEANS = T0MEANS; pop_LAMBDA = LAMBDA; pop_DRIFT = DRIFT; pop_DIFFUSION = DIFFUSION; pop_MANIFESTVAR = MANIFESTVAR; pop_MANIFESTMEANS = MANIFESTMEANS; pop_CINT = CINT; pop_T0VAR = T0VAR; pop_TDPREDEFFECT = TDPREDEFFECT; pop_DIFFUSIONcov = DIFFUSIONcov; pop_MANIFESTcov = MANIFESTcov; pop_T0cov = T0cov; pop_asymCINT = asymCINT; pop_asymDIFFUSIONcov = asymDIFFUSIONcov; 
  }

  
  if(si > 0 && dosmoother && (rowi==ndatapoints || subject[rowi+1] != subject[rowi])){ //at subjects last datapoint, smooth
    int sri = rowi;
    while(sri>0 && subject[sri]==si){
      if(sri==rowi) {
        etab[3,sri]=etab[2,sri];
        yb[3,sri]=yb[2,sri];
        etacovb[3,sri]=etacovb[2,sri];
        ycovb[3,sri]=ycovb[2,sri];
      } else{
        matrix[nlatentpop,nlatentpop] smoother;
        smoother = etacovb[2,sri] * eJAxs[sri+1,,]' / makesym(etacovb[1,sri+1],verbose,1);
        etab[3,sri,]= etab[2,sri,] + smoother * (etab[3,sri+1,] - etab[1,sri+1,]);
        etacovb[3,sri]= etacovb[2,sri] + smoother * ( etacovb[3,sri+1] - etacovb[1,sri+1]) * smoother';
        yb[3,sri,] = yb[2,sri,] + Jys[sri,,] * (etab[3,sri,] - etab[2,sri,]);
        ycovb[3,sri] = ycovb[2,sri] + Jys[sri,,] * (etacovb[3,sri] - etacovb[2,sri]) * Jys[sri,,]';
      }
      state=etab[3,sri,]'; //update for t0means saving
      sri += -1;
      while(sri > 0 && dokalmanrows[sri]==0) sri+= -1; //skip rows if requested

      if(savesubjectmatrices && //if getting subj matrices and 
        (sri ? subject[sri] != subject[sri+1] : 1)){ //no more rows, or change of subject
     
    if(sum(whenmat[10,1:5]) > 0 || statedep[10]) subj_PARS[si] = PARS;
    if(sum(whenmat[1,1:5]) > 0 || statedep[1]) subj_T0MEANS[si] = T0MEANS;
    if(sum(whenmat[2,1:5]) > 0 || statedep[2]) subj_LAMBDA[si] = LAMBDA;
    if(sum(whenmat[3,1:5]) > 0 || statedep[3]) subj_DRIFT[si] = DRIFT;
    if(sum(whenmat[4,1:5]) > 0 || statedep[4]) subj_DIFFUSION[si] = DIFFUSION;
    if(sum(whenmat[5,1:5]) > 0 || statedep[5]) subj_MANIFESTVAR[si] = MANIFESTVAR;
    if(sum(whenmat[6,1:5]) > 0 || statedep[6]) subj_MANIFESTMEANS[si] = MANIFESTMEANS;
    if(sum(whenmat[7,1:5]) > 0 || statedep[7]) subj_CINT[si] = CINT;
    if(sum(whenmat[8,1:5]) > 0 || statedep[8]) subj_T0VAR[si] = T0VAR;
    if(sum(whenmat[9,1:5]) > 0 || statedep[9]) subj_TDPREDEFFECT[si] = TDPREDEFFECT;
    if(sum(whenmat[31,1:5]) > 0 || statedep[31]) subj_DIFFUSIONcov[si] = DIFFUSIONcov;
    if(sum(whenmat[32,1:5]) > 0 || statedep[32]) subj_MANIFESTcov[si] = MANIFESTcov;
    if(sum(whenmat[33,1:5]) > 0 || statedep[33]) subj_T0cov[si] = T0cov;
    if(sum(whenmat[21,1:5]) > 0 || statedep[21]) subj_asymCINT[si] = asymCINT;
    if(sum(whenmat[22,1:5]) > 0 || statedep[22]) subj_asymDIFFUSIONcov[si] = asymDIFFUSIONcov;
     if(sum(whenmat[1,1:5]) > 0 || statedep[1]) subj_T0MEANS[si,,1] = state'; //t0means updated, other pars as per final time point
        }
      
    }
  } //end smoother

  
 } // end si loop (includes sub 0)
  
  prevrow = rowi; //update previous row marker only after doing necessary calcs
}//end active rowi

if(savescores){
  ya=yb;
  ycova=ycovb;
  etaa=etab;
  etacova=etacovb;
}
ll+=sum(llrow);


  }
}
      
model{
  real priormod2 = priormod / nsubsets;
  if(intoverpop==0 && nindvarying > 0) target+= multi_normal_cholesky_lpdf(baseindparams | rep_vector(0,nindvarying), IIlatentpop[1:nindvarying,1:nindvarying]);

  if(ntipred > 0){ 
    if(nopriors==0 && laplacetipreds==0) target+= priormod2 * normal_lpdf(tipredeffectparams / tipredeffectscale| 0, 1);
    if(nopriors==0 && laplacetipreds==1) for(i in 1:ntipredeffects) target+= priormod2 * double_exponential_lpdf(pow(fabs(tipredeffectparams[i]),1+.1/((tipredeffectparams[i]*100)^2+.1)) / tipredeffectscale| 0, 1);
    target+= normal_lpdf(tipredsimputed| 0, tipredsimputedscale); //consider better handling of this when using subset approach
  }

  if(nopriors==0){ //if split files over subjects, just compute priors once
    for(i in 1:nparams){
      if(laplaceprior[i]==1) target+= priormod2 * double_exponential_lpdf(pow(fabs(rawpopmeans[i]) ,1+.1/((rawpopmeans[i]*100)^2+.1))|0,1);
    }
  }

  if(nopriors==0 && !laplaceprioronly){ //if split files over subjects, just compute priors once
  for(i in 1:nparams){
    if(laplaceprior[i]==0) target+= priormod2 * normal_lpdf(rawpopmeans[i]|0,1);
  }
  
    if(nindvarying > 0){
      if(nindvarying >1) target+= priormod2 * normal_lpdf(sqrtpcov | 0, 1);
      target+= priormod2 * normal_lpdf(rawpopsdbase | 0,1);
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
        x[ri,] = (rawpopcovchol * 
          to_vector(normal_rng(rawpopmeans[indvaryingindex],rep_vector(1,nindvarying))) )';
      }
    }
    
    for(pi in 1:nparams){
      int found=0;
      int pr1;
      int pr2;
      real rawpoppar = rawpopmeans[pi];
      while(!found){ //currently seems useless, instead just references last match if multiple
        for(ri in 1:size(matsetup)){
          if(matsetup[ri,3]==pi && matsetup[ri,8]<=0) { //if a free parameter 
            pr1 = ri; 
            pr2=ri;// unless intoverpop, pop matrix row reference is simply current row
            found=1;
            if(intoverpop && matsetup[ri,5]) { //check if shifted
              for(ri2 in 1:size(matsetup)){ //check when state reference param of matsetup corresponds to row of t0means in current matsetup row
                if(matsetup[ri2,8]  && matsetup[ri2,3] == matsetup[ri,1] && 
                matsetup[ri2,3] > nlatent && matsetup[ri2,7] < 20) pr2 = ri2; //if param is dynamic and matches row (state ref) and is not in jacobian
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
