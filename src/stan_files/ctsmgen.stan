
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
      if(pd ==1 && mat[coli,coli] < 1e-5){
        out[coli,coli] = 1e-5;// 
      } else out[coli,coli] = mat[coli,coli]; 
      for(rowi in coli:rows(mat)){
        if(rowi > coli) {
          out[rowi,coli] = mat[rowi,coli];
          out[coli,rowi] = mat[rowi,coli];
        }
        if(is_nan(out[rowi,coli])){
          if(verbose > 0) print("nan during makesym row ", rowi, " col ", coli);
          if(rowi==coli) out[rowi,coli] = 99999;
          if(rowi!=coli) {
            out[rowi,coli] = 0;
            out[coli,rowi] = 0;
          }
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
  
  
  matrix mcalc(matrix matin, vector pars, int when, int m, int[,] ms, real[,] mval, int subi, int[] subindices){
    matrix[rows(matin),cols(matin)] matout = matin;

    for(ri in 1:size(ms)){ //for each row of matrix setup
      if(m==ms[ri,7] && ms[ri,8]==when){ // if correct matrix and when -- check this for t0var, maybe restrictive
      
        if(subi ==0 ||  //if population parameter
          ( ms[ri,7] == 8 && subindices[8]) || //or a covariance parameter in an individually varying matrix
          (ms[ri,3] > 0 && (ms[ri,5] > 0 || ms[ri,6] > 0 || ms[ri,8] > 0)) //or there is individual variation
        ){ //otherwise repeated values
        
          if( (when > 0 && ms[ri,8]>0) || 
              ( (when== -1 || when == 0) && ms[ri,8]==0) ){ //if doing statecalcs do them, if doing static calcs do them
              
            if(ms[ri,3] > 0)  matout[ms[ri,1], ms[ri,2] ] = pars[ms[ri,3]]; //already tformed
            
            if(ms[ri,3] < 1) matout[ms[ri,1], ms[ri,2] ] = mval[ri, 1]; //doing this once over all subjects unless covariance matrix -- speed ups possible here, check properly!
            
            if(ms[ri,9] > 0){ //then get par from other row
              for(ri2 in 1:size(ms)){
                if(ms[ri2,9] == -ri) matout[ms[ri,1], ms[ri,2] ] = pars[ms[ri2,3]];
              }
            }
          }
        }
      }
    }
  return(matout);
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
  int drcintoffdiag[nlatent+1];

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
  int whenvecp[2,nparams];
  int whenvecs[4,nlatentpop];
  int matrixdims[10,2];
  int savescores;
  int savesubjectmatrices;
  int fixedsubpars;
  vector[fixedsubpars ? nindvarying : 0] fixedindparams[fixedsubpars ? nsubjects : 0];
  int dokalman;
  int dokalmanrowsdata[ndatapoints];
  real Jstep;
  real dokalmanpriormodifier;
  int intoverpopindvaryingindex[intoverpop ? nindvarying : 0];
  int nsJAxfinite;
  int sJAxfinite[nsJAxfinite];
  int nJyfinite;
  int sJyfinite[nJyfinite];
  int taylorheun;
  int difftype;
  int jacoffdiag[nlatentpop];
  int njacoffdiagindex;
  int jacoffdiagindex[njacoffdiagindex];
  int popcovn;
  int llsinglerow;
  int doonesubject;
}
      
transformed data{
  matrix[nlatent+nindvarying,nlatent+nindvarying] IIlatentpop = diag_matrix(rep_vector(1,nlatent+nindvarying));
  vector[nlatentpop-nlatent] nlpzerovec = rep_vector(0,nlatentpop-nlatent);
  vector[nlatent+1] nlplusonezerovec = rep_vector(0,nlatent+1);
  int nsubjects2 = doonesubject ? 1 : nsubjects;
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
  vector[fixedsubpars ? 0 : (intoverpop ? 0 : nindvarying)] baseindparams[fixedsubpars ? 0 : (intoverpop ? 0 :  (doonesubject ? 1 : nsubjects2) )]; //vector of subject level deviations, on the raw scale
  
  vector[ntipredeffects] tipredeffectparams; // effects of time independent covariates
  vector[nmissingtipreds] tipredsimputed;
  //vector[ (( (ntipredeffects-1) * (1-nopriors) ) > 0) ? 1 : 0] tipredglobalscalepar;
  
  vector[intoverstates ? 0 : nlatentpop*ndatapoints] etaupdbasestates; //sampled latent states posterior
  real onesubject[doonesubject ? doonesubject : 0]; //allows multiple specific
}
      
transformed parameters{
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying, nindvarying] rawpopc[4];



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

  }
}
      
model{
  if(doonesubject==0 ||onesubject[1] > .5){ 
    if(intoverpop==0 && fixedsubpars == 1 && nindvarying > 0) target+= multi_normal_cholesky_lpdf(fixedindparams | rep_vector(0,nindvarying),IIlatentpop[1:nindvarying,1:nindvarying]);
    if(intoverpop==0 && fixedsubpars == 0 && nindvarying > 0) target+= multi_normal_cholesky_lpdf(baseindparams | rep_vector(0,nindvarying), IIlatentpop[1:nindvarying,1:nindvarying]);
  }
  if(doonesubject==0 ||onesubject[1] < .5){ 
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
  }
  
  
  
  
  if(verbose > 0) print("lp = ", target());
}
generated quantities{
  vector[nparams] popmeans;
  vector[nindvarying] popsd; // = rep_vector(0,nparams);
  matrix[nindvarying,nindvarying] popcov;
  matrix[nparams,ntipred] linearTIPREDEFFECT;

  real ll = 0;
  vector[savescores ? ndatapoints : 1] llrow = rep_vector(0.0,savescores ? ndatapoints : 1);
  matrix[nlatentpop,nlatentpop] etacova[3,savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ycova[3,savescores ? ndatapoints : 0];
  vector[nlatentpop] etaa[3,savescores ? ndatapoints : 0];
  vector[nmanifest] ya[3,savescores ? ndatapoints : 0];
  vector[nmanifest] Ygen[ndatapoints];
  //vector[nlatentpop*ndatapoints] etaupdbasestates = to_vector(normal_rng(rep_array(0.0,ndatapoints*nlatentpop),rep_array(1.0,ndatapoints*nlatentpop))); //sampled latent states posterior
     matrix[matrixdims[10, 1], matrixdims[10, 2] ] PARS[subindices[10]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
      matrix[matrixdims[1, 1], matrixdims[1, 2] ] T0MEANS[subindices[1]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
      matrix[matrixdims[2, 1], matrixdims[2, 2] ] LAMBDA[subindices[2]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
      matrix[matrixdims[3, 1], matrixdims[3, 2] ] DRIFT[subindices[3]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
      matrix[matrixdims[4, 1], matrixdims[4, 2] ] DIFFUSION[subindices[4]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
      matrix[matrixdims[5, 1], matrixdims[5, 2] ] MANIFESTVAR[subindices[5]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
      matrix[matrixdims[6, 1], matrixdims[6, 2] ] MANIFESTMEANS[subindices[6]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
      matrix[matrixdims[7, 1], matrixdims[7, 2] ] CINT[subindices[7]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
      matrix[matrixdims[8, 1], matrixdims[8, 2] ] T0VAR[subindices[8]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
      matrix[matrixdims[9, 1], matrixdims[9, 2] ] TDPREDEFFECT[subindices[9]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1];

  matrix[nlatent,nlatent] asymDIFFUSION[(subindices[3] + subindices[4])  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; //stationary latent process variance
  vector[nlatent] asymCINT[(subindices[3] + subindices[7])  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; // latent process asymptotic level
matrix[nlatent, nlatent] DIFFUSIONcov[subindices[4] ? (savesubjectmatrices ? nsubjects2 : 1) : 1];
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

  matrix[nlatent,nlatent] pop_asymDIFFUSION; //stationary latent process variance
  vector[nlatent] pop_asymCINT; // latent process asymptotic level
matrix[nlatent, nlatent] pop_DIFFUSIONcov;

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


{
  vector[nmanifest] Ygenbase[ndatapoints];
  Ygen = rep_array(rep_vector(99999,nmanifest),ndatapoints);
  for(mi in 1:nmanifest){
    if(manifesttype[mi]==0 || manifesttype[mi]==2) {
      Ygenbase[1:ndatapoints,mi] = normal_rng(rep_vector(0,ndatapoints),rep_vector(1,ndatapoints));
    }
    if(manifesttype[mi]==1){
      Ygenbase[1:ndatapoints,mi] =  uniform_rng(rep_vector(0,ndatapoints),rep_vector(1,ndatapoints));
    }
  }
{

  int si = 0;
  int prevrow=0;
  real prevdt=0;
  real dt;
  real dtsmall;
  int dtchange=1;
  int T0check=0;
  int subjectcount = 0;
  int counter = 1;
  matrix[nlatentpop, nlatentpop] etacov; //covariance of latent states

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] syprior;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypriorcov_sqrt; 
  matrix[nmanifest, nmanifest] ycov; 
  
  matrix[nlatentpop,nlatentpop] Je[savescores ? ndatapoints : 1]; //time evolved jacobian, saved for smoother

  vector[nlatentpop] state = rep_vector(-999,nlatentpop); 
  vector[nlatentpop] statetform;
  matrix[nlatentpop,nlatentpop] sJAx; //Jacobian for drift
  matrix[nlatentpop,nlatentpop] sJ0; //Jacobian for t0
  matrix[nlatentpop,nlatentpop] sJtd;//diag_matrix(rep_vector(1),nlatentpop); //Jacobian for nltdpredeffect
  matrix[ nmanifest,nlatentpop] sJy;//Jacobian for measurement 
  
  matrix[taylorheun ? nlatentpop : 0, taylorheun ? nlatentpop : 0] Kth = rep_matrix(0.0,taylorheun ? nlatentpop : 0, taylorheun ? nlatentpop : 0); 
  matrix[taylorheun ? nlatentpop : 0, taylorheun ? nlatentpop : 0] Mth = Kth;

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

  matrix[nlatent,nlatent] sasymDIFFUSION; //stationary latent process variance
  vector[nlatent] sasymCINT; // latent process asymptotic level
matrix[nlatent, nlatent] sDIFFUSIONcov;
  
  int dokalmanrows[ndatapoints] = dokalmanrowsdata;
  
  if(doonesubject==1){
    dokalmanrows=rep_array(0,ndatapoints);
    for(i in 1:ndatapoints){
      for(subi in 1:size(onesubject)){
        if(fabs(subject[i]-onesubject[subi]) < .5){
          dokalmanrows[i] = 1; 
        }
      }
    }
  }

  for(rowi in 1:(dokalman ? ndatapoints :1)){
  if(dokalmanrows[rowi] ==1) { //used for subset selection
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
    
  
    si = subject[rowi];
    if(prevrow != 0) T0check = (si==subject[prevrow]) ? (T0check+1) : 0; //if same subject, add one, else zero
    if(T0check > 0){
      dt = time[rowi] - time[prevrow];
      dtchange = dt!=prevdt; 
      prevdt = dt; //update previous dt store after checking for change
    }
    //if(savescores || prevrow==0) Je[savescores ? rowi : 1] = IIlatentpop[1:nlatentpop,1:nlatentpop]; //elements updated later
    if(savescores && prevrow!=0) Je[rowi,,] = Je[prevrow,,];
       
    if(T0check == 0) { // calculate initial matrices if this is first row for si
  
    
  for(subi in (((prevrow!=0) ? si : 0)):si){
  
  rawindparams=rawpopmeans;
  
  if(subi > 0 && nindvarying > 0 && intoverpop==0) {
    if(fixedsubpars==0) rawindparams[indvaryingindex] += rawpopc[2] * baseindparams[doonesubject ? 1 : subi];
    if(fixedsubpars==1) rawindparams[indvaryingindex] += rawpopc[2] * fixedindparams[doonesubject ? 1 : subi];
  }

  if(subi > 0 &&  ntieffects > 0){
  if(nmissingtipreds > 0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipreds[subi]';
    if(nmissingtipreds==0) rawindparams[tieffectindices[1:ntieffects]] += 
    TIPREDEFFECT[tieffectindices[1:ntieffects]] *  tipredsdata[subi]';
  }

    indparams[whichequals(whenvecp[subi ? 2 : 1], 0, 0)]= 
      parvectform(size(whichequals(whenvecp[subi ? 2 : 1], 0, 0)),rawindparams, 
      0, matsetup, matvalues, subi, subindices, whenvecp[subi ? 2 : 1]);
      
    if(whenmat[1, 5] >= (subi ? 1 : 0)) sT0MEANS = mcalc(sT0MEANS, indparams, 0, 1, matsetup, matvalues, subi, subindices);
    state=sT0MEANS[,1];
    ;
; 

    
      for(wheni in 1:4) {
        statetform[whichequals(whenvecs[wheni],0,0)] = 
          parvectform( size(whichequals(whenvecs[wheni],0,0)),state, 
          wheni, matsetup, matvalues, subi, subindices, whenvecs[wheni]);
          
        if(whenmat[10,wheni])sPARS=mcalc(sPARS,statetform,wheni, 10, matsetup, matvalues, subi, subindices); 
if(whenmat[1,wheni])sT0MEANS=mcalc(sT0MEANS,statetform,wheni, 1, matsetup, matvalues, subi, subindices); 
if(whenmat[2,wheni])sLAMBDA=mcalc(sLAMBDA,statetform,wheni, 2, matsetup, matvalues, subi, subindices); 
if(whenmat[3,wheni])sDRIFT=mcalc(sDRIFT,statetform,wheni, 3, matsetup, matvalues, subi, subindices); 
if(whenmat[4,wheni])sDIFFUSION=mcalc(sDIFFUSION,statetform,wheni, 4, matsetup, matvalues, subi, subindices); 
if(whenmat[5,wheni])sMANIFESTVAR=mcalc(sMANIFESTVAR,statetform,wheni, 5, matsetup, matvalues, subi, subindices); 
if(whenmat[6,wheni])sMANIFESTMEANS=mcalc(sMANIFESTMEANS,statetform,wheni, 6, matsetup, matvalues, subi, subindices); 
if(whenmat[7,wheni])sCINT=mcalc(sCINT,statetform,wheni, 7, matsetup, matvalues, subi, subindices); 
if(whenmat[8,wheni])sT0VAR=mcalc(sT0VAR,statetform,wheni, 8, matsetup, matvalues, subi, subindices); 
if(whenmat[9,wheni])sTDPREDEFFECT=mcalc(sTDPREDEFFECT,statetform,wheni, 9, matsetup, matvalues, subi, subindices); 
if(whenmat[51,wheni])sJ0=mcalc(sJ0,statetform,wheni, 51, matsetup, matvalues, subi, subindices); 
if(whenmat[52,wheni])sJAx=mcalc(sJAx,statetform,wheni, 52, matsetup, matvalues, subi, subindices); 
if(whenmat[53,wheni])sJtd=mcalc(sJtd,statetform,wheni, 53, matsetup, matvalues, subi, subindices); 
if(whenmat[54,wheni])sJy=mcalc(sJy,statetform,wheni, 54, matsetup, matvalues, subi, subindices); 

        } 
        if(whenmat[10,5] || subi==0) sPARS=mcalc(sPARS,indparams,0, 10, matsetup, matvalues, subi, subindices); 
  if(whenmat[1,5] || subi==0) sT0MEANS=mcalc(sT0MEANS,indparams,0, 1, matsetup, matvalues, subi, subindices); 
  if(whenmat[2,5] || subi==0) sLAMBDA=mcalc(sLAMBDA,indparams,0, 2, matsetup, matvalues, subi, subindices); 
  if(whenmat[3,5] || subi==0) sDRIFT=mcalc(sDRIFT,indparams,0, 3, matsetup, matvalues, subi, subindices); 
  if(whenmat[4,5] || subi==0) sDIFFUSION=mcalc(sDIFFUSION,indparams,0, 4, matsetup, matvalues, subi, subindices); 
  if(whenmat[5,5] || subi==0) sMANIFESTVAR=mcalc(sMANIFESTVAR,indparams,0, 5, matsetup, matvalues, subi, subindices); 
  if(whenmat[6,5] || subi==0) sMANIFESTMEANS=mcalc(sMANIFESTMEANS,indparams,0, 6, matsetup, matvalues, subi, subindices); 
  if(whenmat[7,5] || subi==0) sCINT=mcalc(sCINT,indparams,0, 7, matsetup, matvalues, subi, subindices); 
  if(whenmat[8,5] || subi==0) sT0VAR=mcalc(sT0VAR,indparams,0, 8, matsetup, matvalues, subi, subindices); 
  if(whenmat[9,5] || subi==0) sTDPREDEFFECT=mcalc(sTDPREDEFFECT,indparams,0, 9, matsetup, matvalues, subi, subindices); 
  if(whenmat[51,5] || subi==0) sJ0=mcalc(sJ0,indparams,0, 51, matsetup, matvalues, subi, subindices); 
  if(whenmat[52,5] || subi==0) sJAx=mcalc(sJAx,indparams,0, 52, matsetup, matvalues, subi, subindices); 
  if(whenmat[53,5] || subi==0) sJtd=mcalc(sJtd,indparams,0, 53, matsetup, matvalues, subi, subindices); 
  if(whenmat[54,5] || subi==0) sJy=mcalc(sJy,indparams,0, 54, matsetup, matvalues, subi, subindices); 

    //
  
  //print("state = ", state, "   statetform = ",statetform, "   DRIFT= ", sDRIFT);
  //print("rawindparams = ", rawindparams, "    indparams = ", indparams,"   DRIFT= ", sDRIFT);
    

  // perform any whole matrix transformations, nonlinear calcs based on t0 in order to fill matrices
  
      {
  ;  
  
  } 
  ;
;
  ;
;
  ;
;
  ;
;
  
  if(subi <= (subindices[4] ? nsubjects2 : 0)) {
    sDIFFUSIONcov = sdcovsqrt2cov(sDIFFUSION,choleskymats);
  }
  if(subi <= ((subindices[3] + subindices[4])  ? nsubjects2 : 0)) {
    if(ndiffusion < nlatent) sasymDIFFUSION = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

    if(continuoustime==1) sasymDIFFUSION[ derrind, derrind] = to_matrix( 
    -sqkron_sumii(sDRIFT[ derrind, derrind ]) \  to_vector( 
         sDIFFUSIONcov[ derrind, derrind ]), ndiffusion,ndiffusion);

    if(continuoustime==0) sasymDIFFUSION[ derrind, derrind ] = to_matrix( (add_diag( 
      -sqkron_prod(sDRIFT[ derrind, derrind ], sDRIFT[ derrind, derrind ]),1)) \  to_vector(sDIFFUSIONcov[ derrind, derrind ]), ndiffusion, ndiffusion);
  } //end asymdiffusion loops

     if(subi <= (subindices[8] ? nsubjects2 : 0)) {
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
    
    if(subi <= ((subindices[3] + subindices[7])  ? nsubjects2 : 0)){
      if(continuoustime==1) sasymCINT =  -sDRIFT[1:nlatent,1:nlatent] \ sCINT[ ,1 ];
      if(continuoustime==0) sasymCINT =  add_diag(-sDRIFT[1:nlatent,1:nlatent],1) \ sCINT[,1 ];
    }
    
       
    

    
  
  if(subi == 0){
pop_PARS = sPARS; pop_T0MEANS = sT0MEANS; pop_LAMBDA = sLAMBDA; pop_DRIFT = sDRIFT; pop_DIFFUSION = sDIFFUSION; pop_MANIFESTVAR = sMANIFESTVAR; pop_MANIFESTMEANS = sMANIFESTMEANS; pop_CINT = sCINT; pop_T0VAR = sT0VAR; pop_TDPREDEFFECT = sTDPREDEFFECT; pop_DIFFUSIONcov = sDIFFUSIONcov; pop_asymDIFFUSION = sasymDIFFUSION; pop_asymCINT = sasymCINT; 
  }

    
} // end subject matrix creation
  

    etacov =  sT0VAR;

    } //end T0 matrices
if(verbose > 1) print ("below t0 row ", rowi);

      if(T0check>0){
        vector[nlatentpop] base;
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
      statetform[whichequals(whenvecs[2],0,0)] = parvectform(size(whichequals(whenvecs[2],0,0)),state, 2, matsetup, matvalues, si, subindices, whenvecs[2]);
    
    if(whenmat[10,2])sPARS=mcalc(sPARS,statetform,2, 10, matsetup, matvalues, si, subindices); 
if(whenmat[3,2])sDRIFT=mcalc(sDRIFT,statetform,2, 3, matsetup, matvalues, si, subindices); 
if(whenmat[7,2])sCINT=mcalc(sCINT,statetform,2, 7, matsetup, matvalues, si, subindices); 

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
    statetform[whichequals(whenvecs[2],0,0)] = parvectform(size(whichequals(whenvecs[2],0,0)),state, 2, matsetup, matvalues, si, subindices, whenvecs[2]);
    
    if(whenmat[4,2])sDIFFUSION=mcalc(sDIFFUSION,statetform,2, 4, matsetup, matvalues, si, subindices); 
if(whenmat[52,2])sJAx=mcalc(sJAx,statetform,2, 52, matsetup, matvalues, si, subindices); 
    {
  ;  
  
  } 
   
      if(statedep[4]) sDIFFUSIONcov[derrind,derrind] = sdcovsqrt2cov(sDIFFUSION[derrind,derrind],choleskymats);
      
        if(continuoustime){
          if(taylorheun==0){
            if(dtchange==1 || statedep[3]||statedep[4]||statedep[7] || 
              (T0check == 1 && (subindices[3] + subindices[4] + subindices[7]) > 0)){
                
              
              
              if(difftype==0 || (statedep[3]==0 && statedep[4]==0)){
                discreteDRIFT = matrix_exp(append_row(append_col(sDRIFT[1:nlatent, 1:nlatent],sCINT),nlplusonezerovec') * dtsmall);
                Je[savescores ? rowi : 1,,] =  matrix_exp(sJAx * dtsmall);
               
                if(statedep[3]||statedep[4]) sasymDIFFUSION[derrind,derrind] = to_matrix(  -sqkron_sumii(sJAx[derrind,derrind]) \ 
                  to_vector(sDIFFUSIONcov[derrind,derrind]), ndiffusion,ndiffusion);
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
            
          if(taylorheun==1){
            if(dtchange==1 || statedep[3]||statedep[4] || 
              (T0check == 1 && (subindices[3] + subindices[4]) > 0)){
                Kth = (add_diag(-sJAx * (dtsmall /2),1) );
                Mth = Kth \ (add_diag(sJAx * (dtsmall /2),1) );
              }
            state[1:nlatent] += Kth[1:nlatent,1:nlatent] \
              (sDRIFT[1:nlatent,1:nlatent] * state[1:nlatent] + sCINT[1:nlatent,1]) * dtsmall;
            etacov = quad_form_sym((etacov), Mth');
            etacov[derrind,derrind] += (Kth[derrind,derrind] \ sDIFFUSIONcov[derrind,derrind] / Kth[derrind,derrind]') * dtsmall;
          }
            if(intstepi >= (dt-1e-10) && savescores) Je[rowi,,] = matrix_exp(sJAx * dt); //save approximate exponentiated jacobian for smoothing
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
  
    if(T0check==0){ //nl t0
    state = sT0MEANS[,1]; //in case of t0 dependencies, may have missingness
        {
  ;  
  
  } 
  statetform[whichequals(whenvecs[1],0,0)] = parvectform(size(whichequals(whenvecs[1],0,0)),state, 1, matsetup, matvalues, si, subindices, whenvecs[1]);
    
    if(whenmat[1,1])sT0MEANS=mcalc(sT0MEANS,statetform,1, 1, matsetup, matvalues, si, subindices); 
if(whenmat[8,1])sT0VAR=mcalc(sT0VAR,statetform,1, 8, matsetup, matvalues, si, subindices); 
if(whenmat[51,1])sJ0=mcalc(sJ0,statetform,1, 51, matsetup, matvalues, si, subindices); 

    ;
      state = sT0MEANS[,1];
      if(verbose>1) print("nl t0var = ", sT0VAR, "   sJ0 = ", sJ0);
      etacov = quad_form(sT0VAR, sJ0');
    }//end nonlinear t0
    
    if(ntdpred > 0) {
      int nonzerotdpred = 0;
      for(tdi in 1:ntdpred) if(tdpreds[rowi,tdi] != 0.0) nonzerotdpred = 1;
      if(nonzerotdpred){
          {
  ;  
  
  } 
  statetform[whichequals(whenvecs[3],0,0)] = parvectform(size(whichequals(whenvecs[3],0,0)),state, 3, matsetup, matvalues, si, subindices, whenvecs[3]);
    
    if(whenmat[9,3])sTDPREDEFFECT=mcalc(sTDPREDEFFECT,statetform,3, 9, matsetup, matvalues, si, subindices); 
if(whenmat[53,3])sJtd=mcalc(sJtd,statetform,3, 53, matsetup, matvalues, si, subindices); 

      ;
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



 if (nobs_y[rowi] > 0 || savescores) {  // if some observations create right size matrices for missingness and calculate...
    
      
    int zeroint[1];
    vector[nlatentpop] basestate = state;
    zeroint[1] = 0;
    for(statei in append_array(sJyfinite,zeroint)){ //if some finite differences to do, compute these first
      state = basestate;
      if(statei>0 && (savescores + intoverstates) > 0)  state[statei] += Jstep;
          {
  ;  
  
  } 
  statetform[whichequals(whenvecs[4],0,0)] = parvectform(size(whichequals(whenvecs[4],0,0)),state, 4, matsetup, matvalues, si, subindices, whenvecs[4]);
    
    if(whenmat[10,4])sPARS=mcalc(sPARS,statetform,4, 10, matsetup, matvalues, si, subindices); 
if(whenmat[2,4])sLAMBDA=mcalc(sLAMBDA,statetform,4, 2, matsetup, matvalues, si, subindices); 
if(whenmat[6,4])sMANIFESTMEANS=mcalc(sMANIFESTMEANS,statetform,4, 6, matsetup, matvalues, si, subindices); 
if(whenmat[5,4])sMANIFESTVAR=mcalc(sMANIFESTVAR,statetform,4, 5, matsetup, matvalues, si, subindices); 
if(whenmat[54,4])sJy=mcalc(sJy,statetform,4, 54, matsetup, matvalues, si, subindices); 
    {
  ;  
  
  } 
   
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
    

      if(intoverstates==1 || savescores==1) { //classic kalman
        ycov[o,o] = quad_form(etacov, sJy[o,]'); // + sMANIFESTVAR[o,o]; shifted measurement error down
        for(wi in 1:nmanifest){ 
          if(Y[rowi,wi] != 99999 || savescores==1) ycov[wi,wi] += square(sMANIFESTVAR[wi,wi]);
          if(manifesttype[wi]==1 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += fabs((syprior[wi] - 1) .* (syprior[wi]));
          if(manifesttype[wi]==2 && (Y[rowi,wi] != 99999  || savescores==1)) ycov[wi,wi] += square(fabs((syprior[wi] - round(syprior[wi])))); 
        }
      }
        
      if(intoverstates==0) { //sampled states
        if(ncont_y[rowi] > 0) ypriorcov_sqrt[o0,o0] = sMANIFESTVAR[o0,o0]+1e-8;
      }
        
     

  {
  int skipupd = 0;
        for(vi in 1:nobs_y[rowi]){
            if(fabs(syprior[od[vi]]) > 1e10 || is_nan(syprior[od[vi]]) || is_inf(syprior[od[vi]])) {
              skipupd = 1; 
              syprior[od[vi]] =99999;
  if(verbose > 1) print("pp syprior problem! row ", rowi);
            }
          }
        if(skipupd==0){ 
          if(ncont_y[rowi] > 0){
            ypriorcov_sqrt[o0d,o0d]=cholesky_decompose(makesym(ycov[o0d, o0d],verbose,1)); 
            Ygen[ rowi, o0d] = syprior[o0d] + ypriorcov_sqrt[o0d,o0d] * Ygenbase[rowi,o0d];
          }
          if(nbinary_y[rowi] > 0) for(obsi in 1:size(o1d)) Ygen[rowi, o1d[obsi]] = (syprior[o1d[obsi]] > Ygenbase[rowi,o1d[obsi]]) ? 1 : 0; 
          for(vi in 1:nobs_y[rowi]) if(is_nan(Ygen[rowi,od[vi]])) {
            Ygen[rowi,od[vi]] = 99999;
            print("pp ygen problem! row ", rowi);
          }
        err[od] = Ygen[rowi,od] - syprior[od]; // prediction error
        }
}
      
    
      if(intoverstates==1 && size(od) > 0) {
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
        K[,od] = mdivide_right(etacov * sJy[od,]', ycov[od,od]); 
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
  
      
    }//end nobs > 0 section
  
  
  
  } // end dokalmanrows subset selection
  prevrow = rowi; //update previous row marker only after doing necessary calcs
}//end rowi
ll+=sum(llrow);


}}


}
