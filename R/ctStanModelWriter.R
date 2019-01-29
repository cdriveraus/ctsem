ctStanMatricesList <- function(ctstanmodel){
  m <- list()
  m$base <- unique(c(ctstanmodel$pars$matrix,'TDPREDEFFECT'))
  m$dynamic <- c('DRIFT','DIFFUSION','CINT')
  m$tdpred <- 'TDPREDEFFECT'
  m$measurement <- c('LAMBDA','MANIFESTMEANS','MANIFESTVAR')
  m$t0 <- c('T0MEANS','T0VAR')
  return(m)
}

ctStanCalcsList <- function(ctstanmodel){
  #extract any calcs from model and ensure spec is correct
  temp <- ctstanmodel$calcs
  calcs<-list()
  mats<-ctStanMatricesList(ctstanmodel)
  calcindices <- grep('\\]|\\[',ctstanmodel$pars$param)
  if(length(calcindices) > 0){
    for(ci in calcindices){
      temp <- c(temp,paste0(ctstanmodel$pars$matrix[ci],
        '[',ctstanmodel$pars$row[ci], ', ', ctstanmodel$pars$col[ci],'] = ',
        ctstanmodel$pars$param[ci]),';')
    }
    ctstanmodel$pars$value[calcindices] <- 99999
    ctstanmodel$pars$param[calcindices] <- NA
    if(ctstanmodel$n.TIpred > 0) ctstanmodel$pars[calcindices,c('indvarying',paste0(ctstanmodel$TIpredNames,'_effect'))] <- FALSE
  }
  
  for(mati in unique(ctstanmodel$pars$matrix)){
    temp <- gsub(mati,paste0('s',mati),temp)
  }
  
  calcs <- lapply(mats[-1], function(mlist){
    unlist(lapply(mlist, function(y) temp[grep(paste0('s',y),temp)] ))
  })
  
  ctstanmodel$calcs <- calcs
  return(ctstanmodel)
}


ctStanModelWriter <- function(ctstanmodel, gendata, extratforms){
#if arguments change make sure to change ctStanFit !
  
mats <- ctStanMatricesList(ctstanmodel)

#check when / if PARS needs to be computed
for(mlist in names(mats[-1])){
  if(any(unlist(lapply(ctstanmodel$calcs[[mlist]], function(m) grepl('sPARS',m))))) mats[[mlist]]=c(mats[[mlist]],'PARS')
}


   #intoverpop calcs setup
    intoverpopdynamiccalcs <- paste0('
    if(intoverpop==1){
',paste0(
      unlist(lapply(c(mats$dynamic), function(m){
      paste0('
    for(ri in 1:size(',m,'setup)){ //for each row of matrix setup
      if(',m,'setup[ ri,5] > 0 && ( statei == 0 || statei == nlatent + ',m,'setup[ ri,5])){ // if individually varying
        s',m,'[',m,'setup[ ri,1], ',m,'setup[ri,2]] = 
          ','tform(state[nlatent +',m,'setup[ri,5] ], ',m,'setup[ri,4], ',m,'values[ri,2], ',m,'values[ri,3], ',m,'values[ri,4], ',m,'values[ri,6] ); 
      }
    }')
    })),collapse='\n'),'}',collapse='\n')
    
    ctstanmodel$calcs$tdpred <- paste0(ctstanmodel$calcs$tdpred,';
      if(intoverpop==1){
        ',paste0(unlist(lapply('TDPREDEFFECT', function(m){
      paste0('
    for(ri in 1:size(',m,'setup)){
      if(',m,'setup[ ri,5] > 0){ 
         s',m,'[',m,'setup[ ri,1], ',m,'setup[ri,2]] = ',
          'tform(state[nlatent +',m,'setup[ri,5] ], ',m,'setup[ri,4], ',m,'values[ri,2], ',m,'values[ri,3], ',m,'values[ri,4], ',m,'values[ri,6] ); 
      }
    }')
    })),collapse='\n'),'}',collapse='\n')
    
    ctstanmodel$calcs$measurement <- paste0(ctstanmodel$calcs$measurement,';
      if(intoverpop==1){
        ',paste0(unlist(lapply(mats$measurement, function(m){
      paste0('
    for(ri in 1:size(',m,'setup)){
      if(',m,'setup[ ri,5] > 0){ 
         s',m,'[',m,'setup[ ri,1], ',m,'setup[ri,2]] = ',
         'tform(ukfstates[nlatent +',m,'setup[ri,5], statei ], ',m,'setup[ri,4], ',
            m,'values[ri,2], ',m,'values[ri,3], ',m,'values[ri,4], ',m,'values[ri,6] ); 
      }
    }')
    })),collapse='\n'),'}',collapse='\n')
    
    ctstanmodel$calcs$t0 <- paste0(ctstanmodel$calcs$t0,';
      if(intoverpop==1){
        ',paste0(unlist(lapply(mats$t0, function(m){
      paste0('
    for(ri in 1:size(',m,'setup)){
      if(',m,'setup[ ri,5] > 0){ 
         s',m,'[',m,'setup[ ri,1], ',m,'setup[ri,2]] = ',
          'tform(state[nlatent +',m,'setup[ri,5] ], ',m,'setup[ri,4], ',m,'values[ri,2], ',m,'values[ri,3], ',m,'values[ri,4], ',m,'values[ri,6] ); 
      }
    }')
    })),collapse='\n'),'}',collapse='\n')
    
    #end intoverpop setup
    


ukfilterfunc<-function(ppchecking){
  out<-paste0('
  int si;
  int counter = 0;
  vector[nlatentpop] eta; //latent states
  matrix[nlatentpop, nlatentpop] etacov; //covariance of latent states

  //measurement 
  vector[nmanifest] err;
  vector[nmanifest] ypred;
  vector[ukf ? nmanifest : 0] ystate;
  matrix[nmanifest, nmanifest] ypredcov;
  matrix[nlatentpop, nmanifest] K; // kalman gain
  matrix[nmanifest, nmanifest] ypredcov_sqrt; 

  vector[nmanifest+nmanifest+ (savescores ? nmanifest*2+nlatentpop*2 : 0)] kout[ndatapoints];

  matrix[ukf ? nlatentpop :0,ukf ? nlatentpop :0] sigpoints;

  //linear continuous time calcs
  matrix[nlatent+1,nlatent+1] discreteDRIFT;
  matrix[nlatent,nlatent] discreteDIFFUSION;

  //dynamic system matrices
  ',subjectparaminit(pop=FALSE,smats=TRUE),'

  if(nldynamics==0) discreteDIFFUSION = rep_matrix(0,nlatent,nlatent); //in case some elements remain zero due to derrind

  if(savescores) kout = rep_array(rep_vector(99999,rows(kout[1])),ndatapoints);

  for(rowi in 1:ndatapoints){

    matrix[ukf ? nlatentpop : 0, ukffull ? 2*nlatentpop +2 : nlatentpop + 2 ] ukfstates; //sampled states relevant for dynamics
    matrix[ukf ? nmanifest : 0 , ukffull ? 2*nlatentpop +2 : nlatentpop + 2] ukfmeasures; // expected measures based on sampled states
    si=subject[rowi];
    
    if(T0check[rowi] == 1) { // calculate initial matrices if this is first row for si

    ',subjectparscalc2(pop=FALSE,smats=TRUE),'

      if(ukf==1){
        eta = rep_vector(0,nlatentpop); // because some values stay zero
        sigpoints = rep_matrix(0, nlatentpop,nlatentpop);
        if(intoverpop==1) {
          if(ntipred ==0) eta[ (nlatent+1):(nlatentpop)] = rawpopmeans[indvaryingindex];
          if(ntipred >0) eta[ (nlatent+1):(nlatentpop)] = rawpopmeans[indvaryingindex] + TIPREDEFFECT[indvaryingindex] * tipreds[si]\';
        }
      }

      if(ukf==0){
        eta = sT0MEANS[,1]; //prior for initial latent state
        if(ntdpred > 0) eta += sTDPREDEFFECT * tdpreds[rowi];
        etacov =  sT0VAR;
      }

    } //end T0 matrices

    if(nldynamics==0 && ukf==0 && T0check[rowi]==0){ //linear kf time update
    
      if(continuoustime ==1){
        int dtchange = 0;
        if(si==1 && T0check[rowi -1] == 1) {
          dtchange = 1;
        } else if(T0check[rowi-1] == 1 && dT[rowi-2] != dT[rowi]){
          dtchange = 1;
        } else if(T0check[rowi-1] == 0 && dT[rowi-1] != dT[rowi]) dtchange = 1;
        
        if(dtchange==1 || (T0check[rowi-1]==1 && (si <= DRIFTsubindex[nsubjects] || si <= CINTsubindex[nsubjects]))){
          discreteDRIFT = expm2(append_row(append_col(sDRIFT,sCINT),rep_matrix(0,1,nlatent+1)) * dT[rowi]);
        }
    
        if(dtchange==1 || (T0check[rowi-1]==1 && (si <= DIFFUSIONsubindex[nsubjects]|| si <= DRIFTsubindex[nsubjects]))){
          discreteDIFFUSION[derrind, derrind] = sasymDIFFUSION[derrind, derrind] - 
            quad_form( sasymDIFFUSION[derrind, derrind], discreteDRIFT[derrind, derrind]\' );
          if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(makesym(discreteDIFFUSION,verbose,1));
        }
      }
  
      if(continuoustime==0 && T0check[rowi-1] == 1){
        discreteDRIFT=append_row(append_col(sDRIFT,sCINT),rep_matrix(0,1,nlatent+1));
        discreteDRIFT[nlatent+1,nlatent+1] = 1;
        discreteDIFFUSION=sDIFFUSION;
        if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(makesym(discreteDIFFUSION,verbose,1));
      }

      eta = (discreteDRIFT * append_row(eta,1.0))[1:nlatent];
      if(ntdpred > 0) eta += sTDPREDEFFECT * tdpreds[rowi];
      if(intoverstates==1) {
        etacov = quad_form(etacov, discreteDRIFT[1:nlatent,1:nlatent]\');
        if(ndiffusion > 0) etacov += discreteDIFFUSION;
      }
    }//end linear time update


    if(ukf==1){ //ukf time update
      vector[nlatentpop] state;
      if(T0check[rowi]==0){
        matrix[nlatentpop,nlatentpop] J;
        vector[nlatent] base;
        J = rep_matrix(0,nlatentpop,nlatentpop); //dont necessarily need to loop over tdpreds here...
        if(continuoustime==1){
          matrix[nlatentpop,nlatentpop] Je;
          matrix[nlatent*2,nlatent*2] dQi;
          for(stepi in 1:integrationsteps[rowi]){
            for(statei in 0:nlatentpop){
              if(statei>0){
                J[statei,statei] = 1e-6;
                state = eta + J[,statei];
              } else {
                state = eta;
              }
              ',paste0(ctstanmodel$calcs$dynamic,';\n',collapse=' '),';\n', intoverpopdynamiccalcs,' //remove diffusion calcs and do elsewhere
              if(statei== 0) {
                base = sDRIFT * state[1:nlatent] + sCINT[,1];
              }
              if(statei > 0) {
                J[1:nlatent,statei] = (( sDRIFT * state[1:nlatent] + sCINT[,1]) - base)/1e-6;
              }
            }
            ',paste0(ctstanmodel$calcs$dynamic,';\n',collapse=' '),' //find a way to remove this repeat
            Je= expm2(J * dTsmall[rowi]) ;
            discreteDRIFT = expm2(append_row(append_col(sDRIFT,sCINT),rep_vector(0,nlatent+1)\') * dTsmall[rowi]);
            sasymDIFFUSION = to_matrix(  -kronsum(J[1:nlatent,1:nlatent]) \\ to_vector(tcrossprod(sDIFFUSION)), nlatent,nlatent);
            discreteDIFFUSION =  sasymDIFFUSION - quad_form( sasymDIFFUSION, Je[1:nlatent,1:nlatent]\' );
            etacov = quad_form(etacov, Je\');
            etacov[1:nlatent,1:nlatent] += discreteDIFFUSION;
            eta[1:nlatent] = (discreteDRIFT * append_row(eta[1:nlatent],1.0))[1:nlatent];
          }
        }

        if(continuoustime==0){ //need covariance in here
            for(statei in 0:nlatentpop){
              if(statei>0){
                J[statei,statei] = 1e-6;
                state = eta + J[,statei];
              } else {
                state = eta;
              }
              ',paste0(ctstanmodel$calcs$dynamic,';\n',collapse=' '),';\n', intoverpopdynamiccalcs,' //remove diffusion calcs and do elsewhere
              if(statei== 0) {
                base = sDRIFT * state[1:nlatent] + sCINT[,1];
              }
              if(statei > 0) J[1:nlatent,statei] = (( sDRIFT * state[1:nlatent] + sCINT[,1]) - base)/1e-6;
            }
          ',paste0(ctstanmodel$calcs$dynamic,';\n',collapse=' '),'
          discreteDRIFT=append_row(append_col(sDRIFT,sCINT),rep_matrix(0,1,nlatent+1));
          discreteDRIFT[nlatent+1,nlatent+1] = 1;
          etacov = quad_form(etacov, J\');
          etacov[1:nlatent,1:nlatent] += sDIFFUSION;
          if(intoverstates==0) discreteDIFFUSION = cholesky_decompose(makesym(discreteDIFFUSION,verbose,1));
          eta[1:nlatent] = (discreteDRIFT * append_row(eta[1:nlatent],1.0))[1:nlatent];
        }
      } // end of non t0 time update
  
  
    if(nlmeasurement==1 || ntdpred > 0 || T0check[rowi]==1){ //ukf time update
  
      if(T0check[rowi]==1) {
        if(intoverpop==1) sigpoints[(nlatent+1):(nlatentpop), (nlatent+1):(nlatentpop)] = rawpopcovsqrt * sqrtukfadjust;
        sigpoints[1:nlatent,1:nlatent] = sT0VAR * sqrtukfadjust;
      }
      
      if(T0check[rowi]==0)  sigpoints = cholesky_decompose(makesym(etacov,verbose,1)) * sqrtukfadjust;
    
      //configure ukf states
      for(statei in 2:cols(ukfstates) ){ //for each ukf state sample
  
          state = eta; 
          if(statei > (2+nlatentpop)){
            state += -sigpoints[,statei-(2+nlatentpop)];
          } else
          if(statei > 2) state += sigpoints[,statei-2]; 

        if(T0check[rowi]==1){
          ',paste0(ctstanmodel$calcs$t0,';',collapse=' '),'
          state[1:nlatent] += sT0MEANS[,1];
        } 
        ',paste0(ctstanmodel$calcs$tdpred,';',collapse=' '),'
        if(ntdpred > 0) state[1:nlatent] +=   (sTDPREDEFFECT * tdpreds[rowi]); //tdpred effect only influences at observed time point','
        ukfstates[, statei] = state; //now contains time updated state
        if(statei==2 && ukffull==1) ukfstates[, 1] = state; //mean goes in twice for weighting
      }
  
      if(ukffull == 1) {
        eta = colMeans(ukfstates\');
        etacov = cov_of_matrix(ukfstates\') / asquared;
      }
      if(ukffull == 0){
        eta = ukfstates[,2];
        etacov = tcrossprod(ukfstates[,3:(nlatentpop+2)] - rep_matrix(ukfstates[,2],nlatentpop)) /asquared / (nlatentpop+.5);
      }
    } //end ukf if necessary time update
  } // end non linear time update


  if(savescores==1) kout[rowi,(nmanifest*4+1):(nmanifest*4+nlatentpop)] = eta;
if(verbose > 1) print("etaprior = ", eta, " etapriorcov = ",etacov);

    if(intoverstates==0 && nldynamics == 0) {
      if(T0check[rowi]==1) eta += cholesky_decompose(sT0VAR) * etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
      if(T0check[rowi]==0) eta +=  discreteDIFFUSION * etaupdbasestates[(1+(rowi-1)*nlatent):(rowi*nlatent)];
    }

    if (nobs_y[rowi] > 0) {  // if some observations create right size matrices for missingness and calculate...
    
      int o[nobs_y[rowi]]= whichobs_y[rowi,1:nobs_y[rowi]]; //which obs are not missing in this row
      int o1[nbinary_y[rowi]]= whichbinary_y[rowi,1:nbinary_y[rowi]];
      int o0[ncont_y[rowi]]= whichcont_y[rowi,1:ncont_y[rowi]];

      if(nlmeasurement==0){ //non ukf measurement
        if(intoverstates==1) { //classic kalman
          ypred[o] = sMANIFESTMEANS[o,1] + sLAMBDA[o,] * eta[1:nlatent];
          if(nbinary_y[rowi] > 0) ypred[o1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * eta[1:nlatent])));
          ypredcov[o,o] = quad_form(etacov[1:nlatent,1:nlatent], sLAMBDA[o,]\') + sMANIFESTVAR[o,o];
          for(wi in 1:nmanifest){ 
            if(manifesttype[wi]==1 && Y[rowi,wi] != 99999) ypredcov[wi,wi] += fabs((ypred[wi] - 1) .* (ypred[wi]));
            if(manifesttype[wi]==2 && Y[rowi,wi] != 99999) ypredcov[wi,wi] += square(fabs((ypred[wi] - round(ypred[wi])))); 
          }
          K[,o] = mdivide_right(etacov * append_row(sLAMBDA[o,]\',rep_matrix(0,nlatentpop-nlatent,nmanifest)[,o]), ypredcov[o,o]); 
          etacov += -K[,o] * append_col(sLAMBDA[o,],rep_matrix(0,nmanifest,nlatentpop-nlatent)[o,]) * etacov;
        }
        if(intoverstates==0) { //sampled states
          if(ncont_y[rowi] > 0) {
            ypred[o0] = sMANIFESTMEANS[o0,1] + sLAMBDA[o0,] * eta[1:nlatent];
            ypredcov_sqrt[o0,o0] = sMANIFESTVAR[o0,o0];
          }
          if(nbinary_y[rowi] > 0) ypred[o1] = to_vector(inv_logit(to_array_1d(sMANIFESTMEANS[o1,1] +sLAMBDA[o1,] * eta[1:nlatent])));
        }
      } 
  

      if(nlmeasurement==1){ //ukf measurement
        vector[nlatentpop] state; //dynamic portion of current states
        matrix[nmanifest,cols(ukfmeasures)] merrorstates;

        for(statei in 2:cols(ukfmeasures)){
          state = ukfstates[, statei];
          ',paste0(ctstanmodel$calcs$measurement,';\n',collapse=''),'
  
          ukfmeasures[, statei] = sMANIFESTMEANS[,1] + sLAMBDA * state[1:nlatent];
          if(nbinary_y[rowi] > 0) {
            ukfmeasures[o1 , statei] = to_vector(inv_logit(to_array_1d(ukfmeasures[o1 , statei])));
          }
        
        merrorstates[,statei] = diagonal(sMANIFESTVAR);
        for(wi in 1:nmanifest){ 
          if(manifesttype[wi]==1 && Y[rowi,wi] != 99999) merrorstates[wi,statei] = sMANIFESTVAR[wi,wi] + fabs((ukfmeasures[wi,statei] - 1) .* (ukfmeasures[wi,statei])); //sMANIFESTVAR[wi,wi] + (merror[wi] / cols(ukfmeasures) +1e-8);
          if(manifesttype[wi]==2 && Y[rowi,wi] != 99999) merrorstates[wi,statei] = sMANIFESTVAR[wi,wi] + square(fabs((ukfmeasures[wi,statei]- round(ukfmeasures[wi,statei])))); 
        }
          
          if(statei==2 && ukffull == 1) { 
            merrorstates[,1] = merrorstates[,2];
            ukfmeasures[ , 1] = ukfmeasures [,2];
          }
        } 
        if(ukffull == 1) {
          ypred[o] = colMeans(ukfmeasures[o,]\'); 
          ypredcov[o,o] = cov_of_matrix(ukfmeasures[o,]\') /asquared + diag_matrix(colMeans(merrorstates[o,]\')); //
          K[,o] = mdivide_right(crosscov(ukfstates\', ukfmeasures[o,]\') /asquared, ypredcov[o,o]); 
        }
        if(ukffull == 0){
          ypred[o] = ukfmeasures[o,2];
          for(ci in 3:cols(ukfmeasures)) ukfmeasures[o,ci] += -ukfmeasures[o,2];
          for(ci in 3:cols(ukfstates)) ukfstates[,ci] += -ukfstates[,2];
          ypredcov[o,o] = tcrossprod(ukfmeasures[o,3:(nlatentpop+2)]) /asquared / (nlatentpop +.5) + diag_matrix(merrorstates[o,2]);
          K[,o] = mdivide_right(ukfstates[,3:cols(ukfstates)] * ukfmeasures[o,3:cols(ukfmeasures)]\' /asquared / (nlatentpop+.5), ypredcov[o,o]); 
        }
        etacov +=  - quad_form(ypredcov[o,o],  K[,o]\');
      } //end ukf measurement


'    
,if(ppchecking) paste0('
{
int skipupd = 0;
        for(vi in 1:nobs_y[rowi]){
          if(fabs(ypred[o[vi]]) > 1e10 || is_nan(ypred[o[vi]]) || is_inf(ypred[o[vi]])) {
            skipupd = 1; 
            ypred[o[vi]] =99999;
if(verbose > 1) print("pp ypred problem! row ", rowi);
          }
        }
        if(skipupd==0){ 
          if(ncont_y[rowi] > 0) ypredcov_sqrt[o0,o0]=cholesky_decompose(makesym(ypredcov[o0, o0],verbose,1)); 
          if(ncont_y[rowi] > 0) Ygen[ rowi, o0] = ypred[o0] + ypredcov_sqrt[o0,o0] * Ygenbase[rowi,o0]; 
          if(nbinary_y[rowi] > 0) for(obsi in 1:size(o1)) Ygen[rowi, o1[obsi]] = ypred[o1[obsi]] > Ygenbase[rowi,o1[obsi]] ? 1 : 0; 
          for(vi in 1:nobs_y[rowi]) if(is_nan(Ygen[rowi,o[vi]])) {
            Ygen[rowi,o[vi]] = 99999;
print("pp ygen problem! row ", rowi);
          }
        err[o] = Ygen[rowi,o] - ypred[o]; // prediction error
        }
if(verbose > 1) {
print("rowi ",rowi, "  si ", si, 
          "  eta ",eta,"  etacov ",etacov,"  ypred ",ypred,"  ypredcov ",ypredcov, "  K ",K,
          "  sDRIFT ", sDRIFT, " sDIFFUSION ", sDIFFUSION, " sCINT ", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
          "  sT0VAR", sT0VAR, " sT0MEANS ", sT0MEANS,
          "  rawpopsd ", rawpopsd, "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        print("discreteDRIFT ",discreteDRIFT,  "  discreteDIFFUSION ", discreteDIFFUSION)
}
if(verbose > 2) print("ukfstates ", ukfstates, "  ukfmeasures ", ukfmeasures);
}
      '),
  
      if(!ppchecking) 'err[o] = Y[rowi,o] - ypred[o]; // prediction error','
    
      if(savescores==1) {
        int tmpindex[nobs_y[rowi]] = o;
        for(oi in 1:ncont_y[rowi]) tmpindex[oi] +=  nmanifest*2;
        kout[rowi,tmpindex] = err[o];
        for(oi in 1:ncont_y[rowi]) tmpindex[oi] +=  nmanifest;
        kout[rowi,tmpindex] = ypred[o];
      }
      if(intoverstates==1) eta +=  (K[,o] * err[o]);
  
      ',if(!ppchecking){
        'if(nbinary_y[rowi] > 0) kout[rowi,o1] =  Y[rowi,o1] .* (ypred[o1]) + (1-Y[rowi,o1]) .* (1-ypred[o1]); //intoverstates==0 && 
  
        if(verbose > 1) {
          print("rowi ",rowi, "  si ", si, "  eta ",eta,"  etacov ",etacov,
            "  ypred ",ypred,"  ypredcov ",ypredcov, "  K ",K,
            "  sDRIFT ", sDRIFT, " sDIFFUSION ", sDIFFUSION, " sCINT ", sCINT, "  sMANIFESTVAR ", diagonal(sMANIFESTVAR), "  sMANIFESTMEANS ", sMANIFESTMEANS, 
            "  sT0VAR", sT0VAR,  " sT0MEANS ", sT0MEANS,
            "discreteDRIFT ", discreteDRIFT, "  discreteDIFFUSION ", discreteDIFFUSION, "  sasymDIFFUSION ", sasymDIFFUSION, 
            "  rawpopsd ", rawpopsd,  "  rawpopsdbase ", rawpopsdbase, "  rawpopmeans ", rawpopmeans );
        }
        if(verbose > 2) print("ukfstates ", ukfstates, "  ukfmeasures ", ukfmeasures);
  
        if(size(o0) > 0){
          int tmpindex[ncont_y[rowi]] = o0;
          for(oi in 1:ncont_y[rowi]) tmpindex[oi] +=  nmanifest;
           if(intoverstates==1) ypredcov_sqrt[o0,o0]=cholesky_decompose(makesym(ypredcov[o0,o0],verbose,1));
           kout[rowi,o0] = mdivide_left_tri_low(ypredcov_sqrt[o0,o0], err[o0]); //transform pred errors to standard normal dist and collect
           kout[rowi,tmpindex] = log(diagonal(ypredcov_sqrt[o0,o0])); //account for transformation of scale in loglik
        }
      '},'
    }//end nobs > 0 section
  if(savescores==1) kout[rowi,(nmanifest*4+nlatentpop+1):(nmanifest*4+nlatentpop+nlatentpop)] = eta;','
}//end rowi
')}

  kalmanll <- function(){'
  if(sum(nbinary_y) > 0) {
    vector[sum(nbinary_y)] binaryll;
    counter = 1;
    for(ri in 1:ndatapoints){
      int o1[nbinary_y[ri]] = whichbinary_y[ri,1:nbinary_y[ri]]; //which indicators are observed and binary
      binaryll[counter:(counter + nbinary_y[ri]-1)] = kout[ri,o1];
      counter+= nbinary_y[ri];
    }
    ll+= sum(log(binaryll));
  }

  if(( sum(ncont_y) > 0)) {
    vector[sum(ncont_y)] errtrans[2];
    counter = 1;
    for(ri in 1:ndatapoints){
      int o0[ncont_y[ri]] = whichcont_y[ri,1:ncont_y[ri]]; //which indicators are observed and continuous
      errtrans[1,counter:(counter + ncont_y[ri]-1)] = kout[ri, o0];
      for(oi in 1:ncont_y[ri]) o0[oi] +=  nmanifest; //modify o0 index
      errtrans[2,counter:(counter + ncont_y[ri]-1)] = kout[ri, o0];
      counter+= ncont_y[ri];
    }
    ll += normal_lpdf(errtrans[1]|0,1) - sum(errtrans[2]);
  }
if(savescores) kalman = kout;
'}


subjectparaminit<- function(pop=FALSE,smats=TRUE){
  if(smats && pop) stop('smats and pop cannot both be TRUE!')
  paste0(
   paste0( unlist(lapply(mats$base,function(m){
      paste0('matrix[ ',m,'setup_rowcount ? max(',m,'setup[,1]) : 0, ',m,'setup_rowcount ? max(',m,'setup[,2]) : 0 ] ',
        ifelse(smats,'s',''),ifelse(pop,'pop_',''),m,if(!smats && !pop) paste0('[',m,'subindex[nsubjects]]'),';',collapse='\n')
    })),collapse='\n'),'

  matrix[nlatent,nlatent] ',ifelse(smats,'s',''),ifelse(pop,'pop_',''),'asymDIFFUSION',if(!smats && !pop) '[asymDIFFUSIONsubindex[nsubjects]]','; //stationary latent process variance
  vector[nt0meansstationary ? nlatent : 0] ',ifelse(smats,'s',''),ifelse(pop,'pop_',''),'asymCINT',if(!smats && !pop) '[asymCINTsubindex[nsubjects]]','; // latent process asymptotic level
  ',collapse='\n')
}


subjectparscalc2 <- function(pop=FALSE,smats=TRUE){
  out <- paste0(
 '{
  vector[nparams] rawindparams;
  vector[nparams] tipredaddition;
  vector[nparams] indvaraddition;
  int output = ',as.integer(smats==FALSE && pop==TRUE),';
  
  if(si < 2 || (si > 1 && (nindvarying >0 || ntipred > 0))){ //si of 0 used
    tipredaddition = rep_vector(0,nparams);
    indvaraddition = rep_vector(0,nparams);

    if(si > 0 && nindvarying > 0 && intoverpop==0) {
      if(fixedsubpars==0) indvaraddition[indvaryingindex] = rawpopcovsqrt * baseindparams[si];
      if(fixedsubpars==1) indvaraddition[indvaryingindex] = rawpopcovsqrt * fixedindparams[si];
    }
  
    if(si > 0 &&  ntipred > 0) tipredaddition = TIPREDEFFECT * tipreds[si]\';
  
    rawindparams = rawpopmeans + tipredaddition + indvaraddition;

  }
',

    paste0(unlist(lapply(mats$base, function(m) {
      paste0('
  if(si <= ',m,'subindex[nsubjects]){
    for(ri in 1:size(',m,'setup)){
      if(',m,'setup[ ri,5] > 0 || ',m,'setup[ ri,6] > 0 || si < 2){
        s',m,'[',m,'setup[ ri,1], ',m,'setup[ri,2]] = ',
          m,'setup[ri,3] ? ', 
            'tform(rawindparams[ ',m,'setup[ri,3] ], ',m,'setup[ri,4], ',m,'values[ri,2], ',m,'values[ri,3], ',m,'values[ri,4], ',m,'values[ri,6] ) : ',
            m,'values[ri,1]; //either transformed, scaled and offset free par, or fixed value
      }
    }
  }
      ',collapse='\n')
    })),collapse='\n'),'

  // perform any whole matrix transformations 
    

  if(si <= DIFFUSIONsubindex[nsubjects] &&(output==1|| ukf==0)) sDIFFUSION = sdcovsqrt2cov(sDIFFUSION,ukf);
  if(ukf==0 || output ==1 || nt0varstationary > 0){
    if(si <= asymDIFFUSIONsubindex[nsubjects]) {
      if(ndiffusion < nlatent) sasymDIFFUSION = to_matrix(rep_vector(0,nlatent * nlatent),nlatent,nlatent);

      if(continuoustime==1) sasymDIFFUSION[ derrind, derrind] = to_matrix( 
      -kronsum(sDRIFT[ derrind, derrind ]) \\  to_vector( 
          ( (output==0 && ukf==1) ? tcrossprod(sDIFFUSION[ derrind, derrind ]) : sDIFFUSION[ derrind, derrind ] ) + 
          IIlatent[ derrind, derrind ] * 1e-5), ndiffusion,ndiffusion);

      if(continuoustime==0) sasymDIFFUSION = to_matrix( (IIlatent2 - 
        sqkron_prod(sDRIFT, sDRIFT)) *  to_vector(
          (output==0 && ukf==1) ? tcrossprod(sDIFFUSION[ derrind, derrind ]) : sDIFFUSION[ derrind, derrind ])
        ,ndiffusion, ndiffusion);
    } //end asymdiffusion loops
  }
          
    if(nt0meansstationary > 0){
      if(si <= asymCINTsubindex[nsubjects]){
        if(continuoustime==1) sasymCINT =  -sDRIFT \\ sCINT[ ,1 ];
        if(continuoustime==0) sasymCINT =  (IIlatent - sDRIFT) \\ sCINT[,1 ];
      }
    }

          
    if(binomial==0){
      if(si <= MANIFESTVARsubindex[nsubjects]) {
         for(ri in 1:nmanifest) sMANIFESTVAR[ri,ri] = square(sMANIFESTVAR[ri,ri]);
      }
    }
          
          
    if(si <= T0VARsubindex[nsubjects]) {
      if(ukf==0 || output==1) sT0VAR = sdcovsqrt2cov(sT0VAR,ukf);
      if(nt0varstationary > 0) {
        if(ukf==1 && output ==0) sasymDIFFUSION = cholesky_decompose(makesym(sasymDIFFUSION,verbose,1));
          for(ri in 1:nt0varstationary){ // in case of negative variances
          sT0VAR[t0varstationary[ri,1],t0varstationary[ri,2] ] = 
           t0varstationary[ri,1] == t0varstationary[ri,2] ? 
            sqrt(square(sasymDIFFUSION[t0varstationary[ri,1],t0varstationary[ri,2] ])) : 
            sasymDIFFUSION[t0varstationary[ri,1],t0varstationary[ri,2] ];
        }
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
  ',if(!smats) paste0('
  if(si > 0){
',collectsubmats(pop=FALSE),'
  }',collapse=''),'
',if(pop) paste0('
  if(si == 0){
',collectsubmats(pop=TRUE),'
  }
',collapse=''),'
}
  ',collapse='\n')
  
  return(out)
}

collectsubmats <- function(pop=FALSE,matrices=c(mats$base,'asymDIFFUSION','asymCINT')){
  out<-''
  for(m in matrices){
    if(!pop) out <- paste0(out, m,'[',m,'subindex[si]] = s',m,'; \n')
    if(pop) out <- paste0(out, 'pop_',m,' = s',m,'; \n')
  }
  return(out)
}





  writemodel<-function(){
    paste0('
functions{

   matrix covsqrt2corsqrt(matrix mat, int invert){ //converts from unconstrained lower tri matrix to cor
    matrix[rows(mat),cols(mat)] o;
    vector[rows(mat)] s;
  
    for(i in 1:rows(o)){ //set upper tri to lower
      for(j in min(i+1,rows(mat)):rows(mat)){
        o[j,i] = inv_logit(mat[j,i])*2-1;  // can change cor prior here
        o[i,j] = o[j,i];
      }
      o[i,i]=1; // change to adjust prior for correlations
    }
  
    if(invert==1) o = inverse(o);
  
    for(i in 1:rows(o)){
      s[i] = inv_sqrt(o[i,] * o[,i]);
      if(is_inf(s[i])) s[i]=0;
    }
    return diag_pre_multiply(s,o);
  }

  int[] checkoffdiagzero(matrix M){
    int z[rows(M)];
    for(i in 1:rows(M)){
      z[i] = 0;
      for(j in 1:cols(M)){
        if(i!=j){
          if(M[i,j] != 0){
            z[i] = 1;
            break;
          }
        }
      }
      if(z[i]==0){ //check rows
        for(j in 1:rows(M)){
          if(i!=j){
            if(M[j,i] != 0){
              z[i] = 1;
              break;
            }
          }
        }
      }
    }
  return z;
  }

  matrix expm2(matrix M){
    matrix[rows(M),rows(M)] out;
    int z[rows(out)] = checkoffdiagzero(M);
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
    for(i in 1:size(z0)) out[z0[i],z0[i]] = exp(M[z0[i],z0[i]]);
    return out;
  }
      
  matrix sdcovsqrt2cov(matrix mat, int cholbasis){ //converts from lower partial sd and diag sd to cov or cholesky cov
    if(cholbasis==0)  {
      return(tcrossprod(diag_pre_multiply(diagonal(mat),covsqrt2corsqrt(mat, 0))));
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
  
    da = a - rep_matrix( (colMeans(a))\',rows(a));
    db = b - rep_matrix( (colMeans(b))\',rows(b));
    out = da\' * db ./ (rows(a)-1.0);
    return out;
  }

  matrix chol(matrix a, int verbose){
    matrix[rows(a),rows(a)] l;
    for(j in 1:cols(a)){
      for(i in 1:rows(a)){
        if(j==i) {
          if(j == 1){ 
            if(a[j,j] <=0) { 
              l[j,j] = 1e-5;
              if(verbose > 0) print("Negative variance ", a[j,j], " set to 1e-5 for Cholesky decomp");
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
        ',if(gendata) paste0('if(is_nan(out[rowi,coli])){
          if(verbose > 0) print("nan during makesym row ", rowi, " col ", coli);
          if(rowi==coli) out[rowi,coli] = 99999;
          if(rowi!=coli) {
            out[rowi,coli] = 0;
            out[coli,rowi] = 0;
          }
        }'),'
      }
    }
    return out;
  }

  real tform(real param, int transform, data real multiplier, data real meanscale, data real offset, data real inneroffset){
    real out;
  
    ',tformshapes(),
  if(length(extratforms) > 0) paste0(extratforms,collapse=" \n"),'

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
  int binomial; //binary data only
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

  ',paste0(unlist(lapply(c(mats$base,'asymCINT','asymDIFFUSION'),function(mati) paste0('int ',mati,'subindex[nsubjects];',collapse='\n'))),collapse='\n'),'
  ',paste0(unlist(lapply(mats$base,function(m) paste0('int ',m,'setup_rowcount;',collapse='\n'))),collapse='\n'),'
  ',paste0(unlist(lapply(mats$base,function(m) paste0('int ',m,'setup[',m,'setup_rowcount,6 ];',collapse='\n'))),collapse='\n'),'
  ',paste0(unlist(lapply(mats$base,function(m) paste0('matrix[',m,'setup_rowcount, 6] ',m,'values;'))),collapse='\n'),'
  int TIPREDEFFECTsetup[nparams, ntipred];
  int nmatrixslots;
  int popsetup[nmatrixslots,6];
  real popvalues[nmatrixslots,6];
  int savescores;
  int nlatentpop;
  int fixedhyper;
  int fixedsubpars;
  vector[fixedhyper ? nindvarying : 0] fixedrawpopmeans;
  matrix[fixedhyper ? nindvarying : 0, fixedhyper ? nindvarying : 0] fixedrawpopchol;
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
  vector[fixedhyper ? nparams-nindvarying : nparams] rawpopmeansbase; // population level means \n','
  vector',if(!is.na(ctstanmodel$rawpopsdbaselowerbound)) paste0('<lower=',ctstanmodel$rawpopsdbaselowerbound[1],'>'),'[(nindvarying * (1-fixedhyper))? nindvarying : 0] rawpopsdbase; //population level std dev
  //vector[(nindvarying * (1-fixedhyper)) ? nindvarying : 0] rawpopsdpropsbase; //[(nindvarying * (1-fixedhyper)) ? 1 : 0];
  vector[fixedhyper ? 0 : nindvaryingoffdiagonals] sqrtpcov; // unconstrained basis of correlation parameters
  vector[fixedsubpars ? 0 : (intoverpop ? 0 : nindvarying)] baseindparams[fixedsubpars ? 0 : (intoverpop ? 0 : nsubjects)]; //vector of subject level deviations, on the raw scale
  
  vector[ntipredeffects] tipredeffectparams; // effects of time independent covariates
  vector[nmissingtipreds] tipredsimputed;
  
  vector[intoverstates ? 0 : nlatent*ndatapoints] etaupdbasestates; //sampled latent states posterior
}
      
transformed parameters{
  vector[nparams] rawpopmeans; // population level means \n','
  vector[fixedhyper ? 0 : nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying,nindvarying] rawpopcovsqrt; 
  real ll;
  vector[nmanifest+nmanifest+ (savescores ? nmanifest*2+nlatentpop*2 : 0)] kalman[savescores ? ndatapoints : 0];

  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipreds; //tipred values to fill from data and, when needed, imputation vector
  matrix[nparams, ntipred] TIPREDEFFECT; //design matrix of individual time independent predictor effects

  if(fixedhyper==1) {
    rawpopmeans[indvaryingindex] = fixedrawpopmeans;
    rawpopmeans[notindvaryingindex]=rawpopmeansbase;
}
  if(fixedhyper==0) rawpopmeans = rawpopmeansbase;

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

  if(nindvarying > 0 &&fixedhyper ==0){
    int counter =0;
    rawpopsd = ',ctstanmodel$rawpopsdtransform, '; // sqrts of proportions of total variance
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
      covsqrt2corsqrt(rawpopcovsqrt,0))),verbose,1)); 
  }//end indvarying par setup

  if(nindvarying > 0 &&fixedhyper ==1) rawpopcovsqrt = fixedrawpopchol; 

  ll=0;
if(dokalman==1){',
if(!gendata) ukfilterfunc(ppchecking=FALSE),
if(!gendata) kalmanll(),'
}

}
      
model{

  if(fixedhyper==1) baseindparams ~ multi_normal_cholesky(rep_vector(0,nindvarying),IIindvar);
  if(intoverpop==0 && fixedhyper==0 && fixedsubpars == 1) fixedindparams ~ multi_normal_cholesky(rep_vector(0,nindvarying),IIindvar);

  if(nopriors==0){
   target += normal_lpdf(rawpopmeansbase|0,1);
  
    if(ntipred > 0){ 
      tipredeffectparams ~ normal(0,tipredeffectscale);
      tipredsimputed ~ normal(0,tipredsimputedscale);
    }
    
    if(nindvarying > 0){
      if(nindvarying >1 && fixedhyper==0) {
        sqrtpcov ~ normal(0,1);
      }
      if(intoverpop==0 && fixedhyper==0 && fixedsubpars == 0) baseindparams ~ multi_normal_cholesky(rep_vector(0,nindvarying),IIindvar);
      if(fixedhyper==0) rawpopsdbase ~ ',ctstanmodel$rawpopsdbase,';
    }
  } //end pop priors section
  
  if(intoverstates==0) etaupdbasestates ~ normal(0,1);
  
  target += ll;
  if(verbose > 0) print("lp = ", target());
}
generated quantities{
  vector[nparams] popmeans;
  vector[nparams] popsd;
  matrix[nindvarying,nindvarying] rawpopcov;
  matrix[nindvarying,nindvarying] rawpopcorr;
  matrix[nparams,ntipred] linearTIPREDEFFECT;
  ',if(gendata) 'vector[nmanifest] Ygen[ndatapoints];','

  ',subjectparaminit(pop=TRUE,smats=FALSE),'

  ',subjectparaminit(pop=FALSE,smats=FALSE),'

  {
  ',subjectparaminit(pop=FALSE,smats=TRUE),'

  for(si in 0:subject[ndatapoints]){
    ',subjectparscalc2(pop=TRUE,smats=FALSE),'
    }
  }

rawpopcov = tcrossprod(rawpopcovsqrt);
rawpopcorr = quad_form_diag(rawpopcov,inv_sqrt(diagonal(rawpopcov)));

popsd = rep_vector(0,nparams);
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

',if(gendata) paste0('
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
',ukfilterfunc(ppchecking=TRUE),'

}}
',collapse=";\n"),'

}
')
  }
m <- writemodel()
return(m)
}
