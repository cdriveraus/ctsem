#' apply discrete latent autoregressive moving average SEM model
larmaModel<-function(datawide,Tpoints,n.manifest,n.latent,LAMBDA,n.TDpred=0,
                     n.TIpred=0,merror=FALSE,merrorcov=FALSE,traits=FALSE,
                     reasonable=FALSE,causalpredictors=TRUE,QAR=FALSE,
                     Qterm=TRUE,randomeffects=FALSE,tconstrainedAR=TRUE,
                     tconstrainedQ=TRUE,AR=TRUE,ARx=FALSE,lagxorder=2,MA=FALSE,
                     MAcross=FALSE,detrend=FALSE,inits=NULL,
                     withinconstraints=FALSE,freemanifestmeans=FALSE,
                     betweenconstraints=FALSE,continuous=FALSE,omxExp=TRUE,
                     CINT=TRUE){
  require(OpenMx)
  
  colnames(datawide)    <- 
    c(paste0("V",1:(Tpoints * n.manifest)), paste0("I",1:(Tpoints - 1))) #manifests
#   colnames(datawide)[(Tpoints * n.manifest + (Tpoints-1)*n.TDpred + 1):(Tpoints * n.manifest + (Tpoints-1)*n.TDpred + Tpoints-1)]    <- 
#      # intervals  
  
  if(n.TDpred > 0){ 
    for(i in 1:n.TDpred){
      colnames(datawide)[(Tpoints * n.manifest + (i - 1) * (Tpoints - 1) +1): (Tpoints * n.manifest  + i * (Tpoints - 1))]    <- 
        c(paste0("P",rep(i,each = (Tpoints - 1)),"_",1:(Tpoints - 1)))
    }
  }
  
  if(n.TIpred > 0) {
    colnames(datawide)[(Tpoints * n.manifest+ (Tpoints - 1) * n.TDpred  + (Tpoints - 1) + 1):
        (Tpoints * n.manifest+(Tpoints - 1) * n.TDpred + (Tpoints - 1) +n.TIpred)]    <- 
      c(paste0("Z",1:(n.TIpred)))
  }
  
  
  manifests<-c()
  predictors<-c()
  latents<-c()  
  
  if(nrow(datawide)==1) withinconstraints<-TRUE #can't estimate initial variance
  
  if(n.manifest>n.latent && all(freemanifestmeans==FALSE)) {
    print("Warning:  Mean structure may be over-constrained - consider freeing
          some manifest means")
  }
  
  for(i in 1:n.manifest){ #for every manifest process create indicator and residual variables
    assign(paste0("manifestV",i),paste0("V",seq(i,Tpoints*n.manifest,n.manifest))) #indicator
    if(merror==TRUE){ #then create manifest residuals
      assign(paste0("manifestV",i,"e"),paste0("V",seq(i,Tpoints*n.manifest,n.manifest),"e")) 
      latents<-c(latents,get(paste0("manifestV",i,"e"))) #also to the latent list
    }
  }
  manifests<-paste0("V",1:(Tpoints*n.manifest)) #create manifest list
  manifests<-sort(manifests)
  
  
  
  Tinit<-2
  if(MA==TRUE) Tinit<-1
  for(i in 1:n.latent){
    assign(paste0("latentV",i),paste0("L",i,"_",seq(1,Tpoints,1)))
    latents<-c(latents,get(paste0("latentV",i)))
    #     assign(paste0("latentV",i,"d"),paste0("L",i,"_",seq(2,Tpoints,1),"d"))
    #     assign(paste0("Q",i),paste0("Q",i,"_",seq(1,Tpoints,1)))
    #     assign(paste0("latentV",i,"dt"),paste0("L",i,"_",seq(1,Tpoints,1),"dt"))
    if(CINT==TRUE) latents<-c(latents,paste0("i",i),paste0("ci",i))
    if(Qterm==TRUE) latents<-c(latents,paste0("Q",i,"_",seq(1,Tpoints,1)))
    if(traits==TRUE) latents<-c(latents,paste0("trait",i),paste0("traitint",i))
    if(QAR==TRUE) latents<-c(latents,paste0("QARi",i))
  }
  latents<-sort(latents)
  
  
  if(n.TDpred>0){
    for(i in 1:n.TDpred){
      assign(paste0("TDpredProcess",i),paste0("P",i,"_",seq(1,Tpoints-1,1)))
      manifests<-c(manifests,get(paste0("TDpredProcess",i)))
      predictors<-c(predictors,get(paste0("TDpredProcess",i)))
    }
  }
  
  if(n.TIpred>0){
    for(i in 1:n.TIpred){
      manifests<-c(manifests,paste0("Z",i))
      predictors<-c(predictors,paste0("Z",i))
    }
  }
  
  processInputMatrix<-function(x,symmetric=FALSE,diagadd=0,addCharacters=NULL,...){
    inputm<-eval(parse(text=x))
    inputmFixed<-suppressWarnings(matrix(as.numeric(inputm),nrow=nrow(inputm),ncol=ncol(inputm))) #calc fixed only matrix with f's appended
    #inputmFixed[!is.na(inputmFixed)]<-matrix(paste0(addCharacters,inputm[!is.na(inputmFixed)]),nrow=nrow(inputm))
    inputmFixed[!is.na(inputmFixed)]<-paste0(addCharacters,inputm[!is.na(inputmFixed)])
    inputmCharacter<-inputm
    inputmCharacter[!is.na(inputmFixed)]<-NA #calc character only matrix
    assign(paste0(x,"fixed"),inputmFixed,pos=parent.frame(1)) #output numeric only matrix (kept to determine fixed params)
    assign(paste0(x,"character"),inputmCharacter,pos=parent.frame(1)) #output character only matrix
    
    #starting values
    #tempx<-matrix(runif(nrow(inputm)*ncol(inputm),0.01,.05),nrow=nrow(inputm),ncol=ncol(inputm)) #generate random matrix
    tempx<-matrix(.2,nrow=nrow(inputm),ncol=ncol(inputm)) #generate fixed low correlation matrix
    if(diagadd!=0)  {tempx<-tempx+chardiag(ncol(inputm))*diagadd} #increase diagonals if necessary
    if(symmetric==TRUE)  {tempx[row(tempx)>col(tempx)]<-t(tempx)[col(tempx)<row(tempx)]} #set symmetric if necessary
    tempx[!is.na(inputmFixed)]<-inputmFixed[!is.na(inputmFixed)] #input fixed values
    
    if(!is.null(inits)){ #if there are some inits specified,check each part of x
      for(i in 1:length(inputmCharacter)){ #for every cell in x
        if(!is.na(inputmCharacter[i])){ #if it is not NA
          if(any(inputmCharacter[i]==inits[,1])){  #and any cells match column 1 of labelinits
            print(paste0("Processing inits for ",x)) 
            tempx[i]<-as.numeric(inits[which(inputmCharacter[i]==inits[,1],arr.ind=TRUE)[1],2]) #for labels with init specified,insert init into value matrix
            #eval(parse(text=paste0(x,"inits[i]<<-inits[which(inits==get(x)[i]),1]"))) #set xinits[i] to the labelname
            #eval(parse(text=paste0(x,"[i]<<-as.numeric(inits[which(inits==get(x)[i]),2])")))  #and x[i] to the init value for the label
          }}}}
    assign(x,tempx,pos=parent.frame(1))
  }
  
  processInputMatrix("LAMBDA",rows=n.manifest,cols=n.latent,symmetric=FALSE,diagadd=0,addCharacters="f")
  
  LAMBDALabels<-matrix(paste0("lambda",1:n.manifest,rep(1:n.latent,each=n.manifest)),nrow=n.manifest)
  if(any(!is.na(LAMBDAcharacter))) LAMBDALabels[!is.na(LAMBDAcharacter)]<-LAMBDAcharacter[!is.na(LAMBDAcharacter)]
  
  #remove F's from value matrices and coerce to numeric
  removeF<-function(matrices){
    for(i in 1:length(matrices)){
      x<-eval(parse(text=matrices[i]))
      assign(matrices[i],matrix(as.numeric(gsub("f","",x)),nrow=nrow(x)),pos=sys.frame(-1))
    }
  }
  
  removeF(c("LAMBDA"))
  
  
  
  larmaModel <-OpenMx::mxModel("larmaModel",type="RAM",
                       mxData(datawide,type="raw"),
                       manifestVars=manifests,
                       latentVars=latents)
  
  
  #-------------------MEASUREMENT ERROR----------------------#  
  if(merror==TRUE){ #fix if varying number of manifests
    if(reasonable==TRUE) merrorLbound<-0
    if(reasonable==FALSE) merrorLbound<-NA
    
    larmaModel<-OpenMx::mxModel(larmaModel,
                        
                        #measurement error loading
                        mxPath(from=paste0("V",1:(n.manifest*Tpoints),"e"),
                               to=paste0("V",1:(n.manifest*Tpoints)),
                               arrows=1,free=FALSE,values=1,connect="single"),
                        
                        #measurement error variance
                        mxPath(from=paste0("V",1:(n.manifest*Tpoints),"e"),
                               arrows=2,free=TRUE,values=.1,lbound=merrorLbound,
                               labels=paste0("theta",1:n.manifest))      
    )
    
    if(any(merrorcov!=FALSE)){ #measurement error structure over time
      for(t in 1:(Tpoints-1)){
        larmaModel<-OpenMx::mxModel(larmaModel,
                            mxPath(from=paste0("V",1:(n.manifest*(Tpoints-1))),
                                   to=paste0("V",(n.manifest+1):(n.manifest*Tpoints)),
                                   arrows=2,free=merrorcov,values=0.05,lbound=merrorLbound,
                                   labels=paste0("thetacov",1:n.manifest))
        )
      }#end tpoints loop
    }#end merrorcov    
  } #end merror
  
  
  
  
  
  
  ########### MANIFEST MEANS
  if(any(freemanifestmeans!=FALSE)){
    
    larmaModel<-OpenMx::mxModel(larmaModel,
                        mxPath(from="one",
                               to=paste0("V",1:(n.manifest*Tpoints)),
                               arrows=1,free=freemanifestmeans,values=0,
                               labels=paste0("manint",1:n.manifest))
    )
  }
  
  
  #latent to manifest loadings
  for(i in 1:n.latent){
    for(j in 1:Tpoints){
      larmaModel<-OpenMx::mxModel(larmaModel,         
                          
                          mxPath(paste0("L",i,"_",j),
                                 to=paste0("V",((j-1)*n.manifest+1):((j-1)*n.manifest+n.manifest)),
                                 connect="all.bivariate",arrows=1,free=!is.na(LAMBDAcharacter[,i]),values=LAMBDA[,i],
                                 labels=LAMBDALabels[,i])
      )
    }
  }
  
  
  
  
  
  
  #-------------------CINT----------------------# 
  if(CINT==TRUE){
    if(continuous==FALSE) {
      CINTloadings<-NA
      CINTfree<-T
      CINTlabels<-paste0("cint",1:n.latent)
    }
    
    if(continuous==TRUE) { ######Continuous CINT setup
      CINTloadings<-paste0("intd",1:(Tpoints-1),"[",rep(1:n.latent,each=(Tpoints-1)),",1]")
      CINTfree<-F
      CINTlabels<-NA
      
      INTalgs<-list()
      for(i in 1:(Tpoints-1)){
        algName<-paste0("intd",i)
        defcall<-paste0("datawide.I",i)
        fullAlgString<-paste0("solve(DRIFT)%*%((EVEC %*% (exp(",defcall," %x% EVAL) - tempa)%*% solve(EVEC))-II)%*%t(CINT)")
        if(omxExp==TRUE) fullAlgString<-paste0("solve(DRIFT)%*%(expm(DRIFT %x%",defcall,")-II)%*%t(CINT)")
        INTalgs[i]<-eval(substitute(mxAlgebra(theExpression,name = algName),
                                    list(theExpression = parse(text=fullAlgString)[[1]])))
        
        larmaModel<-OpenMx::mxModel(larmaModel,
                            INTalgs,
                            mxMatrix(type="Iden",nrow=n.latent,ncol=n.latent,name="II"),
                            mxMatrix(type="Full",labels=paste0("cint",1:n.latent),values=.5,free=TRUE,nrow=1,ncol=n.latent,name="CINT")
        )
      }
    }#####END continuous CINT setup
    
    for(i in 1:n.latent){
      larmaModel<-OpenMx::mxModel(larmaModel,
                          
                          mxPath(to=paste0("L",i,"_",2:(Tpoints)),#continuous intercept loadings
                                 from=paste0("ci",i),
                                 values=1,free=FALSE,arrows=1,connect="single",labels=CINTloadings[((i-1)*(Tpoints-1)+1):(i*(Tpoints-1))]),
                          
                          # continuous intercept mean
                          mxPath(from="one",to=paste0("ci",i),arrows=1,free=CINTfree,values=1,labels=CINTlabels[i]),
                          
                          mxPath(to=paste0("L",i,"_",1),#initial intercept loadings
                                 from=paste0("i",i),
                                 values=1,free=FALSE,arrows=1,connect="single")
      )
      #initial intercept means
      if(withinconstraints==FALSE) {
        int.labels<-paste0("m1_",i)
      }
      if(withinconstraints==TRUE) {
#         int.labels<-paste0("m1Alg[",i,",]")
        int.labels<-paste0("m1_",i)
#         larmaModel<-OpenMx::mxModel(larmaModel,
#                             mxAlgebra(CINT%x%DRIFT,name="m1Alg")
      }
        larmaModel<-OpenMx::mxModel(larmaModel,     
                            mxPath(from="one",
                                   to=paste0("i",i),
                                   arrows=1,free=TRUE,values=.05,
                                   labels=paste0("m1_",i))
        )
      
    }
  }#end CINT and latent loop
  
  
  
  
  #####INNOVATIONS
  if(Qterm==TRUE){
    
    if(continuous==FALSE) {
      Qlabels<-rep(indexMatrix(dimension=n.latent,starttext=paste0("Q["),sep=",",
                               endtext="]",unique=FALSE,symmetrical=TRUE),times=(Tpoints-1))
    }
    
    if(continuous==TRUE){ #####Continuous Q algebras and labels
      Qlabels<-paste0("Qd",rep(1:(Tpoints-1),each=n.latent^2),"[",1:(n.latent^2),",1]")
      
      Qdalgs<-list()
      for(i in 1:(Tpoints-1)){
        defcall<-paste0("datawide.I",i)
        algName<-paste0("Qd",i)
        fullAlgString<-paste0("solve(DRIFTHATCH)%*%((EVECH%*%(exp(",defcall,"%x%EVALH)-tempb)%*%solve(EVECH))-(II%x%II))%*%rvectorize(Q)")
        if(omxExp==TRUE) fullAlgString<-paste0("solve(DRIFTHATCH)%*%((expm(DRIFTHATCH %x% ",defcall,"))-(II%x%II))%*%rvectorize(Q)") 
        Qdalgs[i]<-eval(substitute(mxAlgebra(theExpression,name = algName),list(theExpression = parse(text=fullAlgString)[[1]])))
      }
      larmaModel<-OpenMx::mxModel(larmaModel,
                          Qdalgs,
                          mxAlgebra(DRIFT%x%II + II%x%DRIFT,name = "DRIFTHATCH"),
                          mxMatrix("Full",values=(matrix(1,n.latent**2,n.latent**2)-diag(n.latent**2)),name = "tempb")
      )
    }### end continuous Q section    
    
    
    for(t in Tinit:Tpoints){  # Q add latent disturbances to latents
      for(i in 1:n.latent){
        larmaModel<-OpenMx::mxModel(larmaModel,     
                            mxPath(to=paste0("L",i,"_",t),from=paste0("Q",i,"_",t),
                                   connect="single",free=FALSE,arrows=1,values=1)
        )
      }
    } #end tpoint loop and Q loadings
    
    #disturbance variance and covariance  
    if(reasonable==TRUE) {
      Qbounds<-matrix(NA,nrow=n.latent,ncol=n.latent)
      Qbounds[row(Qbounds)==col(Qbounds)]<-0
    }
    if(reasonable==FALSE) Qbounds<-NA
    
    for(i in 2:Tpoints){
      larmaModel<-OpenMx::mxModel(larmaModel,          
                          mxPath(from=paste0("Q",1:n.latent,"_",i),
                                 arrows=2,free=FALSE,connect="unique.pairs",
                                 values=(diag(rnorm(n.latent,.5,.1))+.1)[lower.tri(diag(n.latent),diag=TRUE)],       
                                 labels=Qlabels [((i-2)*n.latent^2+1):((i-1)*n.latent^2)] [lower.tri(diag(n.latent),diag=TRUE)] )        
      )
    }
    #Q matrix
    larmaModel<-OpenMx::mxModel(larmaModel, 
                        mxMatrix(values=diag(.5,n.latent)+.1,free=TRUE,nrow=n.latent,ncol=n.latent,name="Q",
                                 labels=c(indexMatrix(dimension=n.latent,symmetric=TRUE,starttext=paste0("Q"))),
                                 lbound=c(Qbounds))
    )
    
  }#end if Qterm==T
  
  
  
  
  
  
  ###########initial covariance of processes if multiple individuals
    PHIbounds<-NA
    if(reasonable==TRUE) {
      PHIbounds<-matrix(NA,nrow=n.latent,ncol=n.latent)
      PHIbounds[row(PHIbounds)==col(PHIbounds)]<-0
    }
    
    if(withinconstraints==TRUE){
      
      if(continuous==TRUE){ #add discrete matrices
        larmaModel<-OpenMx::mxModel(larmaModel,
                            mxAlgebra(name="dQalg",
                                      solve(DRIFTHATCH)%*%((expm(DRIFTHATCH))-(II%x%II))%*%rvectorize(Q)),
                            mxMatrix(name="dQ",values=NA,labels=paste0("dQalg[",1:(n.latent^2),",1]"),nrow=n.latent,ncol=n.latent),
                            mxAlgebra(name="discreteDRIFT",expm(DRIFT))
        )}
      
      philabels<-paste0("withinphi[",matrix(1:(n.latent^2),nrow=n.latent)[lower.tri(diag(n.latent),diag=TRUE)],",1]")
      
      larmaModel<-OpenMx::mxModel(larmaModel,
                          mxMatrix(name="bigI",values=diag(n.latent^2),nrow=n.latent^2,ncol=n.latent^2,free=FALSE),
                          mxAlgebra(solve(bigI-DRIFT%x%DRIFT)%*%cvectorize(Q),name="withinphi")
      )     
    }
    
    if(withinconstraints==FALSE) {
      philabels<-indexMatrix(dimension=n.latent,starttext="withinphi",unique=TRUE,symmetric=TRUE)
    }
    
    #asymptotic within person variance
    if(Qterm==TRUE){
      larmaModel<-OpenMx::mxModel(larmaModel,
                          mxPath(from=paste0("Q",1:n.latent,"_1"),
                                 to=paste0("L",1:n.latent,"_1"),
                                 values=1,free=FALSE,connect="single"),
                          
                          mxPath(from=paste0("Q",1:n.latent,"_1"),
                                 to=paste0("Q",1:n.latent,"_1"),
                                 connect="unique.pairs",arrows=2,
                                 free=ifelse(withinconstraints==TRUE,F,T),
                                 lbound=PHIbounds[lower.tri(diag(n.latent),
                                                            diag=TRUE)],
                                 values=(diag(rnorm(n.latent,.5,.1),n.latent)+.1)[which(lower.tri(diag(n.latent),diag=TRUE))],
                                 labels=philabels)
      )
      
      if(traits==TRUE && withinconstraints==FALSE) #remove asymptotic individual variance
        larmaModel<-OpenMx::mxModel(larmaModel,
                            mxPath(from=paste0("Q",1:n.latent,"_1"),
                                   to=paste0("Q",1:n.latent,"_1"),
                                   connect="unique.pairs",arrows=2,
                                   free=FALSE,
                                   values=0)
        )}#end if Qterm==T
  
  
  if(ARx==TRUE){
    for(i in 2:lagxorder){
      larmaModel<-OpenMx::mxModel(larmaModel,          
        mxPath(from=paste0("Q",1:n.latent,"_",i),
          arrows=2,free=TRUE,connect="unique.pairs")         
      )
    }
      
      for(i in 1:n.latent){
      OpenMx::mxModel(larmaModel,        
        mxPath(to=paste0("L",i,"_",2:(Tpoints)),#continuous intercept loadings
          from=paste0("ci",i),
          values=1,free=TRUE,arrows=1,connect="single")
      )
    }
}
    
    
  
  
  
  
  ####autoregressive structure
  if(AR==TRUE){
    if(continuous==FALSE) {
      if(reasonable==TRUE) {
        ARlbounds<-matrix(-1,nrow=n.latent,ncol=n.latent)
        ARlbounds[row(ARlbounds)==col(ARlbounds)]<-0      
        ARubounds<-1
      }
      DRIFTlabels<-rep(c(indexMatrix(dimension=n.latent,starttext="DRIFT[",sep=",",endtext="]")),times=(Tpoints-1))
      DRIFTinits<-c((diag(.5,n.latent)+.1))
    }
    
    if(continuous==TRUE){ ########Continuous DRIFT labels and algebras
      if(reasonable==TRUE) {
        ARlbounds<-matrix(-20,nrow=n.latent,ncol=n.latent)
        ARlbounds[row(ARlbounds)==col(ARlbounds)]<--2.6      
        ARubounds<-matrix(3,nrow=n.latent,ncol=n.latent)
        ARubounds[row(ARubounds)==col(ARubounds)]<-0
      }
      DRIFTlabels<-paste0("EXPd",rep(1:(Tpoints-1),each=n.latent^2),
                          c(indexMatrix(dimension=n.latent,starttext="[",sep=",",endtext="]")))
      
      DRIFTinits<-c((diag(-.5,n.latent)+.1))
      
      EXPalgs<-list()
      for(i in 1:(Tpoints-1)){
        defcall<-paste0("datawide.I",i)
        algName<-paste0("EXPd",i)
        fullAlgString<-paste0("EVEC %*% (exp(",defcall," %x% EVAL) - tempa) %*% solve(EVEC)")
        if(omxExp==TRUE) fullAlgString<-paste0("expm(DRIFT %x%",defcall,")")
        EXPalgs[i]<-eval(substitute(mxAlgebra(theExpression,name = algName),
                                    list(theExpression = parse(text=fullAlgString)[[1]])))
      }
      larmaModel<-OpenMx::mxModel(larmaModel,
                          EXPalgs,
                          mxAlgebra(eigenvec(DRIFT),name="EVEC"),
                          mxAlgebra(vec2diag(eigenval(DRIFT)),name="EVAL"),
                          mxAlgebra(vec2diag(eigenval(DRIFTHATCH)),name = "EVALH"),
                          mxAlgebra(eigenvec(DRIFTHATCH),name = "EVECH")
      )
    }###End continuous DRIFT section
    
    
    
    if(reasonable==FALSE) {
      ARubounds<-NA
      ARlbounds<-NA
    }
    for(j in 1:(Tpoints-1)){
      
      larmaModel<-OpenMx::mxModel(larmaModel,   
                          
                          mxPath(from=paste0("L",1:n.latent,"_",j),
                                 to=paste0("L",1:n.latent,"_",j+1),
                                 arrows=1,free=FALSE,
                                 values=c((diag(.5,n.latent)+.1)),
                                 connect="all.bivariate",
                                 labels=DRIFTlabels [((j-1)*n.latent^2+1):((j)*n.latent^2)])
      )
    }
    #AR matrix
    larmaModel<-OpenMx::mxModel(larmaModel,
                        mxMatrix(values=DRIFTinits,lbound=ARlbounds,ubound=ARubounds,
                                 free=TRUE,nrow=n.latent,ncol=n.latent,name="DRIFT",
                                 labels=c(indexMatrix(dimension=n.latent,starttext="DRIFT")))
    )
  }
  
  if(ARx==TRUE){
    for(j in 1:(Tpoints-lagxorder)){
      
      larmaModel<-OpenMx::mxModel(larmaModel,   
                          
                          mxPath(from=paste0("L",1:n.latent,"_",j),
                                 to=paste0("L",1:n.latent,"_",j+lagxorder),
                                 arrows=1,free=TRUE,
                                 values=c((diag(.3,n.latent)+.05)),
                                 connect="all.bivariate",
                                 labels=paste0("lag",lagxorder,"_A",which(diag(n.latent)==diag(n.latent),arr.ind=TRUE)[,2,drop=F],
                                               which(diag(n.latent)==diag(n.latent),arr.ind=TRUE)[,1,drop=F]))
      )
      
      
      
    }
  }
  
  
  
  
  
  
  if(traits==TRUE){ #########TRAIT VARIANCE  
    
    
    for(t in 2:Tpoints){ #trait latent loadings
      larmaModel<-OpenMx::mxModel(larmaModel,
                          mxPath(from=paste0("trait",1:n.latent),
                                 to=paste0("L",1:n.latent,"_",t),
                                 connect="single",values=1,arrows=1,free=FALSE)
      )}
    
    if(reasonable==TRUE) { #set lower bounds
      PHITRAITbounds<-matrix(NA,nrow=n.latent,ncol=n.latent)
      PHITRAITbounds[row(PHITRAITbounds)==col(PHITRAITbounds)]<-0
    }    
    if(reasonable==FALSE) PHITRAITbounds<-NA
    
    if(betweenconstraints==TRUE){ #algebra constraint for between person overall variance
      betweenphilabels<-paste0("betweenphi[",matrix(1:(n.latent^2),nrow=n.latent)[lower.tri(diag(n.latent),diag=TRUE)],",1]") 
      #       betweenphilabels<-c(indexMatrix(starttext="betweenphi[",sep=",",endtext="]",dimension=n.latent,symmetric=FALSE,unique=FALSE))
      larmaModel<-OpenMx::mxModel(larmaModel,
                          
                          mxMatrix(name="smallI",values=diag(n.latent),nrow=n.latent,ncol=n.latent,free=FALSE),
                          mxAlgebra(solve((smallI-DRIFT)%x%(smallI-DRIFT))%*%cvectorize(PHITRAIT),name="betweenphi"),
                          mxMatrix(name="betweenphimatrix",values=NA,free=FALSE,nrow=n.latent,ncol=n.latent,
                                   labels=paste0("betweenphi[",1:(n.latent^2),",1]")),
                          
                          mxAlgebra((smallI-DRIFT)%*%betweenphimatrix,name="traitT1cov"),
                          
                          mxPath(from=paste0("traitint",1:n.latent),#initial trait variance
                                 to=paste0("L",1:n.latent,"_",1),
                                 connect="single",values=1,arrows=1,free=FALSE),#,labels=betweenphilabels)        
                          
                          # trait intercept and trait covariance (PHIT1TRAIT in ctsem)                
                          mxPath(from=paste0("traitint",1:n.latent),
                                 to=paste0("trait",1:n.latent),
                                 arrows=2,free=FALSE,values=.1,connect="unique.bivariate",
                                 labels=c(indexMatrix(starttext="traitT1cov[",sep=",",endtext="]",
                                                      dimension=n.latent,symmetric=FALSE,unique=FALSE))),
                          
                          mxPath(from=paste0("traitint",1:n.latent),#trait intercept variance
                                 arrows=2,free=FALSE,
                                 values=c(diag(.2,n.latent)+.05)[lower.tri(diag(n.latent),diag=TRUE)],
                                 connect="unique.pairs",labels=betweenphilabels)
      )
    }
    
    if(betweenconstraints==FALSE) {
      betweenphilabels<-c(indexMatrix(starttext="betweenphi",dimension=n.latent,symmetric=TRUE,unique=TRUE))
      
      larmaModel<-OpenMx::mxModel(larmaModel,#traitint latent loadings
                          mxPath(from=paste0("traitint",1:n.latent),
                                 to=paste0("L",1:n.latent,"_1"),
                                 connect="single",values=1,arrows=1,free=FALSE),
                          
                          ##algebras for testing purposes
                          mxMatrix(name="smallI",values=diag(n.latent),nrow=n.latent,ncol=n.latent,free=FALSE),
                          mxAlgebra(solve((smallI-DRIFT)%x%(smallI-DRIFT))%*%cvectorize(PHITRAIT),name="betweenphi"),
                          
                          mxMatrix(name="betweenphimatrix",values=NA,free=FALSE,nrow=n.latent,ncol=n.latent,
                                   labels=paste0("betweenphi[",1:(n.latent^2),",1]")),
                          
                          mxAlgebra((smallI-DRIFT)%*%betweenphimatrix,name="traitT1cov"),
                          
                          # trait intercept and trait covariance (PHIT1TRAIT in ctsem)                
                          mxPath(from=paste0("traitint",1:n.latent),
                                 to=paste0("trait",1:n.latent),
                                 arrows=2,free=TRUE,values=.1,connect="unique.bivariate",
                                 labels=c(indexMatrix(starttext="phiT1trait",dimension=n.latent,symmetric=FALSE,unique=TRUE))),
                          
                          mxPath(from=paste0("traitint",1:n.latent),#trait intercept variance
                                 arrows=2,free=ifelse(betweenconstraints==TRUE,F,T),
                                 values=c(diag(.2,n.latent)+.05)[lower.tri(diag(n.latent),diag=TRUE)],
                                 connect="unique.pairs",
                                 lbound=PHITRAITbounds[lower.tri(diag(n.latent),diag=TRUE)],
                                 labels=betweenphilabels)
      )
    }      
    
    larmaModel<-OpenMx::mxModel(larmaModel,#trait variance
                        mxPath(from=paste0("trait",1:n.latent),
                               arrows=2,free=FALSE,values=c(diag(.2,n.latent)+.05)[lower.tri(diag(n.latent),diag=TRUE)],connect="unique.pairs",lbound=PHITRAITbounds[lower.tri(diag(n.latent),diag=TRUE)],
                               #         labels=indexMatrix(starttext="phitrait",dimension=n.latent,symmetric=TRUE,unique=TRUE))
                               labels=indexMatrix(starttext="PHITRAIT[",sep=",",endtext="]",dimension=n.latent,symmetric=TRUE,unique=TRUE)),
                        
                        mxMatrix(name="PHITRAIT",labels=indexMatrix(starttext="phitrait",dimension=n.latent,symmetric=TRUE,unique=FALSE),
                                 free=TRUE,lbound=PHITRAITbounds,values=c(diag(.2,n.latent)+.05),nrow=n.latent,ncol=n.latent)
    )
  }
  
  
  
  
  
  
  
  #PREDICTORS 
  if(n.TDpred>0){
    for(i in 1:n.TDpred){
      for(j in 1:n.latent){
        
        if(causalpredictors==FALSE){
          larmaModel<-OpenMx::mxModel(larmaModel,
                              
                              mxPath(to=get(paste0("latentV",j))[-1],
                                     from=get(paste0("TDpredProcess",i)),
                                     connect="single",arrows=2,free=TRUE,values=.05,labels=paste0("covTDpred",j,"_",i))
          )
        }
        
        if(causalpredictors==TRUE){
          larmaModel<-OpenMx::mxModel(larmaModel,
                              
                              mxPath(to=get(paste0("latentV",j))[-1],
                                     from=get(paste0("TDpredProcess",i)),
                                     connect="single",arrows=1,free=TRUE,values=.5,labels=paste0("TDpred",j,"_",i))
          )
        }
        
        larmaModel<-OpenMx::mxModel(larmaModel,
                            mxPath(from="one",
                                   to=get(paste0("TDpredProcess",i)),
                                   connect="single",arrows=1,free=TRUE,values=.5,labels=paste0("MTDpred",i,"_",1:(Tpoints-1)))
        )
      }
    } #end npred and nlatent loops,remain in TDpred
    larmaModel<-OpenMx::mxModel(larmaModel,
                        mxPath(from=paste0("P",rep(1:n.TDpred,each=(Tpoints-1)),"_",1:(Tpoints-1)),
                               to=paste0("i",rep(1:n.latent,each=(Tpoints-1)*n.latent)),
                               connect="single",arrows=2,free=TRUE,values=.05,
                               labels=paste0("phiTDpred",rep(1:n.TDpred,each=(Tpoints-1)),"_","T",1:(Tpoints-1),"_T1_",rep(1:n.latent,each=((Tpoints-1)*n.latent))))
    )
    if(randomeffects==FALSE && traits==TRUE){
      larmaModel<-OpenMx::mxModel(larmaModel,
                          mxPath(from=paste0("P",rep(1:n.TDpred,each=((Tpoints-1)*n.latent)),"_",1:(Tpoints-1)),
                                 to=paste0("ci",rep(1:n.latent,each=(Tpoints-1))),
                                 connect="single",arrows=2,free=TRUE,values=.05,
                                 labels=paste0("phiTDpred",rep(1:n.TDpred,each=((Tpoints-1)*n.latent)),"T",1:(Tpoints-1),"_trait",rep(1:n.latent,each=(Tpoints-1))))
      )
    }
  } #end TDpred
  
  if(n.TIpred>0){
    for(i in 1:n.TIpred){
      for(j in 1:n.latent){
        larmaModel<-OpenMx::mxModel(larmaModel,
                            
                            mxPath(to=get(paste0("latentV",j))[-1],
                                   from=paste0("Z",i),
                                   connect="single",arrows=1,free=TRUE,values=.5,labels=paste0("TIpred",j,"_",i)),
                            
                            mxPath(from="one",
                                   to=paste0("Z",i),
                                   connect="single",arrows=1,free=TRUE,values=.5,labels=paste0("MTIpred",i))
        )
      }
    }
    larmaModel<-OpenMx::mxModel(larmaModel,
                        mxPath(from=paste0("Z",rep(1:n.TIpred,times=n.latent)),
                               to=paste0("i",rep(1:n.latent,each=n.TIpred)),
                               connect="single",arrows=2,free=TRUE,values=.05,
                               labels=paste0("phiTIpred_",rep(1:n.TIpred,times=n.latent),"T1_",rep(1:n.latent,each=n.TIpred)))
    )
    if(randomeffects==FALSE && traits==TRUE){
      larmaModel<-OpenMx::mxModel(larmaModel,
                          mxPath(from=paste0("Z",rep(1:n.TIpred,times=n.latent)),
                                 to=paste0("ci",rep(1:n.latent,each=n.TIpred)),
                                 connect="single",arrows=2,free=TRUE,values=.05,
                                 labels=paste0("phiTIpred_",rep(1:n.TIpred,times=n.latent),"trait",rep(1:n.latent,each=n.TIpred)))
      )
    }
  }
  
  if (n.TDpred+n.TIpred>0){
    larmaModel<-OpenMx::mxModel(larmaModel,
                        
                        mxPath(from=get(paste0("predictors")),
                               connect="unique.pairs",arrows=2,free=TRUE,
                               values=(diag(n.TDpred*(Tpoints-1)+n.TIpred)+.2)[lower.tri(diag(n.TDpred*(Tpoints-1)+n.TIpred),diag=TRUE)],
                               labels=paste0(indexMatrix(dimension=n.TDpred*(Tpoints-1)+n.TIpred,unique=TRUE,symmetrical=TRUE,starttext="phipred")))
    )
  }
  
  
  
  
  
  
  
  if(MA==TRUE){
    if(reasonable==TRUE) MAbounds<-c(-1,1)
    if(reasonable==FALSE) MAbounds<-c(-99,99)
    
    #MA1 cov with intercept
    larmaModel<-OpenMx::mxModel(larmaModel,
                        mxPath(from=paste0("Q",1:n.latent,"_1"),
                               to=paste0("i",1:n.latent),
                               arrows=2,free=TRUE,
                               values=.05,
                               connect="unique.bivariate",
                               labels=indexMatrix(dimension=n.latent,starttext=paste0("phiQ_"),unique=TRUE,symmetrical=TRUE))
    )
    
    for(j in (Tinit-1):(Tpoints-1)){
      larmaModel<-OpenMx::mxModel(larmaModel,   
                          
                          mxPath(from=paste0("Q",1:n.latent,"_",j),
                                 to=paste0("L",1:n.latent,"_",j+1),
                                 arrows=1,free=TRUE,
                                 values=.05,
                                 lbound=MAbounds[1],ubound=MAbounds[2],
                                 connect="single",
                                 labels=paste0("MA",1:n.latent,1:n.latent))
      )
    } 
  }
  #     if(MAcross==TRUE){
  #       for(j in 1:(Tpoints-1)){
  #         
  #         larmaModel<-OpenMx::mxModel(larmaModel,   
  #           
  #           mxPath(from=paste0("Q",1:n.latent,"_",j),
  #             to=paste0("L",1:n.latent,"_",j+1),
  #             arrows=1,free=TRUE,
  #             values=.05,
  #             connect="single",
  #             labels=paste0("MA",1:n.latent,1:n.latent))
  #         )
  #       } 
  #     }
  
  if(QAR==TRUE){
    print("QAR inits may not be working well")
    
    # QAR intercept to initialise
    larmaModel<-OpenMx::mxModel(larmaModel,   
                        
                        mxPath(from=paste0("QARi",1:n.latent),
                               to=paste0("Q",1:n.latent,"_",Tinit),
                               arrows=1,free=TRUE,
                               values=.05,
                               connect="single",
                               labels=paste0("QAR",1:n.latent)),
                        
                        #       #QAR var/cov
                        #       mxPath(from=paste0("QARi",1:n.latent),
                        #         arrows=2,free=TRUE,
                        #         values=.05,
                        #         connect="unique.pairs",
                        #         labels=indexMatrix(dimension=n.latent,starttext=paste0("QARi"),unique=TRUE,symmetrical=TRUE)),
                        
                        #QAR cov with initial within variance
                        mxPath(from=paste0("QARi",1:n.latent),
                               to=paste0("Q",1:n.latent,"_1"),
                               arrows=2,free=TRUE,
                               values=.05,
                               connect="unique.bivariate",
                               labels=indexMatrix(dimension=n.latent,starttext=paste0("phiQARi"),unique=TRUE,symmetrical=TRUE))
#                         labels=Qlabels [((i-2)*n.latent^2+1):((i-1)*n.latent^2)] )
#                                [lower.tri(diag(n.latent),diag=TRUE)] )
    )
    
    #QAR effect between time points
    for(j in Tinit:(Tpoints-1)){
      larmaModel<-OpenMx::mxModel(larmaModel,   
                          
                          mxPath(from=paste0("Q",1:n.latent,"_",j),
                                 to=paste0("Q",1:n.latent,"_",j+1),
                                 arrows=1,free=TRUE,
                                 values=.05,
                                 connect="single",
                                 labels=paste0("QAR",1:n.latent))
      )
    } 
  }
  
  
  if(detrend==TRUE){
    larmaModel<-OpenMx::mxModel(larmaModel,
                        #-----------------DETREND-----------------------#
                        
                        #random trend means
                        mxPath(from="one",
                               to=paste0("L",1:n.latent,"_",2:Tpoints),
                               arrows=1,free=TRUE,values=.1,
                               labels=c(paste0("meanL",1:n.latent,"_",2:Tpoints)))
    )
  }
  ##objective
  #   objectiveIlength<-nrow(larmaModel$A@values)
  #   larmaModel<-OpenMx::mxModel(larmaModel,
  #     mxAlgebra(F%*%solve(objectiveI-A)%*%S%*%t(solve(objectiveI-A))%*%t(F),name="expcov"),
  #     mxAlgebra(t(F%*%(solve(objectiveI-A))%*%t(M)),name="expmeans"),
  #     mxFIMLObjective(covariance="expcov",means="expmeans",dimnames=paste0("V",1:(n.manifest*Tpoints)),vector=FALSE),
  #     mxMatrix(type="Iden",nrow=objectiveIlength,ncol=objectiveIlength,name="objectiveI")
  #   ) 
  
  #new objective
  #   objectiveIlength<-nrow(larmaModel$A@values)
  #   larmaModel<-OpenMx::mxModel(larmaModel,
  #     mxAlgebra(F%*%solve(objectiveI-A)%*%S%*%t(solve(objectiveI-A))%*%t(F),name="expcov"),
  #     mxAlgebra(t(F%*%(solve(objectiveI-A))%*%t(M)),name="expmeans"),
  #      mxExpectationNormal(covariance="expcov",means="expmeans",dimnames=paste0("V",1:(n.manifest*Tpoints))),
  #     mxMatrix(type="Iden",nrow=objectiveIlength,ncol=objectiveIlength,name="objectiveI"),
  #      mxFitFunctionML()
  #   )
  
  #custom objective
  #   objectiveIlength<-nrow(larmaModel$A@values)
  #   larmaModel<-OpenMx::mxModel(larmaModel,type="default",
  #         mxAlgebra(F%*%solve(objectiveI-A)%*%S%*%t(solve(objectiveI-A))%*%t(F),name="expcov"),
  #         mxAlgebra(t(F%*%(solve(objectiveI-A))%*%t(M)),name="expmeans"),
  #     mxMatrix(type="Iden",nrow=objectiveIlength,ncol=objectiveIlength,name="objectiveI"),
  #   mxFitFunctionR(CondLogLik))
  
  #   #row objective
  #   larmaModel<-addmxrowobjective(larmaModel,paste0("V",1:(n.manifest*Tpoints)))
  #   larmaModel<-OpenMx::mxModel(larmaModel,
  #     mxMatrix(type="Iden",nrow=nrow(larmaModel$A@values),ncol=nrow(larmaModel$A@values),name="objectiveI"),
  #       mxAlgebra(F%*%solve(objectiveI-A)%*%S%*%t(solve(objectiveI-A))%*%t(F),name="expCov"),
  #       mxAlgebra(t(F%*%(solve(objectiveI-A))%*%t(M)),name="expMean")
  #   )
  
  
  
  return(larmaModel)
}


