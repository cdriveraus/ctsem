ctPackageComparison<-function(){
output<-matrix(NA,10,6)
colnames(output)<-c('True', 'ctsem','cts', 'PSM','yuima', 'arima')
for(i in 1:nrow(output)){
generatingModel<-ctModel(n.latent=1,n.manifest=1,Tpoints=500,LAMBDA=diag(1),DRIFT=matrix(-.3,nrow=1),
  CINT=matrix(3,1,1),
  MANIFESTVAR=diag(0,1),
  DIFFUSION=t(chol(diag(5,1))))

output[i,1]<-generatingModel$DRIFT

ctsemData<-ctGenerate(generatingModel,n.subjects=1,burnin=300)
longData<-ctWideToLong(ctsemData,Tpoints=500,n.manifest=1)
longData<-ctDeintervalise(longData)

ctsemModel<-ctModel(n.latent=1,n.manifest=1,Tpoints=500,
  MANIFESTVAR=diag(0,1),
  LAMBDA=diag(1))

ctsemFit<-ctFit(ctsemData,ctsemModel,stationary=c('T0VAR'))
# summary(ctsemFit)
output[i,2]<-mxEval(DRIFT,ctsemFit$mxobj)

### fit using CTS package for comparison
ctsData<-longData[,c('AbsTime', 'Y1')]
library(cts)
ctsFit<-car(ctsData,order=1, scale=1)
ctsFit$phi #transformation of DRIFT parameter
output[i,3]<- -1 * (1+ctsFit$phi) / (1-ctsFit$phi)

### fit using PSM package for comparison
psmFit<-ctPSMfit(ctsemData,omxStartValues=omxGetParameters(ctsemFit$mxobj), ctsemModel)
psmFit$PSMfit$opt$par #parameter estimates comparable to raw openmx parameter estimates of ctsem summary.
psmFit$PSMfit$opt$value *2 #multiply PSM likelihood by 2 to compare to ctsem.
output[i,4]<- -exp(psmFit$PSMfit$opt$par[2])


#yuima
library(yuima)
mod <- setModel(drift="drift*x+cint", diffusion="diffusion")
ou <- setYuima(model=mod, data=setData(longData[,'Y1'], delta=1))
mlout<-qmle(ou,start=list(drift=-.3,diffusion=1, cint=1))
# summary(mlout)
output[i,5]<-mlout@coef[2]

#arima (discrete time analysis only)
arfit<-arima(longData[,'Y1'], order=c(1,0,0))
arfit$coef #comparable to discreteDRIFT and asymCINT from ctsem summary
log(arfit$coef[1]) #transform ar1 parameter to continuous drift parameter
output[i,6]<-log(arfit$coef[1])

if(i > 1){
plot(density(output[1:i,2]),xlim=c(-.5,0), lty=6, lwd=2, main='Density of estimates of drift param (true value -0.3) \n using various packages')
points(density(output[1:i,3]),col='red',type='l', lty=2, lwd=2)
points(density(output[1:i,4]),col='blue',type='l', lty=3, lwd=2)
points(density(output[1:i,5]),col='green',type='l', lwd=2, lty=4)
points(density(output[1:i,6]),col='purple',type='l', lwd=2, lty=5)
legend('topleft', legend=c('ctsem','cts', 'PSM','yuima', 'arima'), bty='n',
  text.col=c('black','red','blue','green','purple'))
}

}



}