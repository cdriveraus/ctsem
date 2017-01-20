#' AnomAuth
#'
#' A dataset containing panel data assessments of individuals Anomia and Authoritarianism.
#' @format data frame with 2722 rows, 14 columns. Column Y1 represents anomia, 
#' Y2 Authoritarianism, dTx the time interval for measurement occasion x.
#' @source See \url{http://psycnet.apa.org/journals/met/17/2/176/} for details.
#' @name AnomAuth
NULL


#' Oscillating
#'
#' Simulated example dataset for the ctsem package.
#' @format 200 by 21 matrix containing containing ctsem wide format data. 
#' 11 measurement occasions and 10 measurement intervals for each of 200 individuals
#' @source See \url{http://onlinelibrary.wiley.com/doi/10.1111/j.2044-8317.2012.02043.x/abstract}
#' @name Oscillating
NULL



#' ctExample1
#'
#' Simulated example dataset for the ctsem package
#' @format 100 by 17 matrix containing containing ctsem wide format data. 
#' 6 measurement occasions of leisure time and happiness and 5 measurement intervals for each of 100 individuals.
#' @name ctExample1
NULL


#' ctExample2
#'
#' Simulated example dataset for the ctsem package
#' @format 100 by 18 matrix containing containing ctsem wide format data. 
#' 8 measurement occasions of leisure time and happiness, 
#' 7 measurement occasions of a money intervention dummy,
#' and 7 measurement intervals for each of 50 individuals.
#' @name ctExample2
#' @examples
#' \dontrun{
#' #two process, one time dependent predictor example
#' Tpoints=20
#' manifestNames<-c('LeisureTime','Happiness')
#' TDpredNames<-'MoneyInt'
#' testm<-ctModel(Tpoints=Tpoints,n.latent=3,n.TDpred=1,n.TIpred=0,n.manifest=2,    
#'   LAMBDA=cbind(diag(1,2),0),
#'   MANIFESTVAR=diag(.1,2),
#'   DRIFT=matrix(c(-.3,.12,0,  -.02,-.3,0, 1,-.3,-.0001  ),nrow=3,ncol=3),
#'   TRAITVAR=t(chol(matrix(c(.2,-.1,0,  -.1,.21,0,  0,0,0.00001),ncol=3,nrow=3))),
#'   DIFFUSION=t(chol(diag(c(1.2,.6,0.0001),3))),
#'   CINT=matrix(c(1,.3,0),nrow=3),
#'   T0MEANS=matrix(0,ncol=1,nrow=3),
#'   T0VAR=diag(c(1,1,0),3),
#'   TDPREDEFFECT=matrix(c(.6,.4,1),nrow=3),
#'   TDPREDVAR=diag(c(rep(0,Tpoints)),Tpoints),
#'   TDPREDMEANS=matrix(c(0,0,0,0,0,1,rep(0,Tpoints-6)),ncol=1,nrow=(Tpoints)))
#' testd<-ctGenerate(testm,n.subjects=10,burnin=10) #generate data
#' 
#' ctIndplot(testd,Tpoints=Tpoints,n.manifest=2,n.subjects=10,colourby="variable")
#' 
#' timestokeep=c(0,1,4,5,7,8,16,19)
#' deltaT<-timestokeep[-1] - timestokeep[-8]
#' testd<-testd[,c(paste0('Y',1:2,'_T',rep(timestokeep,each=2)),paste0('TD1_T',timestokeep))]
#' testd<-cbind(testd,matrix(deltaT,nrow=nrow(testd),ncol=length(deltaT),byrow=TRUE))
#' 
#' colnames(testd)<-ctWideNames(n.manifest=2,Tpoints=8,n.TDpred=1,
#' manifestNames=manifestNames,TDpredNames=TDpredNames)
#' ctExample2<-testd
#' save(ctExample2,file=".\\data\\ctExample2.rda") 
#' }
NULL

#' ctExample3
#'
#' Simulated example dataset for the ctsem package
#' @format 1 by 399 matrix containing containing ctsem wide format data. 
#' 100 observations of variables Y1 and Y2 and 199 measurement intervals, for 1 subject.
#' @name ctExample3
NULL

#' ctExample4
#'
#' Simulated example dataset for the ctsem package
#' @format 20 by 79 matrix containing 20 observations of variables 
#' Y1, Y2, Y3, and 19 measurement intervals dTx, for each of 20 individuals.
#' @name ctExample4
NULL

#' ctExample1TIpred
#'
#' Simulated example dataset for the ctsem package
#' @format 100 by 18 matrix containing containing ctsem wide format data. 
#' 6 measurement occasions of leisure time and happiness, 1 measurement of number of friends,
#' and 5 measurement intervals for each of 100 individuals.
#' @name ctExample1TIpred
NULL


#' ctExample2level
#'
#' Simulated example dataset for the ctsem package
#' @format 100 by 18 matrix containing ctsem wide format data. 
#' 8 measurement occasions of leisure time and happiness, 
#' 7 measurement occasions of a money intervention dummy,
#' and 7 measurement intervals for each of 50 individuals.
#' @name ctExample2level
NULL


#' datastructure
#'
#' Simulated example dataset for the ctsem package
#' @format 2 by 15 matrix containing containing ctsem wide format data. 
#' 3 measurement occasions of manifest variables Y1 and Y2, 
#' 2 measurement occasions of time dependent predictor TD1, 
#' 2 measurement intervals dTx, and 2 time independent predictors 
#' TI1 and TI2, for 2 individuals.
#' @name datastructure
#' @examples
#' \dontrun{
#' Tpoints=30
#' testm<-ctModel(Tpoints=Tpoints,n.latent=1,n.TDpred=1,n.TIpred=2,n.manifest=3,    
#'   LAMBDA=matrix(1,ncol=1,nrow=3),
#'   DRIFT=diag(-.3,1),
#'   DIFFUSION=diag(.1,1),
#'   CINT=diag(2,1),
#'   MANIFESTVAR=diag(1,3),
#'   TDPREDEFFECT=diag(.2,1),
#'   TIPREDEFFECT=matrix(.8,nrow=1,ncol=2),
#'   TDPREDVAR=diag(1,1*(Tpoints)),
#'   TIPREDVAR=diag(1,2)
#' )
#' longexample<-round(ctGenerate(testm,n.subjects=2,logdtsd = 1,burnin=3,wide=FALSE)[c(1:3,32:34),],2)
#' longexample[2,c(2,7)]<-NA
#' longexample[4,c(3)]<-NA
#' datastructure <- ctLongToWide(datalong = longexample,id='id',time='time',
#'   manifestNames = testm$manifestNames,TDpredNames = testm$TDpredNames,
#'   TIpredNames=testm$TIpredNames)
#' datastructure<-ctIntervalise(datawide = datastructure,
#'   Tpoints = 3,n.manifest = testm$n.manifest,n.TDpred = testm$n.TDpred,
#'   n.TIpred=testm$n.TIpred)
#' save(datastructure,file='.\\data\\datastructure.rda')
#' }
NULL


#' longexample
#'
#' Simulated example dataset for the ctsem package
#' @format 7 by 8 matrix containing ctsem long format data, for two subjects, 
#' with three manifest variables Y1, Y2, Y3, 
#' one time dependent predictor TD1, two time independent predictors TI1 and TI2, 
#' and absolute timing information Time.
#' @name longexample
#' @examples
#' \dontrun{
#' #long example (using datastructure base)
#' Tpoints=30
#' testm<-ctModel(Tpoints=Tpoints,n.latent=1,n.TDpred=1,n.TIpred=2,n.manifest=3,    
#'   LAMBDA=matrix(1,ncol=1,nrow=3),
#'   DRIFT=diag(-.3,1),
#'   DIFFUSION=diag(.1,1),
#'   CINT=diag(2,1),
#'   MANIFESTVAR=diag(1,3),
#'   TDPREDEFFECT=diag(.2,1),
#'   TIPREDEFFECT=matrix(.8,nrow=1,ncol=2),
#'   TDPREDVAR=diag(1,1*(Tpoints)),
#'   TIPREDVAR=diag(1,2)
#' )
#' longexample<-round(ctGenerate(testm,n.subjects=2,logdtsd = 1,burnin=3,wide=FALSE)[c(1:3,32:35),],2)
#' longexample[2,c(2,7)]<-NA
#' longexample[4,c(3)]<-NA
#' longexample
#' save(longexample,file='.\\data\\longexample.rda')
#' }
NULL

#' ctstantestfit
#' 
#' Minimal output from \code{\link{ctStanFit}} from ctsem package.
#' @format stanfit class.
#' @name ctstantestfit
#' @examples 
#' \dontrun{
#' ### generator for ctstantestfit
#' set.seed(2)
#' Tpoints=50
#' n.manifest=2
#' n.TDpred=0
#' n.TIpred=3
#' n.latent=2
#' n.subjects=3
#' 
#' testm<-ctModel(type='omx',Tpoints=Tpoints,n.latent=n.latent,
#' n.TDpred=n.TDpred,n.TIpred=n.TIpred,n.manifest=n.manifest,
#'   MANIFESTVAR=diag(0.5,2),
#'   TIPREDEFFECT=matrix(c(0,0,0,0,0,0),nrow=2),
#'   TIPREDVAR=matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3),
#'   TDPREDEFFECT=matrix(c(.1,-.2),nrow=2),
#'   TDPREDVAR=matrix(0,nrow=n.TDpred*(Tpoints-1),ncol=n.TDpred*(Tpoints-1)),
#'   TDPREDMEANS=matrix(rnorm(n.TDpred*(Tpoints-1),0,1),nrow=n.TDpred*(Tpoints-1)),
#'   LAMBDA=diag(1,2),
#'   DRIFT=matrix(c(-.3,.2,-.1,-.2),nrow=2),
#'   DIFFUSION=matrix(c(.3,.1,0,.2),2),CINT=matrix(c(0,0),nrow=2),
#'   T0MEANS=matrix(0,ncol=1,nrow=2),
#'   T0VAR=diag(100,2))
#' cd<-ctGenerate(testm,n.subjects=n.subjects,burnin=300,simultdpredeffect=TRUE,wide=FALSE)
#' 
#' checkm<-ctModel(type='stanct',Tpoints=Tpoints,
#'   n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,
#'   n.manifest=n.manifest,LAMBDA=diag(2))
#'   
#' checkm$pars$indvarying[-1:-2] <- FALSE
#' 
#' ctstantestfit<-ctStanFit(cd,checkm,iter=20,chains=1,initwithoptim=TRUE)
#' save(ctstantestfit,file='.\\data\\ctstantestfit.rda')
#' }
NULL




#' ctstantestdat
#' 
#' Generated dataset for testing \code{\link{ctStanFit}} from ctsem package.
#' @format matrix
#' @name ctstantestdat
#' @examples 
#' \dontrun{
#' Tpoints=50
#' n.manifest=2
#' n.TDpred=1
#' n.TIpred=3
#' n.latent=2
#' n.subjects=5
#' gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
#' n.TDpred=n.TDpred,n.TIpred=n.TIpred,n.manifest=n.manifest,
#'   MANIFESTVAR=diag(0.5,2),
#'   TIPREDEFFECT=matrix(c(.5,0,0,-.5,0,0),nrow=2),
#'   TIPREDVAR=matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3),
#'   TDPREDEFFECT=matrix(c(.1,-.2),nrow=2),
#'   TDPREDVAR=matrix(0,nrow=n.TDpred*(Tpoints-1),ncol=n.TDpred*(Tpoints-1)),
#'   TDPREDMEANS=matrix(rnorm(n.TDpred*(Tpoints-1),0,1),
#'    nrow=n.TDpred*(Tpoints-1)),
#'   LAMBDA=diag(1,2),
#'   DRIFT=matrix(c(-.3,.2,-.1,-.2),nrow=2),
#'   TRAITVAR=t(chol(matrix(c(4,3,3,4),nrow=2))),
#'   DIFFUSION=matrix(c(.3,.1,0,.2),2),CINT=matrix(c(0,0),nrow=2),
#'   T0MEANS=matrix(0,ncol=1,nrow=2),
#'   T0VAR=diag(100,2))
#' 
#' ctstantestdat<-ctGenerate(gm,n.subjects=n.subjects,burnin=30,
#' wide=FALSE, simultdpredeffect = TRUE)
#' save(ctstantestdat,file='.\\data\\ctstantestdat.rda')
#' paths <- sort(Sys.glob(c("data/*.rda", "data/*.RData")))
#' library(tools)
#' resaveRdaFiles(paths)
#' }
NULL

