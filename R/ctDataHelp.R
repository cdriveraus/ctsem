#' AnomAuth
#'
#' A dataset containing panel data assessments of individuals Anomia and Authoritarianism.
#' @format data frame with 2722 rows, 14 columns. Column Y1 represents anomia, Y2 Authoritarianism, dTx the time interval for measurement occasion x.
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
#' 8 measurement occasions of leisure time and happiness, 7 measurement occasions of a money intervention dummy,
#' and 7 measurement intervals for each of 50 individuals.
#' @name ctExample2
NULL

#' ctExample3
#'
#' Simulated example dataset for the ctsem package
#' @format 1 by 399 matrix containing containing ctsem wide format data. 100 observations of variables Y1 and Y2 and 199 measurement intervals, for 1 subject.
#' @name ctExample3
NULL

#' ctExample4
#'
#' Simulated example dataset for the ctsem package
#' @format 20 by 79 matrix containing 20 observations of variables Y1, Y2, Y3, and 19 measurement intervals dTx, for each of 20 individuals.
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
#' @format 100 by 18 matrix containing ctsem wide format data. 8 measurement occasions of leisure time and happiness, 7 measurement occasions of a money intervention dummy,
#' and 7 measurement intervals for each of 50 individuals.
#' @name ctExample2level
NULL


#' datastructure
#'
#' Simulated example dataset for the ctsem package
#' @format 2 by 15 matrix containing containing ctsem wide format data. 
#' 3 measurement occasions of manifest variables Y1 and Y2, 
#' 2 measurement occasions of time dependent predictor TD1, 
#' 2 measurement intervals dTx, and 2 time independent predictors TI1 and TI2, for 2 individuals.
#' @name datastructure
NULL


#' longexample
#'
#' Simulated example dataset for the ctsem package
#' @format 7 by 8 matrix containing ctsem long format data, for two subjects, 
#' with three manifest variables Y1, Y2, Y3, 
#' one time dependent predictor TD1, two time independent predictors TI1 and TI2, 
#' and absolute timing information Time.
#' @name longexample
NULL

#' ctstantestfit
#' 
#' Minimal output from \link{\code{ctStanFit}} from ctsem package.
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
#' testm<-ctModel(type='omx',Tpoints=Tpoints,n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,n.manifest=n.manifest,
#'   MANIFESTVAR=diag(0.5,2),
#'   TIPREDEFFECT=matrix(c(0,0,0,0,0,0),nrow=2),
#'   TIPREDVAR=matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3),
#'   TDPREDEFFECT=matrix(c(.1,-.2),nrow=2),
#'   TDPREDVAR=matrix(0,nrow=n.TDpred*(Tpoints-1),ncol=n.TDpred*(Tpoints-1)),
#'   TDPREDMEANS=matrix(rnorm(n.TDpred*(Tpoints-1),0,1),nrow=n.TDpred*(Tpoints-1)),
#'   LAMBDA=diag(1,2),
#'   DRIFT=matrix(c(-.3,.2,-.1,-.2),nrow=2),
#'   DIFFUSION=matrix(c(.3,.1,0,.2),2),CINT=matrix(c(0,0),nrow=2),T0MEANS=matrix(0,ncol=1,nrow=2),
#'   T0VAR=diag(100,2))
#' cd<-ctGenerate(testm,n.subjects=n.subjects,burnin=300,simulTDpredeffect=TRUE)
#' 
#' dlong<-ctWideToLong(cd,Tpoints,n.manifest=n.manifest, n.TDpred = n.TDpred, n.TIpred = n.TIpred)
#' dlong<-ctDeintervalise(dlong)
#' dlong<-dlong[c(1:10,18:40,48:80,90:120,140:150),]
#' 
#' checkm<-ctModel(type='stanct',Tpoints=Tpoints,
#'   n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,
#'   n.manifest=n.manifest,LAMBDA=diag(2))
#'   
#' checkm$parameters$indvarying[-1:-2] <- FALSE
#' 
#' ctstantestfit<-ctStanFit(dlong,checkm,iter=20,chains=1,initwithoptim=TRUE)
#' save(ctstantestfit,file='.\\data\\ctstantestfit.rda')
#' }
NULL


