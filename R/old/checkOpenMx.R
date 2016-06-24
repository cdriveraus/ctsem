checkOpenMx <- function(callingFunc){
  if("OpenMx" %in% installed.packages(noCache=TRUE)[,'Package'] == TRUE){ 
    if(packageVersion('OpenMx') >= '2.0.0.0') openmx <- 'installed'
    
    if(packageVersion('OpenMx') < '2.0.0.0'){
    openmx<-readline("ctsem requires OpenMx 2.0.0.0 or greater, do you want to update OpenMx? y/n \n")
  }
  }
    
    
    if("OpenMx" %in% installed.packages(noCache=TRUE)[,'Package'] == FALSE){
        openmx<-readline("ctsem requires OpenMx 2.0.0.0 or greater, do you want to install OpenMx? y/n \n")
      }
      
    
  if(openmx=='y' | openmx=='Y') {
    source('http://openmx.psyc.virginia.edu/getOpenMx.R')
    openmx<-'installed'
  }
  
  if(openmx !='installed') stop(paste0(callingFunc,' cannot proceed without OpenMx 2.x'))
}