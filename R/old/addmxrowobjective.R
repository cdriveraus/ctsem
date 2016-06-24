#'adds a row objective FIML objective function to an openmx model
#'@param inputmodel model to add objective function to
#'@param dimlabels vector of labels for observed data

addmxrowobjective<-function(inputmodel,dimlabels){
  rowobjmodel <- OpenMx::mxModel(model=inputmodel,
#   mxAlgebra(expression=omxSelectRowsAndCols(expcov, existenceVector),name="filteredExpCov"),
#   mxAlgebra(expression=omxSelectCols(expmeans, existenceVector),name="filteredExpMean"),
#   mxAlgebra(expression=sum(existenceVector),name="obscount"),
#   mxAlgebra(expression=(log(2*pi))*obscount+log(det(filteredExpCov))+
#       (filteredDataRow - filteredExpMean)%&%solve(filteredExpCov),name="rowlikelihood"),
#   mxAlgebra(expression=sum(rowResults),name="reduceAlgebra"),
    
    
    mxMatrix("Full", 1, 1, values = log(2*pi), name = "log2pi"),
    mxAlgebra(
     	expression=omxSelectRowsAndCols(expCov, existenceVector),
     	name="filteredExpCov",
     	),
    	mxAlgebra(
     	expression=omxSelectCols(expMean, existenceVector),
     	name="filteredExpMean",
     	),
     	mxAlgebra(
     	expression=log2pi %*% 2 + log(det(filteredExpCov)),
     	name ="firstHalfCalc",
     	),
    	mxAlgebra(
     	expression=(filteredDataRow - filteredExpMean) %&% solve(filteredExpCov),
     	name = "secondHalfCalc",
     	),
    	mxAlgebra(
     	expression=(firstHalfCalc + secondHalfCalc),
     	name="rowAlgebra",
     	),
     	mxAlgebra(
     	expression=sum(rowResults),
     	name = "reduceAlgebra",
     	),
    	mxFitFunctionRow(
     	rowAlgebra='rowAlgebra',
     	reduceAlgebra='reduceAlgebra',
     	dimnames=dimlabels)
    
#   mxFitFunctionRow(rowAlgebra='rowlikelihood',reduceAlgebra='reduceAlgebra',dimnames=dimlabels)
  )
  return(rowobjmodel)
}

# k*log(2*pi)+log(detexpcovfilt)+t(x-expmeansfilt)%*%solve(expcovfilt)%*%(x-expmeansfilt)