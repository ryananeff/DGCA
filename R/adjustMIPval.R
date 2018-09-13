#' @title Adjusts MI p-values via permutation testing to empirical p-values.
#' @description Wraps around bigEmpPVals to format the p-values returned from the mutual information calculation step properly as input.
#' @param pVals Matrix of p-values.
#' @param corPermMat Matrix of permuted p-values.
#' @param secondMat Whether two matrices are being compared for differential correlation. Default = FALSE.
#' @return A numeric vector of p-values that have been adjusted empirically. 
#' @export
adjustMIPval <- function(corrs, corPermMat, secondMat=FALSE){

  if(!secondMat) corrs = corrs[upper.tri(corrs)]
  corrs = as.numeric(corrs)
  
  if(!secondMat) corPermMat = corPermMat[apply(corPermMat, 3, upper.tri)]
  corrs0 = as.numeric(corPermMat)
  if(is.matrix(corrs0)) corrs0 = as.vector(corrs0)

  pvalues = bigEmpPVals(stat=corrs,stat0=corrs0,increasing=TRUE)

  return(pvalues)
}
