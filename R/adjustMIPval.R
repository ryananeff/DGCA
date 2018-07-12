adjustMIPval <- function(corrs, corPermMat, secondMat=FALSE){

  if(!secondMat) corrs = corrs[upper.tri(corrs)]
  corrs = as.numeric(corrs)
  
  if(!secondMat) corPermMat = corPermMat[apply(corPermMat, 3, upper.tri)]
  corrs0 = as.numeric(corPermMat)
  if(is.matrix(corrs0)) corrs0 = as.vector(corrs0)

  pvalues = bigEmpPVals(stat=corrs,stat0=corrs0,increasing=TRUE)

  return(pvalues)
}