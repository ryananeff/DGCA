adjustMIPVal <- function(corrs, corPermMat, secondMat=FALSE){

  if(!secondMat) corrs = corrs[upper.tri(corrs)]
  corrs = as.numeric(corrs)
  
  if(!secondMat) corPermMat = corPermMat[apply(corPermMat, 3, upper.tri)]
  corrs0 = as.numeric(corPermMat)
  if(is.matrix(corrs0)) corrs0 = as.vector(corrs0)

  corrs_gt_0 = corrs[corrs>0]
  corrs0_gt_0 = corrs0[corrs0>0]
  corrs_lt_0 = abs(corrs[corrs<0]) #correct them so they are increasing
  corrs0_lt_0 = abs(corrs0[corrs0<0]) #correct them so they are increasing

  pvalues_upper = bigEmpPVals(stat=corrs_gt_0,stat0=corrs0_gt_0,increasing=TRUE)
  pvalues_lower = bigEmpPVals(stat=corrs_lt_0,stat0=corrs0_lt_0,increasing=TRUE)

  pvalues = rep(NA,length(corrs))
  pvalues[corrs>0] = pvalues_upper
  pvalues[corrs<0] = pvalues_lower
  pvalues[corrs==0] = 1 #no MI therefore p-value of 1

  return(pvalues)
}