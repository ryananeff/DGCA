#' @title Calculate correlation matrix p-values.
#' @description Calculate two-sided p-values from a pairwise correlations matrix and a corresponding "number of samples" matrix.
#' @param corrs Computed correlation matrix.
#' @param nsamp Computed number of samples used per call in the correlation matrix.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @return A matrix of p-values.
#' @references HMisc R package https://cran.r-project.org/web/packages/Hmisc/index.html
#' @export
matCorSig <- function(corrs, nsamp, secondMat = FALSE){
	if(!secondMat){
		pvals = matrix(NA, nrow = nrow(corrs), ncol = ncol(corrs))
		#pvals_before <<- pvals
		#corrs_before <<- corrs
		corrs = corrs[upper.tri(corrs)] #computed correlation matrix
		nsamp = nsamp[upper.tri(nsamp)] #number of samples used per call in the correlation matrix
		deg_freedom = nsamp - 2 #degrees of freedom
		deg_freedom = deg_freedom[1]
		in_pt = (abs(corrs) * sqrt(deg_freedom) / sqrt(1 - corrs^2))
		in_pt <<- in_pt
		deg_freedom <<- deg_freedom
		pvals_upper = 2 * (1 - pt(in_pt, deg_freedom))
		#### pvals_upper = 2 * (1 - pt(in_pt, df)) # pt gives the distribution function
		pvals[upper.tri(pvals)] = pvals_upper 
		pvals[abs(corrs) == 1] = 0

	} else {

		df = nsamp - 2
		pvals = matrix(2 * (1 - pt(abs(corrs) * sqrt(df) /
			sqrt(1 - corrs^2), df)), nrow = nrow(corrs))
		pvals[abs(corrs) == 1] = 0

	}
	return(pvals)
}

#' @title Calculate correlation matrix p-values, adding a small noise term to the result to ensure normality.
#' @description Calculate two-sided p-values from a pairwise correlations matrix and a corresponding "number of samples" matrix, adding a small random noise term which is defined as epsilon * rand(0,1).
#' @param corrs Computed correlation matrix.
#' @param nsamp Computed number of samples used per call in the correlation matrix.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @param epsilon A number representing the maximum noise to add to each p-value. Default is 1e-6.
#' @return A matrix of p-values.
#' @references HMisc R package https://cran.r-project.org/web/packages/Hmisc/index.html; Yang et al., A Fuzzy Permutation Method for False Discovery Rate Control, Nat Sci Rep (2016), https://www.nature.com/articles/srep28507
#' @export
matCorSig.fuzz <- function(corrs, nsamp, secondMat = FALSE, epsilon=1e-6){
	if(!secondMat){
		pvals = matrix(NA, nrow = nrow(corrs), ncol = ncol(corrs))
		#pvals_before <<- pvals
		#corrs_before <<- corrs
		corrs = corrs[upper.tri(corrs)] #computed correlation matrix
		nsamp = nsamp[upper.tri(nsamp)] #number of samples used per call in the correlation matrix
		deg_freedom = nsamp - 2 #degrees of freedom
		deg_freedom = deg_freedom[1]
		in_pt = (abs(corrs) * sqrt(deg_freedom) / sqrt(1 - corrs^2))
		in_pt <<- in_pt
		deg_freedom <<- deg_freedom
		pvals_upper = 2 * (1 - pt(in_pt, deg_freedom))
		#### pvals_upper = 2 * (1 - pt(in_pt, df)) # pt gives the distribution function
		pvals[upper.tri(pvals)] = pvals_upper 
		pvals[abs(corrs) == 1] = 0

	} else {

		df = nsamp - 2
		pvals = matrix(2 * (1 - pt(abs(corrs) * sqrt(df) /
			sqrt(1 - corrs^2), df)), nrow = nrow(corrs))
		pvals[abs(corrs) == 1] = 0

	}

	## fuzzy implementation of correlation function
	## see Supp Info S3 at URL: https://www.nature.com/articles/srep28507
	pvalsf = pvals + runif(length(pvals),0,epsilon) #add the uniform fuzz to the correlation matrix
	pvalsfm = structure(sapply(pvalsf,function(x){max(min(x,1),0)}),dim=dim(pvalsf)) #ensure P is still in range [0,1]
	rownames(pvalsfm) = rownames(pvals)
	colnames(pvalsfm) = colnames(pvals)

	return(pvalsfm)
}

###

