#' @title Differential correlation between two conditions.
#' @description Z-transforms correlation coefficients (Pearson or Spearman) and then calculates the difference in z-scores between the two conditions divided by the square root of the standard errors (which is inversely proportion to the sample sizes used to calculate the correlations).
#' @param rho1 Numeric vector of correlation coefficients in condition 1.
#' @param rho2 Numeric vector of correlation coefficients in condition 2.
#' @param n1 Numeric vector of the number of samples used in the correlation calculations in condition 1.
#' @param n2 Numeric vector of the number of samples used in the correlation calculations in condition 2.
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman".
#' @param pval1 When running in mutual information mode, do you have empirical p-values for matrix 1? If so, use them to recalculate the z-score difference. Default = "none" 
#' @param pval2 When running in mutual information mode, do you have empirical p-values for matrix 2? If so, use them to recalculate the z-score difference. Default = "none" 
#' @references Tests For Rank Correlation Coefficients, I. http://biomet.oxfordjournals.org/content/44/3-4/470.full.pdf+html
#' @return Numeric vector with scaled difference in z-scores of correlations between the two conditions.
#' @examples
#' rho1 = runif(100, -1, 1); rho2 = runif(100, -1, 1)
#' n1 = rep(100, 100); n2 = rep(110, 100)
#' dcorrs_res = dCorrs(rho1, n1, rho2, n2)
#' @export
dCorrs <- function(rho1, n1, rho2, n2, corrType = "pearson", pval1="none", pval2="none"){

	if(!(all.equal(length(rho1), length(n1), length(rho2), length(n2)))) stop("All of the input vectors must be the same length.")
	if(!identical(pval1,"none")&!identical(pval2,"none")){
		if(!(all.equal(length(rho1), length(n1), length(rho2), length(n2), length(pval1),length(pval2)))) stop("All of the input vectors must be the same length.")
	}

	if(corrType == "pearson"){
		#http://biomet.oxfordjournals.org/content/44/3-4/470.full.pdf+html
		zr1 = atanh(rho1)
		zr2 = atanh(rho2)
		diff12 = (zr2 - zr1)/sqrt((1/(n1 - 3)) + (1/(n2 - 3)))
	}
	if(corrType == "spearman"){
		#http://biomet.oxfordjournals.org/content/44/3-4/470.full.pdf+html
		zr1 = atanh(rho1)
		zr2 = atanh(rho2)
		diff12 = (zr2 - zr1)/sqrt((1.06/(n1 - 3)) + (1.06/(n2 - 3)))
	}
	if(corrType == "bicor"){
		#http://biomet.oxfordjournals.org/content/44/3-4/470.full.pdf+html
		zr1 = atanh(rho1)
		zr2 = atanh(rho2)
		diff12 = (zr2 - zr1)/sqrt((1/(n1 - 3)) + (1/(n2 - 3)))
	}
	if(corrType == "mutualinformation"){
		#https://en.wikipedia.org/wiki/Mutual_information#Linear_correlation
		if(identical(pval1,"none")|identical(pval2,"none")){
			#convert MI to pearson correlation coefficients under bivariate normal assumptions
			message("Calculating first-pass MI z-score difference using bivariate normal assumptions")
			rho1 = sqrt(1-exp(1)**(-2*rho1)) 
			rho2 = sqrt(1-exp(1)**(-2*rho2))
			zr1 = atanh(rho1)
			zr2 = atanh(rho2)
			diff12 = (zr2 - zr1)/sqrt((1/(n1 - 3)) + (1/(n2 - 3)))
		} else {
			message("Recalculating MI z-score difference using empirical p-values")
			pval1 = pval1[upper.tri(pval1)]
			pval2 = pval2[upper.tri(pval2)]
			zr1 = pval1/abs(pval1)*qnorm(1-(pval1/2)) #reverse z-score calculation from empirical p-values
			zr2 = pval2/abs(pval2)*qnorm(1-(pval2/2)) #reverse z-score calculation from empirical p-values
			diff12 = (zr2 - zr1)/sqrt((1/(n1 - 3)) + (1/(n2 - 3)))
		}
	}

	return(diff12)

}
