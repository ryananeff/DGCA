#' @title Adjusts a numeric vector of p-values.
#' @description Wraps around the R base implementation of p.adjust, the methods used in the fdr tool package, and the qvalue package.
#' @param pVals Numeric vector of p-values to adjust.
#' @param adjust Allows for resulting p-values to be corrected for multiple hypothesis tests. Optional and some non-default choices require the "fdrtool" package or the "qvalue" package. Default = "none", which means that no p-value adjustment is performed. Other options include methods in ?p.adjust (i.e., "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), methods in ?fdrtool (i.e., "fndr", "pct0", "locfdr"), and methods in ?qvalue (i.e. "qvalue").
#' @param plotFdr Allows for plotting of fdrtool or qvalue p-value adjustment result, if this is chosen. Requires fdrtool package or qvalue package. Default = FALSE.
#' @param verbose If TRUE, the function prints out the general method used for multiple correction analysis. Default = FALSE.
#' @return A numeric vector of p-values (or q-values) that have been adjusted according to the specified method.
#' @examples
#' pvals = runif(100, 0, 1)
#' adj_pvals_bh = adjustPVals(pvals, adjust = "BH")
#' adj_pvals_hommel = adjustPVals(pvals, adjust = "hommel")
#' adj_pvals_qval = adjustPVals(pvals, adjust = "qvalue")
#' @export
adjustPVals <- function(pVals, adjust = "none", plotFdr = FALSE,
	verbose = FALSE){

	p_adjust_methods = c("holm", "hochberg", "hommel", "bonferroni", "BH",
		"BY", "fdr")
	fdrtool_methods = c("fndr", "pct0", "locfdr")
	qvalue_methods = c("qvalue")

	if(!(adjust %in% c("none", p_adjust_methods, fdrtool_methods))){
		stop("Adjust method is one of the available methods.")
	}

	pVals = as.numeric(pVals)

	if(adjust == "none"){
		if(verbose){
			message("No adjustment of p-values for multiple comparisons.")
		}
		pvals_adj = pVals
	}

	if(adjust %in% qvalue_methods){
		if(verbose){
			message(paste0("Adjusting p-values for multiple comparisons using ", adjust, " method using qvalue library."))
		}
		qobj = tryCatch(
		    {
		      qvalue::qvalue(p = pvalues, lambda = seq(0.05, 0.95, 0.05))
		    }, error=function(cond) {
		      message("Here's the original error message:")
		      message(cond)
		      cat("\n")
		      message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
		      return(NA)
		      #if the qvalue computation returned without error, then its format should be a list; if not, there was an error.
		      qobj = tryCatch(
		        {
		          qvalue::qvalue(p = pvalues, lambda = seq(0.1, 0.9, 0.05),cl=cl)
		        }, error=function(cond) {
		          message("Here's the original error message:")
		          message(cond)
		          cat("\n")
		          message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
		          return(NA)
		          qobj = tryCatch(
		            {
		              qvalue::qvalue(p = pvalues, lambda = seq(0.2, 0.8, 0.05),cl=cl)
		            }, error=function(cond) {
		              message("Here's the original error message:")
		              message(cond)
		              cat("\n")
		              message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
		              return(NA)
		              qobj = tryCatch(
		                {
		                  qvalue::qvalue(p = pvalues, lambda = seq(0.3, 0.7, 0.05),cl=cl)
		                }, error=function(cond) {
		                  message("Here's the original error message:")
		                  message(cond)
		                  cat("\n")
		                  message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence... if this doesn't work, will report the the adjusted q-values as NA values to indicate that q-value adjustment did not work.")
		                  qobj$qvalues = rep(NA, length(pvalues))
		                  return(qobj)
		              })
		          })
		      })
		    }) ## close tryCatch block

		pvals_adj = qobj$qvalues
		if (plotFdr){
			tryCatch(
			    {
				  pdf(file="qvalue_adjustment.pdf")
	              plot(qobj)
	              dev.off()
	            }, error=function(cond) {
	              message("Warning while plotting q-values. Here's the original error message:")
	              message(cond)
	              cat("\n")
	            })
		}
	} ## close qvalue_methods block

	if(adjust %in% p_adjust_methods){
		if(verbose){
			message(paste0("Adjusting p-values for multiple comparisons using ", adjust, " method in p.adjust."))
		}
		pvals_adj = p.adjust(pVals, method = adjust)
	}

	if(adjust %in% fdrtool_methods){
		if (!requireNamespace("fdrtool", quietly = TRUE)) {
			stop("fdrtool is needed for the p-value adjustment technique you chose to work. Please install it.",
				call. = FALSE)
		}
		if(verbose){
			message("Adjusting p-vals for multiple comparisons using fdrtool, and selecting the local FDR result.")
		}
		pvals_adj = suppressMessages(fdrtool::fdrtool(pVals,
			statistic = "pvalue", plot = plotFdr, cutoff.method = adjust, verbose = verbose))
		pvals_adj = pvals_adj$lfdr
	}

	return(pvals_adj)

}
