#' @title Calculate q-values from DGCA class objects based on permutation-based empirical null statistics.
#' @description First, estimate empirical p-values based on a comparison of the actual and permuted test statistics. Next, estimate the proportion of true null hypotheses using the qvalue package as well as qvalues from the empirical p-values, using this value. If the estimated pi0 <= 0, then sequentially recalculates using increasingly conservative set of lambda values, until lambda = 0.5.
#' @param dcObject The original S4 class object containing the test statistics to be extracted.
#' @param permObject The array of matrices containing the null test statistics.
#' @param secondMat Logical, indicating whether a second matrix was used in the construction of this dcObject and permObject. If FALSE, the upper.tri of both are extracted to avoid double counting test statistics.
#' @param testSlot The slot of the dcObject to be removed for use as the actual test statistic.
#' @param verbose Whether summaries of the q-value operations should be reported.
#' @param plotFdr Allows for plotting of fdrtool p-value adjustment result OR empirical FDR q-value adjustment technique, if either of these are chosen. Requires fdrtool package OR qvalue package. Default = FALSE.
#' @param empOnly Whether or not we don't want to calculate qvalues for the empirical pvalues. Default = FALSE.
#' @return A list containing a vectof of empirical p-values and a vector of q-values, both of the same length as the original actual test statistics.
#' @export
permQValue <- function(dcObject, permObject, secondMat, testSlot,
  verbose = FALSE, plotFdr = FALSE,empOnly=FALSE){

  test_stat_actual = slot(dcObject, testSlot)
  if(!secondMat) test_stat_actual = test_stat_actual[upper.tri(test_stat_actual)]
  test_stat_actual = as.numeric(test_stat_actual)

  if(!secondMat) permObject = permObject[apply(permObject, 3, upper.tri)]
  perm_stats = as.numeric(permObject)

  message("Calculating empirical p-values using the permutation sample statistics.")

  pvalues = bigEmpPVals(stat = abs(test_stat_actual), stat0 = abs(perm_stats))

  message("Calculating qvalues from the empirical p-values.")
  if(!empOnly){

    qobj = getQValue(pvalues)

    if(verbose) summary(qobj)
    if (plotFdr){
      tryCatch(
          {
                grDevices::pdf(file="qvalue_adjustment.pdf")
                graphics::plot(qobj)
                grDevices::dev.off()
              }, error=function(cond) {
                message("Warning while plotting q-values. Here's the original error message:")
                message(cond)
                cat("\n")
              })
    }
    return(list(empPVals = pvalues, pValDiffAdj = qobj$qvalues))
  } else {
      message("Skipping q-value calculation step (perhaps to calculate it at a later step?)")
      return(list(empPVals = pvalues, pValDiffAdj = rep(NA,length.out=length(pvalues))))
  }

}