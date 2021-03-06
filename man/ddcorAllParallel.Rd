% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ddcorAllParallel.R
\name{ddcorAllParallel}
\alias{ddcorAllParallel}
\title{Calls the DGCA pairwise pipeline, splitting input matrix into multiple batch jobs on an HPC.}
\usage{
ddcorAllParallel(inputMat, design, compare, outputFile, sigOutput = FALSE,
  impute = FALSE, corrType = "pearson", nPairs = Inf,
  sortBy = "zScoreDiff", adjust = "perm", nPerms = 10, classify = TRUE,
  sigThresh = 1, corSigThresh = 0.05, corPower = 2, verbose = FALSE,
  corr_cutoff = 0.99, getDCorAvg = FALSE, dCorAvgType = "gene_average",
  dCorAvgMethod = "median", signType = "none", oneSidedPVal = FALSE,
  perBatch = 10, coresPerJob = 2, timePerJob = 60, memPerJob = 2000,
  batchConfig = system.file("inst/config/batchConfig_Zhang.R"),
  batchDir = "batchRegistry", batchWarningLevel = 0, batchSeed = 12345,
  maxRetries = 3, testJob = FALSE)
}
\arguments{
\item{inputMat}{The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your compare.}

\item{design}{A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.}

\item{compare}{Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.}

\item{outputFile}{Location to save the output. This is necessary when the number of comparisons >= 2^31-1. Required.}

\item{sigOutput}{Should we save the significant results in a separate file? Default = FALSE.}

\item{impute}{A binary variable specifying whether values should be imputed if there are missing values. Note that the imputation is performed in the full input matrix (i.e., prior to subsetting) and uses k-nearest neighbors.}

\item{corrType}{The correlation type of the analysis, limited to "pearson" or "spearman". Default = "pearson".}

\item{nPairs}{Either a number, specifying the number of top differentially correlated identifier pairs to display in the resulting table, or a the string "all" specifying that all of the pairs should be returned. If splitSet is specified, this is reset to the number of non-splitSet identifiers in the input matrix, and therefore will not be evaluated.}

\item{sortBy}{Character string specifying the way by which you'd like to sort the resulting table. This will happen at the system shell level (requires a UNIX shell).}

\item{adjust}{Allows for resulting p-values to be corrected for multiple hypothesis tests, optional. Some non-default choices require the "fdrtool" package or the "qvalue". Default = "none", which means that no p-value adjustment is performed. Other options include "perm" to use permutation samples, methods in ?p.adjust (i.e., "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), and methods in ?fdrtool (i.e., "fndr", "pct0", "locfdr").}

\item{nPerms}{Number of permutations to generate. If NULL, permutation testing will not be performed. Default = "10".}

\item{classify}{Binary value specifying whether the correlation values in each condition and differential correlation scores should be used to classifying the resulting identifiers into compare. Default = TRUE}

\item{sigThresh}{If classify = TRUE, this numeric value specifies the p-value threshold at which a differential correlation p-value is deemed significant for differential correlation class calculation. Default = 1, as investigators may use different cutoff thresholds; however, this can be lowered to establish significant classes as desired.}

\item{corSigThresh}{If classify = TRUE, this numeric value specifies the p-value threshold at which a correlation p-value is deemed significant. Default = 0.05.}

\item{corPower}{The power to raise the correlations to before plotting the classic heatmap. Larger correlation powers emphasize larger correlation values relatively more compared to smaller correlation values.}

\item{verbose}{Option indicating whether the program should give more frequent updates about its operations. Default = FALSE.}

\item{corr_cutoff}{Cutoff specifying correlation values beyond which will be truncated to this value, to reduce the effect of outlier correlation values when using small sample sizes. Note that this does NOT affect the underlying correlation values, but does affect the z-score difference of correlation calculation in the dcTopPairs table. Default = 0.99}

\item{getDCorAvg}{Logical, specifying whether the average difference in correlation between compare should be calculated. Default = FALSE}

\item{dCorAvgType}{Character vector specifying the type of average differential correlation calculation that should be performed. Only evaluated if dCorAge is TRUE. Types = c("gene_average", "total_average", "both"). gene_average calculates whether each genes' differential correlation with all others is more than expected via permutation samples (and empirical FDR adjustment, in the case of > 1 gene), while total_average calculates whether the total average differential correlation is higher than expected via permutation samples. "both" performs both of these. If splitSet is specified, then only genes in the splitSet have their average gene differential correlation calculated if gene_average is chosen.}

\item{dCorAvgMethod}{Character vector specifying the method for calculating the "average" differential correlation calculation that should be used. Options = "median", "mean".}

\item{signType}{Coerce all correlation coefficients to be either positive (via "positive"), negative (via "negative"), or none (via "none") prior to calculating differential correlation. This could be used if, e.g., you think that going from a positive to a negative correlation is unlikely to occur biologically and is more likely to be due to noise, and you want to ignore these effects. Note that this does NOT affect the reported underlying correlation values, but does affect the z-score difference of correlation calculation. Default = "none", for no coercing.}

\item{oneSidedPVal}{If the dCorAvgType test is total_average, this option specifies whether a one-sided p-value should be reported, as opposed to a two-sided p-value. That is, if the average difference of z-scores is greater than zero, test whether the permutation average difference of z-scores are less than that average to get the p-value, and vice versa for the case that the average difference of z-scores is less than 0. Otherwise, test whether the absolute value of the average difference in z-scores is greater than the absolute values of the permutation average difference in z-scores. Default = FALSE.}

\item{perBatch}{Number of times to split the features of the input data into separate batches. A higher number creates a larger number of jobs, but may be less uniform. Default = 10.}

\item{coresPerJob}{Number of cores to use on each batch job run. Default = 2.}

\item{timePerJob}{Walltime to request for each batch job (e.g. in a HPC cluster), in minutes. Default = 60}

\item{memPerJob}{Memory to request for each batch job (e.g. in a HPC cluster), in MB. Default = 2000}

\item{batchConfig}{Location of the batchtools configuration file (e.g. to configure this tool to work with your HPC cluster). Defaults to one used at inst/config/batchConfig_Zhang.R.}

\item{batchDir}{Location to store temporary files, logs, and results of the batch run. This is the registry for the batchtools R package. Default = batchRegistry/}

\item{batchWarningLevel}{Warning level on remote nodes during DGCA calculation (equivalent to setting options(warn=batchWarningLevel). Default = 0.}

\item{batchSeed}{Random seed to use on all batch jobs. Default = 12345.}

\item{maxRetries}{Number of times to re-submit jobs that failed. This is helpful for jobs that failed due to transient errors on an HPC. Default = 3}

\item{testJob}{Test one job before running it? Default = FALSE}
}
\value{
Returns whether all jobs successfully executed or not. Output is in the output file.
}
\description{
Runs the full discovery of differential correlation (ddcor) section for comparing pairwise correlations across conditions in the Differential Gene Correlation Analysis (DGCA) package.
}
