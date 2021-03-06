% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dcTopPairs.R
\name{dcTopPairs}
\alias{dcTopPairs}
\title{Creates a data frame for the top differentially correlated gene pairs in your data set.}
\usage{
dcTopPairs(dcObject, nPairs, adjust = "none", plotFdr = FALSE,
  classify = TRUE, sigThresh = 1, corSigThresh = 0.05,
  zScorePerm = NULL, verbose = FALSE, compare = NULL, secondMat = FALSE,
  empOnly = FALSE)
}
\arguments{
\item{dcObject}{The dcPair class object which you'd like to convert into a table.}

\item{nPairs}{The number of gene pairs to display in the resulting table.}

\item{adjust}{Allows for resulting p-values to be corrected for multiple hypothesis tests, optional. Some non-default choices require the "fdrtool" package or the "qvalue". Default = "none", which means that no p-value adjustment is performed. Other options include "perm" to use permutation samples, methods in ?p.adjust (i.e., "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), and methods in ?fdrtool (i.e., "fndr", "pct0", "locfdr").}

\item{plotFdr}{Allows for plotting of p-value adjustment result, if this is chosen. Requires fdrtool or qvalue package. Default = FALSE.}

\item{classify}{Binary value specifying whether the correlation values in each condition and differential correlation scores should be used to classifying the resulting identifiers into groups. Default = TRUE}

\item{sigThresh}{If classify = TRUE, this numeric value specifies the p-value threshold at which a differential correlation p-value is deemed significant for differential correlation class calculation. Default = 1, as investigators may use different cutoff thresholds; however, this can be lowered to establish significant classes as desired.}

\item{corSigThresh}{If classify = TRUE, this numeric value specifies the p-value threshold at which a correlation p-value is deemed significant. Default = 0.05.}

\item{zScorePerm}{A matrix of values with z-scores from permutation tests to be used to generate empirical p-values. Default = NULL.}

\item{verbose}{Whether summaries of the operations should be reported.}

\item{compare}{Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.}

\item{secondMat}{Logical indicator of whether there is a second matrix in the comparison or not.}

\item{empOnly}{Whether or not we don't want to calculate qvalues for the empirical pvalues. Default = FALSE.}
}
\value{
A table containing columns for each name in the considered gene pair (the order of which is arbitrary), correlation values in each condition, differences in z-score of the correlation, and p-values for that z-score difference.
}
\description{
Reads in a dcPair object and outputs a table of all gene pairs (or just the top n pairs), sorted by their unadjusted differential correlation p-value.
}
