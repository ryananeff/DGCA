#' @title Calculate a correlation matrix.
#' @description This function takes one or two input matrices and calculates a correlation matrix from it using the speed-optimized correlation function from WGCNA.
#' @param matA Input data matrix with numeric entries.
#' @param corrType The type of correlation to be performed. Either "pearson" or "spearman".
#' @param matB Optional input data matrix with which the comparison with matA will be made.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @param use The "use" method for performing the correlation calculation. See ?cor for more information. Default = "pairwise.complete.obs" (which is one of the speed-optimized versions; see ?WGCNA::cor for more).
#' @param k If method="mutualinformation", the number of nearest neighbors to consider to estimate the mutual information. Must be less than the number of columns of mat.
#' @param noise If method="mutualinformation, the magnitude of the random noise added to break ties.
#' @return A correlation matrix.
#' data(darmanis); darmanis_subset = darmanis[1:30, ]
#' matcor_res = matCorr(matA = darmanis_subset, corrType = "pearson")
#' @export

matCorr <- function(matA, corrType, use = "pairwise.complete.obs", matB = NULL, secondMat = FALSE,
                    k=3,noise=1e-12){

	if(!secondMat){
		if(corrType %in% "pearson"){
			corrs = cor(matA,use=use,method="pearson")
		}
		if(corrType %in% "spearman"){
			corrs = cor(matA,use=use,method="spearman")
		}
		if(corrType %in% "mutualinformation"){
			#EXPERIMENTAL

			### OMP_NUM_THREADS controls the number of threads knnmi uses ###
			# we have to set it to 1 core because it may compete with other
			#corr calculations going on at the same time
			prev_num_cores = Sys.getenv("OMP_NUM_THREADS")
			Sys.setenv(OMP_NUM_THREADS=1) 

			corrs = parmigene::knnmi.all(matA,k=k,noise=noise)

			Sys.setenv(OMP_NUM_THREADS=prev_num_cores)
		}
	}

	if(secondMat){
		if(corrType %in% "pearson"){
			corrs = cor(matA,matB,use=use,method="pearson")
		}
		if(corrType %in% "spearman"){
		  corrs = cor(matA,matB,use=use,method="spearman")
		}
		if(corrType %in% "mutualinformation"){
			#EXPERIMENTAL

			### OMP_NUM_THREADS controls the number of threads knnmi uses ###
			# we have to set it to 1 core because it may compete with other
			#corr calculations going on at the same time
			prev_num_cores = Sys.getenv("OMP_NUM_THREADS")
			Sys.setenv(OMP_NUM_THREADS=1) 

			corrs = parmigene::knnmi.cross(matA,matB,k=k,noise=noise)

			Sys.setenv(OMP_NUM_THREADS=prev_num_cores)
		}
	}

	return(corrs)

}

###