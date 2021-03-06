#' @title Calculate a correlation matrix.
#' @description This function takes one or two input matrices and calculates a correlation matrix from it using the speed-optimized correlation function from WGCNA.
#' @param matA Input data matrix with numeric entries.
#' @param corrType The type of correlation to be performed. Either "pearson" or "spearman".
#' @param matB Optional input data matrix with which the comparison with matA will be made.
#' @param secondMat Logical indicator of whether there is a second matrix in the comparison or not.
#' @param use The "use" method for performing the correlation calculation. See ?cor for more information. Default = "pairwise.complete.obs" (which is one of the speed-optimized versions; see ?WGCNA::cor for more).
#' @param k If method="mutualinformation", the number of nearest neighbors to consider to estimate the mutual information. Must be less than the number of columns of mat.
#' @param noise If method="mutualinformation, the magnitude of the random noise added to break ties.
#' @param k When running in MI mode, the number of intervals to discretize the data into before calculating mutual information. Default = 5.
#' @param k_iter_max When running in MI mode, the number of iterations to determine the k-clusters for discretization before calculating mutual information. Default = 10. 
#' @return A correlation matrix.
#' data(darmanis); darmanis_subset = darmanis[1:30, ]
#' matcor_res = matCorr(matA = darmanis_subset, corrType = "pearson")
#' @export

matCorr <- function(matA, corrType, use = "pairwise.complete.obs", matB = NULL, secondMat = FALSE,
                    k=5,k_iter_max=10,lib.loc=NULL){
	
	if (is.null(lib.loc)==FALSE){ #why
		.libPaths(lib.loc)
		library("DGCA",lib.loc=lib.loc)
		library("arules",lib.loc=lib.loc)
		print("loaded special libraries")
	}
	
	if(!secondMat){
		if(corrType %in% "pearson"){
			corrs = tryCatch({
                WGCNA::cor(matA, use = use, method="pearson")
            }, error = function(e){
                message("Large matrix detected, running slower 64-bit cor()")
                cor(matA, use=use, method="pearson")
            })
		}
		if(corrType %in% "spearman"){

			#speed hack 10/23/18
			colnames_matA = colnames(matA)
			matA = t(apply(matA,1,function(x) rank(x,ties="random")))
			colnames(matA) = colnames_matA

			corrs = tryCatch({
                WGCNA::cor(matA, use = use, method="pearson")
            }, error = function(e){
                message("Large matrix detected, running slower 64-bit cor()")
                cor(matA, use=use, method="pearson")
            })

		}
		if(corrType %in% "mutualinformation"){
			##EXPERIMENTAL - assume already discrete
			corrs = infotheo::mutinformation(data.frame(matA, stringsAsFactors=TRUE)) ##assume it is already a discrete matrix
			
		}
		if(corrType %in% "bicor"){
			#EXPERIMENTAL
			corrs = WGCNA::bicor(x=matA,use=use,quick=1)
			
		}
	}

	if(secondMat){
		if(corrType %in% "pearson"){
           corrs = tryCatch({
                WGCNA::cor(matA,matB,use=use,method="pearson")
            }, error = function(e){
                message("Large matrix detected, running slower 64-bit cor()")
                cor(matA,matB,use=use,method="pearson")
            })

		}
		if(corrType %in% "spearman"){

			#speed hack 10/23/18
			matA = t(matA)
			matA = apply(matA,1,function(x) rank(x,ties="random"))
			
			matB = t(matB)
			matB = apply(matB,1,function(x) rank(x,ties="random"))
            
            corrs = tryCatch({
                WGCNA::cor(matA,matB,use=use,method="pearson")
            }, error = function(e){
                message("Large matrix detected, running slower 64-bit cor()")
                cor(matA,matB,use=use,method="pearson")
            })


		}
		if(corrType %in% "mutualinformation"){
			#EXPERIMENTAL - assume already discrete
			df1 = data.frame(matA, stringsAsFactors=TRUE) ## WAO, there is a use for this after all...
			df2 = data.frame(matB, stringsAsFactors=TRUE) ## WAOWAO

			mi_df_2mat = matrix(data=NA,nrow=ncol(df1),ncol=ncol(df2))
			rownames(mi_df_2mat) = colnames(df1)
			colnames(mi_df_2mat) = colnames(df2)

			df1_out <<-df1
			df2_out <<-df2
			mi_df_2mat_out <<- mi_df_2mat

			for(col1 in 1:length(colnames(df1))){
			  for(col2 in col1:length(colnames(df2))){
			    mi_df_2mat[col1,col2] = infotheo::mutinformation(df1[,col1],df2[,col2])
			    #mi_df_2mat[col2,col1] = mi_df_2mat[col1,col2]
			  }
			}

			corrs = mi_df_2mat
			
		}
		if(corrType %in% "bicor"){
			#EXPERIMENTAL
			corrs = WGCNA::bicor(x=matA,y=matB,use=use,quick=1)
			
		}
	}

	return(corrs)

}

###
