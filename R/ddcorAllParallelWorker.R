#' @title Worker thread for ddcorAllParallel (internal function).
#' @description Runs the discovery of differential correlation (ddcor) section for comparing pairwise correlations across conditions in the Differential Gene Correlation Analysis (DGCA) package on a worker thread for two blocks of data (A and B).
#' @param job A batchtools job instance.
#' @param data A named list containing the program kwargs.
#' @param instance Required by batchtools, but not used currently.
#' @return Typically, the returned object is a data frame of the table of differential correlations between conditions. In the case that dCorAvg is calculated, the returned object is instead a list containing that table as well as the object summarizing the difference in average correlation for the specified portion of the data set.
#' @export
ddcorAllParallelWorker <- function(job,data,instance){

	# Required inputs:
	# no_cores
	# batchWarningsAsErrors
	# matA
	# matB
	# design_mat
	# groups
	# corrType
	# nPerms
	# verbose
	# seed

	for(i in names(data)){assign(i,data[[i]],pos=1)} # assign the data to variables on the worker

	if(verbose){
		print(attributes(data))
	}

	
	options(warn=batchWarningLevel)

	library(methods)
	library(foreach)
	library(doParallel)
	library(qvalue)
	library(DGCA)

	get_qvalues <- function(pvalues){
	  try_lambda = 20
	  qobj = tryCatch(
	      {
            if ( max(0.01,round(min(pvalues)+0.01,2)) >= min(round(max(pvalues)-0.01,2),0.95) ) {
                 qobj = list()
                 qobj$qvalues = rep(NA, length.out=length(pvalues))
                 return(qobj)
            }else{
            rangevals = (max(pvalues)-min(pvalues))
	        qvalue::qvalue(p = pvalues, lambda = seq(max(0.01,round(min(pvalues)+0.01,2)), 
	                                                 min(round(max(pvalues)-0.01,2),0.95), 
	                                                 min(0.01,round(rangevals/try_lambda,4))), lfdr.out=FALSE)
            }
	      }, error=function(cond) {
	        message("Here's the original error message:")
	        #message(cond)
	        cat("\n")
	        message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
	        #if the qvalue computation returned without error, then its format should be a list; if not, there was an error.
	        qobj = tryCatch(
	          {
	            qvalue::qvalue(p = pvalues, lambda = seq(0.1, 0.9, 0.01),lfdr.out=FALSE)
	          }, error=function(cond) {
	            message("Here's the original error message:")
	            #message(cond)
	            cat("\n")
	            message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
	            qobj = tryCatch(
	              {
	                qvalue::qvalue(p = pvalues, lambda = seq(0.2, 0.8, 0.01),lfdr.out=FALSE)
	              }, error=function(cond) {
	                message("Here's the original error message:")
	                #message(cond)
	                cat("\n")
	                message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
	                qobj = tryCatch(
	                  {
	                    qvalue::qvalue(p = pvalues, lambda = seq(0.3, 0.7, 0.01),lfdr.out=FALSE)
	                  }, error=function(cond) {
	                    message("Here's the original error message:")
	                    #message(cond)
	                    cat("\n")
	                    message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence... if this doesn't work, will report the empirical p-values and the adjusted q-values as NA values to indicate that q-value adjustment did not work.")
	                    qobj = list()
	                    qobj$qvalues = rep(NA, length.out=length(pvalues))
	                    return(qobj)
	                })
	            })
	        })
	      })
	  return(as.matrix(qobj$qvalues))
	}

	cl<-makeCluster(n.cores)
	registerDoParallel(cl)

	clusterEvalQ(cl=cl,eval(parse(file.path(path.package("DGCA"), "R/ddcorClasses.R"))))

	set.seed(seed) #random seed for reproducibility

	nPairs = (nrow(matA)*nrow(matA)-nrow(matA))/2+(nrow(matB)*nrow(matB)-nrow(matB))/2

	if(verbose){
		cat("Starting run now...\n")
		cat(paste0(Sys.time(),"\n"))
	}
	ddcor_res = ddcorAll(nPerms = nPerms, nPairs = nPairs, inputMat = matA, inputMatB=matB, design = design_mat,
	                   compare = groups, cl=cl, corrType = corrType,empOnly=TRUE,classify=FALSE)
	#remove NAs caused by ??
	ddcor_res = ddcor_res[!is.na(ddcor_res$pValDiff),]
	ddcor_res = ddcor_res[!is.na(ddcor_res$empPVals),]
	ddcor_res[,"pValDiff_adj"] <- NULL #delete column we're not going to use anymore (b/c confusing name)

	#recalculate the q-values two ways
	if(verbose){
		message("calculating qvalues now.")
		cat(paste0(Sys.time(),"\n"))
	}
	ddcor_res[,"pValDiff_adj"] <- NULL
	ddcor_res[,"qValDiff"]=get_qvalues(ddcor_res$pValDiff)
	ddcor_res[,"qValDiff_emp"]=get_qvalues(ddcor_res$empPVals)
	classes = dCorClass(ddcor_res[,3],ddcor_res[,4],ddcor_res[,5],ddcor_res[,6],ddcor_res[,"qValDiff"],convertClasses=T,corSigThresh=0.05)
    ddcor_res[,"Classes_qval"] = classes
    classes = dCorClass(ddcor_res[,3],ddcor_res[,4],ddcor_res[,5],ddcor_res[,6],ddcor_res[,"qValDiff_emp"],convertClasses=T,corSigThresh=0.05)
    ddcor_res[,"Classes_qval_emp"] = classes
    if(verbose){
		cat("Completed run\n")
		cat(paste0(Sys.time(),"\n"))
	}
	return(ddcor_res) #algorithm
}