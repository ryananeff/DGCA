#' @title Worker thread for ddcorAllParallel (internal function).
#' @description Runs the discovery of differential correlation (ddcor) section for comparing pairwise correlations across conditions in the Differential Gene Correlation Analysis (DGCA) package on a worker thread for two blocks of data (A and B).
#' @param job A batchtools job instance.
#' @param data A named list containing the program kwargs.
#' @param instance Required by batchtools, but not used currently.
#' @return Typically, the returned object is a data frame of the table of differential correlations between conditions. In the case that dCorAvg is calculated, the returned object is instead a list containing that table as well as the object summarizing the difference in average correlation for the specified portion of the data set.
#' @export
ddcorAllParallelWorker <- function(job,data,instance){

	# Required fields in data:
	# no_cores
	# batchWarningLevel
	# matA
	# matB
	# design_mat
	# groups
	# corrType
	# nPerms
	# verbose
	# seed
	print(data$libloc)
	library("DGCA",lib.loc=data$libloc)

	if(data$verbose){
		print(attributes(data))
	}
	
	options(warn=data$batchWarningLevel)

	cl<-parallel::makeCluster(data$n.cores)
	doParallel::registerDoParallel(cl)

	filepath = system.file("config/ddcorClasses.R",package="DGCA",lib.loc=data$libloc)
	parallel::clusterExport(cl=cl,c("filepath"),envir=environment())
	parallel::clusterEvalQ(cl=cl,eval(parse(filepath)))

	set.seed(data$seed) #random seed for reproducibility

	nPairs = (nrow(data$matA)*nrow(data$matA)-nrow(data$matA))/2+(nrow(data$matB)*nrow(data$matB)-nrow(data$matB))/2

	if(data$verbose){
		cat("Starting run now...\n")
		cat(paste0(Sys.time(),"\n"))
	}
	ddcor_res = ddcorAll(nPerms = data$nPerms, nPairs = nPairs, inputMat = data$matA, inputMatB=data$matB, design = data$design,
	                   compare = data$compare, cl=cl, corrType = data$corrType,empOnly=TRUE,classify=FALSE)
	#remove NAs caused by ??
	ddcor_res = ddcor_res[!is.na(ddcor_res$pValDiff),]
	ddcor_res = ddcor_res[!is.na(ddcor_res$empPVals),]
	ddcor_res[,"pValDiff_adj"] <- NULL #delete column we're not going to use anymore (b/c confusing name)

	#recalculate the q-values two ways
	if(data$verbose){
		message("calculating qvalues now.")
		cat(paste0(Sys.time(),"\n"))
	}
	ddcor_res[,"pValDiff_adj"] <- NULL
	ddcor_res[,"qValDiff"]=as.matrix(getQValue(ddcor_res$pValDiff)$qvalues)
	ddcor_res[,"qValDiff_emp"]=as.matrix(getQValue(ddcor_res$empPVals)$qvalues)
	classes = dCorClass(ddcor_res[,3],ddcor_res[,4],ddcor_res[,5],ddcor_res[,6],ddcor_res[,"qValDiff"],convertClasses=T,corSigThresh=0.05)
    ddcor_res[,"Classes_qval"] = classes
    classes = dCorClass(ddcor_res[,3],ddcor_res[,4],ddcor_res[,5],ddcor_res[,6],ddcor_res[,"qValDiff_emp"],convertClasses=T,corSigThresh=0.05)
    ddcor_res[,"Classes_qval_emp"] = classes
    if(data$verbose){
		cat("Completed run\n")
		cat(paste0(Sys.time(),"\n"))
	}
	return(ddcor_res) #algorithm
}