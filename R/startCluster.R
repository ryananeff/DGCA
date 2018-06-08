#' @title DGCA cluster setup function.
#' @description Creates a cluster object for multithreaded single-node execution of the the DGCA pipeline.
#' @param n.cores The number of cores to use for execution. Default = 2.
#' @return A object of type parallel::cluster that is a cluster instance on the local machine with the number of cores specified.
#' @examples
#' n.cores = 8
#' cl = startCluster(n.cores)
#' @export 

startCluster <- function(n.cores=2){
	if(n.cores<2){
		message("Running in single-threaded mode (n.cores < 2).")
		return(FALSE) #if only one core (or less?) is specified
	}
	cl<-parallel::makeCluster(n.cores)
	doParallel::registerDoParallel(cl)
	filepath = system.file("config/ddcorClasses.R",package="DGCA")
	parallel::clusterExport(cl=cl,c("filepath"),envir=environment())
	parallel::clusterEvalQ(cl=cl,eval(parse(filepath)))
	return(cl)
}
