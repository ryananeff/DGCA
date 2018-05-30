## For configuration instructions and examples, see the batchtools R package reference on CRAN. 

cluster.functions = makeClusterFunctionsLSF("inst/config/minerva_LSF.tmpl",scheduler.latency = 1)
default.resources = list(queue="premium", walltime=60, memory=2000, cores=1, project="acc_zhangb03a")
raise.warnings = FALSE
staged.queries = TRUE