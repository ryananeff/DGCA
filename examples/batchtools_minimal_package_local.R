## module load R/3.4.3
##install
if(!dir.exists("~/.RlibDGCA")){dir.create("~/.RlibDGCA")}
source("https://bioconductor.org/biocLite.R")
biocLite(suppressUpdates=T)
withr::with_libpaths("~/.RlibDGCA",devtools::install_github("nosarcasm/DGCA"))

##load
.libPaths(new="~/.RlibDGCA") ## create this directory and install DGCA to it
library("DGCA",lib="~/.RlibDGCA")
options(stringsAsFactors = FALSE)

message("DGCA batch tools package! hooray!")
message("loading data")
message(Sys.time())

data(darmanis)
data(design_mat)

######################

set.seed(12345)
split = 2					#number of times to split the input. This means (split**2-split)/2 comparisons will be run (see below)
nPerms = 15																	#permutations to run 
num_cores = 8																#cores per worker
outputfile = "dgca_test_2split_15perm_8c_output.batchtools.txt"				#output file (TSV)
input_data = darmanis														#input dataframe
design_mat = design_mat 													#input design matrix
groups = c("oligodendrocyte", "neuron")														#groups			
corrType="spearman"							#correlation type (pearson, spearman, bicor, mutualinformation)

walltime = 90 								#walltime in minutes for batch jobs
memory = 3800								#memory per core (MB) for batch jobs

######################

##for multithreaded operation on single machine
cl = parallel::makeCluster(num_cores) 
doParallel::registerDoParallel(cl)

##single machine mode - multithreaded
ddcor_res = ddcorAll(input_data, design_mat, groups,
                 corrType="spearman",adjust="perm",
                 nPerms=nPerms,cl=cl)

## for batch parallel submission mode - writes output file automatically
ddcorAllParallel(input_data, design_mat, groups, outputfile,
                 sigOutput=TRUE,corrType="spearman",adjust="perm",
                 nPerms=nPerms,verbose=TRUE,perBatch=split,
                 coresPerJob=num_cores,testJob=F,
                 memPerJob = memory,timePerJob=walltime, batchDir="batchminimal",
                 batchConfig = "config/batchConfig_Local.R")
