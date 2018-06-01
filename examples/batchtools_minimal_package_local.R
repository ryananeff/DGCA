## module load R/3.5.0
rm(list = ls(all.names = TRUE))
set.seed(12345)

#.libPaths(new="~/.Rlib35DGCA")
library("DGCA")#,lib="~/.Rlib35DGCA")

options(stringsAsFactors = FALSE)

message("DGCA batch tools package! hooray!")
message("loading data")
message(Sys.time())

data(darmanis)
data(design_mat)

setwd("/Volumes/My Book/dgca_tn_er/")

######################

split = 2					#number of times to split the input. This means (split**2-split)/2 comparisons will be run (see below)
nPerms = 15																	#perms to run 
num_cores = 8																#cores per worker
outputfile = "dgca_tn_er_chunks_20split_15perm_7c_output.batchtools.txt"				#output file (TSV)
input_data = darmanis														#input dataframe
design_mat = design_mat 													#input design matrix
groups = c("oligodendrocyte", "neuron")														#groups			
corrType="spearman"							#correlation type (pearson, spearman)

walltime = 90 								#walltime in minutes
memory = 3800								#memory per core (MB)

######################

cl = parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

ddcorAll(input_data, design_mat, groups,
                 corrType="spearman",adjust="perm",
                 nPerms=nPerms,cl=cl)

ddcorAllParallel(input_data, design_mat, groups, outputfile,
                 sigOutput=TRUE,corrType="spearman",adjust="perm",
                 nPerms=nPerms,verbose=TRUE,perBatch=split,
                 coresPerJob=num_cores,testJob=F,
                 memPerJob = memory,timePerJob=walltime, batchDir="batchminimal",
                 batchConfig = "config/batchConfig_Local.R")
