rm(list = ls(all.names = TRUE))
set.seed(12345)

library("DGCA",lib="~/.Rlib")

options(stringsAsFactors = FALSE)

message("loading data")
message(Sys.time())

data(darmanis) #loads daramanis
data(design_mat) #loads design_mat

######################

split = 5									#number of times to split the input. This means split**2 comparisons will be run (see below)
nPerms = 15									#perms to run
outputfile = "output.test.txt" 				#output file (TSV)
input_data = darmanis 						#input dataframe
design_mat = design_mat 					#input design matrix
groups = c("oligodendrocyte", "neuron")		#groups
num_cores = 2								#cores per worker
corrType="pearson"							#correlation type (pearson, spearman)

######################

ddcorAllParallel(input_data, design_mat, groups, outputfile,
                 sigOutput=TRUE,corrType=corrType,adjust="perm",
                 nPerms=nPerms,verbose=TRUE,perBatch=split,
                 coresPerJob=num_cores,batchDir="minimal",testJob=TRUE,
                 memPerJob = 10000,timePerJob=60)
