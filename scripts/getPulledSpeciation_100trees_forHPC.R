# Code to calculate pulled speciation rates
#########
# For Fig 5 of new conception Current Biology
#######
# NS Upham; last modification: 27 Mar 2021
# ===================================================


# load packages
library(ape); library(castor); library(phytools); 
library(foreach);library(doSNOW)

# directory and source
dirname = "/gpfs/loomis/home.grace/nu35/"
#dirname = "/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/"
setwd(dirname)

cl = makeCluster(40, type = 'SOCK', outfile="")
registerDoSNOW(cl)

# start parallel loop
ntrees=100
foreach(i=81:ntrees, .packages=c('castor','ape','phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

	# load tree i
	mamPhy_wOut<-read.tree(file=paste0("final_100trees/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_100trees/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_",i,"_newick.tre"))
	mamPhy<-drop.tip(mamPhy_wOut,"_Anolis_carolinensis")

      # trim the first 1 My off the tree
        #get root
        totalRoot<-max(branching.times(mamPhy))
        #do slice
        mamPhy_sliced<-treeSlice(tree=mamPhy, slice=totalRoot-1, orientation="rootwards")

	# set resDir
	resDir<-"pulledSpeciation_100trees_new/"
	setwd(resDir)

	# get root age
	#root<-round(max(branching.times(mamPhy)),0)
	root<-110 # use 150 to divide the tree at the very root.

	res_PSR<-fit_hbd_psr_on_grid(  tree = mamPhy_sliced, 
                           oldest_age            = root,
                           age0                  = 0,
                           age_grid              = seq(0, root, by=5),
                           min_PSR               = 0,
                           max_PSR               = 10,
                           guess_PSR             = NULL,
                           fixed_PSR             = NULL,
                           splines_degree        = 1,
                           condition             = "auto",
                           relative_dt           = 1e-3,
                           Ntrials               = 1,
                           Nbootstraps           = 10,
                           Ntrials_per_bootstrap = NULL,
                           Nthreads              = 1,
                           max_model_runtime     = NULL,
                           fit_control           = list(),
                           verbose               = FALSE,
                           verbose_prefix        = "")

	save(res_PSR, file=paste0("PSR_5911species_NDexp_tree",i,"_new.Rda"))

#	res_PDR<-fit_hbd_pdr_on_grid(  tree = mamPhy, 
#                           oldest_age            = root,
#                           age0                  = 0,
#                           age_grid              = seq(0, root, by=5),
#                           min_PDR               = -100,
#                           max_PDR               = +100,
#                           min_rholambda0        = 1e-10,
#                           max_rholambda0        = +100,
#                           guess_PDR             = NULL,
#                           guess_rholambda0      = NULL,
#                           fixed_PDR             = NULL,
#                           fixed_rholambda0      = NULL,
#                           splines_degree        = 1,
#                           condition             = "auto",
#                           relative_dt           = 1e-3,
#                           Ntrials               = 1,
#                           Nbootstraps           = 10,
#                           Ntrials_per_bootstrap = NULL,
#                           Nthreads              = 1,
#                           max_model_runtime     = NULL,
#                           fit_control           = list(),
#                           verbose               = FALSE,
#                           verbose_prefix        = "")
#
#	save(res_PDR, file=paste0("PDR_5911species_NDexp_tree",i,".Rda"))

} # end 100 trees

stopCluster(cl)

q()

n
