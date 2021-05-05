# Code to calculate pulled speciation rates
#########
# For Fig 3 of new conception Current Biology
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


cl = makeCluster(100, type = 'SOCK', outfile="")
registerDoSNOW(cl)

# start parallel loop
#ntrees=100
whichTrees<-c(1:100)#[-16]
foreach(i=whichTrees, .packages=c('castor','ape','phytools'), .verbose=TRUE) %dopar% {

	# load tree i
	mamPhy_wOut<-read.tree(file=paste0("final_100trees/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_100trees/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_",i,"_newick.tre"))
	mamPhy<-drop.tip(mamPhy_wOut,"_Anolis_carolinensis")

      # trim the first 1 My off the tree
        #get root
        totalRoot<-max(branching.times(mamPhy))
        #do slice
        mamPhy_sliced<-treeSlice(tree=mamPhy, slice=totalRoot-1, orientation="rootwards")

  # Get fossil-estimated extinction rates -- 5 Ma time bins
  #########
  #perCapRates_5Ma <- read.csv(file="fossilGeneraThroughTime_PBDB_perCapitaRates_bins5Ma.csv", header=TRUE)
  #     # 'impute' the missing rate data from 70-65 Ma, taking the mean of the before / after values
  #     perCapRates_5Ma[13,"pRate"]<-mean(c(perCapRates_5Ma[12,"pRate"],perCapRates_5Ma[14,"pRate"]))
  #     perCapRates_5Ma[13,"qRate"]<-mean(c(perCapRates_5Ma[12,"qRate"],perCapRates_5Ma[14,"qRate"]))
   mamDivDyn<-read.csv(file="mamDivDyn_ON_pbdb_Mammalia_wPiresEtAl2018_binned1to131_wHiLow.csv", header=TRUE)
   # get rate:
      #fossilExtinctionRates <- perCapRates_5Ma[26:4,"qRate"]
      fossilExtinctionRates_HIGH <- mamDivDyn[1:23,"highExt"]
      fossilExtinctionRates_LOW <- mamDivDyn[1:23,"lowExt"]

   # set resDir
   resDir<-"estHBD_fixedMu_100trees_highLow/"
   setwd(resDir)

   # SET root age
   #root<-round(max(branching.times(mamPhy)),0)
   #root<-150 # use 150 to divide the tree at the very root.
   root<-110 # use 111 since that's all that I'd display anyway

    res_HBD_HIGH<- fit_hbd_model_on_grid(tree = mamPhy, 
                           oldest_age        = root,
                           age0              = 0,
                           age_grid          = seq(0, root, by=5),
                           min_lambda        = 0,
                           max_lambda        = 10,
                           min_mu            = 0,
                           max_mu            = 10,
                           min_rho0          = 1e-10,
                           max_rho0          = 1,
                           guess_lambda      = NULL,
                           guess_mu          = NULL,
                           guess_rho0        = 1,
                           fixed_lambda      = NULL,
                           fixed_mu          = fossilExtinctionRates_HIGH,
                           fixed_rho0        = 1.0,
                           const_lambda      = FALSE,
                           const_mu          = FALSE,
                           splines_degree    = 1,
                           condition         = "auto",
                           relative_dt       = 1e-3,
                           Ntrials           = 1,
                           Nthreads          = 1,
                           max_model_runtime = NULL,
                           fit_control       = list())

	save(res_HBD_HIGH, file=paste0("estHBD_fixedMu_5911species_NDexp_tree",i,"_1trial_HIGH.Rda"))

  res_HBD_LOW<- fit_hbd_model_on_grid(tree = mamPhy, 
                           oldest_age        = root,
                           age0              = 0,
                           age_grid          = seq(0, root, by=5),
                           min_lambda        = 0,
                           max_lambda        = 10,
                           min_mu            = 0,
                           max_mu            = 10,
                           min_rho0          = 1e-10,
                           max_rho0          = 1,
                           guess_lambda      = NULL,
                           guess_mu          = NULL,
                           guess_rho0        = 1,
                           fixed_lambda      = NULL,
                           fixed_mu          = fossilExtinctionRates_LOW,
                           fixed_rho0        = 1.0,
                           const_lambda      = FALSE,
                           const_mu          = FALSE,
                           splines_degree    = 1,
                           condition         = "auto",
                           relative_dt       = 1e-3,
                           Ntrials           = 1,
                           Nthreads          = 1,
                           max_model_runtime = NULL,
                           fit_control       = list())

  save(res_HBD_LOW, file=paste0("estHBD_fixedMu_5911species_NDexp_tree",i,"_1trial_LOW.Rda"))


} # end 100 trees

stopCluster(cl)

q()

n
