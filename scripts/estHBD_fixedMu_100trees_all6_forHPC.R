# Code to calculate pulled speciation rates
#########
# For Fig 3 of new conception Current Biology
#######
# NS Upham; last modification: 7 Apr 2021
# ===================================================

#srun --pty -t 4:00:00 --mem=10G -p interactive bash
#module load  R/3.6.1-foss-2018b-X11-20180604

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
#ntrees=100
whichTrees<-c(21:60)#c(1:100)#[-16]
foreach(i=whichTrees, .packages=c('castor','ape','phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

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
      #  fossilExtinctionRates_HIGH <- mamDivDyn[1:23,"highExt"]
      #  fossilExtinctionRates_LOW <- mamDivDyn[1:23,"lowExt"]
   ###
   # Make the 71-66 Ma time bin for all extRates == 1.0 from the E2f3 rate (aiming to make pushed rate estimation possible for all metrics)
   ####
   extinctionRateNames<-c('extPC', 'ext3t', 'extC3t', 'extGF', 'E2f3', 'ext2f3')
   
   mamDivDyn[14,extinctionRateNames]<-c(1,1,1,1,1,1)

   # use lambda rates as initial guesses
    originationRateNames<-c('oriPC', 'ori3t', 'oriC3t', 'oriGF', 'O2f3', 'ori2f3')


   # set resDir
   resDir<-"estHBD_fixedMu_100trees_all6/"
   setwd(resDir)

   # SET root age
   #root<-round(max(branching.times(mamPhy)),0)
   #root<-150 # use 150 to divide the tree at the very root.
   root<-110 # use 110 since the tree was sliced== 111 in reality

    extinctionRateNames<-c('extPC', 'ext3t', 'extC3t', 'extGF', 'E2f3', 'ext2f3')
    for(j in 1:length(extinctionRateNames)){
  
    fossilExtinctionRates_j<-mamDivDyn[1:23,extinctionRateNames[j]]
    fossilOriginationRates_j<-mamDivDyn[1:23,originationRateNames[j]]

    res_HBD<- fit_hbd_model_on_grid(tree = mamPhy_sliced, 
                           oldest_age        = root,
                           age0              = 0,
                           age_grid          = seq(0, root, by=5),
                           min_lambda        = 0,
                           max_lambda        = +Inf, # omitting upper bounds on lambda
                           min_mu            = 0,    # mu is fixed
                           max_mu            = +Inf, # mu is fixed
                           min_rho0          = 1e-10,
                           max_rho0          = 1,
                           guess_lambda      = fossilOriginationRates_j, # using fossil origination rates as initial guesses for subsequent fitting
                           guess_mu          = NULL, # mu is fixed
                           guess_rho0        = 0.80, # rho is guessed from 0.80: assumption of 0.20 species/Ma median extinction rate toward the present
                           fixed_lambda      = NULL, # lambda is estimated at all points from guess starting point
                           fixed_mu          = fossilExtinctionRates_j, # mu is fixed
                           fixed_rho0        = NULL, # rho is estimated from guess starting point
                           const_lambda      = FALSE,
                           const_mu          = FALSE,
                           splines_degree    = 1,
                           condition         = "auto",
                           relative_dt       = 1e-3,
                           Ntrials           = 10,
                           Nthreads          = 1,
                           max_model_runtime = 1,
                           fit_control       = list())

	save(res_HBD, file=paste0("estHBD_fixedMu_5911species_NDexp_tree",i,"_1trial_all6_",extinctionRateNames[j],".Rda"))

  } # end 6 extinction rates

} # end 100 trees

stopCluster(cl)

q()

n
