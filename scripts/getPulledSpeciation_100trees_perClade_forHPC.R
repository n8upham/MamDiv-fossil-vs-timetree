# Code to calculate pulled speciation rates
#########
# For Fig 5 of new conception RSPB
#######
# NS Upham; last modification: 18 May 2020
# ===================================================


# load packages
library(ape); library(castor) #library(phangorn); 
library(foreach);library(doSNOW)

# directory and source
dirname = "/gpfs/loomis/home.grace/nu35/"
#dirname = "/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/"
setwd(dirname)

cl = makeCluster(100, type = 'SOCK', outfile="")
registerDoSNOW(cl)

# start parallel loop
ntrees=100
foreach(i=1:ntrees, .packages=c('castor','ape'), .verbose=TRUE) %dopar% {

	# load tree i
	mamPhy_wOut<-read.tree(file=paste0("final_100trees/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_100trees/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_",i,"_newick.tre"))
	mamPhy<-drop.tip(mamPhy_wOut,"_Anolis_carolinensis")
	
	# load taxonomy
	mamTax<-read.csv(file=paste0("final_100trees/taxonomy_mamPhy_5911species.csv"), header=TRUE)

	# get clades
		# get clade names
		cladeNames_all<-names(table(mamTax$clade))
		cladeNames<-c("Didelphimorphia", "Australidelphia",
					"Atlantogenata","Talpidae", "Erinaceidae", #<< exclude 24 species
					"Soricidae",
					"Feliformes","Caniformes",
					"Suina", #<< exclude 21 species
					"Whippomorpha","Ruminantia",
					"Yinpterochiroptera","Yangochiroptera",
					"Strepsirrhini", "Catarrhini", "Platyrrhini", 
					"Lagomorpha", "Guinea_pig-related", "Squirrel-related", "Mouse-related")

		australidelphiaNames<-c("Diprodontia", "Dasyuromorphia","Peramelemorphia", "Notoryctemorphia", "Microbiotheria")
		atlantogenataNames<-c("Afrotheria","Xenarthra")

		# get clade species
		cladeTipNames<-list()
		for(j in 1:length(cladeNames)){
			if(j==2){
				cladeTipNames[[j]]<-as.character(mamTax[which( 	mamTax$clade==australidelphiaNames[1] | 
												 	mamTax$clade==australidelphiaNames[2] |
												 	mamTax$clade==australidelphiaNames[3] |
												 	mamTax$clade==australidelphiaNames[4] |
												 	mamTax$clade==australidelphiaNames[5]),"tiplabel"])
			} else if(j==3){
				cladeTipNames[[j]]<-as.character(mamTax[which( 	mamTax$clade==atlantogenataNames[1] |
													mamTax$clade==atlantogenataNames[2]),"tiplabel"])
			} else {
				cladeTipNames[[j]]<-as.character(mamTax[which(mamTax$clade==cladeNames[j]),"tiplabel"])				
			}
		}

		# prune tree to each clade
		cladeTrees<-list()
		for(j in 1:length(cladeTipNames)){
			tips<-cladeTipNames[[j]]
			toDrop<-setdiff(mamPhy$tip.label, tips)
			cladeTrees[[j]]<-drop.tip(mamPhy, toDrop)
		}
		#	btimes<-lapply(cladeTrees,branching.times)
		#	bMaxs<-lapply(btimes,max)
		#	min(unlist(bMaxs))
		#	[1] 17.15449 << so that is the youngest root of the clades (Platyrrhini)

	# RUN analyses
	###

	# set resDir
	resDir<-"pulledSpeciation_perClade_100trees/"
	setwd(resDir)

	# set max age
	maxToEst<-15 # use 5 Mya to divide the trees at 1 Mya intervals
	
	res_PSR<-list()
	for(j in 1:length(cladeTrees)){
	# set tree
	cladePhy<-cladeTrees[[j]]
		#root<-max(branching.times(cladePhy))

	res_PSR[[j]]<-fit_hbd_psr_on_grid(  tree = cladePhy, 
                           oldest_age            = maxToEst,
                           age0                  = 0,
                           age_grid              = seq(0, maxToEst, by=1),
                           min_PSR               = 0.01,
                           max_PSR               = 100,
                           guess_PSR             = NULL,
                           fixed_PSR             = NULL,
                           splines_degree        = 1,
                           condition             = "auto",
                           relative_dt           = 1e-3,
                           Ntrials               = 10,
                           Nbootstraps           = 10,
                           Ntrials_per_bootstrap = 10,
                           Nthreads              = 1,
                           max_model_runtime     = NULL,
                           fit_control           = list(),
                           verbose               = FALSE,
                           verbose_prefix        = "")
	}
	save(res_PSR, file=paste0("PSR_5911species_NDexp_tree",i,"_20clades.Rda"))

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
