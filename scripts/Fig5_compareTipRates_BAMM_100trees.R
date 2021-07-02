# Code to compare the tipRates to each other and BAMM shifts
#########
# For Fig 4 of new conception RSPB
#######
# NS Upham; last modification: 17 Nov 2020
# ===================================================


# FIG 4
###########
	dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses"
	setwd(dirname)

	library(viridis); library(plotrix); library(phyloch); library(phytools); library(moments) #library(vioplot)
	# net diversification
	# WITH the BAMM shifts
	####
	intCol<-"black"#"darkorchid4"#grey(0.5)
	medCol<-"white"#"mediumorchid1"#"plum1"#grey(0.1)
	intervalDims<-seq(from=0.025,0.975,by=0.01)
	bbone<-"NDexp"

	# BAMM lineage-level shifts 
		# SUMMARY get data
		agesAllClades_origALL<-read.table(file=paste("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/divTime_24BAMMnodes_",bbone,"_MCC_target_DETAILS_marsupOK.txt", sep=""), header=TRUE)
		firstShift<-c(1,3,5,7,8,10,11,13,15,17:29,31,32) #5 == keeping Placentalia
		shiftDivDat1<-agesAllClades_origALL[firstShift,c("CLADE_label","ID_label","avgFactor","avgIncDec","mean","lower","upper","avgCladeDiv","cladeDiv_low95","cladeDiv_high95","nodeMCC","numTaxa", "from","to", "Indep")]
		shiftDivDat<-shiftDivDat1[order(shiftDivDat1$avgFactor),]

		# SUMMARY prepare plotting symbols
	    pchs<-rep(21,length(shiftDivDat[,1]))
	    alphaCol<-0.8
	    ptCols<-rep(rgb(1,0,0,alpha=alphaCol),length(shiftDivDat[,1]))
	    ptCols[which(shiftDivDat$avgIncDec=="downShift")]<-rgb(0,0,1,alpha=alphaCol)
	    ptCols[which(shiftDivDat$avgIncDec=="up-or-down")]<-grey(0.5,alpha=alphaCol)
	
	    yVar<-shiftDivDat$avgFactor
	    xVar<-(-shiftDivDat$mean)
	    xLow<-(-shiftDivDat$lower)
	    xHigh<-(-shiftDivDat$upper)

		#    pdf(file=paste("plotCI_cladeBAMMfactor-vs-cladeAge_timeBars_agesOK.pdf",sep=""),onefile=TRUE, width=6,height=4)
		#
		#    plotCI(x=xVar, y=yVar,ui=xLow,li=xHigh, xlim=c(min(na.omit(xHigh)),max(na.omit(xLow))), ylim=c(min(yVar),max(yVar)), 
		#    	xlab="", ylab="", sfrac=0,err="x", lwd=2,col="black",pt.bg=ptCols, scol="black",pch=pchs,font.lab=2,cex.axis=0.95,cex.lab=1, cex=1.5,xaxt="n", yaxt="n", bty="n")
		#    axis(side=1,line=2)
		#	axis(side=4, at=c(1:4))
		#
		#    #rug(xVar[which(shiftDivDat$avgIncDec=="upShift")], col=rgb(1,0,0,alpha=alphaCol), side=1, lwd=2, line=1)
		#    #rug(xVar[which(shiftDivDat$avgIncDec=="downShift")], col=rgb(0,0,1,alpha=alphaCol), side=1, lwd=2, line=2)
		#    #rug(xVar[which(shiftDivDat$avgIncDec=="up-or-down")], col=grey(0.5,alpha=alphaCol), side=1, lwd=2, line=2)
		#
		#    dev.off()

	# Compare BAMM Shifts to tip DR distributuon for same clades.
	####
		# get the TIP DR VALUES per shift
		####
			# get the species per shift
			###
			# load MCC tree
			mamMCC<-drop.tip(read.nexus(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"), "_Anolis_carolinensis")

			tipsPerShift<-vector("list",length(shiftDivDat[,1]))
			for(j in 1:length(shiftDivDat[,1])){
				tipFrom <- as.character(shiftDivDat[j,"from"])
				tipTo <- as.character(shiftDivDat[j,"to"])
				nodeToGet <- getMRCA(mamMCC, tip=c(tipFrom,tipTo))#, type="node")
				allNums <- getDescendants(mamMCC, node=nodeToGet)#-length(mamMCC$tip.label))
				#allNums <- getDescendants(mamMCC, node=shiftDivDat[j,"nodeMCC"])#(nodeToGet-length(mamMCC$tip.label)))
				#tipNums <- allNums[allNums-length(mamMCC$tip.label) > 0]

				tipsPerShift[[j]] <- as.character(na.omit(mamMCC$tip.label[allNums]))
			}

			# load the tip DR values, 10k tree summary
			tipDR10k <- read.table(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL/MammaliaTrees_SepOct2019/node-dated_backbone_NDexp_10ktrees_tipDR/DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_expanded.txt", header=TRUE)
			tipDR_perShift_95CI<-vector("list",length(shiftDivDat[,1]))
			tipDR_perShift_skew<-c()
			for(j in 1:length(shiftDivDat[,1])){
				tips_j <- tipsPerShift[[j]]
				tipDR_j <- tipDR10k[ rownames(tipDR10k) %in% tips_j, "harmMeans"]
				tipDR_perShift_95CI[[j]] <- quantile(tipDR_j, c(0.5, 0.025, 0.975))
				tipDR_perShift_skew[j] <- skewness(tipDR_j)
			}
			tipDR_perShift_95CI_ALL<-cbind.data.frame( do.call(rbind,tipDR_perShift_95CI), tipDR_perShift_skew)
			rownames(tipDR_perShift_95CI_ALL)<-shiftDivDat$CLADE_label
			colnames(tipDR_perShift_95CI_ALL)<-c("tipDR_median","tipDR_low95","tipDR_up95", "tipDR_skew")

			allDat<-cbind.data.frame(shiftDivDat, tipDR_perShift_95CI_ALL)
			write.csv(allDat, file=paste0("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/divTime_24BAMMnodes_",bbone,"_MCC_target_DETAILS_marsupOK_wTipDR.txt"))
    
	# read back in
	allDat_all24<-read.csv(file=paste0("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/divTime_24BAMMnodes_",bbone,"_MCC_target_DETAILS_marsupOK_wTipDR.txt"), header=TRUE)
		# subset to the 18 independent shifts (excludes A, C, I, K, S, U)
		allDat <- allDat_all24[which(allDat_all24$Indep==1),]
			cladeNum<-c(13, 17, 20, NA, 10, 13, 20, 20, 12, NA, 14, 20, 11, 6, 2, 20, 12, 18)
		allDat_BAMM<-cbind.data.frame(allDat, cladeNum)

    alphaLevel=0.7
	rodCol<-viridis(10,alpha=alphaLevel)[2]
	batCol<-viridis(10,alpha=alphaLevel)[4]
	shrewCol<-viridis(10,alpha=alphaLevel)[6]
	priCol<-viridis(10,alpha=alphaLevel)[8] #"chartreuse3"
	artCol<-viridis(10,alpha=alphaLevel)[10] #"gold1"
	oCol<-grey(0.5,alpha=alphaLevel)

	shiftCols<-c(batCol,oCol,rodCol,oCol,artCol,batCol,
				rodCol,rodCol,batCol,
				priCol,priCol,rodCol,artCol,shrewCol,
				oCol,rodCol,batCol,rodCol)

    alphaLevel=1.0
	rodCol<-viridis(10,alpha=alphaLevel)[2]
	batCol<-viridis(10,alpha=alphaLevel)[4]
	shrewCol<-viridis(10,alpha=alphaLevel)[6]
	priCol<-viridis(10,alpha=alphaLevel)[8] #"chartreuse3"
	artCol<-viridis(10,alpha=alphaLevel)[10] #"gold1"
	oCol<-grey(0.5,alpha=alphaLevel)

	borderCols<-c(batCol,oCol,rodCol,oCol,artCol,batCol,
				rodCol,rodCol,batCol,
				priCol,priCol,rodCol,artCol,shrewCol,
				oCol,rodCol,batCol,rodCol)

	 #pdf(file=paste0("plotCI_cladeTipDRSkew-vs-cladeTipDR.pdf"),onefile=TRUE, width=6,height=6)
	 #pdf(file=paste0("plotCI_cladeBAMMdiv-vs-cladeTipDR.pdf"),onefile=TRUE, width=6,height=6)
	 pdf(file=paste0("plotCI_cladeBAMMdiv-vs-cladeTipDR-mean_labBigger_18Shifts.pdf"),onefile=TRUE, width=6,height=6)
	
	 plotCI(x=allDat$tipDR_median, y=allDat$avgCladeDiv,ui=allDat$tipDR_up95,li=allDat$tipDR_low95, 
	 	#xlim=c(min(allDat$tipDR_low95),max(allDat$tipDR_up95)), ylim=c(min(allDat$tipDR_low95),max(allDat$tipDR_up95)), #c(min(allDat$cladeDiv_low95),max(allDat$cladeDiv_high95)), 
	 	xlim=c(0.1,0.6), ylim=c(0.1,0.6), #c(min(allDat$cladeDiv_low95),max(allDat$cladeDiv_high95)), 
	 	xlab="", ylab="", sfrac=0,err="x", xpd=NA, lwd=2, cex=1, col="black",pt.bg=shiftCols, scol=grey(0.65),pch=NA,font.lab=2,cex.axis=0.95,cex.lab=1,xaxt="n", yaxt="n", bty="n")

	 plotCI(x=allDat$tipDR_median, y=allDat$avgCladeDiv,ui=allDat$cladeDiv_high95,li=allDat$cladeDiv_low95, 
	 	add=TRUE, sfrac=0,err="y", xpd=NA, lwd=2,col=borderCols,cex=allDat$avgFactor*1.2, pt.bg=shiftCols, scol=grey(0.65),pch=21,font.lab=2,cex.axis=0.95,cex.lab=1)

	 axis(side=1,line=2, cex.axis=1.3)
	 axis(side=2, line=2, cex.axis=1.3)

	 datToModel<-allDat[which(allDat$Indep==1),]
	 modSum<-summary(lm(datToModel$avgCladeDiv ~ datToModel$tipDR_median))
	abline(a= modSum$coef[1], b=modSum$coef[2], lty=2, lwd=2, col="black")
	#abline(a= 0, b=1, lty=1, lwd=2, col="black")

	mtext(side=3, text=paste0(round(modSum$coef[1],2), " + ", round(modSum$coef[2],2), "x ; R2 = ",round(modSum$r.squared,2)))
	text(x= allDat$tipDR_median, y= allDat$avgCladeDiv, labels=allDat$ID_label, cex=1.5, font=2, adj=c(-0.8, -0.8))
	 #rug(xVar[which(shiftDivDat$avgIncDec=="upShift")], col=rgb(1,0,0,alpha=alphaCol), side=1, lwd=2, line=1)
	 #rug(xVar[which(shiftDivDat$avgIncDec=="downShift")], col=rgb(0,0,1,alpha=alphaCol), side=1, lwd=2, line=2)
	 #rug(xVar[which(shiftDivDat$avgIncDec=="up-or-down")], col=grey(0.5,alpha=alphaCol), side=1, lwd=2, line=2)
	
	 dev.off()


	# Pulled speciation rates ~ tip DR mean
	#############

	   	# EMPIRICAL -- get the TIP DR VALUES per SUBCLADE
	   	####
	   		# get the species per SUBCLADE
			###
	 		mamTax<-read.csv(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL/MammaliaTrees_SepOct2019/node-dated_backbone_NDexp_10ktrees_tipDR/taxonomy_mamPhy_5911species.csv", header=TRUE)

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

			library(data.table); library(psych)
			# load the tip DR values, 10k tree summary
			tipDR10k <- read.table(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL/MammaliaTrees_SepOct2019/node-dated_backbone_NDexp_10ktrees_tipDR/DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_expanded.txt", header=TRUE)
			tipDR_perClade_95CI<-vector("list",length(cladeTipNames))
			tipDR_perClade_skew<-c()
			for(j in 1:length(cladeTipNames)){
				tips_j <- cladeTipNames[[j]]
				tipDR_j <- tipDR10k[ rownames(tipDR10k) %in% tips_j, "harmMeans"]
				tipDR_perClade_95CI[[j]] <- quantile(tipDR_j, c(0.5, 0.025, 0.975))
				tipDR_perClade_skew[j] <- skewness(tipDR_j)
			}
			tipDR_perClade_95CI_ALL<-cbind.data.frame( do.call(rbind,tipDR_perClade_95CI), tipDR_perClade_skew)
			rownames(tipDR_perClade_95CI_ALL)<-cladeNames
				colnames(tipDR_perClade_95CI_ALL)<-c("tipDR_median","tipDR_low95","tipDR_up95", "tipDR_skew")

			# load the tipDR values ALL 10k trees
			tipDR10k_allTrees <- fread(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL/MammaliaTrees_SepOct2019/node-dated_backbone_NDexp_10ktrees_tipDR/DR-matrix_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2.txt")
				# now for each tree, get those data
				tipDR_i_j_95CI_ALL<-list()
				for(i in 2:10000){
		
					tipDR_i_j_harm<-c(); tipDR_i_j_skew<-c()
					for(j in 1:length(cladeTipNames)){
						tips_j <- cladeTipNames[[j]]
						tipNums_j <- as.numeric(c(1:5912)[as.character(unlist(tipDR10k_allTrees[,1])) %in% tips_j])

						select_col<-paste0("V",i)
						tipDR_i_j <- unlist(tipDR10k_allTrees[ tipNums_j, ..select_col])

						tipDR_i_j_harm[j] <- harmonic.mean(tipDR_i_j)
						tipDR_i_j_skew[j] <- skewness(tipDR_i_j)
					}
					tipDR_i_j_95CI_ALL[[i]]<-cbind.data.frame(clade=cladeNames, tipDR_i_j_harm, tipDR_i_j_skew, tree=i-1)
				}
				tipDR_i_j_95CI_ALL_10k<-do.call(rbind,tipDR_i_j_95CI_ALL)
				write.csv(tipDR_i_j_95CI_ALL_10k, file="tipDR_20clades_harmMean_skew_all10ktrees_all.csv")

				# get medians and 95% CI of those tipDR distribution MOMENTS
				cladeDR<-list(); cladeSkew<-list()
				for(j in 1:length(cladeNames)){
					cladeDat<-tipDR_i_j_95CI_ALL_10k[which(tipDR_i_j_95CI_ALL_10k$clade==cladeNames[j]),]
					cladeDR[[j]]<-quantile(cladeDat$tipDR_i_j_harm, c(0.5, 0.025, 0.975))
					cladeSkew[[j]]<-quantile(cladeDat$tipDR_i_j_skew, c(0.5, 0.025, 0.975))
				}
				cladeDR_ALL<-do.call(rbind, cladeDR)
				cladeSkew_ALL<-do.call(rbind, cladeSkew)
				clade_tipDR_sum<-cbind.data.frame(clade=cladeNames, 
									tipDR_harm_med=cladeDR_ALL[,1], tipDR_harm_low95=cladeDR_ALL[,2], tipDR_harm_high95=cladeDR_ALL[,3],
									tipDR_skew_med=cladeSkew_ALL[,1], tipDR_skew_low95=cladeSkew_ALL[,2], tipDR_skew_high95=cladeSkew_ALL[,3])

				write.csv(clade_tipDR_sum, file="tipDR_20clades_harmMean_skew_all10ktrees_SUMMARY.csv")


		   	# SIMULATED -- get the TIP DR VALUES per SUBCLADE
		   	####
				setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/tipDR_clades_compareVersusSIMS")
				library(psych); library(moments)

		   		# get the species per SUBCLADE
				###
				cladeNames<-c("DIDELPHIMORPHIA", #"Australidelphia",
							#"Atlantogenata",
							"Talpidae",# "Erinaceidae", #<< exclude 24 species
							"Soricidae",
							"Feliformes","Caniformes",
							#"Suina", #<< exclude 21 species
							"Whippomorpha","Ruminantia",
							"Yinpterochiroptera","Yangochiroptera",
							"Strepsirrhini", "Catarrhini", "Platyrrhini", 
							"LAGOMORPHA", "Guinea_pig-related", "Squirrel-related", "Mouse-related")
				ntrees<-100

				# load the TREES
				load("MamPhy_SIMS_mamPhyE_NDexp_CLADE_27-all-loaded.Rda")
				load("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_CLADES_27-all-loaded.Rda")

				# load the CALCULATED DR values -- PER CLADE
				simDR<-list(); obsDR<-list()
				tipDR_i_j_skew_SIM<-list(); tipDR_i_j_skew_OBS<-list()
				mannWhitneyResult<-c(); pVal<-c(); fullTest<-list()

				for(j in 1:length(cladeNames)){
					# per clade, 100 trees tip DR values
					simDR[[j]]<-read.table(paste0("DRtips_MamPhy_SIMS_mamPhyE_NDexp_",cladeNames[j],"_all100trees.txt"), header=TRUE)
					obsDR[[j]]<-read.table(paste0("DRtips_MamPhy_OBS_NDexp_",cladeNames[j],"_all100trees.txt"), header=TRUE)

					# per tree, get the harmonic mean & skewness of tip DR
					tipDR_skew_SIM<-c(); tipDR_skew_OBS<-c()
					for(i in 1:ntrees){
						perTree_SIM<-simDR[[j]][,i]
						perTree_OBS<-obsDR[[j]][,i]

						tipDR_skew_SIM[i] <- skewness(perTree_SIM)
						tipDR_skew_OBS[i] <- skewness(perTree_OBS)

					}
						fullTest[[j]]<- wilcox.test(x=tipDR_skew_OBS, y=tipDR_skew_SIM, alternative="two.sided")
							mannWhitneyResult[j]<-if(fullTest[[j]]$p.value < 0.05){ "different" } else { "not different" }
							pVal[j]<-fullTest[[j]]$p.value
					
					tipDR_i_j_skew_SIM[[j]]<-quantile(tipDR_skew_SIM, c(0.5, 0.025, 0.975))
					tipDR_i_j_skew_OBS[[j]]<-quantile(tipDR_skew_OBS, c(0.5, 0.025, 0.975))
				}

				tipDR_skew_ALL<-cbind.data.frame(clade=cladeNames, do.call(rbind,tipDR_i_j_skew_OBS), do.call(rbind,tipDR_i_j_skew_SIM), mannWhitneyResult, pVal)
				colnames(tipDR_skew_ALL)<-c("clade",paste0("skew_OBS_",c("med","low95","high95")), paste0("skew_SIM_",c("med","low95","high95")), "mannWhitneyResult", "pVal")
				
				write.csv(tipDR_skew_ALL, file="tipDR_16clades_SIMS_summary.csv")




		# load the PULLED SPECIATION RATES per SUBCLADE (already calculated on Grace HPC)
		##########
			# get the time=0 (instant present) pulled rate per clade, per tree
			ntrees<-100	
			library(miscTools)

			PSR_0_res_allClades<-list()
			for(i in 1:ntrees){
				xx<-load(file=paste0("pulledSpeciation_perClade_100trees/PSR_5911species_NDexp_tree",i,"_20clades.Rda"))
				psrClades_i<-get(xx)

				PSR_0_res<-list()
				for(j in 1:length(cladeNames)){
					psrClade <- psrClades_i[[j]]
						PSR_0_est <- psrClade$fitted_PSR[1]
						PSR_0_low95 <- psrClade$CI95lower[1]
						PSR_0_high95 <- psrClade$CI95upper[1]
						clade <- cladeNames[j]
						tree <- i
					PSR_0_res[[j]] <- cbind(clade, PSR_0_est, PSR_0_low95, PSR_0_high95, tree)
				}
				PSR_0_res_allClades[[i]]<- as.data.frame(do.call(rbind,PSR_0_res))
			}
			PSR_0_res_allClades_ALL<-do.call(rbind,PSR_0_res_allClades)

			write.csv(PSR_0_res_allClades_ALL, file="PSR0_20clades_100trees_all.csv")

			# get the medians of those values across 100 trees
			psrDat_sum<-list()
			for(j in 1:length(cladeNames)){
				psrDat<-PSR_0_res_allClades_ALL[which(PSR_0_res_allClades_ALL$clade==cladeNames[j]),]
				psrDat_sum[[j]]<-cbind.data.frame(clade=cladeNames[j], PSR_0_est=median( as.numeric(as.character(psrDat[,2])) ), PSR_0_low95=median( as.numeric(as.character(psrDat[,3])) ), PSR_0_high95=median( as.numeric(as.character(psrDat[,4])) ) )
			}
			psrDat_sum_ALL<-do.call(rbind,psrDat_sum)

			write.csv(psrDat_sum_ALL, file="PSR0_20clades_100trees_SUMMARY.csv")

		# COMBINE
			allDat_psrDR<-cbind.data.frame(cladeNum=c(1:20), clade_tipDR_sum, psrDat_sum_ALL[,-1])
			write.csv(allDat_psrDR, file="tipDR_20clades_harmMean_skew_all10ktrees_SUMMARY_withPSR0.csv")

	# Read back in:
	#####
	allDat_psrDR<-read.csv(file="tipDR_20clades_harmMean_skew_all10ktrees_SUMMARY_withPSR0.csv")

		# now PLOT
		#######
	    alphaLevel=0.7
		rodCol<-viridis(10,alpha=alphaLevel)[2]
		batCol<-viridis(10,alpha=alphaLevel)[4]
		shrewCol<-viridis(10,alpha=alphaLevel)[6]
		priCol<-viridis(10,alpha=alphaLevel)[8] #"chartreuse3"
		artCol<-viridis(10,alpha=alphaLevel)[10] #"gold1"
		oCol<-grey(0.5,alpha=alphaLevel)

		cladeCols<-c(oCol,oCol,oCol,shrewCol,shrewCol,shrewCol,
					oCol,oCol,artCol,artCol,artCol,
					batCol,batCol,priCol,priCol,priCol,oCol,
					rodCol,rodCol,rodCol)

	    alphaLevel=1.0
		rodCol<-viridis(10,alpha=alphaLevel)[2]
		batCol<-viridis(10,alpha=alphaLevel)[4]
		shrewCol<-viridis(10,alpha=alphaLevel)[6]
		priCol<-viridis(10,alpha=alphaLevel)[8] #"chartreuse3"
		artCol<-"gold2" #viridis(10,alpha=alphaLevel)[10] #
		oCol<-grey(0.5,alpha=alphaLevel)

		borderCols<-c(oCol,oCol,oCol,shrewCol,shrewCol,shrewCol,
					oCol,oCol,artCol,artCol,artCol,
					batCol,batCol,priCol,priCol,priCol,oCol,
					rodCol,rodCol,rodCol)

		pdf(file=paste0("plotCI_PSR0-vs-cladeTipDR-mean_20Clades.pdf"), height=6, width=6, onefile=TRUE)

			dat<-allDat_psrDR
			plotCI(x=dat$tipDR_harm_med, y=dat$PSR_0_est,ui=dat$tipDR_harm_high95,li=dat$tipDR_harm_low95, 
				#xlim=c(min(allDat$tipDR_low95),max(allDat$tipDR_up95)), ylim=c(min(allDat$tipDR_low95),max(allDat$tipDR_up95)), #c(min(allDat$cladeDiv_low95),max(allDat$cladeDiv_high95)), 
				xlim=c(0.05,0.55), ylim=c(0.05,0.55), #c(min(allDat$cladeDiv_low95),max(allDat$cladeDiv_high95)), 
				xlab="", ylab="", sfrac=0,err="x", xpd=NA, lwd=2, cex=1, col="black", scol=grey(0.65),pch=NA,font.lab=2,cex.axis=0.95,cex.lab=1,xaxt="n", yaxt="n", bty="n")

			plotCI(x=dat$tipDR_harm_med, y=dat$PSR_0_est,ui=dat$PSR_0_high95,li=dat$PSR_0_low95, 
				add=TRUE, sfrac=0,err="y", xpd=NA, lwd=2,col=borderCols,cex=2, pt.bg=cladeCols, scol=grey(0.65),pch=21,font.lab=2,cex.axis=0.95,cex.lab=1)
		
			axis(side=1,line=2, cex.axis=1.3)
			axis(side=2, line=2, cex.axis=1.3)

			datToModel<-dat
			modSum<-summary(lm(datToModel$PSR_0_est ~ datToModel$tipDR_harm_med))
			abline(a= modSum$coef[1], b=modSum$coef[2], lty=2, lwd=2, col="black")
			#abline(a= 0, b=1, lty=1, lwd=2, col="black")

			mtext(side=3, text=paste0(round(modSum$coef[1],2), " + ", round(modSum$coef[2],2), "x ; R2 = ",round(modSum$r.squared,2)))
			text(x= dat$tipDR_harm_med, y= dat$PSR_0_est, labels=c(1:20), cex=1.5, col=borderCols, font=2, adj=c(-0.8, -0.8))

		dev.off()

			# Coefficients:
			#                           Estimate Std. Error t value Pr(>|t|)   
			# (Intercept)                0.02323    0.04330   0.536  0.59820   
			# datToModel$tipDR_harm_med  0.78295    0.23025   3.400  0.00319 **
			# ---
			# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
			# 
			# Residual standard error: 0.06883 on 18 degrees of freedom
			# Multiple R-squared:  0.3911,	Adjusted R-squared:  0.3573 
			# F-statistic: 11.56 on 1 and 18 DF,  p-value: 0.003188



	# sumBAMMShiftFactors ~ tip DR skew PER CLADE
	########
		# sum up the BAMM shift factors
		cladesWithShifts<-names(table(allDat_BAMM$cladeNum))
		
		cladeTotalShiftFactors<-c()
		for(q in 1:length(cladesWithShifts)){
			cladeTotalShiftFactors[q] <- sum(allDat_BAMM[which(allDat_BAMM$cladeNum==cladesWithShifts[q]),"avgFactor"])
		}
		names(cladeTotalShiftFactors)<-cladesWithShifts

		# add to other matrix
		allDat_psrDR_wShifts<-cbind.data.frame(allDat_psrDR, totalShiftFactors=rep(0,20))
		allDat_psrDR_wShifts[ (allDat_psrDR_wShifts$cladeNum %in% names(cladeTotalShiftFactors)), "totalShiftFactors" ] <- as.numeric(cladeTotalShiftFactors)
			# make the Strepsirrhini shift negative
			allDat_psrDR_wShifts[14,"totalShiftFactors"]<- -1.7

			write.csv(allDat_psrDR_wShifts, file="tipDR_20clades_harmMean_skew_all10ktrees_SUMMARY_withPSR0_withBAMMshifts.csv")


		# now PLOT
		####
		cladeNumLabs<-c(1:20)
		toUnlabel<-c(1:3,7:8,17)
		cladeNumLabs[toUnlabel]<-""

		pdf(file=paste0("plotCI_cladeBAMMshiftFactor-vs-cladeTipDRSkew_20Clades_negStrep.pdf"), height=6, width=6, onefile=TRUE)
			
			dat<-allDat_psrDR_wShifts
			plotCI(x=dat$tipDR_skew_med, y=dat$totalShiftFactors,ui=dat$tipDR_skew_high95,li=dat$tipDR_skew_low95, 
				#xlim=c(min(allDat$tipDR_low95),max(allDat$tipDR_up95)), ylim=c(min(allDat$tipDR_low95),max(allDat$tipDR_up95)), #c(min(allDat$cladeDiv_low95),max(allDat$cladeDiv_high95)), 
				xlim=c(-0.8,1.5), ylim=c(-2,8.5), #c(min(allDat$cladeDiv_low95),max(allDat$cladeDiv_high95)), 
				xlab="", ylab="", sfrac=0,err="x", xpd=NA, lwd=2, col=borderCols,cex=2, pt.bg=cladeCols, scol=grey(0.65),pch=21,font.lab=2,cex.axis=0.95,cex.lab=1,xaxt="n", yaxt="n", bty="n")

			#plotCI(x=dat$tipDR_harm_med, y=dat$PSR_0_est,ui=dat$PSR_0_high95,li=dat$PSR_0_low95, 
			#	add=TRUE, sfrac=0,err="y", xpd=NA, lwd=2,col=borderCols,cex=2, pt.bg=cladeCols, scol=grey(0.65),pch=21,font.lab=2,cex.axis=0.95,cex.lab=1)
		
			axis(side=1,line=2, cex.axis=1.3)
			axis(side=2, line=2, cex.axis=1.3)

			datToModel<-dat
			modSum<-summary(lm(datToModel$totalShiftFactors ~ datToModel$tipDR_skew_med))
			abline(a= modSum$coef[1], b=modSum$coef[2], lty=2, lwd=2, col="black")
			#abline(a= 0, b=1, lty=1, lwd=2, col="black")

			mtext(side=3, text=paste0(round(modSum$coef[1],2), " + ", round(modSum$coef[2],2), "x ; R2 = ",round(modSum$r.squared,2)))
			text(x= dat$tipDR_skew_med, y= dat$totalShiftFactors, labels=cladeNumLabs, cex=1.5, col=borderCols, font=2, adj=c(-0.4, -0.4))

		dev.off()

			# Coefficients:
			#                           Estimate Std. Error t value Pr(>|t|)    
			# (Intercept)                -0.1541     0.5194  -0.297 0.770056    
			# datToModel$tipDR_skew_med   3.2776     0.8105   4.044 0.000762 ***
			# ---
			# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
			# 
			# Residual standard error: 1.671 on 18 degrees of freedom
			# Multiple R-squared:  0.476,	Adjusted R-squared:  0.4469 
			# F-statistic: 16.35 on 1 and 18 DF,  p-value: 0.0007619



# SIMILATE THOSE PHYLOGENIES::

# Plot -- Part E (PART 1) -- Likelihood model fitting RC vs RV models in RPANDA 
########
# packages
library(ape); library(phytools); library(moments)

# directory
setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

# specify the backbone to use - FBD or NDexp
bbone <- "NDexp" # "FBD"


# ML models of lineage diversification
#######

# Load in data
#mamFBD_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")
mamNDexp_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees")

cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)
sorted<-cladesDR[order(cladesDR$tiplabel),]

taxonomyVars<-c("tiplabel","gen", "fam", "ord", "higher", "clade", "genes","extinct.")
taxFinal<-sorted[,taxonomyVars]

cladesDR<-taxFinal

ordNames<-names(which(table(cladesDR$ord) >= 25))
cladeNames_1<-names( which(table(cladesDR$clade) >= 25 ))
cladeNames<-cladeNames_1[ which(cladeNames_1!="Dasyuromorphia" & cladeNames_1!="Didelphimorphia" & cladeNames_1 !="Diprodontia" & cladeNames_1!="Lagomorpha") ]

allOrdsClades_names<-c(ordNames,cladeNames)


# ORDERS breakout
ordTipNames<-vector("list",length(ordNames))
for (i in 1:length(ordNames)){
    ordTipNames[[i]]<-cladesDR[which(cladesDR$ord==ordNames[i]),"tiplabel"]
}

mamNDexp_100_ord<-vector("list",length(mamNDexp_100))
for (k in 1:length(ordTipNames)){
    toDrop<-setdiff(mamNDexp_100[[1]]$tip.label,ordTipNames[[k]])
    
    for (i in 1:length(mamNDexp_100)){
        mamNDexp_100_ord[[i]]<-drop.tip(mamNDexp_100[[i]],toDrop)
    }
    write.nexus(mamNDexp_100_ord, file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames[k],".trees",sep=""))
}

# CLADES breakout
cladeTipNames<-vector("list",length(cladeNames))
for (i in 1:length(cladeNames)){
    cladeTipNames[[i]]<-cladesDR[which(cladesDR$clade==cladeNames[i]),"tiplabel"]
}

mamNDexp_100_clade<-vector("list",length(mamNDexp_100))
for (k in 1:length(cladeTipNames)){
    toDrop<-setdiff(mamNDexp_100[[1]]$tip.label,cladeTipNames[[k]])
    
    for (i in 1:length(mamNDexp_100)){
        mamNDexp_100_clade[[i]]<-drop.tip(mamNDexp_100[[i]],toDrop)
    }
    write.nexus(mamNDexp_100_clade, file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_CLADES_",cladeNames[k],".trees",sep=""))
}

setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/__mamPhy_31namedClades_results")
# re-load
for (i in 1:length(ordNames)){
    assign(paste(ordNames[i],".trees",sep=""), read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_ORDS_",ordNames[i],".trees",sep=""))) 
    }
for (i in 1:length(cladeNames)){
    assign(paste(cladeNames[i],".trees",sep=""), read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_CLADES_",cladeNames[i],".trees",sep=""))) 
    }

# put all together = 11 + 16 = 27
allOrdsClades_100trees_first<-list( lapply(paste(ordNames,".trees",sep=""),get), lapply(paste(cladeNames,".trees",sep=""), get) )

allOrdsClades_100trees<-unlist(allOrdsClades_100trees_first, recursive=FALSE)

allOrdsClades_names<-c(ordNames,cladeNames)
save(allOrdsClades_100trees, file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_CLADES_31-all-loaded.Rda")


####
# SIMULATIONS -- need simulated rate-constant trees of the same size and age as each of the 31 clades.
# Make SIMULATIONS of the CLADES for testing these null expectations of RC vs RV models
#===============
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)

ntrees<-100

setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/__mamPhy_31namedClades_results/_sims")

for(j in 1:length(allOrdsClades_100trees)){
    treeSet<-allOrdsClades_100trees[[j]]

    simMam_100<-vector("list",ntrees)
    for(i in 1:ntrees){

        tree_i<-treeSet[[i]]
        root<-max(node.age(tree_i)$ages)

        # First I want to get the ML birth and death rates for each clade
        bd<-birthdeath(tree_i)

        extFrac<-bd$para[[1]] # d/b
        netDiv<-bd$para[[2]] # b-d
        lam<-netDiv/(1-extFrac)
        mu<-lam-netDiv
        # Then sim for mamPhy rates, diversity, and re-scale to age
        for (k in 1:100){
        simMam<-pbtree(b=lam,d=mu,n=length(tree_i$tip.label),t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
        if (class(simMam) == "NULL"){
            cat("another sim...")
            simMam<-pbtree(b=lam,d=mu,n=length(tree_i$tip.label),t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
            } else break 
            cat("got one, next simulation...") 
        }
        simMam_100[[i]]<-simMam  
        class(simMam_100)<-"multiPhylo"
    }
    check<-lapply(simMam_100, is.null)
    if( length(check[which(check==FALSE)]) == 100 ){
    write.nexus(simMam_100, file=paste("MamPhy_SIMS_mamPhyE_",bbone,"_CLADE_",allOrdsClades_names[j],"_100trees.trees",sep=""))
    } else break

}


## CALCULATE THE TIP DR RATES::

	# code.... to find



