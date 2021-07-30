# Code to plot the RTT output of 10 trees run in BAMM
#########
# For Fig 3 and 4a of new conception RSPB
#######
# NS Upham; last modification: 11 May 2020
# ===================================================

library(BAMMtools); library(coda); library(phytools); library(ape); library(phangorn)

folders<-c(1:10)
bbone = "NDexp" # "FBD" 

rateMatrices_1Ma<-vector("list",length(folders))
for (i in 1:length(folders)){
	#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",dir,folders[i],sep=""))
	setwd(paste("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_",folders[i],sep=""))

	# load edata
	tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
	edata <- getEventData(tree, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)

	root<-max(branching.times(tree))

	rateMatrices_1Ma[[i]]<-getRateThroughTimeMatrix(edata,nslices=(root/1))

}
#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/",sep=""))
setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_10trees")
# now SAVE that file:
	# save(rateMatrices_1Ma,file="rateThroughTimeMatrix_10trees_1Ma.Rda")
# now LOAD BACK IN that file:
	load(file="rateThroughTimeMatrix_10trees_1Ma.Rda")


# plots as 10 separate RTTs (1 per tree)
# ===========

	# RATES through time...
	pdf(file=paste("rateThroughTime_spec_ext_new_sameAxes_greyScale.1-10.pdf",sep=""), height=5, width=7, onefile=TRUE)
	for(i in 1:length(folders)){
		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="speciation", nBins=100, intervalCol=grey(0.4), avgCol="black", axis.labels=FALSE, ylim=c(0,0.3), xlim=c(150,0)) #start.time=150, 

		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="extinction", nBins=100, add=TRUE, intervalCol=grey(0.4), avgCol="white", axis.labels=FALSE, ylim=c(0,0.3), xlim=c(150,0)) #start.time=150, 
		abline(v= 66, col="grey", lty=2, lwd=2)
		mtext(side=3, text=paste("Tree ",i,sep=""), cex=2, font=2, col="black", adj=0.05, line=-3)

		mtext(side=2, text="Rate", cex=1.2, font=1, col="black", line=3)
		#mtext(side=2, text="Time before present", cex=0.8, font=2, col="black", line=3)

	}
	dev.off()

library(viridis); library(plotrix); library(phyloch); library(phytools); library(moments) #library(vioplot)
source("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/_R-CODE/vioplot2_source.R")

#plot(1:10,cex=3,pch=20,col=viridis(10)[c(2,9)])
# plots as COMBINED RTTs (10-tree variation is represented)
# ===========
greyCol<-grey(0.5)
intervalDims<-seq(from=0.025,0.975,by=0.01)

exCol<-viridis(10)[2]#"royalblue1"##rgb(0,0,1,alpha=0.8)
spCol<-"darkgoldenrod1"#"red2"##rgb(1,0,0,alpha=0.8)

	# Speciation & extinction
	####
	png(file="rateThroughTime_spec_ext_new_sameAxes_yellowPurple.1-10_Combo_110ma_95ci.png", height=5, width=7, units="in", res=300)
#	png(file="rateThroughTime_spec_ext_new_sameAxes_redBlue.1-10_Combo_110ma_95ci.png", height=5, width=7, units="in", res=300)
#	pdf(file=paste("rateThroughTime_spec_ext_new_sameAxes_redBlue.1-10_Combo_110ma_95ci.pdf",sep=""), height=5, width=7, onefile=TRUE)
		#par(bg="white")
#		plotRateThroughTime(rateMatrices_1Ma[[1]], ratetype="speciation", nBins=100, interval=intervalDims, intervalCol=greyCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(110,0)) #start.time=150, 
#		plotRateThroughTime(rateMatrices_1Ma[[1]], ratetype="extinction", nBins=100, add=TRUE, interval=intervalDims, intervalCol=greyCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(110,0)) #start.time=150, 
		plotRateThroughTime(rateMatrices_1Ma[[1]], ratetype="speciation", nBins=100, interval=intervalDims, intervalCol=spCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(110,0), xticks=-1, yticks=-1) #start.time=150, 
		plotRateThroughTime(rateMatrices_1Ma[[1]], ratetype="extinction", nBins=100, add=TRUE, interval=intervalDims, intervalCol=exCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(110,0)) #start.time=150, 
#		abline(v= 66, col="grey", lty=2, lwd=2)
#		mtext(side=2, text="Rate", cex=1.2, font=1, col="black", line=3)
#		mtext(side=1, text="Time before present", cex=1.2, font=1, col="black", line=3)

	# plot CIs
	for(i in 2:length(folders)){
#		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="speciation", nBins=100, add=TRUE, interval=intervalDims, intervalCol=greyCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(150,0)) #start.time=150, 
#		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="extinction", nBins=100, add=TRUE, interval=intervalDims, intervalCol=greyCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(150,0)) #start.time=150, 
		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="speciation", nBins=100, add=TRUE, interval=intervalDims, intervalCol=spCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(110,0)) #start.time=150, 
		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="extinction", nBins=100, add=TRUE, interval=intervalDims, intervalCol=exCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(110,0)) #start.time=150, 
	}
	# medians
	for(i in 1:length(folders)){
#		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="speciation", nBins=100, add=TRUE, intervalCol=NULL, avgCol="black", axis.labels=FALSE, ylim=c(0,0.3), xlim=c(150,0)) #start.time=150, 
#		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="extinction", nBins=100, add=TRUE, intervalCol=NULL, avgCol="white", axis.labels=FALSE, ylim=c(0,0.3), xlim=c(150,0)) #start.time=150, 
		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="speciation", nBins=100, add=TRUE, intervalCol=NULL, avgCol="thistle1", axis.labels=FALSE, ylim=c(0,0.3), xlim=c(110,0)) #start.time=150, 
		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="extinction", nBins=100, add=TRUE, intervalCol=NULL, avgCol="mediumorchid2", axis.labels=FALSE, ylim=c(0,0.3), xlim=c(110,0)) #start.time=150, 
	}
	# remove most recent 2 Ma
	rect(xleft=2, ybottom=0, xright=-2, ytop=0.4, col="white", border=NA)

	dev.off()

	# net diversification
	# ONLY
	####
	intCol<-"black"#"darkorchid4"#grey(0.5)
	medCol<-"white"#"mediumorchid1"#"plum1"#grey(0.1)
	intervalDims<-seq(from=0.025,0.975,by=0.01)

	#pdf(file=paste("rateThroughTime_netdiv_new_sameAxes_greyScale.1-10_Combo_110ma_95ci.pdf",sep=""), height=5, width=7, onefile=TRUE)
	#pdf(file=paste("rateThroughTime_netdiv_new_sameAxes_cyan.1-10_Combo_110ma_95ci.pdf",sep=""), height=5, width=7, onefile=TRUE)
	#pdf(file=paste("rateThroughTime_netdiv_new_sameAxes_white.1-10_Combo_110ma_95ci.pdf",sep=""), height=5, width=7, onefile=TRUE)
	png(file=paste("rateThroughTime_netdiv_new_sameAxes_white.1-10_Combo_110ma_95ci_ylim.png",sep=""), height=5, width=7, units="in", res=300)
	#pdf(file=paste("rateThroughTime_netdiv_new_sameAxes_purple.1-10_Combo_110ma_95ci_ylim.pdf",sep=""), height=5, width=7, onefile=TRUE)
	#pdf(file=paste("rateThroughTime_netdiv_new_sameAxes_gold.1-10_Combo_110ma_95ci.pdf",sep=""), height=5, width=7, onefile=TRUE)
	par(oma=c(1,1,1,3)) #‘c(bottom, left, top, right)’

		plotRateThroughTime(rateMatrices_1Ma[[1]], ratetype="netdiv", nBins=100, interval=intervalDims, intervalCol=greyCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.175), xlim=c(110,0), xtick=-1,ytick=-1) #start.time=150, 
#		abline(v= 66, col="grey", lty=2, lwd=2)
#		mtext(side=2, text="Rate", cex=1.2, font=1, col="black", line=3)
#		mtext(side=1, text="Time before present", cex=1.2, font=1, col="black", line=3)

	# plot CIs
	for(i in 2:length(folders)){
		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="netdiv", nBins=100, add=TRUE, interval=intervalDims, intervalCol=intCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(150,0)) #start.time=150, 
	}
	# medians
	for(i in 1:length(folders)){
		plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="netdiv", nBins=100, add=TRUE, intervalCol=NULL, avgCol=medCol, axis.labels=FALSE, ylim=c(0,0.3), xlim=c(150,0)) #start.time=150, 
	}
	# remove most recent 2 Ma
	rect(xleft=2, ybottom=0, xright=-2, ytop=0.4, col="white", border=NA)

	dev.off()




		# PER 10 TREES get data
			# all 10 trees, 253 shifts
			shift10trees<-read.table(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_10trees_MSC-focus/SUM10_results.MSC.nodeRATES.NDexp.txt",header=TRUE)
			# load MCC tree
			mamMCC<-read.nexus(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")
			allDivTimes<-read.csv(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL/divTime_allTable_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.csv")

			# find node on MCC tree:
			shiftDivs<-vector("list",length(shift10trees[,1]))
			for(j in 1:length(shift10trees[,1])){
				node<-findMRCA(tree=mamMCC, tips=c(as.character(shift10trees[j,"from"]),as.character(shift10trees[j,"to"])), type="node")
				shiftDivs[[j]]<-cbind.data.frame(shift10trees[j,],
								MCC_mean=allDivTimes[which(allDivTimes$node==node),"height"],
								MCC_low95=allDivTimes[which(allDivTimes$node==node),"height_95._HPD_MIN"],
								MCC_up95=allDivTimes[which(allDivTimes$node==node),"height_95._HPD_MAX"] )
			}
			shiftDivs_ALL<-do.call(rbind,shiftDivs)
			write.csv(shiftDivs_ALL,file="SUM10_results.MSC.nodeRATES.NDexp_wMCC-ages.csv")

		allShifts_forRug <- -shiftDivs_ALL[,"MCC_mean"]

	# PLOT TOGETHER
	pdf(file=paste("rateThroughTime_netdiv_new_sameAxes_greyScale.1-10_Combo_110ma_95ci_wBAMMshifts_2part.pdf",sep=""), height=5, width=7, onefile=TRUE)
	par(oma=c(1,1,1,3)) #‘c(bottom, left, top, right)’
	# RTT plot
		plotRateThroughTime(rateMatrices_1Ma[[1]], ratetype="netdiv", nBins=100, interval=intervalDims, intervalCol=greyCol, avgCol=NA, axis.labels=FALSE, ylim=c(0.0,0.175), xlim=c(110,0)) #start.time=150, 
		abline(v= 66, col="grey", lty=2, lwd=2)
		mtext(side=2, text="Rate", cex=1.2, font=1, col="black", line=3)
		#mtext(side=1, text="Time before present", cex=1.2, font=1, col="black", line=3)
		#axis(side=2,at=c(0,0.1,0.2))

		# plot CIs
		for(i in 2:length(folders)){
			plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="netdiv", nBins=100, add=TRUE, interval=intervalDims, intervalCol=intCol, avgCol=NA, axis.labels=FALSE, ylim=c(0,0.3)) #start.time=150, 
		}
		# medians
		for(i in 1:length(folders)){
			plotRateThroughTime(rateMatrices_1Ma[[i]], ratetype="netdiv", nBins=100, add=TRUE, intervalCol=NULL, avgCol=medCol, axis.labels=FALSE, ylim=c(0,0.3)) #start.time=150, 
		}
	# remove most recent 2 Ma
	rect(xleft=2, ybottom=0, xright=-2, ytop=0.4, col="white", border=NA)

	#par(new=TRUE)
	# ADD the BAMM shifts
		plotCI(x=xVar, y=yVar,ui=xLow,li=xHigh, ylim=c(1,5), xlim=c(-110,0),
	    	xlab="", ylab="", sfrac=0,err="x", lwd=2,col="black",pt.bg=ptCols, scol="black",pch=pchs,font.lab=2,cex.axis=0.95,cex.lab=1, cex=1.5,xaxt="n", yaxt="n", bty="n", xpd=NA)
	    axis(side=1,line=3)
		axis(side=4, at=c(1:4))
		text(x=xVar+3,y=yVar-0.13,label=shiftDivDat$ID_label)#, adj=c(0.5,0.5))

	    rug(allShifts_forRug[which(shiftDivs_ALL$incDec=="upShift")], ticksize=0.05, col=rgb(1,0,0,alpha=alphaCol), side=1, lwd=2, line=1.25)
	    rug(allShifts_forRug[which(shiftDivs_ALL$incDec=="downShift")], ticksize=0.05, col=rgb(0,0,1,alpha=alphaCol), side=1, lwd=2, line=2.5)
	    #rug(allShifts_forRug[which(shiftDivs_ALL$incDec=="up-or-down")], col=grey(0.5,alpha=alphaCol), side=1, lwd=2, line=5)
	
	dev.off()

length(allShifts_forRug[which(shiftDivs_ALL$incDec=="downShift")])
# 45
quantile( allShifts_forRug[which(shiftDivs_ALL$incDec=="downShift")], c(0.5, 0.025,0.975, 0.25, 0.75) )
#    50%    2.5%   97.5%     25%     75% 
#-44.820 -78.764  -6.888 -50.920 -27.590 


length(allShifts_forRug[which(shiftDivs_ALL$incDec=="upShift")])
# 208
quantile( allShifts_forRug[which(shiftDivs_ALL$incDec=="upShift")], c(0.5, 0.025,0.975, 0.25, 0.75) )
#     50%     2.5%    97.5%      25%      75% 
#-24.9400 -83.5400  -4.6100 -38.7475 -14.1700 

# PLOT RUGS as VIOLIN PLOTS
########
source("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/_R-CODE/vioplot2_source.R")

vioplot2(x=)

vioplot2()#, 
	at=j, range=1.5, h=hVals[j], col=grey(0.7,alpha=0.5), violBord=grey(0.5), horizontal=FALSE, add=TRUE,lty=1, rectCol=grey(0.7,alpha=0.5), colMed=grey(0.2), pchMed=19, drawRect=TRUE, wex=1, boxwex=0.15, rectBord=grey(0.2))

upDat<-allShifts_forRug[which(shiftDivs_ALL$incDec=="upShift")]
downDat<-allShifts_forRug[which(shiftDivs_ALL$incDec=="downShift")]

library(ggplot2); library(dplyr)
# create a dataset
data2 <- data.frame(
  name=c( rep("up shifts",length(upDat)), rep("down shifts",length(downDat)) ),
  value=c( upDat, downDat )
)

# sample size
sample_size = data2 %>% group_by(name) %>% summarize(num=n())

pdf(file="BAMMshifts_10trees_violinPlots.pdf")
# Plot
data2 %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(name, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=value, fill=name)) +
    geom_violin(width=1.4) +
    geom_boxplot(width=0.1, color="white", alpha=0.2)+#, notch=TRUE, notchwidth = 0.5) +
    scale_fill_manual(values = c("red","blue")) +
    #scale_fill_viridis(discrete = TRUE) +
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    xlab("") + ylab("Time")
dev.off()

    # geom_boxplot
    ###
	#    The lower and upper hinges correspond to the first and third
	#     quartiles (the 25th and 75th percentiles). This differs slightly
	#     from the method used by the ‘boxplot()’ function, and may be
	#     apparent with small samples. See ‘boxplot.stats()’ for for more
	#     information on how hinge positions are calculated for ‘boxplot()’.
	#
	#     The upper whisker extends from the hinge to the largest value no
	#     further than 1.5 * IQR from the hinge (where IQR is the
	#     inter-quartile range, or distance between the first and third
	#     quartiles). The lower whisker extends from the hinge to the
	#     smallest value at most 1.5 * IQR of the hinge. Data beyond the end
	#     of the whiskers are called "outlying" points and are plotted
	#     individually.
	#
	#     In a notched box plot, the notches extend ‘1.58 * IQR / sqrt(n)’.
	#     This gives a roughly 95% confidence interval for comparing
	#     medians. See McGill et al. (1978) for more details.
	#    width width of boxplot
	#
	#    ymin lower whisker = smallest observation greater than or equal to
	#         lower hinge - 1.5 * IQR
	#
	#    lower lower hinge, 25% quantile
	#
	#    notchlower lower edge of notch = median - 1.58 * IQR / sqrt(n)
	#
	#    middle median, 50% quantile
	#
	#    notchupper upper edge of notch = median + 1.58 * IQR / sqrt(n)
	#
	#    upper upper hinge, 75% quantile
	#
	#    ymax upper whisker = largest observation less than or equal to
	#         upper hinge + 1.5 * IQR





# OR, if want to FULL tree data files
#######

edataSUB_perTree<-vector("list",length(folders))
for (i in 1:length(folders)){
	#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",dir,folders[i],sep=""))
	setwd(paste("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_",folders[i],sep=""))

	# load edata
	treeFull <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
	edata <- getEventData(treeFull, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)

	edataSUB_perTree[[i]]<-subtreeBAMM(edata, tips=setdiff(treeFull$tip.label,toDrop))
}
# << OK, this currently has 5675 species, extant and non-marine.

setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_STRAPP_analyses")
# now SAVE that workspace...
	# save(edataSUB_perTree,file="eventData_10trees_BAMMdata_subsetToExtantNonmarine_5675sp.RData")
# now LOAD BACK IN that workspace...
	load(file="eventData_10trees_BAMMdata_subsetToExtantNonmarine_5675sp.RData")


	# MSC shifts on the PHYLORATE plot per tree...
	for(i in 1:length(folders)){
		#i=1
		msc.set <- maximumShiftCredibility(edataSUB_perTree[[i]], maximize='product')
		msc.config <- subsetEventData(edataSUB_perTree[[i]], index = msc.set$sampleindex)

		pdf(file=paste("MSC_ShiftSet_LIN_small_spec_new_legend.",folders[i],".pdf",sep=""),width=5,height=5.5) 

		par(oma=rep(0,4), mar=rep(1,4))
		X<-plot.bammdata(msc.config, legend=FALSE, color.interval=c(0,1), lwd=2,spex="s", par.reset=FALSE)#, direction="upwards")
		addBAMMshifts(msc.config, cex = 2)
		mtext(side=3, text=paste("Tree ",folders[i],sep=""), font=2, cex=0.8, line=0, adj=0.1)
		#mtext("MSC speciation rate")
		axisPhylo()
		addBAMMlegend(X, location=c(0,100,100,50))
		#nodelabels(cex=0.3)
		dev.off()
	}








