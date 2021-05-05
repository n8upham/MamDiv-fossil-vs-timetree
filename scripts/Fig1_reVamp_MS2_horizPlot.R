# new FIG 1 -- MS#2 MamPhy -- For PNAS / EcologyLetters
#######


#==================
# Plotting TRAITS on the MCC tree... 
#==================
# initialize packages
library(ape); library(nlme); library(geiger); library(dplyr); library(phytools)
source("https://raw.githubusercontent.com/liamrevell/phytools/master/R/plotTree.wBars.R")
source("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/_R-CODE/source_functions/circleTree_plotting_functions.R")

# directory
dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/"
setwd(dirname)

# load MCC tree to start with
mamMCC<-drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"),"_Anolis_carolinensis")

# load in clade labels
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_NDexp_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)

# load trait data
datOrig_ALL<-read.table("MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt",header=TRUE)
	# 2737 values
# load RANGE data
latDat<-read.table("latLong_minMax_sortAll_5788species.txt",header=TRUE)

# Vars to plot together
	varNames<-c("sciName","DispDist_kmMAX_final","Activity123_BM_med","Lat_centroid")
	varsToPlot<-left_join(x=datOrig_ALL[,varNames], y=latDat, by="sciName")

	# make Diurnality BINARY
	DiurnalOrNot<-rep(0,length(varsToPlot[,1]))
	DiurnalOrNot[which(varsToPlot$Activity123_BM_med==3)]<-1
	names(DiurnalOrNot)<-rownames(varsToPlot)
	varsToPlot<-cbind.data.frame(varsToPlot,DiurnalOrNot)
	rownames(varsToPlot)<-datOrig_ALL$tiplabel

# Get the tip DR values...
###
	# Load DRs
	cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_NDexp_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
	DR_mamPhy <- cladesDR[,"harmMeans"]
	names(DR_mamPhy) <- cladesDR$tiplabel

	# Make the color scale
	## basing this on the DR tip range, not the reconRate range, which is buffered by the BM ancestral recon...
		range <- quantile(DR_mamPhy, seq(0,1, 0.01))[2:101]
		range[100] <- range[100]+0.1
		#cols<-rev(colorRampPalette(c('red','orange','yellow','green','deepskyblue1','blue','black'))(100))
		#library(viridis)
		#cols<-viridis(100)
		cols<-rich.colors(100)

		x.tick <- quantile(DR_mamPhy, c(0.01,0.5,0.99,1))

	# ready tree for plotting...
		plottree<-ladderize(mamMCC)
		DR<-DR_mamPhy

		DR_ordered <- DR[match(plottree$tip.label,names(DR))]
		DR_anc <- ace(DR_ordered, plottree, method="pic", scaled=TRUE)$ace

		# Match ancestors totree edges
	    match.anc <- match(as.numeric(names(DR_anc)), plottree$edge[,2])[-1]

	    # Assign rates to each internal node
	    reconRate <- vector(mode="numeric", length=length(plottree$edge[,2]))
	    reconRate[match.anc] <- DR_anc[2:(length(plottree$tip.label)-1)] #[2:7237]

	    # Assign rates to tips
	    tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

	    reconRate[match(tip.idx, plottree$edge[,2])] <- DR[sort(names(DR))]

	    # Create colour palette
	    reconColors <- reconRate

		range <- quantile(DR, seq(0,1, 0.01))[2:101]

	    for (i in 1:100) {
	        if (i==100) {range[i] <- range[i]+0.1}
	        if (i==1) {reconColors[which(reconRate <= range[i])] <- cols[i] } 
	        if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- cols[i]}
	        }

	   #pdf(file="treeForNewFig1_wDR_w66Ma-KPg.pdf", width=6, height=12)
	   pdf(file="treeForNewFig1_wDR_w66Ma-KPg_w56Ma-PETM.pdf", width=6, height=12)
	   #pdf(file="treeForNewFig1_wDR_smaller.pdf", width=3, height=6)
	    # Plot the tree 
		obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, tip.color="black",#x.lim=c(0,roots[j]),
			label.offset=0.1, type="phylogram", edge.width=0.8, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors))
		axisPhylo( line=-1)
		abline(v=max(branching.times(plottree))-66, lty=2, lwd=1, col="gray")
		abline(v=max(branching.times(plottree))-56, lty=2, lwd=1, col="gray")
		#mtext(side=3,text=Names[j], cex=0.8, font=2)
		#par(new=TRUE, fig=c(0.0, 0.15,0.05,0.25),mar=rep(1,4)) #c(x1, x2, y1, y2)
		#par(fig=c(0.1, 0.9,0,0.9))#,mar=c(0,0,0,0), mar=c(0,0,0,0)) #c(x1, x2, y1, y2)
	#	plot(density(DR), add=TRUE, col="dark grey", main="", bty="n", xlab="", ylab="",axes=F, xlim=range(DR))
	#	polygon(density(DR), col="light grey", border="black", bty="n",main="")
	#	dens.rate <- density(DR)$y
	#	axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1, tck=-0.05, mgp=c(1,1,0))

	#	x.tick_other <- quantile(DR, c(0.01,0.5,0.99,1))
		dev.off()


	# PLOT tree with the BAMM rate-shifts...
	######
	rootAge=max(node.depth.edgelength(plottree))

	# for the BAMM shifts
	toPlotFIN<-read.table(file=paste0("NDexp_BAMMshifts_24shifts_ReadyToPlot.txt"), header=TRUE)

	up<-rgb(1,0,0,0.5) #red
	down<-rgb(0,0,1,0.5) #blue
	upDown<-grey(0.5, alpha = 0.5)
	bgcols<-rep(NA,length(toPlotFIN$shift))

	bgcols[which(toPlotFIN$shift=="upShift")]<-up
	bgcols[which(toPlotFIN$shift=="downShift")]<-down
	bgcols[which(toPlotFIN$shift=="up-or-down")]<-upDown

	bgsizes<-toPlotFIN$factor

	# for the grey order boxes
	orderOfTips<-obj$xy$yy[seq_len(obj$Ntip)]
	tipDetails1<-cbind.data.frame(plottree$tip.label,do.call(rbind,strsplit(plottree$tip.label,split="_")))
	colnames(tipDetails1)<-c("tiplabel","gen","sp","fam","ord")
	toJoin<-cladesDR[,c("tiplabel","clade")]
	tipDetails<-left_join(tipDetails1,toJoin,by="tiplabel")

	blackOrds<-c(#"LAGOMORPHA", #"SCANDENTIA",
				"PRIMATES","CETARTIODACTYLA", "EULIPOTYPHLA",#"PERISSODACTYLA","PHOLIDOTA","AFROSORICIDA","TUBULIDENTATA","HYRACOIDEA","CINGULATA",
				"DIDELPHIMORPHIA")#,"PAUCITUBERCULATA")
	grayOrds<-c("RODENTIA", "CHIROPTERA", "CARNIVORA",#"DERMOPTERA",
				 #"MACROSCELIDEA","PROBOSCIDEA","SIRENIA","PILOSA",
				"Australidelphia")#"DIPROTODONTIA","DASYUROMORPHIA","NOTORYCTEMORPHIA","PERAMELEMORPHIA","MICROBIOTHERIA")#,"MONOTREMATA")
	

	majorOrds<-c("CHIROPTERA", "CARNIVORA","DIPROTODONTIA","DIDELPHIMORPHIA",
				"RODENTIA", "PRIMATES", "CETARTIODACTYLA","EULIPOTYPHLA", "DASYUROMORPHIA")

	#grayCLADES<-c("Mouse-related","Guinea_pig-related","PRIMATES","Yinpterochiroptera","Whippomorpha","EULIPOTYPHLA","Marsupialia")
	grayCLADES<-c("Mouse-related","Guinea_pig-related","Platyrrhini","Strepsirrhini","Yangochiroptera",
				"Ruminantia","Suina","Caniformes",#"Pholidota",
				"Soricidae", "Talpidae", "Marsupialia")


# ALL PLOTTED TOGETHER
#####
library(viridis)
varCols<-magma(20)[c(2,11,7)]
	#latCols<-rgb(col2rgb(magma(20)[7])[1,][[1]], col2rgb(magma(20)[7])[2,][[1]], col2rgb(magma(20)[7])[3,][[1]])
	fc <- colorRampPalette(c("mistyrose", "darkorchid4"))
	plot(rep(1, 30),col = fc(30), pch = 19, cex = 3)
	latCols<-c(fc(30)[10],fc(30)[27])

labSize<-1.5


			pdf(file="newFig1_wDR_wBAMMshifts_wVag_wDiur_wLat_newLabels_latCols.pdf", width=20, height=14)
			# 	par(oma=c(0,0,2,20, xpd=NA) #‘c(bottom, left, top, right)
			# 	layout(matrix(c(1:2), 1, 2, byrow = TRUE), widths=c(8,12), heights=rep(14,4))
			 	layout(matrix(c(1:4), 1, 4, byrow = TRUE), widths=c(8,4,4,4), heights=rep(14,4))

					# plot
					###
			#	   pdf(file="treeForNewFig1_wDR_wBAMMshifts_wlabel_wLeg.pdf", width=8, height=14)
				    par(oma=c(0,0,2,4), mar=c(1,1,1,1)) #‘c(bottom, left, top, right)

					# plot dummy tree
					obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, tip.color="black", #x.lim=c(0,roots[j]),
						label.offset=0.1, type="phylogram", edge.width=0.8, no.margin=TRUE, root.edge=TRUE, edge.color="white")
					# plot CLADES
					for(j in 1:length(grayCLADES)){
						if(j==11){		
							ordTips<-which(tipDetails$ord=="PAUCITUBERCULATA" | tipDetails$ord=="DIDELPHIMORPHIA" | tipDetails$ord=="MICROBIOTHERIA" | tipDetails$ord=="NOTORYCTEMORPHIA" | tipDetails$ord=="PERAMELEMORPHIA" | tipDetails$ord=="DASYUROMORPHIA" | tipDetails$ord=="DIPROTODONTIA")
						} else {
							ordTips<-which(tipDetails$clade==grayCLADES[j])
						}			
						par(new=T, xpd=NA)
						#rect(xleft=0, ybottom=min(orderOfTips[ordTips]), xright=obj$xy$xx[1]+9, ytop=max(orderOfTips[ordTips]), col=gray(0.5, alpha=0.2), border=gray(0.5, alpha=0.2))
						rect(xleft=obj$xy$xx[1], ybottom=min(orderOfTips[ordTips]), xright=obj$xy$xx[1]+250, ytop=max(orderOfTips[ordTips]), col=gray(0.5, alpha=0.2), border=NA)
					}

					# plot ORDERS
					for(j in 1:length(blackOrds)){
							blackOrdTips<-which(tipDetails$ord==blackOrds[j])
						par(new=T, xpd=NA)
						#rect(xleft=0, ybottom=min(orderOfTips[ordTips]), xright=obj$xy$xx[1]+9, ytop=max(orderOfTips[ordTips]), col=gray(0.5, alpha=0.2), border=gray(0.5, alpha=0.2))
						rect(xleft=obj$xy$xx[1], ybottom=min(orderOfTips[blackOrdTips]), xright=obj$xy$xx[1]+7, ytop=max(orderOfTips[blackOrdTips]), col=gray(0.2), border=NA)
						}
					for(j in 1:length(grayOrds)){
						if(j==4){
							grayOrdTips<-which(tipDetails$ord=="DIPROTODONTIA" | tipDetails$ord=="DASYUROMORPHIA" | tipDetails$ord=="NOTORYCTEMORPHIA" | tipDetails$ord=="PERAMELEMORPHIA" | tipDetails$ord=="MICROBIOTHERIA")
						} else {
							grayOrdTips<-which(tipDetails$ord==grayOrds[j])
						}
						par(new=T, xpd=NA)
						#rect(xleft=0, ybottom=min(orderOfTips[ordTips]), xright=obj$xy$xx[1]+9, ytop=max(orderOfTips[ordTips]), col=gray(0.5, alpha=0.2), border=gray(0.5, alpha=0.2))
						rect(xleft=obj$xy$xx[1], ybottom=min(orderOfTips[grayOrdTips]), xright=obj$xy$xx[1]+7, ytop=max(orderOfTips[grayOrdTips]), col=gray(0.7), border=NA)
						}
					
					par(new=T)
				    # Plot the real tree 
					obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, tip.color="black", #x.lim=c(0,roots[j]),
						label.offset=0.1, type="phylogram", edge.width=0.3, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors))
					axis(side=3, line=-2, at=rootAge-c(0,50,100,150), labels=FALSE, cex.axis=labSize) #c(0,50,100,150)
					text(y=c(6120,6120,6140,6140), x=rootAge-c(0,50,100,150), labels=c(0,50,100,150), cex=labSize, srt=270)

					nodelabels(srt=270,pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= 2*bgsizes)
			#		nodelabels(srt=270,text=toPlotFIN$ID[1], node=toPlotFIN$node[1], adj=c(-0.8,0.8), font=2, frame = "none", col = "black", cex= 1.75)
			#		nodelabels(srt=270,text=toPlotFIN$ID[2:33], node=toPlotFIN$node[2:33], adj=c(1.25,1.25), font=2, frame = "none", col = "black", cex= 1.75)
					#nodelabels(text=toPlotFIN$CLADE, node=toPlotFIN$node, adj=c(1.5,2.5), font=1, frame = "none", col = "black", cex= 1.5)

			#		# BAMM legend
			#		text(x=rootAge-130,y=5650, labels="Diversification\nrate shifts", cex=labSize, font=2)
			#		legend(x=rootAge-150,y=5800, legend = c("4.0x", "2.0x", "1.1x"), title="", title.adj=0.5, x.intersp=-4.75, y.intersp=2.2, bty = "n", box.lwd=2, box.col="black",lwd=2, cex=labSize, pt.cex = c(8,4,2), lty=c(NA,NA,NA), pt.bg = "white", pch = c(21, 21, 21)) #c(up, down, upDown)
			#		legend(x=rootAge-130,y=5800, legend = c("up", "down", "up/down"), title="", title.adj=0, x.intersp=0.1, y.intersp=2.2, bty = "n", box.lwd=2, box.col="black",lwd=2, cex=labSize, pt.cex = c(3,3,3), lty=c(NA,NA,NA), pt.bg = c(up, down, upDown), pch = c(22, 22, 22))

	#				dev.off()


	# SEPARATE
	####
	# SEGMENTS STYLE
	############
	#	plottree<-ladderize(mamMCC)
	#	# COLS FOR TRAITS
	#	library(viridis)
	#	varCols<-viridis(10)[c(2,5,8)]

	# VAGILITY
		VAR<-log10(varsToPlot[,"DispDist_kmMAX_final"]*1000)#-min(log(varsToPlot[,"DispDist_kmMAX_final"]*1000))
			# 0 ==  0.00447 km max dispersal distance == 4.5 meters... 1.001499 == exp(min(log(varsToPlot[,"DispDist_kmMAX_final"]*1000))/1000)
			# Melomys_rubicola_MURIDAE_RODENTIA == min
		names(VAR)<-rownames(varsToPlot)
		VAR[is.na(VAR)] <- 0

		treeDat<-treedata(plottree, as.data.frame(VAR))
		trait<-as.vector(treeDat$data)
		names(trait)<-rownames(treeDat$data)
		
		magma20<-magma(20)
		FACTOR<-15
		#pdf(file="mamMCC_withTraits_Vagility_5911sp_ladder_tipNames_nCol_seg_wNames.pdf", width=8, height=140)
		#pdf(file="mamMCC_withTraits_Vagility_5911sp_ladder_tipNames_nCol_seg_log10_wNames.pdf", width=8, height=140)
	#	pdf(file="mamMCC_withTraits_Vagility_5911sp_ladder_tipNames_nCol_seg_log10.pdf", width=8, height=14)
			#    par(oma=c(0,0,2,20), xpd=NA) #‘c(bottom, left, top, right)

		    # Plot the real tree 
			obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.1, tip.color="black", #x.lim=c(0,roots[j]),
				label.offset=0.1, type="phylogram", edge.width=0.8, no.margin=TRUE, root.edge=TRUE, edge.color=NA)#"black")

			# for the TIP order 
			orderOfTips<-obj$xy$yy[seq_len(obj$Ntip)]

#			eqLine <- obj$xy$xx[1]+40
			# wLESS
			LESS <- -180
			eqLine <- obj$xy$xx[1]+40 +LESS

			ticksAt <- c(0, 1, 2, 3, 4, 5, 6, 7)
			ticksLabels <- c("0", round((10^ticksAt[2])/1000,2), round((10^ticksAt[3])/1000,2), round((10^ticksAt[4])/1000,2), round((10^ticksAt[5])/1000,2), round((10^ticksAt[6])/1000,2), round((10^ticksAt[7])/1000,2), round((10^ticksAt[8])/1000,2))
			#ticksLabels <- c("0", "0.005", "0.025", "13.3", "727.5", "39,721")
			#ticksLabels_log <- c("1.5", "9.5", "17.5")

				# > exp((8+1.497388))/1000
				# [1] 13.32488
				# > exp((16+1.497388))/1000
				# [1] 39720.9

			buffer<-150#10#150
			# plot coord axis
			segments( x0=eqLine, y0=1, x1=eqLine, y1=5911+30, lty=1, col=gray(0.2), lwd=1)
				text( x=eqLine, y=5911+buffer, labels=ticksLabels[1], cex=labSize, srt=270)
			
			for(j in 2:length(ticksAt)){
			mod<-ticksAt[j]*FACTOR
			segments( x0=eqLine+mod, y0=1, x1=eqLine+mod, y1=5911+30, lty=2, col=gray(0.5), linklabelswd=0.5)
				text( x=eqLine+mod, y=5911+buffer, labels=ticksLabels[j],  cex=labSize, srt=270)
			}

			# get coords PER TIP
			for(j in 1:length(orderOfTips)){
				tipInTree <- orderOfTips[j]
				species <- plottree$tip.label[j]
				if(varsToPlot[species,"DispDist_kmMAX_final"] <= 1 & varsToPlot[species,"DispDist_kmMAX_final"] > 0.1){
					vagCol<-magma20[17]
				} else if(varsToPlot[species,"DispDist_kmMAX_final"] <= 0.1 ){
					vagCol<-magma20[20]
				} else {
					vagCol<-magma20[13]
				}
				sppDat <- as.data.frame(trait[which(names(trait)==species)])
				valX <- eqLine + sppDat[,1]*FACTOR #"lat_min"
				segments( x0=eqLine, y0=tipInTree, x1=valX, y1=tipInTree, lty=1, col=vagCol, lwd=0.3)
			}

	#	dev.off()



	# DIURNALITY
		VAR<-(varsToPlot[,"DiurnalOrNot"])
		names(VAR)<-rownames(varsToPlot)
		VAR[is.na(VAR)] <- 0

		treeDat<-treedata(plottree, as.data.frame(VAR))
		trait<-as.vector(treeDat$data)
		names(trait)<-rownames(treeDat$data)
		
		#pdf(file="mamMCC_withTraits_Diurnality_5911sp_ladder_tipNames_nCol_seg_wNames.pdf", width=8, height=140)
	#	pdf(file="mamMCC_withTraits_Diurnality_5911sp_ladder_tipNames_nCol_seg.pdf", width=8, height=14)
	#	    par(oma=c(0,0,2,20), xpd=NA) #‘c(bottom, left, top, right)

		    # Plot the real tree 
			obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.1, tip.color="black", #x.lim=c(0,roots[j]),
				label.offset=0.1, type="phylogram", edge.width=0.8, no.margin=TRUE, root.edge=TRUE, edge.color=NA)#"black")

			# for the TIP order 
			orderOfTips<-obj$xy$yy[seq_len(obj$Ntip)]
	#		eqLine <- obj$xy$xx[1]+10 
	#		# wEXTRA
	#		EXTRA <- 145
	#		eqLine <- obj$xy$xx[1]+10 +EXTRA
			# wLESS
			LESS <- -210
			eqLine <- obj$xy$xx[1]+40 +LESS

			#buffer<-150#10#150
			## plot coord axis
			#segments( x0=eqLine, y0=1, x1=eqLine, y1=5911+30, lty=1, col=gray(0.2), lwd=1)
				text( x=eqLine, y=5911+buffer, labels="0", cex=labSize, srt=270)
			#segments( x0=eqLine+50, y0=1, x1=eqLine+50, y1=5911+30, lty=1, col=gray(0.2), lwd=0.5)
				text( x=eqLine+30, y=5911+buffer, labels="1", cex=labSize, srt=270)

			# get coords PER TIP
			for(j in 1:length(orderOfTips)){
				tipInTree <- orderOfTips[j]
				species <- plottree$tip.label[j]
				sppDat <- as.data.frame(trait[which(names(trait)==species)])
				valX <- eqLine + sppDat[,1]*30 #"lat_min"
				segments( x0=eqLine, y0=tipInTree, x1=valX, y1=tipInTree, lty=1, col=magma20[11], lwd=0.5) #"darkorange"
			}

	#	dev.off()

# LATITUDE 
	VAR<-varsToPlot[,c("lat_min","lat_max")]
	VAR[is.na(VAR[,1]),1] <- 0
	VAR[is.na(VAR[,1]),2] <- 0

	treeDat<-treedata(plottree, as.data.frame(VAR))
	trait<-treeDat$data

	#pdf(file="mamMCC_withTraits_latMax_5911sp_ladder_tipNames_nCol_wNames.pdf", width=8, height=140)
#	pdf(file="mamMCC_withTraits_latMax_5911sp_ladder_tipNames_nCol_latCol.pdf", width=8, height=14)
#	    par(oma=c(0,0,2,20), xpd=NA) #‘c(bottom, left, top, right)

	    # Plot the real tree 
		obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.1, tip.color="black", #x.lim=c(0,roots[j]),
			label.offset=0.1, type="phylogram", edge.width=0.8, no.margin=TRUE, root.edge=TRUE, edge.color=NA)

		# for the TIP order 
		orderOfTips<-obj$xy$yy[seq_len(obj$Ntip)]
#		eqLine <- obj$xy$xx[1]+90
#			# wEXTRA
#			EXTRA <- 170
#			eqLine <- obj$xy$xx[1]+90 +EXTRA
			# wLESS
			LESS <- -300
			eqLine <- obj$xy$xx[1]+90 +LESS

		buffer<-150#10#150
		# plot coord axis
		segments( x0=eqLine, y0=1, x1=eqLine, y1=5911+30, lty=1, col=gray(0.2), lwd=1)
			text( x=eqLine, y=5911+buffer, labels="0", cex=labSize, srt=270)
		segments( x0=eqLine-23.5, y0=1, x1=eqLine-23.5, y1=5911+30, lty=2, col=gray(0.5), lwd=0.5)
			text( x=eqLine-23.5, y=5911+buffer, labels="-23.5", cex=labSize, srt=270)
		segments( x0=eqLine+23.5, y0=1, x1=eqLine+23.5, y1=5911+30, lty=2, col=gray(0.5), lwd=0.5)
			text( x=eqLine+23.5, y=5911+buffer, labels="23.5", cex=labSize, srt=270)
		segments( x0=eqLine-60, y0=1, x1=eqLine-60, y1=5911+30, lty=2, col=gray(0.5), lwd=0.5)
			text( x=eqLine-60, y=5911+buffer, labels="-60", cex=labSize, srt=270)
		segments( x0=eqLine+60, y0=1, x1=eqLine+60, y1=5911+30, lty=2, col=gray(0.5), lwd=0.5)
			text( x=eqLine+60, y=5911+buffer, labels="60", cex=labSize, srt=270)

		# get coords PER TIP
		for(j in 1:length(orderOfTips)){
			tipInTree <- orderOfTips[j]
			species <- plottree$tip.label[j]

			if(is.na(varsToPlot[species,"lat_min"])){
				next
			} else if(varsToPlot[species,"lat_min"] > -23.5 & varsToPlot[species,"lat_max"] < 23.5){
				col<-"violet" #"magenta" #latCols[1] #"plum1" #magma(20, alpha=0.3)[7]
			} else {
				col<-magma(20)[7] #"darkorchid3" #latCols[2] #"darkorchid3" #magma(20)[7]
			}
			sppDat <- as.data.frame(trait[which(rownames(trait)==species),])
			minX <- eqLine + sppDat[,1][1] #"lat_min"
			maxX <- eqLine + sppDat[,1][2] #"lat_max"
			segments( x0=minX, y0=tipInTree, x1=maxX, y1=tipInTree, lty=1, col=col, lwd=0.5)
		}

	dev.off()


plot(1:20, col=magma(20), pch=20, cex=5)

plot(1:10, col=viridis(10), pch=20, cex=5)




# BARPLOT STYLE
############

# COLS FOR TRAITS
library(viridis)
varCols<-viridis(10)[c(2,5,8)]

# VAGILITY
	VAR<-log(varsToPlot[,"DispDist_kmMAX_final"]*1000)-min(VAR)
	names(VAR)<-rownames(varsToPlot)
	VAR[is.na(VAR)] <- 0

	treeDat<-treedata(plottree, as.data.frame(VAR))
	trait<-as.vector(treeDat$data)
	names(trait)<-rownames(treeDat$data)
	
	pdf(file="mamMCC_withTraits_Vagility_5911sp_ladder_tipNames_nCol.pdf", width=8, height=14)
	plotTree.wBars(tree=(treeDat$phy), x=trait, edge.color="black",
	               scale=4, type="phylogram", method="plotTree", 
	               tip.labels=FALSE, fsize=0.05, border=varCols[1], col=varCols[1], lwd=0.4)
	dev.off()


# DIURNALITY
	VAR<-(varsToPlot[,"DiurnalOrNot"])
	names(VAR)<-rownames(varsToPlot)
	VAR[is.na(VAR)] <- 0

	treeDat<-treedata(plottree, as.data.frame(VAR))
	trait<-as.vector(treeDat$data)
	names(trait)<-rownames(treeDat$data)
	
	pdf(file="mamMCC_withTraits_Diurnality_5911sp_ladder_tipNames_nCol.pdf", width=8, height=14)
	plotTree.wBars(tree=(treeDat$phy), x=trait, 
	               scale=4, type="phylogram", method="plotTree", 
	               tip.labels=FALSE, fsize=0.05, border=varCols[2], col=varCols[2], lwd=0.4)
	dev.off()





