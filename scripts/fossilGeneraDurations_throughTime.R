# New Fossil-genus-through-time analysis
#####
library(plotrix)
setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-geo-trait-analyses/PBDB_FOSSIL_DATA")

# load data
allDat <- read.csv("pbdb_Mammalia_occ-lumpByGenus_regularTaxa_Aug2018_mod.csv", header=TRUE, skip=21)

# subset data
vars<-c("accepted_name","accepted_attr", "accepted_rank", "accepted_no", "early_interval",	"late_interval",	"max_ma",	"min_ma",	"ref_author",	"ref_pubyr",	"reference_no",	
	"phylum",	"class",	"order",	"family",	"genus", "lng",	"lat", "paleolng",	"paleolat","taxon_environment",	"motility",	"life_habit",	"vision",	"diet",	"reproduction")

readyDat <- allDat[which(allDat$max_ma <= 130),vars]

# COLLAPSE TO GENERA (as native, is indiv occurrences of a genus with an error bar for geological age)
###
genNames<-names(table(readyDat$genus))
datGen <- vector("list",length(genNames))

for(j in 1:length(genNames)){
	genDat<-readyDat[which(readyDat$genus==genNames[j]),]
	maxAll<- max(genDat$max_ma)
	minAll<- min(genDat$min_ma)
	intLength <- maxAll - minAll
	midAll <- maxAll - (intLength/2)

	datGen[[j]] <- cbind.data.frame(genDat[1,c("order","family","genus")],maxAll,minAll,midAll,intLength )
}
datGen_ALL<-do.call(rbind,datGen)

# calc K-Pg ending vs. starting vs. crossing
###
# GENERA
ending_KPg<-rep(0,length(datGen_ALL[,1]))
ending_KPg[ datGen_ALL$minAll==66 ]<-1
	# ending_KPg
	#    0    1 
	# 5368   60 
starting_KPg<-rep(0,length(datGen_ALL[,1]))
starting_KPg[ datGen_ALL$maxAll==66 ]<-1
	# starting_KPg
	#    0    1 
	# 5323  105 
crossing_KPg<-rep(0,length(datGen_ALL[,1]))
crossing_KPg[ datGen_ALL$maxAll > 66 & datGen_ALL$minAll < 66]<-1
	# crossing_KPg
	#    0    1 
	# 5420    8
pre_KPg<-rep(0,length(datGen_ALL[,1]))
pre_KPg[ datGen_ALL$maxAll > 66 & datGen_ALL$minAll >= 66]<-1
	# pre_KPg
	#    0    1 
	# 5071  357 
post_KPg<-rep(0,length(datGen_ALL[,1]))
post_KPg[ datGen_ALL$maxAll <= 66 & datGen_ALL$minAll < 66]<-1
	# post_KPg
	#    0    1 
	#  365 5063 
	# 357 + 5063 + 8 == 5428 == all genera in dataset.

# FAMILIES
famNames<-names(table(datGen_ALL$family))[-1]
famCross_KPg<-rep(0,length(famNames))
names(famCross_KPg)<-famNames
for(j in 1:length(famNames)){
	famDat<-datGen_ALL[which(datGen_ALL$family==famNames[j]),]
	maxMax<- max(famDat$maxAll)
	minMin<- min(famDat$minAll)

	if( maxMax > 66 & minMin < 66 ){
		famCross_KPg[j]<-1
	} else { next }
}
	# famCross_KPg
	#   0   1 
	# 539  21 
# get those 21 fam names, translate back to genus records
famNames_thatCrossed<-names(famCross_KPg[famCross_KPg==1])
genFamCross_KPg<-rep(0,length(datGen_ALL[,1]))
genFamCross_KPg[ datGen_ALL$family %in% famNames_thatCrossed ]<-1
	# genFamCross_KPg
	#    0    1 
	# 5264  164 
# ORDERS
ordNames<-names(table(datGen_ALL$order))[-1]
ordCross_KPg<-rep(0,length(ordNames))
names(ordCross_KPg)<-ordNames
for(j in 1:length(ordNames)){
	ordDat<-datGen_ALL[which(datGen_ALL$order==ordNames[j]),]
	maxMax<- max(ordDat$maxAll)
	minMin<- min(ordDat$minAll)

	if( maxMax > 66 & minMin < 66 ){
		ordCross_KPg[j]<-1
	} else { next }
}
	# ordCross_KPg
	#  0  1 
	# 54 9
# get those 9 order names, translate back to genus records
ordNames_thatCrossed<-names(ordCross_KPg[ordCross_KPg==1])
genOrdCross_KPg<-rep(0,length(datGen_ALL[,1]))
genOrdCross_KPg[ datGen_ALL$order %in% ordNames_thatCrossed ]<-1
	# genOrdCross_KPg
	#    0    1 
	# 5087  341 


# join all
###
genALL <- cbind.data.frame(datGen_ALL, ending_KPg, starting_KPg, crossing_KPg, pre_KPg, post_KPg,genFamCross_KPg, genOrdCross_KPg)
write.csv(genALL, file="fossilGen_throughTime_5428gen_Mar2020.csv")

# READ BACK IN (edited to add higher taxa)
###
genALL_wHigher2<-read.csv("fossilGen_throughTime_5428gen_Mar2020_wHigher.csv", header=TRUE)
	# exclude incertaeCedis
	genALL_wHigher <- genALL_wHigher2[which(genALL_wHigher2$HIGHER_1!="incertaeCedis" & genALL_wHigher2$HIGHER_1!="stem Mammalia"),]

	# Clades...
	###
	cladeNames1 <- names(table(genALL_wHigher$HIGHER_1))
	cladeNamesOrdered1 <- cladeNames1[c(6,10,9,5,7,1,11,4,2)]
		cladeNames_index1 <- rep(NA,length(genALL_wHigher[,1]))
		for(j in 1:length(cladeNames_index1)){
			cladeNames_index1[which(genALL_wHigher[,"HIGHER_1"]==cladeNamesOrdered1[j])]<-j
		}

	cladeNames2 <- names(table(genALL_wHigher$HIGHER_2))
	cladeNamesOrdered2 <- cladeNames2[c(3,8,7,2,5,4)]
		cladeNames_index2 <- rep(NA,length(genALL_wHigher[,1]))
		for(j in 1:length(cladeNames_index2)){
			cladeNames_index2[which(genALL_wHigher[,"HIGHER_2"]==cladeNamesOrdered2[j])]<-j
		}

	# extant Orders...
	###
	dat_ordsExtant <- genALL_wHigher[which(genALL_wHigher$'extinctOrder.'==0),]
	ordNames<-names(sort(table(as.vector(dat_ordsExtant$order)), decreasing=TRUE)[c(1:7,11)])
		ordNames_index <- rep(NA,length(genALL_wHigher[,1]))
		for(j in 1:length(ordNames)){
			ordNames_index[which(genALL_wHigher[,"order"]==ordNames[j])]<-j
		}

	# extinct Orders...
	###
	dat_ordsExtinct <- genALL_wHigher[which(genALL_wHigher$'extinctOrder.'==1),]
	ordNames_ext<-names(sort(table(as.vector(dat_ordsExtinct$order)), decreasing=TRUE)[c(2,6)])
		ordNames_index_ext <- rep(NA,length(genALL_wHigher[,1]))
		for(j in 1:length(ordNames_ext)){
			ordNames_index_ext[which(genALL_wHigher[,"order"]==ordNames_ext[j])]<-j
		}

	# join together
	genALL<-cbind.data.frame( cladeNames_index1, cladeNames_index2, ordNames_index, ordNames_index_ext, genALL_wHigher)


# SORT
###
# by timeInt
#index<-order(genALL$intLength, decreasing=TRUE)
# by ending_KPg
#index<-order(genALL$ending_KPg, decreasing=TRUE)
# by timeInt and ending_KPg
#index<-order(genALL$intLength,genALL$ending_KPg, decreasing=TRUE)
# by ending_KPg and timeInt
#index<-order(genALL$ending_KPg,genALL$intLength, decreasing=TRUE)
# by ending_KPg, starting_KPg, and timeInt
#index<-order(genALL$ending_KPg,genALL$starting_KPg, genALL$crossing_KPg, genALL$intLength, decreasing=TRUE)
# by pre_KPg, crossing, then post
#index<-order(genALL$pre_KPg, genALL$crossing_KPg, genALL$post_KPg, genALL$intLength, decreasing=TRUE)
# by pre_KPg, then post
#index<-order(genALL$pre_KPg, genALL$post_KPg, genALL$intLength, decreasing=TRUE)
# by ORDER, pre_KPg, then post
#index<-order(genALL$order, genALL$pre_KPg, genALL$post_KPg, genALL$intLength, decreasing=TRUE)
# by CLADE NAME, pre_KPg, then post, then order
index<-order(genALL$cladeNames_index1, 
	#genALL$pre_KPg, #genALL$post_KPg, 
	genALL$order, genALL$intLength, #genALL$family, 
	decreasing=FALSE)

# plot data
###
library(viridis)

nrows<-length(genALL[,1])
dat<-genALL[index,]

colsToUse<-c("black","grey", viridis(7)[1:6],"goldenrod1") #magma(10)[c(1,3,5,7,9)], 
#colsToUse<-c(viridis(7)[1:6],"goldenrod1")
lineCols<-rep("black",nrows)
lineCols[dat$cladeNames_index1]
for(j in 1:length(cladeNamesOrdered1)){
	lineCols[which(dat[,"HIGHER_1"]==cladeNamesOrdered1[j])]<-colsToUse[j]
}
#labFonts<-c(rep(1,3),rep(2,7))
labFonts<-c(rep(2,9))
#fontColsToUse<-c("black","grey",viridis(7)[1:6],"goldenrod1")
fontColsToUse<-colsToUse

placental_colsToUse<-c(rep(c("black",grey(0.6)),2), grey(0.6), "black", grey(0.6), "black") #viridis(8)[c(1,3,5,7)]


#pdf(file="plot_fossilGen_throughTime_byIntLength-EndingKPg.pdf", height=10, width=6)
#pdf(file="plot_fossilGen_throughTime_byEndingStartingCrossingKPg-IntLength.pdf", height=10, width=6)
#pdf(file="plot_fossilGen_throughTime_byCLADE_PrePostKPg-order.pdf", height=10, width=6)
#pdf(file="plot_fossilGen_throughTime_byCLADE_order-family_wide.pdf", height=6, width=10)
#pdf(file="plot_fossilGen_throughTime_byCLADE_order-family_wLabels.pdf", height=50, width=6)
#pdf(file="plot_fossilGen_throughTime_byCLADE_PrePostKPg_wLabels.pdf", height=10, width=6)
#pdf(file="plot_fossilGen_throughTime_byCLADE_order-family_wColors.pdf", height=10, width=10)
#pdf(file="plot_fossilGen_throughTime_byCLADE_order-family_wColors_7cats_WIDE.pdf", height=8, width=12)
pdf(file="plot_fossilGen_throughTime_byCLADE_order-family_wColors_7cats_byOrd-Int_agesFixed.pdf", height=10, width=8)
#pdf(file="plot_fossilGen_throughTime_byCLADE_PreKPg_wColors.pdf", height=10, width=10)
#pdf(file="plot_fossilGen_throughTime_PrePostKPg.pdf", height=10, width=10)
#quartz(height=50,width=6)
	plot(c(1:nrows)~seq(50,-120,length.out=nrows), col=NA, yaxt="n",xaxt="n",bty="n", xlab="",ylab="", xpd=NA)
	axis(side=1, at=seq(0,-120, by=-20), labels=seq(0,-120, by=-20))
	abline(v= -66, lty=2, col="grey", lwd=4)

	#plot(seq(50,-250,length.out=nrows)~c(1:nrows), col=NA, xaxt="n",yaxt="n",bty="n", xlab="",ylab="")
	#axis(side=2, at=seq(0,-200, by=-50), labels=seq(0,-200, by=-50))
	#abline(h= -66, lty=2, col="grey", lwd=4)

#plotCI(x=c(1:nrows), y = -dat$midAll, ui= -dat$minAll, li= -dat$maxAll, err="y",
plotCI(y=c(1:nrows), x = -dat$midAll, ui= -dat$minAll, li= -dat$maxAll, err="x",
       sfrac=0, gap=0, add=TRUE, col= lineCols, scol=lineCols, pch=20, lwd=1.5, cex=0.0000001)

for(j in 1:length(cladeNamesOrdered1)){
	rowsOfClade <- c(1:nrows)[dat$cladeNames_index1==j]
#	segments(x0=j, x1 = j, y0= min(rowsOfClade), y1 = max(rowsOfClade),
#              col = "blue", lty = 1, lwd = 1)
	if(j %in% c(2, 3, 5, 8:9)){
    text(x= -120, y = min(rowsOfClade)+((max(rowsOfClade)-min(rowsOfClade))/2), srt=0,adj=0, cex=1.1, font=labFonts[j], 
    	labels = cladeNamesOrdered1[j], col=fontColsToUse[j])#, adj = NULL,
	} else {
    text(x= 2, y = min(rowsOfClade)+((max(rowsOfClade)-min(rowsOfClade))/2), srt=0,adj=0, cex=1.1, font=labFonts[j], 
    	labels = cladeNamesOrdered1[j], col=fontColsToUse[j])#, adj = NULL,		
	}

}

# label placental orders
for(j in 1:length(ordNames)){
	rowsOfClade <- c(1:nrows)[which(dat$ordNames_index==j)]
	segments(x0= 2, x1 = 2, y0= min(rowsOfClade)+10, y1 = max(rowsOfClade)-10,
              col = placental_colsToUse[j], lty = 1, lwd = 2)
    text(x= 4, y = min(rowsOfClade)+((max(rowsOfClade)-min(rowsOfClade))/2), srt=0,adj=0, cex=1.1, font=1, 
    	labels = ordNames[j], col=placental_colsToUse[j])#, adj = NULL,
}

# label extinct taxa
for(j in 2:length(ordNames_ext)){
	rowsOfClade <- c(1:nrows)[which(dat$ordNames_index_ext==j)]
	segments(x0= -75, x1 = -75, y0= min(rowsOfClade)+10, y1 = max(rowsOfClade)-10,
              col = grey(0.6), lty = 1, lwd = 2)
    text(x= -77, y = min(rowsOfClade)+((max(rowsOfClade)-min(rowsOfClade))/2), srt=0,adj=1, cex=0.9, font=1, 
    	labels = ordNames_ext[j], col=grey(0.6))#, adj = NULL,
}


dev.off()








# *********************
# IF DEALING WITH ALL OCCURRENCES...
#####

# calc extra vars
###
intLength <- (readyDat$max_ma - readyDat$min_ma)
mid_ma <- readyDat$max_ma - (intLength/2)

# calc K-Pg ending vs. starting
###
ending_KPg<-rep(0,length(readyDat[,1]))
ending_KPg[ readyDat$min_ma==66 ]<-1

starting_KPg<-rep(0,length(readyDat[,1]))
starting_KPg[ readyDat$max_ma==66 ]<-1

# calc K-Pg crossing -- more complex, need to collapse to families...
###
crossing_KPg<-rep(0,length(readyDat[,1]))
crossing_KPg[ readyDat$max_ma > 66 & readyDat$min_ma < 66]<-1
	# NONE == 0, because no genera are known to cross the K-Pg!

famNames<-names(table(readyDat$family))[-1]
famCross_KPg<-rep(0,length(famNames))
names(famCross_KPg)<-famNames
for(j in 1:length(famNames)){
	famDat<-readyDat[which(readyDat$family==famNames[j]),]
	maxMax<- max(famDat$max_ma)
	minMin<- min(famDat$min_ma)

	if( maxMax > 66 & minMin < 66 ){
		famCross_KPg[j]<-1
	} else { next }
}
	# famCross_KPg
	#   0   1 
	# 534  26 
# get those 26 fam names, translate back to genus records
famNames_thatCrossed<-names(famCross_KPg[famCross_KPg==1])
genFamCross_KPg<-rep(0,length(readyDat[,1]))
genFamCross_KPg[ readyDat$family %in% famNames_thatCrossed ]<-1

# same for ORDERS
ordNames<-names(table(readyDat$order))[-1]
ordCross_KPg<-rep(0,length(ordNames))
names(ordCross_KPg)<-ordNames
for(j in 1:length(ordNames)){
	ordDat<-readyDat[which(readyDat$order==ordNames[j]),]
	maxMax<- max(ordDat$max_ma)
	minMin<- min(ordDat$min_ma)

	if( maxMax > 66 & minMin < 66 ){
		ordCross_KPg[j]<-1
	} else { next }
}
	# ordCross_KPg
	#  0  1 
	# 54 10 
# get those 10 ord names, translate back to genus records
ordNames_thatCrossed<-names(ordCross_KPg[ordCross_KPg==1])
genOrdCross_KPg<-rep(0,length(readyDat[,1]))
genOrdCross_KPg[ readyDat$order %in% ordNames_thatCrossed ]<-1

# join all
###
timeTaxoDat<- cbind.data.frame(readyDat[,c(12:16,5:8)], mid_ma, intLength, genFamCross_KPg, genOrdCross_KPg)

# SORT
###
# by timeInt
index<-order(timeTaxoDat$intLength, decreasing=TRUE)

# plot data
###
nrows<-length(timeTaxoDat[,1])
dat<-timeTaxoDat[index,]

quartz(height=10,width=6)
plot(c(1:nrows)~seq(0,-250,length.out=nrows), col=NA, yaxt="n",bty="n", xlab="",ylab="")

plotCI(x = -dat$mid_ma, y=c(1:nrows), ui= -dat$max_ma, li= -dat$min_ma, err="x",
       sfrac=0, gap=0, add=TRUE, col= "black", scol="black", pch=20, lwd=2, cex=0.1)






