# Code to assemble & summarize the fossil record of mammals from the PBDB 
# and from the literature (Pires et al. 2018 / Grossnickle & Newham 2016)
#########
# Fig 2
########
# NS Upham; last modification: 6 April 2021
# ===================================================

setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-geo-trait-analyses/PBDB_FOSSIL_DATA")
source("sqs_function.R")
library(paleotree); library(divDyn)
library(dplyr)

# read in PBDB (lightly) edited data
#########
PBDB_all <- read.csv("pbdb_Mammalia_occ-lumpByGenus_regularTaxa_Aug2018_mod.csv", header=TRUE, skip=21)
    # 75,286 observations

    # subset data
    vars<-c("accepted_name", "max_ma", "min_ma")
    PBDB_dat <- PBDB_all[,vars]

# read in Pires et al. 2018 data on K-Pg mammals
#########
Pires_all<- read.csv("PiresEtAl2018_KPg_fossilMammals.csv", header=TRUE)

    # subset data
    vars<-c("genus", "ma_max", "ma_min")
    Pires_dat <- Pires_all[,vars]
    colnames(Pires_dat)<-c("accepted_name", "max_ma", "min_ma")

# which genera in Pires et al. are also in the PBDB? (want to use Pires et al. ages preferentially)
###
    PiresNames<-names(table(Pires_dat$accepted_name))
    PbdbNames<-names(table(PBDB_dat$accepted_name))
        PiresNames_in_PBDB<-PbdbNames[PbdbNames %in% PiresNames] #as.numeric(na.omit(match(PiresNames, PbdbNames)))]
            # 193 genera in Pires overlap with PBDB (and 96 genera do not!)
        allOther_PbdbNames<-setdiff(PbdbNames,PiresNames_in_PBDB)

    # exclude the PiresNames from the PBDB set
    PBDB_sub<-PBDB_dat[PBDB_dat$accepted_name %in% allOther_PbdbNames,]
        # removes 5068 records of 193 genera that overlap btwn PBDB and Pires et al.
    # sort these
    PBDB_subSorted<-PBDB_sub[order(PBDB_sub$accepted_name),]

    # add Pires et al. set (193 + 96 genera) to this subsetted PBDB set:
    readyDat_all<-rbind(PBDB_subSorted,Pires_dat)
        # now has 72,888 observations (adds the 2670 obs from Pires et al.)
        # now has 5,480 genera (adds 289 genera from Pires et al.)

    # subset to < 145 Ma to focus analyses on the Cretaceous-Recent
    readyDat<-readyDat_all[which(readyDat_all$max_ma < 131),]
        # now has 72,579 observations (excludes 309 observations of 160 genera)
        # now has 5,320 genera

    write.csv(readyDat, file="pbdb_Mammalia_wPiresEtAl2018_72579obs.csv")


    # read in Grossnickle & Newham 2016 data on K-Pg mammals
    #########
    #Gross_all<- read.csv("Grossnickle&Newham2016_FossilOccurrenceDataset_mod.csv", header=TRUE)
        # chose not to use this dataset, aim is to focus on the Pires et al. contributions


# Do manual binning, according to 0/1 in a given time bin.
#######
    readyToBin <- cbind.data.frame(FAD=readyDat$max_ma, LAD=readyDat$min_ma)

    ######
    # 10-Ma bins
    ###
  #  bins_10Ma<- seq(from=0,to=130,by= 10)
  #  binNames_10Ma<-paste0("bin_",bins_10Ma)[-14]
    ###
    # 5-Ma bins (from ZERO Ma)
    ###
  #  bins_5Ma<- seq(from=0,to=130,by= 5)
  #  binNames_5Ma<-paste0("bin_",bins_5Ma)[-27]
    ###
    # 5-Ma bins (from ONE Ma)
    ###
    bins_5Ma<- seq(from=1,to=131,by= 5)
    binNames_5Ma<-paste0("bin_",bins_5Ma,"_",bins_5Ma[-1])[-length(bins_5Ma)]
    # load STAGES
  #  stages<-read.csv(file="stagesThroughTime.csv", header=TRUE)
  #  bins_stages<-stages$max_ma[-c(36:48)]
  #  binNames_stages<-as.character(stages$interval_name)[-c(36:48)]


    bins <- bins_5Ma
    binNames <- binNames_5Ma
  #  bins <- bins_stages
  #  binNames <- binNames_stages

    binDat<-matrix(
        data=rep(0,(length(binNames)*length(readyDat[,1]))),
        nrow=length(readyDat[,1]),
        ncol=length(binNames),
        dimnames=list(1:length(readyDat[,1]),binNames)
        )

    binNames_j_ALL<-list()
    for(j in 1:length(binDat[,1])){

        genus_j<-as.character(readyDat[j,"accepted_name"])
        dat_j<-as.numeric(readyDat[j,c("max_ma","min_ma")])

        binNum_j<-as.numeric(cut(dat_j, breaks=bins, include.lowest=TRUE,  right=FALSE, labels=binNames))
        if( is.na(binNum_j[1]) | is.na(binNum_j[2]) ) { next } else {

        if(length(binNum_j) > 1){

            if( max(binNum_j)-min(binNum_j) == 0 ){
                binNum_all_j <- binNum_j[1]
            } else if(
                max(binNum_j)-min(binNum_j) == 1 ){
                binNum_all_j <- binNum_j
            } else if( 
                max(binNum_j)-min(binNum_j) == 2 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1)
            } else if( 
                max(binNum_j)-min(binNum_j) == 3 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2)
            } else if( 
                max(binNum_j)-min(binNum_j) == 4 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2, max(binNum_j)-3)
            } else if( 
                max(binNum_j)-min(binNum_j) == 5 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2, max(binNum_j)-3, max(binNum_j)-4)
            } else if( 
                max(binNum_j)-min(binNum_j) == 6 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2, max(binNum_j)-3, max(binNum_j)-4, max(binNum_j)-5)
            } else if( 
                max(binNum_j)-min(binNum_j) == 7 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2, max(binNum_j)-3, max(binNum_j)-4, max(binNum_j)-5, max(binNum_j)-6)
            } else if( 
                max(binNum_j)-min(binNum_j) == 8 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2, max(binNum_j)-3, max(binNum_j)-4, max(binNum_j)-5, max(binNum_j)-6, max(binNum_j)-7)
            } else if( 
                max(binNum_j)-min(binNum_j) == 9 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2, max(binNum_j)-3, max(binNum_j)-4, max(binNum_j)-5, max(binNum_j)-6, max(binNum_j)-7, max(binNum_j)-8)
            } else if( 
                max(binNum_j)-min(binNum_j) == 10 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2, max(binNum_j)-3, max(binNum_j)-4, max(binNum_j)-5, max(binNum_j)-6, max(binNum_j)-7, max(binNum_j)-8, max(binNum_j)-9)
            } else if( 
                max(binNum_j)-min(binNum_j) == 11 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2, max(binNum_j)-3, max(binNum_j)-4, max(binNum_j)-5, max(binNum_j)-6, max(binNum_j)-7, max(binNum_j)-8, max(binNum_j)-9, max(binNum_j)-10)
            } else if( 
                max(binNum_j)-min(binNum_j) == 12 ){
                binNum_all_j <- c(binNum_j, max(binNum_j)-1, max(binNum_j)-2, max(binNum_j)-3, max(binNum_j)-4, max(binNum_j)-5, max(binNum_j)-6, max(binNum_j)-7, max(binNum_j)-8, max(binNum_j)-9, max(binNum_j)-10, max(binNum_j)-11)
            }
        } else {
            binNum_all_j <- binNum_j
        }
        }

        binNames_j<-binNames[binNum_all_j]

        binNames_j_list<-list()
        for(k in 1:length(binNames_j)){
            binDat[j, binNum_all_j[k]]<-1
            binNames_j_list[[k]]<-c(genus_j,binNames_j[k], binNum_all_j[k])
        }
        binNames_j_ALL[[j]]<-do.call(rbind,binNames_j_list)
    }
    genBins_ALL<-as.data.frame(do.call(rbind,binNames_j_ALL))
    colnames(genBins_ALL)<-c("genus","binned", "binNum")
        # NEW (w Pires et al. and starting at 1 Ma): yields 90,548 genus occs in a bin

    write.csv(genBins_ALL, file="pbdb_Mammalia_wPiresEtAl2018_binned1to131.csv")


            # HOW MANY of those genera in the 0-5 Ma bin are EXTANT?
            #######
            recentGen<-names(table(genBins_ALL[which(genBins_ALL$binned=="bin_1_6"),1]))
            mamTax<-read.csv("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL/taxonomy_mamPhy_5911species.csv")
                mamTax<-read.csv("/Users/nate/Desktop/ASM_Service/ASM_BiodivComm/MDD_versions/v1.4_11April2021_11April2021_11April2021/MDD/MDD_v1.4_6533species.csv")

            mamGen<-names(table(mamTax$gen))
            table(recentGen %in% mamGen) 
             #   FALSE  TRUE 
             #    1150   660 #<< previous
             #   FALSE  TRUE 
             #   1141   351  #<< new w Pires and 1 Ma start
             #       FALSE  TRUE 
             #        1149   343 #< MDD taxonomy -- new w Pires and 1 Ma start
            table(mamGen %in% recentGen)
             #   FALSE  TRUE 
             #     623   660 #<< old: so then 623 additional extant genera are without fossil records.
             #           # means that 51% (51.44) of extant genera have a fossil record at all...
             #   FALSE  TRUE 
             #     932   351 #<< new: so then 932 additional extant genera are without fossil records.
             #           # means that 27% (27.35) of extant genera have a fossil record older than 1 Million years
             #       FALSE  TRUE 
             #         989   343  #< MDD taxonomy means that 26% (25.78) of extant genera have a fossil record older than 1 Million years

# Reduce to GENUS-LEVEL fossil occurrences
######
# For subsequent plotting / visualization of the whole fossil record
#########
genNames<-names(table(readyDat$accepted_name))
datGen <- vector("list",length(genNames))

for(j in 1:length(genNames)){
    genDat<-readyDat[which(readyDat$accepted_name==genNames[j]),]
    maxAll<- max(genDat$max_ma)
    minAll<- min(genDat$min_ma)
    intLength <- maxAll - minAll
    midAll <- maxAll - (intLength/2)

    datGen[[j]] <- cbind.data.frame(accepted_name=genDat[1,"accepted_name"],maxAll,minAll,midAll,intLength )
}
datGen_ALL<-do.call(rbind,datGen)

write.csv(datGen_ALL, file="fossilGen_throughTime_5320gen_wPiresEtAl_Apr2021.csv")

# Now *match* this broader taxonomy:
####
    # load taxonomy (hand curated based on 16 Aug 2018 download of the PBDB)
    fossilTaxo<-read.csv("fossilGen_throughTime_5428gen_Mar2020_justTaxonomy.csv")

    # merge with newly updated fossil ranges
    library(dplyr)

    datGen_ALL_wTax<-right_join(x=fossilTaxo, y=datGen_ALL, by=c("genus"="accepted_name"))
        unmatched<-anti_join(x=datGen_ALL, y=fossilTaxo, by=c("accepted_name"="genus"))

    write.csv(datGen_ALL_wTax, file="fossilGen_throughTime_5320gen_wPiresEtAl_Apr2021_wTax.csv")

# read back in the * cleaned * taxonomy file
    genALL_wHigher2<-read.csv("fossilGen_throughTime_5320gen_wPiresEtAl_Apr2021_wTax_cleaned.csv", header=TRUE)
    # exclude incertaeCedis
    genALL_wHigher <- genALL_wHigher2[which(genALL_wHigher2$HIGHER_1!="incertaeCedis" & genALL_wHigher2$HIGHER_1!="stem Mammalia"),]

#    # OLD DATA FOR REFERENCE
#    genALL_wHigher2<-read.csv("fossilGen_throughTime_5428gen_Mar2020_wHigher.csv", header=TRUE)
#        # exclude incertaeCedis
#        genALL_wHigher <- genALL_wHigher2[which(genALL_wHigher2$HIGHER_1!="incertaeCedis" & genALL_wHigher2$HIGHER_1!="stem Mammalia"),]


    # Clades...
    ###
    cladeNames1 <- names(table(genALL_wHigher$HIGHER_1))
    #cladeNamesOrdered1 <- cladeNames1[c(6,10,9,5,7,1,11,4,2)]
    cladeNamesOrdered1 <- cladeNames1[c(5,8,7,4,6,1,9,3,2)]
        cladeNames_index1 <- rep(NA,length(genALL_wHigher[,1]))
        for(j in 1:length(cladeNames_index1)){
            cladeNames_index1[which(genALL_wHigher[,"HIGHER_1"]==cladeNamesOrdered1[j])]<-j
        }

    #cladeNames2 <- names(table(genALL_wHigher$HIGHER_2))
    #cladeNamesOrdered2 <- cladeNames2[c(3,8,7,2,5,4)]
    #    cladeNames_index2 <- rep(NA,length(genALL_wHigher[,1]))
    #    for(j in 1:length(cladeNames_index2)){
    #        cladeNames_index2[which(genALL_wHigher[,"HIGHER_2"]==cladeNamesOrdered2[j])]<-j
    #    }

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
    ordNames_ext<-names(sort(table(as.vector(dat_ordsExtinct$order)), decreasing=TRUE)[c(1,4)])
        ordNames_index_ext <- rep(NA,length(genALL_wHigher[,1]))
        for(j in 1:length(ordNames_ext)){
            ordNames_index_ext[which(genALL_wHigher[,"order"]==ordNames_ext[j])]<-j
        }

    # join together
    genALL<-cbind.data.frame( cladeNames_index1, #cladeNames_index2, 
                                ordNames_index, ordNames_index_ext, genALL_wHigher)


# by CLADE NAME, pre_KPg, then post, then order
index<-order(genALL$cladeNames_index1, 
    #genALL$pre_KPg, #genALL$post_KPg, 
    genALL$order, genALL$intLength, #genALL$family, 
    decreasing=FALSE)

    # VISUALIZE-- plot data
    ###
    library(viridis); library(plotrix)

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

    pdf(file="plot_fossilGen_throughTime_byCLADE_order-family_wColors_7cats_byOrd-Int_agesFixed_wPiresEtAl.pdf", height=10, width=8)
    #quartz(height=50,width=6)
        plot(c(1:nrows)~seq(50,-120,length.out=nrows), col=NA, yaxt="n",xaxt="n",bty="n", xlab="",ylab="", xpd=NA)
        axis(side=1, at=seq(0,-120, by=-20), labels=seq(0,-120, by=-20))
        #add K-Pg
        abline(v= -66, lty=2, col="grey", lwd=3)
        #add PETM
        abline(v= -56, lty=2, col="grey", lwd=3)

        #plot(seq(50,-250,length.out=nrows)~c(1:nrows), col=NA, xaxt="n",yaxt="n",bty="n", xlab="",ylab="")
        #axis(side=2, at=seq(0,-200, by=-50), labels=seq(0,-200, by=-50))
        #abline(h= -66, lty=2, col="grey", lwd=4)

    #plotCI(x=c(1:nrows), y = -dat$midAll, ui= -dat$minAll, li= -dat$maxAll, err="y",
    plotCI(y=c(1:nrows), x = -dat$midAll, ui= -dat$minAll, li= -dat$maxAll, err="x",
           sfrac=0, gap=0, add=TRUE, col= lineCols, scol=lineCols, pch=20, lwd=1.5, cex=0.0000001)

    for(j in 1:length(cladeNamesOrdered1)){
        rowsOfClade <- c(1:nrows)[dat$cladeNames_index1==j]
    #   segments(x0=j, x1 = j, y0= min(rowsOfClade), y1 = max(rowsOfClade),
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







# Do SQS subsamplinging of the fossil record-- look at unbiased diversity through time
##########

SQStoRun<-seq(0.1,0.9,by=0.1)
for(q in 1:length(SQStoRun)){
    # do the SQS subsampling
    sqsGenusDiv_manual<-vector("list",length(binNames))
    for(i in 1:length(binNames)){
        bin<-binNames[i]
        t = table(genBins_ALL[which(genBins_ALL$binned==bin),"genus"])
        sqsGenusDiv_manual[[i]]<-sqs(ab=t[t != 0],q=SQStoRun[q],trials=1000, ignore.singletons=FALSE)
    }
    sqsGenusDiv_manual_ALL<-do.call(rbind,sqsGenusDiv_manual)
    rownames(sqsGenusDiv_manual_ALL)<-binNames

    sqsGenusDiv_manual_ALL[,1]<-bins[-length(bins)]
    colnames(sqsGenusDiv_manual_ALL)<-c("Time",colnames(sqsGenusDiv_manual_ALL)[-1])

  #  write.csv(sqsGenusDiv_manual_ALL, file="SQSresults_PBDB_quorum0p8_bins10Ma.csv")
  #  write.csv(sqsGenusDiv_manual_ALL, file="SQSresults_PBDB_quorum0p5_bins10Ma.csv")
  #  write.csv(sqsGenusDiv_manual_ALL, file=paste0("SQSresults_PBDB_quorum0p",q,"_binsStages.csv"))
    write.csv(sqsGenusDiv_manual_ALL, file=paste0("SQSresults_PBDB_quorum0p",q,"_bins5Ma_withPiresFrom1.csv"))
  #  write.csv(sqsGenusDiv_manual_ALL, file=paste0("SQSresults_PBDB_quorum0p",q,"_bins5Ma.csv"))

}  # end SQS runs


    # Do the plot of SQS subsampled diversity now.
    ######
    library(viridis)

    Lwds<-2
    q <- 5
    SQSdat <- read.csv(file=paste0("tables_DivTT_RateTT/SQSresults_PBDB_quorum0p",q,"_bins5Ma_withPiresFrom1.csv"))
    SQStime <- seq(1,126,by=5)+2.5
    COLS<-rep(c(gray(0.95), "white"), 13)

    pdf(file="fossilGeneraThroughTime_PBDB_SQS-0p1-3-5_5Ma-fullDat_manual_withPiresFrom1_wBins.pdf",height=5, width=7)

    par(oma = c(2,2,2,2)+0.1, mar = rep(2,4) + 0.1)

    plot(SQSdat[,"raw.richness"] ~ SQStime, type="n", xlim=c(115,0), log="y", ylim=c(1,2000), yaxt="n",bty="n",ylab="", xlab="")
    axis(side=4, at=c(5,20,200,2000))
    mtext(side=4,text="Fossil genera (log)",line=2.5, adj=0.85)
    mtext(side=1,text="Million years before present (Ma)",line=2.5, adj=0.5)

    for(k in 1: length(SQSdat[,"Time"])){
        rect(xleft=SQSdat[k,"Time"]+5, xright=SQSdat[k,"Time"], ybottom=1, ytop=3000, log="y", lty=1, lwd=1, col=COLS[k], border=NA)
    }

    # draw K-Pg
    abline(v=66, lty=1, lwd=3, col="white")
    abline(v=66, lty=2, lwd=3, col=magma(4)[2])
    # draw PETM
    abline(v=56, lty=1, lwd=3, col="white")
    abline(v=56, lty=2, lwd=3, col=magma(4)[3])
    points(SQSdat[,"raw.richness"] ~ SQStime, type="l", lwd=Lwds)#, xlim=c(110,0), log="y", ylim=c(5,2000))
    
    resSQS<-list()
    #for(q in 1:length(SQStoRun)){
    for(q in c(1,3,5)){
        resSQS[[q]]<-read.csv(file=paste0("tables_DivTT_RateTT/SQSresults_PBDB_quorum0p",q,"_bins5Ma_withPiresFrom1.csv"))
    
        points(resSQS[[q]][,"subsampled.richness"] ~ SQStime, type="l", lty=2, lwd=Lwds, col="black")#, xlim=c(110,0), log="y", ylim=c(5,200))
        text(x=5, y=resSQS[[q]][1,"subsampled.richness"], label=paste0("q=0.",q), adj=c(1,1))
    }
        text(x=56-6, y= 1, label="PETM", col=magma(4)[3])
        text(x=66+6, y= 1, label="K-Pg", col=magma(4)[2])

#    par(new=TRUE)
#    plot(rateDat[-c(1:exclVal),"origination"] ~ rateDat[-c(1:exclVal),"max_ma"], type="l", lwd=Lwds, col="blue", xlim=c(110,0), ylim=c(-0,1.5), bty="n",yaxt="n",xaxt="n",ylab="",xlab="")
#    points(rateDat[-c(1:exclVal),"extinction"] ~ rateDat[-c(1:exclVal),"max_ma"], type="l", lwd=Lwds, col="red")#, xlim=c(110,0), log="y", ylim=c(5,200))
#    axis(side=4, at=c(0,0.2,0.4))
#    mtext(side=4,text="Rate (genera / stage)",line=2.5, adj=-0.1)

   # legend(x=110.5, y=1.5, legend=c("total diversity","subsampled diversity","origination rate","extinction rate"),cex=0.8,lty=c(1,2,1,1),col=c(1,1,"blue","red"),lwd=rep(Lwds,4))

    dev.off()



##########
# Calculate the RATES now
##########
# using divDyn package
####
 #   divDyn(x, tax, bin = NULL, age = NULL, revtime = FALSE,
 #      breaks = NULL, coll = NULL, ref = NULL, om = NULL,
 #      noNAStart = FALSE, data.frame = TRUE, filterNA = FALSE)
library(viridis)

# Make the binNum into a factor
    genBins_ALL$binNum<-factor(genBins_ALL$binNum, levels=1:26)

# calculate metrics of diversity dynamics
    mamDivDyn <- divDyn(genBins_ALL, tax="genus", bin="binNum")
    write.csv(mamDivDyn, file="mamDivDyn_ON_pbdb_Mammalia_wPiresEtAl2018_binned1to131.csv")

# load back
    mamDivDyn<-read.csv(file="mamDivDyn_ON_pbdb_Mammalia_wPiresEtAl2018_binned1to131.csv", header=TRUE)


    # PLOTTING
    ######    
    # Get the 5 Ma bin details...
        bins_5Ma<- seq(from=1,to=131,by= 5)
        binNames_5Ma<-paste0("bin_",bins_5Ma,"_",bins_5Ma[-1])[-27]

        binDeets<-cbind.data.frame(bin=binNames_5Ma, 
            bottom=bins_5Ma[-1], top=bins_5Ma[-27], mid=bins_5Ma[-27]+2.5,
            period=c(rep("N",4),rep("Pg",9), rep("K",13)),
            epoch=c("Po",rep("Mi",3), rep("Og",3),rep("Eo",4), rep("Pe",2), rep("Late",7), rep("Early",6))
            )
            #stages[which(stages$bottom < 131),]

    # test plot of richness    
        tsplot(binDeets, shading="bin", boxes=c("period","epoch"), labels.args= list( cex=c(0.8)),
            xlim=c(131,0), ylab="range-through diversity (genera)", ylim=c(0,2000))
           lines(binDeets$mid, mamDivDyn$divRT, lwd=2)

    # just get the extinctionRate metrics:
        extinctionRateNames<-c('extPC', 'ext3t', 'extC3t', 'extGF', 'E2f3', 'ext2f3')
        mamExtRates<-mamDivDyn[,c("binNum",extinctionRateNames)]
    # just get the originationRate metrics:
        originationRateNames<-c('oriPC', 'ori3t', 'oriC3t', 'oriGF', 'O2f3', 'ori2f3')
        mamOriRates<-mamDivDyn[,c("binNum",originationRateNames)]

        # for each bin, get the high and low extinctionRate estimate
        extRes<-list(); oriRes<-list()
        for(k in 2:(length(mamExtRates[,1])-2) ){
            
            highExt<-max( na.omit(as.numeric(mamExtRates[k,2:7])) )
            lowExt<-min( na.omit(as.numeric(mamExtRates[k,2:7])) )
            extRes[[k]]<-cbind.data.frame(highExt, lowExt)

            highOri<-max( na.omit(as.numeric(mamOriRates[k,2:7])) )
            lowOri<-min( na.omit(as.numeric(mamOriRates[k,2:7])) )
            oriRes[[k]]<-cbind.data.frame(highOri, lowOri)
            
        }
        extRes_ALL<-rbind.data.frame(rep(NA,4), do.call(rbind,extRes),rep(NA,4),rep(NA,4))
        oriRes_ALL<-rbind.data.frame(rep(NA,4), do.call(rbind,oriRes),rep(NA,4),rep(NA,4))

    # join with the full rateMatrix and binDeets
        mamDivDyn_wbins_wHiLow<-cbind.data.frame(binDeets, mamDivDyn, extRes_ALL, oriRes_ALL )

        write.csv(mamDivDyn_wbins_wHiLow, file="mamDivDyn_ON_pbdb_Mammalia_wPiresEtAl2018_binned1to131_wHiLow.csv")

    # reload rates now
    mamDivDyn_wbins_wHiLow<-read.csv(file="mamDivDyn_ON_pbdb_Mammalia_wPiresEtAl2018_binned1to131_wHiLow.csv", header=TRUE)


    # FINAL PLOT of rates as high & low polygon
        rowsWithoutNA<-2:24
        oriIntCol<-gray(0.8,alpha=0.6)
        extIntCol<-gray(0.1,alpha=0.6)
        lwdWidths<-0.5
        ltyTypes<-1
        rateCols<-"white"
        LOG<-0.013

       #pdf(file="mamDivDyn_extinctionOrigination_poly_shaded_logY.pdf", height=5, width=7)
       pdf(file="mamDivDyn_extinctionOrigination_poly_shaded.pdf", height=5, width=7)
       #pdf(file="mamDivDyn_extinctionOrigination_poly_shaded-wide.pdf", width=10, height=5)
            tsplot(mamDivDyn_wbins_wHiLow, #plot.args=list(log="y"), 
                shading="bin", shading.col=c("white",gray(0.95)), boxes=c("period","epoch"), labels.args = list(cex=0.8),
                xlim=c(131,0), ylim=c(0.01,3), ylab="Fossil rates")
            mtext(side=3, text="Extinction and Origination")

                # draw K-Pg
                abline(v=66, lty=2, lwd=3, col=gray(0.5))
                    text(x=66+4.5, y= LOG, label="K-Pg", col=gray(0.5), font=2, cex=0.8)

                # draw PETM
                abline(v=56, lty=2, lwd=3, col=gray(0.5))
                    text(x=56-4.5, y= LOG, label="PETM", col=gray(0.5), font=2, cex=0.8)

                #extinction
                polygon(x=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ),
                        y=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"highExt"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"lowExt"])
                            )+LOG, 
                        lwd=2, border="black",col=extIntCol)#rateCols[q])

                    ##All 6 rates as thin line
                    #for(q in 1:length(extinctionRateNames)){
                    #    lines(mamDivDyn_wbins_wHiLow$mid, mamDivDyn_wbins_wHiLow[,extinctionRateNames[q]], lwd=lwdWidths, lty=ltyTypes, col=rateCols)#"black")#
                    #}
                    

                #origination
                polygon(x=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ),
                        y=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"highOri"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"lowOri"])
                            )+LOG, 
                        lwd=2, border="black",col=oriIntCol)#rateCols[q])

            legend("topleft", legend=c("origination rate", "extinction rate"), 
              fill=c(oriIntCol,extIntCol), border="black", bg="white")
        dev.off()



    # PLOT of individual rates
        rateCols<-viridis(6)
        lwdWidths<-seq(from=8, to=2, length.out=6)
        ltyTypes<-c(1, 1, 1, 6, 6, 6)

        # extinction--all 6 rates
        pdf(file="mamDivDyn_extinction_6rates.pdf", width=6, height=5)
            tsplot(mamDivDyn_wbins_wHiLow, shading="bin", boxes="cat", xlim=c(131,0), 
                 ylab="Fossil rates", ylim=c(0,3))
            abline(v=66, lty=3, lwd=2)
            mtext(side=3, text="Extinction")
            for(q in 1:length(extinctionRateNames)){
                lines(mamDivDyn_wbins_wHiLow$mid, mamDivDyn_wbins_wHiLow[,extinctionRateNames[q]], lwd=lwdWidths[q], lty=ltyTypes[q], col=rateCols[q])#"black")#
            }
            legend("topleft", legend=extinctionRateNames, 
              col=rateCols, lwd=lwdWidths, lty=ltyTypes, bg="white")
        dev.off()
        # origination--all 6 rates
        pdf(file="mamDivDyn_origination_6rates.pdf", width=6, height=5)
            tsplot(mamDivDyn_wbins_wHiLow, shading="bin", boxes="cat", xlim=c(131,0), 
                 ylab="Fossil rates", ylim=c(0,3))
            abline(v=66, lty=3, lwd=2)
            mtext(side=3, text="Origination")
            for(q in 1:length(originationRateNames)){
                lines(mamDivDyn_wbins_wHiLow$mid, mamDivDyn_wbins_wHiLow[,originationRateNames[q]], lwd=lwdWidths[q], lty=ltyTypes[q], col=rateCols[q])#"black")#
            }
            legend("topleft", legend=originationRateNames, 
              col=rateCols, lwd=lwdWidths, lty=ltyTypes, bg="white")
        dev.off()


    # ADD POLYGON - Min & Max of rates
        rowsWithoutNA<-2:24
        intCol<-gray(0.2,alpha=0.5)
        #extinction
            pdf(file="mamDivDyn_extinction_minMax_poly_all6.pdf", width=6, height=5)
            #pdf(file="mamDivDyn_extinction_minMax_poly.pdf", width=6, height=5)

                tsplot(binDeets, shading="bin", boxes="cat", xlim=c(131,0), 
                     ylab="Fossil rates", ylim=c(0,3))
                abline(v=66, lty=3, lwd=2)
                mtext(side=3, text="Extinction")
        
                polygon(x=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ),
                        y=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"highExt"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"lowExt"])
                            ), 
                        lwd=1, border=intCol,col=intCol)#rateCols[q])

                for(q in 1:length(extinctionRateNames)){
                    lines(mamDivDyn_wbins_wHiLow$mid, mamDivDyn_wbins_wHiLow[,extinctionRateNames[q]], lwd=lwdWidths[q], lty=ltyTypes[q], col=rateCols[q])#"black")#
                }
                legend("topleft", legend=extinctionRateNames, 
                  col=rateCols, lwd=lwdWidths, lty=ltyTypes, bg="white")
            dev.off()

        #originaton
            pdf(file="mamDivDyn_originaton_minMax_poly_all6.pdf", width=6, height=5)

                tsplot(mamDivDyn_wbins_wHiLow, shading="bin", boxes="cat", xlim=c(131,0), 
                     ylab="Fossil rates", ylim=c(0,3))
                abline(v=66, lty=3, lwd=2)
                mtext(side=3, text="Originaton")
        
                polygon(x=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ),
                        y=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"highOri"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"lowOri"])
                            ), 
                        lwd=1, border=intCol,col=intCol)#rateCols[q])

                for(q in 1:length(originationRateNames)){
                    lines(mamDivDyn_wbins_wHiLow$mid, mamDivDyn_wbins_wHiLow[,originationRateNames[q]], lwd=lwdWidths[q], lty=ltyTypes[q], col=rateCols[q])#"black")#
                }
                legend("topleft", legend=originationRateNames, 
                  col=rateCols, lwd=lwdWidths, lty=ltyTypes, bg="white")
            dev.off()

    # COMBINED extinction & origination
        extRateCols<-viridis(6)
        oriRateCols<-magma(6)

        lwdWidths<-seq(from=8, to=2, length.out=6)
        ltyTypes<-c(1, 1, 1, 6, 6, 6)

        # both--all 6 rates
        pdf(file="mamDivDyn_extinctionOrigination_6rates.pdf", width=6, height=5)
            tsplot(mamDivDyn_wbins_wHiLow, shading="bin", boxes="cat", xlim=c(131,0), 
                 ylab="Fossil rates", ylim=c(0,3))
            abline(v=66, lty=3, lwd=2)
            mtext(side=3, text="Extinction and Origination")
            for(q in 1:length(extinctionRateNames)){
                lines(mamDivDyn_wbins_wHiLow$mid, mamDivDyn_wbins_wHiLow[,extinctionRateNames[q]], lwd=lwdWidths[q], lty=ltyTypes[q], col=extRateCols[q])#"black")#
            }        
            for(q in 1:length(originationRateNames)){
                lines(mamDivDyn_wbins_wHiLow$mid, mamDivDyn_wbins_wHiLow[,originationRateNames[q]], lwd=lwdWidths[q], lty=ltyTypes[q], col=oriRateCols[q])#"black")#
            }

            legend("topright", legend=originationRateNames, 
              col=oriRateCols, lwd=lwdWidths, lty=ltyTypes, bg="white")

            legend("topleft", legend=extinctionRateNames, 
              col=extRateCols, lwd=lwdWidths, lty=ltyTypes, bg="white")
        dev.off()


        # both--as POLYGONS -- black & white
        rowsWithoutNA<-2:24
    #    extIntCol<-rgb(1,0,0,alpha=0.5)
    #    oriIntCol<-rgb(0,0,1,alpha=0.5)
        extIntCol<-gray(0.5,alpha=0.5)
        oriIntCol<-gray(0.99,alpha=0.8)

        #pdf(file="mamDivDyn_extinctionOrigination_poly.pdf", width=6, height=5)
        pdf(file="mamDivDyn_extinctionOrigination_poly_bw.pdf", width=6, height=5)
            tsplot(mamDivDyn_wbins_wHiLow, shading="bin", shading.col=c(gray(0.95),"white"), boxes="cat", xlim=c(131,0), 
                 ylab="Fossil rates", ylim=c(0,3))
            abline(v=66, lty=3, lwd=2)
            mtext(side=3, text="Extinction and Origination")

                #extinction
                polygon(x=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ),
                        y=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"highExt"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"lowExt"])
                            ), 
                        lwd=3, border="black",col=extIntCol)#rateCols[q])

                #origination
                polygon(x=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ),
                        y=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"highOri"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"lowOri"])
                            ), 
                        lwd=3, border=gray(0.5),col=oriIntCol)#rateCols[q])

            legend("topleft", legend=c("origination rate", "extinction rate"), 
              col=c(gray(0.5),"black"), lwd=10, lty=1, bg="white")
        dev.off()

        # both--as POLYGONS -- red & blue
        rowsWithoutNA<-2:24
        oriIntCol<-rgb(1,0,0,alpha=0.5)
        extIntCol<-rgb(0,0,1,alpha=0.5)

        pdf(file="mamDivDyn_extinctionOrigination_poly_redBlue-wide.pdf", width=10, height=5)
            tsplot(mamDivDyn_wbins_wHiLow, shading="bin", shading.col=c(gray(0.95),"white"), boxes="cat", xlim=c(131,0), 
                 ylab="Fossil rates", ylim=c(0,3))
            abline(v=66, lty=3, lwd=2)
            mtext(side=3, text="Extinction and Origination")

                #extinction
                polygon(x=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ),
                        y=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"highExt"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"lowExt"])
                            ), 
                        lwd=0.1, border="blue",col=extIntCol)#rateCols[q])

                #origination
                polygon(x=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"mid"])
                            ),
                        y=c((mamDivDyn_wbins_wHiLow[rowsWithoutNA,"highOri"])
                            ,rev(mamDivDyn_wbins_wHiLow[rowsWithoutNA,"lowOri"])
                            ), 
                        lwd=0.1, border="red",col=oriIntCol)#rateCols[q])

          #  legend("topleft", legend=c("origination rate", "extinction rate"), 
          #    col=c(oriIntCol,extIntCol), lwd=10, lty=1, bg="white")
        dev.off()



