# Code to summarized the calculated pulled and pushed speciation rates
#########
# For Fig 4B of new conception Current Biology
#######
# NS Upham; last modification: 6 April 2021
# ===================================================

library(castor); library(ape); library(phytools); library(phangorn)

# set wd
####
dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses"
setwd(dirname)

	# # load in 100 trees
	# ######
	# mamPhy100_wOut<-read.nexus("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL/MammaliaTrees_SepOct2019/node-dated_backbone_NDexp_10ktrees_tipDR/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_nexus.trees")
	# ntrees<-length(mamPhy100_wOut)
	# 
	# mamPhy100<-list()
	# for(i in 1:ntrees){
	# 	mamPhy100[[i]]<-drop.tip(mamPhy100_wOut[[i]],"_Anolis_carolinensis")
	# }
	# class(mamPhy100)<-"multiPhylo"

# Load in the *breakpoints* from the fossil ages.
######
age_grid <- seq(1, 111, by=5)

# re-load in the PSR results, from the HPC
######
whichTrees<-c(1:100)#[-c(67)]

res_PSR_ALL<-vector("list",100)
for(i in whichTrees){

	res_PSR_ALL[[i]] <- get(load(file=paste0("pulledSpeciation_100trees_new/PSR_5911species_NDexp_tree",i,"_new.Rda")))

}



# Load in the fossil diversity
#####
  fossilDir<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-geo-trait-analyses/PBDB_FOSSIL_DATA"
  
  # Function to calculate the REAL LAMBDA
  getLambda <- function(PSR, E){
    result <- PSR / (1-E)
    return(result)
  }

  # # 5 Ma time bins
  # #########
  # perCapRates_5Ma <- read.csv(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-geo-trait-analyses/PBDB_FOSSIL_DATA/tables_DivTT_RateTT/fossilGeneraThroughTime_PBDB_perCapitaRates_bins5Ma.csv", header=TRUE)
  #   # 'impute' the missing rate data from 70-65 Ma, taking the mean of the before / after values
  #   perCapRates_5Ma[13,"pRate"]<-mean(c(perCapRates_5Ma[12,"pRate"],perCapRates_5Ma[14,"pRate"]))
  #   perCapRates_5Ma[13,"qRate"]<-mean(c(perCapRates_5Ma[12,"qRate"],perCapRates_5Ma[14,"qRate"]))

  # Calculate the CORRECTED lambda (and pulled speciation rates), across multiple SUBSAMPLES
  # ======================
  # do FOR EACH OF 100 TREES
  LAMBDA_bySQS_ALL100<-vector("list",100)
  allPSR<-vector("list",100)
  allLTT<-vector("list",100)
  for(i in whichTrees){

    allPSR[[i]]<-res_PSR_ALL[[i]]$fitted_PSR
    allLTT[[i]]<-res_PSR_ALL[[i]]$fitted_LTT

    LTT_i<-res_PSR_ALL[[i]]$fitted_LTT
    names(LTT_i)<-age_grid

    PSR_i<-cbind.data.frame(time=age_grid, est=res_PSR_ALL[[i]]$fitted_PSR, low95=res_PSR_ALL[[i]]$CI95lower, high95=res_PSR_ALL[[i]]$CI95upper)
    rownames(PSR_i)<-age_grid
    
    #LAMBDA_bySQS_q<-list()
    #for(q in 1:length(SQStoRun)){
    #for(q in c(1,3,5,7)){
    q=3
        # get SQS fossil data:
        resSQS<-read.csv(file=paste0(fossilDir,"/tables_DivTT_RateTT/SQSresults_PBDB_quorum0p",q,"_bins5Ma_withPiresFrom1.csv"))

      # # IF "E" is based on EXTINCTION RATE
      # ########
      # extinctionRate <- perCapRates_5Ma[-c(1,26),"qRate"]
      # totalExtantLineages <- rev(allLTT_ALL_Med)[7:30]
      # totalMissingLineages <- totalExtantLineages * extinctionRate
      # ALL_lineages <- totalExtantLineages + totalMissingLineages
      # 
      # fractionMissingLineages <- totalMissingLineages / ALL_lineages
      # 
      # LAMBDA <- getLambda(PSR = allPSR_ALL_Med[-c(1, 26:31)], E = rev(fractionMissingLineages))
      # timesToPlot<-rev(perCapRates_5Ma[-c(1,26),"X.2"])

        # IF "E" is based on LTT lineages and subsampled fossils...
        ########
        # per 100 trees
        totalExtantLineages_i <- LTT_i[1:23]    
          # if shift by 5 Ma
         # totalMissingLineages_i <- c(0,resSQS[c(1:22),"subsampled.richness"])
          # if no shift (take as given)
          totalMissingLineages_i <- resSQS[c(1:23),"subsampled.richness"]

        ALL_lineages_i <- totalExtantLineages_i + totalMissingLineages_i
        fractionMissingLineages_i <- totalMissingLineages_i / ALL_lineages_i
        
        LAMBDA_bySQS_est    <- getLambda(PSR = (PSR_i[1:23,"est"]), E = fractionMissingLineages_i)#, q=as.numeric(paste0("0.",q)) )
        LAMBDA_bySQS_low95  <- getLambda(PSR = (PSR_i[1:23,"low95"]), E = fractionMissingLineages_i)#, q=as.numeric(paste0("0.",q)) )
        LAMBDA_bySQS_high95 <- getLambda(PSR = (PSR_i[1:23,"high95"]), E = fractionMissingLineages_i)#, q=as.numeric(paste0("0.",q)) )
        time <- PSR_i[1:23,"time"]
  #}
    LAMBDA_bySQS_ALL100[[i]]<-cbind.data.frame(time, LAMBDA_bySQS_est, LAMBDA_bySQS_low95, LAMBDA_bySQS_high95) 
  }

allPSR_ALL100<-do.call(rbind, allPSR)[,1:23]
allLTT_ALL100<-do.call(rbind, allLTT)[,1:23]

# just get the EST from all trees
##
toCombine<-vector("list",100)
for(i in whichTrees){
  res4<-LAMBDA_bySQS_ALL100[[i]]
  toCombine[[i]]<-res4[,2]
}
LAMBDA_bySQS_ALL100_EST<-do.call(rbind,toCombine)

  # get the 95% CIs and medians of the EST only
  pushedRates_quant_j<-vector("list", length(age_grid))
  pulledRates_quant_j<-vector("list", length(age_grid))
  LTT_quant_j<-vector("list", length(age_grid))
  for(j in 1:length(age_grid)){
    res1_j<-LAMBDA_bySQS_ALL100_EST[,j]
    pushedRates_quant_j[[j]]<-quantile(res1_j, c(0.025,0.975,0.5))

    res2_j<-allPSR_ALL100[,j]
    pulledRates_quant_j[[j]]<-quantile(res2_j, c(0.025,0.975,0.5))

    res3_j<-allLTT_ALL100[,j]
    LTT_quant_j[[j]]<-quantile(res3_j, c(0.025,0.975,0.5))

  }

  # NOW LOAD IN THE HBD (homogeneous birth-death) ESTIMATES WITH FIXED EXTINCTION
  #########
  # I have 76 trees only from this, for some reason. (others didn't run.... maybe something with breakpoints in the phylogeny at the 5Ma intervals)
  #whichTrees<-c(1, 2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 17, 19, 21, 22, 23, 25, 26, 27, 28, 29, 30, 34, 36, 37, 38, 39, 40, 42, 43, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 
  #            63, 64, 66, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 87, 90, 92, 93, 94, 95, 96, 98, 100)
  ######
  # I now have 100 trees from this, across all but tree90 and 91 for 6 types of extinction rate...
  whichTrees<-c(1:100)

extinctionRateNames<-c('extPC', 'ext3t', 'extC3t', 'extGF', 'E2f3', 'ext2f3')
estHBD_Lam_quant_ALL<-list()

for(k in 1:length(extinctionRateNames)){

  estHBD_all<-vector("list",100)
  for(i in whichTrees){
    estHBD_all[[i]]<-get(load(paste0("estHBD_fixedMu_100trees_all6/estHBD_fixedMu_5911species_NDexp_tree",i,"_1trial_all6_",extinctionRateNames[k],".Rda")))
  }

  # get all the estimated LAMBDAs
  allLam<-vector("list",100)
  for(i in whichTrees){
    allLam[[i]]<-estHBD_all[[i]]$fitted_lambda
  }
  allLam_ALL<-do.call(rbind,allLam)

  # get the 95% CI and MED of those TRUE LAMBDAs
  estHBD_Lam_quant_j<-list()
  for(j in 1:length(age_grid)){
    resLam<-allLam_ALL[,j]
    estHBD_Lam_quant_j[[j]]<-quantile(resLam, c(0.025, 0.975, 0.5))
  }
  estHBD_Lam_quant_ALL[[k]]<-do.call(rbind,estHBD_Lam_quant_j)
  colnames(estHBD_Lam_quant_ALL[[k]])<-c( paste0(extinctionRateNames[k],"_", c("Lam_low95","Lam_high95","Lam_med")) )
  
}
estHBD_Lam_quant_ALL_6rates<-do.call(cbind, estHBD_Lam_quant_ALL)


  # Now get the 95% CI and med of ALL SIX fixMu Lam rates simultaneously
  #####
  allLam_all6<-list()
  for(k in 1:length(extinctionRateNames)){

    allLam<-vector("list",100)
    for(i in whichTrees){
      res<-get(load(paste0("estHBD_fixedMu_100trees_all6/estHBD_fixedMu_5911species_NDexp_tree",i,"_1trial_all6_",extinctionRateNames[k],".Rda")))
      allLam[[i]]<-res$fitted_lambda
    }
    allLam_all6[[k]]<-do.call(rbind,allLam)
  }
  allLam_all6_ALL<-do.call(rbind, allLam_all6)

      # get the 95% CI and MED of those 6rate estimates of TRUE LAMBDAs
      allLam_all6_ALL_95perGrid<-list()
      for(j in 1:length(age_grid)){
        resLam<-allLam_all6_ALL[,j]
        allLam_all6_ALL_95perGrid[[j]]<-quantile(resLam, c(0.025, 0.975, 0.5))
      }
      allLam_all6_ALL_95perGrid_ALL<-do.call(rbind,allLam_all6_ALL_95perGrid)
      colnames(allLam_all6_ALL_95perGrid_ALL)<-paste0("all6Mu_", c("Lam_low95","Lam_high95","Lam_med")) 
      

# COMBINE this into a cool SUMMARY TABLE
###
pushedRates_quant_ALL<-cbind.data.frame(time=age_grid, 
                                        mid=age_grid+2.5, 
                                        do.call(rbind,LTT_quant_j),
                                        do.call(rbind,pulledRates_quant_j), 
                                        do.call(rbind,pushedRates_quant_j), 
                                     #   rbind(rep(0,dim(resSQS)[2]),resSQS[c(1:22),])  
                                     #  resSQS[c(1:23),],
                                       estHBD_Lam_quant_ALL_6rates,
                                       allLam_all6_ALL_95perGrid_ALL
                                        )

colnames(pushedRates_quant_ALL)<-c("time", "mid",
                                    paste0("LTT_",c("low95","high95","med")),
                                    paste0("pulled_",c("low95","high95","med")), 
                                    paste0("E(t)_Lam_",c("low95","high95","med")), 
                                    colnames(estHBD_Lam_quant_ALL_6rates),
                                    colnames(allLam_all6_ALL_95perGrid_ALL)
                                    )

write.csv(pushedRates_quant_ALL, file="allRates_withFossilData_1to111Ma_bins5Ma_wPires.csv")

# LOAD BACK IN THE SUMMARY TABLE
pushedRates_quant_ALL<-read.csv(file="allRates_withFossilData_1to111Ma_bins5Ma_wPires.csv")

#       # 5 Ma time bins
#       #########
#       perCapRates_5Ma <- read.csv(file=paste0(fossilDir,"/tables_DivTT_RateTT/fossilGeneraThroughTime_PBDB_perCapitaRates_bins5Ma.csv"), header=TRUE)
#         # 'impute' the missing rate data from 70-65 Ma, taking the mean of the before / after values
#       #  perCapRates_5Ma[13,"pRate"]<-mean(c(perCapRates_5Ma[12,"pRate"],perCapRates_5Ma[14,"pRate"]))
#       #  perCapRates_5Ma[13,"qRate"]<-mean(c(perCapRates_5Ma[12,"qRate"],perCapRates_5Ma[14,"qRate"]))
#
#    allRates_quant_ALL<- cbind.data.frame(pushedRates_quant_ALL, perCapRates_5Ma[rev(4:26),])



# PLOT ALL -- simplified-- Only the PULLED and PUSHED (fixedExtinction) -- 95% CIs and MEDIANS of each
##########
library(viridis); library(RColorBrewer)

#lamCols<-viridis(8)
#lamCols<-YlOrRd(9)[c(9,7,5)]
lamCols<-magma(10, alpha=0.5)
#lamCols<-rep("darkorchid1",10)

pdf(file="fittedRates_PSR_wLambda_EstHBDfixedMu_111to1Ma_95CI-all_6rates_medAll_logY_max.pdf", height=7, width=7)
#pdf(file="fittedRates_PSR_wLambda_EstHBDfixedMu_111to1Ma_95CI-all_6rates_medAll.pdf", height=7, width=7)
#pdf(file="fittedRates_PSR_wLambda_EstHBDfixedMu_111to1Ma_95CI-summary6rates.pdf", height=7, width=7)

YMAX<- 65#2.2
BASE<- 0#0.8*2
XMIN<- 115
COLS<-rep(c(gray(0.95), "white"), 13)
LOG<-0.001

    # draw the BASE
    ##############
    plot(    x     = pushedRates_quant_ALL$mid,
             y     = pushedRates_quant_ALL$pulled_high95,
             main  = '',
             xlab  = '',
             ylab  = '',
             type  = 'n',
             xlim  = c(XMIN,0),
             ylim  = c(0.01,YMAX),
             bty="n",
             log="y",
             xaxt="n",
             #yaxt="n",
             xpd=FALSE) 

    for(q in 1: length(pushedRates_quant_ALL[,"time"])){
        rect(xleft=pushedRates_quant_ALL[q,"time"]+5, xright=pushedRates_quant_ALL[q,"time"], ybottom=LOG, ytop=65, lty=1, lwd=1, col=COLS[q], border=NA)
    }
    
    #axis(side=2, at=c(0, 0.2, 0.4))
    axis(side=1)

    # draw the 95% CI of 100-tree PUSHED Lam results 
    ##############
      #  SUBSET<-13 # 60 Ma
      #  SUBSET<-16 # 75 Ma
      #  SUBSET<-15 # 70 Ma
         SUBSET<-23 # 110 Ma


      # all 6 rates -- fixedFossilExtinction
      #######
       extinctionRateNames<-c('extPC', 'ext3t', 'extC3t', 'extGF', 'E2f3', 'ext2f3')

      for(k in 1:length(extinctionRateNames)){
      # get target data
      pushedRates_fixedMu<-pushedRates_quant_ALL[,c("mid",paste0(extinctionRateNames[k],"_", c("Lam_low95","Lam_high95","Lam_med"))) ]

        # 95% CI of 100 trees
        toPlot<-pushedRates_fixedMu[1:SUBSET,]
        polygon(x = c(toPlot$mid, rev(toPlot$mid)), 
              y= c(toPlot[,2], rev(toPlot[,3]))+LOG, 
            col=rgb(1,0.3,0.3,alpha=0.5),#magma(10,alpha=0.5)[6],#, 
            border=NA)
        # median of 100 trees
        #points(x= toPlot$mid, y=toPlot[,4], type="l", lwd=3, col="black")#rgb(0,0,1))
      }
        points(x= pushedRates_quant_ALL$mid, y=pushedRates_quant_ALL[,"all6Mu_Lam_med"]+LOG, type="l", lwd=3, col="white")#rgb(0,0,1))


    # draw the 95% CI of 100-tree PULLED Lam results 
    ##############
      # get target data
      pulledRates<-pushedRates_quant_ALL[,c(2,6:8)]

        # 95% CI of 100 trees
        toPlot<-pulledRates[1:SUBSET,]
        polygon(x = c(toPlot$mid, rev(toPlot$mid)), 
              y= c(toPlot$pulled_low95, rev(toPlot$pulled_high95))+LOG, 
            col=magma(10,alpha=0.7)[4],#rgb(1,0,0,alpha=0.3), 
            border=NA)
        # median of 100 trees
        points(x= toPlot$mid, y=toPlot$pulled_med+LOG, type="l", lwd=3, col="white")#rgb(1,0,0))

  # general labels
  ####
    # draw K-Pg
        abline(v=66, lty=2, lwd=3, col=gray(0.5))
        text(x=66+6, y= 0.01, label="K-Pg", col=gray(0.5), font=2)
    # draw PETM
        abline(v=56, lty=2, lwd=3, col=gray(0.5))
        text(x=56-6, y= 0.01, label="PETM", col=gray(0.5), font=2)

dev.off()
      

pdf(file="colors_magma10.pdf")
  plot(1:10,pch=20,col=lamCols,cex=5)
dev.off()


# PLOT ALL -- simplified-- just the 95% CIs and MEDIANS of each
##########

pdf(file="fittedRates_PSR_wLambda_fromSQS0p3_andEstHBDfixedMu_111to1Ma_95CI_all3_new_wGrid.pdf", height=10, width=7)

BASE<- 0.8*2
XMIN<- 115
COLS<-rep(c(gray(0.95), "white"), 13)

    # draw the BASE
    ##############
    plot(    x     = pushedRates_quant_ALL$time,
             y     = pushedRates_quant_ALL$pulled_high95,
             main  = '',
             xlab  = '',
             ylab  = '',
             type  = 'n',
             xlim  = c(XMIN,0),
             ylim  = c(-BASE,0.5),
             #log="y",
             bty="n",
             xaxt="n",
             yaxt="n",
             xpd=FALSE) 
    
        for(q in 1: length(pushedRates_quant_ALL[,"time"])){
        rect(xleft=pushedRates_quant_ALL[q,"time"]+5, xright=pushedRates_quant_ALL[q,"time"], ybottom=-2, ytop=1, lty=1, lwd=1, col=COLS[q], border=NA)
        }

    axis(side=2, at=c(0,0.2, 0.4))
    axis(side=1)
    # draw K-Pg
    abline(v=66, lwd=2, lty=2)
    # draw PETM
    abline(v=56, lwd=2, lty=2)

    # draw the 95% CI of 100-tree PULLED Lam results 
    ##############
      #  SUBSET<-13 # 60 Ma
      #  SUBSET<-16 # 75 Ma
      #  SUBSET<-15 # 70 Ma
         SUBSET<-23 # 110 Ma
      # get target data
      pulledRates<-pushedRates_quant_ALL[,c(3,7:9)]

        # 95% CI of 100 trees
        toPlot<-pulledRates[1:SUBSET,]
        polygon(x = c(toPlot$mid, rev(toPlot$mid)), 
              y= c(toPlot$pulled_low95, rev(toPlot$pulled_high95)), 
            col=rgb(1,0,0,alpha=0.3), border=NA)
        # median of 100 trees
        points(x= toPlot$mid, y=toPlot$pulled_med, type="l", lwd=3, col=rgb(1,0,0))


    # draw the 95% CI of 100-tree PUSHED Lam results 
    ##############
      # from (fixedFossilLineages)
      #####
      # get target data
      pushedRates_fixedLineages<-pushedRates_quant_ALL[,c("mid","E.t._Lam_low95", "E.t._Lam_high95", "E.t._Lam_med")]

        # 95% CI of 100 trees
        DOWN<- BASE/2

        toPlot<-pushedRates_fixedLineages[1:SUBSET,]
        polygon(x = c(toPlot$mid, rev(toPlot$mid)), 
              y= c(toPlot$E.t._Lam_low95, rev(toPlot$E.t._Lam_high95))-DOWN, 
            col=rgb(0,1,0,alpha=0.3), border=NA)
        # median of 100 trees
        points(x= toPlot$mid, y=toPlot$E.t._Lam_med-DOWN, type="l", lwd=3, col=rgb(0,1,0))

      axis(side=2, at=c(-DOWN, -DOWN+0.2, -DOWN+0.4), labels=c(0, 0.2, 0.4))

      # from (fixedFossilExtinction)
      #######
      # get target data
      pushedRates_fixedMu<-pushedRates_quant_ALL[,c("mid","all6Mu_Lam_low95","all6Mu_Lam_high95","all6Mu_Lam_med")]

        # 95% CI of 100 trees
        DOWN<- BASE

        toPlot<-pushedRates_fixedMu[1:SUBSET,]
        polygon(x = c(toPlot$mid, rev(toPlot$mid)), 
              y= c(toPlot$all6Mu_Lam_low95, rev(toPlot$all6Mu_Lam_high95))-DOWN, 
            col=rgb(0,0,1,alpha=0.3), border=NA)
        # median of 100 trees
        points(x= toPlot$mid, y=toPlot$all6Mu_Lam_med-DOWN, type="l", lwd=3, col=rgb(0,0,1))

      axis(side=2, at=c(-DOWN, -DOWN+0.2, -DOWN+0.4), labels=c(0, 0.2, 0.4))

dev.off()
      






# PLOT ALL as overlaying
######
library(miscTools); library(viridis); library(RColorBrewer)

#lamCols<-viridis(8)
#lamCols<-YlOrRd(9)[c(9,7,5)]
lamCols<-magma(10)
#lamCols<-rep("darkorchid1",10)
allPSR<-list(); allLTT<-list()
#png(file="fittedRates_PSR_99trees_w95CI_wLambda_fromSQS-n_multiple_justGrey2.png", height=5, width=7, res=300, units="in")

pdf(file="fittedRates_PSR_100trees_w95CI_wLambda_fromSQS0p3_andEstHBDfixedMu_111to1Ma_.pdf", height=10, width=7)
#pdf(file="fittedRates_PSR_99trees_w95CI_wLambda_fromSQS0p3_70to0Ma_100treeALL_noShift.pdf", height=5, width=7)

#png(file="fittedRates_PSR_99trees_w95CI_wLambda_fromSQS-n_multiple_justGrey2.png", height=5, width=7, res=300, units="in")
#pdf(file="fittedRates_PSR_99trees_w95CI_wLambda_fromSQS-n_multiple.pdf", height=5, width=7)
#pdf(file="fittedRates_PSR_99trees_w95CI_wLambda_fromExRate.pdf", height=5, width=5)
#pdf(file="fittedRates_PSR_w95CI.pdf", height=5, width=5)

BASE<- 0.8*2
XMIN<- 115
  # plot the PULLED RATES, all 100 trees with 95% CIs
    plot( x     = res_PSR_ALL[[1]]$age_grid,
             y     = res_PSR_ALL[[1]]$fitted_PSR,
             main  = '',
             xlab  = '',
             ylab  = '',
             type  = 'n',
             xlim  = c(XMIN,0),
             ylim  = c(-BASE,0.5),
             bty="n",
             xaxt="n",
             yaxt="n",
             xpd=FALSE) 
    axis(side=2, at=c(0,0.2, 0.4))
    axis(side=1)

    for(i in whichTrees){
    res<-res_PSR_ALL[[i]]
    polygon(x = c((res$age_grid[1:24]), rev(res$age_grid[1:24])), y= c(res$CI95lower[1:24], rev(res$CI95upper[1:24])), 
 #  polygon(x = c((res$age_grid), rev(res$age_grid)), y= c(res$CI95lower, rev(res$CI95upper)), 
    		col=grey(0.5,alpha=0.1), border=NA)
    }

    for(i in whichTrees){    
     allPSR[[i]]<-res_PSR_ALL[[i]]$fitted_PSR
     allLTT[[i]]<-res_PSR_ALL[[i]]$fitted_LTT

      points( x     = res_PSR_ALL[[i]]$age_grid,
             y     = res_PSR_ALL[[i]]$fitted_PSR,
             main  = 'Fitted PSR',
             xlab  = 'age',
             ylab  = 'PSR',
             type  = 'l',
             lwd=0.2,
             #pch=".",
             col=grey(0.9,alpha=1)#0.5)#,
            ) 
    }
  #	# calc the MEDIAN of PSR
 # 	allPSR_ALL<-do.call(rbind, allPSR)
 # 	colnames(allPSR_ALL)<-res_PSR_ALL[[i]]$age_grid
 # 	allPSR_ALL_Med<-colMedians(allPSR_ALL)

  	# draw K-Pg
  	abline(v=66, lwd=2, lty=2)

  #	# draw the MEDIAN of PSR -- and PUSHED Lam of MEDIAN
  # ##################
  #	 points(x= res_PSR_ALL[[i]]$age_grid, y=allPSR_ALL_Med, type="l", lwd=3)
  #    
  #  	# draw the REAL LAMBDA
  #      for(q in c(1,3,5)){
  #  		points(x= timesToPlot, y=LAMBDA_bySQS[[q]], type="l", lty=1, lwd=3, col=lamCols[q+2])
  #  		text(x= 20, y= LAMBDA_bySQS[[q]][9]+0.2, labels=paste0("q=0.",q), pos=3, col=lamCols[q+2])
  #  	}
  #  		text(x= 20, y= LAMBDA_bySQS[[1]][9]+0.15, labels="pulled", pos=3, col="black")

    # draw the 95% CI of 100-tree PUSHED Lam results 
    ##############
      # from (fixedFossilLineages)
      #####
      pushedRates_ONLY<-allRates_quant_ALL[,1:4]

        # 95% CI of 100 trees
        DOWN<- BASE/2
      #  SUBSET<-13 # 60 Ma
      #  SUBSET<-16 # 75 Ma
      #   SUBSET<-15 # 70 Ma
        SUBSET<-23 # 110 Ma
        toPlot<-pushedRates_ONLY[1:SUBSET,]
        polygon(x = c(toPlot$time, rev(toPlot$time)), 
              y= c(toPlot$pushed_low95, rev(toPlot$pushed_high95))-DOWN, 
            col=rgb(1,0,0,alpha=0.3), border=NA)
        # median of 100 trees
        points(x= toPlot$time, y=toPlot$pushed_med-DOWN, type="l", lwd=3, col=rgb(1,0,0))

      axis(side=2, at=c(-DOWN, -DOWN+0.2, -DOWN+0.4), labels=c(0, 0.2, 0.4))

      # from (fixedFossilExtinction)
      #######
        # 95% CI of 100 trees
        DOWN<- BASE
      #  SUBSET<-13 # 60 Ma
      #  SUBSET<-16 # 75 Ma
      #  SUBSET<-15 # 70 Ma
      #  SUBSET<-23 # 110 Ma
        toPlot<-estHBD_Lam_quant_ALL[1:SUBSET,]
        polygon(x = c(toPlot$time, rev(toPlot$time)), 
              y= c(toPlot$estHBD_Lam_low95, rev(toPlot$estHBD_Lam_high95))-DOWN, 
            col=rgb(0,0,1,alpha=0.3), border=NA)
        # median of 100 trees
        points(x= toPlot$time, y=toPlot$estHBD_Lam_med-DOWN, type="l", lwd=3, col=rgb(0,0,1))

      axis(side=2, at=c(-DOWN, -DOWN+0.2, -DOWN+0.4), labels=c(0, 0.2, 0.4))

dev.off()






      


    # draw the individual 100-tree PUSHED Lam results
    ##############
    for(i in whichTrees){
      pushedRates_i<-LAMBDA_bySQS_ALL100[[i]]

      # 95% CI for EACH of 100 trees
      DOWN<- BASE
    #  SUBSET<-13 # 60 Ma
    #  SUBSET<-16 # 75 Ma
       SUBSET<-15 # 70 Ma
    #  SUBSET<-23 # 110 Ma
      toPlot<-pushedRates_i[1:SUBSET,]
      polygon(x = c(toPlot$time, rev(toPlot$time)), 
            y= c(toPlot$LAMBDA_bySQS_low95, rev(toPlot$LAMBDA_bySQS_high95))-DOWN, 
          col=rgb(1,0,0,alpha=0.1), border=NA)
      # median of 100 trees
      points(x= toPlot$time, y=toPlot$LAMBDA_bySQS_est-DOWN, type="l", lwd=0.6, col=grey(1,alpha=1))#rgb(1,0,0, alpha=0.2))
    }
    axis(side=2, at=c(-BASE, -BASE+0.2, -BASE+0.4), labels=c(0, 0.2, 0.4))

dev.off()




 