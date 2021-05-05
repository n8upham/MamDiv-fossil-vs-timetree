# Calculate tip DR across the 10k trees
###
# NS Upham; 8 Sept 2020
# ===================== 

library(ape)

# From Dan Greenberg
#this produces a caic-like structure from a newick string
readORDER<-function(file) {
	tree <- file
	tpc <- unlist(strsplit(tree, "[\\(\\),;]"))
	tpc=tpc[grep(":",tpc)]

	# find the tip and edge labels
	tiplabels=tpc[grep(".:",tpc)]
	edgelabels=tpc[-grep(".:",tpc)]
	edgelabels=c(edgelabels,":0")

	# locate the clusters and edges
	tree <- unlist(strsplit(tree, NULL))
	x=which(tree=="(")
	y=which(tree==")")
	v=which(tree==":")
	#these are the locations of the tips
	w=setdiff(v,y+1) 

	##Pass through string from left to right locating the paired parenthesis, thus
	## allowing for easy assigment of edge to subtending tips
	#initialize objects (M is the actual clade matrix, while E is vector with associated edge weights)
	j=2 
	k=length(w)+1
	M=matrix(0,length(w)*2-1,length(w))
	E=as.vector(matrix(0,1,length(w)*2-1))
	x=c(x,y[length(y)]+1)
	# Main Pass
	while (length(x)>1)
	{
		if (x[j+1]<y[1])
	{j=j+1
		} else {
		M[k,which(x[j]<w & w<y[1])]=j-1
		E[k]=strsplit(edgelabels[k-length(w)],"[:]")[[1]][2]
		k=k+1
		y=y[-1]
		x=x[-j]
		j=j-1}
	}

	# Assign branch lengths and finished tip names to the tips
	for (i in 1:length(w))
	{
		tmp=strsplit(tiplabels[i],"[:]")[[1]]
		E[i]=tmp[2]
		tiplabels[i]=tmp[1]
		tmp50=abs(M[,i]-max(M[,i])-1)
		M[,i]=(M[,i]>0)*tmp50
	}
		M=list(M,as.numeric(E),tiplabels)
}


# from Dan Greenberg::
#####
ES_v2<- function (clade.matrix,edge.length,tip.label) {
# second way
	caicmatrix=clade.matrix
	caicmatrix[clade.matrix>1]=1
	caicmatrix[(1:dim(caicmatrix)[2]),(1:dim(caicmatrix)[2])]=diag(dim(caicmatrix)[2])
	lambda=edge.length*caicmatrix
	rm(caicmatrix)
	tmp=lambda/(2^clade.matrix)

	ESS=as.matrix(colSums(tmp))
	rownames(ESS)=tip.label
	return(ESS)
}


# TEST these on a 10-tip tree
########
#options(digits = 10)
#library(phytools); library(picante)
#simtree<-rtree(10)
#write.tree(simtree,"testTree_10tips.tre")
#
#trees1=scan("testTree_10tips.tre", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
#trees=strsplit(trees1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree
#treeOrig<-read.tree(text=trees1)
#
#	caicTree_ES<-readORDER(trees[[1]])
##
#	tree<-caicTree_ES
#	res_ES = ES_v2(clade.matrix=tree[[1]],edge.length=tree[[2]],tip.label=tree[[3]])  #test.new gives the tip name as the rowname and ED measure in column 1.
#	ES<-as.matrix(res_ES[match(simtree$tip.label, rownames(res_ES)),])	
#
#	# with picante functions:
#	res_ED_es <- evol.distinct(treeOrig, type = "equal.splits",# "fair.proportion"),
#         scale = FALSE, use.branch.lengths = TRUE)
#	ED_es<-as.matrix(res_ED_es[match(simtree$tip.label, as.character(res_ED_es$Species)),"w"])	
##
#	# compare:
#	compare<-cbind.data.frame(ES=ES,ED_es=ED_es)
#
#		#>> OK great, this all checks out.
