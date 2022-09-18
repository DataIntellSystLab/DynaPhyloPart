rm(list=ls(all=TRUE))

library(Rcpp)
library(ape)
library(castor)
library(phytools)

setwd("C:/Users/m.prosperi/Downloads/PhyloPart_v2.1")
tree=read.tree("nextstrain_ncov_open_global_all-time_tree.nwk")

#tree=rtree(333)

minpercthresh=0.25
minnumclus=3
dist_limit=25000

nallnodes=Ntip(tree)+Nnode(tree)
if (nallnodes*(nallnodes-1)/2<dist_limit) {dist_limit=nallnodes*(nallnodes-1)/2}

dist_mat=NULL
samplenodes1=sample(1:nallnodes)
samplenodes2=sample(1:nallnodes)
samplenodes1=1:nallnodes
samplenodes2=1:nallnodes
lim=0
for (i in samplenodes1)
{
	for (j in samplenodes2)
	{
		if ( i!=j && ( !(length(which(dist_mat[,1]==i))>0 && length(which(dist_mat[,2]==j))>0) || !(length(which(dist_mat[,1]==j))>0 && length(which(dist_mat[,2]==i))>0) ) )
		{
			dab=get_pairwise_distances(tree, i, j)
			dist_mat=rbind(dist_mat,c(i,j,dab))
			lim=lim+1
			if (lim%%6666==0) {cat(round(100*lim/dist_limit,1),"%\n",sep="")}
			if (lim==dist_limit) break;
		}
	}
	if (lim==dist_limit) break;
}
dist_mat=na.omit(dist_mat)
thres=quantile(dist_mat[,3],minpercthresh)


innernodes=unique(tree$edge[,1])
clusts=list()
c=0
for (i in innernodes)
{
	found=F
	if (length(clusts)>0) {for (j in 1:length(clusts)) {if (i %in% clusts[[j]]) {found=T; break;}}}
	if (!found)
	{
		desc=getDescendants(tree,i);
		clu=i
		dis=NULL
		for (k in desc)
		{
			ds=get_pairwise_distances(tree,i,k)
			if (!is.na(ds)) {if (ds<=thres) {clu=c(clu,k)}}
		}
		c=c+1
		clusts[[c]]=clu
		if (i%%10==0) cat("node",i,"\n")
	}
}

nclus=0
for (i in 1:length(clusts)) {if (length(clusts[[i]])>=minnumclus) nclus=nclus+1}

clusts
nclus

palett=rainbow(length(clusts)); #rainbow; topo.colors; heat.colors; cm.colors; terrain.colors; colors()
palett=sample(palett)

edgecols=rep("grey",length(tree$edge[,2]))
for (i in 1:length(edgecols))
{
	for (j in 1:length(clusts)) {if (length(clusts[[j]])>=minnumclus && tree$edge[i,2] %in% clusts[[j]]) {edgecols[i]=palett[j];break;}}
}
nodecols=rep("grey",length(unique(tree$edge[,1])))
for (i in 1:length(nodecols))
{
	for (j in 1:length(clusts)) {if (length(clusts[[j]])>=minnumclus && unique(tree$edge[,1])[i] %in% clusts[[j]]) {nodecols[i]=palett[j];break;}}
}
tipcols=rep("grey",Ntip(tree))
for (i in 1:length(tipcols))
{
	for (j in 1:length(clusts)) {if (length(clusts[[j]])>=minnumclus && i %in% clusts[[j]]) {tipcols[i]=palett[j];break;}}
}

plot(tree,show.tip.label=F,type="f",edge.color=edgecols)
#tiplabels(frame="none",col=tipcols,cex=0.8)
#nodelabels(frame="none",col=nodecols,cex=0.8)
