#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)
library(RootsExtremaInflections)

library(ggplot2)
library(Toolbox)


dendogramCluster = function(ddes, d_cluster, daff, prout){
  
  
  #calibrate affinity for color
  
  daff = daff[,c("genotox", "E2up",  "P4up" , "ER")]
  
  daff = as.data.frame(daff)
  #daff = cbind(rownames(daff), daff)
  
  dcluster = d_cluster[rownames(ddes),]
  #cluster = rep("1", dim(ddes)[1])
  #dcluster = cbind(rownames(ddes), cluster)
  #dcluster = as.data.frame(dcluster)
  #rownames(dcluster) = as.character(rownames(ddes))
  
  
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  #tupgma2 <- groupOTU(tupgma2, 35)
  
  pfilout = paste(prout, "dendo_cluster_name.png", sep = "")
  
  t1 <- ggtree(tupgma2, layout="circular", size=0.8, aes(color=cluster))
  t1 <- t1 %<+% dcluster + geom_text(aes(color=cluster, label = label, angle=angle,  fontface="bold"), hjust=-0.10, size=1.2) +
    geom_tippoint(aes(color=cluster), alpha=0.75, size=1)+
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
    geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  #daff = daff[,-1]
  t2 = gheatmap(t1, daff, offset=5, width=.2, colnames_angle=95, colnames_offset_y = .25) +
    scale_fill_manual(
      values = c("blue", "gray80", "orange"),
      labels = c("NEG", "NT", "POS")
    )
    open_tree(t2, 10) %>% rotate_tree(10)
  
  ggsave(pfilout, dpi=300, height = 11, width = 11)
}



clusterChem = function(desc, prout){
  
  # scale data in input
  din = scale (desc)
  
  # hclust with ward D2 and a gap stat approach 
  p = fviz_nbclust(din, hcut, hcut_metho = "ward.D2",  method = "gap_stat", k.max = 50)
  
  
  dcluster = as.matrix(p$data)
  d = inflexi(as.double(dcluster[,5]),-1*as.double(dcluster[,6]),1,length(dcluster[,1]),3,3,plots=FALSE)
  nboptimal = d$finfl[1]
  
  outclust = hcut(din, k = nboptimal, hc_method = "ward.D2")
  
  # define cluster
  dcluster2 = cbind(names(outclust$cluster),outclust$cluster)
  colnames(dcluster2) = c("names", "cluster")
  dcluster2 = as.data.frame(dcluster2)
  rownames(dcluster2) = names(outclust$cluster)
  dcluster2 = dcluster2[rownames(din),]
  
  # save cluster
  write.csv(dcluster2, paste(prout, "cluster_hclust_ward2_gapstat.csv", sep = ""), row.names = TRUE)
  
  
  return(dcluster2)
  
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_prop = args[1]
p_desc = args[2]
pr_out = args[3]

val_cor = 0.9
max_q = 90


p_prop = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/clusterMC/all/dataset.csv"
p_desc = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/clusterMC/desc_1D2D.csv"
pr_out = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/clusterMC/"


d_prop = read.csv(p_prop, sep = "\t", row.names = 1)
d_prop = d_prop[which(d_prop$MC == 1),]

drops = c("Group", "Aff")
d_prop = d_prop[ , !(names(d_prop) %in% drops)]

# for E2up
d_prop$E2up[which(is.na(d_prop$E2up))] = "NT"
d_prop$E2up[which(d_prop$E2up == "ns effect")] = "NEG"
d_prop$E2up[which(d_prop$E2up != "NEG" & d_prop$E2up != "NT")] = "POS"

#for P4up
d_prop$P4up[which(is.na(d_prop$P4up))] = "NT"
d_prop$P4up[which(d_prop$P4up == "ns effect")] = "NEG"
d_prop$P4up[which(d_prop$P4up != "NEG" & d_prop$P4up != "NT")] = "POS"

#for ER
d_prop$ER[which(d_prop$ER == "")] = "NT"
d_prop$ER[which(d_prop$ER == "agonist")] = "POS"

#genotox
d_prop$genotox[which(d_prop$genotox == "not tested" | d_prop$genotox == "predicted genotoxic")] = "NT"
d_prop$genotox[which(d_prop$genotox == "not genotoxic")] = "NEG"
d_prop$genotox[which(d_prop$genotox != "NEG" & d_prop$genotox != "NT")] = "POS"


# open and reduce descriptor
l_d_desc = openDataVexcluded(p_desc, val_cor, pr_out,c(1,2))
d_desc = l_d_desc[[1]]
rownames(d_desc) = d_desc[,1]
l_casrn = rownames(d_desc)
d_desc = d_desc[,-c(1,2)]
d_desc = apply(d_desc, 2, as.double)
d_desc = delnohomogeniousdistribution(d_desc, max_q)
rownames(d_desc) = l_casrn

l_casrn = intersect(rownames(d_prop), rownames(d_desc))
d_prop = d_prop[l_casrn,]
d_desc = d_desc[l_casrn,]


#cluster chem for representation
d_cluster = clusterChem(d_desc, pr_out)

# dendogram
dendogramCluster(d_desc, d_cluster, d_prop, pr_out)

