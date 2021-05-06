#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)
library(RootsExtremaInflections)
library(NbClust)

library(ggplot2)
library(Toolbox)


dendogramCluster = function(ddes, d_cluster, daff, prout){
  
  
  #calibrate affinity for color
  
  daff = daff[,c("MC", "genotox", "E2up",  "P4up" , "ER")]
  
  daff = as.data.frame(daff)
  #daff = cbind(rownames(daff), daff)
  
  dcluster = cbind(rownames(ddes), d_cluster[rownames(ddes),])
  #cluster = rep("1", dim(ddes)[1])
  #dcluster = cbind(rownames(ddes), cluster)
  dcluster = as.data.frame(dcluster)
  rownames(dcluster) = as.character(rownames(ddes))
  colnames(dcluster) = c("casrn", "cluster")
  
  
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
  t2 = gheatmap(t1, daff, offset=10, width=.2, colnames_angle=95, colnames_offset_y = .25) +
    scale_fill_manual(
      values = c("blue", "gray80", "orange"),
      labels = c("NEG", "NT", "POS")
    )
    open_tree(t2, 10) %>% rotate_tree(10)
  
  ggsave(pfilout, dpi=300, height = 11, width = 11)
}



clusterChem = function(d_in, prout){
  
  # scale data in input
  d_dis = 1- as.dist(d_in)

  # use nbclusster
  din = NbClust(data=NULL, diss = d_dis, min.nc = 2, max.nc = 15, method = "ward.D2", distance = NULL, index="silhouette")
  
  
  dcluster = as.data.frame(din$Best.partition)
  
  # save cluster
  write.csv(dcluster, paste(prout, "cluster_hclust_ward2_gapstat.csv", sep = ""), row.names = TRUE)
  
  
  return(dcluster)
  
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_prop = args[1]
p_desc = args[2]
pr_out = args[3]

#p_prop = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/setOfChemicals/MC.csv"
#p_desc = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/tanimoto_by_list/MC.csv"
#pr_out = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/Analysis_MC/ToxPrint/hclustDendo/"

d_prop = read.csv(p_prop, sep = "\t", row.names = 1)
drops = c("Group", "Aff", "Chemical.name", "SMILES")
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
d_prop$ER[which(d_prop$ER == "tested")] = "NEG"
d_prop$ER[which(d_prop$ER == "agonist" | d_prop$ER == "antagonist")] = "POS"

#genotox
d_prop$genotox[which(d_prop$genotox == "not tested" | d_prop$genotox == "predicted genotoxic" | d_prop$genotox == "" | is.na(d_prop$genotox))] = "NT"
d_prop$genotox[which(d_prop$genotox == "not genotoxic")] = "NEG"
d_prop$genotox[which(d_prop$genotox != "NEG" & d_prop$genotox != "NT")] = "POS"


#forMC
d_prop$MC[which(d_prop$MC == "" | is.na(d_prop$MC))] = "NT"
d_prop$MC[which(d_prop$MC == "0")] = "NEG"
d_prop$MC[which(d_prop$MC == "1")] = "POS"


# open matrix #
################
d_matrix = read.csv(p_desc, sep = "\t", row.names = 1)
colnames(d_matrix) = rownames(d_matrix)

l_casrn = intersect(rownames(d_matrix), rownames(d_prop))

d_prop = d_prop[l_casrn,]
d_matrix = d_matrix[l_casrn,l_casrn]

#cluster chem for representation
d_cluster = clusterChem(d_matrix, pr_out)

# dendogram
dendogramCluster(d_matrix, d_cluster, d_prop, pr_out)

