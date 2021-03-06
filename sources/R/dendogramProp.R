#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)
library(RootsExtremaInflections)

library(ggplot2)
library(Toolbox)


dendogramCluster = function(ddes, d_cluster, daff, dh, prout){
  
  
  #calibrate affinity for color
  
  daff = daff[,c("MC", "E2up",  "P4up")]
  
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
  
  t3 = t2 + ggnewscale::new_scale_fill()
  
  t4 = gheatmap(t3, dh, offset=25, width=.2, colnames_angle=95, colnames_offset_y = .25) +
    scale_fill_viridis_b(option = "A")
   
  
  open_tree(t4, 10) %>% rotate_tree(10)
  ggsave(pfilout, dpi=300, height = 11, width = 11)
}




clusterChem = function(desc, prout){
  
  # scale data in input
  din = scale (desc)

  print(dim(din))
  
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
p_hormone = args[3]
pr_out = args[4]
val_cor = as.double(args[5])
max_q = as.integer(args[6])

#val_cor = 0.9
#max_q = 90

#p_prop = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/setOfChemicals/H295R.csv"
#p_desc = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/Analysis_H295R/rdkit/rdkit.csv"
#p_hormone = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/similarityHormone/matrix_MACCS-Tanimoto.csv"
#pr_out = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/Analysis_H295R/rdkit/hclustDendo/"


d_in = read.csv(p_prop, sep = "\t", row.names = 1)
drops = c("Group", "Aff", "Chemical.name", "SMILES")
d_prop = d_in[ , !(names(d_in) %in% drops)]

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


#hormone
d_hormone = read.csv(p_hormone, sep = "\t", header = FALSE)
rownames(d_hormone) = d_hormone[,1]
d_hormone = d_hormone[,-1]
colnames(d_hormone) = d_hormone[1,]

# reduce with eastradiol and progesterone
d_hormone = d_hormone[,c("57-83-0", "50-28-2")]

# open and reduce descriptor
l_d_desc = openDataVexcluded(p_desc, val_cor, pr_out,c(1,2))
d_desc = l_d_desc[[1]]
rownames(d_desc) = d_desc[,1]
l_casrn = rownames(d_desc)
d_desc = d_desc[,-c(1,2)]
d_desc = apply(d_desc, 2, as.double)
d_desc = delnohomogeniousdistribution(d_desc, max_q)
rownames(d_desc) = l_casrn
d_desc = as.data.frame(d_desc)

l_casrn = intersect(rownames(d_desc), rownames(d_prop))

d_prop = d_prop[l_casrn,]
d_desc = d_desc[l_casrn,]
d_hormone = d_hormone[l_casrn,]
d_hormone = sapply(d_hormone, as.double, 1)
rownames(d_hormone) = l_casrn

#cluster chem for representation
d_cluster = clusterChem(d_desc, pr_out)
chemical_name = d_in[rownames(d_cluster), "Chemical.name"]
d_w = cbind(d_cluster, chemical_name)
write.csv(d_w, paste(pr_out, "cluster_hclust_ward2_gapstat_name.csv", sep = ""), row.names = FALSE)


# dendogram
dendogramCluster(d_desc, d_cluster, d_prop, d_hormone, pr_out)

