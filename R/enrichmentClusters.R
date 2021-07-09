#!/usr/bin/env Rscript
library(ggpubr) 
source("loadPropList.R")
library(ggplot2)


enrich_by_prop = function(d_prop_in){
  
  l_prop = colnames(d_prop_in)
  nb_chem = dim(d_prop_in)[1]
  d_out = NULL
  for(prop in l_prop){
    nb_pos = length(which(d_prop_in[,prop] == "POS")) / nb_chem
    d_out = cbind(d_out, nb_pos)
  }
  d_out = cbind(d_out, nb_chem)
  colnames(d_out) = append(c(l_prop), "NB chem")
  return(d_out)
}


splitBalloonPlot = function(d_out, nb_page, pr_out){

  
  v_size_scaled = rep(table(d_clusters$cluster),dim(d_out)[2])/8
    
  nb_by_page = dim(d_out)[1] / nb_page
  for(i in 1:nb_page){
    
    d_page = d_out[(nb_by_page*(i-1)+1):(nb_by_page*i),]
    
    ggballoonplot(d_page, fill = "value", size = rep(v_size_scaled[(nb_by_page*(i-1)+1):(nb_by_page*i)],2))+
      scale_fill_viridis_c("% active", option = "C", limits=c(min(d_out), max(d_out)))
    ggsave(paste(pr_out, "cluster_proB_ballon", i, ".png",sep = ""), height = dim(d_page)[1]/8, width = 4, dpi = 300)
    
    

  }
  d_size = cbind(table(d_clusters$cluster), table(d_clusters$cluster))
  ggballoonplot(d_size, legend.title = "Size", fill="white")+
    guides(size=guide_legend("Size"))
  
  ggsave(paste(pr_out, "cluster_proB_ballon_legend.png",sep = ""), height = dim(d_page)[1]/8, width = 4, dpi = 300)
  
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_prop = args[1]
p_clusters = args[2]
pr_out = args[3]


#p_prop = "/mnt/c/Users/AlexandreBorrel/research/SSI/breast_carcinogen/results/setOfChemicals/H295R.csv"
# p_clusters = "/mnt/c/Users/AlexandreBorrel/research/SSI/breast_carcinogen/results/Analysis_H295R/rdkit-OPERA/SOM/SOM_Clusters"
#pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/breast_carcinogen/results/Analysis_H295R/rdkit-OPERA/SOM/"


#p_prop = "./../../results/setOfChemicals/H295R.csv"
#p_clusters = "./../../results/Analysis_H295R/rdkit/hclustDendo/cluster_hclust_ward2_gapstat.csv"
#pr_out = "./../../results/Analysis_H295R/rdkit/hclustDendo/"


prop_in = substr(p_prop, 1, nchar(p_prop)-4)
prop_in = tail(strsplit(prop_in, "/")[[1]], n=1)
d_clusters = read.csv(p_clusters, sep = ",", row.names = 1)


# open prop file and format
########
d_prop = loadProp(p_prop)



l_clusters = unique(d_clusters$cluster)

d_out = NULL
for(cluster in l_clusters){
  l_chem = d_clusters$names[which(d_clusters$cluster == cluster)]
  d_prop_temp = d_prop[l_chem,]
  d_out = rbind(d_out, enrich_by_prop(d_prop_temp))
}



rownames(d_out) = l_clusters
write.csv(d_out, paste(pr_out, "prob_by_clusters.csv", sep = ""))

# keep only E2up abd p4 up
drops = c("NB chem", "MC", "ER", "genotox", prop_in)
d_out = d_out[ , !(colnames(d_out) %in% drops)]
d_out = d_out[order(as.double(rownames(d_out))),]





# split ballon plot
nb_page = 2
splitBalloonPlot(d_out, nb_page, pr_out)
