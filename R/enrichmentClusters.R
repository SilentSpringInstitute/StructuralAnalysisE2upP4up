#!/usr/bin/env Rscript
library(ggpubr) 


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





################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_prop = args[1]
p_clusters = args[2]
pr_out = args[3]


#p_prop = "./../../results/setOfChemicals/H295R.csv"
#p_clusters = "./../../results/Analysis_H295R/rdkit/hclustDendo/cluster_hclust_ward2_gapstat.csv"
#pr_out = "./../../results/Analysis_H295R/rdkit/hclustDendo/"
#d_prop_temp = "H295R"


prop_in = substr(p_prop, 1, nchar(p_prop)-4)
prop_in = tail(strsplit(prop_in, "/")[[1]], n=1)
d_clusters = read.csv(p_clusters, sep = ",", row.names = 1)

d_in = read.csv(p_prop, sep = "\t", row.names = 1)
drops = c("Group", "Aff", "Chemical.name", "SMILES")
d_prop = d_in[ , !(names(d_in) %in% drops)]

###########
### prep d_prop
#
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


#forH295R
d_prop$MC[which(d_prop$MC == "" | is.na(d_prop$H295R))] = "NT"
d_prop$MC[which(d_prop$MC == "0")] = "NEG"
d_prop$MC[which(d_prop$MC == "1")] = "POS"


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


v_size = rep(table(d_clusters$cluster),dim(d_out)[2])/8


png(paste(pr_out, "cluster_proB_ballon.png",sep = ""), height = dim(d_out)[1]*10, width = 300)
ggballoonplot(d_out, fill = "value", size = v_size)+
  scale_fill_viridis_c(option = "C")
dev.off()
