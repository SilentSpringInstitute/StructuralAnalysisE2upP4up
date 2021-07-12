#!/usr/bin/env Rscript
library(kohonen)
library(ggplot2)


coolBlueHotRed <- function(n, alpha = 1) {
  rainbow(n, end=4/6, alpha=alpha)[n:1]
}



applySOM = function(som_model, d_prop, hormone, pr_out, svg_plot = 0){
  
  #write cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(som_model$data[[1]])
  
  xdim = som_model$grid$xdim
  ydim = som_model$grid$ydim
  
  d_sim_clust = NULL
  for(clust in seq(1, xdim*ydim)){
    l_casrn = names(which(dclust==clust))
    l_sim = d_prop[l_casrn]
    d_sim_clust = append(d_sim_clust, mean(as.double(l_sim)))
  }
  
  names(d_sim_clust) = seq(1, xdim*ydim)
  
  # plot with proba #
  ###################
  png(paste(pr_out, "_SOM_sim.png", sep = ""))
  plot(som_model, type = "property", property=d_sim_clust, palette.name=coolBlueHotRed, main = paste("Similarity", hormone), dpi=300, height = 500, width = 500, bg = "transparent")
  dev.off()  
  
  if(svg_plot == 1){
    svg(paste(pr_out, "_SOM_sim.png", sep = ""))
    plot(som_model, type = "property", property=d_sim_clust, palette.name=coolBlueHotRed, main = paste("Similarity", hormone), dpi=300, height = 500, width = 500, bg = "transparent")
    dev.off()  
  }
  
  return(d_sim_clust)
  
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_SOM_model = args[1]
p_hormone_sim = args[2]
pr_out = args[3]


#p_SOM_model = "/mnt/c/Users/AlexandreBorrel/research/SSI/breast_carcinogen/results/Analysis_H295R/rdkit-OPERA/SOM/SOM_model.RData"
#p_hormone_sim = "/mnt/c/Users/AlexandreBorrel/research/SSI/breast_carcinogen/results/similarityHormone//matrix_MACCS-Tanimoto.csv"
#pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/breast_carcinogen/results/Analysis_H295R/rdkit-OPERA/SOM_hormones/"


#p_SOM_model = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/Analysis_H295R/rdkit-OPERA/SOM/SOM_model.RData"
#p_hormone_sim = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/similarityHormone//matrix_MACCS-Tanimoto.csv"
#pr_out = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/Analysis_H295R/rdkit-OPERA/SOM_hormones/"

# load SOM model
load(p_SOM_model)

# define filout
f_hormone = basename(p_hormone_sim)
f_hormone = sub('\\..*$', '', basename(f_hormone))
pr_out = paste(pr_out, strsplit(f_hormone, "_")[[1]][2])

d_sim = read.csv(p_hormone_sim, sep = "\t", header = FALSE)
colnames(d_sim) = d_sim[1,]
rownames(d_sim) =d_sim[,1]
d_sim = d_sim[,-1]
d_sim = d_sim[-1,]

d_w = NULL
for(hormone in colnames(d_sim)){
  d_prop = d_sim[,hormone]
  names(d_prop) = rownames(d_sim)
  d_hor = applySOM(som_model, d_prop, hormone, paste(pr_out, hormone, sep = "_"))
  d_w = cbind(d_w, d_hor)
}
colnames(d_w) = colnames(d_sim) 
write.csv(d_w, paste(pr_out, "sim_matrix.csv", sep = ""))