#!/usr/bin/env Rscript
require(kohonen)
library(ggplot2)
source("loadPropList.R")

#################
# Color reverse #
#################

colors <- function(n, alpha = 1) {
  rev(heat.colors(n, alpha))
}

coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}



applySOM = function(som_model, d_AC50, prop, pr_out, svg_plot = 0){
  
  #write cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(som_model$data[[1]])
  
  xdim = som_model$grid$xdim
  ydim = som_model$grid$ydim
  
  
  lAct = rep(0, xdim*ydim)
  names(lAct) = seq(1, xdim*ydim)
  
  # remove inactive
  if(d_AC50 != "0"){
    d_AC50 = na.omit(d_AC50)
    d_AC50 = d_AC50[-which(d_AC50 == "NEG")]
    d_AC50 = d_AC50[-which(d_AC50 == "NT")]
    dclust = dclust[names(d_AC50)]# take only active chemical
    dclust = cbind(names(dclust), dclust)
    ltable = table(dclust[,2])
    lAct[names(ltable)] = ltable
    colnames(dclust) = c("CASRN", "Cluster")
  }else{
    return() # to not run because no classification
  }
  
  
  ltabinit = table(som_model$unit.classif)
  linitial = rep(0, xdim*ydim)
  names(linitial) = seq(1, xdim*ydim)
  linitial[names(ltabinit)] = ltabinit
  
  lprob = lAct / linitial
  write.csv(dclust, paste(pr_out, "SOM_Clusters_", prop, ".csv", sep = ""))
  
  # count of active #
  ###################
  png(paste(pr_out, "SOM_count_", prop, ".png", sep = ""))
  plot(som_model, type = "property", property=lAct, palette.name=coolBlueHotRed, main = "", dpi=300, height = 500, width = 500, bg = "transparent")
  dev.off()
  
  if(svg_plot == 1){
    svg(paste(pr_out, "SOM_count_", prop, ".svg", sep = ""))
    plot(som_model, type = "property", property=lAct, palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
    dev.off()
  }
  
  
  # count of active calibrate #
  ##############################
  lAct[which(lAct == max(lAct))] = max(table(som_model$unit.classif))# have to calibrate based on the max of the original SOM
  if(svg_plot == 1){
    svg(paste(pr_out, "SOM_count_", prop, "_calibrate.svg", sep = ""))
    plot(som_model, type = "property", property=lAct, palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
    dev.off()
  }
  
  # plot with proba #
  ###################
  png(paste(pr_out, "SOM_prob_", prop, ".png", sep = ""))
  plot(som_model, type = "property", property=lprob, palette.name=coolBlueHotRed, main = "Prob active", dpi=300, height = 500, width = 500, bg = "transparent")
  dev.off()  
  
  if(svg_plot == 1){
    svg(paste(pr_out, "SOM_prob_", prop, ".svg", sep = ""))
    plot(som_model, type = "property", property=lprob, palette.name=coolBlueHotRed, main = "Prob active", dpi=300, height = 20, width = 20, bg = "transparent")
    dev.off()  
  }
  
  write.csv(lprob, paste(pr_out, "SOM_Clusters_prob_", prop, sep = ""))
  
}






################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_SOM_model = args[1]
p_prop = args[2]
pr_out = args[3]


#p_SOM_model = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/Analysis_H295R/rdkit-OPERA/SOM/SOM_model.RData"
#p_prop = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/setOfChemicals/H295R.csv"
#pr_out = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/Analysis_H295R/rdkit-OPERA/SOM/"

# load prop
d_prop = loadProp(p_prop)


# load SOM model
load(p_SOM_model)

l_props = c("E2up", "P4up")

for(prop in l_props){
  d_apply = d_prop[,prop]
  names(d_apply) = rownames(d_prop)
  applySOM(som_model, d_apply, prop, pr_out)
}
