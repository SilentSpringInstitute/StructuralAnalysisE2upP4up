#!/usr/bin/env Rscript
require(kohonen)
library(ggplot2)
source("loadPropList.R")
library(RColorBrewer)
library(viridis)
library(Toolbox)
#################
# Color reverse #
#################


coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}

colBrBG <- function(n){brewer.pal(n,"Blues")}

colViridis <- function(n){viridis(n)}

applySOM = function(som_model, d_prop, prop, pr_out, svg_plot = 0){
  
  #write cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(som_model$data[[1]])
  
  xdim = som_model$grid$xdim
  ydim = som_model$grid$ydim
  
  
  lAct = rep(0, xdim*ydim)
  names(lAct) = seq(1, xdim*ydim)
  
  # remove inactive
  if(d_prop != "0"){
    d_AC50 = as.vector(d_prop[,prop])
    names(d_AC50) = rownames(d_prop)
    d_AC50 = na.omit(d_AC50)
    del_nt = as.vector(which(d_AC50 == "NT"))
    if(is.integer0(del_nt) == FALSE){
      d_AC50 = d_AC50[-del_nt]  
    }
    del_neg = as.vector(which(d_AC50 == "NEG"))
    if(is.integer0(del_neg) == FALSE){
      d_AC50 = d_AC50[-del_neg]  
    }
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
  chemical.name = d_prop[rownames(dclust), c("Chemical.name")]
  dclust = cbind(chemical.name, dclust)
  dclust = dclust[order(as.double(as.character(dclust[,"Cluster"]))),]
  write.csv(dclust, paste(pr_out, "SOM_Clusters_", prop, ".csv", sep = ""))
  
  # count of active #
  ###################
  png(paste(pr_out, "SOM_count_", prop, ".png", sep = ""))
  plot(som_model, type = "property", property=lAct, palette.name=colViridis, main = "", dpi=300, height = 500, width = 500, bg = "transparent")
  dev.off()
  
  if(svg_plot == 1){
    svg(paste(pr_out, "SOM_count_", prop, ".svg", sep = ""))
    plot(som_model, type = "property", property=lAct, palette.name=colViridis, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
    dev.off()
  }
  
  
  # count of active calibrate #
  ##############################
  lAct[which(lAct == max(lAct))] = max(table(som_model$unit.classif))# have to calibrate based on the max of the original SOM
  if(svg_plot == 1){
    svg(paste(pr_out, "SOM_count_", prop, "_calibrate.svg", sep = ""))
    plot(som_model, type = "property", property=lAct, palette.name=colViridis, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
    dev.off()
  }
  
  # plot with proba #
  ###################
  png(paste(pr_out, "SOM_prob_", prop, ".png", sep = ""))
  plot(som_model, type = "property", property=lprob, palette.name=colViridis, main = "Prob active", dpi=300, height = 500, width = 500, bg = "transparent")
  dev.off()  
  
  if(svg_plot == 1){
    svg(paste(pr_out, "SOM_prob_", prop, ".svg", sep = ""))
    plot(som_model, type = "property", property=lprob, palette.name=colViridis, main = "Prob active", dpi=300, height = 20, width = 20, bg = "transparent")
    dev.off()  
  }
  
  
  write.csv(lprob, paste(pr_out, "SOM_Clusters_prob_", prop, ".csv", sep = ""))
  
}






################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_SOM_model = args[1]
p_prop = args[2]
pr_out = args[3]


#p_SOM_model = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/Analysis_H295R/rdkit-OPERA/SOM/SOM_model.RData"
#p_prop = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up//results/setOfChemicals/H295R.csv"
#pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/Analysis_H295R/rdkit-OPERA/SOM/"


# load prop
d_prop = loadProp(p_prop)


# load SOM model
load(p_SOM_model)

l_props = c("E2up", "P4up", "MC")

for(prop in l_props){
  d_apply = d_prop[,c(prop, "Chemical.name")]
  rownames(d_apply) = rownames(d_prop)
  applySOM(som_model, d_apply, prop, pr_out)
}



# take the overlap e2up and p4up
##############

E2up_P4up = rep("NT", dim(d_prop)[1])
E2up_P4up[which(d_prop$E2up == "POS" & d_prop$P4up == "POS")] = "POS" 
E2up_P4up[which(d_prop$E2up == "NEG" & d_prop$P4up == "NEG")] = "NEG" 
d_apply = cbind(d_prop$Chemical.name, E2up_P4up)
rownames(d_apply) = rownames(d_prop)
colnames(d_apply) = c("Chemical.name", "E2up_P4up")
d_apply = as.data.frame(d_apply)
applySOM(som_model, d_apply, "E2up_P4up", pr_out)



# take overlap E2up or P4up and MC
#######
E2up_P4up_MC = rep("NT", dim(d_prop)[1])
E2up_P4up_MC[which(d_prop$MC == "POS" & d_prop$P4up == "POS")] = "POS"
E2up_P4up_MC[which(d_prop$MC == "POS" & d_prop$E2up == "POS")] = "POS"
E2up_P4up_MC[which(d_prop$E2up == "NEG" & d_prop$P4up == "NEG" & d_prop$MC == "NEG")] = "NEG" 

d_apply = cbind(d_prop$Chemical.name, E2up_P4up_MC)
rownames(d_apply) = rownames(d_prop)
colnames(d_apply) = c("Chemical.name", "E2up_P4up_MC")
d_apply = as.data.frame(d_apply)
applySOM(som_model, d_apply, "E2up_P4up_MC", pr_out)

