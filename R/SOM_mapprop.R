#!/usr/bin/env Rscript
require(kohonen)
library(ggplot2)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_SOM_model = args[1]
p_prop = args[2]
pr_out = args[3]


p_SOM_model = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/clusterMC/analysis_desc/rdkit/SOM/SOM_model.RData"
p_prop = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/clusterMC/MC/dataset.csv"

d_prop = read.csv(p_prop, sep = "\t", row.names = 1)
d_prop = d_prop[which(d_prop$MC == 1),]



# load SOM model
load(p_SOM_model)

