#!/usr/bin/env Rscript
library(VennDiagram)


venPlot = function(xin, colfill, prout){
  
  venn.diagram(xin, 
               filename = paste(prout, ".tiff", sep = ""),
               col = "black",
               lty = "dotted",
               lwd = 4,
               fill = colfill,
               alpha = 0.50,
               cex = 1.5,
               fontfamily = "serif",
               fontface = "bold",
               cat.col = colfill,
               cat.cex = 0.5,
               cat.fontfamily = "serif")
}

################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_filin = args[1]

#p_filin = "c:/Users/aborr/research/Silent_Spring/breast_carcinogen/results/Overlap_ER-carcinogen/upset_matrix"
d_in = read.csv(p_filin, sep = "\t", row.names = 1)

# format data for venn diagragram
d_plot = NULL
for(i in 1:dim(d_in)[2]){
  d_plot[[i]] = which(d_in[,i] == 1) 
  
}
names(d_plot) = colnames(d_in)


venPlot(d_plot, "blue", paste(p_filin, "_vennPlot", sep = ""))
