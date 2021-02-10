#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)
library(randomcoloR)

library(Toolbox)


dendogramCircleClass = function(ddes, daff, prout){
  #calibrate affinity for color
  
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  legend.title = "Weight"
  
  
  # all chemicals
  pfilout = paste(prout, "HClust_dendo_class.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=class_chem, label=label, angle=angle), hjust=-0.5, size=2.5) +
    geom_tippoint(aes(color=class_chem), alpha=0.75, size=0.5)+
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    guides(col=guide_legend("List chemicals"))+
    geom_treescale(x = 1, y = 1, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  ggsave(pfilout, dpi=300, height = 22, width = 28)
  
  
  
  pfilout = paste(prout, "HClust_dendo_group.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=Aff, label=Aff, angle=angle), vjust=1,hjust=0,  size=3) +
    geom_tippoint(aes(color=Aff), alpha=0.75, size=0.5)+
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    guides( color = FALSE)+
    geom_treescale(x = 1, y = 1, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  
  ggsave(pfilout, dpi=300, height = 38, width = 45, units="cm")
  
}



PCAplot = function (din, dcolor, prout){
  coord = generatePCAcoords(din)
  data_plot = coord[[1]]
  var_cap = coord[[2]]
  cp = coord[[3]]
  
  col.desc = "black"
  palette <- distinctColorPalette(length(unique(dcolor$Aff)))
  dcolor$Aff = as.integer(as.character(dcolor$Aff))
  
  png (paste (prout, "PCA_text.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  
  color_arrow = col.desc[rownames(cp)]
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = paste("PC1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("PC2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  text (data_plot[,1],data_plot[,2], label = rownames (din), cex = 1.2)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  
  # PCA color  
  png (paste (prout, "PCA_color.png", sep = ""),  1700, 1500, res = 100)
  factor = factorACP (data_plot, cp)
  color_arrow =col.desc
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, col = palette[dcolor$Aff], xlab = paste("PC1: ", round (var_cap[1], 2), "%", sep = ""), ylab = paste("PC2: ", round(var_cap[2], 2), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
  legend(12, -8, legend=unique(dcolor$class_chem), col=palette[unique(dcolor$Aff)], pch=19, cex=0.8)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  
  png (paste (prout, "PCA_descriptor.png", sep = ""), 1700, 1500)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], xlab = paste("PC1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("PC2: ", signif (var_cap[2], 4), "%", sep = ""), pch=20, cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 4 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 2.5)
  dev.off()
  
}





################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]

valcor = args[2]
maxquantile = as.double(args[3])

valcor = 0.90
maxquantile = 90
pdesc = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/carcinogen_breast_012721-H295_E2-up-H295_P4-up-Judson2015_ER_agonist-Judson2015_ER_antagonist/desc1D2D.csv"
prout = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/carcinogen_breast_012721-H295_E2-up-H295_P4-up-Judson2015_ER_agonist-Judson2015_ER_antagonist/"


dglobal = openDataVexcluded(pdesc, valcor, prout, c(1,2, 577))
dglobal = dglobal[[1]]
rownames(dglobal) = dglobal[,1]
dglobal = as.data.frame(dglobal)
ddes = dglobal[,c(-1,-2,-3)]
ddes <- sapply(ddes,as.numeric)
rownames(ddes) = rownames(dglobal)

Aff = as.numeric(as.factor(dglobal$dataset))
class_chem= dglobal$dataset
id = rownames(ddes)
daff = cbind(id, Aff, class_chem)
daff = as.data.frame(daff)

# hclust
dendogramCircleClass (ddes, daff, prout)

# PCA
PCAplot(ddes, daff, prout)

