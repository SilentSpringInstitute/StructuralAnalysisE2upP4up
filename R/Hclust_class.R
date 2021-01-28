#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)

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
  
  t4
  
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


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]

valcor = args[2]
maxquantile = as.double(args[3])

valcor = 0.9
maxquantile = 90
pdesc = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/carcinogen_breast_011321-H295_E2-up-H295_P4-up-nonCarci_list-Judson2015_ER_agonist-Judson2015_ER_antagonist/desc1D2D.csv"
prout = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/carcinogen_breast_011321-H295_E2-up-H295_P4-up-nonCarci_list-Judson2015_ER_agonist-Judson2015_ER_antagonist/"


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

dendogramCircleClass (ddes, daff, prout)
