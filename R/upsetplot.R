#!/usr/bin/env Rscript
library(UpSetR)
library(reshape2)
library(magrittr)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pgene = args[1]

#pgene = "../../results/PFAS/Upset/upset_gene"

dgene = read.csv(pgene, sep = "\t", header = TRUE)
rownames(dgene) = dgene[,1]
dgene = dgene[,-1]
dgene = dgene[,order(colnames(dgene))]

nb_gene = dim(dgene)[1]

png(paste(pgene, ".png", sep = ""), res = 300, 5000, 3000)
upset(dgene, nsets = nb_gene, matrix.color = "#DC267F", 
      main.bar.color = "#648FFF", sets.bar.color = "#FE6100")
dev.off()