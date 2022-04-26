#!/usr/bin/env Rscript
library(UpSetR)
library(reshape2)
library(magrittr)
library(eulerr)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pgene = args[1]

#pgene = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/OverlapList/MC-E2up-P4up-H295R/upset_matrix"

dgene = read.csv(pgene, sep = "\t", header = TRUE)
rownames(dgene) = dgene[,1]
dgene = dgene[,-1]
dgene = dgene[,order(colnames(dgene))]

nb_gene = dim(dgene)[1]

png(paste(pgene, ".png", sep = ""), res = 300, 5000, 3000)
upset(dgene, nsets = nb_gene, matrix.color = "#DC267F", 
      main.bar.color = "#648FFF", sets.bar.color = "#FE6100")
dev.off()


png(paste(pgene, "_venn.png", sep = ""), res = 300, 1000, 1000)
plot(euler(dgene), quantities = TRUE)
dev.off()

svg(paste(pgene, "_venn.svg", sep = ""), 6, 6)
plot(euler(dgene), quantities = TRUE)
dev.off()
