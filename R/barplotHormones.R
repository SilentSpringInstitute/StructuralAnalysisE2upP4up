#!/usr/bin/env Rscript
library(ggplot2)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_matrixIn = args[1]

#p_matrixIn = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/steroidogenesis/single_hitc_matrix.csv"

d_h = read.csv(p_matrixIn, sep = "\t")
d_h = as.data.frame(d_h)

ggplot(d_h, aes(x=factor(SUM_HORMONE_CHANGED)))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  xlab("Number of hormones altered")+
  theme_minimal()

ggsave(paste(substr(p_matrixIn, 1, nchar(p_matrixIn)-4), "_barplot.png", sep = ""),  width = 10, height = 10, dpi = 300, bg="transparent")
