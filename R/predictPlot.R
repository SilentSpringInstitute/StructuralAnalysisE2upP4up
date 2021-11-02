#!/usr/bin/env Rscript
library(ggplot2)
library(hrbrthemes)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_pred = args[1]
pr_out = args[2]


p_pred = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/predMC_E2/ToxPrint_QSAR_filtered_AD_0.75_ToxPrint_3_QSAR_0.5.csv"
pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/predMC_E2/"

d_pred = read.csv(p_pred, sep = "\t")

d_pred$AD.distance.to.first.neighbord = as.double(as.character(d_pred$AD.distance.to.first.neighbord))
d_pred$nb.Toxprint...... = as.double(as.character(d_pred$nb.Toxprint......))
d_pred$Pred.RF.balanced = as.double(as.character(d_pred$Pred.RF.balanced))

ggplot(d_pred, aes(x=nb.Toxprint......, y=Pred.RF.balanced, color=AD.distance.to.first.neighbord)) + 
  geom_point(size=6) +
  labs(colour = "AD", x="Nb ToxPrint(***)", y = "Prediction with balanced RF")

ggsave(paste(pr_out, "pred_plot.png" ,sep = ""), dpi = 300, width = 8)


ggplot(d_pred, aes(x=nb.Toxprint......, y=Pred.RF.balanced, color=AD.distance.to.first.neighbord, label = CASRN)) + 
  geom_text() +
  labs(colour = "AD", x="Nb ToxPrint(***)", y = "Prediction with balanced RF")

ggsave(paste(pr_out, "pred_text_plot.png" ,sep = ""), dpi = 300, width = 8)
