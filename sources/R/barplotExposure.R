#!/usr/bin/env Rscript
library(ggplot2)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_exposure = args[1]
  
d_exp = read.csv(p_exposure, sep = "\t", header = TRUE)
rownames(d_exp) = d_exp$Board.exposure
l_class = d_exp$Board.exposure


d_exp$Count = as.numeric(as.character(d_exp$Count))
d_exp = d_exp[order(d_exp$Count),]

# Barplot basique
p<-ggplot(data=d_exp, aes(x=reorder(Board.exposure, Count), y=Count )) +
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("Board exposure")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Barplot horizontal
ggsave(paste(p_exposure, "_barplot.png", sep = ""),  width = 5, height = 5, dpi = 300, bg="transparent")
