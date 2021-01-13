#!/usr/bin/env Rscript
library(ggplot2)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_toxprint = args[1]
#p_toxprint = "c://Users/aborr/research/Silent_Spring/jenny_carcinogen/results/Carcinogen_list/ToxPrint/count_toxprint"



d_toxprint = read.csv(p_toxprint, sep = "\t", header = TRUE)
rownames(d_toxprint) = d_toxprint$Toxprint
l_class = d_toxprint$Toxprint
l_class_remove = l_class[which(d_toxprint$count < 5)]


for(class_chem in l_class_remove){
  d_toxprint = d_toxprint[-which(d_toxprint$Toxprint == class_chem),]
}

d_toxprint$count = as.numeric(as.character(d_toxprint$count))

d_toxprint = d_toxprint[order(d_toxprint$count),]



# Barplot basique
p<-ggplot(data=d_toxprint, aes(x=reorder(Toxprint, -count), y=count )) +
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("ToxPrint")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Barplot horizontal
ggsave(paste(p_toxprint, "_barplot.png", sep = ""),  width = 8, height = 14, dpi = 300, bg="transparent")
