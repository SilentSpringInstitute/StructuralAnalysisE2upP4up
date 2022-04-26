#!/usr/bin/env Rscript
library(ggplot2)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_chemlist = args[1]
#p_chemlist = "c://Users/aborr/research/Silent_Spring/jenny_carcinogen/results/Carcinogen_list/chem_list/count_chemlist"



d_chemlist = read.csv(p_chemlist, sep = "\t", header = TRUE)
rownames(d_chemlist) = d_chemlist$chem_list
l_class = d_chemlist$chem_list
l_class_remove = l_class[which(d_chemlist$count < 20)]


for(class_chem in l_class_remove){
  d_chemlist = d_chemlist[-which(d_chemlist$chem_list == class_chem),]
}

d_chemlist$count = as.numeric(as.character(d_chemlist$count))

d_chemlist = d_chemlist[order(d_chemlist$count),]



# Barplot basique
p<-ggplot(data=d_chemlist, aes(x=reorder(chem_list, -count), y=count )) +
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("ToxPrint")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Barplot horizontal
ggsave(paste(p_chemlist, "_barplot.png", sep = ""),  width = 8, height = 14, dpi = 300, bg="transparent")
