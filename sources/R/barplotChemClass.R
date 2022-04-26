#!/usr/bin/env Rscript
library(ggplot2)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_filin = args[1]
#p_filin = "/mnt/c/Users/aborr/research/Silent_Spring/breast_carcinogen/results/ChemClassInMC/count_chemical_class_MC.csv"



d_in = read.csv(p_filin, sep = "\t", header = TRUE)
rownames(d_in) = d_in$Chem_class
l_class = d_in$Chem_class
#l_class_remove = l_class[which(d_in$count < 5)]


#for(class_chem in l_class_remove){
#  d_in = d_in[-which(d_in$Chem_class == class_chem),]
#}

d_in$count = as.numeric(as.character(d_in$count))

d_in = d_in[order(d_in$count),]



# Barplot basique
p<-ggplot(data=d_in, aes(x=reorder(Chem_class, -count), y=count )) +
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("Chemical class")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Barplot horizontal
ggsave(paste(substr(p_filin, 1, nchar(p_filin)-4), "_barplot.png", sep = ""),  width = 5, height = 8, dpi = 300, bg="transparent")
