#!/usr/bin/env Rscript
library(ggplot2)




################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc = args[1]
p_desc_cleaned = args[2] 
pr_out = args[3]


#p_desc = "../../results/Phthalates/DESC/desc_1D2D.csv"
#p_desc_cleaned = "../../results/Phthalates/RDKIT_desc/Cleaned_Data/desc1D2D_cleaned.csv"
#pr_out = "./../../results/Phthalates/RDKIT_desc/histDesc/"


d_desc = read.csv(p_desc, sep = "\t")
rownames(d_desc) = d_desc[,1]
d_desc = d_desc[,-1]

d_desc_cleaned = read.csv(p_desc_cleaned, sep = ",")
rownames(d_desc_cleaned) = d_desc_cleaned[,1]
d_desc_cleaned = d_desc_cleaned[,-1]



for(desc in colnames(d_desc_cleaned)){
  
  d_plot = cbind(rownames(d_desc_cleaned), rownames(d_desc_cleaned))
  d_plot = cbind(d_plot, d_desc_cleaned[,desc])
  colnames(d_plot) = c("ID", "ID2", "desc")
  
  d_plot = as.data.frame(d_plot)
  d_plot$desc = as.double(as.character(d_plot$desc))
  
  
  ggplot(d_plot, aes(x=desc)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.20, color="blue", fill="blue")+
    geom_density(alpha=0.3, fill="blue")+
    theme(text = element_text(size=19))+
    scale_color_manual(values=c("#cde2ff")) + 
    labs(title=desc, x=desc, y = "Density")
  
  ggsave(paste(pr_out, "hist_", desc, ".png", sep = ""),  width = 6, height = 7, dpi = 300, bg="transparent")
  
}


l_desc_extra = c("MolWt", "NumAllatoms", "MolLogP")

for(desc in l_desc_extra){
  
  d_plot = cbind(rownames(d_desc), rownames(d_desc))
  if(desc %in% colnames(d_desc) == TRUE){
    d_plot = cbind(d_plot, d_desc[,desc])
    colnames(d_plot) = c("ID", "ID2", "desc")
    
    d_plot = as.data.frame(d_plot)
    d_plot$desc = as.double(as.character(d_plot$desc))
    
    
    ggplot(d_plot, aes(x=desc)) +
      geom_histogram(aes(y=..density..), position="identity", alpha=0.20, color="blue", fill="blue")+
      geom_density(alpha=0.3, fill="blue")+
      theme(text = element_text(size=19))+
      scale_color_manual(values=c("#cde2ff")) + 
      labs(title=desc, x=desc, y = "Density")
    
    ggsave(paste(pr_out, "hist_", desc, ".png", sep = ""),  width = 6, height = 7, dpi = 300, bg="transparent")
    
  }
}

