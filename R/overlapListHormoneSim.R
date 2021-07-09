#!/usr/bin/env Rscript
library(ggpubr)
library(cowplot)


make_hist = function(w_data, pr_out){
  
  # 1. Create the histogram plot
  phist <- gghistogram(
    w_data, x = "Similarity", 
    add = "mean", rug = TRUE,
    fill = "Dataset", palette = c("#00AFBB", "#E7B800")
  )
  
  # 2. Create the density plot with y-axis on the right
  # Remove x axis elements
  pdensity <- ggdensity(
    w_data, x = "Similarity", 
    color= "Dataset", palette = c("#00AFBB", "#E7B800"),
    alpha = 0
  ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")  +
    theme_half_open(11, rel_small = 1) +
    rremove("x.axis")+
    rremove("xlab") +
    rremove("x.text") +
    rremove("x.ticks") +
    rremove("legend")
  
  # 3. Align the two plots and then overlay them.
  aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
  ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
  
  ggsave(paste(pr_out, ".png", sep = ""))
}




################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_upset = args[1]
p_hormone_similarity = args[2]
pr_out = args[3]


#p_upset = "/mnt/c/Users/AlexandreBorrel/research/SSI/breast_carcinogen/results/OverlapList/E2-up-P4-up/upset_matrix"
#p_hormone_similarity = "/mnt/c/Users/AlexandreBorrel/research/SSI/breast_carcinogen/results/similarityHormone/matrix_MACCS-Tanimoto.csv"
#pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/breast_carcinogen/results/OverlapList/E2-up-P4-up/"


d_upset = read.csv(p_upset, sep = "\t", row.names = 1)
d_dataset1 = d_upset[which(dim(d_upset)[2] == rowSums(d_upset)),]
d_dataset2 = d_upset[which(dim(d_upset)[2] != rowSums(d_upset)),]


d_hormone = read.csv(p_hormone_similarity, sep = "\t", row.names = 1, header = FALSE)
colnames(d_hormone) = d_hormone[1,]

l_hormone = colnames(d_hormone)

name_dataset1 = paste("Inter", paste(colnames(d_upset), collapse = "-"), sep = " ")
name_dataset2 = paste("Union", paste(colnames(d_upset), collapse = "-"), sep = " ")

for(hormone in l_hormone){
  w_data1 = cbind(na.omit(d_hormone[rownames(d_dataset1) ,hormone]), name_dataset1)
  w_data2 = cbind(na.omit(d_hormone[rownames(d_dataset2) ,hormone]), name_dataset2)
  w_data = rbind(w_data1, w_data2)
  colnames(w_data) = c("Similarity", "Dataset")
  w_data = as.data.frame(w_data)
  w_data$Similarity = as.double(w_data$Similarity)
  make_hist(w_data, paste(pr_out, hormone, sep = ""))
}


