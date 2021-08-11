#!/usr/bin/env Rscript
library(ggpubr)
library(cowplot)
library(tidyr)

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


#p_upset = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/OverlapList/E2-up-P4-up-H295R/upset_matrix"
#p_hormone_similarity = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/similarityHormone/matrix_MACCS-Tanimoto.csv"
#pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/OverlapList/E2-up-P4-up-H295R/"


d_upset = read.csv(p_upset, sep = "\t", row.names = 1)
d_dataset1 = d_upset[which(dim(d_upset)[2] == rowSums(d_upset)),]
d_dataset2 = d_upset[which(dim(d_upset)[2] != rowSums(d_upset)),]


d_hormone = read.csv(p_hormone_similarity, sep = "\t", row.names = 1, header = FALSE)
colnames(d_hormone) = d_hormone[1,]

l_hormone = colnames(d_hormone)

#name_dataset1 = paste("Inter", paste(colnames(d_upset), collapse = "-"), sep = " ")
#name_dataset2 = paste("Union", paste(colnames(d_upset), collapse = "-"), sep = " ")

for(hormone in l_hormone){
  #w_data1 = cbind(na.omit(d_hormone[rownames(d_dataset1) ,hormone]), name_dataset1)
  #w_data2 = cbind(na.omit(d_hormone[rownames(d_dataset2) ,hormone]), name_dataset2)
  #w_data = rbind(w_data1, w_data2)
  #colnames(w_data) = c("Similarity", "Dataset")
  #w_data = as.data.frame(w_data)
  #w_data$Similarity = as.double(w_data$Similarity)
  #make_hist(w_data, paste(pr_out, hormone, sep = ""))
  
  d_forplot = NULL
  for (i_dataset in colnames(d_upset)){
    d_temp = d_upset[which(d_upset[,i_dataset] == 1),]
    Dataset = rep(i_dataset, dim(d_temp)[1])
    Similarity = d_hormone[rownames(d_temp), hormone]
    d_temp_out = cbind(Dataset, Similarity)
    rownames(d_temp_out) = rownames(d_temp)
    print(dim(d_temp))
    
    d_forplot = rbind(d_forplot, d_temp_out)
  }
  d_forplot = d_forplot[,c("Similarity", "Dataset")]
  d_forplot = as.data.frame(d_forplot)
  d_forplot$Similarity = as.double(as.character(d_forplot$Similarity))
  
  ### add significant on boxplot
  # Visualize: Specify the comparisons you want
  my_comparisons <- as.list(as.data.frame(combn(unique(d_forplot$Dataset), 2)))
  print(my_comparisons)
  ggviolin(d_forplot, x = "Dataset", y = "Similarity",
            color = "Dataset", palette = "jco", ylab = "Similarity", xlab = "Dataset", combine = TRUE, add = "boxplot")+ 
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif")+
    stat_summary(fun.data = function(x) data.frame(y=0.9, label = paste("Mean=",round(mean(x),2))), geom="text")+
    stat_summary(fun.data = function(x) data.frame(y=0.85, label = paste("Median=",round(median(x),2))), geom="text")
    #+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 0.8, method = "anova")     # Add global p-value
  
  ggsave(paste(pr_out, hormone, "_boxplot.png" ,sep = ""))
}


