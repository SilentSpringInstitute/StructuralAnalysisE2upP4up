#!/usr/bin/env Rscript
library(ggpubr)
library(cowplot)
library(dplyr)


make_hist = function(w_data, pr_out){
  
  # 1. Create the histogram plot
  phist <- gghistogram(
    w_data, x = "Similarity", 
    add = "mean", rug = TRUE,
  )
  
  # 2. Create the density plot with y-axis on the right
  # Remove x axis elements
  pdensity <- ggdensity(
    w_data, x = "Similarity", 
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


make_hist_eff = function(w_data, pr_out){
  
  # 1. Create the histogram plot
  phist <- gghistogram(
    w_data, x = "Similarity", 
    add = "mean", rug = TRUE,
    fill = "Efficacy.potency", palette = "jco", binwidth=0.1
  )
  
  # 2. Create the density plot with y-axis on the right
  # Remove x axis elements
  pdensity <- ggdensity(
    w_data, x = "Similarity", 
    color= "Efficacy.potency", palette = "jco",
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
p_cor = args[1]
dataset = args[2]

#p_cor = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/corrSimHormActiveClass/E2up_cor_50-28-2"
#dataset = "E2up"

d_cor = read.csv(p_cor, sep = "\t")

w_data = as.data.frame(d_cor)
w_data$Similarity = as.double(w_data$Similarity)

if("Efficacy.potency" %in% colnames(w_data)){
  make_hist_eff(w_data, p_cor)
  
  # Look similarity by class
  dt = table(d_cor$Efficacy.potency)
  dt = cbind(dt, aggregate(w_data[,2], list(w_data$Efficacy.potency), mean))
  
  # Boxplot with similarity
  sum_data = group_by(w_data, Efficacy.potency) %>%
    summarise(
      count = n(),
      mean = mean(Similarity, na.rm = TRUE),
      sd = sd(Similarity, na.rm = TRUE)
    )
  
  # Add all chemicals
  All = c("All", dim(w_data)[1], mean(w_data$Similarity), sd(w_data$Similarity))
  sum_data = rbind(sum_data, All)
  
  write.csv(sum_data, paste(p_cor, "_sum.csv", sep = ""))
  
  ### add significant on boxplot
  # Visualize: Specify the comparisons you want
  #compare_means(Similarity ~ Efficacy.potency,  data = w_data, method = "t.test")
  my_comparisons <- list( c("medium", "lower"), c("medium", "higher"), c("lower", "higher") )
  ggboxplot(w_data, x = "Efficacy.potency", y = "Similarity",
            color = "Efficacy.potency", palette = "jco", ylab = "Similarity", xlab = "Efficacy/potency", order = c("lower", "medium", "higher"))+ 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")#+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 0.8, method = "anova")     # Add global p-value
  
  ggsave(paste(p_cor, "_boxplot.png", sep = ""))
}else{
  make_hist(w_data, p_cor)
  
  # Boxplot with similarity
  sum_data = group_by(w_data) %>%
    summarise(
      count = n(),
      mean = mean(Similarity, na.rm = TRUE),
      sd = sd(Similarity, na.rm = TRUE)
    )
  
  write.csv(sum_data, paste(p_cor, "_sum.csv", sep = ""))
}


