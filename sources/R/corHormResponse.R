#!/usr/bin/env Rscript
library(ggpubr)
library(cowplot)
library(dplyr)
library(grid)

make_hist = function(d_cor, pr_out){
  

  png(paste(pr_out, "hist_horm.png", sep = ""), 500, 1000)
  pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 2)))
  
  # A helper function to define a region on the layout
  define_region <- function(row, col){
    viewport(layout.pos.row = row, layout.pos.col = col)
  } 
    
  i_col = 1
  i_row = 1
  
  #par(mfrow=c(5,2))
  
  
  for (hormone in colnames(d_cor)){
    d_cor[, c(hormone)] = as.double(as.character(d_cor[, c(hormone)]))
    
    # 1. Create the histogram plot
    phist <- gghistogram(
      d_cor, x = hormone, 
      add = "mean", rug = TRUE,
      fill = "lightgrey", palette = "jco", binwidth=0.5,
      title = hormone
    )
    
    
    # 2. Create the density plot with y-axis on the right
    # Remove x axis elements
    pdensity <- ggdensity(
      d_cor, x = hormone, 
      palette = "jco",
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
    print(ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]]), vp = define_region(row = i_row, col = i_col))
    
    if(i_col == 2){
      i_col = 1
      i_row = i_row + 1
    }else{
      i_col = 2
    }
    
  }
  dev.off()
  
}



make_corplot = function(d_cor, pr_out){
  
  png(paste(pr_out, "cor_horm.png", sep = ""), 1000, 2000)
  par(mfrow = c(10,5))
  
  
  i_col = 1
  i_row = 1

  i = 1
  imax = length(colnames(d_cor))
  while(i <= imax){
    j = i + 1
    while(j <= imax){
     cor_val = cor(d_cor[,i], d_cor[,j])
     plot(d_cor[,i], d_cor[,j], pch = 19, xlab = colnames(d_cor)[i], ylab = colnames(d_cor)[j])
      text(paste("r=", round(cor_val, 2), sep = ""), x=2, y=4,size = 4)
     
      if(i_col == 5){
        i_col = 1
        i_row = i_row + 1
      }else{
        i_col = i_col + 1
      }
      j = j + 1
    }
    i = i + 1
  } 
  
  dev.off()
  
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_cor = args[1]
pr_out = args[2]

#p_cor = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/steroidogenesis/CorHorm/horm_resp.csv"
#pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/steroidogenesis/CorHorm/"

d_cor = read.csv(p_cor, sep = "\t", row.names = 1)
d_cor = as.data.frame(d_cor)
for (hormone in colnames(d_cor)){
  d_cor[, c(hormone)] = as.double(as.character(d_cor[, c(hormone)]))
}

# histogram
make_hist(d_cor, pr_out)


# correlation plot
make_corplot(d_cor, pr_out)


