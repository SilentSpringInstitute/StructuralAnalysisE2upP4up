#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)



dendogramCircle = function(ddes, daff, col_desc, prout){
  
  #calibrate affinity for color
  daff = as.data.frame(daff)
  
  
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  legend.title = "Weight"
  
  
  # all chemicals
  pfilout = paste(prout, "HClust_dendo_", col_desc, ".png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=Aff, label=label, angle=angle), hjust=-0.5, size=3) +
    geom_tippoint(aes(color=Aff), alpha=0.75, size=0.5)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    guides(col=guide_legend(col_desc))+
    geom_treescale(x = 1, y = 1, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  
  ggsave(pfilout, dpi=300, height = 10, width = 14)
  
  
}

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc = args[1]
p_pred = args[2]
pr_out = args[3]

p_desc = "../../results/Cleaned_Data/desc1D2D_cleaned.csv"
p_pred = "../../results/DESC/desc_OPERA.csv"
pr_out = "../../results/HClust/"


# open files
ddesc = read.csv(p_desc, sep = ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

daff = read.csv(p_pred, sep = ",")
rownames(daff) = daff[,1]
daff = as.data.frame(daff)
daff = daff[rownames(ddesc),]


l_aff_to_plot = colnames(daff)
for(aff_to_plot in l_aff_to_plot){
  
  if(!is.integer0(grep("_pred", aff_to_plot, fixed = TRUE))){
    dtemp = cbind(rownames(daff), daff[,aff_to_plot])
    rownames(dtemp) = rownames(daff)
    colnames(dtemp) = c("ID", "Aff")
    dtemp = as.data.frame(dtemp)
    
    dtemp$Aff = as.double(as.character(dtemp$Aff))
    
    dendogramCircle(ddesc, dtemp, aff_to_plot, pr_out)
      
  }
  
}

#dendogramCircle(ddesc, daff, pr_out)
