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


dendogramCircleClass = function(ddes, daff, col_desc, prout){
  
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
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    guides(col=guide_legend(col_desc))+
    geom_treescale(x = 1, y = 1, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  
  ggsave(pfilout, dpi=300, height = 10, width = 18)
  
  
  
  pfilout = paste(prout, "HClust_dendo_group_", col_desc, ".png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=Aff, label=Aff, angle=angle), vjust=1,hjust=0,  size=3) +
    geom_tippoint(aes(color=Aff), alpha=0.75, size=0.5)+
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    guides( color = FALSE)+
    geom_treescale(x = 1, y = 1, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  
  ggsave(pfilout, dpi=300, height = 38, width = 45, units="cm")
  
  
  
  
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
p_dataset = args[2]
p_pred = args[3]
pr_out = args[4]



#p_desc =  '../../ILS/results/Phthalates/rdkit/Cleaned_Data/desc1D2D_cleaned.csv'
#p_dataset = "../../ILS/data/Phthalates.csv"
#p_pred = '../../ILS/results/Phthalates/DESC/desc_OPERA.csv'
#pr_out = '../../ILS/results/Phthalates/rdkit/HClustCircular/'


# open files
ddesc = read.csv(p_desc, sep = ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

daff = read.csv(p_pred, sep = ",")
rownames(daff) = daff[,1]
daff = as.data.frame(daff)
daff = daff[rownames(ddesc),]

# dataset
d_dataset = read.csv(p_dataset, sep = "\t")
if(dim(d_dataset)[2] == 1){
  d_dataset = read.csv(p_dataset, sep = ",")
}


if(length(unique(d_dataset[,1] == "--" | d_dataset[,1] == "")) != 1){
  d_dataset = d_dataset[-which(d_dataset[,1] == "--" | d_dataset[,1] == ""),]
}

# remove duplicate
if(length(which(duplicated(d_dataset$CASRN))) != 0){
  d_dataset = d_dataset[-which(duplicated(d_dataset$CASRN)),]
}

# rowname
rownames(d_dataset) = d_dataset[,1]


d_dataset = d_dataset[rownames(ddesc),]
if(exists("Group", d_dataset)){
  dtemp = cbind(rownames(ddesc), d_dataset$Group)
  rownames(dtemp) = rownames(ddesc)
  colnames(dtemp) = c("ID", "Aff")
  dtemp = as.data.frame(dtemp)
  dendogramCircleClass(ddesc, dtemp, "Group", pr_out)
}



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
