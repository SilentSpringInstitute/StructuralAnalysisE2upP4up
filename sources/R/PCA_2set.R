#!/usr/bin/env Rscript
library(ggplot2)
library(factoextra)
library(Toolbox)


generatePCAcoords = function(din){
  
  dinScale = scale(din)
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  
  return(list(data_plot, var_cap, cp))
}



addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc1 = args[1]
pdesc2 = args[2]
prout = args[3]


#pdesc1 = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/QSAR_E2_H295R_nosampling_nosingledosecheck_noborderline/rdkit-OPERA-toxprint_0.9-0/Cleaned_Data/desc_cleaned.csv"
#pdesc2 = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/predMC_E2/desc_global.csv"
#prout = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/predMC_E2/"

ddesc1 = read.csv(pdesc1, sep = "," ,row.names = 1)
ddesc2 = read.csv(pdesc2, sep = "\t", row.names = 1)


ddesc1 = delSDnull(ddesc1)
i_ddesc1 = elimcor_sansY(ddesc1, 0.85)[[1]]
ddesc1 = ddesc1[, i_ddesc1]

l_inter = intersect(colnames(ddesc1), colnames(ddesc2))
ddesc1 = ddesc1[,l_inter]
ddesc2 = ddesc2[,l_inter]

# remove SMILES
#ddesc2 = ddesc2[,-1]

res.pca <- prcomp(ddesc1, scale = TRUE)
var_cap = generatePCAcoords(ddesc1)[[2]]


add_coord = predict(res.pca, newdata = ddesc2)

vcolor = rep(addTrans("#595959", 60), dim(ddesc1)[1])

vcolor_add = rep("#90EE90", dim(ddesc2)[1])
vshape = rep(18, dim(ddesc2)[1])
#vshape[union(which(daff2$Aff == 0), which(is.na(daff2$Aff)))] = 8

png(paste(prout, "PCA_color.png", sep = ""), 1700, 1500)
par(mar=c(8,8,8,8))
plot(rbind(add_coord[,1],res.pca$x[,1]) ,rbind(add_coord[,2],res.pca$x[,2]) , pch=19, col = "white", xlab = paste("PC1: ", round (var_cap[1], 2), "%", sep = ""), ylab = paste("PC2: ", round(var_cap[2], 2), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 2.5)
points(res.pca$x[,1], res.pca$x[,2], pch=19, col = vcolor, cex = 2.5)
#points(res.pca$x[which(vcolor != "darkgray"),1], res.pca$x[which(vcolor != "darkgray"),2], pch=19, col = vcolor[which(vcolor != "darkgray")], cex = 2.5)
points(res.pca$x[which(vcolor != "#5959593C"),1], res.pca$x[which(vcolor != "#5959593C"),2], pch=19, col = vcolor[which(vcolor != "#5959593C")], cex = 2.5)
points(add_coord[,1], add_coord[,2], col = vcolor_add, pch = vshape, cex = 3)
abline(h=0,v=0)
warnings ()
legend(min(res.pca$x[,1]), max(res.pca$x[,2]), legend=c("Train set", "Test set"),
       col=c("#595959", "#90EE90"),pch = 19, cex=2)

dev.off()

#points(data_plot[which(vcolor != "gray90"),1], data_plot[which(vcolor != "gray90"),2], pch=20, col = vcolor[which(vcolor != "gray90")], cex = 4)
#plot(res.pca$x[,1], res.pca$x[,2], xlim = c(-40, 20), ylim = c(-20, 20))

#points(add_coord[,1], add_coord[,2], col = "blue")
#points(res.pca$x[,1], res.pca$x[,2], col = "red")