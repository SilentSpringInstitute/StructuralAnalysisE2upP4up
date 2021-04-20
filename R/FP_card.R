#!/usr/bin/env Rscript
library(ggplot2)
library(factoextra)


generateColor= function(nb_color){
  
  if (nb_color == 2){
    return (c("#FF0000", "#66FF77"))
  }
  else if(nb_color == 3){
    return (c("#66FF77","#FF0000", "#338800"))
  }
  else if(nb_color == 4){
    #return (c("#FF0000","#FFCCCC","#338800", "#66FF77"))
    return (c("#66FF77","#66FF77","#FF0000", "#66FF77"))
  }
  else if(nb_color == 5){
    return (c("#FFFFFF", "#FFE3E3", "#FFCCCC", "#FF0000","#8B2323"))
  }
  else if (nb_color == 6){
    return (c("#003300","#338800","#66FF77","#FFFFFF","#FFCCCC", "#FF0000","#8B2323"))
  }
  else if (nb_color == 10){
    return (c("#6600CC","#0000CC","#3366FF","#33CCFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFF66","#FFCC33","#FF6633","#CC0000"))
    
  }else if(nb_color == "red"){
    return(c("#FFFAFA", "#F4C2C2", "#FF6961", "#FF5C5C", "#FF1C00", "#FF0800", "#FF0000", "#CD5C5C", "#E34234", "#D73B3E", "#CE1620", "#CC0000", "#B22222", "#B31B1B", "#A40000", "#800000", "#701C1C", "#3C1414", "#321414"))
    
  }
}

generateLegend = function (length_tot, cut){
  
  list_out = seq(0,length_tot, cut)
  if (list_out[length(list_out)] != length_tot){
    list_out = append(list_out, length_tot)
  }
  return (list_out)
}

generatePosition = function (list_legend, value_ecart){
  
  list_out = c()
  for (element in list_legend){
    list_out = append(list_out,element * value_ecart)
  }
  return (list_out)
}


cardMatrix = function(matrixIN, name_file, nb_color){
  
  nb_col = dim(matrixIN)[2]
  nb_line = dim(matrixIN)[1]
  
  list_color = generateColor(nb_color)
  
  png (file = paste (name_file, ".png", sep = ""), width = 1600, height = 2500)
  
  matrixText = matrixIN
  matrixIN[matrixIN > 15] = 15
  
  par( mar=c(25,25,1.5,1.5))
  image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", col = list_color)
  grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
  box()
  
  nbcol = dim(matrixIN)[2]
  nbline = dim(matrixIN)[1]
  for (i in seq(0,nbline-1)){
    for (j in seq(0, nbcol-1)){
      text((1/(nbline-1))*i,(1/(nbcol-1))*j, labels = round(matrixText[i+1,j+1], digits = 1), cex = 1.5)
    }
  }
  
  # place les petites barres
  axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
  axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)
  
  # place les positions en fonction du cut
  ecart1 = 1/(nb_line-1)
  ecart2 = 1/(nb_col-1)
  list_L1 = generateLegend (nb_line,1)
  list_L2 = generateLegend (nb_col,1)
  
  # place les legendes
  posX = generatePosition(list_L1, ecart1)
  posY = generatePosition(list_L2, ecart2)
  axis(1,seq(0,1,(1/(nb_line-1))),rownames (matrixIN), cex.axis = 1.75, las = 2)
  axis(2,seq(0,1,(1/(nb_col-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
  dev.off()
  
  svg (file = paste (name_file, ".svg", sep = ""), 20, 20)
  par( mar=c(5,17,0.5,0.5))
  
  image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", col = list_color)
  grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
  box()
  # place les petites barres
  axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
  axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)
  
  # place les positions en fonction du cut
  ecart1 = 1/(nb_line-1)
  ecart2 = 1/(nb_col-1)
  list_L1 = generateLegend (nb_line,1)
  list_L2 = generateLegend (nb_col,1)
  
  # place les legendes
  posX = generatePosition(list_L1, ecart1)
  posY = generatePosition(list_L2, ecart2)
  axis(1,seq(0,1,(1/(nb_line-1))),rownames (matrixIN), cex.axis = 1.75, las = 2)
  axis(2,seq(0,1,(1/(nb_col-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
  dev.off()
}




################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_matrixIn = args[1]


#p_matrixIn = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/steroidogenesis/card_FP/single_hit_matrix.csv"

d_in = read.csv(p_matrixIn, sep = "\t", header = TRUE)
rownames(d_in)= d_in[,1]
d_in = d_in[,-1]
l_MC = d_in$MC
drops = c("MC")
d_in = d_in[ , !(names(d_in) %in% drops)]
d_in = as.matrix(d_in)


# transform matrix in binary matrix
d_bin = d_in
d_bin[which(d_bin < 0)] = -1
d_bin[which(d_bin > 0)] = 1


# Extract MC
d_MC = d_bin[which(l_MC == 1),]
cardMatrix(t(d_MC), paste(substr(p_matrixIn, 1, nchar(p_matrixIn)-4), "_MC", sep = ""), 3)


# tanimoto in d_MC
dist_MC = dist(d_MC, method = "binary")
dist_all = dist(d_bin, method = "binary")

M_all = mean(dist_all)
M_MC = mean(dist_MC)
sd_all = sd(dist_all)
sd_MC = sd(dist_MC)
p_val = t.test(dist_all, dist_MC)$p.value
l_out = c(M_all, sd_all, M_MC, sd_MC, p_val)
names(l_out) = c("M_all", "sd_all", "M_MC", "sd_MC", "p_val")

write.csv(l_out, paste(substr(p_matrixIn, 1, nchar(p_matrixIn)-4), ".sum", sep = ""))
