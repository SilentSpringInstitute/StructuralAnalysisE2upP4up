#!/usr/bin/env Rscript


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

cardMatrixCor = function(matrixIN, name_file, nb_color){
  
  #print(matrixIN)
  
  nb_col = dim(matrixIN)[2]
  nb_line = dim(matrixIN)[1]
  
  y_names = seq(0,nb_col,1)
  x_names =  seq(0,nb_line,1)
  
  
  list_color = generateColor(nb_color)
  
  #postscript (file = paste (name_file, ".ps", sep = ""), width = 0, height = 0)
  
  #par( mar=c(25,25,1.5,1.5))
  #image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", col = list_color)
  #grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
  #box()
  # place les petites barres
  #axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
  #axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)
  
  # place les positions en fonction du cut
  # ecart1 = 1/(nb_line-1)
  # ecart2 = 1/(nb_col-1)
  # list_L1 = generateLegend (nb_line,1)
  # list_L2 = generateLegend (nb_col,1)
  
  # place les legendes
  # posX = generatePosition(list_L1, ecart1)
  # posY = generatePosition(list_L2, ecart2)
  # axis(1,seq(0,1,(1/(nb_line-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
  # axis(2,seq(0,1,(1/(nb_col-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
  # dev.off()
  
  
  
  
  # PNG file
  
  # list descriptor
  
  # fix breaks
  x = 0
  while ((length (list_color)+1) != length (seq(-1,1, 2/(nb_color+x)))){
    x = x + 1
  }
  
  
  png (file = paste (name_file, ".png", sep = ""), width = 2200, height = 2200)
  
  par( mar=c(20,20,1.5,1.5))
  image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", col = list_color, breaks = seq(-1,1, 2/(nb_color+x)))
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
  axis(1,seq(0,1,(1/(nb_line-1))),colnames (matrixIN), cex.axis = 1, las = 2)
  axis(2,seq(0,1,(1/(nb_col-1))),colnames (matrixIN), cex.axis = 1, las = 2)
  dev.off()
  
  
}

cardMatrixCorpdf = function(matrixIN, nb_color, title){
  
  #print(matrixIN)
  
  nb_col = dim(matrixIN)[2]
  nb_line = dim(matrixIN)[1]
  
  y_names = seq(0,nb_col,1)
  x_names =  seq(0,nb_line,1)
  
  
  list_color = generateColor(nb_color)
  
  x = 0
  while ((length (list_color)+1) != length (seq(-1,1, 2/(nb_color+x)))){
    x = x + 1
  }
  
  par( mar=c(20,20,5,1.5))
  image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", col = list_color, breaks = seq(-1,1, 2/(nb_color+x)), main = title)
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
  axis(1,seq(0,1,(1/(nb_line-1))),colnames (matrixIN), cex.axis = 1, las = 2)
  axis(2,seq(0,1,(1/(nb_col-1))),colnames (matrixIN), cex.axis = 1, las = 2)
  
}

cardMatrix = function(matrixIN, name_file, nb_color){
  
  nb_col = dim(matrixIN)[2]
  nb_line = dim(matrixIN)[1]
  
  list_color = generateColor(nb_color)
  
  png (file = paste (name_file, ".png", sep = ""), width = 2500, height = 1600)
  
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
pFP = args[1]

#pFP = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/similarityInterHormone/matrix_MACCS-Dice.csv"

dFP = read.csv(pFP, sep = "\t", header = TRUE)
rownames(dFP) = dFP[,1]
dFP = dFP[,-1]
colnames(dFP) = rownames(dFP)


cardMatrix(dFP, pFP, 5)

