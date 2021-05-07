#!/usr/bin/env Rscript
library(ggplot2)
library(gplots)
library(graphics)
library(corrplot)


library (FactoMineR)

factorAFC = function (xplot, yplot, xplotdata, yplotdata ){
  # max window
  print (yplotdata)
  
  max_window_x = max (abs (xplotdata))
  max_window_y = max (abs (yplotdata))
  
  # max coordinate
  max_x = max (abs (xplot))
  max_y = max (abs (yplot))
  
  
  factor = 1
  while ( max_x < max_window_x  && max_y < max_window_y ){
    factor = factor + 0.1
    max_x = max_x * factor
    max_y = max_y * factor
    
  }
  
  
  if (factor == 1 ){
    while ( max_x > max_window_x || max_y > max_window_y ){
      
      #print (max_x)
      #print (max_y)
      #print (max_window_x)
      #print (max_window_y)
      #print ("********")
      
      factor = factor - 0.2
      max_x = max_x * factor
      max_y = max_y * factor
      
    }
  }
  return (factor)	
}



AFC = function (d, path_file){
  
  # genere afc coords
  r = CA (d, graph = FALSE)
  
  # plot
  xplot = c(r$row$coord[,1], r$col$coord[,1])
  yplot = c(r$row$coord[,2], r$col$coord[,2])
  l_pch = c(rep (20,length (r$row$coord[,1])), rep (18,length (r$col$coord[,2])))
  l_col_all = c(rep ("#000000", length (r$row$coord[,1])), "#03fc73")
  l_col = "#03fc73"
  
  lim_x = max (abs (xplot))
  lim_y = max (abs (yplot))
  
  
  svg (file = paste (path_file, "_AFC_text.svg", sep = ""), bg = "transparent" ,100, 100)
  par(mar=c(5,5,2,2))
  
  plot (xplot, yplot, xlab = paste("DIM 1 : ", round(r$eig[1,2],1), "%", sep = ""), ylab = paste("DIM 2 : ", round(r$eig[2,2],1), "%", sep = ""), cex.lab = 2.2, cex.axis = 1.8, xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y), pch = l_pch, col = l_col_all, cex = 2.5)
  text (r$col$coord[,1], r$col$coord[,2], label = names(r$col$coord[,1]), col = l_col, cex = 2, pos = 4)
  text (r$row$coord[,1], r$row$coord[,2], col = "black", label = names (r$row$coord[,1]), cex = 2.5, pos = 4)
  abline(h=0,v=0, lwd = 2)
  
  dev.off()
  
  svg (file = paste (path_file, "_AFC_point.svg", sep = ""), bg = "transparent", 100, 100)
  par(mar=c(5,5,2,2))
  plot (xplot, yplot, type = "n", xlab = paste("DIM 1 : ", round(r$eig[1,2],1), "%", sep = ""), ylab = paste("DIM 2 : ", round(r$eig[2,2],1), "%", sep = ""), cex.lab = 2.2, cex.axis = 1.8, xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y))
  points (r$col$coord[,1], r$col$coord[,2], col = l_col, cex = 4, pch = 16)
  text (r$row$coord[,1], r$row$coord[,2], col = "black", label = names (r$row$coord[,1]), cex = 2.5, pos = 4)
  abline(h=0,v=0, lwd = 2)
  
  dev.off()
  
  
  
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_count = args[1]

#p_count = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/ToxPrintComparisonCount/MC-all/count.csv"

housetasks <- read.delim(p_count, row.names = 1)


# 1. convert the data as a table
dt <- as.table(as.matrix(housetasks))

# reduce
dt = dt[-which(rowSums(dt)<=(5*dim(dt)[2])),]



# 2. Graph - ballon
png(paste(substr(p_count, 1, nchar(p_count)-4), "_ballon.png",sep = ""), height = 2500, width = 1000)
balloonplot(t(dt), main ="ToxPrint", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)

dev.off()



# 2. Graph - mosaic
png(paste(substr(p_count, 1, nchar(p_count)-4), "_mosaic.png",sep = ""), height = 2500, width = 4000 , res = 300)
mosaicplot(dt, shade = TRUE, las=2,
           main = "ToxPrint")

dev.off()


# apply X2
chisq <- chisq.test(housetasks)
pval = chisq$p.value
res = chisq$residuals[-which(is.na(chisq$residuals[,1])),]

png(paste(substr(p_count, 1, nchar(p_count)-4), "_corplot.png",sep = ""), height = 6000, width = 5000)
corrplot(res, is.cor = FALSE)
dev.off()

# need at least # dimension to do the plot 2D
if(dim(dt)[2] > 2){
  AFC(dt, substr(p_count, 1, nchar(p_count)-4))
}

