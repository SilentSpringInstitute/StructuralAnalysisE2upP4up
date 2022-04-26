#!/usr/bin/env Rscript
require(kohonen)
library(ggplot2)



#################
# Color reverse #
#################

colors <- function(n, alpha = 1) {
  rev(heat.colors(n, alpha))
}

coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}

#######
# SOM #
#######

generateSOM = function(ddesc, xdim, ydim, pr_out){
  
  ddesc = as.matrix(scale(ddesc))
  som_grid <- somgrid(xdim=xdim, ydim=ydim, topo="hexagonal")
  som_model <- som(ddesc, 
                   grid=som_grid, 
                   rlen=100, 
                   alpha=c(0.05,0.01), 
                   keep.data = TRUE) 
  
  save(som_model, file = paste(pr_out, "SOM_model.RData", sep = ""))
  
  
  # draw SOM
  svg(paste(pr_out, "SOM_model_count.svg", sep = ""))
  plot(som_model, type = "count", palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  # draw in png
  png(paste(pr_out, "SOM_model_count.png", sep = ""))
  plot(som_model, type = "count", palette.name=coolBlueHotRed, main = "", dpi=300, height = 500, width = 500, bg = "transparent")
  dev.off()
  
  # cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(som_model$data[[1]])
  dclust = dclust[rownames(ddesc)]# take only active chemical
  dclust = cbind(names(dclust), dclust)
  colnames(dclust) = c("names", "cluster")
  write.table(dclust, paste(pr_out, "SOM_Clusters", sep = ""), sep = ",", row.names = TRUE)
  
  return(som_model)
}



applySOM = function(som_model, d_AC50, pr_out){
  
  #write cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(som_model$data[[1]])
  
  xdim = som_model$grid$xdim
  ydim = som_model$grid$ydim
  

  lAct = rep(0, xdim*ydim)
  names(lAct) = seq(1, xdim*ydim)
  
  # remove inactive
  if(d_AC50 != "0"){
    d_AC50 = na.omit(d_AC50)
    dclust = dclust[rownames(d_AC50)]# take only active chemical
    dclust = cbind(names(dclust), dclust)
    ltable = table(dclust[,2])
    lAct[names(ltable)] = ltable
    colnames(dclust) = c("CASRN", "Cluster")
  }else{
    return() # to not run because no classification
  }
  
  
  ltabinit = table(som_model$unit.classif)
  linitial = rep(0, xdim*ydim)
  names(linitial) = seq(1, xdim*ydim)
  linitial[names(ltabinit)] = ltabinit
  
  lprob = lAct / linitial
  
  
  write.csv(dclust, paste(pr_out, "SOM_Clusters_act", sep = ""))
  
  # count of active
  svg(paste(pr_out, "SOM_count_act.svg", sep = ""))
  plot(som_model, type = "property", property=lAct, palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  # count of active calibrate
  lAct[which(lAct == max(lAct))] = max(table(som_model$unit.classif))# have to calibrate based on the max of the original SOM
  svg(paste(pr_out, "SOM_count_act_calibrate.svg", sep = ""))
  plot(som_model, type = "property", property=lAct, palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  # plot with proba
  svg(paste(pr_out, "SOM_prob_act.svg", sep = ""))
  plot(som_model, type = "property", property=lprob, palette.name=coolBlueHotRed, main = "Prob active", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
    
  write.csv(lprob, paste(pr_out, "SOM_Clusters_prob", sep = ""))
  
  # enrichment with a fisher score
  nbInact = length(which(is.na(dAC50)))
  nbAct = dim(dAC50)[1] - nbInact
  lpval = NULL
  nbclust = ydim*xdim
  for(clust in seq(1,nbclust)){
    dclust = dAC50[which(som_model$unit.classif == clust),]
    ninactClust = length(which(is.na(dclust[,2])))
    nactClust = dim(dclust)[1] - ninactClust
    FisherMat = matrix(c(nbAct, nbInact, nactClust, ninactClust),
                        nrow = 2,
                        dimnames = list(Pop = c("Act", "Inact"),
                                        Clust = c("Act", "Inact")))
    p = fisher.test(FisherMat, alternative = "greater")
    pval = p$p.value
    lpval = append(lpval, p$p.value)
  }
  svg(paste(pr_out, "SOM_enrich_act.svg", sep = ""))
  plot(som_model, type = "property", property=log10(lpval), palette.name=coolBlueHotRed, main = "Enrichment log(pvalue(Fisher test))", dpi=300, height = 20, width = 20)
  dev.off()
}


analyzeSOMSize = function(ddesc, l_sizes, pr_out){
  
  
  ddesc = as.matrix(scale(ddesc))
  
  M_chem_cluster = NULL
  l_emptyCluster = NULL
  
  for(size in l_sizes){
    M_chem_cluster_temp = NULL
    l_emptyCluster_temp = NULL
    for(i in seq(1,5)){
      som_grid <- somgrid(xdim=size, ydim=size, topo="hexagonal")
      som_model <- som(ddesc, 
                      grid=som_grid, 
                      rlen=100, 
                      alpha=c(0.05,0.01), 
                      keep.data = TRUE) 
      
      ltabinit = table(som_model$unit.classif)
      M_chem_cluster_temp = append(M_chem_cluster_temp, mean(ltabinit))
      l_emptyCluster_temp = append(l_emptyCluster_temp, size * size - length(ltabinit))
    }
    M_chem_cluster = append(M_chem_cluster, mean(M_chem_cluster_temp))
    l_emptyCluster = append(l_emptyCluster, size * size - length(l_emptyCluster_temp))
  }
  
  d_out = cbind(M_chem_cluster, l_emptyCluster)
  d_out = cbind(d_out, l_sizes)
  rownames(d_out) = l_sizes
  
  write.table(d_out, paste(pr_out, "size_analysis.csv", sep = ""), sep = ",")
  
  d_out = as.data.frame(d_out)
  
  scaleFactor <- max(d_out$M_chem_cluster) / max(d_out$l_emptyCluster)
  
  ggplot(d_out, aes(x=l_sizes)) +
    geom_smooth( aes(y=M_chem_cluster), col="blue") + 
    geom_smooth( aes(y=l_emptyCluster * scaleFactor), col="red") + 
    scale_y_continuous(
      # Features of the first axis
      name = "Av. chemical by cluster",
      # Add a second axis and specify its features
      sec.axis = sec_axis( ~./scaleFactor, name="Nb. empty cluster")
    )+
    theme(
      axis.title.y.left=element_text(color="blue"),
      axis.text.y.left=element_text(color="blue"),
      axis.title.y.right=element_text(color="red"),
      axis.text.y.right=element_text(color="red")
    )
  ggsave(paste(pr_out, "SOM_size.png", sep = ""), dpi=300, height = 8, width = 8)
  
  
}

################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc = args[1]
p_AC50 = args[2]
pr_out = args[3]
size = as.integer(args[4])


#p_desc = "../../results/Analysis_H295R/rdkit/Cleaned_Data/desc1D2D_cleaned.csv"
#p_AC50 = "0"
#pr_out = "../../results/Analysis_H295R/rdkit/SOM/"
#size = 12

#p_desc = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/clusterMC/MC/rdkit/Cleaned_Data/desc1D2D_cleaned.csv"
#pr_out = "c://Users/aborr/research/Silent_Spring/breast_carcinogen/results/clusterMC/MC/rdkit/"


# open files
ddesc = read.csv(p_desc, sep = ",", header = TRUE)
if(dim(ddesc)[2] == 1){
  ddesc = read.csv(p_desc, sep = "\t", header = TRUE)
}
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

if (p_AC50 != "0"){
    dAC50 = read.csv(p_AC50, sep = ",", header = TRUE)
    rownames(dAC50) = dAC50[,1]
    dAC50 = dAC50[,-1]
}else{
    dAC50 = "0"
}

#som_model = generateSOM(ddesc, size, size, pr_out)

# define model
p_model = paste(pr_out , "SOM_model.RData", sep = "")
if(!file.exists(p_model)){
  # if size = 0 => need to analyse size of SOM
  if(size == 0){
    l_sizes = c(3,4,5,6,7,8,9,10,11,12,13,14,15)
    analyzeSOMSize (ddesc, l_sizes, pr_out)
  }else{
    som_model = generateSOM(ddesc, size, size, pr_out)
    applySOM(som_model, dAC50, pr_out)
  }
}else{
  load(p_model)
  applySOM(som_model, dAC50, pr_out)
}
