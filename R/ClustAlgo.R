#!/usr/bin/env Rscript
library(factoextra)
library(ggplot2)
require(cluster)
library(RootsExtremaInflections)
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)
library(Toolbox)


signifDescCluster = function(d_desc, d_opera, d_cluster, pr_out){
  
  # check chemicals are on the same order
  d_opera = d_opera[rownames(d_desc),]
  
  
  l_cluster = seq(1, max(d_cluster$cluster))
  
  d_global = cbind(d_desc, d_opera)
  l_desc = colnames(d_global)
  
  m_out = NULL
  
  for (desc in l_desc){
    l_out = NULL
    for(cluster in l_cluster){
      i_chem = which(d_cluster$cluster == cluster)
      len_i_chem = length(i_chem)
      len_j_chem = dim(d_global)[1]
      v_desc = d_global[,desc]
      v_cluster = v_desc[i_chem]
      v_nocluster = v_desc[-i_chem]
      if(len_i_chem <= 3){
        l_out = append(l_out, "-")
      }else if(len_j_chem - len_i_chem <= 3){
        l_out = append(l_out, "-")
      }else{
        parametric = conditionTtest(v_cluster, v_nocluster)
        if(parametric == 1){
          pval = comparisonTest (v_cluster, v_nocluster, "parametric")
        }else{
          pval = comparisonTest (v_cluster, v_nocluster, "no-parametric")
        }
        if (is.na(pval) == TRUE){
          signif = "-"
        }else{
          signif = signifPvalue(pval)  
        }
        
        
        l_out = append(l_out, signif)
        #print("====")
        #print(length(v_desc))
        #print(length(v_cluster))
        #print(length(v_nocluster))
      }
    }
    m_out = rbind(m_out, l_out)
  }
  
  
  N_CHEM = NULL
  for(cluster in l_cluster){
    i_chem = which(d_cluster$cluster == cluster)
    N_CHEM = append(N_CHEM, length(i_chem))  
  }
  
  rownames(m_out) = l_desc
  colnames(m_out) = l_cluster
  
  m_out = rbind(m_out, N_CHEM)
  
  write.csv(m_out, file=paste(pr_out, "clusters_signif_descriptors.csv", sep = ""))
  
}



optimalCluters = function (din, prout, metcluster, metOptNB, metagregation){
  
  
  # scale data in input
  din = scale (din)
  kmax = dim(din)[1]-1
  if (metcluster == "hclust"){
    if(metOptNB == "gap_stat"){
      p = fviz_nbclust(din, kmeans, method = metOptNB, k.max = kmax, nboot = 50)
      ggsave(paste(prout, metcluster, "_" , metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 8, width = 8)
    }else{
      p = fviz_nbclust(din, FUN = hcut, metho = metagregation,  method = metOptNB, k.max = kmax, nboot = 50)
      ggsave(paste(prout, metcluster, "_" , metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 8, width = 8)
    }
  }else if(metcluster == "kmeans"){
    p = fviz_nbclust(din, kmeans, method = metOptNB, k.max = kmax)
    ggsave(paste(prout, metcluster, "_" , metOptNB, ".png", sep = ""), dpi=300, height = 8, width = 15)    
  }
  
  
  if(metOptNB == "wss"){
    dcluster = as.matrix(p$data)
    d = inflexi(as.double(dcluster[,1]),as.double(dcluster[,2]),1,length(dcluster[,1]),3,3,plots=FALSE)
    nboptimal = d$finfl[1]
    print(nboptimal)
    
  }else if (metOptNB ==   "silhouette"){
    nboptimal = which(p$data[,2] == max(p$data[,2]))


  }else if (metOptNB == "gap_stat"){
    dcluster = as.matrix(p$data)
    #distorigin = abs(scale(as.double(dcluster[,5]), 0)-scale(-1*as.double(dcluster[,6]), 0))
    d = inflexi(as.double(dcluster[,5]),-1*as.double(dcluster[,6]),1,length(dcluster[,1]),3,3,plots=FALSE)
    nboptimal = d$finfl[1]
  }
  
  
  if (metcluster == "hclust"){
    outclust = hcut(din, k = nboptimal, hc_method = metagregation)
  }else if(metcluster == "kmeans"){
    outclust = hkmeans(din, nboptimal)
  }
  
  # PCA with clusters
  fviz_cluster(outclust, labelsize = 5)
  ggsave(paste(prout, "PCA_", metcluster, "_",  metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 12, width = 12)
  
  
  # dendogram fviz
  if(metcluster == "hclust"){
    fviz_dend(outclust, show_labels = FALSE, type = "circular")
    ggsave(paste(prout, "dendov1_", metcluster, "_",  metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 12, width = 13)
  }
  
  # dendogram old
  dcluster2 = cbind(names(outclust$cluster),outclust$cluster)
  colnames(dcluster2) = c("names", "cluster")
  dcluster2 = as.data.frame(dcluster2)
  rownames(dcluster2) = names(outclust$cluster)
  dcluster2 = dcluster2[rownames(din),]
  
  # save cluster
  write.csv(dcluster2, paste(prout, "clusters.csv", sep = ""), row.names = FALSE)
  
  
  
  # PLOT Dendogram !!
  ###################
  d <- dist(din, method = "euc")
  tupgma2 <- upgma(d, method = metagregation)
  #tupgma2 = groupOTU(tupgma2, nboptimal)
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% dcluster2 + geom_text(aes(color=cluster, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=2) +
    geom_tippoint(aes(color=cluster), alpha=0.75, size=1)+
    #scale_color_continuous(low='red', high='lightgreen') +
    #scale_color_manual(values=c("grey","red")) +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 5, y = 5, width = 10, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  #print(t4)
  ggsave(paste(prout, "dendo_", metcluster, "_",  metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 8, width = 15)
  
  
  d <- dist(din, method = "euc")
  tupgma2 <- upgma(d, method = metagregation)
  tupgma2 <- groupOTU(tupgma2, nboptimal)
  t4 <- ggtree(tupgma2, layout="circular", size=1, aes(color=cluster))
  t4 <- t4 %<+% dcluster2 + geom_text(aes(color=cluster, label=cluster, angle=angle, fontface="bold"), hjust=-0.15, size=2) +
    geom_tippoint(aes(color=cluster), alpha=0.75, size=1)+
    #scale_color_continuous(low='red', high='lightgreen') +
    #scale_color_manual(values=c("grey","red")) +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 5, y = 5, width = 10, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  #print(t4)
  ggsave(paste(prout, "dendo_cluster", metcluster, "_",  metagregation, "_", metOptNB, ".png", sep = ""), dpi=300, height = 8, width = 15)
  return(dcluster2)
  
}




###########
# MAIN    #
###########
args <- commandArgs(TRUE)
p_desc = args[1]
p_opera = args[2]
prout = args[3]


#p_desc = "../../ILS/results/Phthalates/rdkit/Cleaned_Data/desc1D2D_cleaned.csv"
#p_opera = "../../ILS/results/Phthalates/DESC/desc_OPERA.csv"
#prout = "../../ILS/results/Phthalates/rdkit/clustering/"



# d -> desc
ddesc = read.csv(p_desc, sep = ",", header = TRUE, row.names = 1)

# prediction OPERA
d_opera = read.csv(p_opera, sep = ",", header = TRUE, row.names = 1)
d_opera = d_opera[rownames(ddesc),]

# reduce OPERA desc
l_opera_pred = NULL
for(opera_desc in colnames(d_opera)){
  if(!is.integer0(grep("_pred", opera_desc, fixed = TRUE))){
    if(is.integer0(grep("_predRange", opera_desc, fixed = TRUE)) && is.integer0(grep("pKa_b_pred", opera_desc, fixed = TRUE)) && is.integer0(grep("pKa_a_pred", opera_desc, fixed = TRUE))){
      l_opera_pred = append(l_opera_pred, opera_desc)
    }
  }
}

d_opera = d_opera[,l_opera_pred]


#lmetclustering = c("hclust", "kmeans")
lmetclustering = c("hclust") # => to speed up the process
#lmetagregation = c("ward.D2", "complete", "single", "average") # - ward.D
lmetagregation = c("ward.D2") # - ward.D
#lmetoptimal = c("silhouette", "wss", "gap_stat")
lmetoptimal = c("gap_stat")

for(metclustering in lmetclustering){
  if (metclustering == "kmeans"){
    lmetagregation = c("ward.D2")
  }else{
    lmetagregation = c("ward.D2")
  }
  for (metagregation in lmetagregation){
    for(metoptimal in lmetoptimal){
      pclust = paste(prout, metclustering, "-", metagregation, "-", metoptimal, "/", sep = "")
      dir.create(pclust)
      dclust = optimalCluters(ddesc, pclust, metclustering, metoptimal, metagregation)
      print(dim(dclust))
      signifDescCluster(ddesc, d_opera, dclust, pclust)
      # write by chemical 
      dclust = dclust[rownames(d_opera),]
      cluster = dclust$cluster
      d_w = cbind(d_opera, cluster)
      write.csv(d_w, paste(pclust, "clusters_OPERA.csv", sep = ""), row.names = TRUE)
    }
  }
}
      




