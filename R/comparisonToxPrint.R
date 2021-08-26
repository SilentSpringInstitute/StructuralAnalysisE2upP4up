#!/usr/bin/env Rscript
library(Toolbox)
library(igraph)
library(ggraph)
library(scales)

drawHist = function(d_toxprint, max_pval, min_prob, p_out){
  
  # reduce d_toxprint
  d_toxprint = d_toxprint[which(d_toxprint[, "Pval"] < max_pval  | d_toxprint[, 6] > min_prob | d_toxprint[, 7] > min_prob),]
  
  d_plot = NULL
  for(id in rownames(d_toxprint)){
    d_plot = rbind(d_plot, c(id, d_toxprint[id, c(1, 2)], "",  d_toxprint[id, c(4, 6)], colnames(d_toxprint)[4]))
    d_plot = rbind(d_plot, c(id, d_toxprint[id, c(1, 2, 3, 5, 7)], colnames(d_toxprint)[5]))
  }
  colnames(d_plot) = c("ID", "ToxPrint", "Pval", "Significative", "nb", "prob", "set")
  d_plot = as.data.frame(d_plot)
  d_plot$ToxPrint = as.character(d_plot$ToxPrint)
  d_plot$set = as.character(d_plot$set)
  d_plot$nb = as.double(as.character(d_plot$nb))
  d_plot$prob = as.double(as.character(d_plot$prob))
  d_plot$Pval = as.double(as.character(d_plot$Pval))
  
  if (colnames(d_toxprint)[4] == "E2up"){
    rev = TRUE
    v_col = c( "#1f78b4", "#a6cee3")
  }else{
    rev = FALSE
    v_col = c("#a6cee3", "#1f78b4")
  }
  
  
  ggplot(data = d_plot, aes(x=reorder(ToxPrint, -Pval), y=prob, fill=set)) +
    geom_bar(stat="identity", position=position_dodge2(reverse = rev))+
    geom_text(aes(label=Significative), vjust=1.6, color="black", position = position_dodge(0.6), size=3.5)+
    scale_fill_manual(values=v_col)+
    coord_flip() + 
    labs(y="Prob", x="Toxprint")
  
  if(dim(d_plot)[1] > 100){
    ggsave(p_out,  width = 10, height = 15, dpi = 300)
  }else{
    ggsave(p_out,  width = 10, height = 10, dpi = 300)  
  }
}




drawNet = function(d_alltoxprint, d_signif, cutoff_signif, p_out){
  
  # reduce with significative toxprint
  l_toxprints_signif = d_signif$Toxprint[which(d_signif$Pval < cutoff_signif & d_signif[,6] > d_signif[,7])]
  d_alltoxprint = d_alltoxprint[, l_toxprints_signif]
  
  # define the node
  d_node = NULL
  l_toxprints = NULL
  i = 1
  for(toxprint in colnames(d_alltoxprint)){
    nb_chem = sum(d_alltoxprint[,toxprint])
    if(nb_chem > 1){
      id = paste("s", i, sep="")
      d_node = rbind(d_node, c(id, toxprint, nb_chem))
      l_toxprints = append(l_toxprints, toxprint)
      i = i + 1 
    }
  }
  
  colnames(d_node) = c("id", "Toxprint", "size")
  d_node = as.data.frame(d_node)
  d_node$size = as.double(as.character(d_node$size))
  names(l_toxprints) = d_node$id
  
  
  
  # define the edge
  d_edge = NULL
  i = 1
  imax = length(l_toxprints)
  while(i <= imax){
    j = i + 1
    while(j < imax){
      weight = length(which(d_alltoxprint[, l_toxprints[i]] == 1 & d_alltoxprint[, l_toxprints[j]] == 1))
      if(weight > 1){
        d_edge = rbind(d_edge, c(names(l_toxprints[i]), names(l_toxprints[j]), "mention" ,weight))
      }
      j = j + 1
    }
    i = i + 1
  }
  
  colnames(d_edge) = c("from", "to", "type", "weight")
  d_edge = as.data.frame(d_edge)
  d_edge$weight = as.numeric(as.character(d_edge$weight))
  
  net <- graph_from_data_frame(d=d_edge, vertices=d_node, directed=T)
  
  
 ggraph(net, layout = 'linear') + 
    geom_edge_arc(aes(color = weight, alpha=weight), width = 0.7) +
    scale_edge_alpha(name = "Count chem. sharing")+
    scale_edge_color_gradientn(name = "Count chem. sharing",colors = c("cyan", "blue"), values = rescale(c(10,20,42))) +
    geom_node_point(aes(size = size, color = size)) +
    scale_color_gradientn(name = "Node (chem. count)", colors = c("cyan", "darkblue"), values = rescale(c(10,35,60))) +
    scale_size_continuous(name = "Node (chem. count)")+ 
    geom_node_text(aes(label = d_node$Toxprint), hjust=1, vjust=1, angle=45, nudge_y = -0.2,  color="black") +
    theme_void()+
    ylim(-4.5, NA)+xlim(-4.5, NA)
   
 ggsave(p_out,  width = 16, height = 10, dpi = 300)
 
 
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_toxprint1 = args[1]
p_toxprint2 = args[2]
pr_out = args[3]


#p_toxprint1 = "/mnt/c/Users/AlexandreBorrel/research/SSI/E2up_P4up/results/comparisonDescToxprint_E2up-H295R/E2up_toxprint.csv"
#p_toxprint2 = "/mnt/c/Users/AlexandreBorrel/research/SSI/E2up_P4up/results/comparisonDescToxprint_E2up-H295R/H295R_toxprint.csv"
#pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/E2up_P4up/results/comparisonDescToxprint_E2up-H295R/Toxprint/"

d_toxprint1 = read.csv(p_toxprint1, sep = "\t", row.names = 1)


# check the openning
if (dim(d_toxprint1)[2] == 0){
  d_toxprint1 = read.csv(p_toxprint1, sep = ",", row.names = 1)
}

d_toxprint2 = read.csv(p_toxprint2, sep = "\t",  row.names = 1)
if (dim(d_toxprint2)[2] == 0){
  d_toxprint2 = read.csv(p_toxprint2, sep = ",", row.names = 1)
}

#Drop SMILES col
if("SMILES" %in% colnames(d_toxprint1)){
  drops = c("SMILES")
  d_toxprint1 = d_toxprint1[ , !(names(d_toxprint1) %in% drops)]
}


l_toxprints = intersect(colnames(d_toxprint1), colnames(d_toxprint2))
n_d_toxprint1 = dim(d_toxprint1)[1]
n_d_toxprint2 = dim(d_toxprint2)[1]


d_out = NULL
for(toxprint in l_toxprints){
  
  v_toxprint1 = d_toxprint1[,toxprint]
  v_toxprint2 = d_toxprint2[,toxprint]
  
  n_toxprint1 = sum(v_toxprint1)
  n_toxprint2 = sum(v_toxprint2)
  
  # compute pval on a Zscore
  res <- prop.test(x = c(n_toxprint1, n_toxprint2), n = c(n_d_toxprint1, n_d_toxprint2))
  
  # Printing the results
  pval = res$p.val 
  if (is.na(pval) == FALSE){
    signif = signifPvalue(pval)
    d_out = rbind(d_out, c(toxprint, round(pval, 4), signif, n_toxprint1, n_toxprint2, round(res$estimate[1], 2), round(res$estimate[2], 2)))
  }
}


colnames(d_out) = c("Toxprint", "Pval", "significatif", strsplit(basename(p_toxprint1), "_")[[1]][1], strsplit(basename(p_toxprint2), "_")[[1]][1], paste(strsplit(basename(p_toxprint1), "_")[[1]][1], " prob", sep = ""), paste(strsplit(basename(p_toxprint2), "_")[[1]][1], " prob", sep = ""))
d_out = as.data.frame(d_out)
d_out$Pval = as.double(d_out$Pval)
d_out = d_out[order(d_out$Pval), ]
write.csv(d_out, paste(pr_out, "signif.csv", sep = ""))

# make plot 
# only  > ** significant and more than 20% of the set
drawHist(d_out, 0.01, 0.2, paste(pr_out, "toxprint_signif.png", sep = ""))

###
# make network with cutoff of significativity
drawNet (d_toxprint1, d_out, 0.01, paste(pr_out, "network_001.png", sep = ""))
drawNet (d_toxprint1, d_out, 0.001, paste(pr_out, "network_0001.png", sep = ""))


### compute the average number of Toxprint by chemicals
d_sum1 = apply(d_toxprint1, 1, "sum")
d_sum2 = apply(d_toxprint2, 1, "sum")

Av1 = mean(d_sum1)
sd1 = sd(d_sum1)

Av2 = mean(d_sum2)
sd2 = sd(d_sum2)

d_tox_signif = d_toxprint1[,d_out[which(d_out$Pval < 0.01),"Toxprint"]]
d_sum_signif = apply(d_tox_signif, 1, "sum")
Av_signif = mean(d_sum_signif)

d_sum_out = c(Av1, sd1, Av2, sd2, Av_signif) 
names(d_sum_out) = c("Avg nb toxprint 1", "SD nb toxprint 1", "Avg nb toxprint 2", "SD nb toxprint 2", "Avg signif Toxprint")
write.table(d_sum_out, paste(pr_out, "countToxprint.sum", sep =""))

#### combination of toxprint
######

# remove toxprint with a expected prob < 0.1
l_toxprint_tocombine = d_out[which(d_out[,6] != "NA" & d_out[,7] != "NA"),1]

# do combination of toxprint
########

d_combine = NULL
l_combine = NULL
i = 1
imax = length(l_toxprint_tocombine)
while(i <= imax){
  j = i + 1
  while(j <= imax){
    toxprint_combine = paste(l_toxprint_tocombine[i], l_toxprint_tocombine[j], sep = "+")
    l_combine = append(l_combine, toxprint_combine)
    
    # apply test
    v_toxprint1 = length(which(d_toxprint1[,l_toxprint_tocombine[i]] >=1 & d_toxprint1[,l_toxprint_tocombine[j]] >=1))
    v_toxprint2 = length(which(d_toxprint2[,l_toxprint_tocombine[i]] >=1 & d_toxprint2[,l_toxprint_tocombine[j]] >=1))
    
    n_toxprint1 = sum(v_toxprint1)
    n_toxprint2 = sum(v_toxprint2)
    
    # compute pval on a Zscore
    res <- prop.test(x = c(n_toxprint1, n_toxprint2), n = c(n_d_toxprint1, n_d_toxprint2))
    
    # only take when pval is not NA => count > 1
    pval = res$p.val 
    if (is.na(pval) == FALSE){
      signif = signifPvalue(pval)
      d_combine = rbind(d_combine, c(toxprint_combine, round(pval, 4), signif, n_toxprint1, n_toxprint2, round(res$estimate[1], 2), round(res$estimate[2], 2)))
    }
    j = j + 1
  }
  i = i + 1
}


colnames(d_combine) = c("Toxprint", "Pval", "significatif", strsplit(basename(p_toxprint1), "_")[[1]][1], strsplit(basename(p_toxprint2), "_")[[1]][1], paste(strsplit(basename(p_toxprint1), "_")[[1]][1], " prob", sep = ""), paste(strsplit(basename(p_toxprint2), "_")[[1]][1], " prob", sep = ""))
d_combine = as.data.frame(d_combine)
d_combine$Pval = as.double(d_combine$Pval)
d_combine = d_combine[order(d_combine$Pval), ]
write.csv(d_combine, paste(pr_out, "combined_signif.csv", sep = ""))


drawHist(d_combine, 0.01, 0.20, paste(pr_out, "combine_signif.png", sep = ""))

