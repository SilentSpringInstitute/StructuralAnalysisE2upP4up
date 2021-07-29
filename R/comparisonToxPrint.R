#!/usr/bin/env Rscript
library(Toolbox)




################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_toxprint1 = args[1]
p_toxprint2 = args[2]
pr_out = args[3]


#p_toxprint1 = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/comparisonDesc_E2up-H295R/E2up_toxprint.csv"
#p_toxprint2 = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/comparisonDesc_E2up-H295R/H295R_toxprint.csv"
#pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/comparisonDesc_E2up-H295R/toxprint"

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
  if (is.na(pval) == TRUE){
    d_out = rbind(d_out, c(toxprint, "NA", "-", n_toxprint1, n_toxprint2, "NA", "NA"))
  
  }else{
    signif = signifPvalue(pval)
    d_out = rbind(d_out, c(toxprint, round(pval, 4), signif, n_toxprint1, n_toxprint2, round(res$estimate[1], 2), round(res$estimate[2], 2)))
  }
}


colnames(d_out) = c("Toxprint", "Pval", "significatif", paste("N", strsplit(basename(p_toxprint1), "_")[[1]][1], sep = " "), paste("N", strsplit(basename(p_toxprint2), "_")[[1]][1], sep = " "), "Estimate prob1", "Estimate prob2")
d_out = as.data.frame(d_out)
d_out$Pval = as.double(d_out$Pval)
d_out = d_out[order(d_out$Pval), ]
write.csv(d_out, paste(pr_out, "_signif.csv"))
