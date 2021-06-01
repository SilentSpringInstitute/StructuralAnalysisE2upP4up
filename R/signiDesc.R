#!/usr/bin/env Rscript
library(Toolbox)





################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc1 = args[1]
p_desc2 = args[2]
pr_out = args[3]


#p_desc1 = "./../../results/desc_by_list/E2-up/desc_1D2D.csv"
#p_desc2 = "./../../results/desc_by_list/P4-up/desc_1D2D.csv"
#pr_out = "./../../results/comparisonDesc_E2-up-P4-up/"

d_desc1 = read.csv(p_desc1, sep = "\t", row.names = 1)

if (dim(d_desc1)[2] == 0){
  d_desc1 = read.csv(p_desc1, sep = ",", row.names = 1)
}
d_desc2 = read.csv(p_desc2, sep = "\t",  row.names = 1)
if (dim(d_desc2)[2] == 0){
  d_desc2 = read.csv(p_desc2, sep = ",", row.names = 1)
}

if("SMILES" %in% colnames(d_desc1)){
  drops = c("SMILES")
  d_desc1 = d_desc1[ , !(names(d_desc1) %in% drops)]
}


l_desc = intersect(colnames(d_desc1), colnames(d_desc2))


d_out = NULL
for(desc in l_desc){
  
  v_desc1 = d_desc1[,desc]
  v_desc2 = d_desc2[,desc]
  
  isparametric = conditionTtest(v_desc1, v_desc2)
  if(isparametric == 1){
    pval = comparisonTest (v_desc1, v_desc2, "parametric")
  }else{
    pval = comparisonTest (v_desc1, v_desc2, "no-parametric")
  }
  if (is.na(pval) == TRUE){
    d_out = rbind(d_out, c(desc, "NA", "-", "NA"))
  }else{
    signif = signifPvalue(pval)
    d_out = rbind(d_out, c(desc, round(pval, 4), signif, isparametric))
  }

}

colnames(d_out) = c("Desc", "Pval", "significatif", "isparametric")
d_out = d_out[order(d_out[,2]), ]

write.csv(d_out, paste(pr_out, "_signif_desc.csv"))
