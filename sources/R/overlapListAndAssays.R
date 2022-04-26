#!/usr/bin/env Rscript
library(UpSetR)
library(reshape2)
library(magrittr)
library(eulerr)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_assays = args[1]
p_chemList = args[2]
pr_out = args[3]


#p_assays = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/overlapToxCast/CYP19A1.csv"
#p_chemList = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/OverlapList/E2up-P4up/upset_matrix"
#pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/overlapToxCast/"


# load Assays
d_assays = read.csv(p_assays, sep = "\t", row.names = 1)
l_assays = colnames(d_assays)[seq(2,dim(d_assays)[2],2)]

# load chemList
d_chemList = read.csv(p_chemList, sep = "\t", row.names = 1)
l_chem_all = union(rownames(d_assays), rownames(d_chemList))

# write table for upset plot
d_out = NULL

for(chem in l_chem_all){
  l_rows = NULL
  
  # for E2 and P4 
  if(chem %in% rownames(d_chemList) == TRUE){
    l_rows = append(l_rows, d_chemList[chem, "E2up"])
    l_rows = append(l_rows, d_chemList[chem, "P4up"])
  }else{
    l_rows = append(l_rows, 0)
    l_rows = append(l_rows, 0)
  }
  
  # for assays
  if(chem %in% rownames(d_assays)){
    for (assay in l_assays){
      add = d_assays[chem,assay]
      if(add != "Active"){
        l_rows = append(l_rows, 0)
      }else{
        l_rows = append(l_rows, 1)  
      }
    }
  }else{
    for (assay in l_assays){
      l_rows = append(l_rows, 0)
    }
  }
  d_out = rbind(d_out, l_rows)
}

rownames(d_out) = l_chem_all
colnames(d_out) = c(colnames(d_chemList), l_assays)



png(paste(pr_out, "all_assays_venn.png", sep = ""), res = 300, 1000, 1000)
p = plot(euler(d_out), quantities = TRUE)
print(p)
dev.off()


# remove NVS assays
d_onlyTox21 = d_out[,-c(3, 4)]
png(paste(pr_out, "onlyTox21_venn.png", sep = ""), res = 300, 1000, 1000)
p = plot(euler(d_onlyTox21), quantities = TRUE)
print(p)
dev.off()
write.csv(d_onlyTox21, paste(pr_out, "upset_Tox21overlap.csv", sep = ""))


# remove Tox21 assays
d_onlyNVS = d_out[,-c(4, 5)]
png(paste(pr_out, "onlyNVS_venn.png", sep = ""), res = 300, 1000, 1000)
p = plot(euler(d_onlyNVS), quantities = TRUE)
print(p)
dev.off()

d_onlyNVS = cbind("chemical.name" = d_assays[rownames(d_onlyNVS), 1], d_onlyNVS)
write.csv(d_onlyNVS, paste(pr_out, "upset_NVSoverlap.csv", sep = ""))


# write d_out
d_out = cbind("chemical.name" = d_assays[l_chem_all, 1], d_out)
write.csv(d_out, paste(pr_out, "upset_all.csv", sep = ""))
