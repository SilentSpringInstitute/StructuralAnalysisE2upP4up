#!/usr/bin/env Rscript
source ("../../../../ILS/development/R_toolbox/dataManager.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
prout = args[2]
valcor = args[3]
maxquantile = as.double(args[4])

#pdesc = "../../results/DESC/desc_1D2D.csv"
#prout = "../../results/Cleaned_Data/"
#valcor = 0.9
#maxquantile = 90



##############################
# Process descriptors matrix #
##############################

dglobal = openDataVexcluded(pdesc, valcor, prout, c(1,2))
dglobal = dglobal[[1]]
lCAS = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]

# format in double
dglobal =apply(dglobal, 2, as.numeric)
rownames(dglobal) = lCAS


print("==== Preprocessing ====")
print(paste("Data initial: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

##########
# filter #
##########

dglobal = delnohomogeniousdistribution(dglobal, maxquantile)
print(paste("Data after filtering: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))


pdesout = paste(prout, "desc1D2D_cleaned.csv", sep = "")
write.csv(dglobal, pdesout, col.names = TRUE, row.names = TRUE)