#!/usr/bin/env Rscript
library(hexbin)
library(RColorBrewer)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pcoordMain1D2D = args[1]
pcoordMain3D = args[2]
plistChem = args[3]
name_map = args[4]
prout = args[5]



#pcoordMain1D2D = '../../data/Tox21_coord1D2D.csv'
#pcoordMain3D =  '../../data/Tox21_coord3D.csv'
#plistChem = '../../results/PFAS/Tox21_mapped/list_chem.csv'
#name_map = 'Tox21-PFAS'
#prout = '../../results/PFAS/Tox21_mapped/'



# linux 
#p_listchem, p_1D2D, p_3D, pr_out
# window
#pcoordMain1D2D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord1D2D_light.csv"
#pcoordMain3D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord3D_light.csv"
#plistChem = "/Users/Aborrel/research/NIEHS/ChemMaps/data/NPAHs.csv"
#name_map = "NPAH_test"
#prout = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/NPAHs/"


dcoord1D2D = read.csv(pcoordMain1D2D, sep = ",", header = TRUE)
dcoord1D2D = dcoord1D2D[,c(1,2, 3)]
rownames(dcoord1D2D) = dcoord1D2D[,1]
dcoord1D2D = dcoord1D2D[,-1]

dcoord3D = read.csv(pcoordMain3D, sep = ",", header = TRUE)
dcoord3D = dcoord3D[,c(1,2, 3)]
rownames(dcoord3D) = dcoord3D[,1]
dcoord3D = dcoord3D[,-1]

lchem = intersect(rownames(dcoord1D2D), rownames(dcoord3D))
dcoord3D = dcoord3D[lchem,]
dcoord1D2D = dcoord1D2D[lchem,]

# coord dim1 vs dim3
dcoord_xzD = data.frame(DIM1 = dcoord1D2D$DIM1, DIM3 = dcoord3D$DIM3.1)




# cut the map with list of chemical
dchemIn = read.csv(plistChem, sep = "\t")
lchem = dchemIn$inchikey
lchem = na.omit(lchem)
lchem = unique(lchem)


# color
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(100)



# plot full map #
#################

dcoord_full = dcoord1D2D[,c(1,2)]
dcoord_full = cbind(dcoord_full, dcoord3D[,1])
colnames(dcoord_full) = c("DIM1", "DIM2", "DIM3")


png(paste(prout, "Full_PFAS_xy.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoord_full, xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 50, cex.labels = 0.6, cex.title=0.6 )
dev.off()


png(paste(prout, "Full_PFAS_xz.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM3~DIM1,data=dcoord_full,xlab="DIM1",
           ylab="DIM3",colramp=rf, xbins = 50, cex.labels = 0.6, cex.title=0.6 )
dev.off()


#######################
# define coordinate for list chem
#######################
dcoorMap = cbind(dcoord1D2D[lchem,], dcoord3D[lchem,1])
colnames(dcoorMap) = c("DIM1", "DIM2", "DIM3")


# plot => DIM1 VS DIM2 #
########################

# calibrate the map
dcoorMap = rbind(dcoorMap, c(min(dcoord1D2D$DIM1), min(dcoord1D2D$DIM2), min(dcoord3D$DIM3.1)))
dcoorMap = rbind(dcoorMap, c(max(dcoord1D2D$DIM1), max(dcoord1D2D$DIM2), max(dcoord3D$DIM3.1)))


png(paste(prout, name_map, "_xy.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoorMap,xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 50, cex.labels = 0.6, cex.title=0.6 )
dev.off()




# plot => DIM1 VS DIM3 #
########################
# calibrate the map
png(paste(prout, name_map, "_xz.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM3~DIM1,data=dcoorMap,xlab="DIM1",
           ylab="DIM3",colramp=rf, xbins = 50, cex.labels = 0.6, cex.title=0.6 )
dev.off()


print(lchem)
lchem = append(as.character(lchem), "end")
lchem = append(lchem, "begin")
print(lchem)
rownames(dcoorMap) = lchem

write.csv(x = dcoorMap, file = paste(prout , "coords.csv", sep = ""))

