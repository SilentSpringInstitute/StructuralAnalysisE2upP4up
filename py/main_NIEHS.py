import pathFolder
import dataset
import analysis
from os import path
import toolbox
import mapOnSpace

import PFAS

# Define folder #
#################
PR_ROOT = "../../ILS/"
PR_RESULTS = PR_ROOT + "results/"
PR_DATA = PR_ROOT + "data/"

COR_VAL = 0.90
MAX_QUANTILE = 90

P_COORD_1D2D_PFAS = PR_DATA + "PFAS_coord1D2D.csv"
P_COORD_3D_PFAS = PR_DATA + "PFAS_coord3D.csv"

P_COORD_1D2D_Tox21 = PR_DATA + "Tox21_coord1D2D.csv"
P_COORD_3D_Tox21 = PR_DATA + "Tox21_coord3D.csv"


# Extract data and compute desc #
#################################

########## Phthalates
#####################
#p_dataset = PR_DATA + "Phthalates.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "Phthalates/")
# RDKIT
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

# OPERA
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

# OPERA + rdkit
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit", "opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 0)

###### map on space ===> need some update on the chemmaps 
######
#p_RDKIT = pr_results + "rdkit/rdkit.csv"
#pr_results_mapped = pathFolder.createFolder(pr_results + "mapped_on_PFASMap/")
#c_mapped = mapOnSpace.mapOnSpace(p_RDKIT, P_COORD_1D2D_PFAS, P_COORD_3D_PFAS, pr_results_mapped)
#c_mapped.map("Phthalates")

#pr_results_mapped = pathFolder.createFolder(pr_results + "mapped_on_Tox21Map/")
#c_mapped = mapOnSpace.mapOnSpace(p_RDKIT, P_COORD_1D2D_Tox21, P_COORD_3D_Tox21, pr_results_mapped)
#c_mapped.map("Phthalates")


########## Phthalates alternative
#################################
#p_dataset = PR_DATA + "Phthalates_alternatives.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "Phthalates_alternatives/")
## RDKIT
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

## OPERA
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

## OPERA + rdkit
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit", "opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 0)


###### map on space ===> need some update on the chemmaps 
######
#p_RDKIT = pr_results + "rdkit/rdkit.csv"



############### PFAS
#####################
#p_dataset = PR_DATA + "PFAS.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "PFAS/")
## RDKIT
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

## OPERA
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

# OPERA + rdkit
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit", "opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 0)


###### map on space ===> need some update on the chemmaps 
######
#p_RDKIT = pr_results + "rdkit/rdkit.csv"




################# Combine dataset PFAS + Phthalates
###################################################

#p_dataset1 = PR_DATA + "PFAS.csv"
#p_dataset2 = PR_DATA + "Phthalates.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "PFAS-Phthalates/")

#cPFAS = PFAS.PFAS(p_dataset1, pr_results)
#cPFAS.combineDataset(p_dataset2)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

# OPERA
#cPFAS = PFAS.PFAS(p_dataset1, pr_results)
#cPFAS.combineDataset(p_dataset2)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

# OPERA + rdkit
#cPFAS = PFAS.PFAS(p_dataset1, pr_results)
#cPFAS.combineDataset(p_dataset2)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit", "opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 0)

###### map on space ===> need some update on the chemmaps 
######
#p_RDKIT = pr_results + "rdkit/rdkit.csv"



################# Combine dataset Phthalates + Phthalates alternative
######################################################################

#p_dataset1 = PR_DATA + "Phthalates_alternatives.csv"
#p_dataset2 = PR_DATA + "Phthalates.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "Phthalates_alternatives-Phthalates/")

#cPFAS = PFAS.PFAS(p_dataset1, pr_results)
#cPFAS.combineDataset(p_dataset2)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

## OPERA
#cPFAS = PFAS.PFAS(p_dataset1, pr_results)
#cPFAS.combineDataset(p_dataset2)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

## OPERA + rdkit
#cPFAS = PFAS.PFAS(p_dataset1, pr_results)
#cPFAS.combineDataset(p_dataset2)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit", "opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 0)



################# Combine dataset Phthalates-alternative + PFAS
###################################################

#p_dataset1 = PR_DATA + "PFAS.csv"
#p_dataset2 = PR_DATA + "Phthalates.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "PFAS-Phthalates/")

#cPFAS = PFAS.PFAS(p_dataset1, pr_results)
#cPFAS.combineDataset(p_dataset2)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

# OPERA
#cPFAS = PFAS.PFAS(p_dataset1, pr_results)
#cPFAS.combineDataset(p_dataset2)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

# OPERA + rdkit
#cPFAS = PFAS.PFAS(p_dataset1, pr_results)
#cPFAS.combineDataset(p_dataset2)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit", "opera"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 0)

###### map on space ===> need some update on the chemmaps 
######
#p_RDKIT = pr_results + "rdkit/rdkit.csv"




################# Combine Phthalates + Phthalates alternative + PFAS
####################################################################

p_dataset1 = PR_DATA + "Phthalates_alternatives.csv"
p_dataset2 = PR_DATA + "PFAS_phthalates.csv"
pr_results = pathFolder.createFolder(PR_RESULTS + "Phthalates-Phthalates_alternative-PFAS/")

cPFAS = PFAS.PFAS(p_dataset1, pr_results)
cPFAS.combineDataset(p_dataset2)
cPFAS.computeDesc()
cPFAS.buildDescSet(["rdkit"])
cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

# OPERA
cPFAS = PFAS.PFAS(p_dataset1, pr_results)
cPFAS.combineDataset(p_dataset2)
cPFAS.computeDesc()
cPFAS.buildDescSet(["opera"])
cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 1)

# OPERA + rdkit
cPFAS = PFAS.PFAS(p_dataset1, pr_results)
cPFAS.combineDataset(p_dataset2)
cPFAS.computeDesc()
cPFAS.buildDescSet(["rdkit", "opera"])
cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE,  PCA=1, Hclust=1, clustering=1, histDesc = 0)

###### map on space ===> need some update on the chemmaps 
######
#p_RDKIT = pr_results + "rdkit/rdkit.csv"

