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
p_dataset = PR_DATA + "Phthalates.csv"
pr_results = pathFolder.createFolder(PR_RESULTS + "Phthalates/")
#cPFAS = PFAS.PFAS(p_dataset, pr_results)
#cPFAS.computeDesc()
#cPFAS.buildDescSet(["rdkit"])
#cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE, PCA=1, Hclust=1)

#====> need to be updates
# analysis with RDKIT desc
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["opera"], pr_results_analysis)
#analysisDescBasedData(l_file_desc[0], l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

# analysis with RDKIT + OPERA physico chemical
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT-OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit", "opera"], pr_results_analysis)
#analysisDescBasedData(p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

# analysis with OPERA physico chemical only
#pr_results_analysis = pathFolder.createFolder(pr_results + "OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["opera"], pr_results_analysis)
#analysisDescBasedData(p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)


### map on space
#####
#pr_results_mapped = pathFolder.createFolder(pr_results + "mapped_PFAS/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_PFAS, P_COORD_3D_PFAS, pr_results_mapped)
#c_mapped.map("Phthalates")

#pr_results_mapped = pathFolder.createFolder(pr_results + "mapped_Tox21/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_Tox21, P_COORD_3D_Tox21, pr_results_mapped)
#c_mapped.map("Phthalates")


########## Phthalates alternative
#################################
#p_dataset = PR_DATA + "Phthalates_alternatives.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "Phthalates_alternatives/")
#l_file_desc = computeDesc(p_dataset, pr_results)

## analysis with RDKIT desc
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit"], pr_results_analysis)
#analysisDescBasedData(p_dataset, p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

## analysis with RDKIT + OPERA physico chemical
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT-OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit", "opera"], pr_results_analysis)
#nalysisDescBasedData(p_dataset, p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

## analysis with OPERA physico chemical only
#pr_results_analysis = pathFolder.createFolder(pr_results + "OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["opera"], pr_results_analysis)
#analysisDescBasedData(p_dataset, p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)


### map on space
#####
#pr_results_mapped = pathFolder.createFolder(pr_results + "mappedPFAS/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_PFAS, P_COORD_3D_PFAS, pr_results_mapped)
#c_mapped.map("Phthalates_alternatives")

#pr_results_mapped = pathFolder.createFolder(pr_results + "mappedTox21/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_Tox21, P_COORD_3D_Tox21, pr_results_mapped)
#c_mapped.map("Phthalates_alternatives")


############### PFAS
#####################
#p_dataset = PR_DATA + "PFAS.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "PFAS/")
#l_file_desc = computeDesc(p_dataset, pr_results)

## analysis with RDKIT desc
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit"], pr_results_analysis)
#analysisDescBasedData(l_file_desc[0], l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

## analysis with RDKIT + OPERA physico chemical
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT-OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit", "opera"], pr_results_analysis)
#analysisDescBasedData(p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

## analysis with OPERA physico chemical only
#pr_results_analysis = pathFolder.createFolder(pr_results + "OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["opera"], pr_results_analysis)
#analysisDescBasedData(p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)


#pr_results_mapped = pathFolder.createFolder(pr_results + "PFAS_mapped/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_PFAS, P_COORD_3D_PFAS, pr_results_mapped)
#c_mapped.map("PFAS-EPA-PFAS")


#pr_results_mapped = pathFolder.createFolder(pr_results + "Tox21_mapped/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_Tox21, P_COORD_3D_Tox21, pr_results_mapped)
#c_mapped.map("Tox21_PFAS")



################# Combine dataset PFAS + Phthalates
###################################################

p_dataset1 = PR_DATA + "PFAS.csv"
p_dataset2 = PR_DATA + "Phthalates.csv"
pr_results = pathFolder.createFolder(PR_RESULTS + "PFAS-Phthalates/")

cPFAS = PFAS.PFAS(p_dataset, pr_results)
cPFAS.computeDesc()
cPFAS.combineDataset(p_dataset2)
cPFAS.buildDescSet(["rdkit"])
cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE, PCA=1, Hclust=1)


# =====> need to be update
## analysis with RDKIT + OPERA physico chemical
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT-OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit", "opera"], pr_results_analysis)
#analysisDescBasedData(p_dataset,p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

## analysis with OPERA physico chemical only
#pr_results_analysis = pathFolder.createFolder(pr_results + "OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["opera"], pr_results_analysis)
#analysisDescBasedData(p_dataset, p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)


#pr_results_mapped = pathFolder.createFolder(pr_results + "PFAS_mapped/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_PFAS, P_COORD_3D_PFAS, pr_results_mapped)
#c_mapped.map("PFAS-EPA-PFAS")


#pr_results_mapped = pathFolder.createFolder(pr_results + "Tox21_mapped/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_Tox21, P_COORD_3D_Tox21, pr_results_mapped)
#c_mapped.map("Tox21_PFAS")



################# Combine dataset Phthalates + Phthalates alternative
######################################################################

p_dataset1 = PR_DATA + "Phthalates_alternatives.csv"
p_dataset2 = PR_DATA + "Phthalates.csv"
pr_results = pathFolder.createFolder(PR_RESULTS + "Phthalates_alternatives-Phthalates/")

p_dataset = combineDataset(p_dataset1, p_dataset2, pr_results)
l_file_desc = computeDesc(p_dataset, pr_results)


## analysis with RDKIT desc
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit"], pr_results_analysis)
#analysisDescBasedData(p_dataset, l_file_desc[0], l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

## analysis with RDKIT + OPERA physico chemical
pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT-OPERA_desc/")
p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit", "opera"], pr_results_analysis)
analysisDescBasedData(p_dataset,p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)
dd
## analysis with OPERA physico chemical only
#pr_results_analysis = pathFolder.createFolder(pr_results + "OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["opera"], pr_results_analysis)
#analysisDescBasedData(p_dataset, p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

#pr_results_mapped = pathFolder.createFolder(pr_results + "PFAS_mapped/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_PFAS, P_COORD_3D_PFAS, pr_results_mapped)
#c_mapped.map("PFAS-EPA-PFAS")

#pr_results_mapped = pathFolder.createFolder(pr_results + "Tox21_mapped/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_Tox21, P_COORD_3D_Tox21, pr_results_mapped)
#c_mapped.map("Tox21_PFAS")


################# Combine all
#############################

p_dataset1 = PR_DATA + "Phthalates_alternatives.csv"
p_dataset2 = PR_DATA + "PFAS_phthalates.csv"
pr_results = pathFolder.createFolder(PR_RESULTS + "All/")

p_dataset = combineDataset(p_dataset1, p_dataset2, pr_results)
l_file_desc = computeDesc(p_dataset, pr_results)

# analysis with RDKIT desc
pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT_desc/")
p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit"], pr_results_analysis)
analysisDescBasedData(p_dataset, l_file_desc[0], l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

# analysis with RDKIT + OPERA physico chemical
pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT-OPERA_desc/")
p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit", "opera"], pr_results_analysis)
analysisDescBasedData(p_dataset,p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

# analysis with OPERA physico chemical only
pr_results_analysis = pathFolder.createFolder(pr_results + "OPERA_desc/")
p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["opera"], pr_results_analysis)
analysisDescBasedData(p_dataset, p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)


pr_results_mapped = pathFolder.createFolder(pr_results + "PFAS_mapped/")
c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_PFAS, P_COORD_3D_PFAS, pr_results_mapped)
c_mapped.map("PFAS-EPA-PFAS")


pr_results_mapped = pathFolder.createFolder(pr_results + "Tox21_mapped/")
c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D_Tox21, P_COORD_3D_Tox21, pr_results_mapped)
c_mapped.map("Tox21_PFAS")
