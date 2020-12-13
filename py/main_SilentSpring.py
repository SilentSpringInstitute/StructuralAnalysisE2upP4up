from os import path

import pathFolder
import dataset
import analysis
import Chemicals



# Define folder #
#################
PR_ROOT = path.abspath("../../") + "/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")
PR_DATA = PR_ROOT + "data/"

COR_VAL = 0.9
MAX_QUANTILE = 90

# Extract data and compute desc #
#################################
# PFAS
#p_listChem = PR_DATA + "PFAS_45_updated.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "PFAS-45/")

#p_listChem = PR_DATA + "Master_chemical_List_7-29-20.csv"
#pr_results = pathFolder.createFolder(PR_RESULTS + "Master_list/")

# carcinogen
p_listChem = PR_DATA + "2021_updated_List_MCs_120720.csv"
pr_results = pathFolder.createFolder(PR_RESULTS + "Carcinogen_list/")

# Compute desc
###################
cChem = Chemicals.Chemicals(p_listChem, pr_results)
cChem.computeDesc()
cChem.buildDescSet(["rdkit"])
cChem.analysisDescBasedData(COR_VAL, MAX_QUANTILE, PCA=0, Hclust=0, clustering=0, FP=0, SOM=1)


#Extract only tier 1 and 2
#cChem.getChemTier(["1","2"])
#cChem.analysisDescBasedData(COR_VAL, MAX_QUANTILE, PCA=0, Hclust=0, clustering=0, FP=1)



