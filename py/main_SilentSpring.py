import pathFolder
import dataset
import analysis
import PFAS



# Define folder #
#################
PR_ROOT = "../../Silent_Spring/"
PR_RESULTS = PR_ROOT + "results/"
PR_DATA = PR_ROOT + "data/"

COR_VAL = 0.9
MAX_QUANTILE = 90

# Extract data and compute desc #
#################################
p_sub45set = PR_DATA + "PFAS_45_updated.csv"
p_masterlist = PR_DATA + "Master_chemical_List_7-29-20.csv"


# for master list
###################
pr_results = pathFolder.createFolder(PR_RESULTS + "Master_list/")
cPFAS = PFAS.PFAS(p_masterlist, pr_results)
cPFAS.computeDesc()

cPFAS.buildDescSet(["rdkit"])
cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE, PCA=0, Hclust=0)

#Extract only tier 1 and 2
cPFAS.getChemTier(["1","2"])
cPFAS.analysisDescBasedData(COR_VAL, MAX_QUANTILE, PCA=0, Hclust=0, clustering=0, FP=1)



