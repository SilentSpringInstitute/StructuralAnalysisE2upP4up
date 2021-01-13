from os import path

import pathFolder
import dataset
import analysis
import Chemicals
import comparisonChemicalLists


# Define folder #
#################
PR_ROOT = path.abspath("../../") + "/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")
PR_DATA = PR_ROOT + "data/"

# val
COR_VAL = 0.9
MAX_QUANTILE = 90

def computeDescAnalysisFromList(p_list, p_ToxPrint, p_chemList, PR_RESULTS):

    pr_results = pathFolder.createFolder(PR_RESULTS + p_list.split("/")[-1][0:-4] + "/")

    # Compute desc
    ###################
    cChem = Chemicals.Chemicals(p_list, pr_results)
    cChem.computeDesc()
    #cChem.buildDescSet(["rdkit"])
    #cChem.analysisDescBasedData(COR_VAL, MAX_QUANTILE, PCA=1, Hclust=1, clustering=1, FP=1, SOM=1, histDesc=1)

    cChem.analysisToxPrint(p_ToxPrint)
    cChem.analysisChemList(p_chemList)






# list of chemicals available
#############
p_list_carcinogen = PR_DATA + "carcinogen_breast_122120.csv"
p_list_e2 = PR_DATA + "H295_E2-up.csv"
p_list_p4 = PR_DATA + "H295_P4-up.csv"


# Explore ToxPrint - carcinogen
## carcinogen
p_ToxPrint_carcinogen = PR_ROOT + 'comptox/carcinogen_breast_122120_ToxPrint.csv'
p_chemlist_carcinogen = PR_ROOT + 'comptox/carcinogen_breast_122120_chemlist.csv'

## E2
p_ToxPrint_E2 = PR_ROOT + 'comptox/H295_E2-up_ToxPrint.csv'
p_chemlist_E2 = PR_ROOT + 'comptox/H295_E2-up_chemlist.csv'

## P4
p_ToxPrint_P4 = PR_ROOT + 'comptox/H295_P4-up_ToxPrint.csv'
p_chemlist_P4 = PR_ROOT + 'comptox/H295_P4-up_chemlist.csv'


# compute desc #
#################
# carcinogens - all
#computeDescAnalysisFromList(p_list_carcinogen, p_ToxPrint_carcinogen, p_chemlist_carcinogen, PR_RESULTS)

#E2
#computeDescAnalysisFromList(p_list_e2, p_ToxPrint_E2, p_chemlist_E2, PR_RESULTS)

#P4
#computeDescAnalysisFromList(p_list_p4, p_ToxPrint_P4, p_chemlist_P4, PR_RESULTS)


# ToxPi #
#########
# X2 for lists
cComparison = comparisonChemicalLists.comparisonChemicalLists([p_list_carcinogen, p_list_e2, p_list_p4], PR_RESULTS)
cComparison.X2Comparison()



