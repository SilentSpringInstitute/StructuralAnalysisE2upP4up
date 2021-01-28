from os import path
from copy import deepcopy

import pathFolder
import dataset
import analysis
import Chemicals
import comparisonChemicalLists
import toolbox
import runExternal


# Define folder #
#################
PR_ROOT = path.abspath("../../") + "/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")
PR_DATA = PR_ROOT + "data/"

# val
COR_VAL = 0.9
MAX_QUANTILE = 90

def computeDescAnalysisFromList(p_list, p_ToxPrint, p_chemList, PR_RESULTS):

    pr_results = pathFolder.createFolder(PR_RESULTS + "analysis_individual-dataset/" + p_list.split("/")[-1][0:-4] + "/")
    pr_desc = pathFolder.createFolder(PR_RESULTS + "DESC/")

    # Compute desc
    ###################
    cChem = Chemicals.Chemicals(p_list, pr_results)
    cChem.computeDesc(pr_desc)
    #cChem.buildDescSet(["rdkit"])
    #cChem.analysisDescBasedData(COR_VAL, MAX_QUANTILE, PCA=1, Hclust=1, clustering=1, FP=1, SOM=1, histDesc=1)

    cChem.analysisToxPrint(p_ToxPrint)
    cChem.analysisChemList(p_chemList)

def overlapBetweenList(l_list_chem, pr_results):

    
    d_d_chem = {}
    l_CASRN = []
    for p_list_chem in l_list_chem:
        name_list = p_list_chem.split("/")[-1][0:-4]
        d_chem = toolbox.loadMatrix(p_list_chem, sep = ",")
        l_CASRN_chem = []
        for chem in d_chem.keys():
            try:
                l_CASRN_chem.append(d_chem[chem]["CASRN"])
                l_CASRN.append(d_chem[chem]["CASRN"])
            except:continue
        d_d_chem[name_list] = l_CASRN_chem
    
    p_upset = pr_results + "upset_matrix"
    f_open = open(p_upset, "w")
    f_open.write("\t" + "\t".join(list(d_d_chem.keys())) + "\n")
    l_CASRN = list(set(l_CASRN))
    for CASRN in l_CASRN:
        l_w = []
        for list_chem in d_d_chem:
            if CASRN in d_d_chem[list_chem]:
                l_w.append("1")
            else:
                l_w.append("0")
        f_open.write("%s\t%s\n"%(CASRN, "\t".join(l_w)))
    f_open.close()

    runExternal.upsetPlot(p_upset)

    if len(l_list_chem) <= 4:
        runExternal.vennPlot(p_upset)


    return 

def splitgenotox(p_list, PR_RESULTS):

    pr_results = pathFolder.createFolder(PR_RESULTS + "analysis_individual-dataset/" + p_list.split("/")[-1][0:-4] + "/")

    c_dataset = dataset.dataset(p_list, pr_results)
    l_p_dataset = c_dataset.splitDataset("Genotoxic_CCRIS/QSAR/ToxVal")
    return l_p_dataset

def mergedataset(l_psetchems, PR_RESULTS):
    
    pr_desc = pathFolder.createFolder(PR_RESULTS + "DESC/")

    
    l_dataset = []
    d_desc_1D2D = {}
    d_desc_opera = {}
    for psetchem in l_psetchems:
        name_dataset =  psetchem.split("/")[-1][0:-4]
        pr_results = pathFolder.createFolder(PR_RESULTS + "analysis_individual-dataset/" + name_dataset + "/")       
        l_dataset.append(name_dataset)

        # Compute desc
        ###################
        cChem = Chemicals.Chemicals(psetchem, pr_results)
        cChem.computeDesc(pr_desc)

        d_desc_1D2D_chem = toolbox.loadMatrix(cChem.p_desc_rdkit)
        d_desc_opera_chem = toolbox.loadMatrix(cChem.p_desc_OPERA, sep = ",")
        

        for chem in d_desc_1D2D_chem.keys():
            if not chem in list(d_desc_1D2D.keys()):
                d_desc_1D2D[chem] = deepcopy(d_desc_1D2D_chem[chem])
                d_desc_1D2D[chem]["dataset"] = []
            d_desc_1D2D[chem]["dataset"].append(name_dataset)
        
        for chem in d_desc_opera_chem.keys():
            if not chem in list(d_desc_opera.keys()):
                d_desc_opera[chem] = deepcopy(d_desc_opera_chem[chem])
                d_desc_opera[chem]["dataset"] = []
            d_desc_opera[chem]["dataset"].append(name_dataset)
        
    pr_out = pathFolder.createFolder(PR_RESULTS + "-".join(l_dataset) + "/")
    p_desc2D = pr_out + "desc1D2D.csv"
    l_h = list(d_desc_1D2D[list(d_desc_1D2D.keys())[0]].keys()) 
    l_h.remove("CASRN")
    f_desc2D = open(p_desc2D, "w")
    f_desc2D.write("CASRN\t" + "\t".join(l_h) + "\n")
    for chem in d_desc_1D2D.keys():
        f_desc2D.write(chem)
        for h in l_h:
            if h == "dataset":
                f_desc2D.write("\t%s"%("--".join(d_desc_1D2D[chem][h])))
            else:
                f_desc2D.write("\t%s"%(d_desc_1D2D[chem][h]))
        f_desc2D.write("\n")
    f_desc2D.close()

    return [p_desc2D]

def analysisMultiSets(p_desc1D2D):

    # draw hclust circular and PCA
    runExternal.multiSetsAnalysis(p_desc1D2D)



# list of chemicals available
#############
p_list_carcinogen = PR_DATA + "carcinogen_breast_012721.csv"
p_list_e2 = PR_DATA + "H295_E2-up.csv"
p_list_p4 = PR_DATA + "H295_P4-up.csv"
p_list_no_carci = PR_DATA + "nonCarci_list.csv"
p_list_ER_antagonist = PR_DATA + "Judson2015_ER_antagonist.csv"
p_list_ER_agonist = PR_DATA + "Judson2015_ER_agonist.csv"

# Explore ToxPrint - carcinogen
## carcinogen
p_ToxPrint_carcinogen = PR_ROOT + 'comptox/carcinogen_breast_012721_ToxPrint.csv'
p_chemlist_carcinogen = PR_ROOT + 'comptox/carcinogen_breast_012721_chemlist.csv' # need to update

## E2
p_ToxPrint_E2 = PR_ROOT + 'comptox/H295_E2-up_ToxPrint.csv'
p_chemlist_E2 = PR_ROOT + 'comptox/H295_E2-up_chemlist.csv'

## P4
p_ToxPrint_P4 = PR_ROOT + 'comptox/H295_P4-up_ToxPrint.csv'
p_chemlist_P4 = PR_ROOT + 'comptox/H295_P4-up_chemlist.csv'

## no carci
p_ToxPrint_NoCarci = PR_ROOT + 'comptox/NonCarci-up_ToxPrint.csv'
p_chemlist_NoCarci = PR_ROOT + 'comptox/NonCarci-up_chemlist.csv'

## ER agonist
p_ToxPrint_ER_agnonist = PR_ROOT + 'comptox/Judson2015-ER-agonist_ToxPrint.csv'
p_chemlist_ER_agnonist = PR_ROOT + 'comptox/Judson2015-ER-agonist_chemlist.csv'

## ER antagonist
p_ToxPrint_ER_antagnonist = PR_ROOT + 'comptox/Judson2015-ER-agonist_ToxPrint.csv'
p_chemlist_ER_antagnonist = PR_ROOT + 'comptox/Judson2015-ER-antagonist_chemlist.csv'

# compute desc #
#################
# carcinogens - all
#computeDescAnalysisFromList(p_list_carcinogen, p_ToxPrint_carcinogen, p_chemlist_carcinogen, PR_RESULTS)

# split genotox
l_p_genox = splitgenotox(p_list_carcinogen, PR_RESULTS)
#for p_genox in l_p_genox:
#    computeDescAnalysisFromList(p_genox, p_ToxPrint_carcinogen, p_chemlist_carcinogen, PR_RESULTS)


# overlap between set of chemicals
#####################################
overlapBetweenList([p_list_carcinogen, p_list_e2, p_list_p4, p_list_no_carci, p_list_ER_agonist, p_list_ER_antagonist] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_all_lists_genotox/"))
overlapBetweenList([p_list_carcinogen, p_list_e2, p_list_p4, p_list_no_carci, p_list_ER_agonist, p_list_ER_antagonist], pathFolder.createFolder(PR_RESULTS + "Overlap_all_lists/"))
overlapBetweenList([p_list_carcinogen, p_list_ER_agonist, p_list_ER_antagonist], pathFolder.createFolder(PR_RESULTS + "Overlap_ER-carcinogen/"))
overlapBetweenList([p_list_ER_agonist, p_list_ER_antagonist] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_ER-carcinogen-genotox/"))
overlapBetweenList([p_list_e2, p_list_p4 ] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_P4-E2-carcinogen-genotox/"))
overlapBetweenList([p_list_e2, p_list_p4, p_list_ER_agonist] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_P4-E2-ERago-carcinogen-genotox/"))
overlapBetweenList([p_list_carcinogen] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_carcinogen_genotox/"))


#E2
#computeDescAnalysisFromList(p_list_e2, p_ToxPrint_E2, p_chemlist_E2, PR_RESULTS)

#P4
#computeDescAnalysisFromList(p_list_p4, p_ToxPrint_P4, p_chemlist_P4, PR_RESULTS)

#ER agnonist
#computeDescAnalysisFromList(p_list_ER_agonist, p_ToxPrint_ER_agnonist, p_chemlist_ER_agnonist, PR_RESULTS)

#ER antagonist
#computeDescAnalysisFromList(p_list_ER_antagonist, p_ToxPrint_ER_antagnonist, p_chemlist_ER_antagnonist, PR_RESULTS)

#No carcinogen
#computeDescAnalysisFromList(p_list_no_carci, p_ToxPrint_NoCarci, p_chemlist_NoCarci, PR_RESULTS)


# Merge dataset
#################
#l_pdesc = mergedataset([p_list_carcinogen, p_list_e2, p_list_p4, p_list_ER_agonist, p_list_ER_antagonist], PR_RESULTS)
#analysisMultiSets(l_pdesc[0])

# ToxPi comparison# ==> to rewrite in part
###################
# X2 for lists - carcinogen and E2 - P4
#cComparison = comparisonChemicalLists.comparisonChemicalLists([p_list_carcinogen, p_list_e2, p_list_p4], PR_RESULTS)
#cComparison.X2Comparison()


# X2 for all of lists
#cComparison = comparisonChemicalLists.comparisonChemicalLists([p_list_carcinogen, p_list_e2, p_list_p4, p_list_no_carci, p_list_ER_agonist, p_list_ER_antagonist], PR_RESULTS)
#cComparison.X2Comparison()

# all with division genotox
#cComparison = comparisonChemicalLists.comparisonChemicalLists([p_list_carcinogen, p_list_e2, p_list_p4, p_list_no_carci, p_list_ER_agonist, p_list_ER_antagonist] + l_p_genox, PR_RESULTS)
#cComparison.X2Comparison()
