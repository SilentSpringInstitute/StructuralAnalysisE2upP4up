from os import path


import Mcarcinogen
import pathFolder



# Define folder #
#################
PR_ROOT = path.abspath("../../") + "/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")
PR_DATA = PR_ROOT + "data/"

# val to prepare descriptors table
COR_VAL = 0.9
MAX_QUANTILE = 90


p_crossRef = PR_DATA + "Updated_MC_list_crossref_other_lists_020221.xlsx"

cMC = Mcarcinogen.Mcarcinogen(p_crossRef, COR_VAL, MAX_QUANTILE, PR_ROOT + "comptox/", PR_RESULTS)
cMC.main()

ss




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
#overlapBetweenList([p_list_carcinogen, p_list_e2, p_list_p4, p_list_no_carci, p_list_ER_agonist, p_list_ER_antagonist] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_all_lists_genotox/"))
#overlapBetweenList([p_list_carcinogen, p_list_e2, p_list_p4, p_list_no_carci, p_list_ER_agonist, p_list_ER_antagonist], pathFolder.createFolder(PR_RESULTS + "Overlap_all_lists/"))
#overlapBetweenList([p_list_carcinogen, p_list_ER_agonist, p_list_ER_antagonist], pathFolder.createFolder(PR_RESULTS + "Overlap_ER-carcinogen/"))
#overlapBetweenList([p_list_ER_agonist, p_list_ER_antagonist] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_ER-carcinogen-genotox/"))
#overlapBetweenList([p_list_e2, p_list_p4 ] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_P4-E2-carcinogen-genotox/"))
#overlapBetweenList([p_list_e2, p_list_p4, p_list_ER_agonist] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_P4-E2-ERago-carcinogen-genotox/"))
#overlapBetweenList([p_list_carcinogen] + l_p_genox, pathFolder.createFolder(PR_RESULTS + "Overlap_carcinogen_genotox/"))
#overlapBetweenList([p_list_carcinogen, p_list_e2, p_list_p4, p_list_ER_agonist, l_p_genox[1]], pathFolder.createFolder(PR_RESULTS + "Overlap_carcinogen-E2-P4-ERago-genotoxYes/"))
#overlapBetweenList([p_list_carcinogen, p_list_e2, p_list_p4, p_list_ER_agonist, l_p_genox[2]], pathFolder.createFolder(PR_RESULTS + "Overlap_carcinogen-E2-P4-ERago-genotoxNo/"))
#overlapBetweenList([p_list_carcinogen, p_list_e2, p_list_p4, p_list_ER_agonist, l_p_genox[2]], pathFolder.createFolder(PR_RESULTS + "Overlap_carcinogen-E2-P4-ERago-genotoxNo/"))

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

p_steroid = unionListChem([p_list_e2, p_list_p4], "steroid_up", PR_RESULTS)
overlapBetweenList([l_p_genox[1], p_list_ER_agonist, l_p_genox[2], p_steroid], pathFolder.createFolder(PR_RESULTS + "Overlap_carcinogen-steroid/"))
overlapBetweenList([l_p_genox[1], p_list_ER_agonist, l_p_genox[0], p_steroid], pathFolder.createFolder(PR_RESULTS + "Overlap_carcinogen-genotox-yes-notested-steroid/"))
overlapBetweenList([l_p_genox[1], p_list_ER_agonist, p_steroid, p_list_carcinogen], pathFolder.createFolder(PR_RESULTS + "Overlap_carcinogen-genotox-yes-steroid-ER-agonist/"))
sss

# chemical lists comparison # 
#############################

# X2 for lists - carcinogen and E2 - P4
cComparison = comparisonChemicalLists.comparisonChemicalLists([p_list_carcinogen, p_list_e2, p_list_p4], PR_RESULTS)
cComparison.X2Comparison()

# X2 for stereo
cComparison = comparisonChemicalLists.comparisonChemicalLists([p_list_carcinogen, p_list_e2, p_list_p4, p_list_ER_agonist, p_list_ER_antagonist], PR_RESULTS)
cComparison.X2Comparison()

# x2 only genotox
cComparison = comparisonChemicalLists.comparisonChemicalLists([p_list_carcinogen, l_p_genox[1], l_p_genox[2]], PR_RESULTS)
cComparison.X2Comparison()



