from os import path
import MCcrossref
import steroidogenesis
import pathFolder



# Define folder #
#################
PR_ROOT = path.abspath("../../") + "/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")
PR_DATA = PR_ROOT + "data/"

# Define dataset #
##################
p_crossRef = PR_DATA + "Updated_MC_list_crossref_other_lists_020221.xlsx" # need to update 
# 270 chemicals as input


# Dataset preparation value #
#############################
COR_VAL = 0.9
MAX_QUANTILE = 90


# LOAD Steroidogenesis data from Karmaus2016 #
##############################################

#cStereo = steroidogenesis.Steroidogenesis(PR_DATA, PR_RESULTS)
#cStereo.main()

# LOAD AND RUN MC crossref ANALYSIS #
############################

c_MCcrossref = MCcrossref.MCcrossref(p_crossRef, COR_VAL, MAX_QUANTILE, PR_ROOT + "comptox/", PR_RESULTS)
c_MCcrossref.main()





# MIXTE INFORMATION FROM STEROIDOGENESIS AND MC #
#################################################
#cStereo.PCA_FoldChangeMC(cMC.d_MC)
#cStereo.cardTanimoto(cMC.d_MC)



