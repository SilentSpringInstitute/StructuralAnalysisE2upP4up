from os import path

import MCcrossref_data
import steroidogenesis_data
import merge_MCcrossWithStereo
import pathFolder


# Define folders #
##################
PR_ROOT = path.abspath("../../") + "/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")
PR_DATA = PR_ROOT + "data/"

# Define dataset #
##################
p_crossRef = PR_DATA + "Updated_MC_list_crossref_other_lists_020221.xlsx" # need to update 
p_exposure = PR_DATA + "BCRelExposureSources_P65_051221.csv"
p_hormones = PR_DATA + "hormones.csv"
# 270 chemicals as input


# Dataset preparation value #
#############################
COR_VAL = 0.9
MAX_QUANTILE = 90


# LOAD Steroidogenesis data from Karmaus2016 #
##############################################
c_Stereo = steroidogenesis_data.Steroidogenesis_data(PR_DATA, PR_RESULTS)
#c_Stereo.main()


# LOAD AND RUN MC crossref ANALYSIS #
#####################################
c_MCcrossref = MCcrossref_data.MCcrossref(p_crossRef, p_exposure, p_hormones, COR_VAL, MAX_QUANTILE, PR_ROOT + "comptox/", PR_RESULTS)
c_MCcrossref.main()


# MIXTE INFORMATION FROM STEROIDOGENESIS AND MC #
#################################################
pr_out = pathFolder.createFolder(PR_RESULTS + "MC_stereo/")
c_MC_stereo = merge_MCcrossWithStereo.merge_MCcrossWithStereo(c_MCcrossref, c_Stereo, pr_out)
c_MC_stereo.main()


#cStereo.PCA_FoldChangeMC(cMC.d_MC)
#cStereo.cardTanimoto(cMC.d_MC)



