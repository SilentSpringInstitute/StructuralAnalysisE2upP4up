from os import path

import MCcrossref_data
import steroidogenesis_data
import merge_MCcrossWithStereo
import pathFolder
import buildQSAR

# Define folders #
##################
PR_ROOT = path.abspath("../../") + "/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")
PR_DATA = PR_ROOT + "data/"

# Define dataset #
##################
# cross ref excel including all dataset to process
p_listChem = PR_DATA + "cross_lists_for_analysis_090921.xlsx" # excel created with ehp sup and Jenny list 
# exposure by chemicals
p_exposure = PR_DATA + "BCRelExposureSources_P65_051221.csv"
# list of hormone in SMILES that can be considered for the 
p_hormones = PR_DATA + "hormones.csv"
# 270 chemicals as input


# Dataset preparation value #
#############################
COR_VAL = 0.9
MAX_QUANTILE = 90


# LOAD Steroidogenesis data from Karmaus2016 #
##############################################
c_Stereo = steroidogenesis_data.Steroidogenesis_data(PR_DATA, PR_RESULTS)
c_Stereo.main()


# LOAD AND RUN MC crossref ANALYSIS #
#####################################
c_MCcrossref = MCcrossref_data.MCcrossref(p_listChem, p_exposure, p_hormones, COR_VAL, MAX_QUANTILE, PR_ROOT + "comptox/", PR_ROOT)
c_MCcrossref.main()

# overlap with ToxCast - aromatase assays
##########################################


# MIXTE INFORMATION FROM STEROIDOGENESIS AND cross ref file #
#############################################################
#pr_out = pathFolder.createFolder(PR_RESULTS + "crossRef_Stereo/")
#c_MCcrossref_stereo = merge_MCcrossWithStereo.merge_MCcrossWithStereo(c_MCcrossref, c_Stereo, pr_out)
#c_MCcrossref_stereo.main()


## DEVELOP QSAR FOR E2up AND P4up ##
####################################

## E2up ##
##########

# no undersampling 
##############
#MAX_QUANTILE = 0
#name_QSAR = "QSAR_E2_H295R_nosampling_nosingledosecheck_noborderline"
#c_QSAR_E2up = buildQSAR.buildQSAR(name_QSAR, "E2up", "H295R", c_MCcrossref, PR_RESULTS, COR_VAL, MAX_QUANTILE)
#c_QSAR_E2up.buildDataset(c_Stereo, borderline=0)
#c_QSAR_E2up.buildDescSet(["rdkit", "OPERA", "toxprint"])
#c_QSAR_E2up.prepDesc()
#c_QSAR_E2up.computeSimMatrix()# similarity matrix for the AD
#c_QSAR_E2up.runQSARs()# no sampling add


# undersampling with variable active rate
#####################
#MAX_QUANTILE = 0
#name_QSAR = "QSAR_E2_H295R_variable-sampling_singledosecheck_noborderline"
#c_QSAR_E2up = buildQSAR.buildQSAR(name_QSAR, "E2up", "H295R", c_MCcrossref, PR_RESULTS, COR_VAL, MAX_QUANTILE)
#c_QSAR_E2up.buildDataset(c_Stereo, borderline=0)
#c_QSAR_E2up.buildDescSet(["rdkit", "OPERA", "toxprint"])
#c_QSAR_E2up.prepDesc()
#c_QSAR_E2up.computeSimMatrix()# similarity matrix for the AD
#c_QSAR_E2up.runQSARs([0.10, 0.9])



## P4up ##
##########

# no undersampling 
##############
MAX_QUANTILE = 0
name_QSAR = "QSAR_P4_H295R_nosampling_nosingledosecheck_noborderline"
c_QSAR_P4up = buildQSAR.buildQSAR(name_QSAR, "P4up", "H295R", c_MCcrossref, PR_RESULTS, COR_VAL, MAX_QUANTILE)
c_QSAR_P4up.buildDataset(c_Stereo, borderline=0)
c_QSAR_P4up.buildDescSet(["rdkit", "OPERA", "toxprint"])
c_QSAR_P4up.prepDesc()
c_QSAR_P4up.computeSimMatrix()# similarity matrix for the AD
c_QSAR_P4up.runQSARs()# no sampling add


# undersampling with variable active rate
#######################
MAX_QUANTILE = 0
name_QSAR = "QSAR_P4_H295R_variable-sampling_singledosecheck_noborderline"
c_QSAR_P4up = buildQSAR.buildQSAR(name_QSAR, "P4up", "H295R", c_MCcrossref, PR_RESULTS, COR_VAL, MAX_QUANTILE)
c_QSAR_P4up.buildDataset(c_Stereo, borderline=0)
c_QSAR_P4up.buildDescSet(["rdkit", "OPERA", "toxprint"])
c_QSAR_P4up.prepDesc()
c_QSAR_P4up.computeSimMatrix()# similarity matrix for the AD
c_QSAR_P4up.runQSARs([0.10, 0.9])


