import pathFolder
import dataset
import analysis
from os import path
import toolbox
import mapOnSpace


# selected physico chemical descriptor from OPERa
L_OPERA_DESC = ['LogP_pred', 'MP_pred', 'BP_pred', 'LogVP_pred', 'LogWS_pred', 'LogHL_pred', 'RT_pred', 'LogKOA_pred', 'ionization', 'LogD55_pred', 'LogD74_pred', 'LogOH_pred', 'LogBCF_pred', 'BioDeg_LogHalfLife_pred', 'ReadyBiodeg_pred', 'LogKM_pred', 'LogKoc_pred', 'FUB_pred', 'Clint_pred']


##############
# main functions 
##############

def computeDesc(p_dataset, pr_result):
    #Return list with 1D 2D desc and desc from OPERA


    # load dataset
    c_dataset = dataset.dataset(p_dataset, pr_result)
    c_dataset.loadDataset(loadDb=1)


    # compute desc 2D
    p_desc = c_dataset.computeStructuralDesc()

    # compute PNG
    c_dataset.computePNG()

    # compute OPERA
    p_desc_opera = c_dataset.computeOPERADesc()

    return [p_desc, p_desc_opera]
    



def computeBiotransformation(c_dataset, p_list_mater=""):

    # compute biotransformation 
    c_dataset.predictBiotransformation()

    # extract biostransformation
    c_dataset.extractBioTranformationChemical(p_list_mater)
    c_dataset.searchBiotransformedProduceInDB()

    c_dataset.computeDescProductBiotransformed()




def analysisDescBasedData(p_desc, p_desc_opera, pr_results, cor_val, max_quantile):


    cAnalysis = analysis.analysis(p_desc, pr_results, cor_val, max_quantile)
    cAnalysis.prepDesc()

    # 2.1 PCA
    cAnalysis.PCA_plot()

    # 2.2 Hclust
    cAnalysis.HClust_plot(p_desc_opera)

    # 2.3 Clustering 
    cAnalysis.clustering()

    # 2.4 map on PFAS map
    cAnalysis.mapOnSpace()

    # map on the PFAS map

    # 2.3 SOM
    #size = 15
    #cAnalysis.generate_SOM(15)
    #cAnalysis.signifDescBySOMCluster()
    #cAnalysis.extract_actBySOMCluster(pr_desc + "PNG/") # have to run !!!!


    return 




def buildDescSet(p_rdkit, p_opera, l_type_desc, pr_results):

    """Select from OPERA only physico chem descriptors"""

    p_filout = pr_results + "-".join(l_type_desc) + ".csv"
    if path.exists(p_filout):
        return p_filout

    if len(l_type_desc) == 1 and l_type_desc[0] == "rdkit":
        return p_rdkit

    # open RDKIT use for CARSN and SMILES
    d_rdkit = toolbox.loadMatrix(p_rdkit, sep="\t")

    d_out = {}
    l_header = ["CASRN", "SMILES"]
    if "opera" in l_type_desc:
        d_OPERA = toolbox.loadMatrix(p_opera, sep = ",")
        l_header = l_header + L_OPERA_DESC

        for chem in d_OPERA:
            if not chem in list(d_out.keys()):
                d_out[chem] = {}
                d_out[chem]["CASRN"] = d_rdkit[chem]["CASRN"]
                d_out[chem]["SMILES"] = d_rdkit[chem]["SMILES"]
            for h in L_OPERA_DESC:
                d_out[chem][h] = d_OPERA[chem][h]



    if "rdkit" in l_type_desc:
        l_h_drdkit = list(d_rdkit[list(d_rdkit.keys())[0]].keys())
        l_h_drdkit.remove("CASRN")
        l_h_drdkit.remove("SMILES")
        l_header = l_header + l_h_drdkit

        for chem in d_rdkit.keys():
            if not chem in d_out.keys():
                d_out[chem] = {}
                d_out[chem]["CASRN"] = d_rdkit[chem]["CASRN"]
                d_out[chem]["SMILES"] = d_rdkit[chem]["SMILES"]
            
            for h in l_h_drdkit:
                d_out[chem][h] = d_rdkit[chem][h]

    filout = open(p_filout, "w")
    filout.write("\t".join(l_header) + "\n")
    for chem in d_out.keys():
        filout.write("\t".join([d_out[chem][h] for h in l_header]) + "\n")
    filout.close()


    return p_filout


###############
# main run 
###############



# Define folder #
#################
PR_ROOT = "../../"
PR_RESULTS = PR_ROOT + "results/"
PR_DATA = PR_ROOT + "data/"

COR_VAL = 0.90
MAX_QUANTILE = 90

P_COORD_1D2D = PR_DATA + "coord1D2D.csv"
P_COORD_3D = PR_DATA + "coord3D.csv"

print(PR_DATA)

# Extract data and compute desc #
#################################

########## Phthalates
#####################
p_dataset = PR_DATA + "Phthalates.csv"
pr_results = pathFolder.createFolder(PR_RESULTS + "Phthalates/")
l_file_desc = computeDesc(p_dataset, pr_results)

## analysis with RDKIT desc
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT_desc/")
#analysisDescBasedData(l_file_desc[0], l_file_desc[1], pr_results, COR_VAL, MAX_QUANTILE)

## analysis with RDKIT + OPERA physico chemical
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT-OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit", "opera"], pr_results_analysis)
#analysisDescBasedData(p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

## analysis with OPERA physico chemical only
#pr_results_analysis = pathFolder.createFolder(pr_results + "OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["opera"], pr_results_analysis)
#analysisDescBasedData(p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

### map on space
#####
#pr_results_mapped = pathFolder.createFolder(pr_results + "mapped/")
#c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D, P_COORD_3D, pr_results_mapped)
#c_mapped.map("Phthalates")



########## Phthalates alternative
#################################
p_dataset = PR_DATA + "Phthalates_alternatives.csv"
pr_results = pathFolder.createFolder(PR_RESULTS + "Phthalates_alternatives/")
l_file_desc = computeDesc(p_dataset, pr_results)

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


### map on space
#####
pr_results_mapped = pathFolder.createFolder(pr_results + "mapped/")
c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D, P_COORD_3D, pr_results_mapped)
c_mapped.map("Phthalates_alternatives")


############### PFAS
#####################
p_dataset = PR_DATA + "PFAS.csv"
pr_results = pathFolder.createFolder(PR_RESULTS + "PFAS/")
l_file_desc = computeDesc(p_dataset, pr_results)

# analysis with RDKIT desc
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit"], pr_results_analysis)
#analysisDescBasedData(l_file_desc[0], l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

# analysis with RDKIT + OPERA physico chemical
#pr_results_analysis = pathFolder.createFolder(pr_results + "RDKIT-OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["rdkit", "opera"], pr_results_analysis)
#analysisDescBasedData(p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)

# analysis with OPERA physico chemical only
#pr_results_analysis = pathFolder.createFolder(pr_results + "OPERA_desc/")
#p_desc = buildDescSet(l_file_desc[0], l_file_desc[1], ["opera"], pr_results_analysis)
#analysisDescBasedData(p_desc, l_file_desc[1], pr_results_analysis, COR_VAL, MAX_QUANTILE)


pr_results_mapped = pathFolder.createFolder(pr_results + "mapped/")
c_mapped = mapOnSpace.mapOnSpace(l_file_desc[0], P_COORD_1D2D, P_COORD_3D, pr_results_mapped)
c_mapped.map("PFAS_NTP")
