import pathFolder
import dataset
import analysis



##############
# main functions 
##############

def computeDesc(p_dataset, pr_result):
    """
    Return list with 1D 2D desc and desc from OPERA
    """


    # load dataset
    c_dataset = dataset.dataset(p_dataset, pr_result)
    c_dataset.loadDataset(loadDb=1)


    # compute desc 2D
    p_desc = c_dataset.computeStructuralDesc()

    # compute PNG
    c_dataset.computePNG()

    # compute OPERA
    p_desc_opera = c_dataset.computeOPERADesc()

    



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


    # map on the PFAS map

    # 2.3 SOM
    #size = 15
    #cAnalysis.generate_SOM(15)
    #cAnalysis.signifDescBySOMCluster()
    #cAnalysis.extract_actBySOMCluster(pr_desc + "PNG/") # have to run !!!!


    return 


###############
# main run 
###############



# Define folder #
#################
PR_ROOT = "../../"
PR_RESULTS = PR_ROOT + "results/"
PR_DATA = PR_ROOT + "data/"


# Extract data and compute desc #
#################################
p_dataset = PR_DATA + "PFAS_45_updated.csv"
c_dataset = computeDesc(p_dataset, PR_RESULTS)
computeBiotransformation(c_dataset, PR_DATA + "Master_chemical_List.csv")

dddw

# 2. analysis
#####
COR_VAL = 0.90
MAX_QUANTILE = 90
analysisDescBasedData(c_dataset.p_desc, c_dataset.p_desc_opera, PR_RESULTS, COR_VAL, MAX_QUANTILE)



