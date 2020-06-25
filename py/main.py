import pathFolder
import dataset





# Define folder #
#################
PR_ROOT = "../../"
PR_RESULTS = PR_ROOT + "results/"
PR_DATA = PR_ROOT + "data/"



# Extract data and compute desc #
#################################
p_dataset = PR_DATA + "PFAS_45.csv"

# load dataset
c_dataset = dataset.dataset(p_dataset, PR_RESULTS)
c_dataset.loadDataset()


# compute desc
p_desc = c_dataset.computeDesc()
#c_dataset.computePNG()

ss








# 2. analysis
#####
COR_VAL = 0.90
MAX_QUANTILE = 90

#cAnalysis = analysis.analysis(p_AC50, p_desc, PR_RESULTS, COR_VAL, MAX_QUANTILE)
#cAnalysis.prepDesc()

# 2.1. histogram AC50 and summary
#cAnalysis.sumAC50()

# 2.2 PCA
#cAnalysis.PCA_plot()

# 2.3 SOM
size = 15
#cAnalysis.generate_SOM(15)
#cAnalysis.signifDescBySOMCluster()
#cAnalysis.extract_actBySOMCluster(pr_desc + "PNG/") # have to run !!!!


# 2.4 Hclust
#cAnalysis.HClust_plot()

