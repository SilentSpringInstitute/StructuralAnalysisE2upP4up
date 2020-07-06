import pathFolder
import dataset
import analysis




# Define folder #
#################
PR_ROOT = "../../"
PR_RESULTS = PR_ROOT + "results/"
PR_DATA = PR_ROOT + "data/"



# Extract data and compute desc #
#################################
p_dataset = PR_DATA + "PFAS_45_updated.csv"

# load dataset
c_dataset = dataset.dataset(p_dataset, PR_RESULTS)
c_dataset.loadDataset()


# compute desc 2D
p_desc = c_dataset.computeStructuralDesc()

# compute PNG
#c_dataset.computePNG()

# compute OPERA
#p_desc_opera = c_dataset.computeOPERADesc()

# compute biotransformation 
#c_dataset.predictBiotransformation()

# extract biostransformation
c_dataset.extractBioTranformationChemical(p_check=PR_DATA + "Master_chemical_List.csv")
c_dataset.searchBiotransformedProduceInDB()

#c_dataset.computeDescProductBiotransformed()
ddd

# 2. analysis
#####
COR_VAL = 0.90
MAX_QUANTILE = 90

cAnalysis = analysis.analysis(p_desc, PR_RESULTS, COR_VAL, MAX_QUANTILE)
#cAnalysis.prepDesc()

# 2.1 PCA
#cAnalysis.PCA_plot()

# 2.2 Hclust
cAnalysis.HClust_plot(p_desc_opera)

# 2.3 Clustering 
cAnalysis.clustering()


# map on the PFAS map

wwww

# 2.3 SOM
size = 15
#cAnalysis.generate_SOM(15)
#cAnalysis.signifDescBySOMCluster()
#cAnalysis.extract_actBySOMCluster(pr_desc + "PNG/") # have to run !!!!




