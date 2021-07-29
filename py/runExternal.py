from os import system, path, remove, chdir, getcwd, listdir, name
import subprocess 


P_RSCRIPTS = "../R/"
PR_BIOTRANSFORMER = "C:/Users/Aborrel/research/Silent_Spring/PFAS/BioTransformerJar/biotransformerjar"

R_BIN = "C:\\Program Files\\R\\R-4.0.2\\bin\\Rscript.exe"
P_RQSAR_linux = "/mnt/c/Users/AlexandreBorrel/research/development/QSAR-QSPR"
P_RQSAR_window = "c:/Users/aborr/research/development/QSAR-QSPR/"
######
# Main functions

def runRCMD(cmd, out = 0):

    workdir = getcwd()
    chdir(P_RSCRIPTS)
    if name == "nt":
        l_elem = cmd.split(" ")
        cmd_line = [R_BIN] + l_elem
        print(cmd_line)
        p = subprocess.Popen(cmd_line)
        (output, err) = p.communicate() 
        p.wait()
        print(err)
    else:
        print(cmd)
        system(cmd)
    chdir(workdir)


def runRQSARModeling(cmd):
    """Run external R scripts with QSAR modeling functions"""

    workdir = getcwd()
    
    if name == "nt":
        chdir(P_RQSAR_window)
        #cmd = cmd.replace("/", "\\")
        l_elem = cmd.split(" ")
        cmd_line = [R_BIN] + l_elem
        print(cmd_line)
        p = subprocess.Popen(cmd_line)
        (output, err) = p.communicate() 
        p.wait()
        print(err)
    else:
        chdir(P_RQSAR_linux)
        print(cmd)
        system(cmd)
    chdir(workdir)


############
# Functions for analysis


def preprocData(p_desc, pr_out, cor_val, max_quantile):
    cmd = "./cleanerData.R %s %s %s %s"%(p_desc, pr_out, cor_val, max_quantile)
    runRCMD(cmd)


def PCA(p_desc_cleaned, pr_out):
    cmd = "./PCA_chem.R %s %s"%(p_desc_cleaned, pr_out)
    runRCMD(cmd)


def HClust(p_desc_cleaned, p_dataset, p_opera, pr_out):
    cmd = "./HClust_chem.R %s %s %s %s"%(p_desc_cleaned, p_dataset, p_opera, pr_out)
    runRCMD(cmd)


def Clust(p_desc_cleaned, p_desc_opera, pr_out):
    cmd = "./ClustAlgo.R %s %s %s"%(p_desc_cleaned, p_desc_opera, pr_out)
    runRCMD(cmd)


def projectToSpace(p_listchem, p_1D2D, p_3D, name_map, pr_out):
    cmd = "./densityMapOverlap.R %s %s %s %s %s"%(p_1D2D, p_3D, p_listchem, name_map, pr_out)
    runRCMD(cmd)


def drawHist(p_desc, p_desc_clean, pr_out):

    cmd = "./drawHist.R %s %s %s"%(p_desc, p_desc_clean, pr_out)
    runRCMD(cmd)


def cardSimMatrix(p_filin):

    cmd = "./cardFP.R %s"%(p_filin)
    runRCMD(cmd)
    return 


def SOM(p_desc_cleaned, p_AC50_cleaned, pr_out, grid_size):

    cmd = "./SOM_chem.R %s %s %s %s"%(p_desc_cleaned, p_AC50_cleaned, pr_out, grid_size)
    runRCMD(cmd)


def descSignifByCluster(p_desc, p_cluster, pr_out):

    cmd = "./descSignificantByCluster.R %s %s %s"%(p_desc, p_cluster, pr_out)
    runRCMD(cmd)

def barplotToxPrint(p_filin):
    
    cmd = "./barplotToxPrint.R %s"%(p_filin)
    runRCMD(cmd)


def barplotchemlist(p_filin):
    
    cmd = "./barplotChemList.R %s"%(p_filin)
    runRCMD(cmd)


def plotX2(p_filin):

    cmd = "./comparisonX2.R %s"%(p_filin)
    runRCMD(cmd)

def upsetPlot(p_filin):

    cmd = "./upsetplot.R %s"%(p_filin)
    runRCMD(cmd)

def vennPlot(p_filin):
    cmd = "./vennPlot.R %s"%(p_filin)
    runRCMD(cmd)

def multiSetsAnalysis(p_filin):
    cmd = "./multi_sets_analysis.R %s"%(p_filin)
    runRCMD(cmd)

def barplotHormones(p_filin):
    cmd = "./barplotHormones.R %s"%(p_filin)
    runRCMD(cmd)   

def PCA_SteroiMC(p_matrix_single_hitc):
    cmd = "./PCA_SteroiMC.R %s"%(p_matrix_single_hitc)
    runRCMD(cmd)

def FP_card(p_filin):
    cmd = "./FP_card.R %s"%(p_filin)
    runRCMD(cmd)

def dendogramClusterProp(p_dataset, p_desc, p_hormone, pr_out, cor_val, max_q):

    cmd =  "./dendogramProp.R %s %s %s %s %s %s"%(p_dataset, p_desc, p_hormone, pr_out, cor_val, max_q)
    runRCMD(cmd)

def dendogramClusterTwoProp(p_dataset, p_stereo, p_desc, pr_out, cor_val, max_q):

    cmd =  "./dendogramTowProp.R %s %s %s %s %s %s"%(p_dataset, p_stereo, p_desc, pr_out, cor_val, max_q)
    runRCMD(cmd)

def projectPropInSOM(p_SOM, p_all, pr_out):

    cmd = "./SOM_mapprop.R %s %s %s"%(p_SOM, p_all, pr_out)
    runRCMD(cmd)

def dendogramFPProp(p_dataset, p_FP, pr_out):
    
    cmd =  "./dendogramPropFPTanimoto.R %s %s %s"%(p_dataset, p_FP, pr_out)
    runRCMD(cmd)

def barplotChemClass(p_count):
    cmd = "./barplotChemClass.R %s"%(p_count)
    runRCMD(cmd)   

def enrichmentByCluster(p_prop, p_clusters, pr_out):
    cmd = "./enrichmentClusters.R %s %s %s"%(p_prop, p_clusters, pr_out)
    runRCMD(cmd)

def cardPlotSimilarity(p_filin):
    cmd = "./cardPlotFP.R %s"%(p_filin)
    runRCMD(cmd)

def comparisonDesc(p_desc1, p_desc2, pr_out):
    cmd = "./signiDesc.R %s %s %s"%(p_desc1, p_desc2, pr_out)
    runRCMD(cmd)

def SOM_hormone(p_model_SOM, p_hormone_similarity, pr_out):
    cmd = "./SOM_hormoneSimilarity.R %s %s %s"%(p_model_SOM, p_hormone_similarity, pr_out)
    runRCMD(cmd)


def comparisonWithHormoneSimilarity(p_dataset1, p_dataset2, p_hormone, pr_out):

    cmd = "./comparisonHormoneSimilarity.R %s %s %s %s"%(p_dataset1, p_dataset2, p_hormone, pr_out)
    runRCMD(cmd)


def overlapListHormoneSim(p_upset, p_hormone_similarity, pr_out):

    cmd = "./overlapListHormoneSim.R %s %s %s"%(p_upset, p_hormone_similarity, pr_out)
    runRCMD(cmd)


def overlapAssaysListChem(p_assays, p_upset, pr_out):
    cmd = "./overlapListAndAssays.R %s %s %s"%(p_assays, p_upset, pr_out)
    runRCMD(cmd)


def comparisonToxPrint(p_toxprint1, p_toxprint2, pr_out):
    cmd = "./comparisonToxPrint.R %s %s %s"%(p_toxprint1, p_toxprint2, pr_out)
    runRCMD(cmd)



##############
# run biotransformer tool

def BioTransformer(smi, btType, p_out, nsteps=1):

    p_out = path.abspath(p_out)
    workdir = getcwd()
    chdir(PR_BIOTRANSFORMER)
    cmd = "java.exe -jar ./biotransformer-1.1.5.jar -k pred -b %s -ismi \"%s\" -ocsv %s -s %s"%(btType, smi, p_out, nsteps)
    print("********")
    print(cmd)
    print("***********")
    system(cmd)
    chdir(workdir)


######################
### QSAR modeling 

def combineDesc(p_desc_rdkit, p_desc_opera, p_toxPrint, pr_out):

    cmd = "./mergerDescData.R %s %s %s %s"%(p_desc_rdkit, p_desc_opera, p_toxPrint, pr_out)
    runRCMD(cmd)


def preprocDataQSAR(p_desc, p_AC50, pr_out, cor_val, max_quantile, act=0):

    cmd = "./cleanerDataQSAR.R %s %s %s %s %s %s"%(p_desc, p_AC50, pr_out, cor_val, max_quantile, act)
    runRCMD(cmd)



############
# Function for QSAR

def SplitTrainTest(pdescAc50, prout, splitratio):

    pdescAc50 = path.abspath(pdescAc50)
    prout = path.abspath(prout) + "/"

    cmd = "./prepTrainTestSplitClass.R " + pdescAc50 + " " + str(splitratio) + " " + prout
    runRQSARModeling(cmd)


def prepDataQSARReg(pdesc, paff, prout, corcoef, maxQuantile, valSplit, typeAff="All", logaff=0, nbNA = 10):

    cmd = "./QSARsPrep.R " + str(pdesc) + " " + str(paff) + " " + prout + " " + str(corcoef) + " " + str(
        maxQuantile) + " " + str(valSplit) + " " + str(logaff) + " " + str(typeAff) + " " + str(nbNA)
    runRQSARModeling(cmd) 

def prepDataQSAR(p_desc, p_AC50, rate_active, pr_run):

    p_desc = path.abspath(p_desc)
    p_AC50 = path.abspath(p_AC50)
    pr_run = path.abspath(pr_run) + "/"

    cmd = "./prepClassDataset.R %s %s %s %s"%(p_desc, p_AC50, rate_active, pr_run)
    runRQSARModeling(cmd)


def runRQSAR(p_train, p_test, n_foldCV, pr_run):

    p_train = path.abspath(p_train)
    p_test = path.abspath(p_test)
    pr_run = path.abspath(pr_run) + "/"

    cmd = "./QSARsClass.R " + p_train + " " + p_test + " 0 " + pr_run + " " + str(n_foldCV) + " > " + pr_run + "perf.txt"
    runRQSARModeling(cmd)


def runQSARReg(ptrain, ptest, pcluster, prout, nbfold=10):

    cmd_QSAR = "./QSARsReg.R " + ptrain + " " + ptest + " " + pcluster + " " + prout + " " + str(nbfold) + " 1 >" + prout + "perf.txt"
    runRQSARModeling(cmd_QSAR)
    

def plotAC50VSProb(p_prob):

    cmd = "./plotAC50vsProb.R %s"%(p_prob)
    runRCMD(cmd)

def runImportanceDesc(p_desc, nb):
    
    p_desc = path.abspath(p_desc)

    cmd = "./importancePlot.R " + str(p_desc) + " " + str(nb)
    runRQSARModeling(cmd)


def AD(p_desc_model, p_desc_test, pr_out):

    p_desc_model = path.abspath(p_desc_model)
    p_desc_test = path.abspath(p_desc_test)
    pr_out = path.abspath(pr_out) + "/"

    cmd = "./computeAD.R %s %s %s"%(p_desc_model, p_desc_test, pr_out)
    runRCMD(cmd)


def mergeADs(p_train, p_test, p_desc, pr_out):

    cmd = "./mergeADs.R %s %s %s %s"%(p_train, p_test, p_desc, pr_out)
    runRCMD(cmd)

#def predictDataset(p_desc_test, p_model, ML,  pr_out):

 #   cmd = "./predictTestSet.R %s %s %s %s"%(p_desc_test, p_model, ML, pr_out)
 #   runRCMD(cmd)
