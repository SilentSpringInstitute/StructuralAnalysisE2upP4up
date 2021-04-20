from os import system, path, remove, chdir, getcwd, listdir, name
import subprocess 


P_RSCRIPTS = "../R/"
PR_BIOTRANSFORMER = "C:/Users/Aborrel/research/Silent_Spring/PFAS/BioTransformerJar/biotransformerjar"

R_BIN = "C:\\Program Files\\R\\R-4.0.2\\bin\\Rscript.exe"
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