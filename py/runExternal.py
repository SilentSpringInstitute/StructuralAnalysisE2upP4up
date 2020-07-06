from os import system, path, remove, chdir, getcwd, listdir
import subprocess 


P_RSCRIPTS = "../R/"
PR_BIOTRANSFORMER = "C:/Users/Aborrel/research/Silent_Spring/PFAS/BioTransformerJar/biotransformerjar"

R_BIN = "& 'C:\\Program Files\\R\\R-3.6.2\\bin\\Rscript.exe'"
######
# Main functions

def runRCMD(cmd, out = 0):

    workdir = getcwd()
    chdir(P_RSCRIPTS)
    print(R_BIN + " " + cmd)
    if out == 0:
        system(cmd)
        output = 0
    else:
        import subprocess
        output = subprocess.check_output(cmd, shell=True)
    chdir(workdir)
    return output




############
# Functions for analysis


def preprocData(p_desc, pr_out, cor_val, max_quantile):
    cmd = "./cleanerData.R %s %s %s %s"%(p_desc, pr_out, cor_val, max_quantile)
    runRCMD(cmd)


def PCA(p_desc_cleaned, pr_out):
    cmd = "./PCA_chem.R %s %s"%(p_desc_cleaned, pr_out)
    runRCMD(cmd)


def HClust(p_desc_cleaned, p_opera, pr_out):
    cmd = "./HClust_chem.R %s %s %s"%(p_desc_cleaned, p_opera, pr_out)
    runRCMD(cmd)


def Clust(p_desc_cleaned, pr_out):
    cmd = "./ClustAlgo.R %s %s"%(p_desc_cleaned, pr_out)
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