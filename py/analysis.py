from os import path
from statistics import mean, stdev
from shutil import copyfile

import pathFolder
import runExternal
import toolbox



class analysis:
    def __init__(self, p_dataset, p_desc, pr_out, cor_val, max_quantile):
        self.p_desc = p_desc
        self.p_dataset = p_dataset
        self.pr_out = pr_out
        self.cor_val = cor_val
        self.max_quantile = max_quantile



    def prepDesc(self):
        pr_out = pathFolder.createFolder(self.pr_out + "Cleaned_Data/")

        # add shortcut
        p_desc_cleaned = pr_out + "desc1D2D_cleaned.csv"
        if not path.exists(p_desc_cleaned):
            runExternal.preprocData(self.p_desc, pr_out, self.cor_val, self.max_quantile)
        self.p_desc_cleaned = p_desc_cleaned


    def PCA_plot(self):

        pr_out = pathFolder.createFolder(self.pr_out + "PCA/")
        runExternal.PCA(self.p_desc_cleaned, pr_out)



    def HClust_plot(self, p_opera):

        pr_out = pathFolder.createFolder(self.pr_out + "HClustCircular/")
        runExternal.HClust(self.p_desc_cleaned, self.p_dataset, p_opera, pr_out)

    def clustering(self):

        pr_out = pathFolder.createFolder(self.pr_out + "ClustAlgoTest/")
        runExternal.Clust(self.p_desc_cleaned, pr_out)


    def histDesc(self):
        
        pr_out = pathFolder.createFolder(self.pr_out + "histDesc/")
        runExternal.drawHist(self.p_desc, self.p_desc_cleaned, pr_out)