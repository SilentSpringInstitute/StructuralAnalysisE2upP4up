from os import path
from random import shuffle

import pathFolder
import QSAR
import toolbox
import runExternal

class buildQSAR:

    def __init__(self, name_QSAR, active_dataset, inactive_set, c_dataset, pr_results, cor_val, quantile_desc):
        self.name_QSAR = name_QSAR
        self.pr_results = pr_results
        self.cor_val = cor_val
        self.quantile_desc = quantile_desc
        self.active_dataset = active_dataset
        self.inactive_dataset = inactive_set
        self.c_dataset = c_dataset

        # define exit folder
        self.pr_out = pathFolder.createFolder(pr_results + name_QSAR + "/")

    def buildDataset(self):

        p_filout = self.pr_out + "aff.csv"
        # define a class file
        d_dataset = {}
        d_active_set = toolbox.loadMatrix(self.c_dataset.d_dataset[self.active_dataset])
        for chem in d_active_set.keys():
            d_dataset[chem] = "1"
        
        ## inactive set
        d_inactive = toolbox.loadMatrix(self.c_dataset.d_dataset[self.inactive_dataset])
        for chem in d_inactive.keys():
            if not chem in list(d_dataset.keys()):
                d_dataset[chem] = "NA"
        
        filout = open(p_filout, "w")
        filout.write("CASRN\tAff\n")
        l_casrn = list(d_dataset.keys())
        shuffle(l_casrn)
        for casrn in l_casrn:
            filout.write("%s\t%s\n"%(casrn, d_dataset[casrn]))
        filout.close()
        self.d_dataset = d_dataset
        self.p_aff = p_filout
    
    def buildDescSet(self, l_desc):
        pr_desc = pathFolder.createFolder(self.pr_out + "-".join(l_desc) + "/")

        # toke for all chemicals
        runExternal.combineAndPredDesc(self.c_dataset.c_Desc.d_desc["all"]["rdkit"], self.c_dataset.c_Desc.d_desc["all"]["OPERA"], pr_desc)

        if path.exists(pr_desc + "desc_global.csv"):
            self.p_desc = pr_desc + "desc_global.csv"

    def prepDesc(self):
        pr_out = pathFolder.createFolder(self.pr_out + "Cleaned_Data/")

        # add shortcut
        p_desc_cleaned = pr_out + "desc_cleaned.csv"
        p_AC50_cleaned = pr_out + "AC50_cleaned.csv"
        if not path.exists(p_desc_cleaned) and not path.exists(p_AC50_cleaned):
            runExternal.preprocDataQSAR(self.p_desc, self.p_aff, pr_out, self.cor_val, self.quantile_desc)

        self.p_desc_cleaned = p_desc_cleaned
        self.p_AC50_cleaned = p_AC50_cleaned

    def runQSARs(self):

        pr_out = pathFolder.createFolder(self.pr_out + "classQSAR/")
        self.c_QSAR = QSAR.QSAR(self.p_desc_cleaned, self.p_desc, self.p_AC50_cleaned, self.p_aff, pr_out, 5, 10, 0, 0.15)
        self.c_QSAR.runQSARClassUnderSamplingAllSet()
