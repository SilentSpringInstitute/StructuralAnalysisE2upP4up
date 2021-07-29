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

    def buildDataset(self, c_stereo = 0, borderline=1):
        """
        Build the dataset
        - arg: class setereo to remove single dose on the assays
        """
        p_filout = self.pr_out + "aff.csv"
        
        # redefine E2up and P4up with borderline and not
        if borderline == 0:
            self.c_dataset.defineE2P4active(["higher", "medium", "lower"])
            if self.active_dataset == "E2up":
                l_casr_noBorderline = list(self.c_dataset.d_E2up_active.keys())
            else:
                l_casr_noBorderline = list(self.c_dataset.d_P4up_active.keys())
            
            self.c_dataset.defineE2P4active(["higher", "medium", "lower", "borderline"])
            if self.active_dataset == "E2up":
                l_casr_Borderline = list(self.c_dataset.d_E2up_active.keys())
            else:
                l_casr_Borderline = list(self.c_dataset.d_P4up_active.keys())

            # diff list
            l_Borderline = list(set(l_casr_Borderline) - set(l_casr_noBorderline))



        # define a class file
        d_dataset = {}
        d_active_set = toolbox.loadMatrix(self.c_dataset.d_dataset[self.active_dataset])
        for chem in d_active_set.keys():
            if borderline == 1:
                if chem in l_Borderline:
                    continue
            d_dataset[chem] = "1"
        
        ## inactive set
        d_inactive = toolbox.loadMatrix(self.c_dataset.d_dataset[self.inactive_dataset])
        for chem in d_inactive.keys():
            if not chem in list(d_dataset.keys()):
                
                # check for borderline to remove them - there are included in the dataset folder
                if borderline == 0:
                    if chem in l_Borderline:
                        continue

                # check for single dose remove
                if c_stereo == 0:
                    d_dataset[chem] = "NA"
                else:
                    #remove single dose in the inactive set
                    if self.active_dataset == "E2up":
                        if not chem in list(c_stereo.d_single_hit.keys()) or c_stereo.d_single_hit[chem]["ESTRADIOL"] == "0" :
                            d_dataset[chem] = "NA"
                    if self.active_dataset == "P4up":
                        if not chem in list(c_stereo.d_single_hit.keys()) or c_stereo.d_single_hit[chem]["PROG"] == "0":
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
        pr_desc = pathFolder.createFolder("%s%s_%s-%s/"%(self.pr_out,"-".join(l_desc), self.cor_val, self.quantile_desc))

        # toke for all chemicals
        if not "OPERA" in l_desc:
            p_opera = "0"
        else:
            p_opera = self.c_dataset.c_Desc.d_desc["all"]["OPERA"]

        if not "toxprint" in l_desc:
            p_toxprint = "0"
        else:
            # need to write the toxprint
            p_toxprint = pr_desc + "toxprint.csv"
            l_toxprint = list(self.c_dataset.c_FP.d_toxprint[next(iter(self.c_dataset.c_FP.d_toxprint))].keys())
            l_toxprint.remove("INPUT")
            l_toxprint.remove("PREFERRED_NAME")
            l_toxprint.remove("DTXSID")
            f_toxprint = open(p_toxprint, "w")
            f_toxprint.write("CASRN\t" + "\t".join(l_toxprint) + "\n")
            for chem in self.c_dataset.c_FP.d_toxprint.keys():
                f_toxprint.write("%s\t%s\n"%(chem, "\t".join([str(self.c_dataset.c_FP.d_toxprint[chem][toxprint]) for toxprint in l_toxprint])))
            f_toxprint.close()

        runExternal.combineDesc(self.c_dataset.c_Desc.d_desc["all"]["rdkit"], p_opera, p_toxprint, pr_desc)

        if path.exists(pr_desc + "desc_global.csv"):
            self.p_desc = pr_desc + "desc_global.csv"
        
        # folder once the dataset of descriptors is created
        self.pr_desc_QSAR = pr_desc

    def prepDesc(self):
        pr_out = pathFolder.createFolder(self.pr_desc_QSAR + "Cleaned_Data/")

        # add shortcut
        p_desc_cleaned = pr_out + "desc_cleaned.csv"
        p_AC50_cleaned = pr_out + "AC50_cleaned.csv"
        if not path.exists(p_desc_cleaned) and not path.exists(p_AC50_cleaned):
            runExternal.preprocDataQSAR(self.p_desc, self.p_aff, pr_out, self.cor_val, self.quantile_desc)

        self.p_desc_cleaned = p_desc_cleaned
        self.p_AC50_cleaned = p_AC50_cleaned

    def runQSARs(self, rate_undersamplin=0):

        pr_out = pathFolder.createFolder(self.pr_desc_QSAR + "classQSAR/")
        self.c_QSAR = QSAR.QSAR(self.p_desc_cleaned, self.p_desc, self.p_AC50_cleaned, self.p_aff, pr_out, 10, 5, rate_undersamplin, 0.10)
        self.c_QSAR.runQSARClassUnderSamplingAllSet()
