from random import shuffle
from os import path, listdir
import joblib

import toolbox
import pathFolder
import runExternal


class applyQSAR:
    def __init__(self, c_dataset, pr_QSAR, pr_out):

        self.c_dataset = c_dataset
        self.pr_out = pr_out
        self.pr_models = pr_QSAR

    def buildDataset(self, dataset, aff = ""):
        """
        Build the dataset
        arg: -class setereo to remove single dose on the assays
        """

        p_filout = self.pr_out + "aff.csv"
        p_sum = self.pr_out + "aff.sum"


        d_set = toolbox.loadMatrix(self.c_dataset.d_dataset[dataset])

        d_dataset = {}
        for chem in d_set.keys():
            if aff == 1:
                d_dataset[chem] = "1"
            else:
                # take aff in the dataset
                d_dataset[chem] = d_set[chem]["aff"]
        

        f_sum = open(p_sum, "w")
        f_sum.write("Nb chemicals: %s\n"%(len(list(d_set.keys()))))
        f_sum.write("Aff formated: %s\n"%(aff))
        f_sum.close()

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
        pr_desc = pathFolder.createFolder("%s%s_%s-%s/"%(self.pr_out,"-".join(l_desc)))

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

    def applyAllModel(self):

        l_pr_models = listdir(self.pr_models)

        for pr_model in l_pr_models:
            if pr_model == "AD":
                continue
            if path.isdir(self.pr_models + pr_model + "/"):
                l_model_files = listdir(self.pr_models + pr_model + "/")
                for model_file in l_model_files:
                    if model_file == "model.RData":
                        # make a R prediction
                        p_model = self.pr_models + pr_model + "/" + model_file
                        runExternal.applyQSARModel(self.p_desc, p_model, pr_model, "%sperf_%s.csv"%(self.pr_out, pr_model))

                    elif model_file == "model.joblib":

                        # load model
                        loaded_model = joblib.load(self.pr_models + pr_model + "/" + model_file)

                        # apply model
                        result = loaded_model.score(X_test, Y_test)
                        pass
                        # make a python prediction

        return 
