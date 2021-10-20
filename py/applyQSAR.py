from random import shuffle
from os import path, listdir
import joblib
from re import search
from tensorflow.keras.models import load_model
from statistics import mean

import toolbox
import pathFolder
import runExternal
import ML_toolbox


class applyQSAR:
    def __init__(self, c_dataset, pr_QSAR, pr_out):

        self.c_dataset = c_dataset
        self.pr_out = pr_out
        self.pr_models = pr_QSAR

    def loadDataFromCrossRef(self, dataset, aff = ""):
        """
        Build the dataset
        arg: -class setereo to remove single dose on the assays
        """

        p_filout = self.pr_out + "aff.csv"
        p_sum = self.pr_out + "aff.sum"


        d_set = toolbox.loadMatrix(self.c_dataset.d_dataset[dataset])
        self.d_set = d_set

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
        self.dataset = dataset


    def buildDescSet(self, l_desc):

        # toke for all chemicals
        if not "OPERA" in l_desc:
            p_opera = "0"
        else:
            p_opera = self.c_dataset.c_Desc.d_desc["all"]["OPERA"]

        if not "toxprint" in l_desc:
            p_toxprint = "0"
        else:
            # need to write the toxprint
            p_toxprint = self.pr_out + "toxprint.csv"
            l_toxprint = list(self.c_dataset.c_FP.d_toxprint[next(iter(self.c_dataset.c_FP.d_toxprint))].keys())
            l_toxprint.remove("INPUT")
            l_toxprint.remove("PREFERRED_NAME")
            l_toxprint.remove("DTXSID")
            f_toxprint = open(p_toxprint, "w")
            f_toxprint.write("CASRN\t" + "\t".join(l_toxprint) + "\n")
            for chem in self.c_dataset.c_FP.d_toxprint.keys():
                f_toxprint.write("%s\t%s\n"%(chem, "\t".join([str(self.c_dataset.c_FP.d_toxprint[chem][toxprint]) for toxprint in l_toxprint])))
            f_toxprint.close()

        runExternal.combineDesc(self.c_dataset.c_Desc.d_desc[self.dataset]["rdkit"], p_opera, p_toxprint, self.pr_out)

        if path.exists(self.pr_out + "desc_global.csv"):
            self.p_desc = self.pr_out + "desc_global.csv"
        

    def applyAllModel(self):

        l_pr_models = listdir(self.pr_models)

        # load train set to select
        p_train = self.pr_models + "trainGlobal.csv"
        d_train = ML_toolbox.loadSet(p_train)
        d_topred = ML_toolbox.loadSet(self.p_desc, d_train["features"], "", "\t")

        d_pred = {}

        for pr_model in l_pr_models:
            if pr_model == "AD":
                continue
            if path.isdir(self.pr_models + pr_model + "/"):
                l_model_files = listdir(self.pr_models + pr_model + "/")
                for model_file in l_model_files:
                    if model_file == "model.RData":
                        
                        # make a R prediction
                        p_model = self.pr_models + pr_model + "/" + model_file
                        p_predict = "%sperf_%s.csv"%(self.pr_out, pr_model)
                        runExternal.applyQSARModel(self.p_desc, p_model, pr_model, p_predict)

                        d_pred_model = toolbox.loadMatrix(p_predict, sep = ",")
                        d_pred[pr_model] = [d_pred_model[ID]["Pred"] for ID in d_pred_model.keys()]                        

                    elif model_file == "model.joblib":

                        # load model
                        loaded_model = joblib.load(self.pr_models + pr_model + "/" + model_file)

                        # apply model
                        y_pred = loaded_model.predict_proba(d_topred["dataset"])

                        try:y_pred = [pred[1] for pred in y_pred]
                        except:y_pred = [pred[0] for pred in y_pred]

                        d_pred[pr_model] = y_pred


                    elif search(".h5", model_file):
                        model = load_model(self.pr_models + pr_model + "/" + model_file)
                        y_pred = model.predict(d_topred["dataset"])    
                        try:y_pred = [pred[1] for pred in y_pred]
                        except:y_pred = [pred[0] for pred in y_pred] 
                        d_pred[pr_model] = y_pred 
                
        p_filout = self.pr_out + "predict_all.csv"
        l_ml = list(d_pred.keys())
        
        filout = open(p_filout, "w")
        filout.write("CASRN\tChemical name\t" + "\t".join(l_ml) + "\tAvg\n")
        i = 0
        imax = len(d_topred["id"])
        while i < imax:
            med = [float(d_pred[ml][i]) for ml in l_ml]
            med = mean(med)
            filout.write("%s\t%s\t%s\t%.2f\n"%(d_topred["id"][i], self.d_set[d_topred["id"][i]]["Chemical name"], "\t".join([d_pred[ml][i] for ml in l_ml]), med))
            i = i + 1
        filout.close()


        print(d_pred)                  

        return 
