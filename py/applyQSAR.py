from random import shuffle
from os import path, listdir
import joblib
from re import T, search
from tensorflow.keras.models import load_model
from statistics import mean

import toolbox
import pathFolder
import runExternal
import ML_toolbox
import CompDesc

class applyQSAR:
    def __init__(self, c_dataset, pr_QSAR, pr_out):

        self.c_dataset = c_dataset
        self.pr_out = pr_out
        self.pr_models = pr_QSAR

    def loadDataFromCrossRef(self, dataset_test, l_dataset_train, aff = ""):
        """
        Build the dataset
        arg: -class setereo to remove single dose on the assays
        """

        p_filout = self.pr_out + "aff.csv"
        p_sum = self.pr_out + "aff.sum"

        self.dataset_train_act = l_dataset_train[0]
        self.dataset_train_inact = l_dataset_train[1]

        self.dataset_test = dataset_test

        d_set_test = toolbox.loadMatrix(self.c_dataset.d_dataset[dataset_test])
        self.d_set_test = d_set_test

        d_set = toolbox.loadMatrix(self.c_dataset.d_dataset[l_dataset_train[0]])
        self.d_set_train_act = d_set        

        d_set = toolbox.loadMatrix(self.c_dataset.d_dataset[l_dataset_train[1]])
        self.d_set_train_inact = d_set

        d_dataset = {}
        for chem in d_set_test.keys():
            if aff == 1:
                d_dataset[chem] = "1"
            else:
                # take aff in the dataset
                d_dataset[chem] = d_set_test[chem]["aff"]
        

        f_sum = open(p_sum, "w")
        f_sum.write("Nb chemicals: %s\n"%(len(list(d_set_test.keys()))))
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

        runExternal.combineDesc(self.c_dataset.c_Desc.d_desc[self.dataset_test]["rdkit"], p_opera, p_toxprint, self.pr_out)

        if path.exists(self.pr_out + "desc_global.csv"):
            self.p_desc = self.pr_out + "desc_global.csv"        

    def applyQSARModels(self):

        l_pr_models = listdir(self.pr_models)

        # load train set to select
        p_train = self.pr_models + "trainGlobal.csv"
        p_test = self.pr_models + "test.csv"
        self.d_train = ML_toolbox.loadSet(p_train)
        self.d_test = ML_toolbox.loadSet(p_test)
        d_topred = ML_toolbox.loadSet(self.p_desc, self.d_train["features"], "", "\t")

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

                        # need to scale if ghost opt
                        if search("ghost", pr_model):
                            p_prob_cutoff = self.pr_models + pr_model + "/ghost_threshold.txt"
                            filin = open(p_prob_cutoff, "r")
                            cutoff = float(filin.read())
                            filin.close()
                            y_pred = [cutoff * float(y) / 0.5 for y in y_pred]
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
        filout.write("CASRN\tChemical name\t" + "\t".join(l_ml) + "\tAvg\tAvg RF\tIn train\tIn test\t" + self.dataset_train_act + "\t" + self.dataset_train_inact + "\n")
        i = 0
        imax = len(d_topred["id"])
        l_med_rf = []
        d_pred["Avg"] = {}
        d_pred["Avg RF"] = {}
        d_pred["In train"] = {}
        d_pred["In test"] = {}
        while i < imax:
            med = [float(d_pred[ml][i]) for ml in l_ml]
            med_rf = [float(d_pred[ml][i]) for ml in ["RF_py_ghost", "RF_balanced_ghost"]]
            med = mean(med)
            med_rf = mean(med_rf)
            l_med_rf.append(med_rf)
            d_pred["Avg"][i] = med
            d_pred["Avg RF"][i] = med_rf
            
            # check if included in the train set
            CASRN  = d_topred["id"][i]
            if CASRN in self.d_train["id"]:
                d_pred["In train"][i] = "1"
            else:
                d_pred["In train"][i] = "0"
            
            if CASRN in self.d_test["id"]:
                d_pred["In test"][i] = "1"
            else:
                d_pred["In test"][i] = "0"
            i = i + 1


        l_med_rf.sort(reverse=T)
        l_i = [i for i in range(0, len(d_topred["id"]))]
        for med_rf in l_med_rf:
            for i in l_i:
                if d_pred["Avg RF"][i] == med_rf:
                    try: class_train_act = self.d_set_train_act[d_topred["id"][i]][self.dataset_train_act]
                    except: class_train_act = "-"

                    try: class_train_inact = self.d_set_train_inact[d_topred["id"][i]][self.dataset_train_inact]
                    except: class_train_inact = "-"

                    filout.write("%s\t%s\t%s\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n"%(d_topred["id"][i], self.d_set_test[d_topred["id"][i]]["Chemical name"], "\t".join([str(d_pred[ml][i]) for ml in l_ml]), d_pred["Avg"][i], d_pred["Avg RF"][i], d_pred["In train"][i], d_pred["In test"][i], class_train_act, class_train_inact))
                    l_i.remove(i) 
        filout.close()

        self.p_pred_QSAR = p_filout


    def computeAD(self):

        # need to compute before the prediction
        if not "p_pred_QSAR" in self.__dict__:
            print("Compute first the QSAR prediction")
            return 
        
        d_pred_all = toolbox.loadMatrix(self.p_pred_QSAR)
        pr_AD = pathFolder.createFolder(self.pr_out + "AD/")

        #PCA for applicability model #
        ##############################
        pr_PCA = pathFolder.createFolder(pr_AD + "PCA_desc/")
        # load train set to select
        p_train = self.pr_models + "trainGlobal.csv"
        runExternal.PCA_2sets(p_train, self.p_desc, pr_PCA)

        # AD based on similarity score #
        ################################
        pr_AD_sim = pathFolder.createFolder(pr_AD + "chem_similarity/")
        if not path.exists(pr_AD_sim + "CP_similarity_text.png"):
            # compute similarity matrix in the root folder.
            # create matrix with flag active vs inactive and test vs train
            p_matrix_chem = pr_AD_sim + "chem.csv"
            f_matrix_chem = open(p_matrix_chem, "w")
            f_matrix_chem.write("CASRN\tAff\tset\n")
            
            # set to test
            d_SMILES = {}
            for chem in d_pred_all.keys():
                if chem in list(self.d_set_train_act.keys()) or chem in list(self.d_set_train_inact.keys()):
                    continue
                d_SMILES[chem] = self.d_set_test[chem]["SMILES"]
                f_matrix_chem.write("%s\t%s\t%s\n"%(chem, "1", "test"))
            
            # train active
            for casrn_train in self.d_set_train_act.keys():
                d_SMILES[casrn_train] = self.d_set_train_act[casrn_train]["SMILES"]
                f_matrix_chem.write("%s\t%s\t%s\n"%(casrn_train, "1", "train"))
            
            # train inactive
            for casrn_train in self.d_set_train_inact.keys():
                d_SMILES[casrn_train] = self.d_set_train_inact[casrn_train]["SMILES"]
                f_matrix_chem.write("%s\t%s\t%s\n"%(casrn_train, "0", "train"))
            f_matrix_chem.close()

            # define p_sim_matrix
            p_matrix_sim = pr_AD_sim + "matrix_sim.csv"
            
            if not path.exists(p_matrix_sim):
                d_sim = {}
                l_casrn = list(d_SMILES.keys())
                i = 0
                imax = len(l_casrn)
                while i < imax:
                    smi1 = d_SMILES[l_casrn[i]]
                    d_sim[l_casrn[i]] = {}
                    c_smi1 = CompDesc.CompDesc(smi1, "")
                    c_smi1.prepChem()
                    c_smi1.computeFP("MACCS")

                    j = i + 1
                    while j < imax:
                        smi2 = d_SMILES[l_casrn[j]]
                        c_smi2 = CompDesc.CompDesc(smi2, "")
                        c_smi2.prepChem()
                        c_smi2.computeFP("MACCS")
                        score = c_smi1.computeSimilarityFP(c_smi2, "MACCS", "Tanimoto")
                        d_sim[l_casrn[i]][l_casrn[j]] = str(score)
                        j=j+1

                    i = i + 1
            
                # matrix of similarity
                filout = open(p_matrix_sim, "w")
                filout.write("\t".join(l_casrn) + "\n")
                for casrn in l_casrn:
                    filout.write(casrn)
                    for casrn2 in l_casrn:
                        if casrn == casrn2:
                            filout.write("\t1")
                        else:
                            try:filout.write("\t%s"%(d_sim[casrn][casrn2]))
                            except: filout.write("\t%s"%(d_sim[casrn2][casrn]))
                    filout.write("\n")
                filout.close()

            runExternal.computeADBasedOnSimilarityMatrix(p_matrix_sim, p_matrix_chem, pr_AD_sim)  


    def applyToxPrintSignifcant(self, pr_results):

        # define folder with the list of significant toxprint
        p_toxprint_signif = "%scomparisonDescToxprint_%s-%s/Toxprint/signif.csv"%(pr_results, self.dataset_train_act, self.dataset_train_inact)
        d_toxprint_signif = toolbox.loadMatrix(p_toxprint_signif, sep = ",")

        d_count = {}
        l_sum = []
        for CASRN in self.d_dataset.keys():
            d_count[CASRN] = {}
            d_count[CASRN]["***"] = 0
            d_count[CASRN]["**"] = 0
            d_count[CASRN]["*"] = 0
            d_count[CASRN]["-"] = 0

            if CASRN in self.d_train["id"]:
                d_count[CASRN]["In train"] = "1"
            else:
                d_count[CASRN]["In train"] = "0"
            
            if CASRN in self.d_test["id"]:
                 d_count[CASRN]["In test"] = "1"
            else:
                d_count[CASRN]["In test"] = "0"

            
            try: self.c_dataset.c_FP.d_toxprint[CASRN]
            except: 
                continue

            for toxprint in self.c_dataset.c_FP.d_toxprint[CASRN].keys():
                if self.c_dataset.c_FP.d_toxprint[CASRN][toxprint] == "0":
                    continue
                toxprint_match = toxprint.replace(":", '.').replace("-", ".")

                for id_toxprint in d_toxprint_signif.keys():
                    toxprint_signif = d_toxprint_signif[id_toxprint]["Toxprint"]
                    signif = d_toxprint_signif[id_toxprint]["significatif"]

                    if toxprint_match == toxprint_signif:
                        if signif == "***":
                            d_count[CASRN]["***"] = d_count[CASRN]["***"] + 1
                        elif signif == "**":
                            d_count[CASRN]["**"] = d_count[CASRN]["**"] + 1
                        elif signif == "*":
                            d_count[CASRN]["*"] = d_count[CASRN]["*"] + 1
                        elif signif == "-":
                            d_count[CASRN]["-"] = d_count[CASRN]["-"] + 1
                
            l_sum.append(d_count[CASRN]["***"])



        l_sum = list(set(l_sum))
        l_sum.sort(reverse=T)


        p_filout = self.pr_out + "ToxPrints.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\tChemical name\t***\t**\t*\t-\tIn train\tIn test\t" + self.dataset_train_act + "\t" + self.dataset_train_inact + "\n")
        i = 0
        imax = len(l_sum)
        while i < imax:
            for casrn in d_count.keys():
                if d_count[casrn]["***"] == l_sum[i]:
                    try: class_train_act = self.d_set_train_act[casrn][self.dataset_train_act]
                    except: class_train_act = "-"

                    try: class_train_inact = self.d_set_train_inact[casrn][self.dataset_train_inact]
                    except: class_train_inact = "-"

                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(casrn, self.d_set_test[casrn]["Chemical name"], d_count[casrn]["***"], d_count[casrn]["**"], d_count[casrn]["*"], d_count[casrn]["-"], d_count[casrn]["In train"], d_count[casrn]["In test"], class_train_act, class_train_inact))
            i = i + 1
        filout.close()

        self.p_model_toxprint = p_filout

    
    def mergePredToxPrintQSAR(self):

        d_QSAR_pred = toolbox.loadMatrix(self.p_pred_QSAR)
        d_toxprint_pred = toolbox.loadMatrix(self.p_model_toxprint)

        p_filout = self.pr_out + "ToxPrint_QSAR.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\tChemical name\tnb Toxprint (***)\tnb Toxprint (**)\tnb Toxprint (*)\tPred RF\tPred RF balanced\tAvg. RF\tAD distance to first neighbord\n")

        for casrn in d_toxprint_pred.keys():
            if d_toxprint_pred[casrn]["In train"] == "1" or d_toxprint_pred[casrn]["In test"] == "1":
                continue
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(casrn, d_toxprint_pred[casrn]["Chemical name"], d_toxprint_pred[casrn]["***"], d_toxprint_pred[casrn]["**"], d_toxprint_pred[casrn]["*"], d_QSAR_pred[casrn]["RF_py_ghost"], d_QSAR_pred[casrn]["RF_balanced_ghost"], d_QSAR_pred[casrn]["Avg RF"]))
        filout.close()



        return