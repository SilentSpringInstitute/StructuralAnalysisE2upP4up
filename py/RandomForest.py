import pathFolder
import toolbox
import ghost


from numpy import loadtxt, arange
from re import search
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.preprocessing import StandardScaler
sc = StandardScaler()
import joblib
from os import path, listdir
from copy import deepcopy
import numpy as np


from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV

from imblearn.ensemble import BalancedRandomForestClassifier

class RandomForest:

    def __init__(self, p_train, p_test, p_aff, n_foldCV, type_ML, ghost, pr_out):
        self.p_train = p_train
        self.p_test = p_test
        self.p_aff = p_aff
        self.ghost = ghost
        self.l_ghost_threshold = np.round(np.arange(0.05,0.55,0.05),2)
        self.n_foldCV = n_foldCV
        self.verbose = 0
        self.test = 0
        self.force_run = 0
        self.typeModel = "classification"
        self.type_ML = type_ML
        self.pr_out = pr_out

        self.d_model = {}
        # grid optimization - use a test criteria to reduce the grid size for testing
        self.test = 0
        self.force_run = 0

        self.n_jobs = 4
        self.random_state = 42
        if self.test == 0:
            self.n_estimators =  [10, 50, 100, 300, 500, 800, 1200]
            self.max_depth = [5, 8, 15, 25, 30]
            self.min_samples_leaf = [1, 2, 5, 10] 
            self.min_samples_split = [2, 5, 10, 15, 100]
            self.n_iter = 100
        else:
            self.n_estimators =  [10]
            self.max_depth = [5]
            self.min_samples_leaf = [10] 
            self.min_samples_split = [2]
            self.n_iter = 20

    def loadSet(self):
    
        # descriptor train set 
        with open(self.p_train) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        self.dataset_train = loadtxt(self.p_train, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        self.dataset_train = sc.fit_transform(self.dataset_train)
        self.aff_train = loadtxt(self.p_train, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)
        self.id_train = loadtxt(self.p_train, dtype=str,  delimiter=",",usecols=arange(0, 1), skiprows=1)
        self.nb_desc_input = n_cols - 2


        # descriptor test set 
        with open(self.p_test) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        self.dataset_test = loadtxt(self.p_test, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        self.dataset_test = sc.fit_transform(self.dataset_test)
        self.aff_test = loadtxt(self.p_test, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)
        self.id_test = loadtxt(self.p_test, dtype=str, delimiter=",",usecols=arange(0, 1), skiprows=1)

    def run_RF(self, CV = 0, **kwargs):

        # add ghost here
        if self.ghost == 1:
            type_RF = self.type_ML + "_ghost"
        else:
            type_RF = self.type_ML

        self.pr_results = pathFolder.createFolder(self.pr_out + type_RF + "/")

        if CV == 1:
            pr_out = self.pr_results + "CV_"
        else:
            pr_out = self.pr_results

        # short cut with the model already computed and save
        p_model = pr_out + "model.joblib"
        
        # do not check if CV model exist
        if path.exists(p_model) and self.force_run == 0 and CV == 0:
            self.d_model[type_RF] = joblib.load(p_model)

        else:
            random_grid = {'n_estimators': self.n_estimators,
                'max_depth': self.max_depth,
                'min_samples_split': self.min_samples_split,
                'min_samples_leaf': self.min_samples_leaf}

            if self.type_ML == "RF_balanced":
                rf = BalancedRandomForestClassifier(n_jobs=self.n_jobs)
            else:
                rf = RandomForestClassifier(n_jobs=self.n_jobs)

            rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = self.n_iter, cv = self.n_foldCV, verbose=2, random_state=self.random_state, n_jobs = self.n_jobs)
            rf_random.fit(self.dataset_train, self.aff_train)
            best_random = rf_random.best_estimator_
            
            # save the model
            joblib.dump(best_random, p_model)

            # do not save CV on the save as RF
            if CV == 0:
                self.d_model[type_RF] = joblib.load(p_model)

            else:
                self.d_model["CV_" + type_RF] = joblib.load(p_model)

    def apply_model(self, a_set, id_set, type_RF, name_out = "", w=0): 
       
        y_pred = self.d_model[type_RF].predict_proba(a_set)

        # compute quality performance
        if w == 1:
            p_filout = self.pr_results + name_out
            filout = open(p_filout, "w")
            filout.write("\"\",\"ID\",\"Pred\"\n")
            i = 0
            l_chem = list(id_set)
            imax = len(l_chem)
            
            while i < imax:
                filout.write("\"%s\",\"%s\",\"%s\"\n"%(l_chem[i], l_chem[i], y_pred[i][1]))
                i = i + 1
            filout.close()
            
        return y_pred
    
    def performance(self, y_real, y_pred, th_prob = 0.5, p_filout = ""):
            
        if self.typeModel == "classification":

            # change prob -> value
            y_pred = [1. if pred[1] > th_prob else 0. for pred in y_pred]

            acc = metrics.accuracy_score(y_real, y_pred)
            bacc = metrics.balanced_accuracy_score(y_real, y_pred)
            mcc = metrics.matthews_corrcoef(y_real, y_pred)
            recall = metrics.recall_score(y_real, y_pred)
            roc_auc = metrics.roc_auc_score(y_real, y_pred)
            f1b = metrics.fbeta_score(y_real, y_pred, beta=0.5)

            # specificity & sensitivity 
            tn, fp, fn, tp = metrics.confusion_matrix(y_real, y_pred).ravel()
            specificity = float(tn) / (tn+fp)
            sensitivity = float(tp) / (tp+fn)

            print("======= PERF ======")
            print("Acc: ", acc)
            print("b-Acc: ", bacc)
            print("Sp: ", specificity)
            print("Se: ", sensitivity)
            print("MCC: ", mcc)
            print("Recall: ", recall)
            print("AUC: ", roc_auc)
            print("fb1: ", f1b)

            return {"Acc": acc, "b-Acc": bacc, "MCC": mcc, "Recall": recall, "AUC": roc_auc, "Se": sensitivity, "Sp": specificity, "f1b": f1b}
        
        else:
            MAE = metrics.mean_absolute_error(y_real, y_pred)
            R2 = metrics.r2_score(y_real, y_pred)
            EVS = metrics.explained_variance_score(y_real, y_pred)
            MSE = metrics.mean_squared_error(y_real, y_pred)
            MAXERR = metrics.max_error(y_real, y_pred)
            try:MSE_log = metrics.mean_squared_log_error(y_real, y_pred)
            except: MSE_log = 0.0
            MDAE = metrics.median_absolute_error(y_real, y_pred)
            MTD = metrics.mean_tweedie_deviance(y_real, y_pred)
            try:
                MPD = metrics.mean_poisson_deviance(y_real, y_pred)
                MGD = metrics.mean_gamma_deviance(y_real, y_pred)
            except:
                MPD = 0.0
                MGD = 0.0

            print("======= TRAIN ======")
            print("MAE: ", MAE)
            print("R2: ", R2)
            print("Explain Variance score: ", EVS)
            print("MSE: ", MSE)
            print("Max error: ", MAXERR)
            print("MSE log: ", MSE_log)
            print("Median absolute error: ", MDAE)
            print("Mean tweedie deviance: ", MTD)
            print("Mean poisson deviance: ", MPD)
            print("Mean gamma deviance: ", MGD)

            if p_filout != "":
                filout = open(p_filout, "w")
                filout.write("\tMAE\tR2\tExplain Variance score\tMSE\tMax error\tMSE log\tMedian absolute error\tMean tweedie deviance\tMean poisson deviance\tMean gamma deviance\n")
                filout.write("TEST-DNN\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(MAE, R2, EVS, MSE, MAXERR, MSE_log, MDAE, MTD, MPD, MGD))
                filout.close()


            return {"MAE": MAE, "R2": R2, "EVS": EVS, "MSE": MSE, "MAXERR": MAXERR, "MSE_log": MSE_log , "MDAE": MDAE, "MTD": MTD, "MPD":MPD , "MGD":MGD}

    def evaluateOnTest(self, th_prob=0.5):
    
        y_pred = self.model.predict_proba(self.dataset_test)
        if self.typeModel == "classification":
            return self.performance(self.aff_test, y_pred, self.threshold_ghost)
        else:
            return self.performance(self.aff_test, y_pred, self.threshold_ghost)

    def CrossValidation(self):
        """
        Realize a cross validation in N folds
        args: - type of RF (balanced/umbalanced) 
        """
        
        if self.ghost == 1:
            type_RF = self.type_ML + "_ghost"
        else:
            type_RF = self.type_ML
        self.pr_results = pathFolder.createFolder(self.pr_out + type_RF + "/")

        seed = 7
        d_CV = {}
        if self.typeModel == "classification":
            kfold = StratifiedKFold(n_splits=self.n_foldCV, shuffle=True, random_state=seed)
        else:
            kfold = KFold(n_splits=self.n_foldCV, shuffle=True, random_state=seed)
        # load parameters
        
        d_train = deepcopy(self.dataset_train)
        Y_train = deepcopy(self.aff_train)
        id_train = deepcopy(self.id_train)

        
        for train, test in kfold.split(d_train, Y_train):
            
            # update the class
            self.dataset_train = d_train[train]
            self.aff_train = Y_train[train]
            self.id_train = id_train[train]

            self.dataset_test = d_train[test]
            self.aff_test = Y_train[test]
            self.id_test = id_train[test]

            # optimize model
            self.run_RF(CV=1)

            # apply prediction on the test
            y_pred = self.apply_model(self.dataset_test, self.id_test, "CV_" + type_RF)
            if self.ghost == 1:
                y_pred_train = self.apply_model(self.dataset_train, self.id_train, "CV_" + type_RF)
                y_pred_train_ghost = [pred[1] for pred in y_pred_train]
                self.threshold_ghost = ghost.optimize_threshold_from_predictions(self.aff_train, y_pred_train_ghost, self.l_ghost_threshold, ThOpt_metrics = 'Kappa') 
            else:
                self.threshold_ghost = 0.5

            d_pref_test = self.performance(self.aff_test, y_pred, self.threshold_ghost)
            
            for perf_criteria in d_pref_test.keys():
                if not perf_criteria in list(d_CV.keys()):
                    d_CV[perf_criteria] = []
                d_CV[perf_criteria].append(d_pref_test[perf_criteria])

        l_criteria = list(d_CV.keys())
        l_val = ["%.2f"%(np.mean(d_CV[criteria])) for criteria in l_criteria]
        filout = open(self.pr_results + "CV_%s.perf"%(self.n_foldCV), "w")
        filout.write("%s\n%s"%("\t".join(l_criteria), "\t".join(l_val)))

        # populate d_pref
        if not "d_perf" in self.__dict__:
            self.d_perf = {}
        
        self.d_perf["CV"] = {}
        i = 0
        imax = len(l_criteria)
        while i < imax:
            self.d_perf["CV"][l_criteria[i]] = float(l_val[i])
            i  = i + 1   

    def TrainTestPrediction(self):
        
        # load dataset
        self.loadSet()

        # optimize model
        self.run_RF(CV = 0)
        
        if not "d_perf" in self.__dict__:
            self.d_perf = {}

        if self.typeModel == "classification":
            # trainning set
            # add ghost here
            if self.ghost == 1:
                type_RF = self.type_ML + "_ghost"
            else:
                type_RF = self.type_ML
            y_train_pred = self.apply_model(self.dataset_train, self.id_train, type_RF, "train_pred.csv", w=1)
            if self.ghost == 1:
                y_train_pred_ghost = [pred[1] for pred in y_train_pred]                
                threshold_ghost = ghost.optimize_threshold_from_predictions(self.aff_train, y_train_pred_ghost, self.l_ghost_threshold, ThOpt_metrics = 'Kappa') 
            else:
                threshold_ghost = 0.5

            self.ghost_treshold = threshold_ghost
            self.d_perf["Train"] = self.performance(self.aff_train, y_train_pred, self.ghost_treshold)
            
            # test set
            y_pred = self.apply_model(self.dataset_test, self.id_test, type_RF, "test_pred.csv", w=1)
            self.d_perf["Test"] = self.performance(self.aff_test, y_pred, self.ghost_treshold)
           
    def combineResults(self):
       
        if not "d_perf" in self.__dict__:
            print("No data to write")
            return 

        if self.ghost == 1:
            type_RF = self.type_ML + "_ghost"
        else:
            type_RF = self.type_ML
        
        p_filout = self.pr_results + "combined_perf.csv"

        filout = open(p_filout, "w")
        l_perf = ["CV", "Train", "Test"]
        l_criteria = list(self.d_perf["CV"].keys())
        
        filout.write("\t%s\n"%("\t".join(l_criteria)))
        
        for perf in l_perf:
            filout.write("%s-%s\t%s\n"%(perf, type_RF, "\t".join(["%.2f"%(self.d_perf[perf][criteria]) for criteria in l_criteria])))
        filout.close() 
    
    def run_main(self):

        if self.ghost == 1:
            type_RF = self.type_ML + "_ghost"
        else:
            type_RF = self.type_ML

        p_combined_results = self.pr_out + type_RF + "/" + "combined_perf.csv"
        if path.exists(p_combined_results) and self.force_run == 0:
            return  

        self.TrainTestPrediction()
        self.CrossValidation()
        self.combineResults()





