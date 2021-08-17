import pathFolder
import toolbox

from numpy import loadtxt, arange
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

    def __init__(self, p_train, p_test, p_aff, n_foldCV, ghost, pr_out):
        self.p_train = p_train
        self.p_test = p_test
        self.p_aff = p_aff
        self.ghost = ghost
        self.n_foldCV = n_foldCV
        self.verbose = 0
        self.test = 0
        self.force_run = 0
        self.typeModel = "classification"

        # create folder -> root of QSAR results
        if self.ghost == 1:
            self.pr_out = pathFolder.createFolder(pr_out + "RFGhostpy/")
        else:
            self.pr_out = pathFolder.createFolder(pr_out + "RFpy/")   

        # optimization
        self.n_estimators = 500
        self.max_depth = 15
        self.min_samples_leaf = 2
        self.n_jobs = 4
    
        #n_estimators=500,max_depth=15,min_samples_leaf=2,n_jobs=4,oob_score=True, **kwargs)


    def loadSet(self):
    
        # descriptor train set 
        with open(self.p_train) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        self.dataset_train = loadtxt(self.p_train, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        self.dataset_train = sc.fit_transform(self.dataset_train)
        self.aff_train = loadtxt(self.p_train, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)
        self.nb_desc_input = n_cols - 2

        d_train = toolbox.loadMatrix(self.p_train, sep = ",")
        self.l_train = list(d_train.keys())

        # descriptor test set 
        with open(self.p_test) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        self.dataset_test = loadtxt(self.p_test, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        self.dataset_test = sc.fit_transform(self.dataset_test)
        self.aff_test = loadtxt(self.p_test, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)

        d_test = toolbox.loadMatrix(self.p_test, sep = ",")
        self.l_test = list(d_test.keys())
    

    def run_balancedRF(fps_train, fps_test, labels_train, labels_test, **kwargs):
        # build the classification model:
        cls = BalancedRandomForestClassifier(n_estimators=500,max_depth=15,min_samples_leaf=2,n_jobs=4, **kwargs)
        cls.fit(fps_train, labels_train)
        probs = cls.predict_proba(fps_test)[:,1]
        return cls, probs 
    

    def run_RF(self, **kwargs):
        # RF
        n_estimators = [10, 50, 100, 300, 500, 800, 1200]
        max_depth = [5, 8, 15, 25, 30]
        min_samples_split = [2, 5, 10, 15, 100]
        min_samples_leaf = [1, 2, 5, 10] 

        random_grid = {'n_estimators': n_estimators,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf}

        rf = RandomForestClassifier(n_jobs=4)
        rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 100, cv = 10, verbose=2, random_state=42, n_jobs = 4)
        #rf_random = RandomForestClassifier(random_state = 1, max_depth = 15,  n_estimators = 500, min_samples_split = 2, min_samples_leaf = 1)
        # Fit the random search model
        print(self.aff_train)
        rf_random.fit(self.dataset_train, self.aff_train)
        best_random = rf_random.best_estimator_
        
        self.model = best_random

        y_pred = rf_random.predict_proba(self.dataset_train)
        self.performance(self.aff_train, y_pred)

        # save model
        joblib.dump(best_random, self.pr_out + "RF.joblib")
        self.p_model = self.pr_out + "RF.joblib"
        

    def evaluateModel(self, th_prob = 0.5):
    
        if not "d_perf" in self.__dict__:
            self.d_perf = {}

        if self.typeModel == "classification":
            # trainning set
            y_train_pred = self.apply_model(self.p_train)
            self.d_perf["Train"] = self.performance(self.aff_train, y_train_pred, th_prob)
            
            # test set
            y_pred = self.apply_model(self.p_test)
            self.d_perf["Test"] = self.performance(self.aff_test, y_pred, th_prob)


    def apply_model(self, p_test, p_model="", name_out = "test_pred"): 
        
        # descriptor test set 
        with open(p_test) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        dataset_test = loadtxt(p_test, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        dataset_test = sc.fit_transform(dataset_test)
        
        # remove last col
        d_test = toolbox.loadMatrix(p_test, sep = ",")

        # descriptor test set 
        if p_model != "":
            self.model = joblib.load(p_model)
        else:
            if not "model" in self.__dict__:
                print("ERROR: No model loaded")
                return 

        y_pred = self.model.predict_proba(dataset_test)


        # compute quality performance
        p_filout = self.pr_out + name_out
        filout = open(p_filout, "w")
        filout.write("\"\",\"ID\",\"Pred\"\n")
        i = 0
        l_chem = list(d_test.keys())
        imax = len(l_chem)
        
        while i < imax:
            filout.write("\"\",\"ID\",\"Pred\"\n")
            filout.write("\"%s\",\"%s\",\"%s\"\n"%(l_chem[i], l_chem[i], y_pred[i]))
            i = i + 1
        filout.close()
        
        return y_pred
    
    def performance(self, y_real, y_pred, th_prob = 0.5, p_filout = ""):
            
        if self.typeModel == "classification":
            # change prob -> value
            y_pred = [1. if pred[1] > th_prob else 0. for pred in y_pred]

            print(y_pred)
            print(y_real)
            
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
            return self.performance(self.aff_test, y_pred, th_prob)
        else:
            return self.performance(self.aff_test, y_pred, th_prob)

    def CrossValidation(self, type_RF = "RF"):
        seed = 7
        d_CV = {}
        if self.typeModel == "classification":
            kfold = StratifiedKFold(n_splits=self.n_foldCV, shuffle=True, random_state=seed)
        else:
            kfold = KFold(n_splits=self.n_foldCV, shuffle=True, random_state=seed)
        # load parameters
        l_files = listdir(self.pr_out)
        
        d_train = deepcopy(self.dataset_train)
        Y_train = deepcopy(self.aff_train)

        for train, test in kfold.split(d_train, Y_train):
            self.dataset_train = d_train[train]
            self.aff_train = Y_train[train]

            self.dataset_test = d_train[test]
            self.aff_test = Y_train[test]

            # create model
            if type_RF == "RF":
                self.run_RF()
            d_pref_test = self.evaluateOnTest()
            
            for perf_criteria in d_pref_test.keys():
                if not perf_criteria in list(d_CV.keys()):
                    d_CV[perf_criteria] = []
                d_CV[perf_criteria].append(d_pref_test[perf_criteria])

        l_criteria = list(d_CV.keys())
        l_val = ["%.2f"%(np.mean(d_CV[criteria])) for criteria in l_criteria]
        filout = open(self.pr_out + "CV_%s.perf"%(self.n_foldCV), "w")
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


    def combineResults(self, name_combine):
        if not "d_perf" in self.__dict__:
            print("No data to write")
            return 

        p_filout = self.pr_out + "combined_perf.csv"
        if path.exists(p_filout) and self.force_run == 0:
            return
            
        filout = open(p_filout, "w")
        l_perf = ["CV", "Train", "Test"]
        l_criteria = list(self.d_perf["CV"].keys())
        
        filout.write("\t%s\n"%("\t".join(l_criteria)))
        
        for perf in l_perf:
            filout.write("%s-%s\t%s\n"%(perf, name_combine, "\t".join(["%.2f"%(self.d_perf[perf][criteria]) for criteria in l_criteria])))
        filout.close() 