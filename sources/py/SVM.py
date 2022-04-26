from os import path
import numpy as np
from sklearn.model_selection import StratifiedKFold, KFold, RepeatedStratifiedKFold
from copy import deepcopy
import joblib
import ghost
from sklearn.svm import SVC
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV

import ML_toolbox
import pathFolder

class SVM:
    
    def __init__(self, p_train, p_test, p_aff, n_foldCV, kernel, pr_out, ghost = 0):
        self.p_train = p_train
        self.p_test = p_test
        self.p_aff = p_aff

        # optimization with ghost
        self.ghost = ghost
        self.l_ghost_threshold = np.round(np.arange(0.05,0.55,0.05),2)
        self.n_foldCV = n_foldCV
        self.verbose = 0
        self.test = 0
        self.force_run = 0
        self.typeModel = "classification"
        self.kernel = kernel #['linear', 'poly', 'rbf', 'sigmoid']
        self.pr_out = pr_out

        self.d_model = {}
        # grid optimization - use a test criteria to reduce the grid size for testing
        self.test = 0

        if self.ghost == 1:
            type_SVM = self.kernel + "_ghost"
        else:
            type_SVM = self.kernel

        if self.force_run == 1:
            self.pr_results = pathFolder.createFolder(self.pr_out + "SVM-" + type_SVM + "/", clean = 1)
        else:
            self.pr_results = pathFolder.createFolder(self.pr_out + "SVM-" + type_SVM + "/")

        self.n_jobs = 5
        self.random_state = 1
        # use to limit optimization
        if self.test == 0:
            # define search space
            self.params_space = dict()
            self.params_space['C'] = [0.1, 1, 2, 10, 100]
            self.params_space['gamma'] = [1, 0.1, 0.01, 0.001, 0.0001]
            self.params_space['kernel'] = [kernel]
        else:
            self.params_space = dict()
            self.params_space['C'] = [0.1, 1]
            self.params_space['gamma'] = [1, 0.1]
            self.params_space['kernel'] = [kernel]

    def loadSet(self):
    
        d_train = ML_toolbox.loadSet(self.p_train, variableToPredict="Aff")
        self.dataset_train = d_train["dataset"]
        self.aff_train = d_train["aff"]
        self.id_train = d_train["id"]
        self.nb_desc_input = d_train["nb_desc_input"]
        self.l_features_train = d_train["features"]

        # descriptor test set 
        d_test = ML_toolbox.loadSet(self.p_test, self.l_features_train, variableToPredict = "Aff")
        self.dataset_test = d_test["dataset"]
        self.aff_test = d_test["aff"]
        self.id_test = d_test["id"]
        self.nb_desc_test = d_test["nb_desc_input"]

        # control the number of desc are the same
        if self.nb_desc_input != self.nb_desc_test:
            print("===CHECK SIZE TRAIN AND TEST===")
            STOPHERE # crash script

    def run_SVM(self, CV = 0, **kwargs):
        """
        Run modeling 
        Tune SVM -> https://machinelearningmastery.com/scikit-optimize-for-hyperparameter-tuning-in-machine-learning/
        """
        # add ghost here
        if self.ghost == 1:
            type_SVM = self.kernel + "_ghost"
        else:
            type_SVM = self.kernel

        if CV == 1:
            pr_out = self.pr_results + "CV_"
        else:
            pr_out = self.pr_results

        # short cut with the model already computed and save
        p_model = pr_out + "model.joblib"
        
        # do not check if CV model exist
        if path.exists(p_model) and self.force_run == 0 and CV == 0:
            self.d_model[type_SVM] = joblib.load(p_model)

        else:
            # define evaluation -> repetition 3 times the modeling
            svm = SVC(probability=True)
            # define the search
            svm_random = RandomizedSearchCV(estimator = svm, param_distributions = self.params_space, cv = self.n_foldCV, verbose=2, random_state=self.random_state, n_jobs = self.n_jobs)
            svm_random.fit(self.dataset_train, self.aff_train)
            best_random = svm_random.best_estimator_

            # save the model
            joblib.dump(best_random, p_model)

            # do not save CV on the save as RF
            if CV == 0:
                self.d_model[type_SVM] = joblib.load(p_model)

            else:
                self.d_model["CV_" + type_SVM] = joblib.load(p_model)

    def apply_model(self, a_set, id_set, type_SVM, name_out = "", w=0): 
       
        y_pred = self.d_model[type_SVM].predict_proba(a_set)

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

    def CrossValidation(self):
        """
        Realize a cross validation in N folds
        args: - type of RF (balanced/umbalanced) 
        """

        self.loadSet()

        if self.ghost == 1:
            type_SVM = self.kernel + "_ghost"
        else:
            type_SVM = self.kernel

        d_CV = {}
        if self.typeModel == "classification":
            kfold = StratifiedKFold(n_splits=self.n_foldCV, shuffle=True, random_state=self.random_state)
        else:
            kfold = KFold(n_splits=self.n_foldCV, shuffle=True, random_state=self.random_state)
        # load parameters
        
        d_train = deepcopy(self.dataset_train)
        Y_train = deepcopy(self.aff_train)
        id_train = deepcopy(self.id_train)

        
        for train, test in kfold.split(d_train, Y_train):
            
            # update the class
            self.dataset_train = d_train[train]
            self.aff_train = Y_train[train]
            self.id_train = list(np.array(id_train)[list(train)])
            

            self.dataset_test = d_train[test]
            self.aff_test = Y_train[test]
            self.id_test = list(np.array(id_train)[list(test)])

            # optimize model
            self.run_SVM(CV=1)

            # apply prediction on the test
            y_pred = self.apply_model(self.dataset_test, self.id_test, "CV_" + type_SVM)
            #if self.ghost == 1:
            #    y_pred_train = self.apply_model(self.dataset_train, self.id_train, "CV_" + type_SVM)
            #    y_pred_train_ghost = [pred[1] for pred in y_pred_train]
            #    self.threshold_ghost = ghost.optimize_threshold_from_predictions(self.aff_train, y_pred_train_ghost, self.l_ghost_threshold, ThOpt_metrics = 'Kappa') 
            #else:
            self.threshold_ghost = 0.5

            d_pref_test = ML_toolbox.performance(self.aff_test, y_pred, typeModel = "classification", th_prob=self.threshold_ghost)
            
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
        self.run_SVM(CV = 0)
        
        if not "d_perf" in self.__dict__:
            self.d_perf = {}

        # trainning set
        # add ghost here
        if self.ghost == 1:
            type_SVM = self.kernel + "_ghost"
        else:
            type_SVM = self.kernel
        y_train_pred = self.apply_model(self.dataset_train, self.id_train, type_SVM, "train_pred.csv", w=1)
        if self.ghost == 1:
            y_train_pred_ghost = [pred[1] for pred in y_train_pred]                
            threshold_ghost = ghost.optimize_threshold_from_predictions(self.aff_train, y_train_pred_ghost, self.l_ghost_threshold, ThOpt_metrics = 'Kappa') 
        else:
            threshold_ghost = 0.5

        # train set
        self.ghost_treshold = threshold_ghost
        self.d_perf["Train"] = ML_toolbox.performance(self.aff_train, y_train_pred, typeModel = "classification", th_prob = self.ghost_treshold)
            
        # test set
        y_pred = self.apply_model(self.dataset_test, self.id_test, type_SVM, "test_pred.csv", w=1)
        self.d_perf["Test"] = ML_toolbox.performance(self.aff_test, y_pred, typeModel = "classification", th_prob = self.ghost_treshold) 

    def combineResults(self):
           
        if not "d_perf" in self.__dict__:
            print("No data to write")
            return 

        if self.ghost == 1:
            type_SVM = self.kernel + "_ghost"
        else:
            type_SVM = self.kernel
        
        p_filout = self.pr_results + "combined_perf.csv"

        filout = open(p_filout, "w")
        l_perf = ["CV", "Train", "Test"]
        l_criteria = list(self.d_perf["CV"].keys())
        
        filout.write("\t%s\n"%("\t".join(l_criteria)))
        
        for perf in l_perf:
            filout.write("%s-SVM%s\t%s\n"%(perf, type_SVM, "\t".join(["%.2f"%(self.d_perf[perf][criteria]) for criteria in l_criteria])))
        filout.close() 

    def run_main(self):
    
        if self.ghost == 1:
            type_SVM = self.kernel + "_ghost"
        else:
            type_SVM = self.kernel

        p_combined_results = self.pr_out + "SVM-" + type_SVM + "/" + "combined_perf.csv"
        if path.exists(p_combined_results):
            return  

        self.TrainTestPrediction()
        self.CrossValidation()
        self.combineResults()