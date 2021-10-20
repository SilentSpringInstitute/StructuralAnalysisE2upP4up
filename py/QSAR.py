import RandomForest
import pathFolder
import runExternal
import toolbox
import DNN
import SVM

from os import path, listdir, rename, remove
from re import search
from numpy import mean, std
from copy import deepcopy
from shutil import copyfile
from random import shuffle, uniform


class QSAR:
    def __init__(self, p_desc_clean, p_desc_origin, p_AC50, p_AC50_origin, p_sim_matrix, pr_out, nb_repetition = 1, nb_sample = 0, n_foldCV = 10, rate_active = 0.30, rate_splitTrainTest = 0.15):
        self.p_desc = p_desc_clean
        self.p_desc_orign = p_desc_origin
        self.p_AC50 = p_AC50
        self.p_AC50_orign = p_AC50_origin
        self.p_sim_matrix = p_sim_matrix
        self.pr_out = pr_out
        self.nb_sample = nb_sample
        self.nb_repetition = nb_repetition
        self.n_foldCV = n_foldCV
        self.rate_splitTrainTest = rate_splitTrainTest

        #Add variable sampling 
        self.rate_active = rate_active

        # environement
        self.force_run = 0

    def runQSARClassUnderSamplingAllSet(self):

        l_run = list(range(1, self.nb_repetition + 1))
        shuffle(l_run)
        
        for run in l_run:
            pr_run = pathFolder.createFolder(self.pr_out + str(run) + "/")
            self.pr_run = pr_run
        
            # check applicability model
            pr_AD = pathFolder.createFolder(self.pr_run + "AD/")
            l_sample = list(range(1, self.nb_sample + 1))
            shuffle(l_sample)

            for i in l_sample:
                pr_run = self.pr_run + str(i) + "/"
                pathFolder.createFolder(pr_run)

                # prepare dataset => split train / test
                self.prepTrainSetforUnderSampling(pr_run)

                # check AD
                pr_AD_run = pathFolder.createFolder(pr_AD + str(i) + "/")
                if len(listdir(pr_AD_run)) < 6:
                    runExternal.AD(self.p_train, self.p_test, pr_AD_run)

                # build QSAR
                self.buildQSARs(pr_run)

            # merge results
            self.mergeModelSampled()

    def runQSARClassNoSampling(self):
        l_run = list(range(1, self.nb_repetition + 1))
        shuffle(l_run)

        # reduce number of run
        for run in l_run:
            pr_run = pathFolder.createFolder(self.pr_out + str(run) + "/")
            self.pr_run = pr_run

            # define train and test set here
            self.prepSplitTrainTestSet()
            self.p_train = self.p_trainGlobal

            # build QSAR
            self.buildQSARs(pr_run)

        # combine repetition 
        self.combineRepetitionNoSampled()

    def runQSARClassUnderSamplingTrain(self):

        l_run = list(range(1, self.nb_repetition + 1))
        shuffle(l_run)
        
        for run in l_run:
            pr_run = pathFolder.createFolder(self.pr_out + str(run) + "/")
            self.pr_run = pr_run

            # define train and test set here
            self.prepSplitTrainTestSet()

            l_sample = list(range(1, self.nb_sample + 1))
            shuffle(l_sample)
            
            for sample in l_sample:
                # redefine rate undersampling
                pr_sample = self.pr_run + str(sample) + "/"
                pathFolder.createFolder(pr_sample)
                
                # prepare dataset => split train / test
                self.prepTrainSetforUnderSampling(pr_sample, 0)
                
                # build QSAR
                self.buildQSARs(pr_sample)

            # merge results
            self.mergeModelSampled()
        
        # combine repetition 
        self.combineRepetitionSampling()

    def runQSARReg(self, corcoef, maxQuantile):
        l_run = list(range(1, self.nb_repetition + 1))
        shuffle(l_run)

        for i in l_run:
            pr_run = pathFolder.createFolder(self.pr_out + str(i) + "/")
            if not path.exists(pr_run + "trainSet.csv"):
                runExternal.prepDataQSARReg(self.p_desc, self.p_AC50, pr_run, corcoef, maxQuantile, self.rate_splitTrainTest,  typeAff="All", logaff=0, nbNA = 10)
            
            if not path.exists(pr_run + "perfTest.csv"): # control the perf test file is not present
                runExternal.runQSARReg(pr_run + "trainSet.csv", pr_run + "testSet.csv", "0", pr_run, self.n_foldCV)

    def mergeRegResults(self):

        d_result = {}
        pr_result = pathFolder.createFolder(self.pr_out + "mergeResults/")
        l_pr_run = listdir(self.pr_out)

        for pr_run in l_pr_run:
            if pr_run == "mergeResults" or search("csv", pr_run) or pr_run == "Cleaned_Data":
                continue
            else:
                d_result[pr_run] = {}
                p_perfTrain = self.pr_out + pr_run + "/perfTrain.csv"
                p_perfCV = self.pr_out + pr_run + "/perfCV.csv"
                p_perfTest = self.pr_out + pr_run + "/perfTest.csv"

                d_result[pr_run]["train"] =  toolbox.loadMatrix(p_perfTrain, sep = ",")
                d_result[pr_run]["test"] =  toolbox.loadMatrix(p_perfTest, sep = ",")
                d_result[pr_run]["CV"] =  toolbox.loadMatrix(p_perfCV, sep = ",")

        p_filout = pr_result + "results.csv"
        filout = open(p_filout, "w")
        for run in d_result.keys():
            filout.write("Run: %s\n"%(run))
            filout.write("\tCV\t\t\t\t\t\tTrain\t\t\t\t\t\tTest\t\t\t\n")
            filout.write("ML\tR2\tR0\tMAE\tcor\tRMSEP\t\tR2\tR0\tMAE\tcor\tRMSEP\t\tR2\tR0\tMAE\tcor\tRMSEP\n")
            for ML in d_result[run]["CV"].keys():
                filout.write("%s\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\n"%(ML, d_result[run]["CV"][ML]["R2"], d_result[run]["CV"][ML]["R02"], d_result[run]["CV"][ML]["MAE"], d_result[run]["CV"][ML]["r"], d_result[run]["CV"][ML]["RMSEP"], d_result[run]["train"][ML]["R2"], d_result[run]["train"][ML]["R02"], d_result[run]["train"][ML]["MAE"], d_result[run]["train"][ML]["r"], d_result[run]["train"][ML]["RMSEP"], d_result[run]["test"][ML]["R2"], d_result[run]["test"][ML]["R02"], d_result[run]["test"][ML]["MAE"], d_result[run]["test"][ML]["r"], d_result[run]["test"][ML]["RMSEP"]))   
            filout.write("\n")

    def prepSplitTrainTestSet(self):


        # define train and test set
        p_train = self.pr_run + "trainGlobal.csv"
        p_test = self.pr_run + "test.csv"

        if path.exists(p_train) and path.exists(p_test):
            self.p_trainGlobal = p_train
            self.p_test = p_test
            return 

        # prep descriptor with classes without ratio of active => define train global
        runExternal.prepDataQSAR(self.p_desc, self.p_AC50, 0, self.pr_run)
        runExternal.SplitTrainTest(self.pr_run + "desc_Class.csv", self.pr_run, self.rate_splitTrainTest)
        # rename train in global train
        rename(self.pr_run + "train.csv", self.pr_run + "trainGlobal.csv")
        self.p_trainGlobal = p_train
        self.p_test = p_test

    def prepTrainSetforUnderSampling(self, pr_run, splitRatio=""):

        #Add variable sampling 
        if type(self.rate_active) == list:
            self.sampling_rate = uniform(self.rate_active[0], self.rate_active[1])
        else:
            self.sampling_rate = self.rate_active

        # Case where only the train set is created
        if splitRatio == 0: 
            # define train and test set
            p_train = pr_run + "train.csv"

            if path.exists(p_train):
                self.p_train = p_train
            else:
                # prep descriptor with classes without ratio of active
                runExternal.prepDataQSAR(self.p_trainGlobal, self.p_AC50, self.sampling_rate, pr_run)
                # to apply similar formating and create train.csv
                runExternal.SplitTrainTest(pr_run + "desc_Class.csv", pr_run, 0)
                
                # Rename file in train set file
                self.p_train = p_train

        else:
            # define train and test set
            p_train = pr_run + "train.csv"
            p_test = pr_run + "test.csv"

            if path.exists(p_train) and path.exists(p_test):
                self.p_train = p_train
                self.p_test = p_test
                return 

            # prep descriptor with classes
            if not path.exists(pr_run + "desc_Class.csv"):
                runExternal.prepDataQSAR(self.p_desc, self.p_AC50, self.sampling_rate, pr_run)

            runExternal.SplitTrainTest(pr_run + "desc_Class.csv", pr_run, self.rate_splitTrainTest)
            self.p_train = p_train
            self.p_test = p_test
    
    def buildQSARs(self, pr_run):
        """Build QSAR models
        arg: - path folder results
             (- check self.force_run to rebuild the modeling)
        return: None
            write output in the folder 
        """
        # flush the result folder in case of forcing
        if self.force_run == 1:
            pathFolder.cleanFolder(pr_run, [self.p_train, self.p_test])
        
        # run AD
        self.computeAD(pr_run)
                
        # classic machine learning with R
        # only used for NN singleton and CART
        ############################
        p_perfTrain = pr_run + "perfTrain.csv"
        p_perfCV = pr_run + "perfCV.csv"
        p_perfTest = pr_run + "perfTest.csv"

        if not path.exists(p_perfCV) or not path.exists(p_perfTrain) or not path.exists(p_perfTest):
            runExternal.runRQSAR(self.p_train, self.p_test, self.n_foldCV, pr_run)


        # DNN with keras / tensorflow
        # self, p_train, p_test, p_aff, pr_out, n_foldCV, typeModel
        c_DNN = DNN.DNN(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "classification", 0, pr_run)
        c_DNN.run_main()

        # DNN with a ghost optimization
        c_DNN = DNN.DNN(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "classification", 1, pr_run)
        c_DNN.run_main()

        # use the ghost approach for RF and DNN
        # define a class specific to this approch for testing
        # need to be integrated in the DNN class and make new class for RF and 
        
        # run SVM from python script
        ## kernel to test ['linear', 'poly', 'rbf', 'sigmoid']
        #c_SVM = SVM.SVM(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "linear", pr_run) 
        #c_SVM.run_main()

        #c_SVM = SVM.SVM(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "poly", pr_run) 
        #c_SVM.run_main()

        c_SVM = SVM.SVM(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "rbf", pr_run) 
        c_SVM.run_main()

        #c_SVM = SVM.SVM(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "sigmoid", pr_run) 
        #c_SVM.run_main()
       
        # classic RF
        ##############
        c_RF = RandomForest.RandomForest(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "RF_py", 0, pr_run)
        c_RF.run_main()
        
        # with ghost
        c_RF = RandomForest.RandomForest(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "RF_py", 1, pr_run)
        c_RF.run_main()

        # balanced RF
        c_RFb = RandomForest.RandomForest(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "RF_balanced", 0, pr_run)
        c_RFb.run_main()

        c_RFb = RandomForest.RandomForest(self.p_train, self.p_test, self.p_AC50, self.n_foldCV, "RF_balanced", 1, pr_run)
        c_RFb.run_main()

        # summarize the run
        self.writeSumFile(pr_run)    

    def mergeModelSampled(self):

        pr_QSAR_average = pathFolder.createFolder(self.pr_run + "Merge_results/")
        pr_QSAR_proba = self.pr_run + "Merge_probRF/"# proba of prediction merged
        pr_QSAR_desc_involved = self.pr_run + "Merge_involvedDesc/"
        pr_AD =  pathFolder.createFolder(self.pr_run + "Merge_AD/")

        self.mergeSamplingResults(pr_QSAR_average)
        self.mergeSamplingProbaRF(self.p_AC50_orign, pr_QSAR_proba)
        self.mergeSamplingInvolvedDesc("RF", 10, pr_QSAR_desc_involved)
        self.mergeSamplingInvolvedDesc("LDA", 10, pr_QSAR_desc_involved)
        self.mergeSamplingAD(pr_AD)

    def mergeSamplingResults(self, pr_av):
        """
        Merge results from several ML
        args: - folder for results
        return: None put results in the folder
        """

        p_filout = pr_av + "average_perf.csv"
        if path.exists(p_filout) and self.force_run == 0:
            return 
        l_criteria = ["Acc", "Sp", "Se", "MCC"]
        l_dataset = ["CV", "train", "test"]
        
        # output to write
        d_result = {}
        for dataset in l_dataset:
            d_result[dataset] = {}

        # performance by ML
        d_perf = {}
        for criteria in l_criteria:
            d_perf[criteria] = []

        l_pr_run = listdir(self.pr_run)
        for pr_run in l_pr_run:
            if search("Merge", pr_run) or search("AD", pr_run):
                continue

            p_perfCV = self.pr_run + pr_run + "/perfCV.csv"
            p_perfTrain = self.pr_run + pr_run + "/perfTrain.csv"
            p_perfTest = self.pr_run + pr_run + "/perfTest.csv"

            try:
                M_CV = toolbox.loadMatrix(p_perfCV, sep=",")
                M_train = toolbox.loadMatrix(p_perfTrain, sep=",")
                M_test = toolbox.loadMatrix(p_perfTest, sep=",")
            except:
                continue

            # check if DNN here
            p_dnn = self.pr_run + pr_run + "/DNN/combined_perf.csv"
            if path.exists(p_dnn):
                d_DNN = toolbox.loadMatrix(p_dnn)
                M_CV["DNN"] = {}
                M_train["DNN"] = {}
                M_test["DNN"] = {}
                for criteria in l_criteria:
                    M_CV["DNN"][criteria] = 0.0
                    M_train["DNN"][criteria] = 0.0
                    M_test["DNN"][criteria] = 0.0
                
                M_CV["DNN"]["Acc"] =  d_DNN["CV-DNN"]["Acc"]
                M_CV["DNN"]["MCC"] =  d_DNN["CV-DNN"]["MCC"]
                M_CV["DNN"]["Sp"] =  d_DNN["CV-DNN"]["Sp"]
                M_CV["DNN"]["Se"] =  d_DNN["CV-DNN"]["Se"]
                M_train["DNN"]["Acc"] =  d_DNN["Train-DNN"]["Acc"]
                M_train["DNN"]["MCC"] =  d_DNN["Train-DNN"]["MCC"]
                M_train["DNN"]["Se"] =  d_DNN["Train-DNN"]["Se"]
                M_train["DNN"]["Sp"] =  d_DNN["Train-DNN"]["Sp"]
                M_test["DNN"]["Acc"] =  d_DNN["Test-DNN"]["Acc"]
                M_test["DNN"]["MCC"] =  d_DNN["Test-DNN"]["MCC"]
                M_test["DNN"]["Sp"] =  d_DNN["Test-DNN"]["Sp"]
                M_test["DNN"]["Se"] =  d_DNN["Test-DNN"]["Se"]

            l_ML = list(M_CV.keys())

            #build results
            if not l_ML[0] in list(d_result["CV"].keys()):
                for ML in l_ML:
                    d_result["CV"][ML] = deepcopy(d_perf)
                    d_result["train"][ML] = deepcopy(d_perf)
                    d_result["test"][ML] = deepcopy(d_perf)

            for ML in l_ML:
                for criteria in l_criteria:
                    d_result["CV"][ML][criteria].append(float(M_CV[ML][criteria]))
                    d_result["train"][ML][criteria].append(float(M_train[ML][criteria]))
                    d_result["test"][ML][criteria].append(float(M_test[ML][criteria]))


        d_out = deepcopy(d_result)
        for dataset in l_dataset:
            for ML in l_ML:
                for criteria in l_criteria:
                    AV = round(mean(d_result[dataset][ML][criteria]),3)
                    SD = round(std(d_result[dataset][ML][criteria]),3)
                    d_out[dataset][ML][criteria] = [AV, SD]


        # write result
        filout = open(p_filout, "w")
        for dataset in l_dataset:
            filout.write(str(dataset) + "\n")
            filout.write("\t" + "\t".join(["M-" + str(c) + "\t" + "SD-" + str(c) for c in l_criteria]) + "\n")
            for ML in l_ML:
                filout.write(ML)
                for criteria in l_criteria:
                    filout.write("\t" + str(d_out[dataset][ML][criteria][0]) + "\t" + str(d_out[dataset][ML][criteria][1]))
                filout.write("\n")
            filout.write("\n")
        filout.close()

    def mergeSamplingAD(self, pr_AD_merged):

        # use a short cut
        l_AD_files = listdir(pr_AD_merged) 
        if len(l_AD_files) > 5 and self.force_run == 0:
            return

        l_pr_run = listdir(self.pr_run)

        d_Zscore_train = {}
        d_Zscore_test = {}
        for pr_run in l_pr_run:
            try: int(pr_run)
            except:continue
            p_Zscore_train = self.pr_run + pr_run + "/AD/descPCA_based/AD_Train_zscore.csv"
            p_Zscore_test = self.pr_run + pr_run + "/AD/descPCA_based/AD_Test_zscore.csv"
            
            if not path.exists(p_Zscore_train) and not path.exists(p_Zscore_test):
                print("Check - AD:", p_Zscore_train)
                print("Check - AD:", p_Zscore_test)
                continue
                
            d_train = toolbox.loadMatrix(p_Zscore_train, sep = ",")
            d_test = toolbox.loadMatrix(p_Zscore_test, sep = ",")

            for chem_train in d_train.keys():
                if not chem_train in list(d_Zscore_train.keys()):
                    d_Zscore_train[chem_train] = []    
                d_Zscore_train[chem_train].append(float(d_train[chem_train]["Zscore"]))

            for chem_test in d_test.keys():
                if not chem_test in list(d_Zscore_test.keys()):
                    d_Zscore_test[chem_test] = []
                d_Zscore_test[chem_test].append(float(d_test[chem_test]["Zscore"]))
        
        if d_Zscore_test == {}:
            print("No AD found")
            return 

        # write zscore
        p_train = pr_AD_merged + "Zscores_train.csv"
        ftrain = open(p_train, "w")
        ftrain.write("CASRN\tZscore\n")
        for chem in d_Zscore_train.keys():
            ftrain.write("%s\t%s\n"%(chem, mean(d_Zscore_train[chem])))
        ftrain.close()

        p_test = pr_AD_merged + "Zscores_test.csv"
        ftest = open(p_test, "w")
        ftest.write("CASRN\tZscore\n")
        for chem in d_Zscore_test.keys():
            ftest.write("%s\t%s\n"%(chem, mean(d_Zscore_test[chem])))
        ftest.close()

        # draw histogram for AD and add summary and PCA
        runExternal.mergeADs(p_train, p_test, self.p_desc, pr_AD_merged)

    def mergeSamplingProbaRF(self, p_AC50, pr_prob):
        # need to change R scripts for 

        if path.exists(pr_prob) and len(listdir(pr_prob)) > 10:
            return 

        d_prob = {}
        l_dataset = ["CV", "train", "test"]

        # check if RF folder exist
        l_pr_run = listdir(self.pr_run)
        i = 0
        imax = len(l_pr_run)
        while i < imax: 
            try: int(l_pr_run[i])
            except: 
                del l_pr_run[i]
                imax = imax - 1
                continue

            # check if RF is computed
            if not path.exists("%s%s/RFclass/"%(self.pr_run, l_pr_run[i])):
                del l_pr_run[i]
                imax = imax - 1
                continue
            i = i + 1

        if len(l_pr_run) == 0:
            print("No RF computed")
            return 


        # load AC50
        d_AC50 = toolbox.loadMatrix(p_AC50)

        for pr_run in l_pr_run:

            # check if RF is computed
            if not path.exists("%s%s/RFclass/"%(self.pr_run, pr_run)):
                print("Error run: %s"%(pr_run))
                continue

            d_prob[pr_run] = {}

            # CV
            p_CV = "%s%s/RFclass/PerfRFClassCV%s.txt"%(self.pr_run, pr_run, self.n_foldCV)
            d_CV = toolbox.loadMatrix(p_CV, sep = "\t")
            d_prob[pr_run]["CV"] = d_CV

            # train
            p_train = self.pr_run + pr_run + "/RFclass/classTrain.csv"
            d_train = toolbox.loadMatrix(p_train, sep = ",")
            d_prob[pr_run]["train"] = d_train

            # test
            p_test = self.pr_run + pr_run + "/RFclass/classTest.csv"
            d_test = toolbox.loadMatrix(p_test, sep = ",")
            d_prob[pr_run]["test"] = d_test


        d_av = {}
        d_av["CV"] = {}
        d_av["train"] = {}
        d_av["test"] = {}

        for run in l_pr_run:

            # CV
            for chem in d_prob[run]["CV"].keys():
                if not chem in list(d_av["CV"].keys()):
                    d_av["CV"][chem] = {}
                    try: d_av["CV"][chem]["LogAC50"] = d_AC50[chem]["LogAC50"] 
                    except: d_av["CV"][chem]["LogAC50"] = d_AC50[chem]["Aff"] 
                    d_av["CV"][chem]["predicted"] = []
                d_av["CV"][chem]["predicted"].append(float(d_prob[run]["CV"][chem]["Predict"]))
            
            # train
            for chem in d_prob[run]["train"].keys():
                if not chem in list(d_av["train"].keys()):
                    d_av["train"][chem] = {}
                    try: d_av["train"][chem]["LogAC50"] = d_AC50[chem]["LogAC50"] 
                    except: d_av["train"][chem]["LogAC50"] = d_AC50[chem]["Aff"]
                    d_av["train"][chem]["predicted"] = []
                d_av["train"][chem]["predicted"].append(float(d_prob[run]["train"][chem]["x"]))

            # test
            for chem in d_prob[run]["test"].keys():
                if not chem in list(d_av["test"].keys()):
                    d_av["test"][chem] = {}
                    try:d_av["test"][chem]["LogAC50"] = d_AC50[chem]["LogAC50"] 
                    except:d_av["test"][chem]["LogAC50"] = d_AC50[chem]["Aff"]
                    d_av["test"][chem]["predicted"] = []
                d_av["test"][chem]["predicted"].append(float(d_prob[run]["test"][chem]["x"]))


            
        # train
        if not path.exists(pr_prob):
            pr_prob = pathFolder.createFolder(pr_prob)
        p_prob_train = pr_prob + "Prob_train"
        f_prob_train = open(p_prob_train, "w")
        f_prob_train.write("ID\tMpred\tSDpred\tReal\n")
        for IDtrain in list(d_av["train"].keys()):
            f_prob_train.write("%s\t%.3f\t%.3f\t%s\n"%(IDtrain, mean(d_av["train"][IDtrain]["predicted"]), std(d_av["train"][IDtrain]["predicted"]), d_av["train"][IDtrain]["LogAC50"]))
        f_prob_train.close()
        
        runExternal.plotAC50VSProb(p_prob_train)

        # CV
        p_prob_CV = pr_prob + "Prob_CV"
        f_prob_CV = open(p_prob_CV, "w")
        f_prob_CV.write("ID\tMpred\tSDpred\tReal\n")
        for IDCV in list(d_av["CV"].keys()):
            f_prob_CV.write("%s\t%.3f\t%.3f\t%s\n"%(IDCV, mean(d_av["CV"][IDCV]["predicted"]), std(d_av["CV"][IDCV]["predicted"]), d_av["CV"][IDCV]["LogAC50"]))
        f_prob_CV.close()
        runExternal.plotAC50VSProb(p_prob_CV)

        # test
        p_prob_test = pr_prob + "Prob_test"
        f_prob_test = open(p_prob_test, "w")
        f_prob_test.write("ID\tMpred\tSDpred\tReal\n")
        for IDtest in list(d_av["test"].keys()):
            f_prob_test.write("%s\t%.3f\t%.3f\t%s\n"%(IDtest, mean(d_av["test"][IDtest]["predicted"]), std(d_av["test"][IDtest]["predicted"]), d_av["test"][IDtest]["LogAC50"]))
        f_prob_test.close()
        runExternal.plotAC50VSProb(p_prob_test)

    def mergeSamplingInvolvedDesc(self, ML, nbdesc, pr_involvedDesc):


        d_importance = {}
        pr_out = pr_involvedDesc + ML + "/"
        p_desc_importance = pr_out + "Av_importance"
        if path.exists(p_desc_importance):
            return

        # check if file exist
        l_pr_run = listdir(self.pr_run)
        i = 0
        imax = len(l_pr_run)
        while i < imax: 
            try: int(l_pr_run[i])
            except: 
                del l_pr_run[i]
                imax = imax - 1
                continue

            # check if model is computed
            if not path.exists("%s%s/%sclass/"%(self.pr_run, l_pr_run[i], ML)):
                del l_pr_run[i]
                imax = imax - 1
                continue
            i = i + 1

        if len(l_pr_run) == 0:
            print("No %s computed"%(ML))
            return 


        l_pr_run = listdir(self.pr_run)
        for run in l_pr_run:
            if search("Merge", run) or search("AD", run) or search("Cleaned", run) or search("desc", run) or search("csv", run) or search("models", run) or search("mergeDNN", run):
                continue

            if not run in list(d_importance.keys()):
                d_importance[run] = {}

                p_desc_involved = self.pr_run + run + "/" + str(ML) + "class/ImportanceDesc"
                d_desc_involved = toolbox.loadMatrix(p_desc_involved , sep="\t")
                d_importance[run] = d_desc_involved


        # global importance 
        if not path.exists(pr_out):
            pathFolder.createFolder(pr_out)
        f_desc_importance = open(p_desc_importance, "w")
        f_desc_importance.write("Desc\tRun\tval\n")

        l_desc = list(d_importance["1"].keys())
        
        for desc in l_desc:
            for run in d_importance.keys():
                try: w = desc + "\t" + str(run) + "\t" + str(d_importance[run][desc]["x"]) + "\n"
                except: w = desc + "\t" + str(run) + "\t0.0\n"
                f_desc_importance.write(w)
        f_desc_importance.close()
        
        runExternal.runImportanceDesc(p_desc_importance, nbdesc)

    def extractModels(self, pr_results, ML):

        pr_out = pathFolder.createFolder(pr_results)

        l_pr_run = listdir(self.pr_run)
        l_model = listdir(pr_out)
        if len(l_model) > 2:
            return pr_out

        for pr_run in l_pr_run:
            if search("Merge", pr_run) or search("AD", pr_run) or search("Cleaned", pr_run) or search("desc", pr_run) or search("csv", pr_run) or search("models", pr_run) or search("mergeDNN", pr_run):
                continue

            if search("SVM", ML):
                p_model = "%s%s/SVMclass_%s/model.RData"%(self.pr_run, pr_run, ML.split("-")[-1])
            elif search("DNN", ML):
                l_DNNfile = listdir("%s%s/DNNclass/"%(self.pr_run, pr_run))
                p_model = "out"
                for DNNfile in l_DNNfile:
                    if search(".h5", DNNfile):
                        p_model = "%s%s/DNNclass/%s"%(self.pr_run, pr_run, DNNfile)
                        break
            else:
                p_model = "%s%s/%sclass/model.RData"%(self.pr_run, pr_run, ML)
            if path.exists(p_model) and not path.exists(pr_out + pr_run + ".RData"):
                if p_model[-5:] == "RData":
                    copyfile(p_model, pr_out + pr_run + ".RData")
                else:
                    copyfile(p_model, pr_out + pr_run + ".h5")
        
        return pr_out
        
    def writeSumFile(self, pr_out):
        """Function use to wirte summary for the model"""

        if not "sampling_rate" in self.__dict__ and type(self.rate_active) == float:
            self.sampling_rate = self.rate_active

        # write summary of the QSAR modeling
        p_fsum = pr_out + "QSAR.sum"
        fsum = open(p_fsum, "w")
        fsum.write("Descriptor file prepared: %s\nDescriptor file origine: %s\nClass file prepared: %s\nClass file original: %s\nfolder out: %s\nRepetition: %s\nNb sample: %s\nFold cross_validation: %s\nActive rate: %s\nSplit train-test set: %s\n"%(self.p_desc, self.p_desc_orign, self.p_AC50, self.p_AC50_orign, pr_out, self.nb_repetition, self.nb_sample, self.n_foldCV, self.sampling_rate, self.rate_splitTrainTest))
        fsum.close()

    def computeAD(self, pr_out):
        """Compute applicability model with two ways, PCA with descriptors and similarit on MACCS tanimoto
        - folder with QSAR models
        - create AD with descriptors and similarity matrix
        """
        
        pr_AD = pathFolder.createFolder(pr_out + "AD/")

        # AD based on PCA descriptors
        pr_AD_desc = pathFolder.createFolder(pr_AD + "descPCA_based/")
        if not path.exists(pr_AD_desc + "PCA_color.png"): 
            runExternal.AD(self.p_train, self.p_test, pr_AD_desc)

        
        # AD based on similarity score
        pr_AD_sim = pathFolder.createFolder(pr_AD + "chem_similarity/")
        if not path.exists(pr_AD_sim + "CP_similarity_text.png"):
            # compute similarity matrix in the root folder.
            d_train = toolbox.loadMatrix(self.p_train, sep = ",")
            d_test = toolbox.loadMatrix(self.p_test, sep = ",")

            # create matrix with flag active vs inactive and test vs train
            p_matrix_chem = pr_AD_sim + "chem.csv"
            f_matrix_chem = open(p_matrix_chem, "w")
            f_matrix_chem.write("CASRN\tAff\tset\n")
            for chem in d_train.keys():
                f_matrix_chem.write("%s\t%s\t%s\n"%(chem, d_train[chem]["Aff"], "train"))
            
            for chem in d_test.keys():
                f_matrix_chem.write("%s\t%s\t%s\n"%(chem, d_test[chem]["Aff"], "test"))
            
            f_matrix_chem.close()

            runExternal.computeADBasedOnSimilarityMatrix(self.p_sim_matrix, p_matrix_chem, pr_AD_sim)        

    def combineRepetitionSampling(self):

        pr_out = pathFolder.createFolder(self.pr_out + "mergeRep/")
        p_filout = pr_out + "mergeRep.csv"

        l_rep  = listdir(self.pr_out)
        d_out = {}
        l_datasets = ["CV", "train", "test"] 
        for rep in l_rep:
            try: int(rep)
            except:continue
            p_results = "%s%s/Merge_results/average_perf.csv"%(self.pr_out, rep)
            if path.exists(p_results):
                f_results = open(p_results, "r")
                l_results = f_results.readlines()
                f_results.close()
            else:
                continue

            d_out[rep] = {}
            i = 0
            imax = len(l_results)
            while i < imax:
                line_results = l_results[i].strip()
                if line_results == "":
                    i = i + 1
                    continue
                if line_results in l_datasets:
                    dataset = line_results
                    d_out[rep][dataset] = {}
                    l_criteria = l_results[i + 1].strip().split("\t")
                    i = i + 2
                else:
                    l_val = line_results.strip().split("\t")
                    ML = l_val[0]
                    d_out[rep][dataset][ML] = {}
                    j = 0
                    jmax = len(l_criteria)

                    while j < jmax:
                        d_out[rep][dataset][ML][l_criteria[j]] = l_val[j + 1]
                        j = j + 1
                    i = i + 1


        filout = open(p_filout, "w")
        for run in d_out.keys():
            filout.write("Run: %s\n"%(run))
            filout.write("\tCV\t%s\tTrain\t%s\tTest\t%s\n"%("\t".join(["" for criteria in l_criteria]), "\t".join(["" for criteria in l_criteria]), "\t".join(["" for criteria in l_criteria])))
            filout.write("ML\t%s\t\t%s\t\t%s\n"%("\t".join(l_criteria),"\t".join(l_criteria),"\t".join(l_criteria)))
            for ML in d_out[run]["CV"].keys():
                filout.write("%s\t%s\t\t%s\t\t%s\n"%(ML, "\t".join(["%s"%(d_out[run]["CV"][ML][criteria]) for criteria in l_criteria]), "\t".join(["%s"%(d_out[run]["train"][ML][criteria]) for criteria in l_criteria]), "\t".join(["%s"%(d_out[run]["test"][ML][criteria]) for criteria in l_criteria]),))   
            filout.write("\n")
        filout.close()

        return 

    def combineRepetitionNoSampled(self):
        """
        Combine repetition model for analysis where the model is not a sampling
        """

        pr_out = pathFolder.createFolder(self.pr_out + "mergeRep/")
        p_filout = pr_out + "mergeRep.csv"

        l_rep  = listdir(self.pr_out)
        d_out = {}
        l_datasets = ["CV", "train", "test"] 
        l_criteria = []

        for rep in l_rep:
            try: int(rep)
            except:continue
            
            # format rep
            d_out[rep] = {}
            for dataset in l_datasets :d_out[rep][dataset] = {}

            p_cv = "%s%s/perfCV.csv"%(self.pr_out, rep)
            p_train = "%s%s/perfTrain.csv"%(self.pr_out, rep)
            p_test = "%s%s/perfTest.csv"%(self.pr_out, rep)
            if path.exists(p_cv) and path.exists(p_train) and path.exists(p_test):
                d_out[rep]["CV"] = toolbox.loadMatrix(p_cv, sep = ",")
                d_out[rep]["train"] = toolbox.loadMatrix(p_train, sep = ",")
                d_out[rep]["test"] = toolbox.loadMatrix(p_test, sep = ",")

                for criteria in d_out[rep]["train"][list(d_out[rep]["train"].keys())[0]]:
                    if not criteria in l_criteria and criteria != "ID":
                        l_criteria.append(criteria)

            l_folders = listdir(self.pr_out + rep)
            for folder in l_folders:
                p_combined = "%s%s/%s/combined_perf.csv"%(self.pr_out, rep, folder)
                if path.exists(p_combined):
                    d_temp = toolbox.loadMatrix(p_combined)
                    for k in d_temp.keys():
                        ML = k.split("-")[1]
                        dataset = k.split("-")[0]
                        if dataset != "CV":
                            dataset = dataset.lower()
                        d_out[rep][dataset][ML] = {}
                        d_out[rep][dataset][ML] = d_temp[k]

                        # add criteria parameter 
                        for criteria in d_out[rep][dataset][ML].keys():
                            if not criteria in l_criteria and criteria != "ID":
                                l_criteria.append(criteria)


        filout = open(p_filout, "w")
        for run in d_out.keys():
            filout.write("Run: %s\n"%(run))
            filout.write("\tCV\t%s\tTrain\t%s\tTest\t%s\n"%("\t".join(["" for criteria in l_criteria]), "\t".join(["" for criteria in l_criteria]), "\t".join(["" for criteria in l_criteria])))
            filout.write("ML\t%s\t\t%s\t\t%s\n"%("\t".join(l_criteria),"\t".join(l_criteria),"\t".join(l_criteria)))
            for ML in d_out[run]["CV"].keys():
                for criteria in l_criteria: 
                    if not criteria in list(d_out[run]["CV"][ML].keys()): d_out[run]["CV"][ML][criteria] = "NA" 
                    else: d_out[run]["CV"][ML][criteria] = "%.2f"%(float(d_out[run]["CV"][ML][criteria]))
                    if not criteria in list(d_out[run]["train"][ML].keys()): d_out[run]["train"][ML][criteria] = "NA"
                    else: d_out[run]["train"][ML][criteria] = "%.2f"%(float(d_out[run]["train"][ML][criteria]))
                    if not criteria in list(d_out[run]["test"][ML].keys()): d_out[run]["test"][ML][criteria] = "NA"
                    else: d_out[run]["test"][ML][criteria] = "%.2f"%(float(d_out[run]["test"][ML][criteria]))
                filout.write("%s\t%s\t\t%s\t\t%s\n"%(ML, "\t".join(["%s"%(d_out[run]["CV"][ML][criteria]) for criteria in l_criteria]), "\t".join(["%s"%(d_out[run]["train"][ML][criteria]) for criteria in l_criteria]), "\t".join(["%s"%(d_out[run]["test"][ML][criteria]) for criteria in l_criteria]),))   
            filout.write("\n")
        filout.close()