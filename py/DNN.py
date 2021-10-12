import toolbox
import pathFolder
import runExternal

from numpy import loadtxt, arange
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

#import tensorflow import keras
from tensorflow.keras.datasets import mnist 
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout 
from tensorflow.keras.optimizers import RMSprop 
import tensorflow as tf
import numpy as np 
from os import listdir, remove, path
from re import search
from copy import deepcopy

import ML_toolbox

class DNN:
    def __init__(self, p_train, p_test, p_aff, n_foldCV, typeModel, ghost, pr_root):
        self.p_train = p_train
        self.p_test = p_test
        self.p_aff = p_aff
        self.n_foldCV = n_foldCV
        self.verbose = 0
        self.test=0
        self.force_run = 0
        # ghost approach to optimize the prob threshold
        self.ghost = ghost

        # create folder
        if self.ghost == 1:
            type_opt = "_ghost"
        else:
            type_opt = ""

        pr_out = pathFolder.createFolder(pr_root + "DNN" + type_opt + "/")
        self.pr_out = pr_out

        self.pr_root = pr_root        

        # optimization
        self.typeModel = typeModel
        if typeModel == "classification":
            self.kernel_initializer = "random_normal"
            self.loss = "binary_crossentropy"
        else:
            self.kernel_initializer = "normal"
            self.loss = "mean_squared_error"

        
        self.optimizer = "adam"
        self.l_epochs = [50, 100, 120]
        self.l_batch_size = [32, 64, 128]
        self.l_dense_layer = [3, 4, 5]
        self.l_dense_layer = [3]
        self.l_dense_candidate = [50, 25, 20, 10, 1]
        self.l_activation = ["relu", "selu"]

        if self.test == 1:
            self.l_epochs = [50]
            self.l_batch_size = [32]
            self.l_dense_layer = [3]
            self.l_dense_candidate = [4,3,2,1]
            self.l_activation = ["relu"]

    def loadSet(self):
        
        # descriptor train set 
        d_train  = ML_toolbox.loadSet(self.p_train, variableToPredict = 1)
        self.dataset_train = d_train["dataset"]
        self.aff_train = d_train["aff"]
        self.id_train = d_train["id"]
        self.nb_desc_input = d_train["nb_desc_input"]

        # descriptor test set 
        d_test = ML_toolbox.loadSet(self.p_test, variableToPredict = 1)
        self.dataset_test = d_test["dataset"]
        self.aff_test = d_test["aff"]
        self.id_test = d_test["id"]
        self.nb_desc_test = d_test["nb_desc_input"]

        # control the number of desc are the same
        if self.nb_desc_input != self.nb_desc_test:
            print("===CHECK SIZE TRAIN AND TEST===")
            STOPHERE # crash script

    def GridOptimizeDNN(self, criteria_best_perf):

        # check if model is in the 
        l_filout = listdir(self.pr_out)
        for p_filout in l_filout:
            if search(".h5$", p_filout):
                self.p_model = self.pr_out + p_filout
                return 

        filout = open(self.pr_out + "optimization.txt", "w")
        best_perf = 0.0
        for activation in self.l_activation:
            for dense_layer in self.l_dense_layer:
                nb_layer = dense_layer
                self.initModel(activation, self.l_dense_candidate[-nb_layer])
                nb_layer = nb_layer - 1
                while nb_layer > 1:
                    self.addDense(activation, self.l_dense_candidate[-nb_layer])
                    nb_layer = nb_layer - 1
                for epochs in self.l_epochs:
                    for batch_size in self.l_batch_size:
                        d_pref_train = self.fit_compileModel(epochs, batch_size)
                        d_pref_test = self.evaluateOnTest()

                        filout.write("ACTIVATION: %s; Nb layers: %s; Epochs: %s; Batch size: %s;\n"%(activation, dense_layer, epochs, batch_size))
                        filout.write("TRAINNING SET\n%s\n"%("\n".join(["%s: %.4f"%(criteria, d_pref_train[criteria]) for criteria in d_pref_train.keys()])))
                        filout.write("TEST SET\n%s\n[====================]\n"%("\n".join(["%s: %.4f"%(criteria, d_pref_test[criteria]) for criteria in d_pref_test.keys()])))
                        
                        # save model
                        if  float(d_pref_train[criteria_best_perf]) > best_perf:
                            # remove previous model
                            l_items = listdir(self.pr_out)
                            for item in l_items:
                                if item.endswith(".h5"):
                                    remove(path.join(self.pr_out, item))
                            best_perf = float(d_pref_train[criteria_best_perf])
                            self.p_model = "%sbest_model_%s--%s--%s--%s.h5"%(self.pr_out, activation, dense_layer, epochs, batch_size)
                            self.model.save(self.p_model )

        filout.close()                

    def addDense(self, activation, dense_candidate):
        self.model.add(Dense(int(dense_candidate), activation=activation, kernel_initializer=self.kernel_initializer))
        #self.model.add(Dropout(0.5))

    def initModel(self, activation, dense_candidate):
        self.model = Sequential()
        if dense_candidate > self.nb_desc_input:
            dense_candidate = self.nb_desc_input
        self.model.add(Dense(dense_candidate, input_dim=self.nb_desc_input, activation=activation, kernel_initializer=self.kernel_initializer))
        #self.model.add(Dropout(0.5))
    
    def fit_compileModel(self, epochs, batch_size, l_metric = ["mse"], learning_rate = 0.1):

        # last layer
        if self.typeModel == "classification":
            self.model.add(Dense(1, activation='sigmoid'))
        else:
            self.model.add(Dense(1))

        if self.optimizer == "sgd":
            sgd = tf.keras.optimizers.SGD(learning_rate=learning_rate)
            self.optimizer = sgd
        #elif optimizer == "adam":
        #    adam =  tf.keras.optimizers.Adam(learning_rate=learning_rate)
        #    optimizer = adam
        #    print(adam)

        # compile the keras model
        self.model.compile(loss=self.loss, optimizer=self.optimizer, metrics=l_metric)
        # fit the keras model on the dataset
        self.model.fit(self.dataset_train, self.aff_train, epochs=epochs, batch_size=batch_size), 

        # evaluate the keras model
        self.model.evaluate(self.dataset_train, self.aff_train)
        y_pred = self.model.predict(self.dataset_train)
            
        return ML_toolbox.performance(self.aff_train, y_pred, self.typeModel)
        
    def evaluateOnTest(self):

        y_pred = self.model.predict(self.dataset_test)
        return ML_toolbox.performance(self.aff_train, y_pred, self.typeModel)
    
    def evaluateModel(self):

        if "p_model" in self.__dict__:
            model = load_model(self.p_model)
        else:
            print("No model ready to be loaded")
            return

        if not "d_perf" in self.__dict__:
            self.d_perf = {}

        # trainning set
        y_train_pred = model.predict(self.dataset_train)    
        self.d_perf["Train"] = ML_toolbox.performance(self.aff_train, y_train_pred, self.typeModel, p_filout=self.pr_out + "train_pred.csv")
            
        # test set
        y_pred = model.predict(self.dataset_test)
        self.d_perf["Test"] = ML_toolbox.performance(self.aff_test, y_pred, self.typeModel, p_filout=self.pr_out + "test_pred.csv")

    def combineResults(self):
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
            filout.write("%s-DNN\t%s\n"%(perf, "\t".join(["%.2f"%(self.d_perf[perf][criteria]) for criteria in l_criteria])))
        filout.close()     

    def CrossValidation(self):
        seed = 7
        d_CV = {}
        if self.typeModel == "classification":
            kfold = StratifiedKFold(n_splits=self.n_foldCV, shuffle=True, random_state=seed)
        else:
            kfold = KFold(n_splits=self.n_foldCV, shuffle=True, random_state=seed)
        # load parameters
        l_files = listdir(self.pr_out)
        for name_file in l_files:
            if search("h5$", name_file):
                l_parameter = name_file.split("_")[-1][:-3].split("--")
                act = l_parameter[0]
                nb_layer = int(l_parameter[1])
                epochs = int(l_parameter[2])
                batch_size = int(l_parameter[3])
                break
        

        d_train = deepcopy(self.dataset_train)
        Y_train = deepcopy(self.aff_train)

        for train, test in kfold.split(d_train, Y_train):
            self.dataset_train = d_train[train]
            self.aff_train = Y_train[train]

            self.dataset_test = d_train[test]
            self.aff_test = Y_train[test]

            # create model
            self.initModel(act, self.l_dense_candidate[-nb_layer])
            nb_layer_temp = nb_layer - 1
            while nb_layer_temp > 1:
                self.addDense(act, self.l_dense_candidate[-nb_layer_temp])
                nb_layer_temp = nb_layer_temp - 1

            d_pref_train = self.fit_compileModel(epochs, batch_size)
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

    def predictDNN(self, p_model, p_test, name_out = "test_pred"):
        # descriptor test set 
        with open(p_test) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        dataset_test = loadtxt(p_test, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        dataset_test = sc.fit_transform(dataset_test)
        
        y_test = loadtxt(p_test, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)
        d_test = toolbox.loadMatrix(p_test)

        model = load_model(p_model)
        y_pred = model.predict(dataset_test)

        # compute quality performance
        p_filout_sum = self.pr_out + name_out + ".sum"
        self.performance(y_test, y_pred, p_filout_sum)
        
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


    # create a main in the DNN
    def run_main(self):
        
        ## to shortcurt
        #p_combined_results = self.pr_out + "combined_perf.csv"
        #if path.exists(p_combined_results) and self.force_run == 0:
        #    return  

        self.loadSet()
        self.GridOptimizeDNN("AUC")
        self.evaluateModel()
        self.CrossValidation()
        self.combineResults()
        stophere

### not used in the DNN -> relocated in the QSAR class
    def mergeUndersampling(self, pr_rep):

        d_out = {}
        d_out ["train"] = {"Acc":[], "Se":[], "Sp":[], "MCC":[], "b-Acc":[]}
        d_out ["test"] = {"Acc":[], "Se":[], "Sp":[], "MCC":[], "b-Acc":[]}
        d_out ["CV"] = {"Acc":[], "Se":[], "Sp":[], "MCC":[], "b-Acc":[]}
        pr_out = pathFolder.createFolder(pr_rep + "mergeDNN/")
        for i in range(1, self.nb_repetition + 1):
            p_perf = "%s/%s/DNNclass/combined_perf.csv"%(pr_rep, i)
            d_perf = toolbox.loadMatrix(p_perf)
            d_out["train"]["Acc"].append(float(d_perf["Train-DNN"]["Acc"]))
            d_out["train"]["Se"].append(float(d_perf["Train-DNN"]["Se"]))
            d_out["train"]["Sp"].append(float(d_perf["Train-DNN"]["Sp"]))
            d_out["train"]["MCC"].append(float(d_perf["Train-DNN"]["MCC"]))
            d_out["train"]["b-Acc"].append(float(d_perf["Train-DNN"]["b-Acc"]))

            d_out["test"]["Acc"].append(float(d_perf["Test-DNN"]["Acc"]))
            d_out["test"]["Se"].append(float(d_perf["Test-DNN"]["Se"]))
            d_out["test"]["Sp"].append(float(d_perf["Test-DNN"]["Sp"]))
            d_out["test"]["MCC"].append(float(d_perf["Test-DNN"]["MCC"]))
            d_out["test"]["b-Acc"].append(float(d_perf["Test-DNN"]["b-Acc"]))

            d_out["CV"]["Acc"].append(float(d_perf["CV-DNN"]["Acc"]))
            d_out["CV"]["Se"].append(float(d_perf["CV-DNN"]["Se"]))
            d_out["CV"]["Sp"].append(float(d_perf["CV-DNN"]["Sp"]))
            d_out["CV"]["MCC"].append(float(d_perf["CV-DNN"]["MCC"]))
            d_out["CV"]["b-Acc"].append(float(d_perf["CV-DNN"]["b-Acc"]))
        
        l_perf = ["Acc", "Se", "Sp", "MCC", "b-Acc"]
        p_filout = pr_out + "average_perf.csv"
        filout = open(p_filout, "w")
        filout.write("\t" + "\t".join(["M-%s\tSD-%s"%(perf, perf)for perf in l_perf]) + "\n")

        for k in d_out.keys():
            filout.write("%s\t%s\n"%(k, "\t".join(["%s\t%s"%(np.average(d_out[k][perf]), np.std(d_out[k][perf]))for perf in l_perf])))
        filout.close()






#########
### FOR TESTING

    def buildDNNTest(self, l_metric = ["mse"]):
        

        # define the keras model
        model = Sequential()
        model.add(Dense(12, input_dim=self.nb_desc_input, activation='relu'))
        model.add(Dense(8, activation='relu'))
        model.add(Dense(1, activation='sigmoid'))
        # compile the keras model
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=l_metric)
        # fit the keras model on the dataset
        model.fit(self.dataset_train, self.aff_train, epochs=10, batch_size=10)

        # evaluate the keras model
        _, accuracy = model.evaluate(self.dataset_train, self.aff_train)
        print('ACC: %.2f' % (accuracy))

        y_pred = model.predict(self.dataset_test)
        y_pred = [1. if pred[0] > 0.5 else 0. for pred in y_pred]
        
        acc = metrics.accuracy_score(self.aff_test, y_pred)
        bacc = metrics.balanced_accuracy_score(self.aff_test, y_pred)
        mcc = metrics.matthews_corrcoef(self.aff_test, y_pred)
        recall = metrics.recall_score(self.aff_test, y_pred)

        print("======= TEST ======")
        print("Acc: ", acc)
        print("b-Acc: ", bacc)
        print("MCC: ", mcc)
        print("Recall: ", recall)

    def test(self):

        #dataset = loadtxt('C:/Users/aborr/research/ILS/HERG/data/pima-indians-diabetes.data.csv', delimiter=',')
        dataset = loadtxt('/mnt/c/Users/aborr/research/ILS/HERG/data/pima-indians-diabetes.data.csv', delimiter=',')
        X = dataset[:,0:8]
        y = dataset[:,8]

        print(X)

        # define the keras model
        model = Sequential()
        model.add(Dense(12, input_dim=8, activation='relu'))
        model.add(Dense(8, activation='relu'))
        model.add(Dense(1, activation='sigmoid'))
        # compile the keras model
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
        # fit the keras model on the dataset
        model.fit(X, y, epochs=150, batch_size=10)
        # evaluate the keras model
        _, accuracy = model.evaluate(X, y)
        print('Accuracy: %.2f' % (accuracy*100))



        return 

    def test2(self):

        (x_train, y_train), (x_test, y_test) = mnist.load_data() 

        print(x_train)
        sss

        x_train = x_train.reshape(60000, 784) 
        x_test = x_test.reshape(10000, 784) 
        x_train = x_train.astype('float32') 
        x_test = x_test.astype('float32') 
        x_train /= 255 
        x_test /= 255 

        y_train = keras.utils.to_categorical(y_train, 10) 
        y_test = keras.utils.to_categorical(y_test, 10) 

        model = Sequential() 
        model.add(Dense(512, activation='relu', input_shape = (784,))) 
        model.add(Dropout(0.2)) 
        model.add(Dense(512, activation = 'relu')) 
        model.add(Dropout(0.2)) 
        model.add(Dense(10, activation = 'softmax'))
        model.compile(loss = 'categorical_crossentropy', 
        optimizer = RMSprop(), 
        metrics = ['accuracy']) 

        history = model.fit(x_train, y_train, 
        batch_size = 128, epochs = 20, verbose = 1, validation_data = (x_test, y_test))
