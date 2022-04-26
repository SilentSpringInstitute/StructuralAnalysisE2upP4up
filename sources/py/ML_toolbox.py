from sklearn import metrics
from numpy import loadtxt, arange
import pandas
from sklearn.preprocessing import StandardScaler
from tensorflow.python.module.module import valid_identifier
sc = StandardScaler()

def performance(y_real, y_pred, typeModel, th_prob = 0.5, p_filout = ""):
            
    if typeModel == "classification":

        # change prob -> value
        # case of 2 values predicted
        try:y_pred = [1. if pred[1] > th_prob else 0. for pred in y_pred]
        except:y_pred = [1. if pred[0] > th_prob else 0. for pred in y_pred]

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

def loadSet(p_set, l_col_to_order = [], variableToPredict = "Aff", sep = ","):
    """
    Data in input open first with panda and next convert to numpy format to exit
        - sep = ,
        - aff in last col with the name in input
        - also col selected previously
    """

    d_out = {}

    # open data
    d_out["dataset"] = pandas.read_csv(p_set, sep = sep)
    print(d_out["dataset"])
    
    # take affinity col
    if variableToPredict != "":
        l_aff = d_out["dataset"][variableToPredict]
        l_aff = l_aff.to_numpy()
        d_out["aff"] = l_aff
        d_out["dataset"] = d_out["dataset"].drop(columns = [variableToPredict])
    
    # extra ID 
    if "ID" in list(d_out["dataset"].keys()):
        l_id =  list(d_out["dataset"]["ID"])
    elif "CASRN" in list(d_out["dataset"].keys()):
        l_id =  list(d_out["dataset"]["CASRN"])
    else:
        l_id = []
    d_out["id"] = l_id
    d_out["dataset"] = d_out["dataset"].iloc[:, 1:]
    
    # select specific col
    if l_col_to_order != []:
        d_out["dataset"] = d_out["dataset"][l_col_to_order]

    # list of features
    l_features = list(d_out["dataset"].columns)
    d_out["features"] = l_features
    d_out["nb_desc_input"] = len(l_features)

    # format dataset here
    d_out["dataset"] = d_out["dataset"].to_numpy()

    return d_out