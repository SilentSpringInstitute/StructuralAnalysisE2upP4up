from sklearn import metrics
from numpy import loadtxt, arange
from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

def performance(y_real, y_pred, typeModel, th_prob = 0.5, p_filout = ""):
            
    if typeModel == "classification":

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



def loadSet(p_set, variableToPredict = 0):
    """
    Data in input 
        - sep = ,
        - aff in last col
        - rownames in the first col
    """

    d_out = {}
    # descriptor train set 
    with open(p_set) as f:
        #determining number of columns from the first line of text
        n_cols = len(f.readline().split(","))
    
    if variableToPredict == 1:
        d_out["dataset"] = loadtxt(p_set, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        d_out["dataset"] = sc.fit_transform(d_out["dataset"])
        # take last col
        d_out["aff"] = loadtxt(p_set, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)
        d_out["nb_desc_input"] = n_cols - 2 # remove rowname and the aff
    else:
        d_out["dataset"] = loadtxt(p_set, delimiter=",",usecols=arange(1, n_cols), skiprows=1)
        d_out["dataset"] = sc.fit_transform(d_out["dataset"])
        d_out["nb_desc_input"] = n_cols - 1 # remove only the rownames

    d_out["id"] = loadtxt(p_set, dtype=str,  delimiter=",",usecols=arange(0, 1), skiprows=1)

    return d_out