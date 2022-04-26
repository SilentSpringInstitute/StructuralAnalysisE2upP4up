import pandas
from re import search


def loadExcelSheet(p_excel, name_sheet, k_head):
    """
    TO DO: Add check duplicate rownames
    """

    d_out = {}
    # load MC list
    data_frame = pandas.read_excel(p_excel, name_sheet, engine='openpyxl')
    data_size = data_frame.shape
    nb_row = data_size[0]
    nb_col = data_size[1]
    l_col = data_frame.columns


    i = 0
    while i < nb_row:
        rowname = str(data_frame.iloc[i][k_head])
        rowname = rowname.replace("'", "")
        if not rowname in list(d_out.keys()):
            d_out[rowname] = {}
        j = 0
        while j < nb_col:
            colname = l_col[j]
            if search("Unnamed:", colname):
                j = j + 1
                continue
            else:
                d_out[rowname][colname] = data_frame.iloc[i][l_col[j]]
            j = j + 1
        i = i + 1
    
    return d_out

def loadMatrixToList(pmatrixIn, sep = "\t"):

    filin = open(pmatrixIn, "r", encoding="utf8", errors='ignore')
    llinesMat = filin.readlines()
    filin.close()

    l_out = []
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1) -1):
        lheaders = ["ID"] + lheaders


    i = 0
    while i < len(lheaders):
        if lheaders[i] == "":
            lheaders[i] = "ID"
        i += 1

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        j = 0
        if len(lvalues) != len(lheaders):
            print("ERROR - line: ", i)
            print(lvalues)
            print(lheaders)
        jmax = len(lheaders)
        dtemp = {}
        while j < jmax:
            try:dtemp[lheaders[j]] = lvalues[j]
            except:pass
            j += 1
        l_out.append(dtemp)
        i += 1

    return l_out

def loadMatrix(pmatrixIn, sep = "\t"):

    filin = open(pmatrixIn, "r", encoding="utf-8-sig", errors='ignore')
    llinesMat = filin.readlines()
    filin.close()

    dout = {}
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    if len(line1) == 1:
        sep = ","
        lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1) -1):
        lheaders = ["ID"] + lheaders


    i = 0
    while i < len(lheaders):
        if lheaders[i] == "":
            lheaders[i] = "ID"
        i += 1

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        kin = lvalues[0]
        dout[kin] = {}
        j = 0
        if len(lvalues) != len(lheaders):
            print("ERROR - line: ", i)
            print(lvalues)
            print(lheaders)
        jmax = len(lheaders)
        while j < jmax:
            try:dout[kin][lheaders[j]] = lvalues[j]
            except:pass
            j += 1
        i += 1

    return dout

def formatLine(linein):

    linein = linein.replace("\n", "")
    linenew = ""

    imax = len(linein)
    i = 0
    flagchar = 0
    while i < imax:
        if linein[i] == '"' and flagchar == 0:
            flagchar = 1
        elif linein[i] == '"' and flagchar == 1:
            flagchar = 0

        if flagchar == 1 and linein[i] == ",":
            linenew = linenew + " "
        else:
            linenew = linenew + linein[i]

        i += 1

    linenew = linenew.replace('\"', "")
    return linenew

def writeMatrix(ddesc, pdescAct, sep = "\t"):


    filout = open(pdescAct, "w")
    lheader = list(ddesc[list(ddesc.keys())[0]].keys())

    # put header in first

    if "CAS" in lheader:
        del lheader[lheader.index("CAS")]
        lheader = ["CAS"] + lheader
    elif "CASID" in lheader:
        del lheader[lheader.index("CASID")]
        lheader = ["CASID"] + lheader
    else:
        lheader = ["CASID"] + lheader
        for casID in list(ddesc.keys()):
            ddesc[casID]["CASID"] = casID


    filout.write(sep.join(lheader) + "\n")
    for casID in list(ddesc.keys()):
        lval = [str(ddesc[casID][i]) for i in lheader]
        filout.write(sep.join(lval) + "\n")
    filout.close()

def formatLineDataset(linein):

    linein = linein.replace("\n", "")
    linenew = ""

    imax = len(linein)
    i = 0
    flagchar = 0
    while i < imax:
        if linein[i] == '"' and flagchar == 0:
            flagchar = 1
        elif linein[i] == '"' and flagchar == 1:
            flagchar = 0

        if flagchar == 1 and linein[i] == ",":
            linenew = linenew + " "
        else:
            linenew = linenew + linein[i]
        i += 1

    linenew = linenew.replace('\"', "")
    return linenew

def colNameDict(dict_in):

    k1 = list(dict_in.keys())[0]
    return list(dict_in[k1].keys())

def combineDict(d_main, d_toadd):

    # define list of k that need to be added in the main dictionnary
    l_kadd = []
    for k in d_main.keys():
        if k in list(d_toadd.keys()):
            for k_toadd in d_toadd[k].keys():
                if not k_toadd in list(d_main[k].keys()) and k_toadd not in l_kadd:
                    l_kadd.append(k_toadd)

    for kadd in l_kadd:
        for k in d_main.keys():
            try:
                d_main[k][kadd] = d_toadd[k][kadd]
            except:
                d_main[k][kadd] = "NA"