from os import path
from statistics import mean, stdev
from shutil import copyfile
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys

import pathFolder
import runExternal
import toolbox



class analysis:
    def __init__(self, p_dataset, p_desc, p_opera, pr_out, cor_val, max_quantile):
        self.p_desc = p_desc
        self.p_dataset = p_dataset
        self.p_opera = p_opera
        self.pr_out = pr_out
        self.cor_val = cor_val
        self.max_quantile = max_quantile



    def prepDesc(self):
        pr_out = pathFolder.createFolder(self.pr_out + "Cleaned_Data/")

        # add shortcut
        p_desc_cleaned = pr_out + "desc1D2D_cleaned.csv"
        if not path.exists(p_desc_cleaned):
            runExternal.preprocData(self.p_desc, pr_out, self.cor_val, self.max_quantile)
        self.p_desc_cleaned = p_desc_cleaned


    def PCA_plot(self):

        pr_out = pathFolder.createFolder(self.pr_out + "PCA/")
        runExternal.PCA(self.p_desc_cleaned, pr_out)



    def HClust_plot(self, p_opera):

        pr_out = pathFolder.createFolder(self.pr_out + "HClustCircular/")
        runExternal.HClust(self.p_desc_cleaned, self.p_dataset, p_opera, pr_out)

    def clustering(self):

        pr_out = pathFolder.createFolder(self.pr_out + "Clustering/")
        runExternal.Clust(self.p_desc_cleaned, self.p_opera, pr_out)


    def histDesc(self):
        
        pr_out = pathFolder.createFolder(self.pr_out + "histDesc/")
        runExternal.drawHist(self.p_desc, self.p_desc_cleaned, pr_out)
    

    def FPTanimoto(self, l_typeFP):

        # load SMILES in list
        d_desc = toolbox.loadMatrix(self.p_desc, sep = "\t") # SMILES are already cleaned

        pr_out = pathFolder.createFolder(self.pr_out + "FP/")

        d_chem = {}
        for chem in d_desc.keys():
            d_chem[chem] = {}
            d_chem[chem]["SMILES"] = d_desc[chem]["SMILES"]
            d_chem[chem]["chem"] = Chem.MolFromSmiles(d_chem[chem]["SMILES"])

        l_chem = list(d_chem.keys())
        dFP = {}

        i = 0
        imax = len(l_chem)
        while i < imax:
            CASRN1 = l_chem[i]
            dFP[CASRN1] = {}

            # FP for CASRN1
            if "topo" in l_typeFP:
                fps1_topo = Chem.RDKFingerprint(d_chem[CASRN1]["chem"])
            if "MACCS" in l_typeFP:
                fps1_MACCS = MACCSkeys.GenMACCSKeys(d_chem[CASRN1]["chem"])
            if "Morgan" in l_typeFP:
                fps1_Morgan = Chem.AllChem.GetMorganFingerprint(d_chem[CASRN1]["chem"], 2)


            j = i
            while j < imax:
                CASRN2 = l_chem[j]
                dFP[CASRN1][CASRN2] = {}

                # FP for CASRN1
                if "topo" in l_typeFP:
                    fps2_topo = Chem.RDKFingerprint(d_chem[CASRN2]["chem"])
                    sim = DataStructs.FingerprintSimilarity(fps1_topo, fps2_topo)
                    dFP[CASRN1][CASRN2]["topo"] = sim
                if "MACCS" in l_typeFP:
                    fps2_MACCS = MACCSkeys.GenMACCSKeys(d_chem[CASRN2]["chem"])
                    sim = DataStructs.FingerprintSimilarity(fps1_MACCS, fps2_MACCS)
                    dFP[CASRN1][CASRN2]["MACCS"] = sim
                if "Morgan" in l_typeFP:
                    fps2_Morgan = Chem.AllChem.GetMorganFingerprint(d_chem[CASRN2]["chem"], 2)
                    sim = DataStructs.DiceSimilarity(fps1_Morgan,fps2_Morgan)
                    dFP[CASRN1][CASRN2]["Morgan"] = sim


                j = j + 1
            i = i + 1

        # write matix of sim
        for typeFP in l_typeFP:
            p_filout = pr_out + "sim_" + typeFP
            filout = open(p_filout, "w")
            filout.write("\t" + "\t".join(l_chem) + "\n")
            i = 0
            imax = len(l_chem)
            while i < imax:
                lw = [l_chem[i]]
                j = 0
                while j < imax:
                    try:
                        lw.append(dFP[l_chem[i]][l_chem[j]][typeFP])
                    except:
                        lw.append(dFP[l_chem[j]][l_chem[i]][typeFP])

                    j = j + 1
                filout.write("\t".join([str(w) for w in lw]) + "\n")
                i = i + 1
            
            filout.close()
            runExternal.cardSimMatrix(p_filout)
