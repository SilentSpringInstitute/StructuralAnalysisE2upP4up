from os import path
from statistics import mean, stdev
from shutil import copyfile
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys

import pathFolder
import runExternal
import toolbox



class MakePlots_fromDesc:
    def __init__(self, p_dataset, pr_out, pr_desc,  p_desc="", p_hormone_similarity = "", p_FP="", p_opera_all="", cor_val="", max_quantile=""):
        self.p_desc = p_desc
        self.p_dataset = p_dataset
        self.p_hormone_similarity = p_hormone_similarity
        self.p_opera = p_opera_all
        self.p_FP = p_FP
        self.pr_desc = pr_desc
        self.pr_out = pr_out
        self.cor_val = cor_val
        self.max_quantile = max_quantile

    def hclusterByProp(self):
        """clustering of Hclust and compute enrichement by cluster
        - agreg: ward2
        - clustering: hclust
        - opt number of cluster gap-stat
        """

        pr_out = pathFolder.createFolder(self.pr_out + "hclustDendo/")

        # run dendogram with circle of prop
        p_cluster = pr_out + "cluster_hclust_ward2_gapstat.csv"
        p_enrich = pr_out + "prob_by_clusters.csv"
        if not path.exists(p_enrich) and not path.exists(p_cluster):
            runExternal.dendogramClusterProp(self.p_dataset, self.p_desc, self.p_hormone_similarity, pr_out, self.cor_val, self.max_quantile)
            if path.exists(p_cluster):
                runExternal.enrichmentByCluster(self.p_dataset, p_cluster, pr_out)

        # extract png by cluster
        pr_clusters_png = pathFolder.createFolder(pr_out + "png_by_cluster/")
        d_cluster = toolbox.loadMatrix(p_cluster)

        for casrn in d_cluster.keys():
            cluster = d_cluster[casrn]["cluster"]
            pr_cluster = pathFolder.createFolder(pr_clusters_png + cluster + "/")
            copyfile("%sPNG/%s.png"%(self.pr_desc, casrn), "%s%s.png"%(pr_cluster, casrn))



    def hclusterFromFPByProp(self):

        pr_out = pathFolder.createFolder(self.pr_out + "hclustDendo/")

        # run dendogram with circle of prop
        p_dendo = pr_out + "dendo_cluster_name.png"
        if not path.exists(p_dendo):
            runExternal.dendogramFPProp(self.p_dataset, self.p_FP, pr_out)


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

    def SOMClustering(self, nb_cluster):

        # need to prepare the dataset before the clustering
        if not "p_desc_cleaned" in self.__dict__:
            self.prepDesc()

        pr_out = pathFolder.createFolder(self.pr_out + "SOM/")
        p_model = pr_out + "SOM_model.RData"
        if not path.exists(p_model) or nb_cluster == 0:
            runExternal.SOM(self.p_desc_cleaned, "0", pr_out, nb_cluster)

        p_cluster = pr_out + "SOM_Clusters"
        if path.exists(p_cluster) and nb_cluster != 0:
            runExternal.enrichmentByCluster(self.p_dataset, p_cluster, pr_out)

            # extract png by cluster
            pr_clusters_png = pathFolder.createFolder(pr_out + "png_by_cluster/")
            d_cluster = toolbox.loadMatrix(p_cluster)

            for casrn in d_cluster.keys():
                cluster = d_cluster[casrn]["cluster"]
                pr_cluster = pathFolder.createFolder(pr_clusters_png + cluster + "/")
                copyfile("%sPNG/%s.png"%(self.pr_desc, casrn), "%s%s.png"%(pr_cluster, casrn))

        self.p_model_SOM = p_model

    def SOMMapProp(self):

        if not "p_model_SOM" in self.__dict__:
            print("Generate SOM first")
            return 

        pr_out = pathFolder.createFolder(self.pr_out + "SOM/")
        runExternal.projectPropInSOM(self.p_model_SOM, self.p_dataset, pr_out)

    def SOMHormoneSimilarity(self):
        if not "p_model_SOM" in self.__dict__:
            print("Generate SOM first")
            return 

        pr_out = pathFolder.createFolder(self.pr_out + "SOM_hormones/")
        runExternal.SOM_hormone(self.p_model_SOM, self.p_hormone_similarity, pr_out)


    def signifDescBySOMCluster(self):

        # check if SOM is computed
        p_cluster = self.pr_out + "SOM/SOM_Clusters"
        if not path.exists(p_cluster):
            print("ERROR: no cluster file existed")
            return 
        
        pr_out = pathFolder.createFolder(self.pr_out + "SOM/DescriptorSignif/")
        runExternal.descSignifByCluster(self.p_desc_cleaned, p_cluster, pr_out)

    def extract_actBySOMCluster(self, pr_png):

        p_cluster = pathFolder.createFolder(self.pr_out + "SOM/") + "SOM_Clusters_act"
        if not path.exists(p_cluster):
            return "Error"

        pr_out = pathFolder.createFolder(self.pr_out + "SOM/Active_by_cluster/")

        dcluster = toolbox.loadMatrix(p_cluster, sep = ",")
        for CASRN in dcluster.keys():
            try:cluster = dcluster[CASRN]["Cluster"]
            except:cluster = dcluster[CASRN]["x"]
            pr_cluster = pathFolder.createFolder(self.pr_out + "SOM/Active_by_cluster/" + cluster + "/")
            pathFolder.createFolder(pr_cluster)

            try:copyfile(pr_png + CASRN + ".png", pr_cluster + CASRN + ".png")
            except: pass