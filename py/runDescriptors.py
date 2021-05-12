import dataset
import pathFolder
import toolbox
import runExternal

from os import path
from shutil import copyfile

# selected physico chemical descriptor from OPERa
L_OPERA_DESC = ['LogP_pred', 'MP_pred', 'BP_pred', 'LogVP_pred', 'LogWS_pred', 'LogHL_pred', 'RT_pred', 'LogKOA_pred', 'ionization', 'LogD55_pred', 'LogD74_pred', 'LogOH_pred', 'LogBCF_pred', 'BioDeg_LogHalfLife_pred', 'ReadyBiodeg_pred', 'LogKM_pred', 'LogKoc_pred', 'FUB_pred', 'Clint_pred']


class runDescriptors:
    def __init__(self, d_dataset, pr_desc, pr_out):
        self.d_dataset = d_dataset
        self.pr_out = pr_out
        self.pr_desc = pr_desc
        self.loadDB = 1

    def computeDesc(self, p_dataset, pr_out):
        # load dataset
        c_dataset = dataset.dataset(p_dataset, pr_out)
        c_dataset.loadDataset(loadDb=self.loadDB)
        self.c_dataset = c_dataset

        # compute desc 2D
        p_desc_rdkit = self.c_dataset.computeStructuralDesc(self.pr_desc)

        # compute PNG
        self.c_dataset.computePNG(self.pr_desc)

        # compute OPERA
        p_desc_opera = self.c_dataset.computeOPERADesc(self.pr_desc)

        # filout
        d_out = {}
        d_out["rdkit"] = p_desc_rdkit
        d_out["OPERA"] = p_desc_opera

        return d_out

    def computeAllOperaPred(self, p_dataset, pr_out):
        # load dataset
        if not "c_dataset" in self.__dict__:
            c_dataset = dataset.dataset(p_dataset, pr_out)
            c_dataset.loadDataset(loadDb=self.loadDB)
            self.c_dataset = c_dataset

        p_desc_opera = self.c_dataset.computeAllOPERADesc(self.pr_desc)
        return p_desc_opera

    def buildDescSet(self, dataset, l_type_desc, pr_results):
        """Select from OPERA only physico chem descriptors"""

        if not "d_desc" in self.__dict__:
            self.compute_all()

        p_filout = pr_results + "-".join(l_type_desc) + ".csv"
        if path.exists(p_filout):
            self.p_desc = p_filout

        # need to reformat the 
        #if len(l_type_desc) == 1 and l_type_desc[0] == "rdkit":
        #    self.p_desc = self.d_dataset[dataset]["rdkit"]


        # open RDKIT use for CARSN and SMILES
        d_rdkit = toolbox.loadMatrix(self.d_desc[dataset]["rdkit"], sep="\t")
        
        d_out = {}
        l_header = ["CASRN", "SMILES"]
        if "opera" in l_type_desc:
            d_OPERA = toolbox.loadMatrix(self.d_desc[dataset]["OPERA"], sep = ",")
            l_header = l_header + L_OPERA_DESC

            for chem in d_OPERA.keys():
                if not chem in list(d_out.keys()) and chem in list(d_rdkit.keys()):
                    d_out[chem] = {}
                    d_out[chem]["CASRN"] = d_rdkit[chem]["CASRN"]
                    d_out[chem]["SMILES"] = d_rdkit[chem]["SMILES"]
                else:
                    continue
                for h in L_OPERA_DESC:
                    d_out[chem][h] = d_OPERA[chem][h]



        if "rdkit" in l_type_desc:
            l_h_drdkit = list(d_rdkit[list(d_rdkit.keys())[0]].keys())
            l_h_drdkit.remove("CASRN")
            l_h_drdkit.remove("SMILES")
            l_header = l_header + l_h_drdkit

            for chem in d_rdkit.keys():
                if not chem in d_out.keys():
                    d_out[chem] = {}
                    d_out[chem]["CASRN"] = d_rdkit[chem]["CASRN"]
                    d_out[chem]["SMILES"] = d_rdkit[chem]["SMILES"]
                
                for h in l_h_drdkit:
                    d_out[chem][h] = d_rdkit[chem][h]

        filout = open(p_filout, "w")
        filout.write("\t".join(l_header) + "\n")
        for chem in d_out.keys():
            filout.write("\t".join(["NA" if not h in list(d_out[chem].keys()) else d_out[chem][h] for h in l_header]) + "\n")
        filout.close()

        return p_filout

    def getChemByColValues(self, k, l_vals, pr_out):
        """
        Refine the table of chemicals using values in the on column
        arg
        k: name column to read
        l_vals: values of k that need to be selected

        return
        self.p_desc: refined descriptor table.
        """

        pr_results = pathFolder.createFolder(pr_out + str(k) + "-".join(l_vals) + "/")

        # need to define a new p_desc
        l_CASRN = []
        for CASRN in self.c_dataset.d_dataset.keys():
            if self.c_dataset.d_dataset[CASRN][k] in l_vals:
                l_CASRN.append(CASRN)

        # redefine p_desc
        d_desc = toolbox.loadMatrix(self.p_desc)
        l_h = list(d_desc[list(d_desc.keys())[0]].keys())
        l_h.remove("CASRN")
        l_h.remove("SMILES")

        p_desc_out = pr_results + "desc.csv"
        f_desc_out = open(p_desc_out, "w")
        f_desc_out.write("CASRN\tSMILES\t" + "\t".join(l_h) + "\n")
        for chem in d_desc.keys():
            if chem in l_CASRN:
                f_desc_out.write("%s\t%s\t%s\n"%(chem, d_desc[chem]["SMILES"], "\t".join([d_desc[chem][h] for h in l_h])))
        f_desc_out.close()

        self.p_desc = p_desc_out

    def combineDataset(self, p_dataset2):
        c_dataset = dataset.dataset(self.p_dataset, self.pr_out)
        self.c_dataset = c_dataset
        
        self.c_dataset.combineDataset(p_dataset2)

    def computeBiotransformation(self, p_list_mater=""):

        # compute biotransformation 
        self.c_dataset.predictBiotransformation()

        # extract biostransformation
        self.c_dataset.extractBioTranformationChemical(p_list_mater)
        self.c_dataset.searchBiotransformedProduceInDB()

        self.c_dataset.computeDescProductBiotransformed()

    def analysisDescBasedData(self, cor_val, max_quantile, PCA=0, Hclust=0, clustering=0, SOM=0, SOM_SIZE =10, histDesc =0, FP = 0):
        """
        need to be rewrtie
        """


        if not "p_desc_rdkit" in self.__dict__:
            self.err = 1
            self.log = "load dataset first"
            return 

        cAnalysis = analysis.analysis(self.p_dataset, self.p_desc_rdkit, self.p_desc_OPERA, self.pr_results, cor_val, max_quantile)
        cAnalysis.prepDesc()

        # 2.1 PCA
        if PCA == 1:
            cAnalysis.PCA_plot()

        # 2.2 Hclust
        if Hclust == 1:
            cAnalysis.HClust_plot(self.p_desc_OPERA)

        # 2.3 Clustering
        if clustering == 1: 
            cAnalysis.clustering()

        # 2.4 SOM
        if SOM == 1:
            cAnalysis.generate_SOM(SOM_SIZE)
            cAnalysis.signifDescBySOMCluster()
            cAnalysis.extract_actBySOMCluster(self.pr_desc + "PNG/") # have to run !!!!

        # 2.5 histogram by descriptor
        if histDesc == 1:
            cAnalysis.histDesc()  

        # 2.6 Tanimoto in fingerprint
        if FP == 1:
            cAnalysis.FPTanimoto(["topo", "MACCS", "Morgan"])

    def analysisChemList(self, p_chemlist):

        d_chemList = toolbox.loadMatrix(p_chemlist, ",") 
        pr_out = pathFolder.createFolder(self.pr_out + "chem_list/")

        l_chemlist = list(d_chemList[list(d_chemList.keys())[0]].keys())
        l_chemlist.remove('INPUT')
        l_chemlist.remove('DTXSID')
        l_chemlist.remove('PREFERRED_NAME')
        l_chemlist.remove('FOUND_BY')

        d_out = {}
        for chemlist in l_chemlist:
            d_out[chemlist] = 0

        for chem in d_chemList.keys():
            for chemlist in l_chemlist:
                if d_chemList[chem][chemlist] == "Y":
                    d_out[chemlist] = d_out[chemlist] + 1


        p_filout = pr_out + "count_chemlist"
        filout = open(p_filout, "w")
        filout.write("chem_list\tcount\n")
        for chemlist in d_out.keys():
            if d_out[chemlist] != 0:
                filout.write("%s\t%s\n"%(chemlist, d_out[chemlist]))
        filout.close()

        runExternal.barplotchemlist(p_filout)

        return 
    
    def compute_all(self):
        """
        Compute all descriptors available for the this project
        - can add option to compute only partially the descriptor set
        return: self.d_desc
        """
        d_out = {}

        for dataset in self.d_dataset.keys():
            d_out[dataset] = {}
            pr_chem = pathFolder.createFolder(self.pr_out + dataset + "/")
            d_out[dataset].update(self.computeDesc(self.d_dataset[dataset], pr_chem))
            d_out[dataset]["all OPERA pred"] = self.computeAllOperaPred(self.d_dataset[dataset], pr_chem)
        self.d_desc = d_out

    def png_by_list(self, pr_png):

        for dataset in self.d_dataset.keys():
            pr_png_list = pathFolder.createFolder(pr_png + dataset + "/")
            d_dataset = toolbox.loadMatrix(self.d_dataset[dataset])
            for chem in d_dataset.keys():
                p_png_desc = self.pr_desc + "PNG/" + chem + ".png"
                if path.exists(p_png_desc):
                    copyfile(p_png_desc, pr_png_list + chem + ".png")

                