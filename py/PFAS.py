import dataset
import pathFolder
import analysis
import toolbox

from os import path

# selected physico chemical descriptor from OPERa
L_OPERA_DESC = ['LogP_pred', 'MP_pred', 'BP_pred', 'LogVP_pred', 'LogWS_pred', 'LogHL_pred', 'RT_pred', 'LogKOA_pred', 'ionization', 'LogD55_pred', 'LogD74_pred', 'LogOH_pred', 'LogBCF_pred', 'BioDeg_LogHalfLife_pred', 'ReadyBiodeg_pred', 'LogKM_pred', 'LogKoc_pred', 'FUB_pred', 'Clint_pred']


class PFAS:
    def __init__(self, p_dataset, pr_out):
        self.p_dataset = p_dataset
        self.pr_out = pr_out
        self.loadDB = 1

    def computeDesc(self):
        # load dataset
        c_dataset = dataset.dataset(self.p_dataset, self.pr_out)
        c_dataset.loadDataset(loadDb=self.loadDB)


        # compute desc 2D
        p_desc = c_dataset.computeStructuralDesc()

        # compute PNG
        c_dataset.computePNG()

        # compute OPERA
        p_desc_opera = c_dataset.computeOPERADesc()

        # filout
        self.p_desc_rdkit = p_desc
        self.p_desc_OPERA = p_desc_opera
        self.pr_desc = c_dataset.pr_desc
        self.c_dataset = c_dataset



    def buildDescSet(self, l_type_desc):
        """Select from OPERA only physico chem descriptors"""

        pr_results = pathFolder.createFolder(self.pr_out + "-".join(l_type_desc) + "/")
        self.pr_results = pr_results


        p_filout = pr_results + "-".join(l_type_desc) + ".csv"
        if path.exists(p_filout):
            self.p_desc = p_filout

        if len(l_type_desc) == 1 and l_type_desc[0] == "rdkit":
            self.p_desc = self.p_desc_rdkit


        # open RDKIT use for CARSN and SMILES
        d_rdkit = toolbox.loadMatrix(self.p_desc_rdkit, sep="\t")

        d_out = {}
        l_header = ["CASRN", "SMILES"]
        if "opera" in l_type_desc:
            d_OPERA = toolbox.loadMatrix(self.p_desc_OPERA, sep = ",")
            l_header = l_header + L_OPERA_DESC

            for chem in d_OPERA:
                if not chem in list(d_out.keys()):
                    d_out[chem] = {}
                    d_out[chem]["CASRN"] = d_rdkit[chem]["CASRN"]
                    d_out[chem]["SMILES"] = d_rdkit[chem]["SMILES"]
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
            filout.write("\t".join([d_out[chem][h] for h in l_header]) + "\n")
        filout.close()

        self.p_desc = p_filout



    def getChemTier(self, l_tiers):

        pr_results = pathFolder.createFolder(self.pr_results + "Tier" + "-".join(l_tiers) + "/")
        self.pr_results = pr_results

        # need to define a new p_desc
        l_CASRN = []
        for CASRN in self.c_dataset.d_dataset.keys():
            if self.c_dataset.d_dataset[CASRN]["Tier"] in l_tiers:
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
        self.c_dataset.combineDataset(p_dataset2)




    def computeBiotransformation(self, p_list_mater=""):

        # compute biotransformation 
        self.c_dataset.predictBiotransformation()

        # extract biostransformation
        self.c_dataset.extractBioTranformationChemical(p_list_mater)
        self.c_dataset.searchBiotransformedProduceInDB()

        self.c_dataset.computeDescProductBiotransformed()




    def analysisDescBasedData(self, cor_val, max_quantile, PCA=0, Hclust=0, clustering=0, SOM=0, histDesc =0, FP = 0):

        if not "p_desc" in self.__dict__:
            self.err = 1
            self.log = "load dataset first"

        cAnalysis = analysis.analysis(self.p_dataset, self.p_desc, self.p_desc_OPERA, self.pr_results, cor_val, max_quantile)
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
            size = 15
            cAnalysis.generate_SOM(15)
            cAnalysis.signifDescBySOMCluster()
            cAnalysis.extract_actBySOMCluster(pr_desc + "PNG/") # have to run !!!!

        # 2.5 histogram by descriptor
        if histDesc == 1:
            cAnalysis.histDesc()  

        # 2.6 Tanimoto in fingerprint
        if FP == 1:
            cAnalysis.FPTanimoto(["topo", "MACCS", "Morgan"])
