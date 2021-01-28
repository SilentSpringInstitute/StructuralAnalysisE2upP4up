import dataset
import pathFolder
import analysis
import toolbox
import runExternal

from os import path

# selected physico chemical descriptor from OPERa
L_OPERA_DESC = ['LogP_pred', 'MP_pred', 'BP_pred', 'LogVP_pred', 'LogWS_pred', 'LogHL_pred', 'RT_pred', 'LogKOA_pred', 'ionization', 'LogD55_pred', 'LogD74_pred', 'LogOH_pred', 'LogBCF_pred', 'BioDeg_LogHalfLife_pred', 'ReadyBiodeg_pred', 'LogKM_pred', 'LogKoc_pred', 'FUB_pred', 'Clint_pred']


class Chemicals:
    def __init__(self, p_dataset, pr_out):
        self.p_dataset = p_dataset
        self.pr_out = pr_out
        self.loadDB = 0

    def computeDesc(self, pr_desc):
        # load dataset
        if not "c_dataset" in self.__dict__:
            c_dataset = dataset.dataset(self.p_dataset, self.pr_out)
            c_dataset.loadDataset(loadDb=self.loadDB)
            self.c_dataset = c_dataset


        # compute desc 2D
        p_desc = self.c_dataset.computeStructuralDesc(pr_desc)

        # compute PNG
        self.c_dataset.computePNG(pr_desc)

        # compute OPERA
        p_desc_opera = self.c_dataset.computeOPERADesc(pr_desc)

        # filout
        self.p_desc_rdkit = p_desc
        self.p_desc_OPERA = p_desc_opera
        self.pr_desc = self.c_dataset.pr_desc
       
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
            cAnalysis.generate_SOM(10)
            cAnalysis.signifDescBySOMCluster()
            cAnalysis.extract_actBySOMCluster(self.pr_desc + "PNG/") # have to run !!!!

        # 2.5 histogram by descriptor
        if histDesc == 1:
            cAnalysis.histDesc()  

        # 2.6 Tanimoto in fingerprint
        if FP == 1:
            cAnalysis.FPTanimoto(["topo", "MACCS", "Morgan"])

    def analysisToxPrint(self, p_ToxPrint):
        
        d_toxprint = toolbox.loadMatrix(p_ToxPrint, ",") 
        pr_out = pathFolder.createFolder(self.pr_out + "ToxPrint/")

        l_toxprints = list(d_toxprint[list(d_toxprint.keys())[0]].keys())
        l_toxprints.remove('INPUT')
        l_toxprints.remove('DTXSID')
        l_toxprints.remove('PREFERRED_NAME')

        d_out = {}
        for toxprint in l_toxprints:
            d_out[toxprint] = 0

        for chem in d_toxprint.keys():
            # case of genotoxic in breast carcinogen list
            flag = 0
            for chem_dataset in list(self.c_dataset.d_dataset.keys()):
                if "DTXSID" in list(self.c_dataset.d_dataset[chem_dataset].keys()):
                    dtxsid =  d_toxprint[chem]["DTXSID"]
                    if dtxsid == self.c_dataset.d_dataset[chem_dataset]["DTXSID"]:
                        for toxprint in l_toxprints:
                            if d_toxprint[chem][toxprint] == "1":
                                d_out[toxprint] = d_out[toxprint] + 1
                        flag=1
                        break
            if flag == 0:
                for toxprint in l_toxprints:
                    if d_toxprint[chem][toxprint] == "1":
                        d_out[toxprint] = d_out[toxprint] + 1

        p_filout = pr_out + "count_toxprint"
        filout = open(p_filout, "w")
        filout.write("Toxprint\tcount\n")
        for toxprint in d_out.keys():
            if d_out[toxprint] != 0:
                filout.write("%s\t%s\n"%(toxprint, d_out[toxprint]))
        filout.close()

        runExternal.barplotToxPrint(p_filout)

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