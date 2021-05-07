from os import path, listdir
from copy import deepcopy
from re import search
import math

import pathFolder
import dataset
import MakePlots_fromDesc
import runDescriptors
import comparisonChemicalLists
import toolbox
import runExternal
import runToxPrint



class MCcrossref:

    def __init__(self, p_crossref, COR_VAL, MAX_QUANTILE, pr_ToxPrints, pr_out):

        self.pr_toxprint = pr_ToxPrints
        self.pr_out = pr_out
        self.p_crossref = p_crossref

        self.COR_VAL = COR_VAL
        self.MAX_QUANTILE = MAX_QUANTILE

        self.pr_desc = pathFolder.createFolder(pr_out + "DESC/")
    
    def loadCrossRefExcel(self, rm_radiation=1):

        self.d_MC = toolbox.loadExcelSheet(self.p_crossref, name_sheet='Updated 2021 MClist', k_head = "CASRN")
        if rm_radiation == 1:
            l_chem = list(self.d_MC.keys())
            nb_chem  = len(l_chem)
            i = 0
            while i < nb_chem:
                if search("radiation", str(self.d_MC[l_chem[i]]["chemical_use SEE CPDAT-EXPOCAST SPREADSHEET"]).lower()):
                    del self.d_MC[l_chem[i]]
                i = i + 1
        self.d_E2up = toolbox.loadExcelSheet(self.p_crossref, name_sheet='E2-up', k_head = "CASRN")
        self.d_P4up = toolbox.loadExcelSheet(self.p_crossref, name_sheet='P4-up', k_head = "CASRN")
        self.d_ER = toolbox.loadExcelSheet(self.p_crossref, name_sheet='Judson ER active', k_head = "CASRN")

    def splitMCGenotox(self):

        d_out = {}
        for chem in self.d_MC.keys():
            genotox = self.d_MC[chem]["Genotoxic_CCRIS/QSAR/ToxValDB"]
            try:
                genotox = float(genotox)
            except:pass
            if type(genotox) == float and math.isnan(genotox) == True:
                genotox = "not tested"
            elif genotox == 'genotoxic radiation':
                genotox = "genotoxic"
            if not genotox in list(d_out.keys()):
                d_out[genotox] = {}
            d_out[genotox][chem] = deepcopy(self.d_MC[chem])
        self.d_MCgenotox = d_out

    def defineE2P4active(self, l_class_active):

        self.d_steroid_active = {}
        self.d_steroid = {}

        self.d_E2up_active = {}
        for chem in self.d_E2up.keys():
            active = self.d_E2up[chem]["Efficacy/potency"]
            if active in l_class_active:
                self.d_E2up_active[chem] = deepcopy(self.d_E2up[chem])
                self.d_steroid_active[chem] = deepcopy(self.d_E2up[chem])
            self.d_steroid[chem] = deepcopy(self.d_E2up[chem])
            self.d_steroid[chem]["Efficacy/potency"] = "E2-%s"%(self.d_steroid[chem]["Efficacy/potency"])

        self.d_P4up_active = {}
        for chem in self.d_P4up.keys():
            active = self.d_P4up[chem]["Efficacy/potency"]
            if active in l_class_active:
                self.d_P4up_active[chem] = deepcopy(self.d_P4up[chem])
                if not chem in list(self.d_steroid_active.keys()):
                   self.d_steroid_active[chem] =  deepcopy(self.d_P4up[chem])


            if chem in list(self.d_steroid.keys()):
                self.d_steroid[chem]["Efficacy/potency"] = "%s--P4-%s"%(self.d_steroid[chem]["Efficacy/potency"], self.d_P4up[chem]["Efficacy/potency"])
            else:
                self.d_steroid[chem] =  deepcopy(self.d_P4up[chem])
                self.d_steroid[chem]["Efficacy/potency"] = "P4-%s"%(self.d_steroid[chem]["Efficacy/potency"])

    def defineERactive(self):
        self.d_ERagonist = {}
        self.d_ERantagonist = {}

        for chem in self.d_ER.keys():
            if float(self.d_ER[chem]["AUC.Agonist"]) >= 0.1:
                self.d_ERagonist[chem] = deepcopy(self.d_ER[chem])
        
        for chem in self.d_ER.keys():
            if float(self.d_ER[chem]["AUC.Antagonist"]) >= 0.1:
                self.d_ERantagonist[chem] = deepcopy(self.d_ER[chem])

    def overlapBetweenListChem(self, l_list_chem):

        pr_out = pathFolder.createFolder(self.pr_out + "OverlapList/")
        pr_results = pathFolder.createFolder(pr_out + "-".join(l_list_chem) + "/")

        # check if file exist to not repete
        p_upset = pr_results + "upset_matrix"
        if path.exists(p_upset):
            return 


        d_d_chem = {}
        l_CASRN = []
        for list_chem in l_list_chem:
            d_d_chem[list_chem] = []
            if list_chem == "MC":
                for chem in self.d_MC.keys():
                   d_d_chem[list_chem].append(chem)
                   l_CASRN.append(chem)
            elif list_chem == "genotoxic" or list_chem == 'not genotoxic' or list_chem == 'not tested':
                for chem in self.d_MCgenotox[list_chem].keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "E2-up":
                for chem in self.d_E2up_active.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "P4-up":
                for chem in self.d_P4up_active.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "Steroid-up":
                for chem in self.d_steroid_active.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "ER-agonist":
                for chem in self.d_ERagonist.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "ER-antagonist":
                for chem in self.d_ERantagonist.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "ER":
                for chem in self.d_ER.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "Steroid":
                for chem in self.d_steroid.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            
        
        f_open = open(p_upset, "w")
        f_open.write("\t" + "\t".join(list(d_d_chem.keys())) + "\n")
        l_CASRN = list(set(l_CASRN))
        for CASRN in l_CASRN:
            l_w = []
            for list_chem in d_d_chem:
                if CASRN in d_d_chem[list_chem]:
                    l_w.append("1")
                else:
                    l_w.append("0")
            f_open.write("%s\t%s\n"%(CASRN, "\t".join(l_w)))
        f_open.close()

        runExternal.upsetPlot(p_upset)

        if len(l_list_chem) <= 4:
            runExternal.vennPlot(p_upset)

    def mergeAllSets(self):

        d_all = {}
        for CASRN in self.d_MC.keys():
            if not CASRN in list(d_all.keys()):
                d_all[CASRN] = {}
                d_all[CASRN]["name"] = self.d_MC[CASRN]["name"]
                d_all[CASRN]["SMILES"] = self.d_MC[CASRN]["SMILES"]
                if self.d_MC[CASRN]["Genotoxic_CCRIS/QSAR/ToxValDB"] == "":
                    d_all[CASRN]["genotox"] = "NA"
                else:
                    d_all[CASRN]["genotox"] = self.d_MC[CASRN]["Genotoxic_CCRIS/QSAR/ToxValDB"]
                d_all[CASRN]["MC"] = "1"
                d_all[CASRN]["ER"] = []
                d_all[CASRN]["P4up"] = "NA"
                d_all[CASRN]["E2up"] = "NA"
                d_all[CASRN]["Group"] = ["MC"]

        ## ER lists    
        for CASRN in self.d_ERagonist.keys():
            if not CASRN in list(d_all.keys()):
                d_all[CASRN] = {}
                d_all[CASRN]["name"] = self.d_ERagonist[CASRN]["Name"]
                d_all[CASRN]["SMILES"] = self.d_ERagonist[CASRN]["SMILES"]
                d_all[CASRN]["genotox"] = "NA"
                d_all[CASRN]["ER"] = ["agonist"]
                d_all[CASRN]["P4up"] = "NA"
                d_all[CASRN]["E2up"] = "NA"
                d_all[CASRN]["Group"] = ["ER"]
                d_all[CASRN]["MC"] = "0"
            else:
                d_all[CASRN]["ER"].append("agonist")
                d_all[CASRN]["Group"].append("ER")
        
        for CASRN in self.d_ERantagonist.keys():
            if not CASRN in list(d_all.keys()):
                d_all[CASRN] = {}
                d_all[CASRN]["name"] = self.d_ERantagonist[CASRN]["Name"]
                d_all[CASRN]["SMILES"] = self.d_ERantagonist[CASRN]["SMILES"]
                d_all[CASRN]["genotox"] = "NA"
                d_all[CASRN]["ER"] = ["antagonist"]
                d_all[CASRN]["P4up"] = "NA"
                d_all[CASRN]["E2up"] = "NA"
                d_all[CASRN]["Group"] = ["ER"]
                d_all[CASRN]["MC"] = "0"
            else:
                d_all[CASRN]["ER"].append("antagonist")
                if not "ER" in d_all[CASRN]["Group"]:
                    d_all[CASRN]["Group"].append("ER")

        for CASRN in self.d_ER.keys():
            if not CASRN in list(d_all.keys()):
                d_all[CASRN] = {}
                d_all[CASRN]["name"] = self.d_ER[CASRN]["Name"]
                d_all[CASRN]["SMILES"] = self.d_ER[CASRN]["SMILES"]
                d_all[CASRN]["genotox"] = "NA"
                d_all[CASRN]["ER"] = ["tested"]
                d_all[CASRN]["P4up"] = "NA"
                d_all[CASRN]["E2up"] = "NA"
                d_all[CASRN]["Group"] = ["ER"]
                d_all[CASRN]["MC"] = "0"
            else:
                if not "ER" in d_all[CASRN]["Group"]:
                    d_all[CASRN]["Group"].append("ER")


        for CASRN in self.d_P4up.keys():
            if not CASRN in list(d_all.keys()):
                d_all[CASRN] = {}
                d_all[CASRN]["name"] = self.d_P4up[CASRN]["Chemical name"]
                d_all[CASRN]["SMILES"] = self.d_P4up[CASRN]["SMILES"]
                d_all[CASRN]["genotox"] = "NA"
                d_all[CASRN]["ER"] = []
                d_all[CASRN]["P4up"] = self.d_P4up[CASRN]["Efficacy/potency"]
                d_all[CASRN]["E2up"] = "NA"
                d_all[CASRN]["Group"] = ["P4up"]
                d_all[CASRN]["MC"] = "0"
            else:
                d_all[CASRN]["P4up"] = self.d_P4up[CASRN]["Efficacy/potency"]
                d_all[CASRN]["Group"].append("P4up")

        for CASRN in self.d_E2up.keys():
            if not CASRN in list(d_all.keys()):
                d_all[CASRN] = {}
                d_all[CASRN]["name"] = self.d_E2up[CASRN]["Chemical name"]
                d_all[CASRN]["SMILES"] = self.d_E2up[CASRN]["SMILES"]
                d_all[CASRN]["genotox"] = "NA"
                d_all[CASRN]["ER"] = []
                d_all[CASRN]["P4up"] = "NA"
                d_all[CASRN]["E2up"] = self.d_E2up[CASRN]["Efficacy/potency"]
                d_all[CASRN]["Group"] = ["E2up"]
                d_all[CASRN]["MC"] = "0"
            else:
                d_all[CASRN]["E2up"] = self.d_E2up[CASRN]["Efficacy/potency"]
                d_all[CASRN]["Group"].append("E2up")

        self.d_all = d_all

    def formatSetofChem(self, l_chemsets = ["ER", "MC", "Steroid", "E2", "P4", "all", "Steroid-up", "ER-agonist"]):
        pr_dataset = pathFolder.createFolder(self.pr_out + "setOfChemicals/")

        d_out = {}
        for chemset in l_chemsets:
            p_filout = pr_dataset + chemset + ".csv"
            d_out[chemset] = p_filout
            if path.exists(p_filout):
                continue
            filout = open(p_filout, "w")
            filout.write("CASRN\tSMILES\tChemical name\tMC\tgenotox\tE2up\tP4up\tER\tGroup\tAff\n")

            if chemset == "MC":
                for CASRN in self.d_MC.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))
            
            if chemset == "genotox":
                for CASRN in self.d_MCgenotox["genotox"].keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))
            
            if chemset == "Steroid-up":
                for CASRN in self.d_steroid_active.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))
            
            if chemset == "Steroid":
                for CASRN in self.d_steroid.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))

            if chemset == "ER-agonist":
                for CASRN in self.d_ERagonist.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))

            if chemset == "ER":
                for CASRN in self.d_ER.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))

            if chemset == "E2":
                for CASRN in self.d_E2up.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))

            if chemset == "P4":
                for CASRN in self.d_P4up.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))

            if chemset == "all":
                for CASRN in self.d_all.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))
            filout.close()

        self.d_dataset = d_out

    def prepSets(self):
        self.loadCrossRefExcel()
        self.splitMCGenotox()
        self.defineE2P4active(["higher", "medium", "lower"])
        self.defineERactive()
        self.mergeAllSets()

        # write different chemical sets
        self.formatSetofChem(["ER", "MC", "Steroid", "E2", "P4", "all", "Steroid-up", "ER-agonist"])#["Steroid-up", "Steroid"])#, "ER-agonist", "Steroid"])

    def analysisMDescByDataset(self, dataset, l_desc, cor_val, max_q, hclust=0, SOM=0, SOM_size=8):

        pr_out = pathFolder.createFolder("%sAnalysis_%s/%s/"%(self.pr_out, dataset, "-".join(l_desc)))
        
        # build dataset - based on list of descriptor available
        if not "c_Desc" in self.__dict__:
            return "Error - build descriptor before run clustering"

        # create the descriptor set
        p_desc = self.c_Desc.buildDescSet(dataset, l_desc, pr_out)

        #load the analysis class
        c_MakePlot = MakePlots.MakePlots(p_dataset=self.d_dataset[dataset], p_desc=p_desc, pr_out=pr_out, p_opera_all = self.c_Desc.d_desc[dataset]["all OPERA pred"], cor_val=cor_val, max_quantile=max_q)
        
        if hclust == 1:
            c_MakePlot.hclusterByProp()


        if SOM == 1:
            c_MakePlot.SOMClustering(nb_cluster=0)
            c_MakePlot.SOMClustering(nb_cluster=SOM_size)
            #c_MakePlot.SOMMapProp() # maybe need to be developed to extract by cluster the percentage of MC for example or other prop

    def analysisToxPrintByDataset(self, dataset, hclust=0):

        pr_out = pathFolder.createFolder("%sAnalysis_%s/ToxPrint/"%(self.pr_out, dataset))
        
        # build dataset - based on list of descriptor available
        if not "c_FP" in self.__dict__:
            return "Error - build descriptor before run clustering"

        # create the descriptor set
        p_FPMatrix = self.c_FP.d_FPMatrix[dataset]

        #load the analysis class
        c_MakePlot = MakePlots_fromDesc.MakePlots_fromDesc(p_dataset=self.d_dataset[dataset], p_FP=p_FPMatrix, pr_out=pr_out)
        
        if hclust == 1:
            c_MakePlot.hclusterFromFPByProp()





    # to del but need to put back the other analysis in the new class
    def computeDescAnalysisFromList(p_list, p_ToxPrint, p_chemList, PR_RESULTS):

        """
        TO DELETE
        
        """

        pr_results = pathFolder.createFolder(PR_RESULTS + "analysis_individual-dataset/" + p_list.split("/")[-1][0:-4] + "/")
        pr_desc = pathFolder.createFolder(PR_RESULTS + "DESC/")

        # Compute desc
        ###################
        cChem = Chemicals.Chemicals(p_list, pr_results)
        cChem.computeDesc(pr_desc)
        #cChem.buildDescSet(["rdkit"])
        #cChem.analysisDescBasedData(COR_VAL, MAX_QUANTILE, PCA=1, Hclust=1, clustering=1, FP=1, SOM=1, histDesc=1)

        cChem.analysisToxPrint(p_ToxPrint)
        cChem.analysisChemList(p_chemList)







    def analysisMultiSets(p_desc1D2D):

        # draw hclust circular and PCA
        runExternal.multiSetsAnalysis(p_desc1D2D)

    def unionListChem(l_p_set, name_setchem, PR_RESULTS):
        pr_union = pathFolder.createFolder(PR_RESULTS + "union_sets/")
        l_name_set = []
        d_out = {}
        for p_set in l_p_set:
            name_set = p_set.split("/")[-1][0:-4]
            l_name_set.append(name_set)
            l_chemset = toolbox.loadMatrixToList(p_set, sep = ",")
            for chemset in l_chemset:
                print(chemset)
                CASRN = chemset["CASRN"]
                SMILES = chemset["SMILES"]
                if not CASRN in list(d_out.keys()):
                    d_out[CASRN] = {}
                    d_out[CASRN]["SMILES"] = SMILES
                    d_out[CASRN]["list"] = []
                d_out[CASRN]["list"].append(name_set)

        p_filout = pr_union + name_setchem + ".csv"
        filout = open(p_filout, "w")
        filout.write("CASRN,SMILES,list\n")
        for chem in d_out.keys():
            filout.write("%s,%s,%s\n"%(chem, d_out[chem]["SMILES"], "--".join(d_out[chem]["list"])))
        filout.close()

        return p_filout

    def main(self):

        # prepare set of chemicals - split by list
        self.prepSets()
        
        # compute Venn diagram
        #self.overlapBetweenListChem(["MC", "genotoxic", "Steroid-up", "ER-agonist"])
        #self.overlapBetweenListChem(["MC", "Steroid", "ER"])
        #self.overlapBetweenListChem(["MC", "Steroid", "ER-agonist", "genotoxic"])
        #self.overlapBetweenListChem(["MC", "Steroid-up", "ER-agonist", "genotoxic"])
        
        # Compute and/or load descriptors by set of chemicals
        #pr_desc_by_list = pathFolder.createFolder(self.pr_out + "desc_by_list/")
        #self.c_Desc = runDescriptors.runDescriptors(self.d_dataset, self.pr_desc, pr_desc_by_list)
        #self.c_Desc.compute_all() # here included all of the descriptors for the full set of chemicals

        # analyze by daataset and set of descriptors #
        ##############################################
        
        # MC #
        ######
        #self.analysisMDescByDataset(dataset="MC", l_desc=["rdkit"], hclust=0, SOM=1, cor_val=self.COR_VAL, max_q=self.MAX_QUANTILE) #hclust

        # E2 #
        ######
        #self.analysisMDescByDataset(dataset="E2", l_desc=["rdkit"], hclust=1, SOM=1, cor_val=self.COR_VAL, max_q=self.MAX_QUANTILE) #hclust

        # P4 #
        ######
        #self.analysisMDescByDataset(dataset="P4", l_desc=["rdkit"], hclust=1, SOM=1, cor_val=self.COR_VAL, max_q=self.MAX_QUANTILE) #hclust



        # analysis of the toxprint #
        ############################
        #self.c_FP = runToxPrint.runToxPrint(self.d_dataset, self.pr_toxprint, self.pr_out)
        #self.c_FP.loadToxPrint()
        #self.c_FP.ToxPrintCount()
        #self.c_FP.computeTanimotoMatrix()


        # MC # 
        ######
        #self.analysisToxPrintByDataset(dataset="MC", hclust=1) #hclust


        # Comparison ToxPrint #
        #######################
        #self.c_FP.comparisonToxPrintCount(["MC", "Steroid", "all"])

     