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

    def __init__(self, p_crossref, p_exposure, COR_VAL, MAX_QUANTILE, pr_ToxPrints, pr_out):

        self.pr_toxprint = pr_ToxPrints
        self.pr_out = pr_out
        self.p_crossref = p_crossref
        self.p_exposure = p_exposure

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
        self.d_H295R = toolbox.loadExcelSheet(self.p_crossref, name_sheet='H295R_steroid_synth', k_head = "CASRN")

    def updateExposureForMC(self):
        if not "d_MC" in self.__dict__:
            print("ERROR - load excel file with chemical first")
            return 
        
        d_exposure = toolbox.loadMatrix(self.p_exposure, sep=",")
        for chem in self.d_MC.keys():
            self.d_MC[chem]["exposure"] = []

            if chem in list(d_exposure.keys()):
                if d_exposure[chem]["Diet"] != "NA":
                    self.d_MC[chem]["exposure"].append("Diet")
                if d_exposure[chem]["Pharma"] != "NA":
                    self.d_MC[chem]["exposure"].append("Pharma")
                if d_exposure[chem]["Consumer"] != "NA":
                    self.d_MC[chem]["exposure"].append("Consumer")
                if d_exposure[chem]["Pesticide"] != "NA":
                    self.d_MC[chem]["exposure"].append("Pesticide")
                if d_exposure[chem]["Industrial"] != "NA":
                    self.d_MC[chem]["exposure"].append("Industrial")
        
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
            elif list_chem == "H295R":
                for chem in self.d_H295R.keys():
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
                d_all[CASRN]["H295R"] = "0"

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
                d_all[CASRN]["H295R"] = "0"
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
                d_all[CASRN]["H295R"] = "0"
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
                d_all[CASRN]["H295R"] = "0"
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
                d_all[CASRN]["H295R"] = "0"
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
                d_all[CASRN]["H295R"] = "0"
            else:
                d_all[CASRN]["E2up"] = self.d_E2up[CASRN]["Efficacy/potency"]
                d_all[CASRN]["Group"].append("E2up")

        # d_H295R
        for CASRN in self.d_H295R.keys():
            if not CASRN in list(d_all.keys()):
                d_all[CASRN] = {}
                d_all[CASRN]["name"] = self.d_H295R[CASRN]["PREFERRED_NAME"]
                d_all[CASRN]["SMILES"] = self.d_H295R[CASRN]["SMILES"]
                d_all[CASRN]["genotox"] = "NA"
                d_all[CASRN]["ER"] = []
                d_all[CASRN]["P4up"] = "NA"
                d_all[CASRN]["E2up"] = "NA"
                d_all[CASRN]["H295R"] = "1"
                d_all[CASRN]["Group"] = ["H295R"]
                d_all[CASRN]["MC"] = "0"
            else:
                d_all[CASRN]["H295R"] = "1"
                d_all[CASRN]["Group"].append("H295R")

        self.d_all = d_all

    def formatSetofChem(self, l_chemsets):
        pr_dataset = pathFolder.createFolder(self.pr_out + "setOfChemicals/")

        d_out = {}
        for chemset in l_chemsets:
            p_filout = pr_dataset + chemset + ".csv"
            d_out[chemset] = p_filout
            if path.exists(p_filout):
                continue
            filout = open(p_filout, "w")
            filout.write("CASRN\tSMILES\tChemical name\tMC\tgenotox\tE2up\tP4up\tER\tH295R\tGroup\tAff\n")

            if chemset == "MC":
                l_casrn = list(self.d_MC.keys())
            elif chemset == "genotoxic":
                l_casrn = list(self.d_MCgenotox["genotoxic"])
            elif chemset == "Steroid-up":
                l_casrn = list(self.d_steroid_active.keys())
            elif chemset == "Steroid":
                l_casrn = list(self.d_steroid.keys())
            elif chemset == "ER-agonist":
                l_casrn = list(self.d_ERagonist.keys())
            elif chemset == "ER":
                l_casrn = list(self.d_ER.keys())
            elif chemset == "E2":
                l_casrn = list(self.d_E2up.keys())
            elif chemset == "E2-up":
                l_casrn = list(self.d_E2up_active.keys())
            elif chemset == "P4":
                l_casrn = list(self.d_P4up.keys())
            elif chemset == "P4-up":
                l_casrn = list(self.d_P4up_active.keys())
            elif chemset == "H295R":
                l_casrn = list(self.d_H295R.keys())
            elif chemset == "all":
                l_casrn = list(self.d_all.keys())

            for CASRN in l_casrn:
                filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]), self.d_all[CASRN]["H295R"],  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))
            filout.close()

        self.d_dataset = d_out

    def prepSets(self, l_sets):
        self.loadCrossRefExcel()
        self.splitMCGenotox()
        self.defineE2P4active(["higher", "medium", "lower", "borderline"])  # to change E2 - P4 up defintion // hormone substrate is not included in the active -> NA
        self.defineERactive()
        self.mergeAllSets()
        self.updateExposureForMC()

        # write different chemical sets
        self.formatSetofChem(l_sets)#

    def analysisMDescByDataset(self, dataset, l_desc, cor_val, max_q, hclust=0, SOM=0, SOM_size=8):

        pr_out = pathFolder.createFolder("%sAnalysis_%s/%s/"%(self.pr_out, dataset, "-".join(l_desc)))
        
        # build dataset - based on list of descriptor available
        if not "c_Desc" in self.__dict__:
            return "Error - build descriptor before run clustering"

        # create the descriptor set
        p_desc = self.c_Desc.buildDescSet(dataset, l_desc, pr_out)

        #load the analysis class
        c_MakePlot = MakePlots_fromDesc.MakePlots_fromDesc(p_dataset=self.d_dataset[dataset], p_desc=p_desc, pr_out=pr_out, p_opera_all = self.c_Desc.d_desc[dataset]["all OPERA pred"], cor_val=cor_val, max_quantile=max_q)
        
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

    def ChemClassesByMC(self):

        pr_out = pathFolder.createFolder(self.pr_out + "ChemClassInMC/")

        d_count = {}
        for casrn in self.d_MC.keys():
            lchemClass = self.d_MC[casrn]["exposure"]
            if lchemClass == []:
                lchemClass = ["unclassified"]
            for chemClass in lchemClass:
                if not chemClass in list(d_count.keys()):
                    d_count[chemClass] = 0
                d_count[chemClass] = d_count[chemClass] + 1

        p_count = pr_out + "count_chemical_class_MC.csv"
        filout = open(p_count, "w")
        filout.write("Chem_class\tcount\n")
        for chemClass in d_count.keys():
            filout.write("%s\t%s\n"%(chemClass, str(d_count[chemClass])))
        
        filout.close()
        runExternal.barplotChemClass(p_count)

        stop395




    def main(self):

        # prepare set of chemicals - split by list
        self.prepSets(["ER", "MC", "Steroid", "E2-up", "P4-up", "all", "Steroid-up", "ER-agonist", "genotoxic", "H295R"])

        # compute Venn diagram
        #self.overlapBetweenListChem(["MC", "genotoxic", "Steroid-up", "ER-agonist"])
        #self.overlapBetweenListChem(["MC", "Steroid", "ER"])
        #self.overlapBetweenListChem(["MC", "Steroid", "ER-agonist", "genotoxic"])
        #self.overlapBetweenListChem(["MC", "Steroid-up", "ER-agonist", "genotoxic"])
        #self.overlapBetweenListChem(["Steroid", "E2-up", "P4-up"])
        self.overlapBetweenListChem(["E2-up", "P4-up", "H295R"])
        self.overlapBetweenListChem(["MC", "genotoxic", "H295R"])
        stophere
        # analyse class of chemical by MC
        #self.ChemClassesByMC()
        
        # Compute and/or load descriptors by set of chemicals
        pr_desc_by_list = pathFolder.createFolder(self.pr_out + "desc_by_list/")
        self.c_Desc = runDescriptors.runDescriptors(self.d_dataset, self.pr_desc, pr_desc_by_list)
        self.c_Desc.compute_all() # here included all of the descriptors for the full set of chemicals


        # compute similarity with hormone derivative
        d_hormone = {}

        # put png in a different folder
        #pr_png_by_list = pathFolder.createFolder(self.pr_out + "png_by_list/")
        #self.c_Desc.png_by_list(pr_png_by_list)


        # analyze by daataset and set of descriptors #
        ##############################################
        
        # MC #
        ######
        #self.analysisMDescByDataset(dataset="MC", l_desc=["rdkit"], hclust=0, SOM=1, cor_val=self.COR_VAL, max_q=self.MAX_QUANTILE) #hclust

        # Steroid #
        ######
        #self.analysisMDescByDataset(dataset="Steroid", l_desc=["rdkit"], hclust=1, SOM=1, cor_val=self.COR_VAL, max_q=self.MAX_QUANTILE) #hclust


        # H295R #
        self.analysisMDescByDataset(dataset="H295R", l_desc=["rdkit"], hclust=1, SOM=1, cor_val=self.COR_VAL, max_q=self.MAX_QUANTILE) #hclust
        
        #hormone and precursor similarity


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
        #self.c_FP.comparisonToxPrintCount(["Steroid", "E2-up", "P4-up"])
        #self.c_FP.comparisonToxPrintCount(["H295R", "E2-up", "P4-up"])

        stop485



## not used anymore but can be reintegrate 
#######

## need to reintegrate chemlinst in the sources as do with toxprint ##
#cChem.analysisChemList(p_chemList)


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