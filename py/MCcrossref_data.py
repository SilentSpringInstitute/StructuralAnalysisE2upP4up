from os import path, listdir
from copy import deepcopy
from shutil import copyfile
from re import search
import math

import pathFolder
import dataset
import MakePlots_fromDesc
import runChemStruct
import comparisonChemicalLists
import toolbox
import runExternal
import runToxPrint
import ToxCast


class MCcrossref:

    def __init__(self, p_crossref, p_exposure, p_hormones, COR_VAL, MAX_QUANTILE, pr_ToxPrints, pr_root):

        self.pr_toxprint = pr_ToxPrints
        self.pr_out = pr_root + "results/"
        self.pr_data = pr_root + "data/"
        self.p_crossref = p_crossref
        self.p_exposure = p_exposure
        self.p_hormones = p_hormones

        self.COR_VAL = COR_VAL
        self.MAX_QUANTILE = MAX_QUANTILE

        self.pr_desc = pathFolder.createFolder(self.pr_out + "DESC/")
    
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

        # load E2up need to load 2 sheet to have SMILES        
        # 1- load AC50 from bethsaida EHP
        self.d_E2up = toolbox.loadExcelSheet(self.p_crossref, name_sheet='E2-up', k_head = "CASRN")
        self.d_P4up = toolbox.loadExcelSheet(self.p_crossref, name_sheet='P4-up', k_head = "CASRN")

        # 2- load AC50
        d_E2up_temp = toolbox.loadExcelSheet(self.p_crossref, name_sheet='H295R_E2up', k_head = "CASN_Protect")
        d_P4up_temp = toolbox.loadExcelSheet(self.p_crossref, name_sheet='H295R_P4up', k_head = "CASN_Protect")
        
        toolbox.combineDict(self.d_E2up, d_E2up_temp)
        toolbox.combineDict(self.d_P4up, d_P4up_temp)

        # for ER and h295 stereo only take the sheet with the right name
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

        pr_dataset = pathFolder.createFolder(self.pr_out + "setOfChemicals/")
        
        d_sum = {}
        d_sum["E2up"] = {}
        d_sum["P4up"] = {}

        p_filout_sum = pr_dataset + "E2up_P4up.sum"


        self.d_steroid_active = {}
        self.d_steroid = {}

        self.d_E2up_active = {}
        for chem in self.d_E2up.keys():
            active = self.d_E2up[chem]["Efficacy/potency"]
            if not active in list(d_sum["E2up"].keys()):
                d_sum["E2up"][active] = 0
            d_sum["E2up"][active] = d_sum["E2up"][active] + 1
            if active in l_class_active:
                self.d_E2up_active[chem] = deepcopy(self.d_E2up[chem])
                self.d_steroid_active[chem] = deepcopy(self.d_E2up[chem])
            self.d_steroid[chem] = deepcopy(self.d_E2up[chem])
            self.d_steroid[chem]["Efficacy/potency"] = "E2-%s"%(self.d_steroid[chem]["Efficacy/potency"])

        self.d_P4up_active = {}
        for chem in self.d_P4up.keys():
            active = self.d_P4up[chem]["Efficacy/potency"]
            if not active in list(d_sum["P4up"].keys()):
                d_sum["P4up"][active] = 0
            d_sum["P4up"][active] = d_sum["P4up"][active] + 1
            if active in l_class_active:
                self.d_P4up_active[chem] = deepcopy(self.d_P4up[chem])
                if not chem in list(self.d_steroid_active.keys()):
                   self.d_steroid_active[chem] =  deepcopy(self.d_P4up[chem])


            if chem in list(self.d_steroid.keys()):
                self.d_steroid[chem]["Efficacy/potency"] = "%s--P4-%s"%(self.d_steroid[chem]["Efficacy/potency"], self.d_P4up[chem]["Efficacy/potency"])
            else:
                self.d_steroid[chem] =  deepcopy(self.d_P4up[chem])
                self.d_steroid[chem]["Efficacy/potency"] = "P4-%s"%(self.d_steroid[chem]["Efficacy/potency"])

        # write the sum
        filout = open(p_filout_sum, "w")
        filout.write("Type\tE2up\tP4up\n")
        l_active = list(d_sum[list(d_sum.keys())[0]].keys())
        for active in l_active:
            filout.write("%s\t%s\t%s\n"%(active, d_sum["E2up"][active], d_sum["P4up"][active]))
        filout.close()

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
        """Compute overlap between different list and draw Venn diagram"""

        pr_out = pathFolder.createFolder(self.pr_out + "OverlapList/")
        pr_results = pathFolder.createFolder(pr_out + "-".join(l_list_chem) + "/")

        # check if file exist to not repete
        p_upset = pr_results + "upset_matrix"
        #if path.exists(p_upset):
        #    return p_upset

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
            elif list_chem == "E2up":
                for chem in self.d_E2up_active.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "P4up":
                for chem in self.d_P4up_active.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "Steroid-up":
                for chem in self.d_steroid_active.keys():
                    d_d_chem[list_chem].append(chem)
                    l_CASRN.append(chem)
            elif list_chem == "ERagonist":
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
            for list_chem in d_d_chem.keys():
                if CASRN in d_d_chem[list_chem]:
                    l_w.append("1")
                else:
                    l_w.append("0")
            f_open.write("%s\t%s\n"%(CASRN, "\t".join(l_w)))
        f_open.close()

        runExternal.upsetPlot(p_upset)

        if len(l_list_chem) <= 4:
            runExternal.vennPlot(p_upset)
        
        d_upset = toolbox.loadMatrix(p_upset)

        # save overlap
        pr_overlap = pathFolder.createFolder(pr_results + "overlapChem/")
        d_upset = toolbox.loadMatrix(p_upset, sep = "\t")
        l_list_chem = list(d_upset[list(d_upset.keys())[0]].keys())

        p_filout = pr_overlap + "overlap_chem.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\tChemical.name\n")

        for chem in d_upset.keys():
            flag = 0
            for list_chem in l_list_chem:
                if d_upset[chem][list_chem] == "0":
                    flag = 1
            if flag == 1:
                continue
            filout.write("%s\t%s\n"%(chem, self.d_all[chem]["name"]))
            # copy paste png chemicals
            try:copyfile("%sPNG/%s.png"%(self.pr_desc, chem), "%s%s.png"%(pr_overlap, chem))
            except: pass
        filout.close()

        return p_upset

    def overlapBetweenListChemWithHormoneSim(self, l_list_chem):
        """Overlap with 2 lists max"""

        pr_out = pathFolder.createFolder(self.pr_out + "OverlapList/")
        pr_results = pathFolder.createFolder(pr_out + "-".join(l_list_chem) + "/")

        # check if file exist to not repete
        p_upset = pr_results + "upset_matrix"
        if not path.exists(p_upset):
            self.overlapBetweenListChem(l_list_chem)

        # comparison hormone similarity overlap
        runExternal.overlapListHormoneSim(p_upset, self.c_Desc.p_hormone_similarity, pr_results)

        # extract structure from overlap
        pr_overlap = pathFolder.createFolder(pr_results + "overlapChem/")
        d_upset = toolbox.loadMatrix(p_upset)

        p_filout = pr_overlap + "overlap_chem.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\tChemical.name\n")

        for chem in d_upset.keys():
            flag = 0
            for list_chem in l_list_chem:
                if d_upset[chem][list_chem] == "0":
                    flag = 1
            if flag == 1:
                continue
            filout.write("%s\t%s\n"%(chem, self.d_all[chem]["name"]))
            # copy paste png chemicals
            try:copyfile("%sPNG/%s.png"%(self.pr_desc, chem), "%s%s.png"%(pr_overlap, chem))
            except: pass
        filout.close()

    def overlapListChemE2upP4up(self):
        pr_out = pathFolder.createFolder(self.pr_out + "OverlapList/")
        pr_results = pathFolder.createFolder(pr_out + "E2_P4_potency/")

        d_chem = {}
        for chem in self.d_E2up_active.keys():
            if not chem in list(d_chem.keys()):
                d_chem[chem] = {"E2up_lower": 0, "E2up_medium": 0, "E2up_higher": 0, "P4up_lower": 0, "P4up_medium": 0, "P4up_higher":0}
            if self.d_E2up_active[chem]["Efficacy/potency"] == "medium":
                d_chem[chem]["E2up_medium"] = 1
            elif self.d_E2up_active[chem]["Efficacy/potency"] == "higher":
                d_chem[chem]["E2up_higher"] = 1
            elif self.d_E2up_active[chem]["Efficacy/potency"] == "lower":
                 d_chem[chem]["E2up_lower"] = 1
        
        for chem in self.d_P4up_active.keys():
            if not chem in list(d_chem.keys()):
                d_chem[chem] = {"E2up_lower": 0, "E2up_medium": 0, "E2up_higher": 0, "P4up_lower": 0, "P4up_medium": 0, "P4up_higher":0}
            if self.d_P4up_active[chem]["Efficacy/potency"] == "medium":
                d_chem[chem]["P4up_medium"] = 1
            elif self.d_P4up_active[chem]["Efficacy/potency"] == "higher":
                d_chem[chem]["P4up_higher"] = 1
            elif self.d_P4up_active[chem]["Efficacy/potency"] == "lower":
                d_chem[chem]["P4up_lower"] = 1

        p_filout = pr_results + "upset_plot"
        filout = open(p_filout, "w")
        filout.write("CARSN\tE2up_lower\tE2up_medium\tE2up_higher\tP4up_lower\tP4up_medium\tP4up_higher\n")
        for chem in d_chem.keys():
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(chem, d_chem[chem]["E2up_lower"],d_chem[chem]["E2up_medium"],d_chem[chem]["E2up_higher"], d_chem[chem]["P4up_lower"], d_chem[chem]["P4up_medium"], d_chem[chem]["P4up_higher"]))
        
            # extract png overlap
            if d_chem[chem]["E2up_medium"] == 1 and d_chem[chem]["P4up_medium"] == 1:
                pr_med = pathFolder.createFolder(pr_results + "E2up_medium__P4up_medium/")
                try:copyfile("%sPNG/%s.png"%(self.pr_desc, chem), "%s%s.png"%(pr_med, chem))
                except: pass

            if d_chem[chem]["E2up_higher"] == 1 and d_chem[chem]["P4up_higher"] == 1:
                pr_higher = pathFolder.createFolder(pr_results + "E2up_higher__P4up_higher/")
                try:copyfile("%sPNG/%s.png"%(self.pr_desc, chem), "%s%s.png"%(pr_higher, chem))
                except: pass

        filout.close()

        runExternal.upsetPlot(p_filout)


        # comparison hormone similarity overlap
        runExternal.overlapListHormoneSim(p_filout, self.c_Desc.p_hormone_similarity, pr_results)

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
                d_all[CASRN]["H295R"] = "NA"

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
                d_all[CASRN]["MC"] = "NA"
                d_all[CASRN]["H295R"] = "NA"
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
                d_all[CASRN]["MC"] = "NA"
                d_all[CASRN]["H295R"] = "NA"
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
                d_all[CASRN]["MC"] = "NA"
                d_all[CASRN]["H295R"] = "NA"
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
                d_all[CASRN]["MC"] = "NA"
                d_all[CASRN]["H295R"] = "NA"
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
                d_all[CASRN]["MC"] = "NA"
                d_all[CASRN]["H295R"] = "NA"
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
                d_all[CASRN]["MC"] = "NA"
            else:
                d_all[CASRN]["H295R"] = "1"
                d_all[CASRN]["Group"].append("H295R")

        self.d_all = d_all

    def formatSetofChem(self, l_chemsets):
        """
        ER add in the format list but not only agonist
        """

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
            elif chemset == "ERagonist":
                l_casrn = list(self.d_ERagonist.keys())
            elif chemset == "ER":
                l_casrn = list(self.d_ER.keys())
            elif chemset == "E2":
                l_casrn = list(self.d_E2up.keys())
            elif chemset == "E2up":
                l_casrn = list(self.d_E2up_active.keys())
            elif chemset == "P4":
                l_casrn = list(self.d_P4up.keys())
            elif chemset == "P4up":
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
        """
        Prepare list of chemicals 
        """
        self.loadCrossRefExcel()
        self.splitMCGenotox()
        self.defineE2P4active(["higher", "medium", "lower"])#, "borderline"])  # to change E2 - P4 up defintion // hormone substrate is not included in the active -> NA
        self.defineERactive()
        self.mergeAllSets()
        self.updateExposureForMC()

        # write different chemical sets
        self.formatSetofChem(l_sets)#

    def analysisMDescByDataset(self, dataset, l_desc, cor_val, max_q, hclust=0, SOM=0, SOM_size=12):

        pr_out = pathFolder.createFolder("%sAnalysis_%s/%s/"%(self.pr_out, dataset, "-".join(l_desc)))
        
        # build dataset - based on list of descriptor available
        if not "c_Desc" in self.__dict__:
            return "Error - build descriptor before run clustering"

        # create the descriptor set
        p_desc = self.c_Desc.buildDescSet(dataset, l_desc, pr_out)

        #load the analysis class
        c_MakePlot = MakePlots_fromDesc.MakePlots_fromDesc(p_dataset=self.d_dataset[dataset], p_desc=p_desc, pr_desc=self.pr_desc, p_hormone_similarity = self.c_Desc.p_hormone_similarity, pr_out=pr_out, p_opera_all = self.c_Desc.d_desc[dataset]["all OPERA pred"], cor_val=cor_val, max_quantile=max_q)
        
        if hclust == 1:
            c_MakePlot.hclusterByProp()

        if SOM == 1:
            # check if model already exist
            pr_SOM = pathFolder.createFolder(pr_out + "SOM/")
            p_model = pr_SOM + "SOM_model.RData"
            if not path.exists(p_model):
                c_MakePlot.SOMClustering(nb_cluster=0)
            
            c_MakePlot.SOMClustering(nb_cluster=SOM_size)
            c_MakePlot.SOMMapProp(self.d_dataset[dataset]) # maybe need to be developed to extract by cluster the percentage of MC for example or other prop
            
            c_MakePlot.SOMHormoneSimilarity()

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

    def ComparisonTwoChemicalLists(self, l_datasets):

        pr_out = pathFolder.createFolder(self.pr_out + "comparisonDescToxprint_" + "-".join(l_datasets) + "/")

        
        # by descriptors
        if not "c_Desc" in self.__dict__:
            print("Compute Descriptor first")
            return 
        
        # with molecular descriptors - quantitative with 
        pr_rdkit = pathFolder.createFolder(pr_out + "rdkit/")
        pr_opera = pathFolder.createFolder(pr_out + "opera/")
        runExternal.comparisonDesc(self.c_Desc.d_desc[l_datasets[0]]["rdkit"], self.c_Desc.d_desc[l_datasets[1]]["rdkit"], pr_rdkit)
        runExternal.comparisonDesc(self.c_Desc.d_desc[l_datasets[0]]["OPERA"], self.c_Desc.d_desc[l_datasets[1]]["OPERA"], pr_opera)

        # toxprint as descriptor - comparison qualitative proportion
        pr_toxprint = pathFolder.createFolder(pr_out + "Toxprint/")
        p_toxprint1 = self.c_FP.writeToxPrintMatrix(list(toolbox.loadMatrix(self.d_dataset[l_datasets[0]]).keys()), pr_out + l_datasets[0] + "_toxprint.csv")
        p_toxprint2 = self.c_FP.writeToxPrintMatrix(list(toolbox.loadMatrix(self.d_dataset[l_datasets[1]]).keys()), pr_out + l_datasets[1] + "_toxprint.csv")
        runExternal.comparisonToxPrint(p_toxprint1, p_toxprint2, pr_toxprint)

        # with the similarity with the hormone
        pr_hormone = pathFolder.createFolder(pr_out + "similarity_hormone/")
        runExternal.comparisonWithHormoneSimilarity(self.d_dataset[l_datasets[0]], self.d_dataset[l_datasets[1]], self.c_Desc.p_hormone_similarity, pr_hormone)

    def crossToxCastAssays(self, l_genes = ["CYP19A1"]):
        pr_results = pathFolder.createFolder(self.pr_out + "overlapToxCast/")
        
        p_assays = pr_results + "".join(l_genes) + ".csv"
        if not path.exists(p_assays):
            c_ToxCast = ToxCast.ToxCast(l_genes, pr_results)
            p_assays = c_ToxCast.loadAssaysMatrix()
        

        # compte overlap E2 - P4 up
        p_overlap = self.overlapBetweenListChem(["E2up", "P4up"])
        
        # overlap
        runExternal.overlapAssaysListChem(p_assays, p_overlap, pr_results)

        p_upset = pr_results + "upset_all.csv"

        # extract structure from overlap
        pr_overlap = pathFolder.createFolder(pr_results + "overlapChem/")
        d_upset = toolbox.loadMatrix(p_upset, sep = ",")
        l_list_chem = list(d_upset[list(d_upset.keys())[0]].keys())

        p_filout = pr_overlap + "overlap_chem.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\tChemical.name\n")

        for chem in d_upset.keys():
            flag = 0
            for list_chem in l_list_chem:
                if d_upset[chem][list_chem] == "1":
                    flag = 1
            if flag != 1:
                continue
            try:filout.write("%s\t%s\n"%(chem, self.d_all[chem]["name"]))
            except: pass
            # copy paste png chemicals
            try:copyfile("%sPNG/%s.png"%(self.pr_desc, chem), "%s%s.png"%(pr_overlap, chem))
            except: pass
        filout.close()

    def AC50ByList(self, aenm, l_list, l_class_active = []):
        """
        Make radial plot with AC50 for selected assays and E2up - P4up AC50 or AC10
        """

        pr_out = pathFolder.createFolder(self.pr_out + "Assays_" + aenm + "-".join(l_list) + '/')
        
        # load ac50 for the selected aenm
        p_ac50 = self.pr_data + aenm + "/ac50.csv"
        if path.exists(p_ac50):
            d_ac50 = toolbox.loadMatrix(p_ac50)
        else:
            return
        
        d_out = {}
        for dtxsid in d_ac50.keys():
            casrn = d_ac50[dtxsid]["casn"]
            ac50 = d_ac50[dtxsid]["ac50"]
            name = d_ac50[dtxsid]["name"]
            QC = d_ac50[dtxsid]["new_hitc"]
            if ac50 != "" and casrn in list(self.d_all.keys()) and QC == "1":
                #flag_in = 0
                #for list_chem in l_list:
                #    if self.d_all[casrn][list_chem] == 1 or self.d_all[casrn][list_chem] in l_class_active:
                #        flag_in = 1
                #        break
                
                #if flag_in == 1:
                    d_out[casrn] = {}
                    d_out[casrn]["ac50"] = ac50
                    d_out[casrn]["name"] = name
                    d_out[casrn]["dtxsid"] = dtxsid
                    for list_chem in l_list:
                        if self.d_all[casrn][list_chem] == 1 or self.d_all[casrn][list_chem] in l_class_active:
                            d_out[casrn][list_chem] = 1
                        else:
                            d_out[casrn][list_chem] = 0
                #d_out[casrn] = {}


        # add in it the variation of the hormone level from the h295R
        pr_assays = self.pr_data + "H295R_assays/"
        l_p_assays = listdir(pr_assays)
        l_horm_delta = []
        for p_assay in l_p_assays:
            d_assays = toolbox.loadMatrix(pr_assays + p_assay)
            horm_delta = p_assay[:-4]
            l_horm_delta.append(horm_delta)
            for k_casn in d_out.keys():
                dtxsid = d_out[k_casn]["dtxsid"]
                
                if dtxsid in list(d_assays.keys()) and d_assays[dtxsid]["new_hitc"] == "1":
                    d_out[k_casn][horm_delta] =  d_assays[dtxsid]["ac50"]
                else:
                    d_out[k_casn][horm_delta] =  "NA"


        # add AC50 + AC10 for E2up and P4up

        l_kac = ["AC50_E2up", "AC10_E2up", "AC50_P4up", "AC10_P4up"]
        for k_casn in d_out.keys():
            try:d_out[k_casn]["AC50_E2up"] = self.d_E2up[k_casn]["AC50 (µM)"]
            except: d_out[k_casn]["AC50_E2up"] = "NA"
            try:d_out[k_casn]["AC10_E2up"] = self.d_E2up[k_casn]["AC10 (µM)"]
            except: d_out[k_casn]["AC10_E2up"] = "NA"
            try:d_out[k_casn]["AC50_P4up"] = self.d_P4up[k_casn]["AC50 (µM)"]
            except: d_out[k_casn]["AC50_P4up"] = "NA"
            try:d_out[k_casn]["AC10_P4up"] = self.d_P4up[k_casn]["AC10 (µM)"]
            except: d_out[k_casn]["AC10_P4up"] = "NA"

        p_filout = pr_out + "ac50_list.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\tname\tAC50\t" + "\t".join(l_list) +"\t" + "\t".join(l_horm_delta) + "\t" + "\t".join(l_kac) + "\n")
        for casrn in d_out.keys():
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(casrn, d_out[casrn]["name"], d_out[casrn]["ac50"], "\t".join([str(d_out[casrn][inlist]) for inlist in l_list]), "\t".join([str(d_out[casrn][delta_horm]) for delta_horm in l_horm_delta]), "\t".join(["0" if search("not active",str(d_out[casrn][kac])) else str(d_out[casrn][kac]) for kac in l_kac])))
        filout.close()

        # draw radial plot with AC50
        stophere
        runExternal.drawRadialPlotAC50(p_filout, pr_out)

    def corHormoneSimilarityClassActive(self, dataset):
        """
        only available for E2up and P4up
        """
        pr_out = pathFolder.createFolder(self.pr_out + "corrSimHormActiveClass/")

        if dataset == "E2up":
            p_subfilout = pr_out + "E2up_cor"
            d_chem = self.d_E2up_active
        elif dataset == "P4up":
            p_subfilout = pr_out + "P4up_cor"
            d_chem = self.d_P4up_active
        else:
            return

        d_sim = toolbox.loadMatrix(self.c_Desc.p_hormone_similarity)
        l_horm = ["521-18-6", "53-41-8", "68-96-2", "57-85-2", "57-83-0", "63-05-8", "53-43-0", "58-18-4", "50-28-2", "57-91-0", "53-16-7", "474-86-2", "50-27-1", "57-63-6", "57-88-5", "145-13-1", "387-79-1", "50-23-7","50-22-6","152-58-9","64-85-7"]
        for horm in l_horm:
            p_filout = p_subfilout + "_" + horm 
            filout = open(p_filout, "w")
            filout.write("CASRN\tSimilarity\tEfficacy/potency\n")
            for chem in d_chem.keys():
                if not chem in list(d_sim.keys()):
                    continue
                filout.write("%s\t%s\t%s\n"%(chem, d_sim[chem][horm], d_chem[chem]["Efficacy/potency"]))
            filout.close()

            runExternal.corHormSimClassActive(p_filout,dataset)

    def main(self):

        # prepare set of chemicals - split by list
        self.prepSets(["ER", "MC", "Steroid", "E2up", "P4up", "all", "Steroid-up", "ERagonist", "genotoxic", "H295R"])
        # compute Venn diagram - overlap between list of chemicals
        ######
        #self.overlapBetweenListChem(["MC", "genotoxic", "Steroid-up", "ER-agonist"])
        #self.overlapBetweenListChem(["MC", "Steroid", "ER"])
        #self.overlapBetweenListChem(["MC", "Steroid", "ER-agonist", "genotoxic"])
        #self.overlapBetweenListChem(["MC", "Steroid-up", "ER-agonist", "genotoxic"])
        #self.overlapBetweenListChem(["Steroid", "E2up", "P4up"])
        #self.overlapBetweenListChem(["E2up", "P4up", "H295R"])
        #self.overlapBetweenListChem(["MC", "genotoxic", "H295R"])
        #self.overlapBetweenListChem(["E2up", "P4up"])
        #self.overlapBetweenListChem(["MC", "E2up", "P4up"])
        #self.overlapBetweenListChem(["MC", "E2up", "P4up", "H295R"])

        
        # analyse class of chemical by MC
        ####
        #self.ChemClassesByMC()


        # Compute and/or load descriptors by set of chemicals
        ##################
        self.c_Desc = runChemStruct.runChemStruct(self.d_dataset, self.pr_desc, self.pr_out)
        self.c_Desc.compute_allDesc() # here included all of the descriptors for the full set of chemicals

        # compute similarity with hormone derivative
        ########################
        # FP types: "MACCS", "Morgan", "Mol"
        # Dist types: "Dice", "Tanimoto"
        #self.c_Desc.compute_similarity_inter_hormones(self.p_hormones)
        #self.c_Desc.compute_similarity_with_hormones(self.p_hormones, "MACCS", "Tanimoto")
        
        # correlation similarity with eff/pot class for E2-P4 up
        ###############################
        #self.corHormoneSimilarityClassActive("E2up")
        #self.corHormoneSimilarityClassActive("P4up")
        
        # overlap E2up - P4up with class Efficacy/Efficiency
        ######
        #self.overlapListChemE2upP4up()


        # put png in a different folder
        #############################
        #pr_png_by_list = pathFolder.createFolder(self.pr_out + "png_by_list/")
        #self.c_Desc.png_by_list(pr_png_by_list)


        # load Toxprint by list
        ##########################
        #self.c_FP = runToxPrint.runToxPrint(self.d_dataset, self.pr_toxprint, self.pr_out)
        #self.c_FP.loadToxPrint()

        # overlap with similarity with hormone
        #######
        #self.overlapBetweenListChemWithHormoneSim(["E2up", "P4up", "H295R"])

        # analyze by daataset and set of descriptors #
        ##############################################
        
        # MC #
        ######
        #self.analysisMDescByDataset(dataset="MC", l_desc=["rdkit"], hclust=0, SOM=1, cor_val=self.COR_VAL, max_q=self.MAX_QUANTILE) #hclust

        # Steroid #
        ######
        #self.analysisMDescByDataset(dataset="Steroid-up", l_desc=["rdkit", "OPERA"], hclust=1, SOM=1, cor_val=self.COR_VAL, max_q=self.MAX_QUANTILE) #hclust

        # H295R #
        ######
        self.analysisMDescByDataset(dataset="H295R", l_desc=["rdkit", "OPERA"], hclust=0, SOM=1, cor_val=self.COR_VAL, max_q=self.MAX_QUANTILE, SOM_size=8) #hclust


        # Comparison between two list of chemicals descriptor
        #####
        #self.ComparisonTwoChemicalLists(["E2up", "H295R"])
        #self.ComparisonTwoChemicalLists(["P4up", "H295R"])
        
        #self.ComparisonTwoChemicalLists(["E2up", "P4up"])

        # analysis of the toxprint #
        ############################

        #self.c_FP.ToxPrintCount()
        #self.c_FP.computeTanimotoMatrix()

        # MC # 
        ######
        #self.analysisToxPrintByDataset(dataset="MC", hclust=1) #hclust


        # Comparison ToxPrint #
        #######################
        #self.c_FP.comparisonToxPrintCount(["Steroid", "E2up", "P4up"])
        #self.c_FP.comparisonToxPrintCount(["H295R", "E2up", "P4up"])


        # cross with ToxCast assays #
        #############################

        #self.crossToxCastAssays(l_genes = ["CYP19A1"])
        self.AC50ByList("TOX21_Aromatase_Inhibition", ["E2up", "P4up"], ["higher", "medium", "lower"])
        #self.AC50ByList("NVS_ADME_hCYP19A1", ["E2up", "P4up"], ["higher", "medium", "lower"])


