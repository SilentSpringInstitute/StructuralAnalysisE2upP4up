from os import path, listdir
from copy import deepcopy
from re import search
import math

import pathFolder
import dataset
import analysis
import Chemicals
import comparisonChemicalLists
import toolbox
import runExternal






class Mcarcinogen:

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
        
    def loadToxPrint(self):
        d_toxprint = {}
        l_p_ToxPrint = listdir(self.pr_toxprint)
        for p_ToxPrint in l_p_ToxPrint:
            if search("ToxPrint", p_ToxPrint):
                d_temp = toolbox.loadMatrix(self.pr_toxprint + p_ToxPrint, sep=",")
                d_toxprint.update(d_temp)
        self.d_toxprint = d_toxprint

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
            
        p_upset = pr_results + "upset_matrix"
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
            if chemset == "MC":
                filout.write("CASRN\tDTXSID\tSMILES\tGroup\tAff\n")
                for CASRN in self.d_MC.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_MC[CASRN]["DTXSID"],self.d_MC[CASRN]["SMILES"], self.d_MC[CASRN]["Genotoxic_CCRIS/QSAR/ToxValDB"], self.d_MC[CASRN]["Genotoxic_CCRIS/QSAR/ToxValDB"]))
            
            if chemset == "genotox":
                filout.write("CASRN\tDTXSID\tSMILES\tGroup\tAff\n")
                for CASRN in self.d_MCgenotox["genotox"].keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_MC[CASRN]["DTXSID"],self.d_MC[CASRN]["SMILES"], self.d_MC[CASRN]["Genotoxic_CCRIS/QSAR/ToxValDB"], self.d_MC[CASRN]["Genotoxic_CCRIS/QSAR/ToxValDB"]))
            
            if chemset == "Steroid-up":
                filout.write("CASRN\tChemical name\tSMILES\tGroup\tAff\n")
                for CASRN in self.d_steroid_active.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_steroid_active[CASRN]["Chemical name"],  self.d_steroid_active[CASRN]["SMILES"], self.d_steroid_active[CASRN]["Efficacy/potency"], self.d_steroid_active[CASRN]["LEC (µM)"]))
            
            if chemset == "Steroid":
                filout.write("CASRN\tChemical name\tGroup\tSMILES\tAff\n")
                for CASRN in self.d_steroid.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_steroid[CASRN]["Chemical name"], self.d_steroid[CASRN]["Efficacy/potency"], self.d_steroid[CASRN]["SMILES"], self.d_steroid[CASRN]["LEC (µM)"]))

            if chemset == "ER-agonist":
                filout.write("CASRN\tSMILES\tChemical name\tGroup\tAff\n")
                for CASRN in self.d_ERagonist.keys():
                    filout.write("%s\t%s\t%s\tactive\t%s\n"%(CASRN, self.d_ERagonist[CASRN]["SMILES"], self.d_ERagonist[CASRN]["Name"], self.d_ERagonist[CASRN]["AUC.Agonist"]))

            if chemset == "ER":
                filout.write("CASRN\tSMILES\tChemical name\tGroup\tAff\n")
                for CASRN in self.d_ER.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_ER[CASRN]["SMILES"], self.d_ER[CASRN]["Name"], self.d_ER[CASRN]["structure_category"], self.d_ER[CASRN]["AUC.Agonist"]))

            if chemset == "E2":
                filout.write("CASRN\tSMILES\tChemical name\tGroup\tAff\n")
                for CASRN in self.d_E2up.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_E2up[CASRN]["SMILES"], self.d_E2up[CASRN]["Chemical name"], self.d_E2up[CASRN]["Efficacy/potency"], self.d_E2up[CASRN]["LEC (µM)"]))

            if chemset == "P4":
                filout.write("CASRN\tSMILES\tChemical name\tGroup\tAff\n")
                for CASRN in self.d_P4up.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_P4up[CASRN]["SMILES"], self.d_P4up[CASRN]["Chemical name"], self.d_P4up[CASRN]["Efficacy/potency"], self.d_P4up[CASRN]["LEC (µM)"]))

            if chemset == "all":
                filout.write("CASRN\tSMILES\tChemical name\tMC\tgenotox\tE2up\tP4up\tER\tGroup\tAff\n")
                for CASRN in self.d_all.keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(CASRN, self.d_all[CASRN]["SMILES"], self.d_all[CASRN]["name"], self.d_all[CASRN]["MC"], self.d_all[CASRN]["genotox"], self.d_all[CASRN]["E2up"], self.d_all[CASRN]["P4up"], "-".join(self.d_all[CASRN]["ER"]),  "-".join(self.d_all[CASRN]["Group"]),  "-".join(self.d_all[CASRN]["Group"])))
            filout.close()

        self.d_dataset = d_out

    def ToxPrintCount(self, l_chemsets):

        pr_out = pathFolder.createFolder(self.pr_out + "ToxPrintCount/")
        l_toxprints = list(self.d_toxprint[list(self.d_toxprint.keys())[0]].keys())
        l_toxprints.remove('INPUT')
        l_toxprints.remove('DTXSID')
        l_toxprints.remove('PREFERRED_NAME')

        for chemset in l_chemsets:
            d_dataset = toolbox.loadMatrix(self.d_dataset[chemset])
            d_temp = {}
            for toxprint in l_toxprints:
                d_temp[toxprint] = 0

            for chem in d_dataset.keys():
                if not chem in list(self.d_toxprint.keys()):
                    continue
                for toxprint in l_toxprints:
                    if self.d_toxprint[chem][toxprint] == "1":
                        d_temp[toxprint] = d_temp[toxprint] + 1

            p_filout = pr_out + chemset
            filout = open(p_filout, "w")
            filout.write("Toxprint\tcount\n")
            for toxprint in l_toxprints:
                filout.write("%s\t%s\n"%(toxprint, d_temp[toxprint]))
            filout.close()
            runExternal.barplotToxPrint(p_filout)

    def comparisonToxPrintCount(self, l_chemsets):
        pr_out = pathFolder.createFolder(self.pr_out + "ComparisonToxPrintCount/")
        l_toxprints = list(self.d_toxprint[list(self.d_toxprint.keys())[0]].keys())
        l_toxprints.remove('INPUT')
        l_toxprints.remove('DTXSID')
        l_toxprints.remove('PREFERRED_NAME')

        pr_countToxPrints = pathFolder.createFolder(self.pr_out + "ToxPrintCount/")

        d_out = {}
        for chemset in l_chemsets:
            d_count = toolbox.loadMatrix(pr_countToxPrints + chemset)
            d_out[chemset] = d_count

        p_filout = pr_out + "-".join(l_chemsets) + ""
        filout = open(p_filout, "w")
        filout.write("Toxprint\t"+"\t".join(l_chemsets) + "\n")
        for toxprint in l_toxprints:
            filout.write("%s\t%s\n"%(toxprint, "\t".join([str(d_out[chemset][toxprint]['count']) for chemset in l_chemsets])))
        filout.close()

        runExternal.plotX2(p_filout)

    def prepSets(self):
        self.loadCrossRefExcel()
        self.splitMCGenotox()
        self.defineE2P4active(["higher", "medium", "lower"])
        self.defineERactive()
        self.mergeAllSets()

        # write different chemical sets
        self.formatSetofChem(["ER", "MC", "Steroid", "E2", "P4", "all", "Steroid-up", "ER-agonist"])#["Steroid-up", "Steroid"])#, "ER-agonist", "Steroid"])
       
    def clusterMC(self):

        pr_out = pathFolder.createFolder(self.pr_out + "clusterMC/")
        p_all = self.d_dataset["all"]

        # extract descriptor
        c_chems = Chemicals.Chemicals(p_all, pr_out)
        c_chems.computeDesc(self.pr_desc)

        # run dendogram with circle of prop
        runExternal.dendogramProp(p_all, c_chems.p_desc_rdkit) 
        STOP402

        return

    def main(self):

        # prepare set of chemicals
        self.prepSets()
        
        STOPMC393

        # Compute descriptor
        pr_analysis = pathFolder.createFolder(self.pr_out + "desc_by_list/")
        for dataset in self.d_dataset.keys():
            pr_chems = pathFolder.createFolder(pr_analysis + dataset + "/")
            c_Chems = Chemicals.Chemicals(self.d_dataset[dataset], pr_chems)
            c_Chems.computeDesc(self.pr_desc)
            c_Chems.computeAllOperaPred(self.pr_desc)
            c_Chems.buildDescSet(["rdkit"])
            c_Chems.analysisDescBasedData(self.COR_VAL, self.MAX_QUANTILE, PCA=1, histDesc=1, SOM=1, clustering=0, Hclust=1)



        # overlap between dataset
        self.overlapBetweenListChem(["MC", "genotoxic", "Steroid-up", "ER-agonist"])
        self.overlapBetweenListChem(["MC", "Steroid", "ER"])
        self.overlapBetweenListChem(["MC", "Steroid", "ER-agonist", "genotoxic"])
        self.overlapBetweenListChem(["MC", "Steroid-up", "ER-agonist", "genotoxic"])
        # load ToxPrint
        self.loadToxPrint()

        # Write datasets
        self.formatSetofChem(["ER", "MC", "Steroid", "E2", "P4", "all", "Steroid-up", "ER-agonist"])#["Steroid-up", "Steroid"])#, "ER-agonist", "Steroid"])
        # plot toxplint barplot
        self.ToxPrintCount(["MC", "Steroid-up", "Steroid", "ER-agonist", "all"])

        # X2 for toxprint
        self.comparisonToxPrintCount(["MC", "Steroid-up", "Steroid"])

    def computeDescAnalysisFromList(p_list, p_ToxPrint, p_chemList, PR_RESULTS):

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

    def mergedataset(l_psetchems, PR_RESULTS):
        
        pr_desc = pathFolder.createFolder(PR_RESULTS + "DESC/")
        l_dataset = []
        d_desc_1D2D = {}
        d_desc_opera = {}
        for psetchem in l_psetchems:
            name_dataset =  psetchem.split("/")[-1][0:-4]
            pr_results = pathFolder.createFolder(PR_RESULTS + "analysis_individual-dataset/" + name_dataset + "/")       
            l_dataset.append(name_dataset)

            # Compute desc
            ###################
            cChem = Chemicals.Chemicals(psetchem, pr_results)
            cChem.computeDesc(pr_desc)

            d_desc_1D2D_chem = toolbox.loadMatrix(cChem.p_desc_rdkit)
            d_desc_opera_chem = toolbox.loadMatrix(cChem.p_desc_OPERA, sep = ",")
            

            for chem in d_desc_1D2D_chem.keys():
                if not chem in list(d_desc_1D2D.keys()):
                    d_desc_1D2D[chem] = deepcopy(d_desc_1D2D_chem[chem])
                    d_desc_1D2D[chem]["dataset"] = []
                d_desc_1D2D[chem]["dataset"].append(name_dataset)
            
            for chem in d_desc_opera_chem.keys():
                if not chem in list(d_desc_opera.keys()):
                    d_desc_opera[chem] = deepcopy(d_desc_opera_chem[chem])
                    d_desc_opera[chem]["dataset"] = []
                d_desc_opera[chem]["dataset"].append(name_dataset)
            
        pr_out = pathFolder.createFolder(PR_RESULTS + "-".join(l_dataset) + "/")
        p_desc2D = pr_out + "desc1D2D.csv"
        l_h = list(d_desc_1D2D[list(d_desc_1D2D.keys())[0]].keys()) 
        l_h.remove("CASRN")
        f_desc2D = open(p_desc2D, "w")
        f_desc2D.write("CASRN\t" + "\t".join(l_h) + "\n")
        for chem in d_desc_1D2D.keys():
            f_desc2D.write(chem)
            for h in l_h:
                if h == "dataset":
                    f_desc2D.write("\t%s"%("--".join(d_desc_1D2D[chem][h])))
                else:
                    f_desc2D.write("\t%s"%(d_desc_1D2D[chem][h]))
            f_desc2D.write("\n")
        f_desc2D.close()

        return [p_desc2D]

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



