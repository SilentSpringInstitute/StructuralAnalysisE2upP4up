from os import path, listdir
from re import search

import FPToolbox
import toolbox
import pathFolder
import runExternal

class runToxPrint:

    def __init__(self, d_dataset, pr_toxprint, pr_out):

        self.d_dataset = d_dataset
        self.pr_toxprint = pr_toxprint
        self.pr_out = pr_out

    def loadToxPrint(self):
        """Load from inputs/comptox the toxprint by chemicals
        """
        d_toxprint = {}
        l_p_ToxPrint = listdir(self.pr_toxprint)
        for p_ToxPrint in l_p_ToxPrint:
            if search("ToxPrint", p_ToxPrint):
                d_temp = toolbox.loadMatrix(self.pr_toxprint + p_ToxPrint, sep=",")
                d_toxprint.update(d_temp)
        
        # clean 
        l_chem = list(d_toxprint.keys())
        i = 0
        imax = len(l_chem)
        while i < imax:
            if d_toxprint[l_chem[i]]["atom:element_main_group"] == "-":
                del d_toxprint[l_chem[i]]
            i = i + 1

        self.d_toxprint = d_toxprint

    def computeTanimotoMatrix(self):

        pr_out = pathFolder.createFolder(self.pr_out + "ToxPrintTanimoto_by_list/")

        all_computed = 1
        d_out ={}
        for dataset in self.d_dataset.keys():
            p_matrix_out = "%s%s.csv"%(pr_out, dataset)
            d_out[dataset] = p_matrix_out
            if not path.exists (p_matrix_out):
                all_computed = 0
        
        self.d_FPMatrix = d_out
        if all_computed == 1:
            return

        d_Tanimoto = {}
        l_casrn = list(self.d_toxprint.keys())

        l_toxprints = list(self.d_toxprint[l_casrn[0]].keys())
        if "INPUT" in l_toxprints:
            l_toxprints.remove("INPUT")
        if "DTXSID" in l_toxprints:
            l_toxprints.remove("DTXSID")
        if "PREFERRED_NAME" in l_toxprints:
            l_toxprints.remove("PREFERRED_NAME")

        # compute all matrix of tanimoto
        
        i = 0
        imax = len(l_casrn)
        while i < imax:
            d_Tanimoto[l_casrn[i]] = {}
            c_FP_i = FPToolbox.FingerPrint(self.d_toxprint[l_casrn[i]])
            c_FP_i.prepFP()

            j = i
            while j < imax:
                c_FP_j = FPToolbox.FingerPrint(self.d_toxprint[l_casrn[j]])
                c_FP_j.prepFP()

                # compute score
                scoreT = c_FP_i.jaccardScore(c_FP_j.bits, exclude0=0)
                #scoreT = c_FP_i.DICEScore(c_FP_j.bits)
                d_Tanimoto[l_casrn[i]][l_casrn[j]] = scoreT
                j = j + 1
            i = i + 1 
                

        for dataset in self.d_dataset.keys():
            p_dataset = self.d_dataset[dataset]
            d_dataset = toolbox.loadMatrix(p_dataset)

            l_casrn_dataset = list(set(list(d_Tanimoto.keys())) & set(list(d_dataset.keys())))

            p_matrix_out = "%s%s.csv"%(pr_out, dataset)

            #write matrix
            filout = open(p_matrix_out, "w")
            filout.write("\t" + "\t".join(l_casrn_dataset) + "\n")
            for casrn in l_casrn_dataset:
                filout.write("%s\t%s\n"%(casrn, "\t".join([str(d_Tanimoto[casrn][c]) if c in list(d_Tanimoto[casrn].keys()) else str(d_Tanimoto[c][casrn]) for c in l_casrn_dataset])))
            filout.close()

    def ToxPrintCount(self):
        
        pr_out = pathFolder.createFolder(self.pr_out + "ToxPrintCount_by_list/")

        all_computed = 1
        d_out ={}
        for dataset in self.d_dataset.keys():
            p_count_out = "%s%s.csv"%(pr_out, dataset)
            d_out[dataset] = p_count_out
            if not path.exists (p_count_out):
                all_computed = 0
        
        
        if all_computed == 1:
            self.d_pcount = d_out
            return


        # compute all of the count
        d_out = {}
        l_toxprints = list(self.d_toxprint[list(self.d_toxprint.keys())[0]].keys())
        l_toxprints.remove('INPUT')
        l_toxprints.remove('DTXSID')
        l_toxprints.remove('PREFERRED_NAME')

        for dataset in self.d_dataset.keys():
            p_dataset = self.d_dataset[dataset]
            d_dataset = toolbox.loadMatrix(p_dataset)

            # initialize the count
            d_temp = {}
            for toxprint in l_toxprints:
                d_temp[toxprint] = 0

            for chem in d_dataset.keys():
                if not chem in list(self.d_toxprint.keys()):
                    continue
                for toxprint in l_toxprints:
                    if self.d_toxprint[chem][toxprint] == "1":
                        d_temp[toxprint] = d_temp[toxprint] + 1
            
            p_count_out = "%s%s.csv"%(pr_out, dataset)
            filout = open(p_count_out, "w")
            filout.write("Toxprint\tcount\n")
            for toxprint in l_toxprints:
                filout.write("%s\t%s\n"%(toxprint, d_temp[toxprint]))
            filout.close()
            
            d_out[dataset] = p_count_out

            # run barplot
            runExternal.barplotToxPrint(p_count_out)
        
        self.d_pcount = d_out
           
    def comparisonToxPrintCount(self, l_chemsets):
        
        if not "d_pcount" in self.__dict__:
            self.ToxPrintCount()
        
        pr_out = pathFolder.createFolder("%sToxPrintComparisonCount/%s"%(self.pr_out, "-".join(l_chemsets))) + "/"

        # define list of toxprint
        l_toxprints = list(self.d_toxprint[list(self.d_toxprint.keys())[0]].keys())
        l_toxprints.remove('INPUT')
        l_toxprints.remove('DTXSID')
        l_toxprints.remove('PREFERRED_NAME')

        d_out = {}
        for chemset in l_chemsets:
            d_count = toolbox.loadMatrix(self.d_pcount[chemset])
            d_out[chemset] = d_count

        p_filout = pr_out + "count.csv"
        filout = open(p_filout, "w")
        filout.write("Toxprint\t"+"\t".join(l_chemsets) + "\n")
        for toxprint in l_toxprints:
            filout.write("%s\t%s\n"%(toxprint, "\t".join([str(d_out[chemset][toxprint]['count']) for chemset in l_chemsets])))
        filout.close()

        runExternal.plotX2(p_filout)


    def writeToxPrintMatrix(self, l_casrn, p_filout):

        filout = open(p_filout, "w")
        
        # define list of toxprint
        l_toxprints = list(self.d_toxprint[list(self.d_toxprint.keys())[0]].keys())
        l_toxprints.remove('INPUT')
        l_toxprints.remove('DTXSID')
        l_toxprints.remove('PREFERRED_NAME')

        filout.write("CASRN\t" + "\t".join(l_toxprints) + "\n")

        for casrn in l_casrn:
            try:filout.write("%s\t%s\n"%(casrn, "\t".join([str(self.d_toxprint[casrn][col]) for col in l_toxprints])))
            except: pass
        filout.close()

        return p_filout

