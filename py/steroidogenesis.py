from os import path
import toolbox
import pathFolder
import runExternal

class Steroidogenesis:
    """
    Class use to load and parse table from Kaumaus et al 2006, study develop from H295R
    """

    def __init__(self, pr_data, pr_results):

        pr_SI = pr_data + "Karmaus2016_SI/kfw002_Supplementary_Data/"
        pr_results = pathFolder.createFolder(pr_results + "steroidogenesis/")
        self.pr_results = pr_results

        self.pf_chem_MTT = pr_SI + "toxsci-15-0570-File006.csv"
        self.pf_concentration_response = pr_SI + "toxsci-15-0570-File006.csv"
        self.pf_hormone_resp = pr_SI + "toxsci-15-0570-File008.csv"
        self.pf_hormone_hitc = pr_SI + "toxsci-15-0570-File009.csv"
        self.pf_hormone_CR = pr_SI + "toxsci-15-0570-File010.csv"
        self.pf_raw_hormone_data = pr_SI + "toxsci-15-0570-File011.csv"
        self.pf_hormone_change = pr_SI + "toxsci-15-0570-File014.csv"
        
    def loadsinglehitc(self):

        l_lines_hitc = toolbox.loadMatrixToList(self.pf_hormone_hitc, sep = ",")
        d_hitc = {}
        l_hormones = []
        if not "d_chem" in self.__dict__:
            self.loadChemicalsMapping()
        
        for line_hitc in l_lines_hitc:
            chid = line_hitc["chid"]
            if chid == "NA":
                continue
            aenm = line_hitc["aenm"]
            hormone = aenm.split("_")[2]
            up_down = aenm.split("_")[-1]
            if not hormone in l_hormones:
                l_hormones.append(hormone)
            if not chid in list(d_hitc.keys()):
                d_hitc[chid] = {}
            if line_hitc["hitc"] != "0":
                if up_down == "up":
                    d_hitc[chid][aenm] = float(line_hitc["max_med"])
                else:
                    d_hitc[chid][aenm] = -float(line_hitc["max_med"])
        
        self.d_hitc = d_hitc
        self.l_hormones = l_hormones

    def summaryHitc(self):

        p_filout = self.pr_results + "single_hitc_matrix.csv"
        #if path.exists(p_filout):
        #    return 
        
        if not "d_hitc" in self.__dict__:
            self.loadsinglehitc()
            
        d_out = {}
        for chem in self.d_hitc.keys():
            if not chem in list(d_out.keys()):
                d_out[chem] = {}
                for h in self.l_hormones:
                    d_out[chem][h] = 0
            
            for aenm in self.d_hitc[chem].keys():
                hormone = aenm.split("_")[2]
                hitc = self.d_hitc[chem][aenm]
                if hitc != 0:
                    d_out[chem][hormone] = hitc
        
        filout = open(p_filout, "w")
        filout.write("CASRN\t%s\tSUM_HORMONE_CHANGED\n"%("\t".join(self.l_hormones)))
        for chem in d_out.keys():
            casrn = self.d_chem_mapping[chem]["casn"]
            nb_h_changed = 0
            for h in self.l_hormones: 
                if d_out[chem][h] != 0: 
                    nb_h_changed = nb_h_changed +1 
            filout.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(d_out[chem][h]) for h in self.l_hormones]), nb_h_changed))
        filout.close()
    
        # write a summary
        runExternal.barplotHormones(p_filout)

        self.d_single_hit = d_out

    def loadChemicalsMapping(self):

        d_chem = {}
        l_lines = toolbox.loadMatrixToList(self.pf_hormone_resp, sep = ",")

        for line_file in l_lines:
            chid = line_file["chid"]
            spid = line_file["spid"]
            casn = line_file["casn"]
            chnm = line_file["chnm"]

            if not chid in list(d_chem.keys()) and chid != "NA":
                d_chem[chid] = {}
                d_chem[chid]["chid"] = chid
                d_chem[chid]["spid"] = spid
                d_chem[chid]["casn"] = casn
                d_chem[chid]["chnm"] = chnm
        
        self.d_chem_mapping = d_chem

    def loadChemicalCR(self):

        if not "d_chem_mapping" in self.__dict__:
            self.loadChemicalsMapping

        if not "l_hormones" in self.__dict__:
            self.loadsinglehitc

        d_out = {}
        
        l_lines_CR = toolbox.loadMatrixToList(self.pf_hormone_change, sep = ',')
        for line_CR in l_lines_CR:
            chid = line_CR["chid"]
            casn = self.d_chem_mapping[chid]["casn"]
            d_out[casn] = {}
            for h in self.l_hormones:
                d_out[casn][h] = line_CR[h]
        
        self.d_CR_hit = d_out

        p_filout = self.pr_results + "CR_hitc_matrix.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\t%s\n"%("\t".join(self.l_hormones)))
        for casn in d_out.keys():
            filout.write("%s\t%s\n"%(casn, "\t".join([d_out[casn][h] for h in self.l_hormones])))
        
        filout.close()

    def main(self):

        self.loadChemicalsMapping()
        self.loadsinglehitc()
        self.summaryHitc()
        self.loadChemicalCR()

    def PCA_FoldChangeMC(self, d_MC):

        pr_out = pathFolder.createFolder(self.pr_results + "PCA_FoldChange_MC/")
        
        ## single dose
        p_matrix_single_hitc = pr_out + "single_hit_matrix.csv"
        f_matrix_single_hitc = open(p_matrix_single_hitc, "w")
        f_matrix_single_hitc.write("CASRN\t%s\tMC\n"%("\t".join(self.l_hormones)))
        for chid in self.d_single_hit.keys():
            casrn = self.d_chem_mapping[chid]["casn"]
            if casrn in list(d_MC.keys()):
                MC = 1
            else:
                MC = 0
            
            f_matrix_single_hitc.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(self.d_single_hit[chid][h]) for h in self.l_hormones]), MC))
        f_matrix_single_hitc.close()
        runExternal.PCA_SteroiMC(p_matrix_single_hitc)


        ## CR response
        p_matrix_CR_hitc = pr_out + "CR_hit_matrix.csv"
        f_matrix_CR_hitc = open(p_matrix_CR_hitc, "w")
        f_matrix_CR_hitc.write("CASRN\t%s\tMC\n"%("\t".join(self.l_hormones)))
        for casrn in self.d_CR_hit.keys():
            if casrn in list(d_MC.keys()):
                MC = 1
            else:
                MC = 0
            
            f_matrix_CR_hitc.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(self.d_CR_hit[casrn][h]) for h in self.l_hormones]), MC))
        f_matrix_CR_hitc.close()
        runExternal.PCA_SteroiMC(p_matrix_CR_hitc)

    def cardTanimoto(self, d_MC):

        pr_out = pathFolder.createFolder(self.pr_results + "card_FP/")
        
        ## single dose
        p_matrix_single_hitc = pr_out + "single_hit_matrix.csv"
        f_matrix_single_hitc = open(p_matrix_single_hitc, "w")
        f_matrix_single_hitc.write("CASRN\t%s\tMC\n"%("\t".join(self.l_hormones)))
        for chid in self.d_single_hit.keys():
            casrn = self.d_chem_mapping[chid]["casn"]
            if casrn in list(d_MC.keys()):
                MC = 1
            else:
                MC = 0
            f_matrix_single_hitc.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(self.d_single_hit[chid][h]) for h in self.l_hormones]), MC))
        f_matrix_single_hitc.close()
        runExternal.FP_card(p_matrix_single_hitc)

        ## CR response
        p_matrix_CR_hitc = pr_out + "CR_hit_matrix.csv"
        f_matrix_CR_hitc = open(p_matrix_CR_hitc, "w")
        f_matrix_CR_hitc.write("CASRN\t%s\tMC\n"%("\t".join(self.l_hormones)))
        for casrn in self.d_CR_hit.keys():
            if casrn in list(d_MC.keys()):
                MC = 1
            else:
                MC = 0
            
            f_matrix_CR_hitc.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(self.d_CR_hit[casrn][h]) for h in self.l_hormones]), MC))
        f_matrix_CR_hitc.close()
        runExternal.FP_card(p_matrix_CR_hitc)