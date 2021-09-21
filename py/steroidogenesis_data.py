from os import path
import toolbox
import pathFolder
import runExternal

class Steroidogenesis_data:
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
        
        # define hormones for the analysis
        self.l_hormones = ["OHPREG", "PROG", "OHPROG", "DOC", "CORTISOL", "11DCORT", "ANDR", "TESTO", "ESTRONE", "ESTRADIOL"]

    def loadsinglehitc(self):

        l_lines_hitc = toolbox.loadMatrixToList(self.pf_hormone_hitc, sep = ",")
        d_hitc = {}
        if not "d_chem" in self.__dict__:
            self.loadChemicalsMapping()
        
        for line_hitc in l_lines_hitc:
            chid = line_hitc["chid"]
            if chid == "NA":
                continue
            aenm = line_hitc["aenm"]
            hormone = aenm.split("_")[2]
            up_down = aenm.split("_")[-1]
            if not chid in list(d_hitc.keys()):
                d_hitc[chid] = {}
            if line_hitc["hitc"] != "0":
                if up_down == "up":
                    d_hitc[chid][aenm] = float(line_hitc["max_med"])
                else:
                    d_hitc[chid][aenm] = -float(line_hitc["max_med"])
        
        self.d_hitc = d_hitc

    def summaryHitc(self):

        p_filout = self.pr_results + "single_hitc_matrix.csv"
        if path.exists(p_filout):
            d_out = toolbox.loadMatrix(p_filout)
            self.d_single_hit = d_out
            return 
        
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

        d_out = toolbox.loadMatrix(p_filout)
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

        p_filout = self.pr_results + "CR_hitc_matrix.csv"
        if path.exists(p_filout):
            self.d_CR_hit = toolbox.loadMatrix(p_filout)
            return
        
        if not "d_chem_mapping" in self.__dict__:
            self.loadChemicalsMapping

        d_out = {}
        
        l_lines_CR = toolbox.loadMatrixToList(self.pf_hormone_change, sep = ',')
        for line_CR in l_lines_CR:
            chid = line_CR["chid"]
            casn = self.d_chem_mapping[chid]["casn"]
            d_out[casn] = {}
            for h in self.l_hormones:
                d_out[casn][h] = line_CR[h]
        
        
        filout = open(p_filout, "w")
        filout.write("CASRN\t%s\tSUM_HORMONE_CHANGED\n"%("\t".join(self.l_hormones)))
        for casrn in d_out.keys():
            nb_h_changed = 0
            for h in self.l_hormones: 
                if d_out[casrn][h] != "0": 
                    nb_h_changed = nb_h_changed +1 
            filout.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(d_out[casrn][h]) for h in self.l_hormones]), nb_h_changed))
        filout.close()
        runExternal.barplotHormones(p_filout)
        self.d_CR_hit = toolbox.loadMatrix(p_filout)


    def corHormoneEndpoint(self):

        pr_out = pathFolder.createFolder(self.pr_results + "CorHorm/")
        p_filout = pr_out + "horm_resp.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\t%s\n"%("\t".join(self.l_hormones)))


        # plot correlation for 2 hormones
        for chem in self.d_CR_hit.keys():
            filout.write("%s\t%s\n"%(chem, "\t".join([self.d_CR_hit[chem][h] for h in self.l_hormones])))
        filout.close()     

        # plot correlation
        runExternal.corHormResponse(p_filout, pr_out)


    def main(self):
        """
        Load data from stereogenesis 
        """

        self.loadChemicalsMapping()
        self.loadsinglehitc()
        self.summaryHitc()
        self.loadChemicalCR()

        #self.corHormoneEndpoint()

