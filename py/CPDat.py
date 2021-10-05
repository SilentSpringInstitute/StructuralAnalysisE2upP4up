from numpy.lib.npyio import load
import toolbox
from re import search


class CPDat:
    def __init__(self, pr_database="/mnt/d/database/CPDat/CPDatRelease20201216/", version="20201216"):

        self.pr_database = pr_database
        self.version = version

        # load dataset
        self.p_chem_dict = pr_database + "chemical_dictionary_" + self.version + ".csv"
        self.p_presence_data = pr_database + "list_presence_data_" + self.version + ".csv"
        self.p_HHE_data = pr_database + "HHE_data_" + self.version + ".csv"
        self.p_doc_dict = pr_database + "document_dictionary_" + self.version + ".csv"
        self.p_list_presence = pr_database + "list_presence_dictionary_" + self.version + ".csv"
        self.p_PUC_dict = pr_database + "PUC_dictionary_" + self.version + ".csv"
        self.p_funct_used = pr_database + "functional_use_data_" + self.version + ".csv"
        self.p_product_compo = pr_database + "product_composition_data_" + self.version + ".csv"
        self.p_QSUR_data = pr_database + "QSUR_data_" + self.version + ".csv"
        self.p_fuctional_use_dict = pr_database + "functional_use_dictionary_" + self.version + ".csv" 

    def loadMapping(self):

        self.d_chem_mapping = toolbox.loadMatrix(self.p_chem_dict, sep = ",")
        l_d_functionnal_used = toolbox.loadMatrixToList(self.p_fuctional_use_dict, sep = ",")
        d_functional_used_out = {}
        for d_functional_used in l_d_functionnal_used:
            chemical_id = d_functional_used["chemical_id"]
            if not chemical_id in list(d_functional_used_out.keys()):
                d_functional_used_out[chemical_id] = {}
                d_functional_used_out[chemical_id]["functional_use_id"] = []
                d_functional_used_out[chemical_id]["report_funcuse"] = []
                d_functional_used_out[chemical_id]["oecd_function"] = []
            d_functional_used_out[chemical_id]["functional_use_id"].append(d_functional_used["functional_use_id"])
            d_functional_used_out[chemical_id]["report_funcuse"].append(d_functional_used["report_funcuse"])
            d_functional_used_out[chemical_id]["oecd_function"].append(d_functional_used["oecd_function"])
        self.d_functional_used =  d_functional_used_out  

        # presence
        self.d_list_presence = toolbox.loadMatrix(self.p_list_presence, sep = ",")
        l_presence_data = toolbox.loadMatrixToList(self.p_presence_data, sep = ",")

        self.d_presence_map = {}
        for d_presence_data in l_presence_data:
            chemical_id = d_presence_data["chemical_id"]
            document_id = d_presence_data["document_id"]
            presence_id = d_presence_data["list_presence_id"]
            if not chemical_id in list(self.d_presence_map.keys()):
                self.d_presence_map[chemical_id] = {}
                self.d_presence_map[chemical_id]["document_id"] = []
                self.d_presence_map[chemical_id]["list_presence_id"] = []
            if not document_id in self.d_presence_map[chemical_id]["document_id"]:
                self.d_presence_map[chemical_id]["document_id"].append(document_id)
            if not presence_id in d_presence_data["list_presence_id"]:
                self.d_presence_map[chemical_id]["list_presence_id"].append(presence_id)
        
        # load mapping from CASRN
        d_casrn = {}
        for chem in self.d_chem_mapping.keys():
            try:casrn = self.d_chem_mapping[chem]["preferred_casrn"]
            except: continue
            if casrn != "NA" and casrn != "":
                try:d_casrn[casrn].append(chem)
                except: d_casrn[casrn] = [chem]

        self.d_cas_mapping = d_casrn

    
    def casrnToFunctions(self, casrn):
        """
        Use the chemical disctionnary but only pull together data from report used and oecd
        """
        # find the chem_id
        try:l_chem_id = self.d_cas_mapping[casrn]
        except: return {"l_chem_id":[], "funcuse":[], "oecd_function":[]}

        d_out = {"l_chem_id":[],"funcuse":[], "oecd_function":[]}
        for chem_id in l_chem_id:
            d_out["l_chem_id"].append(chem_id)
            try:
                l_funuse = self.d_functional_used[chem_id]["report_funcuse"]
                l_oecd = self.d_functional_used[chem_id]["oecd_function"]
            except: continue

            for funuse in l_funuse:
                if not funuse.lower() in d_out["funcuse"]:
                    d_out["funcuse"].append(funuse.lower())

            for oecd in l_oecd:
                if not oecd.lower() in d_out["oecd_function"]:
                    d_out["oecd_function"].append(oecd.lower())

        return d_out


    def listCasToFunct(self, l_casrn, p_filout = ""):
        """
        Build a dictionnary with the CASRN as key and value from the cpdat
        """
        d_out = {}

        for CASRN in l_casrn:
            d_exposure = self.casrnToFunctions(CASRN)
            d_out[CASRN] = d_exposure

        if p_filout != "":
            filout = open(p_filout, "w")
            filout.write("CASRN\tUse\toecd\tlist_def\n")
            for casn in d_out.keys():
                
                # search presence in list
                l_def_list = []
                for chemical_id in d_out[casn]["l_chem_id"]:
                    if chemical_id in list(self.d_presence_map.keys()):
                        l_id_presence = self.d_presence_map[chemical_id]["list_presence_id"]
                        for id_presence in l_id_presence:
                            l_def_list.append(self.d_list_presence[id_presence]["definition"])
                filout.write("%s\t%s\t%s\n"%(casn, ",".join(d_out[casn]["funcuse"]), ",".join(d_out[casn]["oecd_function"]), ",".join(l_def_list)))
            filout.close()
        
        self.d_casrn_mapped = d_out
        return d_out

    def extractBoardExposure(self, p_filout = ""):
        """
        From Cardona 2021 script
        List of function [Pesticides, Industrial, Consumer products, Diet, Pharmaceuticals, No data]
        """
        
        d_out = {}
        if not "d_casrn_mapped" in self.__dict__:
            print("== Load CPDAT first with the list of CASRN ==")
            return 
            
        if p_filout != "":
            filout = open(p_filout, "w")
            filout.write("CASRN\tchemical_id\tfunction_used\toecd\tpresence_name\tpresence_definition\n")

        for casrn in self.d_casrn_mapped.keys():
            if self.d_casrn_mapped[casrn]["l_chem_id"] == []:
                if p_filout != "":
                    filout.write("%s\tNA\tNA\tNA\tNA\tNA\n"%(casrn))
                continue
            
            # define as a string
            self.d_casrn_mapped[casrn]["funcuse"] = ",".join(self.d_casrn_mapped[casrn]["funcuse"])
            self.d_casrn_mapped[casrn]["oecd_function"] = ",".join(self.d_casrn_mapped[casrn]["oecd_function"])

            l_funct_from_funcuse = self.searchBoardExposureInString(self.d_casrn_mapped[casrn]["funcuse"])
            l_funct_from_oecd = self.searchBoardExposureInString(self.d_casrn_mapped[casrn]["oecd_function"])

            l_presence = []
            l_def_presence = []
            
            for chemical_id in self.d_casrn_mapped[casrn]["l_chem_id"]:
                try:l_presence_id = self.d_presence_map[chemical_id]["list_presence_id"]
                except:continue

                for presence_id in l_presence_id:
                    l_presence.append(self.d_list_presence[presence_id]["name"])
                    l_def_presence.append(self.d_list_presence[presence_id]["definition"])

            # remove duplication
            l_presence = list(set(l_presence))
            l_def_presence = list(set(l_def_presence))

            l_funct_from_presence_name = self.searchBoardExposureInString(",".join(l_presence))
            l_funct_from_presence_def = self.searchBoardExposureInString(",".join(l_def_presence))

            if p_filout != "":
                filout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(casrn, "-".join(self.d_casrn_mapped[casrn]["l_chem_id"]), "-".join(l_funct_from_funcuse), "-".join(l_funct_from_oecd), "-".join(l_funct_from_presence_name), "-".join(l_funct_from_presence_def)))

        if p_filout != "":
            filout.close()


        self.d_board_exp = d_out
        return d_out
                               
    def searchBoardExposureInString(self, str_in):
        l_funct = []
        str_in = str_in.lower()
        if search("pesticide", str_in):
            l_funct.append("Pesticides")
        if search("antimicrobial", str_in):
            l_funct.append("Pesticides")
        if search("fungicide", str_in):
            l_funct.append("Pesticides")
        if search("extermination", str_in): 
            l_funct.append("Pesticides")
        if search("herbicide", str_in):
            l_funct.append("Pesticides")
        if search("herbicdisinfectantide", str_in):
            l_funct.append("Pesticides")

        if search("industrial", str_in): 
            l_funct.append("Industrial")
        if search("NACE", str_in):
            l_funct.append("Industrial")
            
        if search("raw_material", str_in):
            l_funct.append("Industrial")
        if search("industrial_manufacturing", str_in):
            l_funct.append("Industrial")
        if search("industrial_fluid", str_in):
            l_funct.append("Industrial")
        if search("mining", str_in):
            l_funct.append("Industrial")
        if search("manufacturing", str_in):
            l_funct.append("Industrial")
        if search("resource_extraction", str_in):
            l_funct.append("Industrial")
        if search("rubber_processing", str_in):
            l_funct.append("Industrial")
        if search("plasticizer", str_in):
            l_funct.append("Industrial")
        if search("plasticisers", str_in):
            l_funct.append("Industrial")
        if search("catalyst", str_in):
            l_funct.append("Industrial")

        if search("food", str_in):
            l_funct.append("Diet")
        if search("beverage", str_in):
            l_funct.append("Diet")
        if search("Retail", str_in):
            l_funct.append("Diet")
        if search("drinking", str_in):
            l_funct.append("Diet")
        if search("preservative", str_in):
            l_funct.append("Diet")
        if search("flavouring", str_in):
            l_funct.append("Diet")
            
        if search("drug", str_in):
            l_funct.append("Pharmaceuticals")

        if search("apparel", str_in):
            l_funct.append("Consumer products")
        if search("personal_care", str_in):
            l_funct.append("Consumer products")
        if search("arts_crafts", str_in):
            l_funct.append("Consumer products")
        if search("furniture", str_in): 
            l_funct.append("Consumer products")
        if search("child_use", str_in):
            l_funct.append("Consumer products")
        if search("decor", str_in):
            l_funct.append("Consumer products")
        if search("toy", str_in):
            l_funct.append("Consumer products")
        if search("electronics", str_in):
            l_funct.append("Consumer products")
        if search("lawn_garden", str_in):
            l_funct.append("Consumer products")
        if search("sports_equipment", str_in):
            l_funct.append("Consumer products")
        if search("baby_use", str_in):
            l_funct.append("Consumer products")
        if search("pet", str_in):
            l_funct.append("Consumer products")
        if search("dogs", str_in):
            l_funct.append("Consumer products")
        if search("cats", str_in):
            l_funct.append("Consumer products")

        if search("tools", str_in):
            l_funct.append("Consumer products")
        if search("dental", str_in):
            l_funct.append("Consumer products")
        if search("toothbrush", str_in):
            l_funct.append("Consumer products")
        if search("cleaning_washing", str_in):
            l_funct.append("Consumer products")
        if search("soap", str_in):
            l_funct.append("Consumer products")
        if search("automotive_care", str_in):
            l_funct.append("Consumer products")
        if search("hair dyeing", str_in):
            l_funct.append("Consumer products")
        if search("skin-care", str_in):
            l_funct.append("Consumer products")
        if search("hair conditioning", str_in):
            l_funct.append("Consumer products")

        return list(set(l_funct))

## testing
#pr_database = "/mnt/d/database/CPDat/CPDatRelease20201216/"
#c_CPDAT = CPDat(pr_database)
#c_CPDAT.loadMapping()
#out = c_CPDAT.casrnToFunctions("7235-40-7")
#print(out)