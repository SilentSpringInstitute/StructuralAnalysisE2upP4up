import toolbox
from re import search


class CPDat:
    def __init__(self, pr_database="/mnt/d/database/CPDat/CPDatRelease20201216/", version="20201216"):

        self.pr_database = pr_database
        self.version = version

        # load dataset
        self.p_chem_dict = pr_database + "chemical_dictionary_" + self.version + ".csv"
        self.p_l_presence = pr_database + "list_presence_data_" + self.version + ".csv"
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
        self.d_functional_used = toolbox.loadMatrix(self.p_fuctional_use_dict, sep = ",")

        print(len(list(self.d_chem_mapping.keys())))

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

        # find the chem_id
        try:l_chem_id = self.d_cas_mapping[casrn]
        except: return {"funcuse":[], "oecd_function":[]}

        d_out = {"funcuse":[], "oecd_function":[]}
        for chem_id in l_chem_id:
            try:
                funuse = self.d_functional_used[chem_id]["report_funcuse"]
                oecd = self.d_functional_used[chem_id]["oecd_function"]
            except: continue

            if not funuse.lower() in d_out["funcuse"]:
                d_out["funcuse"].append(funuse.lower())

            if not oecd in d_out["oecd_function"]:
                d_out["oecd_function"].append(oecd)

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
            filout.write("CASRN\tUse\toecd\n")
            for casn in d_out.keys():
                filout.write("%s\t%s\t%s\n"%(casn, ",".join(d_out[casn]["funcuse"]), ",".join(d_out[casn]["oecd_function"])))
            filout.close()
        
        self.d_casrn_mapped = d_out
        return d_out

    def extractBoardExposure(self):
        """
        From Cardona 2021 script
        List of function [Pesticides, Industrial, Consumer products, Diet, Pharmaceuticals, No data]

        ====> need to add a list mapping

        """
        
        d_out = {}
        if not "d_casrn_mapped" in self.__dict__:
            print("== Load CPDAT first with the list of CASRN ==")
            
        for casrn in self.d_casrn_mapped.keys():
            
            self.d_casrn_mapped[casrn]["funcuse"] = ",".join(self.d_casrn_mapped[casrn]["funcuse"])
            self.d_casrn_mapped[casrn]["oecd_function"] = ",".join(self.d_casrn_mapped[casrn]["oecd_function"])

            l_funct = []
            if search("pesticide", self.d_casrn_mapped[casrn]["funcuse"]) or search("pesticide", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Pesticides")
            if search("antimicrobial", self.d_casrn_mapped[casrn]["funcuse"]) or search("antimicrobial", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Pesticides")
            if search("fungicide", self.d_casrn_mapped[casrn]["funcuse"]) or search("fungicide", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Pesticides")
            if search("extermination", self.d_casrn_mapped[casrn]["funcuse"]) or search("extermination", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Pesticides")
            if search("herbicide", self.d_casrn_mapped[casrn]["funcuse"]) or search("herbicide", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Pesticides")
            if search("herbicdisinfectantide", self.d_casrn_mapped[casrn]["funcuse"]) or search("disinfectant", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Pesticides")

            if search("industrial", self.d_casrn_mapped[casrn]["funcuse"]) or search("industrial", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("NACE", self.d_casrn_mapped[casrn]["funcuse"]) or search("NACE", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("raw_material", self.d_casrn_mapped[casrn]["funcuse"]) or search("raw_material", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("industrial_manufacturing", self.d_casrn_mapped[casrn]["funcuse"]) or search("industrial_manufacturing", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("industrial_fluid", self.d_casrn_mapped[casrn]["funcuse"]) or search("industrial_fluid", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("mining", self.d_casrn_mapped[casrn]["funcuse"]) or search("mining", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("manufacturing", self.d_casrn_mapped[casrn]["funcuse"]) or search("manufacturing", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("resource_extraction", self.d_casrn_mapped[casrn]["funcuse"]) or search("resource_extraction", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("rubber_processing", self.d_casrn_mapped[casrn]["funcuse"]) or search("rubber_processing", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("plasticizer", self.d_casrn_mapped[casrn]["funcuse"]) or search("plasticizer", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("plasticisers", self.d_casrn_mapped[casrn]["funcuse"]) or search("plasticisers", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")
            if search("catalyst", self.d_casrn_mapped[casrn]["funcuse"]) or search("catalyst", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Industrial")

            if search("food", self.d_casrn_mapped[casrn]["funcuse"]) or search("food", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Diet")
            if search("beverage", self.d_casrn_mapped[casrn]["funcuse"]) or search("beverage", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Diet")
            if search("Retail", self.d_casrn_mapped[casrn]["funcuse"]) or search("Retail", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Diet")
            if search("drinking", self.d_casrn_mapped[casrn]["funcuse"]) or search("drinking", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Diet")
            if search("preservative", self.d_casrn_mapped[casrn]["funcuse"]) or search("preservative", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Diet")
            if search("flavouring", self.d_casrn_mapped[casrn]["funcuse"]) or search("flavouring", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Diet")
            
            if search("drug", self.d_casrn_mapped[casrn]["funcuse"]) or search("drug", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Pharmaceuticals")

            if search("apparel", self.d_casrn_mapped[casrn]["funcuse"]) or search("apparel", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("personal_care", self.d_casrn_mapped[casrn]["funcuse"]) or search("personal_care", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("arts_crafts", self.d_casrn_mapped[casrn]["funcuse"]) or search("arts_crafts", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("furniture", self.d_casrn_mapped[casrn]["funcuse"]) or search("furniture", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("child_use", self.d_casrn_mapped[casrn]["funcuse"]) or search("child_use", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("decor", self.d_casrn_mapped[casrn]["funcuse"]) or search("decor", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("toy", self.d_casrn_mapped[casrn]["funcuse"]) or search("toy", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("electronics", self.d_casrn_mapped[casrn]["funcuse"]) or search("electronics", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("lawn_garden", self.d_casrn_mapped[casrn]["funcuse"]) or search("lawn_garden", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("sports_equipment", self.d_casrn_mapped[casrn]["funcuse"]) or search("sports_equipment", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("baby_use", self.d_casrn_mapped[casrn]["funcuse"]) or search("baby_use", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("pet", self.d_casrn_mapped[casrn]["funcuse"]) or search("pet", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("dogs", self.d_casrn_mapped[casrn]["funcuse"]) or search("dogs", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("cats", self.d_casrn_mapped[casrn]["funcuse"]) or search("cats", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")

            if search("tools", self.d_casrn_mapped[casrn]["funcuse"]) or search("tools", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("dental", self.d_casrn_mapped[casrn]["funcuse"]) or search("dental", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("toothbrush", self.d_casrn_mapped[casrn]["funcuse"]) or search("toothbrush", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("cleaning_washing", self.d_casrn_mapped[casrn]["funcuse"]) or search("cleaning_washing", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("soap", self.d_casrn_mapped[casrn]["funcuse"]) or search("soap", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("automotive_care", self.d_casrn_mapped[casrn]["funcuse"]) or search("automotive_care", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("hair dyeing", self.d_casrn_mapped[casrn]["funcuse"]) or search("hair dyeing", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("skin-care", self.d_casrn_mapped[casrn]["funcuse"]) or search("skin-care", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")
            if search("hair conditioning", self.d_casrn_mapped[casrn]["funcuse"]) or search("hair conditioning", self.d_casrn_mapped[casrn]["oecd_function"]):
                l_funct.append("Consumer products")

            d_out[casrn] = list(set(l_funct))

        
        self.d_board_exp = d_out
        return d_out
                               


## testing
#pr_database = "/mnt/d/database/CPDat/CPDatRelease20201216/"
#c_CPDAT = CPDat(pr_database)
#c_CPDAT.loadMapping()
#out = c_CPDAT.casrnToFunctions("7235-40-7")
#print(out)