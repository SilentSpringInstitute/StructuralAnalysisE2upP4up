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

        # PUC
        self.d_PUC = toolbox.loadMatrix(self.p_PUC_dict, sep = ",")

        # document ID with PUC
        l_doc_PUC = toolbox.loadMatrixToList(self.p_product_compo, sep = ",")

        self.d_PUC_map_chemical = {}
        for d_doc_PUC in l_doc_PUC:
            PUC_id = d_doc_PUC["puc_id"]
            try: self.d_PUC_map_chemical[d_doc_PUC["chemical_id"]]
            except:
                self.d_PUC_map_chemical[d_doc_PUC["chemical_id"]] = []

            if not PUC_id in self.d_PUC_map_chemical[d_doc_PUC["chemical_id"]] :
                self.d_PUC_map_chemical[d_doc_PUC["chemical_id"]].append(PUC_id)

        # functionnal used        
        l_d_functionnal_used = toolbox.loadMatrixToList(self.p_fuctional_use_dict, sep = ",")
        d_functional_used_out = {}
        for d_functional_used in l_d_functionnal_used:
            chemical_id = d_functional_used["chemical_id"]
            try: d_functional_used_out[chemical_id]
            except:
                d_functional_used_out[chemical_id] = {}
                d_functional_used_out[chemical_id]["functional_use_id"] = []
                d_functional_used_out[chemical_id]["report_funcuse"] = []
                d_functional_used_out[chemical_id]["oecd_function"] = []
            d_functional_used_out[chemical_id]["functional_use_id"].append(d_functional_used["functional_use_id"])
            d_functional_used_out[chemical_id]["report_funcuse"].append(d_functional_used["report_funcuse"])
            d_functional_used_out[chemical_id]["oecd_function"].append(d_functional_used["oecd_function"])
        self.d_functional_used =  d_functional_used_out  

        # presence in dict
        self.d_list_presence = toolbox.loadMatrix(self.p_list_presence, sep = ",")
        
        l_presence_data = toolbox.loadMatrixToList(self.p_presence_data, sep = ",")
        self.d_presence_map = {}
        for d_presence_data in l_presence_data:
            chemical_id = d_presence_data["chemical_id"]
            document_id = d_presence_data["document_id"]
            presence_id = d_presence_data["list_presence_id"]
            try:self.d_presence_map[chemical_id]
            except: 
                self.d_presence_map[chemical_id] = {}
                self.d_presence_map[chemical_id]["l_document_id"] = []
                self.d_presence_map[chemical_id]["l_presence_id"] = []
            
            if not document_id in self.d_presence_map[chemical_id]["l_document_id"]:
                self.d_presence_map[chemical_id]["l_document_id"].append(document_id)
            if not presence_id in self.d_presence_map[chemical_id]["l_presence_id"]:
                self.d_presence_map[chemical_id]["l_presence_id"].append(presence_id)
        
        # load mapping from CASRN
        d_casrn = {}
        for chem in self.d_chem_mapping.keys():
            try:casrn = self.d_chem_mapping[chem]["preferred_casrn"]
            except: continue
            if casrn != "NA" and casrn != "":
                try:d_casrn[casrn].append(chem)
                except: d_casrn[casrn] = [chem]
        self.d_cas_mapping = d_casrn

        # load document ID
        self.d_document_by_id = toolbox.loadMatrix(self.p_doc_dict, sep = ",")
    
    def casrnToFunctions(self, casrn):
        """
        Use the chemical disctionnary but only pull together data from report used and oecd
        """
        # find the chem_id
        try:l_chem_id = self.d_cas_mapping[casrn]
        except: return {"l_chem_id":[], "funcuse":[], "oecd_function":[], "l_presence_id":[], "l_document_id":[]}

        d_out = {"l_chem_id":[],"funcuse":[], "oecd_function":[], "l_presence_id":[], "l_document_id":[]}
        for chem_id in l_chem_id:
            d_out["l_chem_id"].append(chem_id)
            try: l_funuse = self.d_functional_used[chem_id]["report_funcuse"]
            except: l_funuse = []
            try:l_oecd = self.d_functional_used[chem_id]["oecd_function"]
            except: l_oecd = []
            try: l_presence_id = self.d_presence_map[chem_id]["l_presence_id"]
            except: l_presence_id = []
            try: l_document_id = self.d_presence_map[chem_id]["l_document_id"]
            except: l_document_id = []

            d_out["l_presence_id"] = l_presence_id
            d_out["l_document_id"] = l_document_id

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
            filout.write("CASRN\tl_chem_id\tUse\toecd\tlist_presence_id\tlist_presence_def\tlist_document_id\tlist_document\n")
            for casn in d_out.keys():
                
                # search presence in list
                l_def_list = []
                for id_presence in d_out[casn]["l_presence_id"]:
                    if id_presence == "NA": ## need to check it is inside
                        continue
                    else:
                        l_def_list.append(self.d_list_presence[id_presence]["definition"])
                
                l_doc_list = []
                for id_doc in d_out[casn]["l_document_id"]:
                    l_doc_list.append(self.d_document_by_id[id_doc]["title"])

                filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(casn, ",".join(d_out[casn]["l_chem_id"]), ",".join(d_out[casn]["funcuse"]), ",".join(d_out[casn]["oecd_function"]),  ",".join(d_out[casn]["l_presence_id"]), ",".join(l_def_list), ",".join(d_out[casn]["l_document_id"]), ",".join(l_doc_list)))
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
            filout.write("CASRN\tchemical_id\tfunction_used\toecd\tpresence_name\tpresence_definition\tdocument_title\tclass_combine\n")

        for casrn in self.d_casrn_mapped.keys():
            d_out[casrn] = []
            if self.d_casrn_mapped[casrn]["l_chem_id"] == []:
                if p_filout != "":
                    filout.write("%s\tNA\tNA\tNA\tNA\tNA\tNA\tNo data\n"%(casrn))
                continue
            
            # define as a string
            self.d_casrn_mapped[casrn]["funcuse"] = ",".join(self.d_casrn_mapped[casrn]["funcuse"])
            self.d_casrn_mapped[casrn]["oecd_function"] = ",".join(self.d_casrn_mapped[casrn]["oecd_function"])
            
            l_funct_from_funcuse = self.searchBoardExposureInString(self.d_casrn_mapped[casrn]["funcuse"])
            l_funct_from_oecd = self.searchBoardExposureInString(self.d_casrn_mapped[casrn]["oecd_function"])

            l_presence = []
            l_def_presence = []
            l_doc_title = []
            l_doc_exp = []
            l_cas_exp = self.individualMappingOnBoardExp(casrn)
            l_PUC_exp = []


            #PUC -- when PUC add in Consumer products
            for chemical_id in self.d_casrn_mapped[casrn]["l_chem_id"]:
                try:
                    l_PUC = self.d_PUC_map_chemical[chemical_id]
                    l_PUC_exp.append("Consumer products")
                except:
                    pass

            # presence ID
            for presence_id in self.d_casrn_mapped[casrn]["l_presence_id"]:
                if presence_id == "NA":
                    continue
                l_presence.append(self.d_list_presence[presence_id]["name"])
                l_def_presence.append(self.d_list_presence[presence_id]["definition"])

            # document
            for document_id in self.d_casrn_mapped[casrn]["l_document_id"]:
                if document_id == "NA":
                    continue
                l_doc_title.append(self.d_document_by_id[document_id]["title"])
                l_doc_exp = l_doc_exp + self.mapDocumentToBoardExp(document_id)

            # remove duplication
            l_presence = list(set(l_presence))
            l_def_presence = list(set(l_def_presence))
            l_doc_title = list(set(l_def_presence))

            l_funct_from_presence_name = self.searchBoardExposureInString(",".join(l_presence))
            l_funct_from_presence_def = self.searchBoardExposureInString(",".join(l_def_presence))
            l_funct_from_doc = self.searchBoardExposureInString(",".join(l_doc_title))

            # add in the d-out dictionnary
            d_out[casrn] = d_out[casrn] + l_funct_from_funcuse + l_funct_from_oecd + l_funct_from_presence_name + l_funct_from_presence_def + l_funct_from_doc + l_doc_exp + l_cas_exp + l_PUC_exp
            d_out[casrn] = list(set(d_out[casrn]))

            if d_out[casrn] == []:
                d_out[casrn] = ["No data"]

            if p_filout != "":
                filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(casrn, "-".join(self.d_casrn_mapped[casrn]["l_chem_id"]), "-".join(l_funct_from_funcuse), "-".join(l_funct_from_oecd), "-".join(l_funct_from_presence_name), "-".join(l_funct_from_presence_def), "-".join(l_funct_from_doc), "-".join(d_out[casrn])))


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
        if search("insecticide", str_in):
            l_funct.append("Pesticides")


        if search("industrial", str_in): 
            l_funct.append("Industrial")
        if search("NACE", str_in):
            l_funct.append("Industrial")
        if search("coal tar", str_in):
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
        if search("uv stabilizer", str_in):
            l_funct.append("Industrial")        
        if search("flame retardant", str_in):
            l_funct.append("Industrial")
        if search("colorants", str_in):
            l_funct.append("Industrial") 

        if search("food", str_in) and not search("not for food", str_in) and not search("Nonfood", str_in):
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
        if search("pharma", str_in):
            l_funct.append("Pharmaceuticals")
        if search("sunscreen agent", str_in):
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
        if search("shampoo", str_in):
            l_funct.append("Consumer products") 
        if search("cosmetic", str_in):
            l_funct.append("Consumer products")
        if search("perfuming", str_in):
            l_funct.append("Consumer products")
        if search("perfume", str_in):
            l_funct.append("Consumer products")
        if search("flame retardant", str_in):
            l_funct.append("Consumer products")
        if search("personal care", str_in):
            l_funct.append("Consumer products")
        if search("skin conditioning", str_in):
            l_funct.append("Consumer products")
        if search("sunscreen agent", str_in):
            l_funct.append("Consumer products")
        if search("coal tar", str_in):
            l_funct.append("Consumer products") 
        if search("colorants", str_in):
            l_funct.append("Consumer products") 


        return list(set(l_funct))

    def mapDocumentToBoardExp(self, document_ID):

        ##Air Water INC Fine Chemicals List
        if document_ID == "1373515":
            return ["Industrial"]

        ##Fl@vis Flavouring Substances
        if document_ID == "1513117":
            return ["Diet"]

        ##U.S. FDA list of Indirect Additives used in Food Contact Substances; presence of a substance in this list indicates that only certain intended uses and use conditions are authorized by FDA regulations (version of list updated 10/4/2018)
        if document_ID == "460":
            return ["Diet"]

        ##Indirect Additives used in Food Contact Substances
        if document_ID == "1372213":
            return ["Pesticides"]

        ##Inert Ingredients Permitted for Use in Nonfood Use Pesticide Products
        if document_ID == "1365244":
            return ["Pesticides"]  

        ##2007 Pesticide Residues in Blueberries
        if document_ID == "1374900":
            return ["Pesticides"] 


        ##Harmonized Tariff Schedule of the United States (2019) Intermediate Chemicals for Dyes Appendix
        if document_ID == "1371498":
            return ["Consumer products"]

        ##present on the WA State Department of Ecology - Chemicals of High Concern to Children Reporting List (version of list pulled 4/24/2020),pertaining to  or intended for use specifically by children
        if document_ID == "453478":
            return ["Consumer products", "Industrial"]
        
        ##chemicals measured or identified in environmental media or products,Sources specific to a European country or subset of countries,writing utensils containing liquid or gel ink
        if document_ID == "400407471":
            return ["Consumer products"]

        #substances on the International Fragrance Association's ordered register of all fragrance ingredients used in consumer goods by the fragrance industry's customers worldwide
        if document_ID == "519":
            return ["Consumer products", "Industrial"]

        ##applied to all data sources used in MN DoH chemical screening proof of concept,chemicals measured or identified in environmental media or products,water intended for drinking  or related to drinking water; includes bottled water  finished water from drinking water treatment plants  and untreated water that has been denoted as a drinking source
        if document_ID == "392400422":
            return ["Pesticides", "Diet"]

        ## Relating to pesticides or pesticide usage. Includes specific types of pesticides  e.g. insecticides   herbicides  fungicides  and fumigants; also includes general biocides,chemical residues  typically from drugs or pesticides
        if document_ID == "423446":
            return ["Pesticides", "Pharmaceuticals"]       

        ##chemicals measured or identified in environmental media or products,Relating to pesticides or pesticide usage. Includes specific types of pesticides  e.g. insecticides   herbicides  fungicides  and fumigants; also includes general biocides,general term pertaining to agricultural practices  including the raising and farming of animals and growing of crops,includes fresh  canned and frozen forms  as well as juices and sauces (e.g. applesauce)  excludes forms intended for consumption by young children (i.e. baby foods); includes green beans and peas,chemical residues  typically from drugs or pesticides
        if document_ID == "400423425442446":
            return ["Pesticides", "Pharmaceuticals"] 

        ##Pharmaceutical Appendix (2019) Table 1
        if document_ID == "1372195":
            return ["Pharmaceuticals"]

        ##Pharmaceutical Appendix (2019) Table 3
        if document_ID == "1372197":
            return ["Pharmaceuticals"]        

        return []
    
    def individualMappingOnBoardExp(self, casrn):
        
        if casrn == "87818-31-3":
            return ["Pesticides"]

        if casrn == "4291-63-8":
            return ["Pharmaceuticals"]

        return []

    def comparisonBoardExposureWithCardona2021(self, p_board_exp, p_csv_cardona2021, pr_out):

        l_d_cardona = toolbox.loadMatrixToList(p_csv_cardona2021, sep = ",")
        d_board_exp = toolbox.loadMatrix(p_board_exp)

        d_out = {}
        for casrn in d_board_exp.keys():
            d_out[casrn] = {}
            d_out[casrn]["cardona2021"] = []
            d_out[casrn]["board_exp"] = []
            d_out[casrn]["aggred"]  = "0"


            exposure = d_board_exp[casrn]["class_combine"]
            d_out[casrn]["board_exp"] = exposure.split("-")
            d_out[casrn]["board_exp"] = list(set(d_out[casrn]["board_exp"]))

            for d_cardona in l_d_cardona:
                casrn_cardona = d_cardona["CASN_protect"].replace(" ", "")
                if casrn == casrn_cardona:
                    if d_cardona["Consumer"] == "1" or d_cardona["Consumer"] == "1*":
                        d_out[casrn]["cardona2021"].append("Consumer products")
                    if d_cardona["Diet"] == "1" or d_cardona["Diet"] == "1*":
                        d_out[casrn]["cardona2021"].append("Diet")
                    if d_cardona["Industrial"] == "1" or d_cardona["Industrial"] == "1*":
                        d_out[casrn]["cardona2021"].append("Industrial")
                    if d_cardona["Pest."] == "1" or d_cardona["Pest."] == "1*":
                        d_out[casrn]["cardona2021"].append("Pesticides")
                    if d_cardona["Pharma."] == "1" or d_cardona["Pharma."] == "1*":
                        d_out[casrn]["cardona2021"].append("Pharmaceuticals")
                    if d_cardona["No exposure source data"] == "1" or d_cardona["No exposure source data"] == "1*":
                        d_out[casrn]["cardona2021"].append("No data")
                    
                    d_out[casrn]["cardona2021"] = list(set(d_out[casrn]["cardona2021"]))
                    if d_out[casrn]["cardona2021"] == d_out[casrn]["board_exp"]:
                        d_out[casrn]["aggred"] = "1"
            
        # write in file
        p_filout = pr_out + "overlapCardonaCPDATMap.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\tCardonaList\tMappingList\tAggree\n")

        for casrn in d_out.keys():
            filout.write("%s\t%s\t%s\t%s\n"%(casrn, "-".join(d_out[casrn]["cardona2021"]), "-".join(d_out[casrn]["board_exp"]), d_out[casrn]["aggred"]))
        filout.close()


## testing
pr_database = "/mnt/d/database/CPDat/CPDatRelease20201216/"
c_CPDAT = CPDat(pr_database)
c_CPDAT.comparisonBoardExposureWithCardona2021("/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/ClassCPDAT_E2up/board_exposure.csv", "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/data/E2up_CPDAT_Cardona2021.csv", "/mnt/c/Users/AlexandreBorrel/research/SSI/e2up_p4up/results/ClassCPDAT_E2up/")
#c_CPDAT.loadMapping()
#out = c_CPDAT.casrnToFunctions("7235-40-7")
#print(out)