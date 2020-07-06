from os import path, rename, listdir
from random import shuffle
from re import search

import toolbox
import pathFolder
import runExternal
import searchInComptox

# import descriptor computation scripts => precise folder where descriptor are included
import sys
sys.path.insert(0, "../../../../ILS/development/molecular-descriptors/")
import Chemical



class dataset:
    def __init__(self, p_data, pr_out):
        self.p_data = p_data
        self.pr_out = pr_out
        self.p_SALTS = path.abspath("./Salts.txt")


    def loadDataset(self, loadDb =0):

        # define output
        d_out = toolbox.loadMatrix(self.p_data, sep = ",")
        self.d_dataset = d_out


        if loadDb == 1:
        # load SMILES from comptox
        l_header = list(d_dataset[list(d_dataset.keys())[0]].keys())
        if not "SMILES" in l_header:
            for chem in d_dataset.keys():
                CASRN = d_dataset[chem]["CASRN"]
                print("CASRN:", CASRN, "LOAD chem")
                c_search = searchInComptox.loadComptox(CASRN)
                c_search.searchInDB()
                if c_search.err != 1:
                    self.d_dataset[chem]["SMILES"] = c_search.SMILES
                else:
                    self.d_dataset[chem]["SMILES"] = "--"

                print("Load done:", self.d_dataset[chem]["SMILES"])


    def computeStructuralDesc(self):

        if not "pr_desc" in self.__dict__:
            pr_desc = pathFolder.createFolder(self.pr_out + "DESC/")
            self.pr_desc = pr_desc

        p_filout = self.pr_desc + "desc_1D2D.csv"
        if path.exists(p_filout):
            return p_filout

        # extract descriptor 2D
        l_desc = Chemical.getLdesc("1D2D")

        # open filout
        filout = open(p_filout, "w")
        filout.write("CASRN\tSMILES\t%s\n"%("\t".join(l_desc)))

        # compute descriptor
        for CASRN in self.d_dataset.keys():
            SMILES = self.d_dataset[CASRN]["SMILES"]
            print(SMILES)
            cChem = Chemical.Chemical(SMILES, self.pr_desc, p_salts=self.p_SALTS)
            cChem.prepChem() # prep
            # case error cleaning
            if cChem.err == 1:
                continue
            cChem.computeAll2D() # compute
            cChem.writeMatrix("2D") # write by chem to save time in case of rerun
            if cChem.err == 1:
                continue
            else:
                # write direcly descriptor
                filout.write("%s\t%s\t%s\n"%(CASRN, cChem.smi, "\t".join([str(cChem.all2D[desc]) for desc in l_desc])))

        filout.close()

        return p_filout


    def computeOPERADesc(self):

        if not "pr_desc" in self.__dict__:
            pr_desc = pathFolder.createFolder(self.pr_out + "DESC/")
            self.pr_desc = pr_desc

        p_filout = self.pr_desc + "desc_OPERA.csv"
        if path.exists(p_filout):
            return p_filout


        print("ERROR - OPERA desc no computed")
        return 0


    def computePNG(self):

        if not "pr_desc" in self.__dict__:
            pr_desc = pathFolder.createFolder(self.pr_out + "DESC/")
            self.pr_desc = pr_desc

        pr_png = pathFolder.createFolder(self.pr_desc + "PNG/")

        # compute descriptor
        l_CASRN = list(self.d_dataset.keys())
        shuffle(l_CASRN)
        for CASRN in l_CASRN:
            p_png = pr_png + CASRN + ".png"
            if path.exists(p_png):
                continue
            else:
                SMILES = self.d_dataset[CASRN]["SMILES"]
                cChem = Chemical.Chemical(SMILES, self.pr_desc, OS="Window", p_salts=self.p_SALTS)
                cChem.prepChem() # prep
                p_png_inch = cChem.computePNG()
                if cChem.err == 0:
                    rename(p_png_inch, p_png)
    



    def predictBiotransformation(self):

        pr_out = pathFolder.createFolder(self.pr_out + "bioTransformation/run/")

        l_CASRN = list(self.d_dataset.keys())
        shuffle(l_CASRN)

        for CASRN in l_CASRN:
            p_biotransformation_human = "%s%s_AllHuman.csv"%(pr_out, CASRN)
            p_biotransformation_env = "%s%s_Env.csv"%(pr_out, CASRN)

            if path.exists(p_biotransformation_env) and path.exists(p_biotransformation_human):
                continue 
            else:
                SMILES = self.d_dataset[CASRN]["SMILES"]
                cChem = Chemical.Chemical(SMILES, self.pr_desc, OS="Window", p_salts=self.p_SALTS)
                cChem.prepChem() 
                if cChem.err == 1:
                    continue
                
                # for all humain
                if not path.exists(p_biotransformation_human):
                    runExternal.BioTransformer(cChem.smi, "allHuman", p_biotransformation_human)

                # for env
                if not path.exists(p_biotransformation_env):
                    runExternal.BioTransformer(cChem.smi, "env", p_biotransformation_env)
        
        self.pr_biotransformation_run = pr_out



    def extractBioTranformationChemical(self, p_check = ""):

        pr_out = pathFolder.createFolder(self.pr_out + "bioTransformation/")
        self.pr_biotransformation = pr_out

        # 1 - Extract results from the run file
        p_biotransformed = pr_out + "sumBiotranformedChemical.csv"
        if not path.exists(p_biotransformed):
            d_biotransform = {}
            l_biotransformed = listdir(pr_out + "run/")
            for f_biotransformed in l_biotransformed:
                if search("sumBiotranformedChemical", f_biotransformed):
                    continue
                else:
                    orignCASRN = f_biotransformed.split("_")[0]
                    typeBiotransfo = f_biotransformed.split("_")[1]

                    if not orignCASRN in list(d_biotransform.keys()):
                        d_biotransform[orignCASRN] = {}

                    d_temp = toolbox.loadMatrix(self.pr_biotransformation + "run/" + f_biotransformed, sep = ",")
                    
                    for chem in d_temp.keys():
                        SMILES = d_temp[chem]["SMILES"]
                        if not SMILES in list(d_biotransform[orignCASRN].keys()):
                            d_biotransform[orignCASRN][SMILES] = {}
                            d_biotransform[orignCASRN][SMILES]["Biosystem"] = []
                            d_biotransform[orignCASRN][SMILES]["Reaction"] = []
                            d_biotransform[orignCASRN][SMILES]["Enzyme(s)"] = []
                            d_biotransform[orignCASRN][SMILES]["Precursor SMILES"] = []
                    
                        if not d_temp[chem]["Biosystem"] in d_biotransform[orignCASRN][SMILES]["Biosystem"]:
                            d_biotransform[orignCASRN][SMILES]["Biosystem"].append(d_temp[chem]["Biosystem"])

                        if not d_temp[chem]["Reaction"] in d_biotransform[orignCASRN][SMILES]["Reaction"]:
                            d_biotransform[orignCASRN][SMILES]["Reaction"].append(d_temp[chem]["Reaction"])

                        if not d_temp[chem]["Enzyme(s)"] in d_biotransform[orignCASRN][SMILES]["Enzyme(s)"]:
                            d_biotransform[orignCASRN][SMILES]["Enzyme(s)"].append(d_temp[chem]["Enzyme(s)"])

                        if not d_temp[chem]["Precursor SMILES"] in d_biotransform[orignCASRN][SMILES]["Precursor SMILES"]:
                            d_biotransform[orignCASRN][SMILES]["Precursor SMILES"].append(d_temp[chem]["Precursor SMILES"])
            


            # write summary
            filout = open(p_biotransformed, "w")
            filout.write("SMILES transformed\tBiosystem\tReaction\tEnzyme(s)\tCASRN origin\tSMILES origin\n")
            for CASRN_orgin in d_biotransform.keys():
                for SMILES in d_biotransform[CASRN_orgin].keys():
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(SMILES, ";".join(d_biotransform[CASRN_orgin][SMILES]["Biosystem"]), ";".join(d_biotransform[CASRN_orgin][SMILES]["Reaction"]), ";".join(d_biotransform[CASRN_orgin][SMILES]["Enzyme(s)"]), CASRN_orgin,";".join(d_biotransform[CASRN_orgin][SMILES]["Precursor SMILES"]) ))
            filout.close()
            self.p_biotransformed = p_biotransformed
            l_product_biotransformed = toolbox.loadMatrixToList(p_biotransformed) # reload in list
       
        else:
            self.p_biotransformed = p_biotransformed
            l_product_biotransformed = toolbox.loadMatrixToList(p_biotransformed)

        

        # 2 - clean and prepare biotransformation output
        p_products_cleaned = self.pr_biotransformation + "biotransformed_cleaned.csv"
        if path.exists(p_products_cleaned): 
            self.p_biotransformed_products_cleaned = p_products_cleaned
            self.d_product = toolbox.loadMatrix(p_products_cleaned)
        else:
        
            # load folder to check 
            if p_check != "":
                d_toCheck = toolbox.loadMatrix(p_check, sep = ",")
                # compute clean SMILES to be sure it is the same
                l_SMI_master = []
                for chem in d_toCheck.keys():
                    SMILES = d_toCheck[chem]["SMILES"]
                    cChem = Chemical.Chemical(SMILES, self.pr_biotransformation, OS="Window", p_salts=self.p_SALTS)
                    cChem.prepChem()
                    if cChem.err != 1:
                        d_toCheck[chem]["SMILES_clean"] = cChem.smi
                        if not cChem.smi in l_SMI_master:
                            l_SMI_master.append(cChem.smi)
            else:
                l_SMI_master = []
            
            d_SMILES_out = {}
            ID = 1
            for d_product in l_product_biotransformed:
                SMILES = d_product["SMILES transformed"]
                CASRN_orgin = d_product["CASRN origin"]
                if not SMILES in list(d_SMILES_out.keys()):  
                    cChem = Chemical.Chemical(SMILES, self.pr_biotransformation, OS="Window", p_salts=self.p_SALTS)
                    cChem.prepChem()
                    if cChem.err != 1:
                        cChem.generateInchiKey()
                        d_SMILES_out[SMILES] = {}
                        d_SMILES_out[SMILES]["SMILES"] = cChem.smi
                        d_SMILES_out[SMILES]["INCHIKEY"] = cChem.inchikey
                        d_SMILES_out[SMILES]["biotransform"]= []
                        d_SMILES_out[SMILES]["CASRN origin"]= []
                        d_SMILES_out[SMILES]["ID"] = ID
                        ID = ID + 1
                    else:
                        continue

                l_biosystems = d_product["Biosystem"].split(";")
                for biosystem in l_biosystems:
                    if not biosystem in d_SMILES_out[SMILES]["biotransform"]:
                            d_SMILES_out[SMILES]["biotransform"].append(biosystem)
                    if not CASRN_orgin in list(d_SMILES_out[SMILES]["CASRN origin"]):
                        d_SMILES_out[SMILES]["CASRN origin"].append(CASRN_orgin)
        

            filout = open(p_products_cleaned, "w")
            filout.write("ID\tSMILES\tIn master list\tINCHIKEY\tType biotransformation\tCASRN origin\n")

            for SMILES in d_SMILES_out.keys():
                inlist = 0
                if SMILES in l_SMI_master:
                    inlist = 1
                filout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%( d_SMILES_out[SMILES]["ID"], SMILES, inlist, d_SMILES_out[SMILES]["INCHIKEY"], ",".join(d_SMILES_out[SMILES]["biotransform"]), ",".join(d_SMILES_out[SMILES]["CASRN origin"]) ))     
            filout.close()

            self.p_biotransformed_products_cleaned = p_products_cleaned
            self.d_product = toolbox.loadMatrix(p_products_cleaned)



    def searchBiotransformedProduceInDB(self):

        if not "d_product" in self.__dict__:
            print("RUN BIOTRANSFORM FIRST")
            return 

        p_filout = self.pr_biotransformation + "biotransformed_cleaned_names.csv"
        
        for chem_ID in self.d_product.keys():
            inchikey = self.d_product[chem_ID]["INCHIKEY"]
            c_loadFromComptox = searchInComptox.loadComptox(inchikey)
            c_loadFromComptox.searchInDB()
            if c_loadFromComptox.err == 0:
                self.d_product[chem_ID]["name"] = c_loadFromComptox.nameChem
                self.d_product[chem_ID]["DTXSID"] = c_loadFromComptox.DTXSID
                self.d_product[chem_ID]["CASRN"] = c_loadFromComptox.CASRN
            else:
                self.d_product[chem_ID]["name"] = "--"
                self.d_product[chem_ID]["DTXSID"] = "--"
                self.d_product[chem_ID]["CASRN"] = "--"
        

        filout = open(p_filout, "w")
        filout.write("ID\tName\tCASRN\tSMILES\tDTXSID\tIn master list\tType biotransformation\tCASRN origin\n")

        for chem_ID in self.d_product.keys():
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(self.d_product[chem_ID]["ID"], self.d_product[chem_ID]["name"], self.d_product[chem_ID]["CASRN"], self.d_product[chem_ID]["SMILES"], self.d_product[chem_ID]["DTXSID"], self.d_product[chem_ID]["In master list"], self.d_product[chem_ID]["Type biotransformation"], self.d_product[chem_ID]["CASRN origin"]))     
        filout.close()



    def computeDescProductBiotransformed(self):

        if not "d_biotransformed_products" in self.__dict__:
            print("Compute biotransformation first and cleanning function")
            return 
        
        p_desc2D = self.pr_biotransformation + "desc2D.csv"

        # extract descriptor 2D
        l_desc = Chemical.getLdesc("1D2D")

        # open filout
        filout = open(p_desc2D, "w")
        filout.write("ID\tSMILES\t%s\n"%("\t".join(l_desc)))

        pr_desc = pathFolder.createFolder(self.pr_biotransformation + "DESC/")
        pr_PNG = pathFolder.createFolder(self.pr_biotransformation + "PNG/")

        for SMILES in self.d_biotransformed_products.keys():
            
            cChem = Chemical.Chemical(SMILES, pr_desc, p_salts = self.p_SALTS, OS="Windows")
            cChem.prepChem() # prep
            # case error cleaning
            if cChem.err == 1:
                continue
            cChem.computeAll2D() # compute
            cChem.writeMatrix("2D") # write by chem to save time in case of rerun
            
            # PNG
            p_inch_png = cChem.computePNG(prPNG=pr_PNG)
            p_png = pr_PNG + str(self.d_biotransformed_products[SMILES]["ID"]) + ".png"
            rename(p_inch_png, p_png)

            if cChem.err == 1:
                continue
            else:
                # write direcly descriptor
                filout.write("%s\t%s\t%s\n"%(self.d_biotransformed_products[SMILES]["ID"], cChem.smi, "\t".join([str(cChem.all2D[desc]) for desc in l_desc])))
        filout.close()