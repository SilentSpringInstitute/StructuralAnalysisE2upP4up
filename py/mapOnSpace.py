from os import path

import toolbox
import runExternal

# import descriptor computation scripts => precise folder where descriptor are included
import sys
sys.path.insert(0, "../../../../ILS/development/molecular-descriptors/")
import Chemical


class mapOnSpace:
    def __init__(self, p_chem, p_space_1D2D, p_space_3D, pr_out):
        self.p_chem = p_chem
        self.p_space_1D2D = p_space_1D2D
        self.p_space_3D = p_space_3D
        self.pr_out = pr_out
        
        self.p_SALTS = path.abspath("./Salts.txt")



    def map(self, name_map):

        p_chem_prep = self.pr_out + "list_chem.csv"
        if not path.exists(p_chem_prep):
            # load list of chem
            f_list_chem = open(p_chem_prep, "w")
            f_list_chem.write("inchikey\tCASRN\tSMILES\n")
            d_chem = toolbox.loadMatrix(self.p_chem, sep = ",")
            for chem in d_chem.keys():
                SMILES = d_chem[chem]["SMILES"]
                cchem = Chemical.Chemical(SMILES, "", self.p_SALTS)
                inchkey = cchem.generateInchiKey()
                f_list_chem.write("%s\t%s\t%s\n"%(inchkey, chem, SMILES))
            f_list_chem.close()
        
        runExternal.projectToSpace(p_chem_prep, self.p_space_1D2D, self.p_space_3D, name_map, self.pr_out)




