from os import path, rename
from random import shuffle

import toolbox
import pathFolder

# import descriptor computation scripts => precise folder where descriptor are included
import sys
sys.path.insert(0, "../../../../ILS/development/molecular-descriptors/")
import Chemical



class dataset:
    def __init__(self, p_data, pr_out):
        self.p_data = p_data
        self.pr_out = pr_out
    

    def loadDataset(self):

        # define output
        d_out = toolbox.loadMatrix(self.p_data, sep = ",")
        self.d_dataset = d_out


    def computeDesc(self):

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
            cChem = Chemical.Chemical(SMILES, self.pr_desc, p_salts=path.abspath("./Salts.txt"))
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
                cChem = Chemical.Chemical(SMILES, self.pr_desc, OS="Window", p_salts=path.abspath("./Salts.txt"))
                cChem.prepChem() # prep
                p_png_inch = cChem.computePNG()
                if cChem.err == 0:
                    rename(p_png_inch, p_png)