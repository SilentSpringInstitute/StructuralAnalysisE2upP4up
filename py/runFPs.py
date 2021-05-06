from os import path
import toolbox

import FPToolbox


class runFPs:

    def __init__(self, d_dataset, d_ToxPrint, pr_out):

        self.d_dataset = d_dataset
        self.d_ToxPrint = d_ToxPrint
        self.pr_out = pr_out



    def computeTanimotoMatrix(self):

        all_computed = 1
        d_out ={}
        for dataset in self.d_dataset.keys():
            p_matrix_out = "%s%s.csv"%(self.pr_out, dataset)
            d_out[dataset] = p_matrix_out
            if not path.exists (p_matrix_out):
                all_computed = 0
        
        self.d_FPMatrix = d_out
        if all_computed == 1:
            return


        d_Tanimoto = {}
        l_casrn = list(self.d_ToxPrint.keys())

        l_toxprints = list(self.d_ToxPrint[l_casrn[0]].keys())
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
            c_FP_i = FPToolbox.FingerPrint(self.d_ToxPrint[l_casrn[i]])
            c_FP_i.prepFP()

            j = i
            while j < imax:
                c_FP_j = FPToolbox.FingerPrint(self.d_ToxPrint[l_casrn[j]])
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

            p_matrix_out = "%s%s.csv"%(self.pr_out, dataset)

            #write matrix
            filout = open(p_matrix_out, "w")
            filout.write("\t" + "\t".join(l_casrn_dataset) + "\n")
            for casrn in l_casrn_dataset:
                filout.write("%s\t%s\n"%(casrn, "\t".join([str(d_Tanimoto[casrn][c]) if c in list(d_Tanimoto[casrn].keys()) else str(d_Tanimoto[c][casrn]) for c in l_casrn_dataset])))
            filout.close()

            


