import FPToolbox


class runFPs:

    def __init__(self, d_dataset, d_ToxPrint, pr_out):

        self.d_dataset = d_dataset
        self.d_ToxPrint = d_ToxPrint
        self.pr_out = pr_out



    def computeTanimotoMatrix(self):

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
        imax = 10
        while i < imax:
            d_Tanimoto[l_casrn[i]] = {}
            c_FP_i = FPToolbox.FingerPrint(self.d_ToxPrint[l_casrn[i]])
            c_FP_i.prepFP()

            j = i
            while j < imax:
                c_FP_j = FPToolbox.FingerPrint(self.d_ToxPrint[l_casrn[j]])
                c_FP_j.prepFP()

                # compute score
                scoreT = c_FP_i.jaccardIndex(c_FP_j.bits)
                d_Tanimoto[l_casrn[i]][l_casrn[j]] = scoreT
                j = j + 1
            i = i + 1 
                
        print(d_Tanimoto)
