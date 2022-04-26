import pathFolder
import toolbox
import runExternal

class comparisonChemicalLists:
    def __init__(self, l_p_list, pr_out):
        self.l_p_list = l_p_list
        self.pr_out = pathFolder.createFolder(pr_out + "__".join([p_list.split("/")[-1][0:-4] for p_list in l_p_list]) + "/")
        self.pr_results = pr_out

    def X2Comparison(self):

        p_forX2 = pathFolder.createFolder(self.pr_out + "forX2Comparison/") + "forX2"
        f_forX2 = open(p_forX2, "w")

        # header 
        f_forX2.write("\t%s\n"%("\t".join([p_list.split("/")[-1][0:-4] for p_list in self.l_p_list])))
        
        # open ToxPrint
        d_toxprint = {}
        l_toxprint = []
        for p_list in self.l_p_list:
            k_in = p_list.split("/")[-1][0:-4]
            p_toxprint = self.pr_results + "analysis_individual-dataset/"+ k_in + "/ToxPrint/count_toxprint"
            print(p_toxprint)
            d_toxprint[k_in] = toolbox.loadMatrix(p_toxprint)
            l_toxprint = l_toxprint + list(d_toxprint[k_in].keys())
        
        l_toxprint = list(set(l_toxprint))

        for toxprint in l_toxprint:
            f_forX2.write("%s"%(toxprint))
            for p_list in self.l_p_list:
                k_in = p_list.split("/")[-1][0:-4]
                try:
                    f_forX2.write("\t%s"%(d_toxprint[k_in][toxprint]["count"]))
                except:
                    f_forX2.write("\t0")
            f_forX2.write("\n")
        f_forX2.close()

        runExternal.plotX2(p_forX2)

    def comparisonDescTwoLists(self):

        


        return 