from os import path

import pathFolder
import toolbox
import runExternal

class merge_MCcrossWithStereo:

    def __init__(self, c_MC, c_stereo, pr_out):
        self.pr_out = pr_out
        self.c_MC = c_MC
        self.c_stereo = c_stereo

    def analyzeByDataset(self, dataset, PCA=0, hclust=0):

        self.pr_analysis = pathFolder.createFolder("%sanalysis_%s/"%(self.pr_out, dataset))

        if PCA == 1:
            self.PCA_FoldChangeMC(dataset)

        if hclust == 1:
            self.addStereoOncluster(dataset)


    def main(self):

        # MC
        self.analyzeByDataset("MC", PCA=0, hclust=1)

        # E2
        #self.analyzeByDataset("E2", PCA=0, hclust=1)

        # P4
        #self.analyzeByDataset("P4", PCA=0, hclust=1)



    
    def PCA_FoldChangeMC(self, dataset):

        pr_out = pathFolder.createFolder(self.pr_analysis + "PCA_FoldChange/")
        d_dataset = toolbox.loadMatrix(self.c_MC.d_dataset[dataset])        

        ## single dose
        p_matrix_single_hitc = pr_out + "single_hit.csv"
        f_matrix_single_hitc = open(p_matrix_single_hitc, "w")
        f_matrix_single_hitc.write("CASRN\t%s\tinset\n"%("\t".join(self.c_stereo.l_hormones)))
        for casrn in self.c_stereo.d_single_hit.keys():
            if casrn in list(d_dataset.keys()):
                inset = 1
            else:
                inset = 0
            
            f_matrix_single_hitc.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(self.c_stereo.d_single_hit[casrn][h]) for h in self.c_stereo.l_hormones]), inset))
        f_matrix_single_hitc.close()
        runExternal.PCA_SteroiMC(p_matrix_single_hitc)


        ## CR response
        p_matrix_CR_hitc = pr_out + "CR_hit.csv"
        f_matrix_CR_hitc = open(p_matrix_CR_hitc, "w")
        f_matrix_CR_hitc.write("CASRN\t%s\tinset\n"%("\t".join(self.c_stereo.l_hormones)))
        for casrn in self.c_stereo.d_CR_hit.keys():
            if casrn in list(d_dataset.keys()):
                inset = 1
            else:
                inset = 0
            
            f_matrix_CR_hitc.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(self.c_stereo.d_CR_hit[casrn][h]) for h in self.c_stereo.l_hormones]), inset))
        f_matrix_CR_hitc.close()
        runExternal.PCA_SteroiMC(p_matrix_CR_hitc)

    def addStereoOncluster(self, dataset):

        pr_out = pathFolder.createFolder(self.pr_analysis + "hclust_FoldChange/")

    

        # CR hormone
        p_dendo_CR = pr_out + "dendogram_CR.png"
        if not path.exists(p_dendo_CR):
            # write CR
            p_stereo_CR = pr_out + "stereo_CR.csv"
            f_stereo_CR = open(p_stereo_CR, "w")
            f_stereo_CR.write("CASRN\t%s\n"%("\t".join(self.c_stereo.l_hormones)))
            for casrn in self.c_stereo.d_CR_hit.keys():
                 f_stereo_CR.write("%s\t%s\n"%(casrn, "\t".join([str(self.c_stereo.d_CR_hit[casrn][h]) for h in self.c_stereo.l_hormones])))
            f_stereo_CR.close()

            runExternal.dendogramClusterTwoProp(self.c_MC.d_dataset[dataset], p_stereo_CR, self.c_MC.c_Desc.d_desc[dataset]["rdkit"], pr_out, self.c_MC.COR_VAL, self.c_MC.MAX_QUANTILE)


        # single point hormone
        p_dendo_single = pr_out + "dendogram_single.png"
        if not path.exists(p_dendo_single):
            # write single
            p_stereo_single = pr_out + "stereo_single.csv"
            f_stereo_single = open(p_stereo_single, "w")
            f_stereo_single.write("CASRN\t%s\n"%("\t".join(self.c_stereo.l_hormones)))
            for casrn in self.c_stereo.d_single_hit.keys():
                 f_stereo_single.write("%s\t%s\n"%(casrn, "\t".join([str(self.c_stereo.d_single_hit[casrn][h]) for h in self.c_stereo.l_hormones])))
            f_stereo_single.close()

            runExternal.dendogramClusterTwoProp(self.c_MC.d_dataset[dataset], p_stereo_single, self.c_MC.c_Desc.d_desc[dataset]["rdkit"], pr_out, self.c_MC.COR_VAL, self.c_MC.MAX_QUANTILE)




    def cardTanimoto(self, d_MC):

        pr_out = pathFolder.createFolder(self.pr_results + "card_FP/")
        
        ## single dose
        p_matrix_single_hitc = pr_out + "single_hit_matrix.csv"
        f_matrix_single_hitc = open(p_matrix_single_hitc, "w")
        f_matrix_single_hitc.write("CASRN\t%s\tMC\n"%("\t".join(self.l_hormones)))
        for chid in self.d_single_hit.keys():
            casrn = self.d_chem_mapping[chid]["casn"]
            if casrn in list(d_MC.keys()):
                MC = 1
            else:
                MC = 0
            f_matrix_single_hitc.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(self.d_single_hit[chid][h]) for h in self.l_hormones]), MC))
        f_matrix_single_hitc.close()
        runExternal.FP_card(p_matrix_single_hitc)

        ## CR response
        p_matrix_CR_hitc = pr_out + "CR_hit_matrix.csv"
        f_matrix_CR_hitc = open(p_matrix_CR_hitc, "w")
        f_matrix_CR_hitc.write("CASRN\t%s\tMC\n"%("\t".join(self.l_hormones)))
        for casrn in self.d_CR_hit.keys():
            if casrn in list(d_MC.keys()):
                MC = 1
            else:
                MC = 0
            
            f_matrix_CR_hitc.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(self.d_CR_hit[casrn][h]) for h in self.l_hormones]), MC))
        f_matrix_CR_hitc.close()
        runExternal.FP_card(p_matrix_CR_hitc)