import pathFolder
import ToxCastLib

p_ICE = "/mnt/d/database/invitroDB3.3_7-14-21/ICE_invitroDB3.3_7-14-21/cHTS2021_invitrodb33_20210128.txt"
p_assays_sum = "/mnt/d/database/invitroDB3.3_7-14-21/EPA_invitroDB-3.3_7-14-21/INVITRODB_V3_3_SUMMARY/assay_annotation_information_invitrodb_v3_3.xlsx"
p_gene_mapping = "/mnt/d/database/invitroDB3.3_7-14-21/EPA_invitroDB-3.3_7-14-21/INVITRODB_V3_3_SUMMARY/gene_target_information_invitrodb_v3_3.xlsx"



class ToxCast:
    def __init__(self, l_genes, pr_out):
    
        self.pr_out = pathFolder.createFolder(pr_out)
        self.l_genes = l_genes


    def loadAssaysMatrix(self):


        self.c_ToxCastLib = ToxCastLib.ToxCastLib(p_ICE, p_assays_sum, p_gene_mapping)
        self.p_assays = self.c_ToxCastLib.get_resultTableFromGenes(self.l_genes, self.pr_out)

        return self.p_assays