import toolbox
import PubMed
import pandas as pd


class pubmedSearch:
    def __init__(self, p_tab_chem, l_term_search, pr_out, email):
        
        self.pr_out = pr_out
        self.l_terms = l_term_search
        self.p_tab_chem = p_tab_chem
        self.email = email
        

    def load_tab_chem(self):
        
        d_chem = toolbox.loadExcelSheet(self.p_tab_chem, name_sheet='ToxPrint_QSAR', k_head = "CASRN")
        self.d_chem = d_chem

    def do_search_bychem(self):
        p_sum = self.pr_out + "chem.sum"
        f_sum = open(p_sum, "w")
        f_sum.write("CASRN\tQuery\tNb Result\tlog\n")
        for chem in list(self.d_chem.keys()):
            if pd.isna(self.d_chem[chem]["PubMed search chemical"]) == False:
                c_query = PubMed.searching(self.email, f_sum)
                c_query.buildNewQuery(chem, self.d_chem[chem]["PubMed search chemical"])
                c_query.search()
        f_sum.close()

    def do_search_combine_term(self):
        for term in self.l_terms:
            p_sum = "%schem_%s.sum"%(self.pr_out, term)
            f_sum = open(p_sum, "w")
            f_sum.write("\"CASRN\"\t\"Query\"\t\"Nb Result\"\t\"Article\"\n")    

            for chem in list(self.d_chem.keys()):
                if pd.isna(self.d_chem[chem]["PubMed search chemical"]) == False:
                    c_query = PubMed.searching(self.email, f_sum)
                    c_query.buildNewQuery(chem, self.d_chem[chem]["PubMed search chemical"])
                    c_query.buildingQuery(term, "AND")
                    c_query.search(w=True)
        f_sum.close()


    def do_search(self, l_combine_search, review = "on", date = ""):
        self.load_tab_chem()
        #self.do_search_bychem()
        self.do_search_combine_term()        
        


