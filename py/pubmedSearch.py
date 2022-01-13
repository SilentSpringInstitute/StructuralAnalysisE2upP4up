import toolbox
import PubMed
import pandas as pd
import xlsxwriter

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
        
        p_xlx = "%schem_search.xlsx"%(self.pr_out)
        workbook = xlsxwriter.Workbook(p_xlx)
        worksheet = workbook.add_worksheet()
        l_headers = ["CASRN", "Query", "Nb Result", "Article"]
        l_chem = list(self.d_chem.keys())
        # write header in the xlxs            
        [worksheet.write(0, i_col, l_headers[i_col]) for i_col in range(0, 4)]
            
        row = 1
        imax = len(list(self.d_chem.keys()))
        i = 0
        while i < imax:
            if pd.isna(self.d_chem[l_chem[i]]["PubMed search chemical"]) == True:
                i = i + 1
                continue
                
            c_query = PubMed.searching(self.email, "")
            c_query.buildNewQuery(l_chem[i], self.d_chem[l_chem[i]]["PubMed search chemical"])
            c_query.search(w=False)
                
            # add casrn
            worksheet.write(row, 0, l_chem[i])
            # add query 
            worksheet.write(row, 1, c_query.Query)
            # add nb results
            worksheet.write(row, 2, c_query.count)
            # add list of article
            worksheet.write(row, 3, "\n".join(c_query.l_articles))
                
            row = row + 1
            i = i + 1
            
        workbook.close()
        
        
    def do_search_combine_term(self):
        
        ## write in excel to have several lines in article list
        for term in self.l_terms:
            p_xlx = "%schem_%s.xlsx"%(self.pr_out, term)
            workbook = xlsxwriter.Workbook(p_xlx)
            worksheet = workbook.add_worksheet()
            l_headers = ["CASRN", "Query", "Nb Result", "Article"]
            l_chem = list(self.d_chem.keys())
            # wrtie header in the xlxs            
            [worksheet.write(0, i_col, l_headers[i_col]) for i_col in range(0, 4)]
            
            row = 1
            imax = len(list(self.d_chem.keys()))
            i = 0
            while i < imax:
                if pd.isna(self.d_chem[l_chem[i]]["PubMed search chemical"]) == True:
                    i = i + 1
                    continue
                
                c_query = PubMed.searching(self.email, "")
                c_query.buildNewQuery(l_chem[i], self.d_chem[l_chem[i]]["PubMed search chemical"])
                c_query.buildingQuery(term, "AND")
                c_query.search(w=False)
                
                # add casrn
                worksheet.write(row, 0, l_chem[i])
                # add query 
                worksheet.write(row, 1, c_query.Query)
                # add nb results
                worksheet.write(row, 2, c_query.count)
                # add list of article
                worksheet.write(row, 3, "\n".join(c_query.l_articles))
                
                row = row + 1
                i = i + 1
            
            workbook.close()

    def do_search(self, l_combine_search, review = "on", date = ""):
        self.load_tab_chem()
        self.do_search_bychem()
        self.do_search_combine_term()        
        


