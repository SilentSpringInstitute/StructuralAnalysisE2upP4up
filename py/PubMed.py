from Bio import Entrez
import datetime
import time

class searching:
    def __init__(self, email, fsum):

        Entrez.email = email
        now = datetime.datetime.now()
        self.now = now
        self.fsum = fsum
        self.Qlog = "%s-%s-%s,%s:%s:%s"%(now.month, now.day, now.year, now.hour, now.minute, now.second)
        self.err = 0

    def buildNewQuery(self, name, strin, keyQuery=""):

        self.nameQ = name
        self.Query = strin
        if keyQuery != "":
            self.Query = self.Query + " [" + keyQuery + "]"


    def buildingQuery(self, strin, connect, keyQuery = ""):

        if not "Query" in self.__dict__:
            self.buildNewQuery(strin, keyQuery)
        else:
            if keyQuery == "":
                self.Query = "(" + self.Query + ") " + connect.upper() + " (" + strin + ")"
            else:
                self.Query = "(" + self.Query + ") " + connect.upper() + " (\"" + strin + "\" [" + keyQuery + "]" ")"
                

    def search(self, w=False, limit_extract=500):

        if not "time" in self.__dict__:
            self.time = time.time()
        else:
            # every 1.2s one request
            gaptime = time.time() - self.time
            if gaptime < 1.2:
                time.sleep(1.2 - gaptime)
                self.time = time.time()

        try:
            handle = Entrez.esearch(db="pubmed", term=self.Query, retmax="10000000")
            result = Entrez.read(handle)
        except:
            self.err = 1
            return
        
        # here limit the extraction to not overload the sys
        if int(result["Count"]) > limit_extract:
            l_titles = ["More than %s articles (n=%s) -> Please refine the query"%(limit_extract, result["Count"])]
            l_abstracts = [""]
            l_id = result["IdList"]
        else:
            l_id = result["IdList"]
            l_titles = []
            l_abstracts = []
            for id in l_id:
                pubmed_entry = {}
                while pubmed_entry == {}:
                    try:pubmed_entry = Entrez.efetch(db="pubmed", id=id, retmode="xml")
                    except: time.sleep(2)
                pubmed_results = Entrez.read(pubmed_entry)
                try: article = pubmed_results['PubmedArticle'][0]['MedlineCitation']['Article']
                except:continue
                try:date = article["ArticleDate"][0]["Year"]
                except: 
                    try:date = article["Journal"]["JournalIssue"]["PubDate"]["Year"]
                    except: date = "----"
                
                try:abstract = article["Abstract"]["AbstractText"][0]
                except: abstract=""
                w_article = "%s; %s; %s"%(article["ArticleTitle"], date, id)
                l_titles.append(w_article)
                l_abstracts.append(abstract)
        
        if w == True:
            self.fsum.write("\"%s\"\t\"%s\"\t\"%s\"\t\"%s\"\n"%(self.nameQ, self.Query.replace("\"", "'"), result["Count"], "\n".join(l_titles)))
        self.l_titles = l_titles
        self.l_abstracts = l_abstracts
        self.pubmedID = l_id
        self.count = result["Count"]

