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
                self.Query = "(" + self.Query + ") " + connect.upper() + " (" + strin + " [" + keyQuery + "]" ")"
                

    def search(self, w=False):

        if not "time" in self.__dict__:
            self.time = time.time()
        else:
            # every 4s one request
            gaptime = time.time() - self.time
            if gaptime < 4.0:
                time.sleep(4 - gaptime)
                self.time = time.time()

        handle = Entrez.esearch(db="pubmed", term=self.Query, retmax="10000000")
        result = Entrez.read(handle)
        
        l_id = result["IdList"]
        l_title = []
        for id in l_id:
            pubmed_entry = Entrez.efetch(db="pubmed", id=id, retmode="xml")
            pubmed_results = Entrez.read(pubmed_entry)
            try: article = pubmed_results['PubmedArticle'][0]['MedlineCitation']['Article']
            except:continue
            try:date = article["ArticleDate"][0]["Year"]
            except: 
                try:date = article["Journal"]["JournalIssue"]["PubDate"]["Year"]
                except: date = "----"
                    
            w_article = "%s; %s; %s"%(article["ArticleTitle"], date, id)
            l_title.append(w_article)
        if w == True:
            self.fsum.write("\"%s\"\t\"%s\"\t\"%s\"\t\"%s\"\n"%(self.nameQ, self.Query.replace("\"", "'"), result["Count"], "\n".join(l_title)))
        self.l_articles = l_title
        self.count = result["Count"]

