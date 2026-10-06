"""Fetch primary model sources only. Downloaded content is never executed."""
from pathlib import Path
import urllib.request, hashlib, json, concurrent.futures, datetime

HERE=Path(__file__).resolve().parent
SOURCES={
 'hfip-2022.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9298724/fullTextXML',
 'pmma-graphene-2021.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7962820/fullTextXML',
 'polyply-library.md':'https://raw.githubusercontent.com/marrink-lab/polyply_1.0/master/LIBRARY.md',
 'polyply-tree.json':'https://api.github.com/repos/marrink-lab/polyply_1.0/git/trees/master?recursive=1',
 'hfip-gaff-2023.pdf':'https://flore.unifi.it/retrieve/ddfb0dbc-48b0-409a-8a5a-9e157872802f/Journal%20of%20Peptide%20Science%20-%202023%20-%20Casoria%20-%20Upgrading%20of%20the%20general%20AMBER%20force%20field%202%20for%20fluorinated%20alcohol.pdf'
}

def fetch(item):
 name,url=item
 path=HERE/'sources'/name
 request=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 (local scientific research)'})
 try:
  with urllib.request.urlopen(request,timeout=35) as r:
   data=r.read(70_000_001)
   if len(data)>70_000_000: raise ValueError('source exceeds 70 MB limit')
   actual=r.url
   content_type=r.headers.get('Content-Type','')
  path.write_bytes(data)
  return dict(name=name,url=url,resolved=actual,content_type=content_type,bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
 except Exception as e:
  return dict(name=name,url=url,error=str(e))

if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
  records=list(pool.map(fetch,SOURCES.items()))
 (HERE/'sources/downloads.json').write_text(json.dumps(records,indent=2)+'\n')
 for r in records: print(json.dumps(r))
