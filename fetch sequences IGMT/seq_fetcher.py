
from bs4 import BeautifulSoup
import requests
import re

pdb_codes = ['1A2Y', '1bql']


def get_chain(the_class):
    the_td = soup.find_all('td', class_ = the_class)
    for td in the_td:
        the_text = td.text
        p = re.compile(r'(.+\s+]\s+.+)Sequence\sin\sFASTA')
        the_link = p.finditer(the_text)
        for match in the_link:
            sequence = match.groups()
            chain_seq= str(sequence)
            linked_chain = re.sub(r'[\\n\s+\]\\n]', '',chain_seq)
            light_chain = str(linked_chain).replace('xa0', '')
            return(light_chain)

def get_light_chain(soup):
    return(get_chain("data_l"))

def get_heavy_chain(soup):
    return(get_chain("data_h"))

for pdb in pdb_codes:
    page_link = "http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode="+pdb
    page_response = requests.get(page_link, timeout=5)
    soup = BeautifulSoup(page_response.text, "html.parser")
    print(get_light_chain(soup))
    print(get_heavy_chain(soup))
