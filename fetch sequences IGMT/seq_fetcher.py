

from bs4 import BeautifulSoup
import requests
import re

pdb_codes = ['1A2Y', '1bql']


def get_chain(the_class):
    """
    Fetches amino acid chains the IGMT's webpage
    :param the_class: 
    :return: an amino acid residues chain as a string
    """
    the_td = soup.find_all('td', class_ = the_class) #the td tag has 'data_l' class for light chain and 'data_h' for heavy chain
    for td in the_td:
        the_text = td.text
        p = re.compile(r'(.+\s+]\s+.+)Sequence\sin\sFASTA')
        the_link = p.finditer(the_text)
        for match in the_link:
            sequence = match.groups()
            chain_seq= ''.join(sequence)
            linked_chain = re.sub(r'[\\n\s+\]\\n]', '',chain_seq)
            light_chain = str(linked_chain).replace('xa0', '')
            return(''.join(light_chain))

def get_light_chain(soup):
    return(get_chain("data_l"))

def get_heavy_chain(soup):
    return(get_chain("data_h"))

#Returns a dictionary comprising of pdb codes as keys and their light and heavy chains
chain_seq = {}
for pdb in pdb_codes:
    page_link = "http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode="+pdb
    page_response = requests.get(page_link, timeout=5)
    soup = BeautifulSoup(page_response.text, "html.parser")
    chain_seq[pdb]= get_light_chain(soup),get_heavy_chain(soup)

print(chain_seq['1A2Y'])
