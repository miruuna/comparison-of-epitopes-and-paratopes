from bs4 import BeautifulSoup
import requests
import re
import json

pdb_codes= json.loads(open('pdb_codes_lysozyme.json').read())

def get_chain_annotation(soup, part):
    """
    Fetches antigen chain annotation the IGMT's webpage
    :param soup
    :return: an array of letter chain names
    """
    all_chains=''
    chains_array=[]
    all_text =[]
    species=[]
    chains_list=''
    the_div = soup.find_all('div', class_ = "TableResume")
    for div in the_div:
        #all_text=all_text+'-----'+h3.text
        for tr_tag in div.find_all('tr'):
            for td in tr_tag.find_all('td', class_="cell_milieu"):
                if part in td:
                    for em_tag in tr_tag.find_all('em'):
                        species.append(em_tag.text)
    return species



def fetch_chains_by_pdb(pdb, part):
    chains =[]
    page_link = "http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode="+pdb.upper()+"&Part=Chain"
    page_response = requests.get(page_link, timeout=5)
    soup = BeautifulSoup(page_response.text, "html.parser")
    chains=get_chain_annotation(soup, part)
    return chains


def store_data():
    """
    Stores data abbout the chain annotation os each pdb in a dictionary of dictionaries
    :return:
    """
    pdb_dict={}
    for pdb in pdb_codes:
        pdb_dict[pdb]=''
        antibody=fetch_chains_by_pdb(pdb, "IG")

        for species in antibody:
            if species!='':
                if species not in pdb_dict[pdb]:
                    pdb_dict[pdb]+=species
    return pdb_dict

with open("species_lysozime.json", 'w') as ctr:
    json.dump(store_data(), ctr)
print(store_data())
