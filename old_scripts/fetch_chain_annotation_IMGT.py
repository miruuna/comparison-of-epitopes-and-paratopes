from bs4 import BeautifulSoup
import requests
import re
import json

page_link = "http://www.imgt.org/3Dstructure-DB/cgi/3Dquery.cgi"
page_response = requests.get(page_link)
soup = BeautifulSoup(page_response.text, "html.parser")
pdb_codes=[]

#Find all the pdb codes associated with lysozyme
for tr in soup.find_all("tr"):
    for td in tr.find_all('td'):
        if "lysozyme" in td.text:
            for a_tag in tr.find_all('a'):
                pdb = a_tag.text
                pdb_codes.append(pdb)

pass_list={}
def test(codes):
    for code in codes:
        if code in pdb_codes:
            pass_list[code]=True
        else:
            pass_list[code]=False
    return pass_list

print(len(pdb_codes))

def get_chain_annotation(soup, part):
    """
    Fetches antigen chain annotation the IGMT's webpage
    :param soup
    :return: an array of letter chain names
    """
    all_chains=''
    chains_array=[]
    all_text =[]
    chains_list=''
    the_h3 = soup.find_all('h3', class_ = "barrevertelight")
    for h3 in the_h3:
    #all_text=all_text+'-----'+h3.text
        for strong_tag in h3.find_all('strong'):
            if part in strong_tag.text:
                for em_tag in strong_tag.find_all('em'):
                    the_text = em_tag.next_sibling
                    the_chain_code1 = the_text.replace('[','\\')
                    the_chain_code = the_chain_code1.replace(']','\\')
                    all_text.append(the_chain_code) #returns it as an array

    for series in all_text:
        p = re.compile(r'_(\w+)')
        the_chains = p.findall(series)
        for match in the_chains:
            all_chains = all_chains +' '+ match

    for amino_acid in all_chains.replace('\\',''):
        if (amino_acid!=' ') and (amino_acid not in chains_array):
            chains_array.append(amino_acid)
    for aa in sorted(chains_array):
        chains_list+=aa
    return chains_array


def fetch_chains_by_pdb(pdb, part):
    chains =[]
    page_link = "http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode="+pdb.upper()+"&Part=Epitope"
    page_response = requests.get(page_link, timeout=5)
    soup = BeautifulSoup(page_response.text, "html.parser")
    chains=get_chain_annotation(soup, part)
    return chains

def get_pdb_no_chains():
    """
    Checks if there is no antigen or antibody in the structure and adds it to a list
    :return: pdb_empyt-list of pdb codes which miss either the antibody or the antigen
    """
    pdb_empty=[]
    for pdb in pdb_codes:
        antigen=fetch_chains_by_pdb(pdb,"EC:")
        antibody=fetch_chains_by_pdb(pdb, "_KAPPA")
        if len(antibody)==0 or len(antigen)==0:
            pdb_empty.append(pdb)
    return (pdb_empty)


# Deletes pdbs with no antigen or antibody from the original generated list of pdb codes
for pdb in get_pdb_no_chains():
    if pdb in pdb_codes:
        pdb_codes.remove(pdb)



def store_data():
    """
    Stores data abbout the chain annotation os each pdb in a dictionary of dictionaries
    :return: 
    """
    pdb_dict={}
    for pdb in pdb_codes:
        pdb_dict[pdb]={}
        antigen=fetch_chains_by_pdb(pdb,"EC:")
        antibody=fetch_chains_by_pdb(pdb, "_KAPPA")
        pdb_dict[pdb]["antigen"]=[]
        pdb_dict[pdb]["antibody"]=[]

        for chain in antigen:
            if chain!='':
                pdb_dict[pdb]["antigen"].append(chain)
        for chain in antibody:
            if chain!='':
                pdb_dict[pdb]["antibody"].append(chain)
    return pdb_dict
print(store_data())


#store the annotation in a json file 
with open('data.json', 'w') as fp:
   json.dump(store_data(), fp)

