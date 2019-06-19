#!/usr/bin/env python3

"""
------------------------------------------------------------------------------------------------------------------------
Author: Miruna Serian
Description: Extracts the chain annotations for given pdb structures
------------------------------------------------------------------------------------------------------------------------
"""

from bs4 import BeautifulSoup
import requests
import re
import get_residues as residues
pdb_codes = residues.get_subdir()




def get_antigen_chains(soup):
    """
    Fetches antigen chain annotation the IGMT's webpage
    :param soup
    :return: an array of letter chain names
    """
    all_chains=''
    chains_list=[]
    all_text =[]
    the_h3 = soup.find_all('h3', class_ = "barrebleue")
    for h3 in the_h3:
        #all_text=all_text+'-----'+h3.text
        for strong_tag in h3.find_all('strong'):
            if 'Hemagglutinin' in strong_tag.text:
                for em_tag in h3.find_all('em'):
                    the_text = em_tag.next_sibling
                    the_chain_code1 = the_text.replace('[','\\')
                    the_chain_code = the_chain_code1.replace(']','\\')
                    all_text.append(the_chain_code) #returns it as an array
    for series in all_text:
        p = re.compile(r'_(.+)')
        the_chains = p.findall(series)
        for match in the_chains:
            all_chains = all_chains +' '+ match

    for amino_acid in all_chains.replace('\\',''):
        if (amino_acid!=' ') and (amino_acid not in chains_list):
            chains_list.append(amino_acid)
    return(sorted(chains_list))


def get_antibody_chains(soup):
    """
    Fetches antibody chain annotation the PyIgClassify's webpage
    :param soup
    :return: an array of letter chain names
    """
    the_chains = [ ]
    for tr in soup.find_all('tr', class_= "GridRow"):
        the_string = tr.find_all('td')[0].text #returns the text of the first td tag
        if str(the_string) not in the_chains:
            the_chains.append(str(the_string))
    for tr in soup.find_all('tr', class_= "GridAlternativeItem"):
        the_string2 = tr.find_all('td')[0].#returns the text of the first td tag
        if str(the_string2) not in the_chains:
            the_chains.append(str(the_string2))
    return sorted(the_chains)




def fetch_antigen_chains():
    antigen_chains = {}
    for pdb in pdb_codes:
        page_link = "http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode="+pdb
        page_response = requests.get(page_link, timeout=5)
        soup = BeautifulSoup(page_response.text, "html.parser")
        antigen_chains[pdb]= get_antigen_chains(soup)
    return antigen_chains

def fetch_antibody_chains():
    light_chains ={}
    for pdb in pdb_codes:
        page_link = "http://dunbrack2.fccc.edu/PyIgClassify/Result/ChainCdrClusterInfo.aspx?pdb="+pdb.upper()
        page_response = requests.get(page_link, timeout=5)
        soup = BeautifulSoup(page_response.text, "html.parser")
        light_chains[pdb]= get_antibody_chains(soup)
    return light_chains

print(fetch_antigen_chains())
