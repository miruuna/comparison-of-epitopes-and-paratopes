
#This script contains the functions for computing the equivalent positions of the interating Ag residues.
#Equivalent positions were calculated using the multiple sequence alignment obtained from EBI'S Clustal Omega.
#The equivalebnt positions for different chains or combination of chains are saved in json format to be used in further analyses
#The multiple sequence alignment page is only stored for 2 weeks, thus page_link has to be changed accordingly

from collections import defaultdict
from itertools import chain

import requests
import re
from bs4 import BeautifulSoup
import requests
import os
import get_contacts_list_lysozyme as epitope
rootdir = os.getcwd()
import json



def hasCharacters(inputString):
    return any(char.isalpha() for char in inputString)

pdb_codes= json.loads(open('pdb_codes_lysozyme.json').read())
species=json.loads(open('species_lysozime.json').read())
germ_line=json.loads(open("germ_lines_lysozime.json").read())


unique_species = [k for k,v in species.items() if list(species.values()).count(v)==1]

page_link = 'http://www.ebi.ac.uk/Tools/services/rest/clustalo/result/clustalo-E20190804-101942-0556-43333684-p1m/aln-clustal_num'
page_response = requests.get(page_link, timeout=5)
text = BeautifulSoup(page_response.text, "html.parser")


counter=0 #countning how many different pdbs in the alignment
the_pdb_list =[]
text_1=str(text).replace('PDB|', '>')
alignment=text_1.replace('|', '')
p = re.compile(r'>[a-zA-Z0-9]{4}')
the_pdbs = p.findall(alignment)
for aa in the_pdbs:
    if aa not in the_pdb_list:
        the_pdb_list.append(str(aa))
#print(alignment)
#print(the_amino_acid_list)


#capture the sequences line by line and add it to dictionary having pdb codes as keys and  the sequences(As a astring) as items
def get_sequences_from_alignment():
    sequences={}
    for pdb in the_pdb_list:
        sequences[pdb[1:]] =''
        q=re.compile(r'{}\w\s+(.+)\s+\d'.format(pdb))
        for string in q.findall(alignment):
            sequences[pdb[1:]]+=(str(string))
    return(sequences)



def get_sequences_from_pdb():
    pdb_seq={}
    elem_list =[]
    with open(rootdir+'\\gp120_for_MSA.txt', 'r') as f:
        contents = f.read()
        contents2 = contents.rstrip("\n")
        string1=str(contents2).replace('>PDB', ':')
        string2=string1.replace('|', '')
        string3=string2.replace('\n', '')
        string3+=':'
        p=re.compile(r':(\w+):')
        for item in p.findall(string3):
            elem_list.append(item)
    for elem in elem_list: #separate previous string.first 4 characters are the pdb code
        pdb_seq[elem[:4]]=elem[5:]
    return pdb_seq

def get_gaps():
    gap_counter=0
    pos_dict={}
    pos_list={}
    equiv_pos_dict={}
    pos_pdb={}
    d =get_sequences_from_alignment()
    for k in d.keys():
        equiv_pos_dict[k]={}
        for pos in range(1, len(d[k])+1):
            pos_dict[pos]=d[k][:pos].count('-')
            if pos!='-':
                pos_pdb[pos]=pos-pos_dict[pos]
            else:
                pos_pdb[pos]=0
            equiv_pos_dict[k].update(pos_pdb)

    return(equiv_pos_dict)



def get_contacts():
    all_dict={}
    for pdb in pdb_codes:
        x= epitope.get_contacts("heavy", "h_bonded")[pdb]
        y= epitope.get_contacts("heavy", "non_h_bonded")[pdb]
        z= dict(y, **x)
        all_x_y = defaultdict(list)
        for k,v in chain(x.items(), y.items()):
            all_x_y[k].extend(v)
        a= epitope.get_contacts("light", "h_bonded")[pdb]
        b= epitope.get_contacts("light", "non_h_bonded")[pdb]
        all_a_b = defaultdict(list)
        for k,v in chain(a.items(), b.items()):
            all_a_b[k].extend(v)
        all = defaultdict(list)
        for k,v in chain(all_a_b.items(), all_x_y.items()):
            all[k].extend(v)
        all_dict[pdb]=dict(all)
    return all_dict




def get_equivalent_contacts():
    contacts_dict=get_contacts()
    equiv_contacts={}
    equiv_contacts_dict={}
    equiv_pos_alignment = get_gaps()
    for pdb in pdb_codes:
        #quiv_contacts[pdb]= {}
        eq={}
        for pdb, v in contacts_dict.items():
            eq[pdb]=[]
            for elem, value in v.items():
                for k in equiv_pos_alignment[pdb].keys():
                    if str(equiv_pos_alignment[pdb][k])!='-':
                        if str(elem) == str(equiv_pos_alignment[pdb][k]):
                            eq[pdb].append(k)


            equiv_contacts[pdb]=eq[pdb]
    return equiv_contacts



with open('equiv_contacts_list_all_contacts_by_pdb.json', 'w') as ctr: #save the contact residues in a json file
    json.dump(get_equivalent_contacts(), ctr)
