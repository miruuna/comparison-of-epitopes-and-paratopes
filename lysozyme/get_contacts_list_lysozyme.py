#This script contains the function get_contacts() which returns a nested dictionary of the interaction residues for the antibody/aantigen and their different chains

import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt

annotation = json.loads(open('data_lysozyme.json').read())
pdb_codes1=["1tzg"]
pdb_codes= json.loads(open('pdb_codes_lysozyme.json').read())
species=json.loads(open('species_lysozime.json').read())
contact_residues= json.loads(open('contact_residues_all_imgtt.json').read())

def get_contacts(chain, param):
    contacts={}
    for pdb in pdb_codes:
        a_dict={}
        pdb_light_=[]
        entry_list=[]

        for e in annotation[pdb]['antigen']:
            if e in ['C','Y','L']:
                antigen=e
            else:
                antigen=annotation[pdb]['antigen'][0]


        #antibody= annotation[pdb]['antibody']

        if chain=="light":
            if len(annotation[pdb]['antibody'])>=2:
                if pdb in ['1jto', '1jtt', '1jtp']:
                    non_h_cons=[]

                else:
                    for elem in annotation[pdb]['antibody']:
                        if elem=='L':
                            non_h_cons=contact_residues[pdb][antigen][elem][param]
                        if elem=='A':
                            non_h_cons=contact_residues[pdb][antigen][elem][param]


            else:
                non_h_cons=[]

        elif chain=="heavy":
            if len(annotation[pdb]['antibody'])>=2:
                if pdb in ['1jto', '1jtt', '1jtp']:
                    non_h_cons=contact_residues[pdb]['L']['A'][param]
                    if non_h_cons==[]:
                        non_h_cons=contact_residues[pdb]['M']['B'][param]
                else:
                    for elem in annotation[pdb]['antibody']:
                        if str(elem) =='H':
                            non_h_cons=contact_residues[pdb][antigen][elem][param]
                        elif elem=='B':
                            non_h_cons=contact_residues[pdb][antigen][elem][param]
                        elif elem=='N':
                            non_h_cons=contact_residues[pdb][antigen][elem][param]

            elif len(annotation[pdb]['antibody'])==1:
                antibody= annotation[pdb]['antibody'][0]
                non_h_cons=contact_residues[pdb][antigen][antibody][param]

        else:
            print("no chain given")
        res1_pos={}
        pos1=[]
        for elem in non_h_cons:
            if elem not in pos1:
                pos1.append(elem['res1_pos'])
                pos=sorted(pos1)

        for res in pos:
            res1_pos[res]=[]
            for elem in non_h_cons:
                if res==elem['res1_pos']:
                    res1_pos[res].append(elem['res1_name'])


            no_elem=non_h_cons.count(elem)
        #dict={pos: int(res1_pos.count(pos)) for pos in res1_pos}
        #heavy_chain[pdb]=(res1_pos)
        #print(pdb, res1_pos)g
        res1_pos_unique=[]
        for c in res1_pos:
            if c not in res1_pos_unique:
                res1_pos_unique.append(c)

        contacts[pdb]=res1_pos


    return contacts
print(get_contacts('light', 'h_bonded'))


