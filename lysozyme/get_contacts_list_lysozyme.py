
import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt

annotation = json.loads(open('data_lysozyme.json').read())
pdb_codes1=["1tzg"]
pdb_codes= json.loads(open('pdb_codes_lysozyme.json').read())
contact_residues= json.loads(open('contact_residues_all_imgtt.json').read())

def get_contacts(chain, param):
    contacts={}
    for pdb in pdb_codes:
        a_dict={}
        pdb_light_=[]
        entry_list=[]
        antigen = annotation[pdb]['antigen'][0]
        #antibody= annotation[pdb]['antibody']
        if chain=="light":
            if len(annotation[pdb]['antibody'])>=2:
                if (elem=='L' for elem in annotation[pdb]['antibody']):
                    antibody= annotation[pdb]['antibody'][1]
                    non_h_cons=contact_residues[pdb][antigen][antibody][param]

            else:
                antibody= annotation[pdb]['antibody'][0]
                non_h_cons=[]

        elif chain=="heavy":
            if len(annotation[pdb]['antibody'])>=2:
                if (elem=='H' for elem in annotation[pdb]['antibody']):
                    antibody= annotation[pdb]['antibody'][0]
                    non_h_cons=contact_residues[pdb][antigen][antibody][param]
                elif (elem=='B' for elem in annotation[pdb]['antibody']):
                    antibody= annotation[pdb]['antibody'][1]
                    non_h_cons=contact_residues[pdb][antigen][antibody][param]
            elif len(annotation[pdb]['antibody'])==1:
                antibody= annotation[pdb]['antibody'][0]
                non_h_cons=contact_residues[pdb][antigen][antibody][param]

        else:
            print("no chain given")
        res1_pos=[]
        for elem in non_h_cons:
            res1_pos.append(elem['res1_pos'])
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

print(get_contacts('heavy','non_h_bonded')['1sq2'])
