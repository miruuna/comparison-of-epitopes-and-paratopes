import json
from bs4 import BeautifulSoup
import requests
import re

#read the chain annotation dictionary  and the pdb codes array saved in json format
annotation = json.loads(open('data.json').read())
pdb_codes1= json.loads(open('pdb_codes.json').read())
pdb_codes=['1a2y']



def get_annotation(pdb):
    anno_dict={}
    for antigen_chain in annotation[pdb]['antigen']:
        antibody_dict={}
        for antibody_chain in annotation[pdb]['antibody']:
            page_link = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb="+pdb+"&chain1="+antibody_chain+"&chain2="+antigen_chain
            page_response = requests.get(page_link)
            soup = BeautifulSoup(page_response.text, "html.parser")
            line_list=[]
            for p in soup.find_all('p'):
                for line in soup.text.split("\n"):
                    line_list.append(line) #for each page it generates a list of the text lines in the body of the page


            #search_h_bonded=line_list[line_list.index('Hydrogen bonds')+7:line_list.index('Non-bonded contacts')-1]
            big_dict={}

            if "Hydrogen bonds" not in line_list: #search for line which contains "Hydrogen bonds
                big_dict['h_bonded']=[]

            else:
                search_h_bonded=line_list[line_list.index('Hydrogen bonds')+7:line_list.index('Non-bonded contacts')-1]
                res_anno_bonded=[] #a list containing dictionaries of contact residues with their names and positions
                for line in search_h_bonded:
                    column=line.split()
                    dict_res_annotation1={
                        "res1_name": column[3],
                        "res1_pos": column[4],
                        "res2_name": column[9],
                        "res2_pos":column[10]
                    }

                    res_anno_bonded.append(dict_res_annotation1)
                big_dict["h_bonded"]=(res_anno_bonded)
            res_anno_non_bonded=[]
            if 'Non-bonded contacts' not in line_list:
                big_dict['non_h_bonded']=[]

            elif "Salt bridges" in line_list:
                search_non_bonded_all=line_list[line_list.index('Non-bonded contacts')+7:line_list.index(line_list[-1])]
                search_non_bonded=line_list[:line_list.index("Salt bridges")-1]

            else:
                search_non_bonded=line_list[line_list.index('Non-bonded contacts')+7:line_list.index(line_list[-6])]


                for line in search_non_bonded:
                    column=line.split()
                    dict_res_annotation1={
                        "res1_name": column[3],
                        "res1_pos": column[4],
                        "res2_name": column[9],
                        "res2_pos":column[10]
                    }

                    res_anno_non_bonded.append(dict_res_annotation1)
                big_dict["non_h_bonded"]=(res_anno_non_bonded)



            antibody_dict[antibody_chain]=big_dict
        anno_dict[antigen_chain]=antibody_dict

    return anno_dict


annotation_all_pdbs={}
for pdb in pdb_codes1:
    annotation_all_pdbs[pdb]=get_annotation(pdb)
    #try:
        #annotation_all_pdbs[pdb]=get_annotation(pdb)
    #except:
        #continue
with open('contact_residues_all_imgtt.json', 'w') as ctr: #save the contact residues in a json file
    json.dump(annotation_all_pdbs, ctr)
