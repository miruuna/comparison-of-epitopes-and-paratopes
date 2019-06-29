import json
from bs4 import BeautifulSoup
import requests
import re


annotation = json.loads(open('data.json').read())
pdb_codes=["1a2y", "1c08", "1dqj", "1fdl", "1g7h", "1g7i", "1g7j", "1g7l", "1g7m", "1ic7", "1j1o", "1j1p", "1j1x", "1jhl", "1kip", "1kiq", "1kir", "1mlc", "1nby", "1nbz", "1ndg", "1ndm", "1p2c", "1ua6", "1uac", "1vfb", "1xgp", "1xgq", "1xgr", "1xgt", "1xgu", "1yqv", "2dqc", "2dqd", "2dqe", "2dqf", "2dqg", "2dqh", "2dqi", "2dqj", "2eiz", "2eks", "2hfm", "2iff", "2yss", "3a6b", "3a6c", "3d9a", "3hfm"]

anno_dict={}
def get_annotation(pdb):
    for antigen_chain in annotation[pdb]['antigen']:
        antibody_dict={}
        for antibody_chain in annotation[pdb]['antibody']:
            page_link = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb="+pdb+"&chain1="+antibody_chain+"&chain2="+antigen_chain
            page_response = requests.get(page_link)
            soup = BeautifulSoup(page_response.text, "html.parser")
            line_list=[]
            for p in soup.find_all('p'):
                for line in soup.text.split("\n"):
                    line_list.append(line)
            search_h_bonded=line_list[line_list.index('Hydrogen bonds')+7:line_list.index('Non-bonded contacts')-1]
            search_non_bonded=line_list[line_list.index('Non-bonded contacts')+7:line_list.index(line_list[-4])]
            search_dict={
                "h_bonded": search_h_bonded,
                "non_bonded": search_non_bonded
            }
            counter=0
            big_dict={}
            big_dict["h_bonded"]={}
            res_anno_bonded=[]
            res_anno_non_bonded=[]
            big_dict["non_bonded"]={}
            for line in search_dict["h_bonded"]:
                column=line.split()
                dict_res_annotation={
                    "res1_name": column[3],
                    "res1_pos": column[4],
                    "res2_name": column[9],
                    "res2_pos":column[10]
                    }
                res_anno_bonded.append(dict_res_annotation)
            big_dict["h_bonded"]=(res_anno_bonded)

            for line in search_dict["non_bonded"][:-1]:
                column=line.split()
                dict_res_annotation1={
                    "res1_name": column[3],
                    "res1_pos": column[4],
                    "res2_name": column[9],
                    "res2_pos":column[10]
                }

                res_anno_non_bonded.append(dict_res_annotation1)
            big_dict["non_bonded"]=(res_anno_non_bonded)

            antibody_dict[antibody_chain]=big_dict
        anno_dict[antigen_chain]=antibody_dict
    return anno_dict

#with open('contact_residues_all_imgt.json', 'w') as ctr:
 #   json.dump(anno_dict, ctr)
##annotation_all_pdbs={}
#for pdb in pdb_codes:
    #annotation_all_pdbs[pdb]=get_annotation(pdb)
all_annotation={}
def main_run():
    for pdb in pdb_codes:
        try:
            all_annotation[pdb]=get_annotation(pdb)
        except:
            pass
    return all_annotation
with open('all_annotation.json', 'w') as pdbc:
    json.dump(main_run(), pdbc)
