import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt

annotation = json.loads(open('data.json').read())
pdb_codes1=["2eks"]
pdb_codes= json.loads(open('pdb_codes.json').read())
contact_residues= json.loads(open('contact_residues_all_imgtt.json').read())

def get_contacts(chain, param):
    contacts=[]
    for pdb in pdb_codes:
        a_dict={}
        try:
            pdb_light_=[]
            entry_list=[]
            antigen = annotation[pdb]['antigen'][0]
            #antibody= annotation[pdb]['antibody']
            if chain=="light":
                if (elem=='L' for elem in annotation[pdb]['antibody']):
                    antibody= annotation[pdb]['antibody'][1]
                    non_h_cons=contact_residues[pdb][antigen][antibody][param]

                else:
                    antibody= annotation[pdb]['antibody'][0]
                    non_h_cons=contact_residues[pdb][antigen][antibody][param]

            elif chain=="heavy":
                if (elem=='H' for elem in annotation[pdb]['antibody']):
                    antibody= annotation[pdb]['antibody'][0]
                    non_h_cons=contact_residues[pdb][antigen][antibody][param]
                elif (elem=='B' for elem in annotation[pdb]['antibody']):
                    antibody= annotation[pdb]['antibody'][1]
                    non_h_cons=contact_residues[pdb][antigen][antibody][param]

            else:
                print("no chain given")
            res1_pos=[]
            for elem in non_h_cons:
                res1_pos.append(elem['res1_pos'])
                no_elem=non_h_cons.count(elem)
            dict={pos: int(res1_pos.count(pos)) for pos in res1_pos}
            #heavy_chain[pdb]=(res1_pos)
            #print(pdb, res1_pos)g
            a_dict[pdb]=dict
        except KeyError:
            continue
        contacts.append(a_dict)

    df = pd.DataFrame.from_records([{'pdb': k, **v} for d in contacts for k,v in d.items()])
    df.set_index('pdb',inplace=True)
    df=df[natsorted(df.columns)]

    fig, ax = plt.subplots(figsize=(60,30))
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
    plt.suptitle(chain+" chain: "+param, fontsize = 36) #title
    g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=.1, annot=True, annot_kws={"size": 20}, cmap="YlGnBu", cbar=False)#.set_title(chain+" chain: "+param)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=25)
    g.get_figure().savefig("heatmap_"+chain+"_"+param+".png")

get_contacts('light','h_bonded')
get_contacts('heavy','h_bonded')
get_contacts('light','non_h_bonded')
get_contacts('heavy','non_h_bonded')
