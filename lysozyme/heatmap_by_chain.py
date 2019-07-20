import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt


contact_residues= json.loads(open('equiv_contacts_dict_lsyozyme_light.json').read())
species = json.loads(open("species_lysozime.json").read())

unique_species = [k for k,v in species.items() if list(species.values()).count(v)==1]
def get_contacts( name):
    res1_pos=[]
    contacts={}
    a_dict={}
    counter=0
    for pdb, v in contact_residues.items():
        contacts[pdb]={}
        for pos, names in v.items():
            try:
                if len(name)!=0:
                    dict={pos: len(names)}
                    #heavy_chain[pdb]=(res1_pos)
                    #print(pdb, res1_pos)g


            except KeyError:
                continue
            contacts[pdb].update(dict)

        df =pd.DataFrame.from_dict(contacts,orient='index')
    #df.set_index('pdb',inplace=True)
    df=df[natsorted(df.columns)]
    df.groupby(by=species,axis=1)
    fig, ax = plt.subplots(figsize=(100,30))
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
    plt.suptitle("Light chain", fontsize = 36) #title
    g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=.1,annot=True, annot_kws={"size": 20}, cmap="YlGnBu", cbar=False)#.set_title(chain+" chain: "+param)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=30)
    g.get_figure().savefig("v1_on_antigen_heatmap_"+"_"+name+"_Light.png")



get_contacts('lysozyme')
get_contacts( 'lysozyme')
get_contacts('lysozyme')
get_contacts('lysozyme')
