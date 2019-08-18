# Dependencies

import distance
import matplotlib.pyplot as plt
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from natsort import natsorted


g1=['1c08', '1dqj', '1ic4', '1ic5', '1ic7', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '4i0c']
g2=['1ri8', '1zvh',  '1ihl']
g3=['1mlc','2iff', '1zv5','1a2y','2i25']
g4=['1p2c','1ndg','1ndm','1op9']
g5=['2i26', '1zmy','1zvy','1jto']
species = json.loads(open("species_lysozime.json").read())
pdb_codes= json.loads(open('pdb_codes_lysozyme.json').read())
contacts=json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_all_all.json').read())
non_bonded_dict= json.loads(open('equiv_contacts_dict_lsyozyme_all_non_h_bonded.json').read())
h_bonded_dict=json.loads(open('equiv_contacts_dict_lsyozyme_all_h_bonded.json').read())
unique1=['1a2y', '1bql', '1c08', '1dqj', '1dzb', '1ic4', '1ic5', '1ic7', '1jhl', '1kip', '1kiq', '1mlc', '1ndg', '1ndm', '1p2c', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '2iff', '1jto', '1op9', '1ri8', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '4i0c']
unique=['1sq2','1t6v','2hfm','2i25','2i26','1a2y', '1bql', '1c08', '1dqj', '1dzb', '1ic4', '1ic5', '1ic7', '1jhl', '1kip', '1kiq', '1mlc', '1ndg', '1ndm', '1p2c', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '2iff', '4tsa', '4tsb', '4tsc', '4ttd', '1jto', '1op9', '1ri8', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '4i0c']
c_and_m=['1jto', '1jtp', '1jtt', '1mel', '1op9', '1ri8', '1rjc', '1sq2', '1t6v', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '2i25', '2i26', '3eba', '4i0c']

unique_minus_1=set(unique).difference(set(g1))


def get_pdb_by_group(group):
    species_list=[]
    species_dict={}
    for k, v in species.items():
        if v not in species_list:
            species_list.append(v)
    for sp in species_list:
        species_dict[sp]=[]
        for k,v in species.items():
            if k in group:
                if sp==v:
                    species_dict[sp].append(k)

    return species_dict
def get_jaccard_score(group,name):
    non_bonded_list={}
    h_bonded_list={}
    for pdb in group:
        for k, v in contacts.items():
            if pdb in k:
                non_bonded_list[pdb]=[]
                for pos, res in v.items():
                    non_bonded_list[pdb].append(pos)
    d=non_bonded_list

    #for pdb in pdb_codes:
        #for k, v in h_bonded_dict.items():
            #if pdb in k:
                #h_bonded_list[pdb]=[]
               # for pos, res in v.items():
                 #   h_bonded_list[pdb].append(pos)
    #if type=="h_bonded":
        #d=h_bonded_list
    #if type=="non_h_bonded":
        #d=non_bonded_list

    a={'a': [2 , 3, 4],
       'b': [3,4,5],
       'c':[4,1,3],
       'd':[5,1,2]}

    p={}
    count=0
    key_list=[]
    for key, v in d.items():
        key_list.append(key)

    for k in key_list:
        p[k]=[]
        cos1=0
        for i in key_list:
            if i!=k:
                while(len(d[k])!=len(d[i])):
                    if len(d[k])> len(d[i]):
                        d[i].append(0)
                    elif len(d[k])< len(d[i]):
                        d[k].append(0)
                p[k].append(1-distance.jaccard(set(d[k]),set(d[i])))
        sc_sum=0
        dict_sum={}
        dict_avg={}
        for k, v in p.items():

            dict_sum[k]=0
            for e in v:
                dict_sum[k]+=e
            dict_avg[k]=dict_sum[k]/len(p[k])

    #with open('all_jaccard_scores.json', 'w') as ctr: #save the contact residues in a json file
    #json.dump(dict_avg, ctr)
    #p_sorted={k: v for k, v in sorted(dict_avg.items(), key=lambda x: x[1])}
    name_list=[]
    score_list=[]
    species_list=[]
    for k, v in dict_avg.items():
            name_list.append(k)
            score_list.append(v)
            species_list.append(species[k])

    df= pd.DataFrame({'Jaccard_score':score_list,
                      'PDB_codes':name_list,
                      'species': species_list
                      })


    df.sort_values("species", inplace=True)
    df.index = np.arange(1, len(df)+1)
    sorted=df.sort_values('Jaccard_score')
    g=sns.scatterplot(y="Jaccard_score", x="PDB_codes", data=sorted, hue='species')
    #g.set(xticks=[])
    plt.xticks(rotation=90)
    g.get_figure().savefig(name+"_jaccard_scatter_plot.png")
    df.to_excel(name+"_jaccard_scores.xlsx")
    return df
    #with open('mytable.txt', 'w') as tf:
    #   tf.write(df.to_latex())

def get_pie_chart(group, name):
    species_count={}
    for k, v in get_pdb_by_group(group).items():
        species_count[k]=len(v)
    df=pd.DataFrame.from_dict(species_count, orient="Index")
    df.to_excel("species_distrib_"+name+".xlsx")
    return df



print(get_jaccard_score(c_and_m, "cm"))
#print(unique_minus_1)

#print(get_pie_chart(unique, "unique"))
