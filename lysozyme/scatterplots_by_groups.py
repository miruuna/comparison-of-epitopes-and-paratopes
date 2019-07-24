import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np
group1=['1bql', '1mlc', '1nby', '1nbz', '1ndg', '1ndm', '1p2c', '1ri8', '1rjc', '1yqv', '1zv5', '1zvh', '2iff', '3d9a']
group2=['1fdl', '1jto', '1jtp', '1jtt', '1xfp'] #removed 1zmy
group3a=['1a2y', '1g7h', '1g7i', '1g7j', '1g7l', '1g7m', '1kip', '1kiq', '1kir',  '1vfb']
group3b=['1sq2', '1t6v', '2hfm', '2i25', '2i26', '3eba', '4i0c']
groupiv=['1c08', '1ic5', '1ic7', '1ua6', '1uac', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqd', '2dqe', '2dqf', '2dqg', '2eiz', '2eks', '2yss', '3a6b', '3a6c', '3hfm']
contacts= json.loads(open('equiv_contacts_dict_lsyozyme_heavy.json').read())
species=json.loads(open('species_lysozime.json').read())


d_g1={}
pos=[]
a_dict={}
g1={}
g11={}
for k, v in contacts.items():

    for pdb in group1:
        if pdb in k:
            g1[pdb]={}
            g11[pdb]={}

            for a, b in v.items():

                res1=''
                resi=[]
                d1={}
                for elem in b:
                    if elem not in resi:
                        res1+=elem+"-"
                        resi.append(elem)
                g1[pdb][int(a)]=res1
                if len(res1)==4:
                    g11[pdb][int(a)]=res1
                else:
                    g11[pdb][int(a)]="Multiple pref"

                d1={res : int(resi.count(res)) for res in resi}
                a_dict[a]=d1
"""
df =pd.DataFrame.from_dict(a_dict,orient='index').T
df=df[natsorted(df.columns)]
df = df.replace(np.nan, 0)
fig, ax = plt.subplots(figsize=(100,30))
fig.autofmt_xdate()
plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
plt.suptitle("Amino acid preference ", fontsize = 36) #title
g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=.3,annot=True, annot_kws={"size": 20}, cmap="GnBu_r", cbar=False)#.set_title(chain+" chain: "+param)
ax.tick_params(labelbottom='off',labeltop='on', labelsize=30)
g.get_figure().savefig("aa_preference_group3a_hm.png")

"""

"""
df = df.replace(np.nan, 0)
#df.set_index('pdb',inplace=True)
df=df[natsorted(df.columns)]
fig, ax = plt.subplots(figsize=(100,30))
fig.autofmt_xdate()
plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
plt.suptitle("try", fontsize = 36) #title
g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=.1,annot=True, annot_kws={"size": 20}, cmap="YlGnBu", cbar=False)#.set_title(chain+" chain: "+param)
ax.tick_params(labelbottom='off',labeltop='on', labelsize=30)
g.get_figure().savefig("a1.png")
"""
#plot 2
#f, ax = plt.subplots()
df2=pd.DataFrame.from_dict(g1,orient='index')
#df2 = df2.replace(np.nan, "-")
df2.index.rename('PDB', inplace=True)
stacked = df2.stack().reset_index()

stacked.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
stacked=stacked.sort_values(by=['Position'])
stacked2=stacked.stack().reset_index()

#fig, ax = plt.subplots(figsize=(100,30))

sns.relplot(data=stacked, x='Position', y='PDB', hue='Residues', s=20 )

plt.savefig("aa_preference_group1_H.png")
sns.set_palette("husl")
#snsheatmap(stacked)
plt.show()
print(stacked)
