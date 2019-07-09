import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt

annotation = json.loads(open('data.json').read())
pdb_codes1=["2eks"]
pdb_codes= json.loads(open('pdb_codes.json').read())
contact_residues= json.loads(open('contact_residues_all_imgtt.json').read())

heavy_chain=[]
for pdb in pdb_codes:
    a_dict={}
    try:
        pdb_light_=[]
        entry_list=[]
        antigen = annotation[pdb]['antigen'][0]
        antibody= annotation[pdb]['antibody']

        non_h_cons=contact_residues[pdb][antigen][antibody[1]]['non_h_bonded']
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
    heavy_chain.append(a_dict)

df = pd.DataFrame.from_records([{'pdb': k, **v} for d in heavy_chain for k,v in d.items()])
df.set_index('pdb',inplace=True)
the_columns=list(df.columns)
sorted_col=[]
for elem in the_columns:
    sorted_col.append(int(elem))

df=df[natsorted(df.columns)]

print(df)
print(natsorted(df.columns))
fig, ax = plt.subplots(figsize=(60,20))
g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=.1, annot=True, annot_kws={"size": 20}, cmap="YlGnBu", cbar=False)
ax.tick_params(labelbottom='off',labeltop='on', labelsize=20)
g.get_figure().savefig("heatmap.png")
