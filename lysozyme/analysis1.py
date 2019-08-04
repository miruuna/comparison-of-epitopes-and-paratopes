import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
c2=['1ic5', '1ic7']
group1=['1bql', '1mlc', '1nby', '1nbz', '1ndg', '1ndm', '1p2c', '1ri8', '1rjc', '1yqv', '1zv5', '1zvh', '2iff', '3d9a']
group2=['1fdl', '1jto', '1jtp', '1jtt', '1xfp'] #removed 1zmy
group3a=['1a2y', '1g7h', '1g7i', '1g7j', '1g7l', '1g7m', '1kip', '1kiq', '1kir',  '1vfb']
group3b=['1sq2', '1t6v', '2hfm', '2i25', '2i26', '3eba', '4i0c']
group4=['1c08', '1ic5', '1ic7', '1ua6', '1uac', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqd', '2dqe', '2dqf', '2dqg', '2eiz', '2eks', '2yss', '3a6b', '3a6c', '3hfm']
group5=['1sq2', '1t6v',  '2i25', '2i26', '3eba', '4i0c','1c08', '1ic5', '1ic7', '1ua6', '1uac', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqd', '2dqe', '2dqf', '2dqg', '2eiz', '2eks', '2yss', '3a6b', '3a6c', '3hfm']
#removed 2hfm
contacts= json.loads(open('equiv_contacts_dict_lsyozyme_heavy.json').read())
species=json.loads(open('species_lysozime.json').read())
pdb_codes1= json.loads(open('pdb_codes_lysozyme.json').read())
germ_line=json.loads(open("germ_lines_lysozime.json").read())
def get_species(pdb):
    sp=''
    if "musculus" in species[pdb]:
        sp+='M'
    if "drom" in species[pdb]:
        sp+='C'
    elif "cirratum" in species[pdb]:
        sp+="S"
    elif "human" in species[pdb]:
        sp+="H"
    return sp
def get_json_dict(group, name):
    g1={}
    g11={}
    for k, v in contacts.items():
        for pdb in group:
            if pdb in k:
                g1[pdb]={}
                for a, b in v.items():
                    res1=''
                    resi=[]

                    for elem in b:
                        if elem not in resi:
                            res1+=elem+"-"
                            resi.append(elem)
                    g1[pdb][int(a)]=resi


    df2=pd.DataFrame.from_dict(g1,orient='index')
    #df2 = df2.replace(np.nan, "-")
    df2.index.rename('PDB', inplace=True)
    stacked = df2.stack().reset_index()

    stacked.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked=stacked.sort_values(by=['Position'])
    stacked2=stacked.stack().reset_index()

    def flatten_column(df, column_name):
        repeat_lens = [len(item) if item is not np.nan else 1 for item in df[column_name]]
        df_columns = list(df.columns)
        df_columns.remove(column_name)
        expanded_df = pd.DataFrame(np.repeat(df.drop(column_name, axis=1).values, repeat_lens, axis=0), columns=df_columns)
        flat_column_values = np.hstack(df[column_name].values)
        expanded_df[column_name] = flat_column_values
        expanded_df[column_name].replace('nan', np.nan, inplace=True)
        return expanded_df
    stackedd=flatten_column(stacked, 'Residues')

    common_pos={}
    for position, df_position in stackedd.groupby('Position'):
        ct=df_position['Residues'].value_counts()[:1].index.tolist()
        common_pos[position]=''
        for elem in ct:
            common_pos[position]+=elem

    with open('common_pos_'+name+'.json', 'w') as ctr: #save the contact residues in a json file
        json.dump(common_pos, ctr)



def get_heatmap_101(group, name):
    common_pos=json.loads(open("common_pos_"+name+".json").read())
    d_g1={}
    a_dict={}
    g1={}
    g11={}
    dif_pos={}
    my_list = [tuple()]
    for k, v in contacts.items():
        for pdb in group:
            d1={}
            if pdb in k:
                g1[pdb]={}
                g11[pdb]={}
                dif_pos[pdb]={}
                a_dict[pdb]={}

                for a, b in v.items():
                    res1=''
                    resi=[]
                    d1[a]=''
                    for elem in b:
                        if elem not in resi:
                            res1+=elem+"-"
                            resi.append(elem)
                    g1[pdb][int(a)]=resi
                    if len(g1[pdb][int(a)])>1:
                        if common_pos[a] in g1[pdb][int(a)]:
                            d1[a]=0

                        else:
                            d1[a]=-1
                            x=(common_pos[a],g1[pdb][int(a)])
                            my_list.append(x)

                    elif len(g1[pdb][int(a)])==1:
                        if common_pos[a] in g1[pdb][int(a)]:
                            d1[a]=1

                        else:
                            d1[a]=-1
                            x=(common_pos[a],g1[pdb][int(a)])
                            my_list.append(x)

                dif_pos[get_species(pdb)+":"+pdb]=d1
    #plot 2
    #f, ax = plt.subplots()
    df2=pd.DataFrame.from_dict(dif_pos,orient='index')
    for d in df2.columns:
        if d is not 'index':
            df2[d] = pd.to_numeric(df2[d])
    df2=df2[natsorted(df2.columns)]
    df2.groupby(by=species,axis=1)

    fig, ax = plt.subplots(figsize=(25,15))
    fig.autofmt_xdate()
    g = sns.heatmap(df2, xticklabels=True, yticklabels=True,  linewidths=.1,annot=False, annot_kws={"size": 20}, cmap="Paired", cbar=True,ax=ax)#.set_title(chain+" chain: "+param)

    g.get_figure().savefig("ep_aa_preference_res_"+name+"_hm.png")
    #------------------------------------Generate relplot-----------------
    df3=pd.DataFrame.from_dict(g1,orient='index')

    df3.index.rename('PDB', inplace=True)
    stacked = df3.stack().reset_index()
    stacked.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked=stacked.sort_values(by=['Position'])
    stacked2=stacked.stack().reset_index()

    #fig, ax = plt.subplots(figsize=(100,30))
    def flatten_column(df, column_name):
        repeat_lens = [len(item) if item is not np.nan else 1 for item in df[column_name]]
        df_columns = list(df.columns)
        df_columns.remove(column_name)
        expanded_df = pd.DataFrame(np.repeat(df.drop(column_name, axis=1).values, repeat_lens, axis=0), columns=df_columns)
        flat_column_values = np.hstack(df[column_name].values)
        expanded_df[column_name] = flat_column_values
        expanded_df[column_name].replace('nan', np.nan, inplace=True)
        return expanded_df
    stackedd=flatten_column(stacked, 'Residues')
"""
    fig, ax = plt.subplots(figsize=(12,7))
    fig.autofmt_xdate()
    b=sns.scatterplot(data=stackedd, x='Position', y='PDB', hue='Residues', ax=ax)
    sns.set_palette("RdBu")
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    b.set(yticks=[])
    b.get_figure().savefig("ep_aa_preference_res_"+name+".png")
"""

    #return my_list

def run_code(group, name):
    get_json_dict(group, name)
    get_heatmap_101(group, name)

run_code(pdb_codes1, "all_codes")

#=get_heatmap_101(group5, "group3+4_heavy")
dd=defaultdict(list)
listt=[('THR', ['SER']), ('SER', ['ARG', 'TYR', 'ALA']), ('TYR', ['ARG']), ('THR', ['TYR', 'ASP']), ('SER', ['GLY', 'TYR']), ('SER', ['VAL', 'TYR']), ('TRP', ['ALA']), ('SER', ['ARG', 'TYR', 'ALA']), ('TYR', ['ARG']), ('THR', ['TYR', 'ASP']), ('SER', ['GLY', 'TYR']), ('TYR', ['ARG']), ('SER', ['VAL', 'TYR']), ('SER', ['TYR', 'PHE']), ('TYR', ['PHE']), ('SER', ['TYR']), ('SER', ['TYR']), ('SER', ['TYR']), ('SER', ['ARG']), ('ARG', ['SER', 'TYR', 'ASP']), ('ARG', ['ASP']), ('ALA', ['ASP', 'ASN']), ('THR', ['GLU', 'TYR', 'ARG']), ('SER', ['ARG']), ('SER', ['VAL']), ('ARG', ['GLY', 'SER']), ('SER', ['ARG']), ('TYR', ['GLY']), ('SER', ['ARG']), ('SER', ['ALA']), ('ARG', ['SER', 'ASP']), ('ARG', ['ASN', 'ASP']), ('ARG', ['GLY', 'SER']), ('THR', ['TYR', 'GLU']), ('SER', ['ARG']), ('ASP', ['ARG']), ('TYR', ['SER']), ('THR', ['ARG']), ('SER', ['ARG']), ('SER', ['ARG', 'ASP']), ('TYR', ['ASP', 'GLY', 'ILE']), ('TRP', ['ALA', 'VAL', 'GLY']), ('TRP', ['VAL']), ('SER', ['GLY']), ('ASP', ['ARG', 'GLY', 'TRP']), ('TRP', ['GLU', 'VAL']), ('TRP', ['TYR'])]

for (a, b) in listt:

    dd[a].append(b)
dicti=dict(dd)
print(dicti)
"""
dict_pref={}
for a, b in dicti.items():
    #dict_pref[a]={}
    list_pref=[]
    list_2=[]
    for l1 in b:
        for elem in l1:
            list_pref.append(elem)
    dict_pref[a]=list_pref
ndf=pd.DataFrame.from_dict(dict_pref,orient='index')

ndf.index.rename('PDB', inplace=True)
stacked_n = ndf.stack().reset_index()
stacked_n.rename(columns={ 'level_1': 'Position',0: 'Residues'}, inplace=True)
#stacked_n=stacked_n.sort_values(by=['Position'])
stacked2_n=stacked_n.stack().reset_index()
common_pos={}
for position, df_position in stacked_n.groupby('PDB'):
    ct=df_position['Residues'].value_counts().index.tolist()
    common_pos[position]= ct
#print(common_pos)
cpd=pd.DataFrame.from_dict(common_pos,orient="index")
"""







