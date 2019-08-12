import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
#import numbering as numbering
from collections import defaultdict
c2=['1ic5', '1ic7']
s1=['1a2y', '1fdl', '1g7h', '1g7i', '1g7j', '1g7l', '1g7m', '1kir', '1vfb']
s2=['1c08', '1j1o', '1j1p', '1j1x', '2dqi', '2dqj', '2znw', '3a67', '3a6b', '3a6c', '3d9a', '3hfm']
s3=['1dqj', '1nby', '1nbz','1ua6', '1uac']
s4=['2eiz', '2eks', '2yss']
s5=['1jto', '1jtp', '1jtt']
g1a=['1ri8', '1bql', '1ndg', '1ndm', '1p2c', '1zv5', '1mlc', '2iff']
g1b=[ '1a2y', '1kiq', '1kip']
g1c=['1jto', '1xfp']
g2=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg']
g1_mice= g2
g2_mice= g1a
g3_cs=['1jto','1jtp', '1jtt','1rjc','1xfp','1zmy','1t6v','1sq2','2i25','2i26']
unique=['1a2y', '1bql', '1c08', '1dqj', '1dzb', '1ic4', '1ic5', '1ic7', '1jhl', '1kip', '1kiq', '1mlc', '1ndg', '1ndm', '1p2c', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '2iff', '1jto', '1op9', '1ri8', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '4i0c']

#removed 2hfm

pos_paratopes=json.loads(open("equiv_contacts_dict_lsyozyme_pos_by_bonded.json").read())#dict type pos epitope:pos paratope

contacts =json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_heavy_all.json').read())
species=json.loads(open('species_lysozime.json').read())
pdb_codes1= json.loads(open('pdb_codes_lysozyme.json').read())
germ_line=json.loads(open("germ_lines_lysozime.json").read())
pos_paratope=json.loads(open("equiv_contacts_dict_lsyozyme_pos_paratop.json").read())

def flatten_column(df, column_name):
    repeat_lens = [len(item) if item is not np.nan else 1 for item in df[column_name]]
    df_columns = list(df.columns)
    df_columns.remove(column_name)
    expanded_df = pd.DataFrame(np.repeat(df.drop(column_name, axis=1).values, repeat_lens, axis=0), columns=df_columns)
    flat_column_values = np.hstack(df[column_name].values)
    expanded_df[column_name] = flat_column_values
    expanded_df[column_name].replace('nan', np.nan, inplace=True)
    return expanded_df

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

    return stackedd

def get_heatmap_101(group, name, chain,type,size):
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

                dif_pos[pdb]=d1
    #plot 2
    #f, ax = plt.subplots()
    df2=pd.DataFrame.from_dict(dif_pos,orient='index')
    for d in df2.columns:
        if d is not 'index':
            df2[d] = pd.to_numeric(df2[d])
    df2=df2[natsorted(df2.columns)]
    df2.groupby(by=species,axis=1)

    fig, ax = plt.subplots(figsize=(size))
    fig.autofmt_xdate()
    g = sns.heatmap(df2, xticklabels=True, yticklabels=True,  linewidths=5,annot=False, annot_kws={"size": 20}, cmap="Paired", cbar=True,ax=ax)#.set_title(chain+" chain: "+param)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=30)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax.set_aspect("equal")
    g.get_figure().savefig("ep_aa_preference_res_hm_"+name+"_"+chain+"_"+type+".png")

    #------------------------------------Generate relplot-----------------
    df3=pd.DataFrame.from_dict(g1,orient='index')

    df3.index.rename('PDB', inplace=True)
    stacked = df3.stack().reset_index()
    stacked.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked=stacked.sort_values(by=['Position'])
    stacked2=stacked.stack().reset_index()

    #fig, ax = plt.subplots(figsize=(100,30))

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

def get_heatmap_no_cont_para(group, name, chain, type,  size):
    cont_paratopes=json.loads(open("equiv_contacts_dict_lsyozyme_pos_paratope"+chain+".json").read())
    pos={}
    x= cont_paratopes[str(chain)+":"+type]
    #y= epitope.get_contacts("heavy", "non_h_bonded")
    #z= dict(y, **x)
    a= cont_paratopes[str(chain)+":"+type]
    #b= epitope.get_contacts("light", "non_h_bonded")
    #c=dict(b, **a)
    contact_dict = dict(a, **x)



    cont_dict={}
    for k, v in contacts.items():
        for pdb in group:
            d1={}
            if pdb in k:
                cont_dict[pdb]={}
                for pos, names in v.items():
                    try:
                        if len(names)!=0:
                            dict={pos: len(names)}


                    except KeyError:
                        continue
                    cont_dict[pdb].update(dict)
        df =pd.DataFrame.from_dict(cont_dict,orient='index')
    #df.set_index('pdb',inplace=True)
    df=df[natsorted(df.columns)]
    df.groupby(by=species,axis=1)
    fig, ax = plt.subplots(figsize=(size))
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
    #plt.suptitle("Number of contacts - " + name, fontsize = 36) #title
    g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=1, linecolor='black',annot=True, annot_kws={"size": 30}, cmap="YlGnBu", cbar=True)#.set_title(chain+" chain: "+param)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=30)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax.set_aspect("equal")
    g.get_figure().savefig("paratope_no_contacts_"+name+"_"+type+".png")


def get_heatmap_no_cont(group, name,chain, type, size,x):
    #the_contacts=json.loads(open("equiv_contacts_dict_lsyozyme_pos_paratope_"+chain+"_"+ type+".json").read())
    cont_dict={}
    for k, v in contacts.items():
        for pdb in group:
            d1={}
            if pdb in k:
                cont_dict[pdb]={}
                for pos, names in v.items():
                    try:
                        if len(names)!=0:
                            dict={pos: len(list(set(names)))}
                            cont_dict[pdb].update(dict)


                    except KeyError:
                        continue

        df =pd.DataFrame.from_dict(cont_dict,orient='index')
    #df.set_index('pdb',inplace=True)
    df=df[natsorted(df.columns)]
    df.groupby(by=species,axis=1)
    fig, ax = plt.subplots(figsize=size)
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
    #plt.suptitle("Number of contacts - " + name, fontsize = 36) #title
    colours=['#66c2a5','#386cb0','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']
    g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=.1,annot=True, annot_kws={"size": 15,}, cmap=ListedColormap(colours[:x]), cbar=False )#.set_title(chain+" chain: "+param)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=15)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax.set_aspect("equal")

    g.get_figure().savefig("epitope_no_contacts_"+name+"_"+chain+"_"+type+".png")


def get_heatmap_CDR( group, name,chain, type, size, x):
    cont_paratopes=json.loads(open("equiv_contacts_dict_lsyozyme_pos_paratope_"+chain+"_"+type+".json").read())
    #x= cont_paratopes[str(chain)+":"+"h_bonded"]
    #y= epitope.get_contacts("heavy", "non_h_bonded")
    #z= dict(y, **x)
    #a= cont_paratopes[str(chain)+":"+"non_h_bonded"]
    #b= epitope.get_contacts("light", "non_h_bonded")
    #c=dict(b, **a)
    #contact_dict = dict(x, **a)
    #contacts_dict=dict(contact_dict, **a)

    cont_dict={}
    for k, v in cont_paratopes.items():
        for pdb in group:
            d1={}
            if pdb in k:
                cont_dict[pdb]={}
                dicto={}
                for pos, paratope_list in v.items():
                    arg_list=[]
                    dicto[pos]=[]
                    unique_p_list=[]

                    for paratope in paratope_list:
                        if paratope not in unique_p_list:
                            unique_p_list.append(paratope)
                    for paratope in unique_p_list:
                        try:
                            if chain =="light":
                                if int(paratope) in list(range(24,35)):
                                    dicto[pos].append("L1")
                                elif int(paratope) in list(range(50,57)):
                                    dicto[pos].append("L2")
                                elif int(paratope) in list(range(89,97)):
                                    dicto[pos].append("L3")
                                else:
                                    dicto[pos].append("O")
                            if chain=="heavy":
                                if int(paratope) in list(range(26,33)):
                                    dicto[pos].append("H1")
                                elif int(paratope) in list(range(52,57)):
                                    dicto[pos].append("H2")
                                elif int(paratope) in list(range(95,103)):
                                    dicto[pos].append("H3")
                                else:
                                    dicto[pos].append("o")

                        except KeyError:
                            continue
                        cont_dict[pdb].update(dicto)


    df3=pd.DataFrame.from_dict(cont_dict,orient='index')

    df3.index.rename('PDB', inplace=True)
    stacked = df3.stack().reset_index()
    stacked.rename(columns={'level_1': 'Position', 0: 'CDR'}, inplace=True)
    stacked=stacked.sort_values(by=['Position'])
    stacked2=stacked.stack().reset_index()

    c2=flatten_column(stacked, 'CDR')
    df3.to_excel("g3_cs_heavy.xlsx")

    eq_dict_hm={}
    if chain=="heavy":
        chain_type="H"
    elif chain=="light":
        chain_type="L"
    for pdb, pos in cont_dict.items():
        d_pos={}
        for p, cdr in pos.items():
            #eq_dict_hm[pdb][p]=[]
            if len(list(set(cdr)))>1:
                if chain_type+"1" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=40
                if chain_type+"1" and chain_type+"3" in list(set(cdr)):
                    d_pos[p]=50
                if chain_type+"3" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=60
                if chain_type+"3" and chain_type+"2" and chain_type+"1" in list(set(cdr)):
                    d_pos[p]=70
                #if chain_type+"1" and "o" in cdr:
                 #   d_pos[p]=70
                #if chain_type+"2" and "o" in cdr:
                #    d_pos[p]=80
                #if "3" and "o" in cdr:
                #    d_pos[p]=90
            if len(list(set(cdr)))==1:
                if "o" in list(set(cdr)):
                    d_pos[p]=90
                if chain_type+"1" in list(set(cdr)):
                    d_pos[p]=5
                if chain_type+"2" in list(set(cdr)):
                    d_pos[p]=20
                if chain_type+"3"in cdr:
                    d_pos[p]=30
        eq_dict_hm[pdb]=d_pos
    df =pd.DataFrame.from_dict(eq_dict_hm,orient='index')
    #df.set_index('pdb',inplace=True)
    df=df[natsorted(df.columns)]
    df.groupby(by=species,axis=1)
    fig, ax = plt.subplots(figsize=size)
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title

    g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=1, linecolor='black',annot=True,annot_kws={"size": 30}, cmap="Set2", robust=True, cbar=True,cbar_kws={"shrink": .82})#.set_title(chain+" chain: "+param)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=30)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax.set_aspect("equal")
    #plt.show()
    g.get_figure().savefig("CDR_heatmap_"+chain+"_"+name+"_"+type+".png")
    return cont_dict
def get_chain_letter(chain):
    if chain=="heavy":
        letter="H"
    elif chain=="light":
        letter="L"
    return letter

def get_cdr_on_epitope(group,name, chain, type, size, x):

    cdr_range={"H1": list(range(26,33)),
               "H2": list(range(52, 57)),
               "H3": list(range(95,103)),
               "L1": list(range(24,35)),
               "L2": list(range(50,57)),
               "L3": list(range(89,98)),
               "o": [0,1]}

    pos_res_paratopes=json.loads(open("equiv_contacts_dict_lsyozyme_pos_paratope_"+chain+"_all.json").read())
    pos_name_paratope=json.loads(open("equiv_contacts_dict_lsyozyme_pos_name_paratope_"+chain+".json").read())

    pos_res_paratope={}
    for pdb in pdb_codes1:
        for k, v in pos_res_paratopes.items():
            if pdb in k:
                pos_res_paratope[pdb]=pos_res_paratopes[k]
    flipped = {}
    for k, v in get_heatmap_CDR(group, name, chain, type, size, x).items(): # flip the dictionary to get the epitope contacts which bind to each CDR
        flipped_pos={}
        for key, v in v.items():
            for value in v:
                if value not in flipped_pos:
                    flipped_pos[value] = [key]
                else:
                    flipped_pos[value].append(key)
        flipped[k]=flipped_pos
    cdr_dict={}
    epitope_dict={}
    unique_epitope_dict={}

    for pdb,cdrs in flipped.items():
        cdr_dict={}
        unique_cdr_dict={}
        pos_dict={}
        unique_pos_dict={}
        for cdr_type, pos_cdr in cdrs.items():

            for a_pos in pos_cdr:
                pos_dict[a_pos]=[]
                unique_pos_dict[a_pos]=[]
                for a_p in pos_res_paratope[pdb][a_pos]:
                    if a_p in cdr_range[cdr_type]:

                        try:
                            pos_dict[a_pos].append(pos_name_paratope[pdb][a_p][0])
                            if pos_name_paratope[pdb][a_p][0] not in unique_pos_dict[a_pos]:
                                unique_pos_dict[a_pos].append(pos_name_paratope[pdb])


                        except KeyError:
                            continue
            cdr_dict[cdr_type]=pos_dict
            unique_cdr_dict[cdr_type]=unique_pos_dict

        unique_epitope_dict[pdb]=unique_cdr_dict
        epitope_dict[pdb]=cdr_dict

    by_cdr_dict_frq={}
    cdr_type_list=[]
    for cdr in [str(get_chain_letter(chain))+'1', str(get_chain_letter(chain))+'2',str(get_chain_letter(chain))+'3',"O"]:
        by_cdr_dict_frq[cdr]={}
        try:
            for pdb, cdr_type in flipped.items():
                by_cdr_dict_frq[cdr][pdb]={}
                for a_cdr, pos_d in cdr_type.items():
                    if a_cdr==cdr:
                        for pos in pos_d:
                            by_cdr_dict_frq[cdr][pdb][pos]=pos_d.count(pos)
        except KeyError:
            continue
        df =pd.DataFrame.from_dict(by_cdr_dict_frq[cdr],orient='index')

    return unique_pos_dict



def get_heatmap_hydrophobicity(group, name):
    cont_dict={}
    for k, v in contacts.items():
        for pdb in group:
            d1={}
            if pdb in k:
                cont_dict[pdb]={}
                for pos, names in v.items():
                    for res in names:
                        try:
                            dict={pos: len(names)}


                        except KeyError:
                            continue
                        cont_dict[pdb].update(dict)

"""
run_code(pdb_codes1, "all_codes")

#=get_heatmap_101(group5, "group3+4_heavy")
dd=defaultdict(list)
listt=[('THR', ['SER']), ('SER', ['ARG', 'TYR', 'ALA']), ('TYR', ['ARG']), ('THR', ['TYR', 'ASP']), ('SER', ['GLY', 'TYR']), ('SER', ['VAL', 'TYR']), ('TRP', ['ALA']), ('SER', ['ARG', 'TYR', 'ALA']), ('TYR', ['ARG']), ('THR', ['TYR', 'ASP']), ('SER', ['GLY', 'TYR']), ('TYR', ['ARG']), ('SER', ['VAL', 'TYR']), ('SER', ['TYR', 'PHE']), ('TYR', ['PHE']), ('SER', ['TYR']), ('SER', ['TYR']), ('SER', ['TYR']), ('SER', ['ARG']), ('ARG', ['SER', 'TYR', 'ASP']), ('ARG', ['ASP']), ('ALA', ['ASP', 'ASN']), ('THR', ['GLU', 'TYR', 'ARG']), ('SER', ['ARG']), ('SER', ['VAL']), ('ARG', ['GLY', 'SER']), ('SER', ['ARG']), ('TYR', ['GLY']), ('SER', ['ARG']), ('SER', ['ALA']), ('ARG', ['SER', 'ASP']), ('ARG', ['ASN', 'ASP']), ('ARG', ['GLY', 'SER']), ('THR', ['TYR', 'GLU']), ('SER', ['ARG']), ('ASP', ['ARG']), ('TYR', ['SER']), ('THR', ['ARG']), ('SER', ['ARG']), ('SER', ['ARG', 'ASP']), ('TYR', ['ASP', 'GLY', 'ILE']), ('TRP', ['ALA', 'VAL', 'GLY']), ('TRP', ['VAL']), ('SER', ['GLY']), ('ASP', ['ARG', 'GLY', 'TRP']), ('TRP', ['GLU', 'VAL']), ('TRP', ['TYR'])]

for (a, b) in listt:

    dd[a].append(b)
dicti=dict(dd)
print(dicti)

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
"""
groups= [(g1_mice, "g1_mice"), (g2_mice, "g2_mice"), (g3_cs, "g3_cs"), (unique, "unique_all_pdbs")]
for chain in ["light"]:
    for (a, b) in groups:
        x=5
        try:
            #get_json_dict(a, b)
            #get_cdr_on_epitope(a, b,chain, "all", (30,30)
            if  "g1" in b:
                #get_heatmap_no_cont(a, b, chain, "non_h_bonded", (30,30))
                get_heatmap_CDR(a, b, chain, "h_bonded",(30,30),x )
                get_heatmap_CDR(a, b, chain, "non_h_bonded",(30,30),x)
                get_heatmap_CDR(a, b, chain, "all",(30,30),x)
            if "g2" in b:
                #get_heatmap_no_cont(a, b, chain, "non_h_bonded", (45,9))
                get_heatmap_CDR(a, b, chain, "h_bonded",(30,30),x)
                get_heatmap_CDR(a, b, chain, "non_h_bonded",(30,30),x)
                get_heatmap_CDR(a, b, chain, "all",(30,30),x)
            if "g3" in b:
                #get_heatmap_no_cont(a, b, chain, "non_h_bonded", (70,20))
                get_heatmap_CDR(a, b, chain, "h_bonded",(30,30),x)
                get_heatmap_CDR(a, b, chain, "non_h_bonded",(30,30),x)
                get_heatmap_CDR(a, b, chain, "all",(30,30),x)
            if "unique" in b:
                #get_heatmap_no_cont(a, b, chain, "non_h_bonded", (100,40))
                get_heatmap_CDR(a, b, chain, "h_bonded",(30,30),x)
                get_heatmap_CDR(a, b, chain, "non_h_bonded",(30,30),x)
                get_heatmap_CDR(a, b, chain, "all",(30,30),x)


        except ValueError:
            continue

#get_heatmap_no_cont(unique, "unique", "all", "all", (30,30))

#g3-45,15
#g2-45,10
#g1- 30, 15

#h-bonded
#g1- 15,30
#g2 -30,10
#g3 - 30,18

"""




print(get_cdr_on_epitope(g1_mice, "g1", "heavy", "all", (30,30), 5))
