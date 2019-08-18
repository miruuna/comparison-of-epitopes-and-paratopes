import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from sortedcontainers import SortedDict
#import numbering as numbering
from collections import defaultdict
c2=['1ic5', '1ic7']
s1=['1a2y', '1fdl', '1g7h', '1g7i', '1g7j', '1g7l', '1g7m', '1kir', '1vfb']
s2=['1c08', '1j1o', '1j1p', '1j1x', '2dqi', '2dqj', '2znw', '3a67', '3a6b', '3a6c', '3d9a', '3hfm']
s3=['1dqj', '1nby', '1nbz','1ua6', '1uac']
s4=['2eiz', '2eks', '2yss']
s5=['1jto', '1jtp', '1jtt']
g1a=['1ri8', '1bql',  '1zv5', '1mlc', '2iff']
g1b=[ '1a2y', '1kiq', '1kip']
g1c=['1jto', '1xfp']
g22=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg']
g1_mice= g22
g2_mice= g1a
g3_cs=['1jto','1jtp', '1jtt','1rjc','1xfp','1zmy','1t6v','1sq2','2i25','2i26']
g1=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg']

c_and_m=['1jto', '1jtp', '1jtt', '1mel', '1op9', '1ri8', '1rjc', '1sq2', '1t6v', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '2i25', '2i26', '3eba', '4i0c']
m1=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg']
m2=['1bql','1mlc', '2iff']
cm=['1jto', '1jtp', '1jtt','1xfp','1zmy','2i25','2i26']
cm2=['1zvy','1sq2','1t6v']
cm3=['1ri8','1rjc','1zv5']
unique=['1sq2','1t6v','2hfm','2i25','2i26','1a2y', '1bql', '1c08', '1dqj', '1dzb', '1ic4', '1ic5', '1ic7', '1jhl', '1kip', '1kiq', '1mlc', '1ndg', '1ndm', '1p2c', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '2iff', '4tsa', '4tsb', '4tsc', '4ttd', '1jto', '1op9', '1ri8', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '4i0c']
unique1=['1a2y', '1c08', '1dqj', '1dzb', '1ic4', '1ic5', '1ic7', '1jhl', '1kip', '1kiq', '1mlc', '1ndg', '1ndm', '1p2c', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '2iff', '1jto', '1op9', '1ri8', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '4i0c']
#removed  1bql from unique
#removed 2hfm

unique_minus_1=set(unique).difference(set(g1))


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
def get_json_dict(group, name, chain, type):
    contacts =json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_'+type+'.json').read())
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
    contacts=json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_'+type+'.json').read())
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

                dif_pos[get_species(pdb)+"-"+pdb]=d1
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
    cmap=['#66c2a5','#fc8d62','#8da0cb']
    g = sns.heatmap(df2, xticklabels=True, yticklabels=True,  linewidths=0.1, linecolor='black', annot=False, annot_kws={"size": 20}, cmap=ListedColormap(cmap), cbar=True,ax=ax)#.set_title(chain+" chain: "+param)
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
    contacts1 =json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_'+type+'.json').read())


    cont_dict={}
    for k, v in contacts1.items():
        for pdb in group:
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
    the_contacts=json.loads(open("equiv_contacts_dict_lsyozyme_pos_paratope_"+chain+"_"+ type+".json").read())
    cont_dict={}
    for k, v in the_contacts.items():
        for pdb in group:
            d1={}
            if pdb in k:
                cont_dict[get_species(pdb)+"-"+pdb]={}
                for pos, names in v.items():
                    try:
                        if len(names)!=0:
                            dict={pos: len(list(set(names)))}
                            cont_dict[get_species(pdb)+"-"+pdb].update(dict)


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
    #colours=['#66c2a5','#386cb0','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']
    g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=.1, linecolor="black", annot=True, annot_kws={"size": 20,}, cmap="YlOrRd", cbar=False )#.set_title(chain+" chain: "+param)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=20)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax.set_aspect("equal")

    g.get_figure().savefig("new_epitope_no_contacts_"+name+"_"+chain+"_"+type+".png")
    print(cont_dict)

def get_heatmap_fw( group, name,chain, type, size, x):
    cont_paratopes=json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_'+type+'.json').read())
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
                                if int(paratope) in list(range(0,24)):
                                    dicto[pos].append("FW11")
                                elif int(paratope) in list(range(34,50)):
                                    dicto[pos].append("FW2")
                                elif int(paratope) in list(range(56,89)):
                                    dicto[pos].append("FW3")
                                elif int(paratope) >96:
                                    dicto[pos].append("FW4")
                            if chain=="heavy":
                                if int(paratope) in list(range(0,26)):
                                    dicto[pos].append("FW11")
                                elif int(paratope) in list(range(32,52)):
                                    dicto[pos].append("FW2")
                                elif int(paratope) in list(range(56,95)):
                                    dicto[pos].append("FW3")
                                elif int(paratope) >103:
                                    dicto[pos].append("FW4")


                        except KeyError:
                            continue
                        cont_dict[pdb].update(dicto)


    df3=pd.DataFrame.from_dict(cont_dict,orient='index')

    df3.index.rename('PDB', inplace=True)
    stacked = df3.stack().reset_index()
    stacked.rename(columns={'level_1': 'Position', 0: 'FWR'}, inplace=True)
    stacked=stacked.sort_values(by=['Position'])
    stacked2=stacked.stack().reset_index()

    c2=flatten_column(stacked, 'FWR')
    df3.to_excel(name+"_cs_"+type+".xlsx")
    return cont_dict

def get_heatmap_CDR( group, name,chain, type, size, x):
    cont_paratopes=json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_'+type+'.json').read())
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
    return cont_dict

def get_the_heatmap( group, name,chain, type, size, x):

    eq_dict_hm={}
    if chain=="heavy":
        chain_type="FW"
    elif chain=="light":
        chain_type="FW"
    for pdb, pos in get_heatmap_fw( group, name,chain, type, size, x).items():
        d_pos={}
        for p, cdr in pos.items():
            #eq_dict_hm[pdb][p]=[]
            if len(list(set(cdr)))>1:
                if chain_type+"1" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=4
                if chain_type+"1" and chain_type+"3" in list(set(cdr)):
                    d_pos[p]=5
                if chain_type+"3" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=6
                if chain_type+"3" and chain_type+"2" and chain_type+"1" in list(set(cdr)):
                    d_pos[p]=7
                #if chain_type+"1" and "o" in cdr:
                #   d_pos[p]=70
                #if chain_type+"2" and "o" in cdr:
                #    d_pos[p]=80
                #if "3" and "o" in cdr:
                #    d_pos[p]=90
            if len(list(set(cdr)))==1:
                if "o" in list(set(cdr)):
                    d_pos[p]=9
                if chain_type+"1" in list(set(cdr)):
                    d_pos[p]=1
                if chain_type+"2" in list(set(cdr)):
                    d_pos[p]=2
                if chain_type+"3"in cdr:
                    d_pos[p]=3
        eq_dict_hm[get_species(pdb)+"-"+pdb]=d_pos
    df =pd.DataFrame.from_dict(eq_dict_hm,orient='index')
    #df.set_index('pdb',inplace=True)
    df=df[natsorted(df.columns)]
    df.groupby(by=species,axis=1)
    df.index.rename('PDB', inplace=True)
    stacked = df.stack().reset_index()
    stacked.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked=stacked.sort_values(by=['Position'])
    maxim=int(stacked['Residues'].max())
    minim=int(stacked['Residues'].min())
    cmap=['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#bdbdbd']
    the_cmap=cmap[minim-1:maxim-1]
    fig, ax = plt.subplots(figsize=size)
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title

    g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=0.1, linecolor='black',annot_kws={"size": 30}, cmap=ListedColormap(the_cmap), cbar=True,cbar_kws={"shrink": .82})#.set_title(chain+" chain: "+param)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=30)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax.set_aspect("equal")
    #plt.show()
    g.get_figure().savefig("FWR_heatmap_"+chain+"_"+name+"_"+type+".png")


def get_chain_letter(chain):
    if chain=="heavy":
        letter="H"
    elif chain=="light":
        letter="L"
    return letter

def flipped_cdr_dict( group, name,chain, type, size, x):

    flipped = {}
    for k, v in get_heatmap_CDR( group, name,chain, type, size, x).items(): # flip the dictionary to get the epitope contacts which bind to each CDR
        flipped_pos={}
        for key, v in v.items():
            for value in v:
                if value not in flipped_pos:
                    flipped_pos[value] = [key]
                else:
                    flipped_pos[value].append(key)

        flipped[k]=flipped_pos
    sorted_l={}
    for k, v in flipped.items():
        sorted_l[k]={}
        for cdr_t, pos in v.items():
            l1=sorted(list(map(int, pos)))
            sorted_l[k][cdr_t]=list(map(str, l1))
    return sorted_l
def get_pos_freq_cdr( group, name,chain, type, size, x):
    by_cdr_dict_frq={}
    cdr_type_list=[]
    for cdr in [str(get_chain_letter(chain))+'1', str(get_chain_letter(chain))+'2',str(get_chain_letter(chain))+'3',"o"]:
        by_cdr_dict_frq[cdr]={}
        try:
            for pdb, cdr_type in flipped_cdr_dict( group, name,chain, type, size, x).items():
                by_cdr_dict_frq[cdr][pdb]={}
                for a_cdr, pos_d in cdr_type.items():
                    if a_cdr==cdr:
                        for pos in pos_d:
                            by_cdr_dict_frq[cdr][pdb][pos]=pos_d.count(pos)
        except KeyError:
            continue
        df =pd.DataFrame.from_dict(by_cdr_dict_frq[cdr],orient='index')

def sort_dict(d):
    sort_dic = {}

    for i in sorted(d):
        sort_dic.update({i:d[i]})
    return(sort_dic)

def get_aa_residues_cdr(group,name, chain, type, size, x, the_cdr):

    cdr_range={"H1": list(range(26,33)),
               "H2": list(range(52, 57)),
               "H3": list(range(95,103)),
               "L1": list(range(24,35)),
               "L2": list(range(50,57)),
               "L3": list(range(89,98)),
               "o": [0,1]}

    pos_res_paratopes=json.loads(open("equiv_contacts_dict_lsyozyme_pos_paratope_"+chain+"_"+type+".json").read())
    pos_name_paratope=json.loads(open("equiv_contacts_dict_lsyozyme_pos_name_paratope_"+chain+"_"+type+".json").read())

    pos_res_paratope={}
    for pdb in pdb_codes1:
        for k, v in pos_res_paratopes.items():
            if pdb in k:
                pos_res_paratope[pdb]=pos_res_paratopes[k]



    unique_epitope_dict={}
    epitope_dict={}
    for cdr in [str(get_chain_letter(chain))+'1', str(get_chain_letter(chain))+'2',str(get_chain_letter(chain))+'3',"o"]:
        unique_epitope_dict[cdr]={}
        for pdb, cdr_type in flipped_cdr_dict( group, name,chain, type, size, x).items():
            unique_epitope_dict[cdr][pdb]={}
            try:
                for a_cdr, pos_d in cdr_type.items():
                    if a_cdr==cdr:
                        for pos in pos_d:
                            inter_pos=pos_res_paratope[pdb][pos]
                            unique_epitope_dict[cdr][pdb][pos]=[]
                            for a_p in inter_pos:
                                if int(a_p) in cdr_range[cdr]:
                                    if pos_name_paratope[pdb][a_p][0] not in unique_epitope_dict[cdr][pdb][pos]:
                                        unique_epitope_dict[cdr][pdb][pos].append(pos_name_paratope[pdb][a_p][0])
            except KeyError:
                continue

    aa_dict={}
    for cdr, pdb_dict in unique_epitope_dict.items():
        aa_dict[cdr]={}
        for pdb, pos_dict in pdb_dict.items():
            aa_dict[cdr][pdb]={}
            for pos, aa_list in pos_dict.items():
                if len(aa_list)==2:
                    aa_dict[cdr][pdb][pos]=aa_list[0]
                    aa_dict[cdr][pdb][str(pos)+"a"]=aa_list[1]
                elif len(aa_list)==3:
                    aa_dict[cdr][pdb][pos]=aa_list[0]
                    aa_dict[cdr][pdb][str(pos)+"a"]=aa_list[1]
                    aa_dict[cdr][pdb][str(pos)+"b"]=aa_list[2]
                elif len(aa_list)==1:
                    aa_dict[cdr][pdb][pos]=aa_list[0]
                else:
                    aa_dict[cdr][pdb][pos]=[]
    return aa_dict
def get_cdr_on_epitope(group,name, chain, type, size, x, the_cdr):


    aa_type_dict={}
    hm_aa_type_dict={}
    for cdr, pdb_dict in get_aa_residues_cdr(group,name, chain, type, size, x, the_cdr).items():
        aa_type_dict[cdr]={}
        hm_aa_type_dict[cdr]={}
        for pdb, pos_dict in pdb_dict.items():
            aa_type_dict[cdr][get_species(pdb)+"-"+pdb]={}
            hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb]={}
            for pos, aa in pos_dict.items():
                if aa in ["ASP", "GLU"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="bright_red"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=1
                elif aa in ["CYS", "MET"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="yellow"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=2
                elif aa in ["LYS", "ARG"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="blue"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=3
                elif aa in ["SER", "THR"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="orange"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=4

                elif aa in ["PHE", "TYR"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="mild_blue"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=5
                elif aa in ["ASN", "GLN"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="cyan"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=6
                elif aa in ["Leu", "VAL", "ILE"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="green"
                elif aa in ["GLY"]:
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=7
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="light_gray"
                elif aa in ["ALA"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="dark_gray"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=8
                elif aa in ["TRP"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="pink"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=9
                elif aa in ["HIS"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="pale_blue"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=10
                elif aa in ["PRO"]:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="flesh"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=11
                else:
                    aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]="black"
                    hm_aa_type_dict[cdr][get_species(pdb)+"-"+pdb][pos]=12

    df_aa=pd.DataFrame.from_dict(aa_type_dict[the_cdr], orient="Index")
    df_aa.index.rename('PDB', inplace=True)
    stacked = df_aa.stack().reset_index()
    stacked.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked=stacked.sort_values(by=['Position'])
    stacked2=stacked.stack().reset_index()

    #b=sns.scatterplot(data=stacked, x='Position', y='PDB', hue='Residues')
    #plt.show()
    cmap2=[(230,10,10),(230,230,0),(20,90,255), (250,150,0),(50,50,170),(0,220,220),(15,130,15), (235,235,235),(200,200,200),(180,90,180),(130,130,210),(220,150,130)]
    cmap=['#ff0000', '#ffd700', '#6897bb', '#ff9a00', '#03396c', '#00e5e5', '#20b2aa', '#b6b4b9', '#635d5a', '#ff99aa', '#c9caff', '#f0c39e','#ffeffa']

    df_aa_hm=pd.DataFrame.from_dict(hm_aa_type_dict[the_cdr], orient="Index")
    df_aa=df_aa[natsorted(df_aa.columns)]
    stacked_hm = df_aa_hm.stack().reset_index()
    stacked_hm.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    #, cmap=ListedColormap(cmap[:z]
    df_aa_hm=df_aa_hm[natsorted(df_aa_hm.columns)]
    fig, ax = plt.subplots(figsize=size)
    fig.autofmt_xdate()
    maxim=int(stacked_hm['Residues'].max())
    minim=int(stacked_hm['Residues'].min())
    the_cmap=cmap[minim-1:maxim-1]
    a=sns.heatmap(df_aa_hm, linewidths=2, linecolor='black', cbar=True, cmap=ListedColormap(the_cmap), ax=ax)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=30)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax.set_aspect("equal")

    a.get_figure().savefig("aatype_"+chain+"_"+name+"_"+type+"_"+the_cdr+".png")

    #plt.show()
    print(hm_aa_type_dict)



def get_heatmap_aa_type_cdr(group,name, chain, type, size, x, the_cdr):
    hm_aa_type_dict=get_aa_residues_cdr(group,name, chain, type, size, x, the_cdr)

    df_aa=pd.DataFrame.from_dict(hm_aa_type_dict[the_cdr], orient="Index")
    df_aa.index.rename('PDB', inplace=True)
    stacked = df_aa.stack().reset_index()
    stacked.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked=stacked.sort_values(by=['Position'])
    stacked2=stacked.stack().reset_index()

    #b=sns.scatterplot(data=stacked, x='Position', y='PDB', hue='Residues')
    #plt.show()
    cmap2=[(230,10,10),(230,230,0),(20,90,255), (250,150,0),(50,50,170),(0,220,220),(15,130,15), (235,235,235),(200,200,200),(180,90,180),(130,130,210),(220,150,130)]
    cmap=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17']
    df_aa=df_aa[natsorted(df_aa.columns)]
    stacked_hm = df_aa.stack().reset_index()
    stacked_hm.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    counts =df_aa.apply(pd.Series.value_counts, axis=1)

    return counts
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
#        for cdr in [str(get_chain_letter(chain))+'1', str(get_chain_letter(chain))+'2',str(get_chain_letter(chain))+'3',"o"]:
for chain in ["light"]:
    for (a, b) in groups:
        x=5
        try:
            #get_json_dict(a, b)
            #get_cdr_on_epitope(a, b,chain, "all", (30,30)
            if  "g1" in b:
                #get_heatmap_no_cont(a, b, chain, "non_h_bonded", (30,30))
                get_the_heatmap(a, b, chain, "h_bonded",(30,30),x )
                get_the_heatmap(a, b, chain, "non_h_bonded",(30,30),x)
                get_the_heatmap(a, b, chain, "all",(30,30),x)
            if "g2" in b:
                #get_heatmap_no_cont(a, b, chain, "non_h_bonded", (45,9))
                get_the_heatmap(a, b, chain, "h_bonded",(30,30),x)
                get_the_heatmap(a, b, chain, "non_h_bonded",(30,30),x)
                get_the_heatmap(a, b, chain, "all",(30,30),x)
            if "g3" in b:
                #get_heatmap_no_cont(a, b, chain, "non_h_bonded", (70,20))
                get_the_heatmap(a, b, chain, "h_bonded",(30,30),x)
                get_the_heatmap(a, b, chain, "non_h_bonded",(30,30),x)
                get_the_heatmap(a, b, chain, "all",(30,30),x)
            if "unique" in b:
                #get_heatmap_no_cont(a, b, chain, "non_h_bonded", (100,40))
                get_the_heatmap(a, b, chain, "h_bonded",(30,30),x)
                get_the_heatmap(a, b, chain, "non_h_bonded",(30,30),x)
                get_the_heatmap(a, b, chain, "all",(30,30),x)


        except ValueError:
            continue


      for cdr in [str(get_chain_letter(chain))+'1', str(get_chain_letter(chain))+'2',str(get_chain_letter(chain))+'3',"o"]:
            try:
                get_the_heatmap(a, b, chain, "all", (30,30), 5)
            except IndexError:
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

groups= [(cm,"cm1"),(cm2,"cm2"),(cm3, "cm3"), (m1, "m1"),(m2,"m2")]
for chain in ["heavy"]:
    for cdr in [str(get_chain_letter(chain))+'1', str(get_chain_letter(chain))+'2',str(get_chain_letter(chain))+'3',"o"]:
        for (a,b) in groups:
            try:
                get_cdr_on_epitope(a, b, chain, "all",(20,20), 5, cdr)
            except IndexError:
                continue
"""               
with open('aa_type_frequency_heavy.json', 'w') as ctr: #save the contact residues in a json file
    json.dump(get_heatmap_aa_type_cdr(unique, "unique_all_pdbs", "heavy", "all", (5,5), 5, "H3"), ctr)
"""

