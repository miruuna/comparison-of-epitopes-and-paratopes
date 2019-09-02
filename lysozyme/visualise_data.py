
#This script contains the functions to plot different heatmaps 



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


c_s=['1jto', '1jtp', '1jtt','1xfp','1zmy','2i25','2i26','1zvy','1sq2','1t6v','1ri8','1rjc','1zv5','c']
m1=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg','c']
m2=['1bql','1mlc', '2iff','c']
cm=['1jto', '1jtp', '1jtt','1xfp','1zmy','2i25','2i26','c']
cm2=['1zvy','1sq2','1t6v','c']
cm3=['1ri8','1rjc','1zv5','c']
unique=['1sq2','1t6v','2hfm','2i25','2i26','1a2y', '1bql', '1c08', '1dqj', '1dzb', '1ic4', '1ic5', '1ic7', '1jhl', '1kip', '1kiq', '1mlc', '1ndg', '1ndm', '1p2c', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '2iff', '4tsa', '4tsb', '4tsc', '4ttd', '1jto', '1op9', '1ri8', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '4i0c']
unique1=['1a2y', '1c08', '1dqj', '1dzb', '1ic4', '1ic5', '1ic7', '1jhl', '1kip', '1kiq', '1mlc', '1ndg', '1ndm', '1p2c', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '2iff', '1jto', '1op9', '1ri8', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '4i0c']
#removed  1bql from unique
#removed 2hfm
all_groups=m1+m2+cm+cm2+cm3

mice=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg','1bql','1mlc', '2iff']
camel=['1jto', '1jtp', '1jtt','1xfp','1zmy','1zvy','1ri8','1rjc','1zv5']
shark=['1sq2','1t6v','2i25','2i26']

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
    #contacts1 =json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_'+type+'.json').read())
    contacts_h_bonded=json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_h_bonded.json').read())
    contacts_all =json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_all.json').read())


    cont_dict_h_bonded={}
    for k, v in contacts_h_bonded.items():
        for pdb in group:
            if pdb in k:
                cont_dict_h_bonded[get_species(pdb)+"-"+pdb]={}
                for pos, names in v.items():
                    try:
                        if len(names)!=0:
                            dict={pos: len(names)}


                    except KeyError:
                        continue
                    cont_dict_h_bonded[get_species(pdb)+"-"+pdb].update(dict)

    cont_dict_all={}
    for k, v in contacts_all.items():
        for pdb in group:
            if pdb in k:
                cont_dict_all[get_species(pdb)+"-"+pdb]={}
                for pos, names in v.items():
                    try:
                        if len(names)!=0:
                            dict={pos: len(names)}


                    except KeyError:
                        continue
                    cont_dict_all[get_species(pdb)+"-"+pdb].update(dict)
    cont_h_bonded=cont_dict_h_bonded
    cont_all=cont_dict_all
    list_pos_all=[]
    for y, b in cont_dict_all.items():
        for pos2, v2 in b.items():
            list_pos_all.append(pos2)
    list_pos_h=[]
    for x, a in cont_dict_h_bonded.items():
        for pos1, v1 in a.items():
            list_pos_h.append(pos1)
    for db in group:
        for pos in list_pos_all:
            for posh in list_pos_h:
                if pos not in posh:
                    cont_h_bonded[get_species(pdb)+"-"+pdb][pos]=0
                if posh !=pos:
                    cont_all[get_species(pdb)+"-"+pdb][posh]=0
    df_all =pd.DataFrame.from_dict(cont_all,orient='index')
    #df.set_index('pdb',inplace=True)
    df_all=df_all[natsorted(df_all.columns)]
    df_all.groupby(by=species,axis=1)


    df_h =pd.DataFrame.from_dict(cont_h_bonded,orient='index')
    #df.set_index('pdb',inplace=True)
    df_h=df_h[natsorted(df_h.columns)]
    df_h.groupby(by=species,axis=1)



    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(size))
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
    #plt.suptitle("Number of contacts - " + name, fontsize = 36) #title
    sns.heatmap(df_all, ax=ax1,xticklabels=True, yticklabels=True, linewidths=0.1, linecolor='black',annot=True, annot_kws={"size": 8}, cmap="YlOrRd", cbar=True)#.set_title("All contacts",y=1.1)
    ax1.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    ax1.set_yticklabels(ax1.get_yticklabels(),rotation=0)
    ax1.set_xticklabels(ax1.get_xticklabels(),rotation=90)
    ax1.set_aspect("equal")
    sns.heatmap(df_h, ax=ax2,xticklabels=True, yticklabels=True, linewidths=0.1, linecolor='black',annot=True, annot_kws={"size": 8}, cmap="YlOrRd", cbar=True)#.set_title("H bonded",y=1.1)
    ax2.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    ax2.set_yticklabels(ax2.get_yticklabels(),rotation=0)
    ax2.set_xticklabels(ax2.get_xticklabels(),rotation=90)
    ax2.set_aspect("equal")
    plt.show()

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
                                    dicto[pos].append("FW1")
                                elif int(paratope) in list(range(35,50)):
                                    dicto[pos].append("FW2")
                                elif int(paratope) in list(range(57,89)):
                                    dicto[pos].append("FW3")
                                elif int(paratope) >96:
                                    dicto[pos].append("FW4")
                            if chain=="heavy":
                                if int(paratope) in list(range(0,31)):
                                    dicto[pos].append("FW1")
                                elif int(paratope) in list(range(36,50)):
                                    dicto[pos].append("FW2")
                                elif int(paratope) in list(range(66,95)):
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
    pos_antibody={}
    for k, v in cont_paratopes.items():
        for pdb in group:
            d1={}
            if pdb in k:
                cont_dict[pdb]={}
                pos_antibody[pdb]={}
                dicto={}
                ant_list={}
                for pos, paratope_list in v.items():

                    dicto[pos]=[]
                    ant_list[pos]=[]
                    unique_p_list=[]

                    for paratope in paratope_list:
                        if paratope not in unique_p_list:
                            unique_p_list.append(paratope)
                    for paratope in unique_p_list:
                        try:
                            if chain =="light":
                                if int(paratope) in list(range(24,35)):
                                    dicto[pos].append("L1")
                                    ant_list[pos].append(int(paratope))
                                elif int(paratope) in list(range(50,57)):
                                    dicto[pos].append("L2")

                                elif int(paratope) in list(range(89,98)):
                                    dicto[pos].append("L3")

                                else:
                                    dicto[pos].append("O")
                            if chain=="heavy":
                                if int(paratope) in list(range(31,36)):
                                    dicto[pos].append("H1")

                                elif int(paratope) in list(range(50,66)):
                                    dicto[pos].append("H2")

                                elif int(paratope) in list(range(95,103)):
                                    dicto[pos].append("H3")

                                else:
                                    dicto[pos].append("o")
                                    ant_list[pos].append(int(paratope))

                        except KeyError:
                            continue
                        cont_dict[pdb].update(dicto)
                        pos_antibody[pdb].update(ant_list)


    df3=pd.DataFrame.from_dict(cont_dict,orient='index')

    df3.index.rename('PDB', inplace=True)
    stacked = df3.stack().reset_index()
    stacked.rename(columns={'level_1': 'Position', 0: 'CDR'}, inplace=True)
    stacked=stacked.sort_values(by=['Position'])
    stacked2=stacked.stack().reset_index()

    c2=flatten_column(stacked, 'CDR')
    df3.to_excel("g3_cs_heavy.xlsx")
    #return cont_dict
    the_list={}
    for k, v in pos_antibody.items():
        the_list[k]=[]
        for pos, e in v.items():
            if len(e)>0:
                for elem in e:
                    the_list[k].append(elem)
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
                    d_pos[p]=8
                if chain_type+"1" in list(set(cdr)):
                    d_pos[p]=1
                if chain_type+"2" in list(set(cdr)):
                    d_pos[p]=2
                if chain_type+"3"in cdr:
                    d_pos[p]=3
                if chain_type+"4"in cdr:
                    d_pos[p]=4
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
    the_cmap=cmap[minim-1:maxim ]
    fig, ax = plt.subplots(figsize=size)
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title

    g = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax, linewidths=0.1, linecolor='black',annot_kws={"size": 30}, cmap=ListedColormap(the_cmap), cbar=True,cbar_kws={"shrink": .82})#.set_title(chain+" chain: "+param)
    ax.tick_params(labelbottom='off',labeltop='on', labelsize=30)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax.set_aspect("equal")
    plt.show()
   # g.get_figure().savefig("FWR_heatmap_"+chain+"_"+name+"_"+type+".png")

def get_heatmap_cdr_bonding( group, name,chain,type,  size, x):
    no_cdr_b={}
    cdr_by_species={}
    eq_dict_all={}
    if chain=="heavy":
        chain_type="H"
    elif chain=="light":
        chain_type="L"
    cdr_by_species[chain_type+"3"]=0
    cdr_by_species[chain_type+"1"]=0
    cdr_by_species[chain_type+"2"]=0
    for pdb, pos in get_heatmap_CDR( group, name,chain, type, size, x).items():
        d_pos={}
        no_cdr_b[get_species(pdb)+"-"+pdb]={}
        for p, cdr in pos.items():
            #eq_dict_hm[pdb][p]=[]
            no_cdr_b[get_species(pdb)+"-"+pdb][chain_type+"1"]=0
            no_cdr_b[get_species(pdb)+"-"+pdb][chain_type+"2"]=0
            no_cdr_b[get_species(pdb)+"-"+pdb][chain_type+"3"]=0
            if len(list(set(cdr)))>1:
                if chain_type+"1" in list(set(cdr)):
                    no_cdr_b[get_species(pdb)+"-"+pdb][chain_type+"1"]+=1
                    cdr_by_species[chain_type+"1"]+=1
                if chain_type+"2" in list(set(cdr)):
                    no_cdr_b[get_species(pdb)+"-"+pdb][chain_type+"2"]+=1
                    cdr_by_species[chain_type+"2"]+=1
                if chain_type+"3" in list(set(cdr)):
                    no_cdr_b[get_species(pdb)+"-"+pdb][chain_type+"3"]+=1
                    cdr_by_species[chain_type+"3"]+=1

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
                    d_pos[p]=8
                if chain_type+"1" in list(set(cdr)):
                    d_pos[p]=1
                    no_cdr_b[get_species(pdb)+"-"+pdb][chain_type+"1"]+=1
                    cdr_by_species[chain_type+"1"]+=1
                if chain_type+"2" in list(set(cdr)):
                    no_cdr_b[get_species(pdb)+"-"+pdb][chain_type+"2"]+=1
                    d_pos[p]=2
                    cdr_by_species[chain_type+"2"]+=1
                if chain_type+"3"in cdr:
                    d_pos[p]=3
                    no_cdr_b[get_species(pdb)+"-"+pdb][chain_type+"3"]+=1
                    cdr_by_species[chain_type+"3"]+=1

        eq_dict_all[pdb]=d_pos
    eq_dict_h_bonded={}
    for pdb, pos in get_heatmap_CDR( group, name,chain, "h_bonded", size, x).items():
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
                    d_pos[p]=8
                if chain_type+"1" in list(set(cdr)):
                    d_pos[p]=1
                if chain_type+"2" in list(set(cdr)):
                    d_pos[p]=2
                if chain_type+"3"in cdr:
                    d_pos[p]=3
        eq_dict_h_bonded[get_species(pdb)+"-"+pdb]=d_pos

    all_list=[]
    h_list=[]
    for k, v in eq_dict_all.items():
        for pos1, v1 in v.items():
            if pos1 not in all_list:
                all_list.append(pos1)
    for k, v in eq_dict_h_bonded.items():
        for pos1, v1 in v.items():
            if pos1 not in h_list:
                h_list.append(pos1)
    eq_dict_h=eq_dict_h_bonded
    for pdb in group:
        for pos_h in h_list:
            for pos_a in all_list:
                if pos_a not in h_list:
                    eq_dict_h_bonded[get_species(pdb)+"-"+pdb][pos_a]=9

    df_all =pd.DataFrame.from_dict(eq_dict_all,orient='index')
    #df.set_index('pdb',inplace=True)
    df_all=df_all[natsorted(df_all.columns)]
    df_all.groupby(by=species,axis=1)
    df_all.index.rename('PDB', inplace=True)
    df_h =pd.DataFrame.from_dict(eq_dict_h_bonded,orient='index')
    #df.set_index('pdb',inplace=True)
    df_h=df_h[natsorted(df_h.columns)]
    df_h.groupby(by=species,axis=1)
    df_h.index.rename('PDB', inplace=True)
    stacked_h = df_h.stack().reset_index()
    stacked_h.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked_h=stacked_h.sort_values(by=['Position'])
    maxim_h=int(stacked_h['Residues'].max())
    minim_h=int(stacked_h['Residues'].min())
    cmap_h=['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#bdbdbd','#FFFAFA','#faed27','#ffffff','#faed27']
    the_cmap_h=cmap_h[minim_h-1:maxim_h]
    stacked_all = df_all.stack().reset_index()
    stacked_all.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked_all=stacked_all.sort_values(by=['Position'])
    maxim_all=int(stacked_all['Residues'].max())
    minim_all=int(stacked_all['Residues'].min())
    cmap_all=['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#bdbdbd','#FFFAFA','#faed27','#ffffff']
    the_cmap_all=cmap_h[minim_all-1:maxim_all]

    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(size))
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
    #plt.suptitle("Number of contacts - " + name, fontsize = 36) #title
    sns.heatmap(df_all, xticklabels=True, yticklabels=True, ax=ax1, linewidths=0.1, linecolor='black',annot=False,annot_kws={"size": 10}, cmap=ListedColormap(the_cmap_all), cbar=True,cbar_kws={"shrink": .82})#.set_title("All contacts", y=1.12)
    ax1.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    ax1.set_yticklabels(ax1.get_yticklabels(),rotation=0)
    ax1.set_xticklabels(ax1.get_xticklabels(),rotation=90)
    ax1.set_aspect("equal")
    sns.heatmap(df_h, xticklabels=True, yticklabels=True, ax=ax2, linewidths=0.1, linecolor='black',annot=False,annot_kws={"size": 10}, cmap=ListedColormap(the_cmap_h), cbar=True,cbar_kws={"shrink": .82})#.set_title("H-bonded contacts", y=1.12)
    ax2.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    ax2.set_yticklabels(ax2.get_yticklabels(),rotation=0)
    ax2.set_xticklabels(ax2.get_xticklabels(),rotation=90)
    ax2.set_aspect("equal")
    #plt.show()



    print(cdr_by_species)
    df=pd.DataFrame.from_dict(cdr_by_species, orient="index")
    df.index.rename('CDR', inplace=True)
    df.rename(columns={ 0: type}, inplace=True)
    #with open("cdr_freq_"+name+"_"+chain+"_"+type+".json", 'w') as ctr: #save the contact residues in a json file
        #json.dump(cdr_by_species, ctr)
    return(eq_dict_all)


def get_the_heatmap_fw(group, name,chain,  size, x):

    eq_dict_all={}
    if chain=="heavy":
        chain_type="FW"
    elif chain=="light":
        chain_type="FW"
    for pdb, pos in get_heatmap_fw( group, name,chain, "all", size, x).items():
        d_pos={}
        for p, cdr in pos.items():
            #eq_dict_hm[pdb][p]=[]
            if len(list(set(cdr)))>1:
                if chain_type+"1" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=7
                if chain_type+"1" and chain_type+"3" in list(set(cdr)):
                    d_pos[p]=5
                if chain_type+"3" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=6
                if chain_type+"3" and chain_type+"2" and chain_type+"1" in list(set(cdr)):
                    d_pos[p]=8
                if chain_type+"1" and chain_type+"4" in list(set(cdr)):
                    d_pos[p]=9
                if chain_type+"4" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=10
                if chain_type+"3" and chain_type+"4" in list(set(cdr)):
                    d_pos[p]=11
                if chain_type+"3" and chain_type+"2" and chain_type+"1" and chain_type+"4" in list(set(cdr)):
                    d_pos[p]=12

                #if chain_type+"1" and "o" in cdr:
                #   d_pos[p]=70
                #if chain_type+"2" and "o" in cdr:
                #    d_pos[p]=80
                #if "3" and "o" in cdr:
                #    d_pos[p]=90
            if len(list(set(cdr)))==1:
                if "o" in list(set(cdr)):
                    d_pos[p]=8
                if chain_type+"1" in list(set(cdr)):
                    d_pos[p]=1
                if chain_type+"2" in list(set(cdr)):
                    d_pos[p]=2
                if chain_type+"3"in cdr:
                    d_pos[p]=3
                if chain_type+"4"in cdr:
                    d_pos[p]=4
        eq_dict_all[get_species(pdb)+"-"+pdb]=d_pos
    eq_dict_h_bonded={}
    for pdb, pos in get_heatmap_fw( group, name,chain, "h_bonded", size, x).items():
        d_pos={}
        for p, cdr in pos.items():
            #eq_dict_hm[pdb][p]=[]
            if len(list(set(cdr)))>1:
                if chain_type+"1" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=7
                if chain_type+"1" and chain_type+"3" in list(set(cdr)):
                    d_pos[p]=5
                if chain_type+"3" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=6
                if chain_type+"3" and chain_type+"2" and chain_type+"1" in list(set(cdr)):
                    d_pos[p]=8
                if chain_type+"1" and chain_type+"4" in list(set(cdr)):
                    d_pos[p]=9
                if chain_type+"4" and chain_type+"2" in list(set(cdr)):
                    d_pos[p]=10
                if chain_type+"3" and chain_type+"4" in list(set(cdr)):
                    d_pos[p]=11
                if chain_type+"3" and chain_type+"2" and chain_type+"1" and chain_type+"4" in list(set(cdr)):
                    d_pos[p]=12
            if len(list(set(cdr)))==1:
                if "o" in list(set(cdr)):
                    d_pos[p]=8
                if chain_type+"1" in list(set(cdr)):
                    d_pos[p]=1
                if chain_type+"2" in list(set(cdr)):
                    d_pos[p]=2
                if chain_type+"3"in cdr:
                    d_pos[p]=3
        eq_dict_h_bonded[get_species(pdb)+"-"+pdb]=d_pos

    all_list=[]
    h_list=[]
    for k, v in eq_dict_all.items():
        for pos1, v1 in v.items():
            if pos1 not in all_list:
                all_list.append(pos1)
    for k, v in eq_dict_h_bonded.items():
        for pos1, v1 in v.items():
            if pos1 not in h_list:
                h_list.append(pos1)
    eq_dict_h=eq_dict_h_bonded
    for pdb in group:
        for pos_h in h_list:
            for pos_a in all_list:
                if pos_a not in h_list:
                    eq_dict_h_bonded[get_species(pdb)+"-"+pdb][pos_a]=13

    df_all =pd.DataFrame.from_dict(eq_dict_all,orient='index')
    #df.set_index('pdb',inplace=True)
    df_all=df_all[natsorted(df_all.columns)]
    df_all.groupby(by=species,axis=1)
    df_all.index.rename('PDB', inplace=True)
    df_h =pd.DataFrame.from_dict(eq_dict_h_bonded,orient='index')
    #df.set_index('pdb',inplace=True)
    df_h=df_h[natsorted(df_h.columns)]
    df_h.groupby(by=species,axis=1)
    df_h.index.rename('PDB', inplace=True)
    stacked_h = df_h.stack().reset_index()
    stacked_h.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked_h=stacked_h.sort_values(by=['Position'])
    maxim_h=int(stacked_h['Residues'].max())
    minim_h=int(stacked_h['Residues'].min())
    cmap_h=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#ffffff']
    the_cmap_h=cmap_h[minim_h-1:maxim_h]
    stacked_all = df_all.stack().reset_index()
    stacked_all.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked_all=stacked_all.sort_values(by=['Position'])
    maxim_all=int(stacked_all['Residues'].max())
    minim_all=int(stacked_all['Residues'].min())
    cmap_all=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#ffffff']
    the_cmap_all=cmap_h[minim_all-1:maxim_all]

    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(size))
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
    #plt.suptitle("Number of contacts - " + name, fontsize = 36) #title
    sns.heatmap(df_all, xticklabels=True, yticklabels=True, ax=ax1, linewidths=0.1, linecolor='black',annot=False,annot_kws={"size": 10}, cmap=ListedColormap(the_cmap_all), cbar=True,cbar_kws={"shrink": .82})#.set_title("All contacts", y=1.12)
    ax1.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    ax1.set_yticklabels(ax1.get_yticklabels(),rotation=0)
    ax1.set_xticklabels(ax1.get_xticklabels(),rotation=90)
    ax1.set_aspect("equal")
    sns.heatmap(df_h, xticklabels=True, yticklabels=True, ax=ax2, linewidths=0.1, linecolor='black',annot=False,annot_kws={"size": 10}, cmap=ListedColormap(the_cmap_h), cbar=True,cbar_kws={"shrink": .82})#.set_title("H-bonded contacts", y=1.12)
    ax2.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    ax2.set_yticklabels(ax2.get_yticklabels(),rotation=0)
    ax2.set_xticklabels(ax2.get_xticklabels(),rotation=90)
    ax2.set_aspect("equal")
    plt.show()


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
                elif aa in ["LEU", "VAL", "ILE"]:
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
def get_heatmap_ep_shapely(group, name, chain, type,  size):
    #contacts1 =json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_'+type+'.json').read())
    contacts_h_bonded=json.loads(open('equiv_contacts_dict_lsyozyme_pos_name_paratope_'+chain+'_h_bonded.json').read())
    contacts_all =json.loads(open('equiv_contacts_dict_lsyozyme_pos_name_paratope_'+chain+'_all.json').read())
    aa_type={"1": ["ASP", "GLU"],
             "2":["CYS","MET"],
             "3":["LYS","ARG"],
             "4":["SER","THR"],
             "5":["PHE","TYR"],
             "6":["ASN","GLN"],
             "7":["LEU","VA","ILE"],
             "8":["GLY"],
             "9":["ALA"],
             "10":["TRP"],
             "11":["HIS"],
             "12":["PRO"]}

    cont_dict_h_bonded={}
    for k, v in contacts_h_bonded.items():
        for pdb in group:
            if pdb in k:
                cont_dict_h_bonded[get_species(pdb)+"-"+pdb]={}
                for pos, names in v.items():
                    try:
                        if all(aa in ["ASP", "GLU"]for aa in names) is True:

                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=1
                        elif all(aa in ["CYS", "MET"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=2
                        elif all(aa in ["LYS", "ARG"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=3

                        elif all(aa in ["SER", "THR"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=4

                        elif all(aa in ["PHE", "TYR"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=5
                        elif all(aa in ["ASN", "GLN"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=6
                        elif all(aa in ["LEU", "VAL", "ILE"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=7
                        elif all(aa in ["GLY"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=8
                        elif all(aa in ["ALA"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=9
                        elif all(aa in ["TRP"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=10
                        elif all(aa in ["HIS"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=11
                        elif all(aa in ["PRO"] for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=12
                        else:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=13

                    except KeyError:
                        continue
    cont_dict_all={}
    for k, v in contacts_all.items():
        for pdb in group:
            if pdb in k:
                cont_dict_all[get_species(pdb)+"-"+pdb]={}
                for pos, names in v.items():
                    try:
                        if all(aa in ["ASP", "GLU"]for aa in names) is True:

                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=1
                        elif all(aa in ["CYS", "MET"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=2
                        elif all(aa in ["LYS", "ARG"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=3

                        elif all(aa in ["SER", "THR"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=4

                        elif all(aa in ["PHE", "TYR"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=5
                        elif all(aa in ["ASN", "GLN"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=6
                        elif all(aa in ["LEU", "VAL", "ILE"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=7
                        elif all(aa in ["GLY"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=8
                        elif all(aa in ["ALA"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=9
                        elif all(aa in ["TRP"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=10
                        elif all(aa in ["HIS"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=11
                        elif all(aa in ["PRO"] for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=12
                        else:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=13

                    except KeyError:
                        continue
    cont_h_bonded=cont_dict_h_bonded
    cont_all=cont_dict_all
    list_pos_all=[]
    for y, b in cont_dict_all.items():
        for pos2, v2 in b.items():
            list_pos_all.append(pos2)
    list_pos_h=[]
    for x, a in cont_dict_h_bonded.items():
        for pos1, v1 in a.items():
            list_pos_h.append(pos1)
    for db in group:
        for pos in list_pos_all:
            for posh in list_pos_h:
                if pos not in posh:
                    cont_h_bonded[get_species(pdb)+"-"+pdb][pos]=14
                if posh !=pos:
                    cont_all[get_species(pdb)+"-"+pdb][posh]=14
    df_all =pd.DataFrame.from_dict(cont_all,orient='index')
    #df.set_index('pdb',inplace=True)
    df_all=df_all[natsorted(df_all.columns)]
    df_all.groupby(by=species,axis=1)


    df_h =pd.DataFrame.from_dict(cont_h_bonded,orient='index')
    #df.set_index('pdb',inplace=True)
    df_h=df_h[natsorted(df_h.columns)]
    df_h.groupby(by=species,axis=1)
    stacked_h = df_h.stack().reset_index()
    stacked_h.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked_h=stacked_h.sort_values(by=['Position'])
    maxim_h=int(stacked_h['Residues'].max())
    minim_h=int(stacked_h['Residues'].min())
    cmap_h=['#ff0000', '#ffd700', '#6897bb', '#ff9a00', '#03396c', '#00e5e5', '#20b2aa', '#b6b4b9', '#635d5a', '#ff99aa', '#c9caff', '#f0c39e','#c8e6c8' ]
    the_cmap_h=cmap_h[minim_h-1:maxim_h-1]
    stacked_all = df_all.stack().reset_index()
    stacked_all.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked_all=stacked_all.sort_values(by=['Position'])
    maxim_all=int(stacked_all['Residues'].max())
    minim_all=int(stacked_all['Residues'].min())
    cmap_all=['#ff0000', '#ffd700', '#6897bb', '#ff9a00', '#03396c', '#00e5e5', '#20b2aa', '#b6b4b9', '#635d5a', '#ff99aa', '#c9caff', '#f0c39e','#c8e6c8']
    the_cmap_all=cmap_h[minim_all-1:maxim_all-1]



    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(size))
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
    #plt.suptitle("Number of contacts - " + name, fontsize = 36) #title
    sns.heatmap(df_all, ax=ax1,xticklabels=True, yticklabels=True, linewidths=0.1, linecolor='black',annot=False, annot_kws={"size": 8}, cmap=(ListedColormap(the_cmap_all)), cbar=True)#.set_title("All contacts",y=1.1)
    ax1.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    ax1.set_yticklabels(ax1.get_yticklabels(),rotation=0)
    ax1.set_xticklabels(ax1.get_xticklabels(),rotation=90)
    ax1.set_aspect("equal")
    sns.heatmap(df_h, ax=ax2,xticklabels=True, yticklabels=True, linewidths=0.1, linecolor='black',annot=False, annot_kws={"size": 8}, cmap=(ListedColormap(the_cmap_h)), cbar=True)#.set_title("H bonded",y=1.1)
    ax2.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    ax2.set_yticklabels(ax2.get_yticklabels(),rotation=0)
    ax2.set_xticklabels(ax2.get_xticklabels(),rotation=90)
    ax2.set_aspect("equal")
    plt.show()


"""
groups= [(cm,"cm1"),(cm2,"cm2"),(cm3, "cm3"), (m1, "m1"),(m2,"m2")]
for chain in ["heavy"]:
    for cdr in [str(get_chain_letter(chain))+'1', str(get_chain_letter(chain))+'2',str(get_chain_letter(chain))+'3',"o"]:
        for (a,b) in groups:
            try:
                get_cdr_on_epitope(a, b, chain,"all", (15,15),5, cdr)
            except IndexError:
                continue
                
"""
def get_heatmap_ep_lesk(group, name, chain, type,  size):
    #contacts1 =json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_'+type+'.json').read())
    contacts_h_bonded=json.loads(open('equiv_contacts_dict_lsyozyme_pos_name_paratope_'+chain+'_h_bonded.json').read())
    contacts_all =json.loads(open('equiv_contacts_dict_lsyozyme_pos_name_paratope_'+chain+'_all.json').read())
    aa_type={"1": ["ASP", "GLU"],
             "2":["CYS","MET"],
             "3":["LYS","ARG"],
             "4":["SER","THR"],
             "5":["PHE","TYR"],
             "6":["ASN","GLN"],
             "7":["LEU","VA","ILE"],
             "8":["GLY"],
             "9":["ALA"],
             "10":["TRP"],
             "11":["HIS"],
             "12":["PRO"]}

    cont_dict_h_bonded={}
    for k, v in contacts_h_bonded.items():
        for pdb in group:
            if pdb in k:
                cont_dict_h_bonded[get_species(pdb)+"-"+pdb]={}
                for pos, names in v.items():
                    try:
                        if all(aa in ["GLY","ALA","SER","THR"]for aa in names) is True:

                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=1
                        elif all(aa in ["CYS", "MET","VAL","ILE","LEU","PRO","PHE","TYR","TRP"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=2
                        elif all(aa in ["ASN", "GLN","HIS"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=3

                        elif all(aa in ["ASP", "GLU"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=4

                        elif all(aa in ["LYS", "ARG"]for aa in names) is True:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=5

                        else:
                            cont_dict_h_bonded[get_species(pdb)+"-"+pdb][pos]=6


                    except KeyError:
                        continue
    cont_dict_all={}
    for k, v in contacts_all.items():
        for pdb in group:
            if pdb in k:
                cont_dict_all[get_species(pdb)+"-"+pdb]={}
                for pos, names in v.items():
                    try:
                        if all(aa in ["GLY","ALA","SER","THR"]for aa in names) is True:

                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=1
                        elif all(aa in ["CYS", "MET","VAL","ILE","LEU","PRO","PHE","TYR","TRP"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=2
                        elif all(aa in ["ASN", "GLN","HIS"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=3

                        elif all(aa in ["ASP", "GLU"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=4

                        elif all(aa in ["LYS", "ARG"]for aa in names) is True:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=5

                        else:
                            cont_dict_all[get_species(pdb)+"-"+pdb][pos]=6


                    except KeyError:
                        continue
    cont_h_bonded=cont_dict_h_bonded
    cont_all=cont_dict_all
    list_pos_all=[]
    for y, b in cont_dict_all.items():
        for pos2, v2 in b.items():
            list_pos_all.append(pos2)
    list_pos_h=[]
    for x, a in cont_dict_h_bonded.items():
        for pos1, v1 in a.items():
            list_pos_h.append(pos1)
    for db in group:
        for pos in list_pos_all:
            for posh in list_pos_h:
                if pos not in posh:
                    cont_h_bonded[get_species(pdb)+"-"+pdb][pos]=0
                if posh !=pos:
                    cont_all[get_species(pdb)+"-"+pdb][posh]=0
    df_all =pd.DataFrame.from_dict(cont_all,orient='index')
    #df.set_index('pdb',inplace=True)
    df_all=df_all[natsorted(df_all.columns)]
    df_all.groupby(by=species,axis=1)


    df_h =pd.DataFrame.from_dict(cont_h_bonded,orient='index')
    #df.set_index('pdb',inplace=True)
    df_h=df_h[natsorted(df_h.columns)]
    df_h.groupby(by=species,axis=1)
    stacked_h = df_h.stack().reset_index()
    stacked_h.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked_h=stacked_h.sort_values(by=['Position'])
    maxim_h=int(stacked_h['Residues'].max())
    minim_h=int(stacked_h['Residues'].min())
    cmap_h=['#ff0000', '#ffd700', '#6897bb', '#ff9a00', '#03396c', '#00e5e5', '#20b2aa', '#b6b4b9', '#635d5a', '#ff99aa', '#c9caff', '#f0c39e','#ffeffa']
    the_cmap_h=cmap_h[minim_h:maxim_h-1]
    stacked_all = df_all.stack().reset_index()
    stacked_all.rename(columns={'level_1': 'Position', 0: 'Residues'}, inplace=True)
    stacked_all=stacked_all.sort_values(by=['Position'])
    maxim_all=int(stacked_all['Residues'].max())
    minim_all=int(stacked_all['Residues'].min())
    cmap_all=['#ff0000', '#ffd700', '#6897bb', '#ff9a00', '#03396c', '#00e5e5', '#20b2aa', '#b6b4b9', '#635d5a', '#ff99aa', '#c9caff', '#f0c39e','#ffeffa']
    the_cmap_all=cmap_h[minim_all:maxim_all-1]



    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(size))
    fig.autofmt_xdate()
    plt.subplots_adjust(top=0.9) #adjust sublots labels to make room for title
    #plt.suptitle("Number of contacts - " + name, fontsize = 36) #title
    sns.heatmap(df_all, ax=ax1,xticklabels=True, yticklabels=True, linewidths=0.1, linecolor='black',annot=True, annot_kws={"size": 8}, cmap=(ListedColormap(the_cmap_all)), cbar=True).set_title("All contacts",y=1.1)
    ax1.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax1.set_aspect("equal")
    sns.heatmap(df_h, ax=ax2,xticklabels=True, yticklabels=True, linewidths=0.1, linecolor='black',annot=True, annot_kws={"size": 8}, cmap=(ListedColormap(the_cmap_h)), cbar=True).set_title("H bonded",y=1.1)
    ax2.tick_params(labelbottom='off',labeltop='on', labelsize=12)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax2.set_aspect("equal")
    plt.show()



groups= [(camel,"camel"),(mice, "mice"),(shark,"shark")]
for (a,b) in groups:
    #print(get_heatmap_cdr_bonding(a, b,"heavy","all" ,(10,10),5))
    print(len(a))




"""               
with open('aa_type_frequency_heavy.json', 'w') as ctr: #save the contact residues in a json file
    json.dump(get_heatmap_aa_type_cdr(unique, "unique_all_pdbs", "heavy", "all", (5,5), 5, "H3"), ctr)
"""

pos_epitope=json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_pos_all.json').read())
#print(pos_epitope['1c08']['102'])
dict_pos={}
def get_cdr_pattenrs():
    cds=get_heatmap_fw(pdb_codes1, "all","heavy","all" ,(10,10), 5)
    for k, v in cds.items():
        list_cdr=[]
        for pos, cd in v.items():
            if len(cd) >0:
                list_cdr.append(cd[0])
        dict_pos[k]=list_cdr
    df_all=pd.DataFrame.from_dict(dict_pos, orient="index")
    #resturn(df_all)
    df_all.to_excel("cdr_binding.xlsx")
    return(dict_pos)




with open('fw_binding.json', 'w') as ctr: #save the contact residues in a json file
    json.dump(get_cdr_pattenrs(), ctr)
