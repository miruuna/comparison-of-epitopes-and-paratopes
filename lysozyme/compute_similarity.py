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
from scipy import stats


species = json.loads(open("species_lysozime.json").read())
pdb_codes= json.loads(open('pdb_codes_lysozyme.json').read())
the_contacts=json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_all_all.json').read())
non_bonded_dict= json.loads(open('equiv_contacts_dict_lsyozyme_all_non_h_bonded.json').read())
h_bonded_dict=json.loads(open('equiv_contacts_dict_lsyozyme_all_h_bonded.json').read())
unique=['1sq2','1t6v','2hfm','2i25','2i26','1a2y', '1bql', '1c08', '1dqj', '1dzb', '1ic4', '1ic5', '1ic7', '1jhl', '1kip', '1kiq', '1mlc', '1ndg', '1ndm', '1p2c', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '2iff', '4tsa', '4tsb', '4tsc', '4ttd', '1jto', '1op9', '1ri8', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '4i0c']
c_and_s=['1jto', '1jtp', '1jtt', '1mel', '1op9', '1ri8', '1rjc', '1sq2', '1t6v', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '2i25', '2i26', '3eba', '4i0c']
m1=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg']
m2=['1bql','1mlc', '2iff']
cm=['1jto', '1jtp', '1jtt','1xfp','1zmy','2i25','2i26']
cm2=['1zvy','1sq2','1t6v']
cma=['2i25','2i26']
cm3=['1ri8','1rjc','1zv5']
unique_minus_1=set(unique).difference(set(g1))
all_unique=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg','1bql','1mlc', '2iff','1jto', '1jtp', '1jtt','1xfp','1zmy','1zvy','1ri8','1rjc','1zv5','1sq2','1t6v','2i25','2i26']
mice=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg','1bql','1mlc', '2iff']
camel=['1jto', '1jtp', '1jtt','1xfp','1zmy','1zvy','1ri8','1rjc','1zv5']
shark=['1sq2','1t6v','2i25','2i26']


# get all pdbs from the 78 complexes that belong to one of the species
all_mice=[]
all_camel=[]
all_shark=[]
for pdb in pdb_codes:
    if "musculus" in species[pdb]:
        all_mice.append(pdb)
    if "drome" in species[pdb]:
        all_camel.append(pdb)
    if "cirrat" in species[pdb]:
        all_shark.append(pdb)


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
def get_jaccard_score(group,name, chain, type):

    contacts= json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_'+chain+'_'+type+'.json').read())
    non_bonded_list={}
    h_bonded_list={}
    for pdb in group:
        for k, v in contacts.items():
            if pdb in k:
                non_bonded_list[pdb]=[]
                for pos, res in v.items():
                    non_bonded_list[pdb].append(pos)
    d=non_bonded_list

    p={}
    key_list=[]
    for key, v in d.items():
        key_list.append(key)

    for k in key_list:
        p[k]=[]
        for i in key_list:
            if i!=k:
                while(len(d[k])!=len(d[i])):
                    if len(d[k])> len(d[i]):
                        d[i].append(0)
                    elif len(d[k])< len(d[i]):
                        d[k].append(0)
                p[k].append(1-distance.jaccard(set(d[k]),set(d[i])))
        dict_sum={}
        dict_avg={}
        for k, v in p.items():
            dict_sum[k]=0
            for e in v:
                dict_sum[k]+=e
            dict_avg[k]=(dict_sum[k]/len(p[k]))

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
    df_groups=pd.DataFrame({"Group 1"})
    df= pd.DataFrame({'Epitope_similarity':score_list,
                      'PDB':name_list
                      })

    df.index = np.arange(1, len(df)+1)
    sorted=df.sort_values('PDB')
    #g=sns.scatterplot(y="Jaccard_score", x="PDB_codes", data=sorted, hue='species')
    #g.set(xticks=[])
    #plt.xticks(rotation=90)
    #plt.show()
    #g.get_figure().savefig(name+"_jaccard_scatter_plot.png")
    #df.to_excel(name+"_"+chain+"_"+type+"_jaccard_scores.xlsx")
    return df
    #with open('mytable.txt', 'w') as tf:
    #   tf.write(df.to_latex())
def get_table():
    """
    Make table with the groups and their associated pdbs. Convert dataframe to excel """
    for x in [m1,m2,cm,cm2,cm3]:
        for y in [m1,m2,cm,cm2,cm3]:
            while len(x)<len(y):
                x.append(" ")

    d_t=pd.DataFrame({"Group 1a":m1,
                      "Group 1b":m2,
                      "Group 2a":cm,
                      "Group 2b":cm2,
                      "Group 2c":cm3})
    d_t.to_excel("new_grouping.xlsx")

def get_pie_chart(group, name):
    """
    Get pie chart of the distribution of species within a group
    """
    species_count={}
    for k, v in get_pdb_by_group(group).items():
        species_count[k]=len(v)
    df=pd.DataFrame.from_dict(species_count, orient="Index")
    #df.to_excel("species_distrib_"+name+".xlsx")
    df.index.rename('Species', inplace=True)
    df_melt=df.melt('Species', var_name='cols',  value_name='Frequency')
    fig1, ax1 = plt.subplots()
    ax1.pie('Species', labels='cols', autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.show()
    return df


def get_amino_acid(group,name):
    """
   Get the paratope amino acids that interact with the epitope.
   Returns dataframe containing amino acids and their occurence.
    """

    all_d= json.loads(open('equiv_contacts_dict_lsyozyme_pos_name_paratope_all_all.json').read())
    ep_all=[]
    aa_type_pos={}
    for k, v in all_d.items():
        if k in group:
            for pos, res in v.items():
                for r in res:
                    ep_all.append(r)
    for elem in ep_all:
        aa_type_pos[elem]=ep_all.count(elem)
    df=pd.DataFrame.from_dict(aa_type_pos, orient="index")
    df.index.rename('Amino Acid', inplace=True)
    df.rename(columns={ 0: name}, inplace=True)
    df.to_excel("aa_pref_paratope_"+name+".xlsx")
    return(ep_all)

def get_barplots_paratope_com():
    """
    Get the 3 barplots of paratope amino acid distribution for each of the three species.
    :return:
    """
    df_mice=pd.read_excel("aa_pref_paratope_mice.xlsx")
    df_mice.sort_values(by="Amino Acid",axis=0)
    df_mice['mouse'] = df_mice['mouse'].divide(19)
    df_camel=pd.read_excel("aa_pref_paratope_camel.xlsx")
    df_camel.sort_values(by="Amino Acid",axis=0)
    df_camel["camel"]=df_camel["camel"].divide(9)
    df_shark=pd.read_excel("aa_pref_paratope_shark.xlsx")
    df_shark.sort_values(by="Amino Acid",axis=0)
    df_shark["shark"]=df_shark["shark"].divide(4)
    df_merge_1=pd.merge(df_mice, df_camel, on="Amino Acid", how="right")
    df_merge_2=pd.merge(df_merge_1,df_shark, on="Amino Acid", how="right")
    df_melt = df_merge_2.melt('Amino Acid', var_name='species',  value_name='frequency')
    g=sns.barplot(x='Amino Acid', y="frequency", hue="species", data=df_melt)
    #fig, (ax1, ax2, ax3)=plt.subplots(1,3,figsize=(50,30))

    #sns.barplot(x="Amino Acid", y="mice", data=df_mice, ax=ax1)
    #sns.barplot(x="Amino Acid", y="camel", data=df_camel, ax=ax2)
    #sns.barplot(x="Amino Acid", y="shark", data=df_shark, ax=ax3)
    plt.show()



def get_stats(group, name):
    """
    Get different statistics computed with Jaccard scores.
    :return:
    """
    pairwise_sc=json.loads(open('pairwise_score'+name+'.txt').read())
    change_pairwise={}
    for k, v in pairwise_sc.items():
        change_pairwise[k]=v/1000000
    df_pairwise=pd.DataFrame.from_dict(change_pairwise,orient='index')
    df_pairwise.index.rename('PDB', inplace=True)
    df_pairwise.rename(columns={"PDB_codes": 'PDB', 0: 'Sequence similarity'}, inplace=True)
    df_jaccard_heavy=get_jaccard_score(group,name, "heavy", "all")

    df_jaccard_heavy.rename(columns={ "PDB_codes": 'PDB', "Jaccard_score": 'Heavy chain CR similarity'}, inplace=True)
    df_jaccard_light=get_jaccard_score(group,name, "light", "all")

    df_jaccard_light.rename(columns={"PDB_codes": 'PDB', "Jaccard_score": 'Light chain CR similarity'}, inplace=True)
    df_merge_2=pd.merge(df_jaccard_light, df_jaccard_heavy, on="PDB", how="left")
    df_merge_3=pd.merge(df_merge_2, df_pairwise, on="PDB", how="left")
    #df_merge_3 = df_merge_2.stack().reset_index()
    #f_melt = df_merge_2.melt('PDB', var_name='cols',  value_name='vals')
    #g = sns.catplot(x="PDB", y="vals", hue='cols', data=df_melt)
    #plt.show()
    df_merge_3.plot(x="PDB", y=["Light chain CR similarity", "Heavy chain CR similarity",'Sequence similarity'])
    plt.show()
    print(df_merge_2)

def get_paratope_sequence(group):
    """
    Get the paratope sequence obtained from the positions of the paratope residues that the epitope interacts with.
    """

    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    paratope=json.loads(open('equiv_paratope_all_all.json').read())
    dict_seq={}
    pos_paratope={}
    for k, v in paratope.items():
        if k in group:
            pos_list=[]
            for pos, a in v.items():
                pos_list.append(int(pos))

            pos_paratope[k]=sorted(pos_list) #find the antibody residue positions that the paratope is interacting with
    for k, pos in pos_paratope.items():
        aa_seq=[]
        for pos1 in pos:
            try:
                aa_seq.append(d[paratope[k][str(pos1)][0]]) #search the name of the antibody residue found at position pos1
            except IndexError:
                continue
        dict_seq[k]=aa_seq
    return dict_seq

def get_jaccard_paratopes(group):
    """
    Compute Jaccard similarity scores for the paratope sequences
    """

    key_list=[]
    d=get_paratope_sequence(group)
    for key, v in d.items():
        key_list.append(key)
    p={}
    for k in key_list:
        p[k]=[]
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
            dict_avg[k]=(dict_sum[k]/len(p[k]))

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
    df_groups=pd.DataFrame({"Group 1"})
    df= pd.DataFrame({'Paratope_similarity':score_list,
                      'PDB':name_list
                      })


    #df.sort_values("species", inplace=True)
    df.index = np.arange(1, len(df)+1)
    sorted=df.sort_values('PDB')
    #g=sns.scatterplot(y="Jaccard_score", x="PDB_codes", data=sorted, hue='species')
    #g.set(xticks=[])
    #plt.xticks(rotation=90)
    #plt.show()
    #g.get_figure().savefig(name+"_jaccard_scatter_plot.png")
    #df.to_excel(name+"_"+chain+"_"+type+"_jaccard_scores.xlsx")
    return df

def get_stats_by_group(group, name):
    pairwise_m1=get_jaccard_paratopes(m1)
    pairwise_m2=get_jaccard_paratopes(m2)
    pairwise_cm2=get_jaccard_paratopes(cm2)
    pairwise_cm1=get_jaccard_paratopes(cm)
    pairwise_cm3=get_jaccard_paratopes(cm3)
    seq_m1=pd.read_excel("m1_scores.xlsx")
    seq_m2=pd.read_excel("m2_scores.xlsx")
    seq_cm1=pd.read_excel("cm1_scores.xlsx")
    seq_cm2=pd.read_excel("cm2_scores.xlsx")
    seq_cm3=pd.read_excel("cm3_scores.xlsx")

    jac_m1=get_jaccard_score(m1, "m1","all","all")
    jac_m2=get_jaccard_score(m2, "m2","all","all")
    jac_cm1=get_jaccard_score(cm, "cm1","all","all")
    jac_cm2=get_jaccard_score(cm2, "cm2","all","all")
    jac_cm3=get_jaccard_score(cm3, "cm3","all","all")

    merge_cm11=pd.merge(seq_cm1, jac_cm1, on="PDB", how="right")
    merge_cm11['Epitope_similarity']=merge_cm11['Epitope_similarity'].multiply(400000000)
    merge_cm21=pd.merge(seq_cm2, jac_cm2, on="PDB", how="right")
    merge_cm21['Epitope_similarity']=merge_cm21['Epitope_similarity'].multiply(400000000)
    merge_cm31=pd.merge(seq_cm3, jac_cm3, on="PDB", how="right")
    merge_cm31['Epitope_similarity']=merge_cm31['Epitope_similarity'].multiply(400000000)
    merge_m11=pd.merge(seq_m1, jac_m1, on="PDB", how="right")
    merge_m11['Epitope_similarity']=merge_m11['Epitope_similarity'].multiply(400000000)
    merge_m21=pd.merge(seq_m2, jac_m2, on="PDB", how="right")
    merge_m21['Epitope_similarity']=merge_m21['Epitope_similarity'].multiply(40000000000)
    df_melt_cm11 = merge_cm11.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_cm21 = merge_cm21.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_cm31 = merge_cm31.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_m11 = merge_m11.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_m21 = merge_m21.melt('PDB', var_name='Legend',  value_name='scores')

    merge_cm1=pd.merge(pairwise_cm1, jac_cm1, on="PDB", how="right")
    merge_cm2=pd.merge(pairwise_cm2, jac_cm2, on="PDB", how="right")
    merge_cm3=pd.merge(pairwise_cm3, jac_cm3, on="PDB", how="right")
    merge_m1=pd.merge(pairwise_m1, jac_m1, on="PDB", how="right")
    merge_m2=pd.merge(pairwise_m2, jac_m2, on="PDB", how="right")
    df_melt_cm1 = merge_cm1.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_cm2 = merge_cm2.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_cm3 = merge_cm3.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_m1 = merge_m1.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_m2 = merge_m2.melt('PDB', var_name='Legend',  value_name='scores')

    fig,((ax11, ax12), (ax21, ax22), (ax31, ax32),(ax41, ax42), (ax51, ax52))=plt.subplots(5,2, figsize=(30,30))
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_cm1, ax=ax31, legend=None)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_cm2, ax=ax41, legend=None)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_cm3, ax=ax51, legend=None)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_m1, ax=ax11,legend=None)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_m2, ax=ax21,legend=None)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_cm11, ax=ax32, marker="D", legend=None)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_cm21, ax=ax42, marker="D",legend=None)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_cm31, ax=ax52, marker="D",legend=None)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_m11, ax=ax12, marker="D",legend=None)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_m21, ax=ax22, marker="D",legend=None)

    plt.show()
    a_cm1=stats.pearsonr(merge_cm1["Epitope_similarity"],merge_cm1["Paratope_similarity"] )
    a_cm2=stats.pearsonr(merge_cm2["Epitope_similarity"],merge_cm2["Paratope_similarity"] )
    a_cm3=stats.pearsonr(merge_cm3["Epitope_similarity"],merge_cm3["Paratope_similarity"] )
    a_m2=stats.pearsonr(merge_m2["Epitope_similarity"],merge_m2["Paratope_similarity"] )
    a_m1=stats.pearsonr(merge_m1["Epitope_similarity"],merge_m1["Paratope_similarity"] )
    a_cm11=stats.pearsonr(merge_cm11["Epitope_similarity"],merge_cm11["Ab_similarity"] )
    a_cm21=stats.pearsonr(merge_cm21["Epitope_similarity"],merge_cm21["Ab_similarity"] )
    a_cm31=stats.pearsonr(merge_cm31["Epitope_similarity"],merge_cm31["Ab_similarity"] )
    a_m21=stats.pearsonr(merge_m21["Epitope_similarity"],merge_m21["Ab_similarity"] )
    a_m11=stats.pearsonr(merge_m11["Epitope_similarity"],merge_m11["Ab_similarity"] )
    dictio={"Group 1a": a_m1,
            "Group 1b": a_m2,
            "Group 2a": a_cm1,
            "Group 2b": a_cm2,
            "Group 2c": a_cm3}
    df=pd.DataFrame.from_dict(data=dictio, orient="index")
    df.to_excel("stats.xlsx")

    dictio1={"Group 1a": a_m11,
             "Group 1b": a_m21,
             "Group 2a": a_cm11,
             "Group 2b": a_cm21,
             "Group 2c": a_cm31}
    df1=pd.DataFrame.from_dict(data=dictio1, orient="index")
    df1.to_excel("stats1.xlsx")

def get_jaccard_cdr(group, name):
    """
    Compute Jaccard similarity scores for CDR sequences.
    """
    key_list=[]
    d=json.loads(open(name+"_binding.json").read())
    for key, v in d.items():
        if key in group:
            key_list.append(key)
    p={}
    for k in key_list:
        p[k]=[]
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
            dict_avg[k]=(dict_sum[k]/len(p[k]))

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
    df_groups=pd.DataFrame({"Group 1"})
    df= pd.DataFrame({name+'_Preference_similarity':score_list,
                      'PDB':name_list
                      })


    #df.sort_values("species", inplace=True)
    df.index = np.arange(1, len(df)+1)
    sorted=df.sort_values('PDB')
    #g=sns.scatterplot(y="Jaccard_score", x="PDB_codes", data=sorted, hue='species')
    #g.set(xticks=[])
    #plt.xticks(rotation=90)
    #plt.show()
    #g.get_figure().savefig(name+"_jaccard_scatter_plot.png")
    #df.to_excel(name+"_"+chain+"_"+type+"_jaccard_scores.xlsx")
    return df

def get_stats_all(group):
    paratope=get_jaccard_paratopes(group)
    epitope=get_jaccard_score(group, "m1","all","all")
    cdr_binding=get_jaccard_cdr(group,"CDR")
    fw_binding=get_jaccard_cdr(group, "FR")
    merge1=pd.merge(epitope,paratope,  on='PDB', how="left")
    merge2=pd.merge(merge1, cdr_binding, on="PDB", how="left")
    merge3=pd.merge( cdr_binding,epitope, on="PDB", how="left")
    merge4=pd.merge( epitope,fw_binding, on="PDB", how="left")
    merge5=pd.merge( paratope,fw_binding, on="PDB", how="left")
    merge6=pd.merge( cdr_binding,fw_binding, on="PDB", how="left")
    df_melt_cm11 = merge3.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_cm12 = merge1.melt('PDB', var_name='Legend',  value_name='scores')
    df_melt_cm13=  merge4.melt('PDB', var_name='Legend',  value_name='scores')
    fig, (ax1, ax2, ax3)=plt.subplots(3,1)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_cm12, ax=ax1)
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_cm11, ax=ax2, palette="Set1")
    sns.lineplot(x='PDB', y='scores', hue='Legend', data=df_melt_cm13, ax=ax3, palette="Set2")

    plt.xticks(rotation=90)
    coeff=stats.pearsonr(merge1["Epitope_similarity"],merge1["Paratope_similarity"])
    coeff1=stats.pearsonr(merge2['CDR_Preference_similarity'], merge1["Epitope_similarity"])
    coeff2=stats.pearsonr(merge4['FR_Preference_similarity'], merge4["Epitope_similarity"])
    coeff3=stats.pearsonr(merge5['FR_Preference_similarity'], merge5["Paratope_similarity"])
    coeff34=stats.pearsonr(merge6['FR_Preference_similarity'], merge6["CDR_Preference_similarity"])
    plt.show()
    print("Ep-Para:",coeff, "Ep-CDR:", coeff1,"Ep-FR:", coeff2, "Para-FR:", coeff3, "CDR-FR:", coeff34)
def get_cdr_freq():
    """
    Get the contribution of each Heavy chain CDR to epitope binding in all of the 3 species
    """

    df_shark_all=pd.read_excel("cdr_freq_shark_heavy_all.xlsx")
    df_shark_all['All Contacts'] = df_shark_all['All Contacts'].divide(4)
    df_shark_all['H-Bonded'] = df_shark_all['H-Bonded'].divide(4)


    df_camel_all=pd.read_excel("cdr_freq_camel_heavy_all.xlsx")

    df_camel_all['All Contacts'] = df_camel_all['All Contacts'].divide(9)
    df_camel_all['H-Bonded'] = df_camel_all['H-Bonded'].divide(9)


    df_mice_all=pd.read_excel("cdr_freq_mice_heavy_all.xlsx")
    df_mice_all['H-Bonded'] = df_mice_all['H-Bonded'].divide(19)
    df_mice_all['All Contacts'] = df_mice_all['All Contacts'].divide(19)
    df_mice_all=df_mice_all[natsorted(df_mice_all.columns)]
    #df_mice_all['All Contacts'] = df_mice_all['All Contacts'].divide(19)

    df_melt_camel = df_camel_all.melt('CDR', var_name='cols',  value_name='Frequency')
    df_melt_shark = df_shark_all.melt('CDR', var_name='cols',  value_name='Frequency')
    df_melt_mice = df_mice_all.melt('CDR', var_name='cols',  value_name='Frequency')
    fig, (ax1,ax2,ax3)=plt.subplots(1,3,figsize=(30,10))
    sns.barplot(x="CDR", y="Frequency", hue='cols', data=df_melt_camel, ax=ax1).set_title("Camel")
    ax1.set_ylim([0,7])
    sns.barplot(x="CDR", y="Frequency", hue='cols', data=df_melt_shark, ax=ax2).set_title("Shark")
    sns.barplot(x="CDR", y="Frequency", hue='cols', data=df_melt_mice, ax=ax3).set_title("Mice")

    plt.show()

    print(df_camel_all)

def get_ep_positions_groups(group,name):
    """
    Get the (equivalent) positions for epitope residues in each group
    """
    all_d= json.loads(open('equiv_contacts_dict_lsyozyme_pos_paratope_all_h_bonded.json').read())
    ep_all=[]
    aa_type_pos={}
    for k, v in all_d.items():
        if k in group:
            for pos, res in v.items():
                ep_all.append(int(pos))
    #return(list(set(ep_all)))
    return(ep_all)

