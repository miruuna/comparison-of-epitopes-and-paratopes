from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.pairwise2 import format_alignment
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt
import seaborn as sns
from bs4 import BeautifulSoup
import requests
import json
from difflib import SequenceMatcher
import pandas as pd
s1=['1a2y', '1fdl', '1g7h', '1g7i', '1g7j', '1g7l', '1g7m', '1kir', '1vfb']
s2=['1c08', '1j1o', '1j1p', '1j1x', '2dqi', '2dqj', '2znw', '3a67', '3a6b', '3a6c', '3d9a', '3hfm']
#s3=['1dqj', '1nby', '1nbz','1ua6', '1uac']
s4=['2eiz', '2eks', '2yss']
s5=['1jto', '1jtp', '1jtt']
g1a=['1ri8', '1bql', '1ndg', '1ndm', '1p2c', '1zv5', '1mlc', '2iff']
g1b=[ '1a2y', '1kiq', '1kip']
g1c=['1jto', '1xfp']
g2=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg']
g1_mice= g2
g2_mice= g1a
g3_cs=['1jto','1jtp', '1jtt','1rjc','1xfp','1zmy','1t6v','1sq2','2i25','2i26']

pdb_codes1= json.loads(open('pdb_codes_lysozyme.json').read())
cdrs= json.loads(open('CDRS_update3.json').read())
cdrs_chain=json.loads(open('CDRs_by_chain.json').read())
jaccard_score=json.loads(open("jaccard_scores.json").read())
species=json.loads(open('species_lysozime.json').read())
pairwise_sc=json.loads(open('pairwise_score.txt').read())
def get_pairwise_score():
    pairwise_dict={}
    for k, v in pairwise_sc.items():
        pairwise_dict[k]=v/71000000
    return pairwise_dict
def flatten_column(df, column_name):
    repeat_lens = [len(item) if item is not np.nan else 1 for item in df[column_name]]
    df_columns = list(df.columns)
    df_columns.remove(column_name)
    expanded_df = pd.DataFrame(np.repeat(df.drop(column_name, axis=1).values, repeat_lens, axis=0), columns=df_columns)
    flat_column_values = np.hstack(df[column_name].values)
    expanded_df[column_name] = flat_column_values
    expanded_df[column_name].replace('nan', np.nan, inplace=True)
    return expanded_df

def get_list(pdb):
    page_link = "http://dunbrack2.fccc.edu/PyIgClassify/Result/ChainCdrClusterInfo.aspx?pdb="+pdb.upper()
    page_response = requests.get(page_link, timeout=5)
    soup = BeautifulSoup(page_response.text, "html.parser")
    td_list=[]
    for tr_tag in soup.find_all('tr'):
        tr_list=[]
        for td in tr_tag.find_all('td'):
            tr_list.append(td.text)
        td_list.append(tr_list)
    return td_list

def get_ind_cdrs(cdr,pdb):
    td_list= get_list(pdb)
    l1_list=[]
    for l1 in td_list:
        if cdr in l1:
            l1_list.append(l1[-2])
    return l1_list[-1]

def get_dict_cdrs(pdb):
    d={}
    for chain in ['L1', 'L2', 'L3', 'H1', 'H2','H3']:
        try:
            d[chain]= get_ind_cdrs(chain, pdb)
        except IndexError:
            continue
    return d
def get_pdb_dict():
    d1=cdrs
    d={}
    for pdb in ['1rjc', '1sq2', '1t6v', '2hfm', '2i25', '2i26']:
        try:
            d[pdb]=get_dict_cdrs(pdb)
        except KeyError:
            continue
    d1.update(d)
    return d1

def get_cdr_by_chain():
    d1={}
    for k, v in cdrs.items():
        d1[k]={'light': '',
               'heavy': ''}
        for c,seq in v.items():
            if "L" in c:
                d1[k]['light']+=seq
            elif "H" in c:
                d1[k]['heavy']+=seq
    return d1

def seq_matcher(chain):
    key_list=[]
    cdr_similarity={}
    score_sum=0
    for k, v in cdrs_chain.items():
        key_list.append(k)
    for k in key_list:
        cdr_similarity[k]=[]
        cos1=0
        for i in key_list:
            if(i!=k):
                s=SequenceMatcher(None, cdrs_chain[k][chain], cdrs_chain[i][chain])
                score=s.ratio()
                score_sum+=score
                cdr_similarity[k].append(score)
    sc_sum={}
    len_sc={}
    for k, v in cdr_similarity.items():
        len_sc[k]=len(v)
        sc_sum[k]=0
        for e in v:
            sc_sum[k]+=e

    avg_sc={}
    for k, v in sc_sum.items():
        avg_sc[k]=v/len(sc_sum)

    return avg_sc

def seq_matcher_individual(chain):
    key_list=[]
    cdr_similarity={}
    score_sum=0
    for k, v in cdrs.items():
        key_list.append(k)
    for k in key_list:
        cdr_similarity[k]=[]
        cos1=0
        for i in key_list:
            if(i!=k):
                s=SequenceMatcher(None, cdrs[k][chain], cdrs[i][chain])
                score=s.ratio()
                score_sum+=score
                cdr_similarity[k].append(score)
    sc_sum={}
    len_sc={}
    for k, v in cdr_similarity.items():
        len_sc[k]=len(v)
        sc_sum[k]=0
        for e in v:
            sc_sum[k]+=e

    avg_sc={}
    for k, v in sc_sum.items():
        avg_sc[k]=v/len(sc_sum)

    return avg_sc
def get_unique_cdrs(chain):


    flipped = {}

    for key, value in seq_matcher_individual(chain).items():
        if value not in flipped:
            flipped[value] = [key]
        else:
            flipped[value].append(key)

    df =pd.DataFrame.from_dict(flipped,orient='index')
    unique_cdrs=[]
    for k, v in flipped.items():
        unique_cdrs.append(v[0])

    df2 =pd.DataFrame.from_dict(flipped,orient='index')
    with open(chain+'_similarity_para.txt', 'w') as tf:
        tf.write(df2.to_latex())
    return unique_cdrs

def get_fasta(group):
    filenames = []
    name=str(group)
    with open(name+'for_MSA.fasta', 'w') as outfile:
        for pdb in group:
            outfile.write(">"+pdb+"\n"+get_cdr_by_chain()[pdb]['heavy']+"\n")

#for group in [s1, s2,s3,s4,s5]:
    #get_fasta(group)

def compare_score():
    pdb=['1sq2','1t6v','1ua6','2dqh','3a6c','2i25','2dqf','2dqe','1j1x','2i26','1zvy','2dqd','3a67','3a6b','2znw','2hfm','1c08','2dqj','2dqi','1rjc','1uac','2dqg','1j1p','1ri8','1jtt','1op9','3eba','1xfp','1zmy','1jto','1jtp','4i0c','1dzb','1fdl','1kir','1a2y','1vfb','1kiq','1kip','1g7h','1g7i','1g7m','1g7j','1g7l','1bql','2iff','1yqv','1mlc','1p2c','3d9a','3hfm','1ndg','1ndm','1xgr','1dqj','1nby','1nbz','1xgp','1xgt','1xgu','1xgq','2yss','2eiz','2eks','1jhl','1zv5','1ic4','2dqc','1zvh','1ic7','1ic5','1j1o']
    d_sc=['0.03739','0.03739','0.03969','0.01384','0.00252','0.01465','0.01465','0.02407','0.00548','0.0044','0.0044','0.0278','0.0077','0.00184','0.07452','0.06543','0.02281','0.00903','0.0001','0.03168','0.03168','0.01028','0.00053','0.05061','0.05061','0.02086','0.02047','0.07435','0.08037','0.00461','0.00063','0.13456','0.27224','0.03685','0.01443','0.00996','0.00507','0.00022','0.00225','0.00224','0.00224','0.00224','0.00193','0.00856','0','0','0.00458','0.03092','0.03069','0.00227','0.00246','0.01913','0.00785','0.00258','0','0','0','0.00118','0.00118','0.00118','0.00114','0.00303','0.00595','0.01042','0.01042','0.02224','0.01097','0.00057','0.03111','0.03111','0.01148','0.0002']
    dist_sc=pd.DataFrame({'PDB':pdb ,
                          'dist_sc': d_sc})
    d_sc=pd.read_excel("cladogram.xlsx")
    pairwise_sc=pd.DataFrame.from_dict(get_pairwise_score(),orient='index')
    pairwise_sc.index.rename('PDB', inplace=True)
    pairwise_sc.rename(columns={ 0: 'pairwse_score'}, inplace=True)
    jaccard_h=pd.read_excel("h_bonded_jaccard_scores.xlsx")
    jaccard_non_h=pd.read_excel("non_h_bonded_jaccard_scores.xlsx")
    heavy_sc=pd.DataFrame.from_dict(seq_matcher('heavy'),orient='index')
    heavy_sc.index.rename('PDB', inplace=True)
    heavy_sc.rename(columns={ 0: 'heavy_score'}, inplace=True)
    species_df=pd.DataFrame.from_dict(species,orient='index')
    species_df.index.rename('PDB', inplace=True)
    species_df.rename(columns={ 0: 'species'}, inplace=True)
    light_sc=pd.DataFrame.from_dict(seq_matcher('light'),orient='index')
    light_sc.index.rename('PDB', inplace=True)
    light_sc.rename(columns={ 0: 'light_score'}, inplace=True)
    jaccard_sc=pd.DataFrame.from_dict(jaccard_score,orient='index')
    jaccard_sc.index.rename('PDB', inplace=True)
    jaccard_sc.rename(columns={ 0: 'jaccard_score'}, inplace=True)
    df_merge_seq_sim=pd.merge(d_sc, pairwise_sc, on="PDB")
    df_merge_col = pd.merge( light_sc, heavy_sc, on='PDB')
    df_merge_2=pd.merge(df_merge_col, pairwise_sc, on="PDB")
    df_merge_1=pd.merge(jaccard_sc, pairwise_sc, on="PDB")
    df_merge_3=pd.merge(df_merge_col, df_merge_seq_sim, on="PDB")
    df_merge_jaccard1=pd.merge(jaccard_h, jaccard_non_h, on="PDB")
    df_merge_jaccard2=pd.merge(df_merge_jaccard1, jaccard_sc, on="PDB")
    df_merge_jaccard_all=pd.merge(df_merge_jaccard1, df_merge_col, on="PDB")
    #df_merge_4=pd.merge(df_merge_3, species_df, on="PDB")
    df_melt = df_merge_jaccard_all.melt('PDB', var_name='cols',  value_name='vals')
    g = sns.catplot(x="PDB", y="vals", hue='cols', data=df_melt)
    #g=sns.lineplot(x='PDB', y=['ligh_score', 'heavy_score'], data=df_merge_2)

    #g.set(xticks=[])
    plt.xticks(rotation=90)
    plt.show()
    g.get_figure().savefig("lineplot_jacc_H_V.png")
    return jaccard_non_h

def flip_dict(dictionary):
    flipped = {}

    for key, value in dictionary.items():
        if value not in flipped:
            flipped[value] = [key]
        else:
            flipped[value].append(key)
    return flipped
def get_species_groups():
    d_gr={"I": g1_mice,
                  "II": g2_mice,
                  "III": g3_cs}
    d_gr_species={}
    for k, v in d_gr.items():
        d_gr_species[k]={}
        for pdb in v:
            d_gr_species[k][pdb]=species[pdb]
    flipped_dict={}
    for grp, v in d_gr_species.items():
        flipped_dict[grp]=flip_dict(v)
    group_composition={}
    for grp, v in flipped_dict.items():
        group_composition[grp]={}
        for s, list_pdb in v.items():
            group_composition[grp][s]=len(list_pdb)


    df =pd.DataFrame.from_dict(d_gr,orient='index').T
    with open('group_table.txt', 'w') as tf:
        tf.write(df.to_latex())


    return
print(compare_score())





#with open('CDRs_by_chain.json', 'w') as ctr: #save the contact residues in a json file
#   json.dump(get_cdr_by_chain(), ctr)

