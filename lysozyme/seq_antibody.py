import subprocess
import json
import re
import os
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
from Bio import pairwise2
import json
import pandas as pd
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np
from Bio.pairwise2 import format_alignment
m1=['2dqd', '1xgr', '1xgt', '1xgq', '1dqj', '2dqc', '1c08', '1xgu', '1ic5', '1ic4', '1ic7', '2dqh', '1ua6', '1xgp', '2dqe', '2dqg']
m2=['1bql','1mlc', '2iff']
cm=['1jto', '1jtp', '1jtt','1xfp','1zmy','2i25','2i26']
cm2=['1zvy','1sq2','1t6v']
cm3=['1ri8','1rjc','1zv5']
the_ag=["4i0c"]

annotation = json.loads(open('data_lysozyme.json').read())
pdb_codes= json.loads(open('pdb_codes_lysozyme.json').read())
unique=['1sq2','1t6v','2hfm','2i25','2i26','1a2y', '1bql', '1c08', '1dqj', '1dzb', '1ic4', '1ic5', '1ic7', '1jhl', '1kip', '1kiq', '1mlc', '1ndg', '1ndm', '1p2c', '1ua6', '1xgp', '1xgq', '1xgr', '1xgt', '1xgu', '2dqc', '2dqd', '2dqe', '2dqf', '2dqg', '2dqh', '2eiz', '2iff', '4tsa', '4tsb', '4tsc', '4ttd', '1jto', '1op9', '1ri8', '1xfp', '1zmy', '1zv5', '1zvh', '1zvy', '4i0c']

def get_annotation_string(pdb):
    """
    :param pdb:
    :return: the antigen chains as string separated by commas as required by the pdb_tofasta command
    """
    antigen_list=''
    if len(annotation[pdb]['antigen'])>2:
        antigen_list=annotation[pdb]['antigen'][0]#+','+annotation[pdb]['antibody'][1]
    else:
        for antigen in annotation[pdb]['antigen']:
            antigen_list=antigen_list+','+antigen
    return antigen_list[1:] #omit the , at the beginning of string

def get_fasta_files(group, name):
    for pdb in group:
        cmd= "pdb_fetch "+pdb+" | pdb_tofasta >> ./fasta1/"+pdb+".fasta"
        subprocess.run(cmd, shell=True)

def process_fasta():
    cmd="cat *.fasta > merged.fasta"
    subprocess.run(cmd, shell=True)




def get_sequences_from_pdb():
    pdb_seq={}
    elem_list =[]
    with open(os.getcwd()+"\\antibody_lysozyme_for_MSA2.txt", 'r') as f:
        contents = f.read()
    contents2 = contents.rstrip("\n")
    string1=str(contents2).replace('>PDB|', ":")
    string12=string1.replace('X','')
    string2=string12.replace('|', '>')
    string3=string2.replace('\n', '')
    p=re.compile(r'(\w+):')
    for item in p.findall(string3):
        elem_list.append(item)
    for i in list(range(0, len(pdb_codes))):
        pdb_seq[pdb_codes[i]]=elem_list[i]
    return pdb_seq

def edit_fasta_for_msa():
    with open(os.getcwd()+"\\fasta_antibody\\merged.fasta", 'r') as f:
        lines_duplicate = []
        lines=[]
        counter=0
        for line in f:
            if line.startswith('>'):
                counter+=1
                lines_duplicate.append(line[:4]+"|"+pdb_codes[counter-1]+line[4:])
            else:
                lines_duplicate.append(line)
        return lines_duplicate

def write_new_file():
    with open('updated_antibody_lysozyme_for_MSA.txt', 'w') as f2:
        for line in edit_fasta_for_msa():
            f2.write(line)


groupa_b=['1xgt', '1xgq']
def get_alignment(group, name):
    matrix = matlist.blosum62
    p={}
    a={}
    b={}
    for k in group:
        b[k]=[]
        p[k]=0
        for i in group:
            if i!=k:
                for al1,al2, score, begin, end in pairwise2.align.globalms(get_sequences_from_pdb()[k], get_sequences_from_pdb()[i], 1, -5, -10, -5):
                    p[k]+=score
        a[k]=abs(p[k])/len(group)

    df=pd.DataFrame.from_dict(data=a, orient="index")
    df.to_excel(name+"_scores.xlsx")
    return a

def get_alignment_matrix(group, name):
    matrix = matlist.blosum62
    p={}
    a={}
    b={}
    for k in group:
        b[k]=[]
        p[k]={}
        for i in group:
            if i!=k:
                p[k][i]=0
                for al1,al2, score, begin, end in pairwise2.align.globalms(get_sequences_from_pdb()[k], get_sequences_from_pdb()[i], 1, -5, -1, -1):
                    p[k][i]+=score
            else:
                p[k][i]=0
        a[k]=abs(p[k])/(len(group)-1)


    df=pd.DataFrame.from_dict(data=p, orient="index")
    df.to_excel(name+"_scores.xlsx")
    cmap=sns.palplot(sns.light_palette("purple"))
    fig, ax = plt.subplots(figsize=(30,30))
    sns.heatmap(df,cmap=cmap,ax=ax)
    plt.show()
    return a

#for (a,b) in [(m1,"m1"),(m2,"m2"),(cm,"cm1"),(cm2,"cm2"),(cm3,"cm3")]:
   # print(get_alignment_matrix(a,b))
print(len(unique))
