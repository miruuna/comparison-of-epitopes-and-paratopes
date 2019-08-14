import subprocess
import json
import re
import os
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

annotation = json.loads(open('data_lysozyme_update.json').read())
pdb_codes= json.loads(open('pdb_codes_lysozyme.json').read())
s1=['1a2y', '1fdl', '1g7h', '1g7i', '1g7j', '1g7l', '1g7m', '1kir', '1vfb']
s2=['1c08', '1j1o', '1j1p', '1j1x', '2dqi', '2dqj', '2znw', '3a67', '3a6b', '3a6c', '3d9a', '3hfm']
s3=['1dqj', '1nby', '1nbz','1ua6', '1uac']
s4=['2eiz', '2eks', '2yss']
s5=['1jto', '1jtp', '1jtt']
sim_gr_0157=['1c08', '1j1o', '1j1p', '1j1x', '2dqi', '2dqj', '2znw', '3a67', '3a6b', '3a6c', '3d9a', '3hfm']
pdb_codes_update=['4tsa', '4tsb', '4tsc','4ttd','4pgj']
def get_annotation_string(pdb):
    """
    :param pdb:
    :return: the antigen chains as string separated by commas as required by the pdb_tofasta command
    """
    antigen_list=''
    if len(annotation[pdb]['antigen'])>2:
        antigen_list=annotation[pdb]['antigen'][0]+','+annotation[pdb]['antigen'][1]
    else:
        for antigen in annotation[pdb]['antigen']:
            antigen_list=antigen_list+','+antigen
    return antigen_list[1:] #omit the , at the beginning of string

def get_fasta_files():
    for pdb in pdb_codes_update:
        cmd= "pdb_fetch "+pdb+" | pdb_selchain -"+get_annotation_string(pdb)+" | pdb_tofasta >> ./fasta_lysozyme/"+pdb+".fasta"
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


def get_alignment():
    matrix = matlist.blosum62
    p={}
    for k in pdb_codes:
        p[k]=0
        for i in pdb_codes:
            if i!=k:
                for al1,al2, score, begin, end in pairwise2.align.globalms(get_sequences_from_pdb()[k], get_sequences_from_pdb()[i], 2, -1, -.5, -.1):
                    p[k]+=score
    return p
with open('pairwise_score.txt', 'w') as ctr: #save the contact residues in a json file
   json.dump(get_alignment(), ctr)