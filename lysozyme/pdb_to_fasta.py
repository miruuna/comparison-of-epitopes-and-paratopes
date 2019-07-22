import subprocess
import json
import re
import os



annotation = json.loads(open('data_gp41.json').read())
pdb_codes= json.loads(open('pdb_codes_gp41.json').read())


def get_annotation_string(pdb):
    """
    :param pdb:
    :return: the antigen chains as string separated by commas as required by the pdb_tofasta command
    """
    antigen_list=''
    for antigen in annotation[pdb]['antigen']:
        antigen_list=antigen_list+','+antigen
    return antigen_list[1:] #omit the , at the beginning of string

def get_fasta_files():
    for pdb in pdb_codes:
        cmd= "pdb_fetch "+pdb+" | pdb_selchain -"+annotation[pdb]['antigen'][0]+" | pdb_tofasta >> ./fasta_gp41/"+pdb+".fasta"
        subprocess.run(cmd, shell=True)

def process_fasta():
    cmd="cat *.fasta > merged.fasta"
    subprocess.run(cmd, shell=True)



def get_sequences_from_pdb():
    pdb_seq={}
    elem_list =[]
    with open(os.getcwd()+"\\fasta_gp41\\merged.fasta", 'r') as f:
        contents = f.read()
    contents2 = contents.rstrip("\n")
    string1=str(contents2).replace('>PDB|', ":")
    string12=string1.replace('X','')
    string2=string12.replace('|', '>')
    string3=string2.replace('\n', '')
    p=re.compile(r'(\w+):')
    for item in p.findall(string3[1:]):
        elem_list.append(item)
    for i in range(0, len(pdb_codes)-1):
        pdb_seq[pdb_codes[i]]=elem_list[i]
    return pdb_seq

def edit_fasta_for_msa():
    with open(os.getcwd()+"\\fasta_lyso\\merged.fasta", 'r') as f:
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
    with open('lysozyme_for_MSA.txt', 'w') as f2:
        for line in edit_fasta_for_msa():
            f2.write(line)

write_new_file()
