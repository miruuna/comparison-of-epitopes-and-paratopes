import matplotlib.pyplot as plt
import os
import csv
import pandas as pd
import seaborn as sns
import plotly.plotly as py
import plotly.graph_objs as go
import numpy as np
import skimage;
from skimage.color import rgb2gray

rootdir = os.getcwd()
def get_subdir():
    """
    :return: array of subdir names, representing the pdb_codes
    """
    subdir_list =[]
    for dirs, subdirs, files in os.walk(rootdir+'\\out_HA1\\'): ##find subdirectories and store them as strings
        for subdir in subdirs:
            a_dir = str((os.path.join(subdir)))
            subdir_list.append(a_dir)
        return(subdir_list)

def store_data(): ##creates a nested dictionary containing the contact residues
    """
    Loops through every directory in the out_/ folder and for each txt file it extracts 4the coordinates of each position
    :return: a dictionary of dictionary of arrays
    """
    file_dict ={}
    for subdir in get_subdir():
        for root, dir, files in os.walk(rootdir+'\\out_HA1\\'+subdir+"\\"):
            file_dict[subdir]={}
            for file in files:
                if file.endswith('.csv'):
                    chains_dict={}
                    chains_dict[file]=[]
                    with open(rootdir+'\\out_HA1\\'+subdir+"\\"+file) as csvfile:
                        reader = csv.reader(csvfile, delimiter=',')
                        for row in reader:
                            chains_dict[file].append(row) ##rename dict key molecule_1 and 2 to antibody and antigen

                        if file == 'molecule_1.txt.csv':
                            chains_dict['antibody'] = chains_dict.pop('molecule_1.txt.csv')
                        else:
                            chains_dict['antigen'] = chains_dict.pop('molecule_2.txt.csv')
                    file_dict[subdir].update(chains_dict)


    return(file_dict)



dict = store_data()

def get_chain():
    """
    Function to extract just the information about light and heavy chains
    :return: dictionary with keys: heavy and light chains
    """
    dict = store_data()
    dict_chains={}

    light={}
    heavy ={}
    for subdir in get_subdir():
        dict_chains[subdir]={}
        for i in range(0, len(dict[subdir]['antibody'])-1):
            #return dict[subdir]['antibody'][i].__getitem__(2)
            if dict[subdir]['antibody'][i][0]=='A':
                position= dict[subdir]['antibody'][i].__getitem__(1)
                contact = dict[subdir]['antibody'][i].__getitem__(2)
                light[position]=contact

            elif dict[subdir]['antibody'][i][0]=='B':
                position= dict[subdir]['antibody'][i].__getitem__(1)
                contact = dict[subdir]['antibody'][i].__getitem__(2)
                heavy[position]=contact
        dict_chains[subdir]['light']=[light]
        dict_chains[subdir]['heavy']=[heavy]
        #light[i].update(contact)
    return dict_chains

def get_epitope():
    """
    Extracts just the coordinates of the epitopes for all the pdb codes
    :return: a dictionary of array
    """
    dict = store_data()
    dict_epitopes={}
    for subdir in get_subdir():
        dict_epitopes[subdir]=[]
        for i in range(0, len(dict[subdir]['antigen'])-1):
            contact = dict[subdir]['antigen'][i].__getitem__(2)
            position= dict[subdir]['antigen'][i].__getitem__(1)
            if contact =='C':
                dict_epitopes[subdir].append(position)
    return dict_epitopes

def get_epitope_list():
    """
    Transforms into an 2D array the dict of arrays containing coordinates of contact residues in the epitopes
    :return: 2D Array
    """
    dict = store_data()
    list_epitopes=[]
    for subdir in get_subdir():
        list_1=[]
        for i in range(0, len(dict[subdir]['antigen'])-1):
            contact = dict[subdir]['antigen'][i].__getitem__(2)
            position= dict[subdir]['antigen'][i].__getitem__(1)
            if contact =='C':
                list_1.append(position)
        list_epitopes.append(list_1)
    return list_epitopes




#print(store_data())
#print(get_chain())
#print(get_epitope_list())

#d = get_epitope_list()


#s = pd.DataFrame(d, index=get_subdir())
#s.fillna(value=pd.np.nan, inplace=True)
#ax = sns.heatmap(s)
#plt.show()
#print(s)

pdb_codes = get_subdir()
a_list =''
for code in pdb_codes:
    a_list=a_list+','+code
print(a_list)