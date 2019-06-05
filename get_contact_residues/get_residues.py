
import os
import csv

rootdir = os.getcwd()
def get_subdir():
    subdir_list =[]
    for dirs, subdirs, files in os.walk(rootdir+'\\out_\\'): ##find subdirectories and store them as strings
        for subdir in subdirs:
            a_dir = str((os.path.join(subdir)))
            subdir_list.append(a_dir)
        return(subdir_list)

def store_data(): ##creates a nested dictionary containing the contact residues
    file_dict ={}
    for subdir in get_subdir():
        for root, dir, files in os.walk(rootdir+'\\out_\\'+subdir+"\\"):
            file_dict[subdir]={}
            for file in files:
                if file.endswith('.csv'):
                    chains_dict={}
                    chains_dict[file]=[]
                    with open(rootdir+'\\out_\\'+subdir+"\\"+file) as csvfile:
                        reader = csv.reader(csvfile, delimiter=',')
                        for row in reader:
                            chains_dict[file].append(row)

                        if file == 'molecule_1.txt.csv':
                            chains_dict['antibody'] = chains_dict.pop('molecule_1.txt.csv')
                        else:
                            chains_dict['antigen'] = chains_dict.pop('molecule_2.txt.csv')
                    file_dict[subdir].update(chains_dict)


    return(file_dict)



dict = store_data()

def get_chain():
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



#print(store_data())
print(get_chain())
