# Dependencies

import distance
import matplotlib.pyplot as plt
import json
import numpy as np
import matplotlib.pyplot as plt


species = json.loads(open("species_lysozime.json").read())
d=json.loads(open('equiv_contacts_list_all_contacts.json').read())


p={}

key_list=[]
for key, v in d.items():
    key_list.append(key)

for k in key_list:
    cos1=0
    for i in key_list:
        if i!=k:
            while(len(d[k])!=len(d[i])):
                if len(d[k])> len(d[i]):
                    d[i].append(0)
                elif len(d[k])< len(d[i]):
                    d[k].append(0)
            cos1+=(1-distance.jaccard(set(d[k]),set(d[i]))) #Jaccard scores = 1-Jaccard distance
    p[k]=cos1/(len(key_list)-1)

p_sorted={k: v for k, v in sorted(p.items(), key=lambda x: x[1])}
name_list=[]
score_list=[]

for k, v in p_sorted.items():
    name_list.append(k)
    score_list.append(v)
   
data=(score_list, name_list)
species_d=species_list
colors = ("red", "green", "blue")


# Create plot
fig = plt.figure()
plt.plot(score_list,name_list, 'o', color="black")
plt.yticks([])

plt.title('Distribution of Jaccard scores')
plt.legend(loc=2)
plt.savefig("scatter.png")
plt.show()
