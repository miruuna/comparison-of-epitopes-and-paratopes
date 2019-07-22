from sklearn.metrics.pairwise import cosine_similarity
import json

print(cosine_similarity([[1, 0, -1]], [[-1,-1, 0]]))
d=json.loads(open('equiv_contacts_list_all_contacts.json').read())

a={'a': [2 , 3, 4],
   'b': [3,4,5],
    'c':[4,1,3],
   'd':[5,1,2]}

p={}
count=0
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
            cos1+=cosine_similarity([d[k]],[d[i]])
    p[k]=cos1/(len(key_list)-1)

for k in key_list:
    print(k, p[k])
