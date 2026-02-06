"""
Alu_oligomer.py, zhipengluchina@gmail.com
python ~/Dropbox/_scripts/Alu_oligomer.py > Alu_list.urls

##### 1. search Dfam for Alu and got 43 familes, from 299nts to 320nts
store in file Alu_list.txt

##### 2. get all the urls
input1 = open("Alu_list.txt",'r').readlines()
Alu_list = [i.split()[0] for i in input1]
s1 = "https://www.dfam.org/api/families/"
s2 = "/sequence?format=fasta&download=true"
for i in Alu_list: print(s1+i+s2)
#python ~/Dropbox/_scripts/Alu_oligomer.py > Alu_list.urls #for this step. 

##### 3. download all the models into a file Alu_models.txt 
touch Alu_models.txt
wget --no-check-certificate -i Alu_list.urls -O - -o /dev/null >> Alu_models.txt

##### 4. analyze the models for oligomers
"""


import re
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

infile = "Alu_models.txt"; f = open("Alu_models.txt",'r')
models = {}; ID = ''
for line in f:
    if line[0] == '>': ID = line[1:-1]; models[ID] = ''
    else: models[ID] += line.strip()
dists = []
for k in models:
    #print(k,models[k])
    matches = re.finditer('AAAAA',models[k])
    l = [m.start() for m in matches]
    positions = [l[0]]+[l[i+1] for i in range(len(l)-1) if l[i+1]-l[i]!=5]
    dists.append([positions[i+1]-positions[i] for i in range(len(positions)-1)])
dists = [j for i in dists for j in i]
counted = Counter(dists)
countedall = [counted[i] if i in counted else 0 for i in range(170)]
plt.bar(list(range(170)),countedall); plt.savefig("a.pdf")
##### extract all A(n) motifs, taking the starting locations. 


