import os
import csv
import math
import matplotlib.pyplot as plt

home_dir = os.getcwd()
input_dir = home_dir + '/Input/'
smiles_dir = home_dir + '/SMILES'
result_dir = home_dir + '/Result/'
analysis_dir = home_dir + '/Analysis'

os.chdir(analysis_dir)
reader = csv.reader(open('Residue Distance (Omitted Zeros).csv', 'r'))
sim_dict = {}
for row in reader:
    if row != ['PROTEIN','LIGAND','ITERATION','AFFINITY','COMPATABILITY SCORE (MAX=9)','WHITELIST','120','124','170','173','292','296','457','458','459','BLACKLIST','268','302']:
        if sim_dict.get(float(row[4])) != None:
            sim_dict[float(row[4])] += 1
        else:
            sim_dict[float(row[4])] = 1
sorted_dict = dict(sorted(sim_dict.items()))
x = list(sorted_dict.keys())
y = list(sorted_dict.values())
'''
x.insert(0, 0)
y.insert(0, 8740)
'''
m = max(y)
i = y.index(max(y))
print(list(zip(x,y)))
plt.plot(x, y, marker='o', linestyle='-')
plt.title('Compatability Distribution')
plt.xlabel('Pair Score')
plt.ylabel('Number of Compounds')
plt.show()