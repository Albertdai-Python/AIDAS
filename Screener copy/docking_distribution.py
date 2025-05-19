import os
import csv
import math
import matplotlib.pyplot as plt
import statistics

home_dir = os.getcwd()
input_dir = home_dir + '/Input/'
smiles_dir = home_dir + '/SMILES'
result_dir = home_dir + '/Result/'
analysis_dir = home_dir + '/Analysis'

os.chdir(analysis_dir)
reader = csv.reader(open('1zum.csv', 'r'))
sim_dict = {}
for row in reader:
    if row != ["PROTEIN","LIGAND","ITERATION","POS1","POS2","POS3","POS4","POS5","POS6","POS7","POS8","POS9"]:
        if sim_dict.get(round(float(row[3]), 1)) != None:
            sim_dict[round(float(row[3]), 1)] += 1
        else:
            sim_dict[round(float(row[3]), 1)] = 1
sorted_dict = dict(sorted(sim_dict.items()))
x = list(sorted_dict.keys())
y = list(sorted_dict.values())
mu = sum(y)/len(y)
sigma = statistics.stdev(y)
m = max(y)
i = y.index(max(y))
print(x[i], m)
print(mu, sigma)
plt.plot(x, y, marker='o', linestyle='-')
plt.title('Affinity Distribution')
plt.xlabel('Affinity')
plt.ylabel('Number of Compounds')
plt.show()