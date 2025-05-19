import os
import csv
import math
import matplotlib.pyplot as plt

home_dir = os.getcwd()
input_dir = home_dir + '/Input/'
smiles_dir = home_dir + '/SMILES'
result_dir = home_dir + '/Result/'
analysis_dir = home_dir + '/Analysis'

os.chdir(smiles_dir)
reader = csv.reader(open('Similarity Ranking.csv', 'r'))
sim_dict = {}
for row in reader:
    if row != ['ZINC_ID', 'SIMILARITY']:
        if sim_dict.get(float(row[1])) != None:
            sim_dict[float(row[1])] += 1
        else:
            sim_dict[float(row[1])] = 1
sorted_dict = dict(sorted(sim_dict.items()))
x = list(sorted_dict.keys())
y = list(sorted_dict.values())
m = max(y)
i = y.index(max(y))
print(x[i], m)
plt.scatter(x, y, marker='.')
plt.title('Similarity Distribution')
plt.xlabel('Similarity')
plt.ylabel('Number of Compounds')
plt.show()