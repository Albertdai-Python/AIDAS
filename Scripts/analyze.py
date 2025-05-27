#IMPORT
import os
import csv

# PARAMETERS
home_dir = os.getcwd().split('Scripts')[0]
vina_dir = home_dir + '/Vina/vina'
data_dir = home_dir + '/Data/'
result_dir = home_dir + '/Result/'
input_dir = home_dir + '/Input/'
logs_dir = home_dir + '/Logs/'
smiles_dir = home_dir + '/SMILES/'
analysis_dir = home_dir + '/Analysis/'

# VARIABLES
ligands = []
data = []
iterations = 1 # HAS TO BE SAME AS MAIN.PY
protein = '1zum'
threshold = 0

# LIGANDS
os.chdir(input_dir)
try:
    with open('ligands.txt', 'r') as f:
        ligands = f.read().split('\n')
        for line in range(len(ligands)):
            zinc_id = ligands[line]
            if 'ZINC' not in zinc_id:
                while len(zinc_id) != 12:
                    zinc_id = '0' + zinc_id
                zinc_id = 'ZINC' + zinc_id
            ligands[line] = zinc_id
        f.close()
except Exception as e:
    e.add_note('ligands.txt does not exist.')
    raise

# DATA FORMATTING
for i in range(iterations * len(ligands)+1):
    data.append(['' for j in range(12)])
data[0] = ['PROTEIN', 'LIGAND', 'ITERATION', 'POS1', 'POS2', 'POS3', 'POS4', 'POS5', 'POS6', 'POS7', 'POS8', 'POS9']

# PRE-PARSING
data[1][0] = protein
for i in range(len(ligands)):
    # data[i*iterations+1][1] = ligands[i]
    # ^^single label^^
    for j in range(iterations):
        data[i*iterations+1+j][2] = str(j+1)
        data[i*iterations+1+j][1] = ligands[i]

# PARSING
os.chdir(logs_dir)
for i in range(len(ligands)):
    for j in range(iterations):
        filename = f'{protein}_{ligands[i]}_{j}.txt'
        try:
            with open(filename, 'r') as f:
                contents = f.readlines()
                if len(contents) >= 47:
                    for k in range(9):
                        data[i*iterations+1+j][3+k] = [x for x in contents[38+k].split(' ') if x][1]
                    f.close()
                else:
                    f.close()
        except Exception:
            pass

# FORMING TABLE
os.chdir(analysis_dir)
with open(f'{protein}.csv', 'w', newline = '') as f:
    writer = csv.writer(f)
    writer.writerow(data.pop(0))
    for i in data:
        if i[3] and float(i[3]) <= threshold:
            writer.writerow(i)
    f.close()
    
