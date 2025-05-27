import os
import csv
import math
# Calculate distances to key residues and implementing a scoring function

home_dir = os.getcwd().split('Scripts')[0]
vina_dir = home_dir + '/Vina/vina'
data_dir = home_dir + '/Data/'
result_dir = home_dir + '/Result/'
input_dir = home_dir + '/Input/'
logs_dir = home_dir + '/Logs/'
smiles_dir = home_dir + '/SMILES/'
analysis_dir = home_dir + '/Analysis/'

def distance(mol1, mol2):
    return math.sqrt((mol1[0]-mol2[0])**2 + (mol1[1]-mol2[1])**2 + (mol1[2]-mol2[2])**2)

class Protein:
    def __init__(self, file, blacklist, whitelist):
        self.file = file
        self.blacklist = {}
        self.whitelist = {}
        for residue in blacklist:
            self.blacklist[residue] = []
        text = open(result_dir+self.file, 'r').read().split('\n')
        for line in text:
            if line.startswith('ATOM'):
                line = [x for x in line.split(' ') if x]
                residue = int(line[5])
                if residue in self.blacklist:
                    self.blacklist[residue].append((float(line[6]), float(line[7]), float(line[8])))
        self.blacklist = dict(sorted(self.blacklist.items()))
        for residue in whitelist:
            self.whitelist[residue] = []
        text = open(result_dir+self.file, 'r').read().split('\n')
        for line in text:
            if line.startswith('ATOM'):
                line = [x for x in line.split(' ') if x]
                residue = int(line[5])
                if residue in self.whitelist:
                    self.whitelist[residue].append((float(line[6]), float(line[7]), float(line[8])))
        self.whitelist = dict(sorted(self.whitelist.items()))

class Ligand:
    def __init__(self, ligand, iteration, protein = '1zum'):
        self.file = protein + '_' + ligand + '_' + str(iteration) + '.pdbqt'
        self.affinity = 0

        self.atoms = []
        text = open(result_dir+self.file, 'r').read().split('MODEL')[1].split('\n')[1:]
        self.affinity = float([x for x in text[0].split(' ') if x][3])
        for line in text:
            if line.startswith('ATOM'):
                line = [x for x in line.split(' ') if x]
                self.atoms.append((float(line[5]), float(line[6]), float(line[7])))
        self.matrix = [protein, ligand, iteration+1, self.affinity, '']
    def score(self, protein):
        self.pair_score = 0
        for white in protein.whitelist:
            min_distance = 0
            distance_list = []
            for atom1 in protein.whitelist[white]:
                for atom2 in self.atoms:
                    distance_list.append(distance(atom1, atom2))
            min_distance = min(distance_list)
            if min_distance <= 4.0:
                self.pair_score += 1
            self.matrix.append(min_distance)
        self.matrix.append('')
        for black in protein.blacklist:
            min_distance = 0
            distance_list = []
            for atom1 in protein.blacklist[black]:
                for atom2 in self.atoms:
                    distance_list.append(distance(atom1, atom2))
            min_distance = min(distance_list)
            if min_distance <= 4.0:
                self.pair_score = 0
            self.matrix.append(min_distance)
        self.matrix.insert(4, self.pair_score)
        if self.pair_score == 0:
            return False
        elif self.pair_score == len(protein.whitelist):
            return self.matrix
        else:
            return False


blacklist = [268, 302]
whitelist = [120, 124, 170, 173, 292, 296, 457, 458, 459]
protein = Protein('1zum.pdbqt', blacklist, whitelist)
os.chdir(input_dir)
ligands = open('ligands.txt', 'r').read().split('\n')
os.chdir(analysis_dir)
writer = csv.writer(open('Residue Distance.csv', 'w'))
first_row = ['PROTEIN', 'LIGAND', 'ITERATION', 'AFFINITY', f'COMPATABILITY SCORE (MAX={len(whitelist)})', 'WHITELIST']
first_row.extend(whitelist)
first_row.append('BLACKLIST')
first_row.extend(blacklist)
writer.writerow(first_row)
for ligand in ligands:
    result = Ligand(ligand, 0, protein='1zum').score(protein)
    if result:
        writer.writerow(result)
