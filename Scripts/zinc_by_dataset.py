#IMPORT
import os
import requests
from openbabel import openbabel
import ssl
import shutil

# INPUT
iterations = 1
protein = '1zum'

# UNVERIFIED CONTEXT
ssl._create_default_https_context = ssl._create_unverified_context

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

def change(string):
    while len(string) < 12:
        string = '0' + string
    string = 'ZINC' + string
    return string
# LIGAND EXTRACTION
os.chdir(input_dir)
with open('ligands.txt', 'r') as f:
    ligands = [change(x) for x in f.read().split('\n')]
    f.close()
os.chdir(data_dir)
library = os.listdir(data_dir)
conv = openbabel.OBConversion()
conv.SetInAndOutFormats('sdf', 'pdbqt')
count = 0
for i in ligands:
    
    if i+'.pdbqt' in library:
        count += 1
        print(f'-> [A] {i} ({count} / {len(ligands)})')
    else:
        try:
            response = requests.get(f'https://zinc.docking.org/substances/{i}.sdf')
            response.raise_for_status()
            mol = openbabel.OBMol()
            conv.ReadString(mol, response.text)
            builder = openbabel.OBBuilder()
            builder.Build(mol)
            mol.AddHydrogens()
            chargeModel = openbabel.OBChargeModel.FindType("gasteiger")
            chargeModel.ComputeCharges(mol)
            pdbqt_text = conv.WriteString(mol)
            with open(f'{i}.pdbqt', 'w') as f:
                f.write(pdbqt_text)
                f.close()
            count += 1
            print(f'-> [Y] {i} ({count} / {len(ligands)})')
            pass
        except Exception as e:
            ligands.remove(i)
            print(f'-> [N] {i}')
            pass
os.chdir(input_dir)
with open('ligands.txt', 'w') as f:
    f.write("\n".join(ligands))
    f.close()

# DOCKING
os.chdir(data_dir)
for l in ligands:
    for i in range(iterations):
        cmd = rf'{vina_dir} --receptor "#{protein}.pdbqt" --ligand {l}.pdbqt --config "#config.txt" --out "{protein}_{l}_{i}.pdbqt"'
        result = os.popen(cmd).read()
        with open(f'{protein}_{l}_{i}.txt', 'w') as f:
            f.write(result)
            f.close()
        with open(f'{protein}_{l}_{i}.txt', 'r') as f:
            content = f.readlines()
            if len(content) >= 39:
                affinity = [x for x in content[38].split(' ') if x][1]
            else:
                affinity = False
            f.close()
        print(f'-> [R] {protein} | {l} ({i+1} / {iterations}) | Affinity: {affinity}')
        if affinity:
            shutil.move(data_dir+f'{protein}_{l}_{i}.txt', logs_dir+f'{protein}_{l}_{i}.txt')
            shutil.move(data_dir+f'{protein}_{l}_{i}.pdbqt', result_dir+f'{protein}_{l}_{i}.pdbqt')
        else:
            os.remove(data_dir+f'{protein}_{l}_{i}.txt')
                
            
