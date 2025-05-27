import requests
import os
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.DataStructs import TanimotoSimilarity
import csv
import sqlite3

home_dir = os.getcwd().split('Scripts')[0]
vina_dir = home_dir + '/Vina/vina'
data_dir = home_dir + '/Data/'
result_dir = home_dir + '/Result/'
input_dir = home_dir + '/Input/'
logs_dir = home_dir + '/Logs/'
smiles_dir = home_dir + '/SMILES/'
analysis_dir = home_dir + '/Analysis/'

alda_smiles = 'C1OC2=C(O1)C=C(C=C2)CNC(=O)C3=C(C=CC=C3Cl)Cl'
alda = Chem.MolFromSmiles(alda_smiles)
alda_fp = MACCSkeys.GenMACCSKeys(alda)

os.chdir(smiles_dir)
db_connection = sqlite3.connect('zinc_compounds.db')
db_cursor = db_connection.cursor()
write_cursor = db_connection.cursor()
db_cursor.execute('''
CREATE TABLE IF NOT EXISTS compounds (
    zinc_id TEXT PRIMARY KEY,
    smiles TEXT NOT NULL
)
''')
db_connection.commit()
db_cursor.execute('''
CREATE TABLE IF NOT EXISTS similarities (
    zinc_id TEXT PRIMARY KEY,
    similarity REAL NOT NULL
)
''')
db_connection.commit()



#'''
os.chdir(input_dir)
url_list = open('ZINC-downloader-3D-smi.txt', 'r').read().split('\n')

os.chdir(smiles_dir)
for i in range(len(url_list)):
    url = url_list[i]
    print(f'Retrieving {url} ({i+1}/{len(url_list)})')
    response = requests.get(url)
    try:
        response.raise_for_status()
        whole_text = response.text.split('\n')[1:][:-1]
        response = False
        print(f'Downloading {len(whole_text)} SMILES structures')
        compound_list = []
        for line in whole_text:
            smiles = line.split(' ')[0]
            zinc_id = line.split(' ')[1]
            # Turn Zinc ID into number format
            if 'ZINC' in zinc_id:
                zinc_id = zinc_id[4:]
                while zinc_id[0] == '0':
                    zinc_id = zinc_id[1:]
            compound_list.append((zinc_id, smiles))
            smiles = False
            zinc_id = False
        db_cursor.executemany(
            'INSERT OR IGNORE INTO compounds (zinc_id, smiles) VALUES (?, ?)',
            compound_list
        )
    except requests.exceptions.HTTPError:
        print(f'Download failed for {url}')
db_connection.commit()
#'''

db_cursor.execute("SELECT COUNT(*) FROM compounds")
num_rows = db_cursor.fetchone()[0]
db_cursor.execute("SELECT zinc_id, smiles FROM compounds")
pointer = 1
for row in db_cursor:
    compound = Chem.MolFromSmiles(row[1])
    if compound is not None:
        fingerprint = MACCSkeys.GenMACCSKeys(compound)
        tanimoto_val = TanimotoSimilarity(alda_fp, fingerprint)
        write_cursor.execute(
            'INSERT OR IGNORE INTO similarities (zinc_id, similarity) VALUES (?, ?)',
            (row[0], tanimoto_val)
        )
    if pointer % 10000 == 0:
        print(f'Successfully downloaded {pointer}/{num_rows} SMILES structures')
    pointer += 1
print(f'Finished downloading {num_rows} SMILES structures')
db_connection.commit()

writer = csv.writer(open('Similarity Ranking.csv', 'w'))
writer.writerow(["ZINC_ID", "SIMILARITY"])
db_cursor.execute("SELECT zinc_id, similarity FROM similarities ORDER BY similarity DESC")
for row in db_cursor:
    writer.writerow(row)