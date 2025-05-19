import os
import plotly.graph_objects as go
from sklearn.cluster import DBSCAN
import numpy as np
import plotly.express as px
import pandas as pd

home_dir = os.getcwd()
input_dir = home_dir + '/Input/'
smiles_dir = home_dir + '/SMILES'
result_dir = home_dir + '/Result/'
analysis_dir = home_dir + '/Analysis'

ligands = open(input_dir + 'ligands.txt', 'r').read().splitlines()
protein = '1zum'
iteration = 0
x = []
y = []
z = []

for ligand in ligands:
    with open(result_dir + f'{protein}_{ligand}_{iteration}.pdbqt', 'r') as f:
        text = f.read().split('MODEL')[1].split('\n')[1:]
        affinity = float([x for x in text[0].split(' ') if x][3])
        for line in text:
            if line.startswith('ATOM'):
                line = [x for x in line.split(' ') if x]
                x.append(float(line[5]))
                y.append(float(line[6]))
                z.append(float(line[7]))

enzyme_x = []
enzyme_y = []
enzyme_z = []

protein_text = open(result_dir + f'{protein}.pdbqt', 'r').read().splitlines()
for line in protein_text:
    if line.startswith('ATOM'):
        line = [x for x in line.split(' ') if x]
        enzyme_x.append(float(line[6]))
        enzyme_y.append(float(line[7]))
        enzyme_z.append(float(line[8]))

coords = np.column_stack((x, y, z))
clustering = DBSCAN(eps=1.0, min_samples=10).fit(coords)
labels = clustering.labels_
df = pd.DataFrame({
    'x': x,
    'y': y,
    'z': z,
    'cluster': labels.astype(str)  # convert to string for categorical coloring
})
df_filtered = df[(df['cluster'] != '-1') & (df['cluster'] != '3')]


fig = px.scatter_3d(df_filtered, x='x', y='y', z='z', color='cluster',
                    color_discrete_sequence=px.colors.qualitative.Set1)

fig.update_traces(marker=dict(size=2))
fig.add_trace(go.Scatter3d(
    x=enzyme_x,
    y=enzyme_y,
    z=enzyme_z,
    mode='markers',
    marker=dict(
        size=2,
        color='gold',
        opacity=1
    )
))


fig.show()

