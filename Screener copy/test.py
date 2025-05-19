import os
import numpy as np
import plotly.graph_objects as go

# Directories (adjust as needed)
home_dir = os.getcwd()
input_dir = os.path.join(home_dir, 'Input')
result_dir = os.path.join(home_dir, 'Result')

# Parameters
protein = '1zum'
iteration = 0

# Read ligand list
with open(os.path.join(input_dir, 'ligands.txt'), 'r') as f:
    ligands = f.read().splitlines()

# Parse ligand docking results: extract centroid coordinates and affinities
coordinates = []
affinities = []

for ligand in ligands:
    pdbqt_path = os.path.join(result_dir, f"{protein}_{ligand}_{iteration}.pdbqt")
    with open(pdbqt_path, 'r') as f:
        content = f.read()
    # Split on MODEL to get first model
    model_text = content.split('MODEL')[1].split('\n')[1:]
    # Extract affinity from first line
    affinity = float([x for x in model_text[0].split(' ') if x][3])
    x_sum = y_sum = z_sum = atom_count = 0
    for line in model_text:
        if line.startswith('ATOM'):
            parts = [p for p in line.split(' ') if p]
            x_sum += float(parts[5])
            y_sum += float(parts[6])
            z_sum += float(parts[7])
            coordinates.append((float(parts[5]), float(parts[6]), float(parts[7])))
            affinities.append(affinity)
            atom_count += 1
    # Compute centroid
    x_cent = x_sum / atom_count
    y_cent = y_sum / atom_count
    z_cent = z_sum / atom_count

coordinates = np.array(coordinates)
affinities = np.array(affinities)

# Normalize affinities so that better (more negative) scores have higher weights
weights = affinities.max() - affinities
weights /= weights.max()

# Define voxel grid bounds (adjust to cover enzyme or region of interest)
xmin, xmax = 15, 85
ymin, ymax = 35, 125
zmin, zmax = 10, 80

# Number of voxels along each axis
Nx, Ny, Nz = 50, 50, 50

# Initialize empty voxel grid
voxel_grid = np.zeros((Nx, Ny, Nz))

# Map each ligand centroid to voxel index and accumulate weighted counts
for (x, y, z), w in zip(coordinates, weights):
    ix = int((x - xmin) / (xmax - xmin) * Nx)
    iy = int((y - ymin) / (ymax - ymin) * Ny)
    iz = int((z - zmin) / (zmax - zmin) * Nz)
    # Boundary check
    if 0 <= ix < Nx and 0 <= iy < Ny and 0 <= iz < Nz:
        voxel_grid[ix, iy, iz] += w

# Create voxel center coordinates
vx = (xmax - xmin) / Nx
vy = (ymax - ymin) / Ny
vz = (zmax - zmin) / Nz

x_centers = np.linspace(xmin + vx/2, xmax - vx/2, Nx)
y_centers = np.linspace(ymin + vy/2, ymax - vy/2, Ny)
z_centers = np.linspace(zmin + vz/2, zmax - vz/2, Nz)

X, Y, Z = np.meshgrid(x_centers, y_centers, z_centers, indexing='ij')

# Plot voxel grid as volume with Plotly
fig = go.Figure(data=go.Volume(
    x=X.ravel(),
    y=Y.ravel(),
    z=Z.ravel(),
    value=voxel_grid.ravel(),
    isomin=voxel_grid.min(),
    isomax=voxel_grid.max(),
    opacity=0.1,          # Adjust transparency
    surface_count=20,     # Number of isosurfaces
    colorscale='Reds',    # Color scale
    colorbar=dict(title='Affinity-weighted count')
))

fig.update_layout(scene=dict(
    xaxis_title='X (Å)',
    yaxis_title='Y (Å)',
    zaxis_title='Z (Å)',
    aspectmode='data'
))

fig.show()