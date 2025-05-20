# CLEAR-FACE
## Computational Ligand Exploration for ALDH2 Restoration – Finding Activators &amp; Catalytic Enhancers
---
### What does each script do?
- `main.py`
  - Downloads *SMILES* and *ZINC_ID* in tranches from `/Input/ZINC-downloader-3D-smi.txt`
  - Writes *ZINC_ID* and *SMILES* into `/SMILES/zinc_compounds.db` with ***sqlite3***
  - Converts compounds from *SMILES* 
  - Uses ***MACCSKeys*** and ***Tanimoto Coefficient*** to compute similarity between compounds and **Alda-1**, a known ALDH2 agonist
  - Outputs `/SMILES/Similarity Ranking.csv`
---
- `zinc_by_dataset.py`
  - ***This is the script for docking!***
  - Grabs compounds listed in `/Input/Zinc_dataset.txt` and check if they are originally in `Data/`
  - Grab the compounds from `https://zinc.docking.org/substances/{compound}.sdf`
  - Adds *Hydrogens* and *Gasteiger Charge* to compound
  - Rewrites `/Input/Zinc_dataset.txt` with the available compounds and excludes failed ones
  - Docking
      - Inputs
        - Vina directory is at `/Vina/vina`, which is an executable for Vina 1.2.0
        - `/Data/#{protein}.pdbqt` as macromolecule
        - `/Data/{ligand}.pdbqt` as ligand
        - `/Data/#config.txt` as config file, which contains docking center, grid box size, energy range, and exhaustiveness
      - Outputs
        - `/Logs/{protein}_{ligand}_{iteration}.txt` is the log file, which contains the affinities for nine poses
        - `/Result/{protein}_{ligand}_{iteration}.pdbqt` is the docking position file, which contains the coordinates of the nine poses
---
- `analyze.py
