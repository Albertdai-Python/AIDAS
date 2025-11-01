# AIDAS
## ALDH2 In silico Docking & Agonist-Screening 
### Source Code for High School Research "Hangover Cure – Searching for ALDH2 E487K Agonists" in Collaboration with China Medical University
**Lead Programmer:** Yu-Cheng Dai\
**Other Research Team Members:** Hong-Wei Chen, Yu-Lun Chen\
**Mentors:** Prof. Shou-Lun Lee, Prof. Yi-Chuan Li
 
![alt text](https://github.com/Albertdai-Python/AIDAS/blob/main/Media/AIDAS_mac_icon_2.png?raw=true)
---
### Project Introduction & Motivation
ALDH2 is an enzyme responsible for human alcohol metabolism, and the E487K East-Asian mutant produces an inactive version, causing harm to the human body. Through molecular docking simulations with [Autodock Vina](https://github.com/ccsb-scripps/AutoDock-Vina), we can acquire the binding position and affinity of a ligand and a given protein (ALDH2), thereby determining its likelihood of affecting enzyme activity. However, there are two major flaws to this approach: 
- First, given the large number of chemical and biological compounds, we would have to hand-pick compounds. 
- Second, it is difficult to determine whether a binding ligand acts as an inhibitor or an agonist to the enzyme solely via molecular docking. 

To combat these issues, our project considers the molecular structure of existing enzyme agonists as a template, performing mass screening on compounds from the [ZINC database](https://zinc.docking.org/) using molecular similarity and ranking them from highest to lowest similarity for docking priority.

For compounds with a molecular similarity (to an existing agonist) above a given threshold, our scripts perform automated molecular docking using an executable version of [Autodock Vina](https://github.com/ccsb-scripps/AutoDock-Vina).

Our scripts then compare enzyme residues near existing agonists to the target ligand, calculating a **Pair Score** to allow easy *in-silico* evaluation of whether a ligand is a potential agonist.

Our pipeline provides automation for both molecular similarity comparison and molecular docking, taking protein code, existing agonists, and other Vina-specified parameters as input, while outputting a CSV file for final evaluation.

---
### Achieved Results
Our project has helped us perform over 560 million similarities comparisons and 9100 molecular docking instances. 

Ultimately, 7 compounds were identified as potential ALDH2 agonists.

---
### Future Directions
- Test script viability with other enzymes
- Design full-fledged Pair Score algorithm for individual residue weights & biases
- Incorporate self-designed neural network with molecular similarity for pre-docking screening
- Integrate ADMET (absorption, distribution, metabolism, excretion, toxicity) prediction after Pair Score calculation
- Build GUI for easier usage


---
### Python Module Prerequisites
- matplotlib
- rdkit
- sqlite3
- plotly
- sklearn
- numpy
- pandas
- openbabel

---
### Detailed Functions of Each Script
- `main.py`
  - ***Docking Preparation and Similarity Calculation***
  - Downloads *SMILES* and *ZINC_ID* in tranches from `/Input/ZINC-downloader-3D-smi.txt`
  - Writes *ZINC_ID* and *SMILES* into `/SMILES/zinc_compounds.db` with ***sqlite3***
  - Converts compounds from *SMILES* 
  - Uses ***MACCSKeys*** and ***Tanimoto Coefficient*** to compute similarity between compounds and **Alda-1**, a known ALDH2 agonist
  - Outputs `/SMILES/Similarity Ranking.csv`
---
- `zinc_by_dataset.py`
  - ***Main docking Script***
  - Grabs compounds listed in `/Input/ligands.txt` and check if they are originally in `Data/`
  - Grab the compounds from `https://zinc.docking.org/substances/{compound}.sdf`
  - Adds *Hydrogens* and *Gasteiger Charge* to compound
  - Rewrites `/Input/ligands.txt` with the available compounds and excludes failed ones
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
- `analyze.py`
  - ***Compound Filtering via Affinity***
  - Accepts a threshold value as input, which can filter out low affinity compounds
  - Reads `/Input/ligands.txt` to get list of compounds
  - Writes the protein, ligand, iteration, affinities of the nine positions into `/Analysis/{protein}({threshold}).csv`
---
- `distance_search.py`
  - ***Calculates Minimum Distances of Compounds to Key Residues; Implements Pair Score to Rank Agonist Possibility***
  - Defines blacklist as residues that the compound should not interact with, likely because it is the substrate or coenzyme binding site
  - Defines whitelist as residues that the compound should interact with
  - Reads `/Result/{protein}.pdbqt` for blacklist and whitelist residues, saving their coordinates as groups for each residue
  - Reads `/Result/{protein}_{ligand}_{iteration}.pdbqt` for the ***FIRST POSE*** of docking results, saving their coordinates in a group
  - Implements minimum distance algorithm for each ligand / residue pair
    - Takes one atom each from ligand / residue
    - Calculate the ***Euclidean Distance*** between the all atom combinations
    - Outputs minimum distance
  - Implements Pair Score algorithm for each ligand
    - Pair Score = numbers of whitelist residues that are within 4.0Å of the compound
    - If any blacklist residues are within 4.0Å of the compound, automatically set Pair Score to 0
    - Writes protein, ligand, iteration, affinity, and distances to whitelist and blacklist residues into `/Analysis/Residue Distance.csv`
---
- `pose_distribution.py`
  - ***3D SCATTERPLOT for Docking Position Clustering***
  - Reads `/Input/ligands.txt` for ligand names
  - Iterates through each ligand: `/Result/{protein}_{ligand}_{iteration}.pdbqt` and read the atom positions of the first position
  - Uses DBSCAN with epsilon value of `1.0` for position clustering
  - Filters out `cluster -1`, which is noise and `cluster 3`, which consists of atoms less than one compound
  - Reads `/Result/{protein}.pdbqt` for atom positions of the protein
  - Adds the protein to the plot with different color
