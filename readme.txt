# FairOmics Complex Networks Project

## Overview
This repository contains scripts and folder structure for analyzing microbial co-occurrence networks from DSMZ data.
Here it ise the link to the folder containing the raw DSMZ data file, BC vs CC .html scatter plot file and 3D spring layout plot of the graph:
https://liveunibo-my.sharepoint.com/:f:/g/personal/simone_baldisserri_studio_unibo_it/IgCzf8URtT6dRZh4DUajXIaWAfIGSKNI-I92HnCa2j0TT-E?e=nz072D 

---

## Repository Structure

complex_fairomics/
|── scripts/
│   ├── dsmz_processing.py       # Preprocessing: parsing DSMZ data and building bacteria-habitat matrix
│   ├── dsmz_matrix.py           # Builds sparse matrices of bacterial similarity from DSMZ data
│   ├── graph.py                 # Constructs weighted similarity graph, computes degree, betweenness, clustering 
|   |                                 coefficient, sparsifies graph, performs 3D layout, Louvain community detection
│   ├── community_analysis.py    # Maps nodes to taxonomic ranks, summarizes communities, exports CSV tables with
|   |                                 node info and community rank distributions
│   └── graphics.py              # Generates plots: 3D community scatter plots, degree histograms, BC vs CC scatterplots,
|                                     taxonomic rank distributions, stacked bar plots
|
├── project_data/
│   ├── initial_data/            # Here must be the DSMZ raw data file, downloadable from the link in the TOP but 
|   |                               contains a describing image of raw data structure
│   └── staging_data/            # Preprocessed and intermediate files
├── results/                     # Generated figures

---

## Requirements

- Python >= 3.12.3 
- Packages (highly suggested to be installed via `pip install -r requirements.txt`, see below):

numpy
scipy
pandas
matplotlib
plotly
networkx
python-louvain
tqdm
ete3
collection
re
json
os
colorsys

**Recommended:** Create a virtual environment before running the scripts:

python -m venv env
source env/bin/activate 
pip install --upgrade pip
pip install -r requirements.txt

---

## Usage
0. **Fundamentals**
Download `DSMZ_Habitat.txt` from the link on TOP and put it in `project_data/initial_data`

1. **Preprocessing and matrix construction**

python scripts/dsmz_matrix.py

- Reads `DSMZ_Habitat.txt` from `project_data/initial_data`
- Runs `dsmz_processing.py`, creating Bacteria dict and filtering taxids and habitats
- Builds sparse bacteria-habitat matrix
- Computes pairwise dissimilarity and similarity between taxa and generates and saves plots:
  - Distribution of pairwise dissimilarities (`distribuzione_dissimilarity_mat22_no_fungi.png`)
  - Distribution of pairwise similarities (`distribuzione_similarity_mat22_no_fungi.png`)
  - Similarity heatmap (`similarity_22_heatmap_no_fungi.html`)
- Saves COO sparse similarity matrix (`similarity_coo_mat22_no_fungi.npz`)

2. **Graph construction and analysis**

python scripts/graph.py

- Loads COO sparse similarity matrix
- Builds weighted NetworkX graph
- Performs and saves intermediate results in `project_data/staging_data`:
  - Degree computation (initial `degrees_i.csv` and after sparsification `degrees_ii.csv`)
  - Betweenness centrality (approximation) and Clustering coefficient (`ncbi_bc_cc(2000-70%)_no_fungi.txt`)
  - 3D spring layout and Louvain community detection (`graph_layout_data_no_fungi.json`)

3. **Community taxonomic analysis**

python scripts/community_analysis.py

- Loads graph and community detection information
- Maps NCBI taxids to taxonomic ranks and scientific names
- Summarizes communities and saves in `project_data/staging_data/`:
  - Full rank distributions per community (`community_rank_distribution.csv`)
  - Node info with community membership (`ncbi_node_info_with_communities.csv`)
  - Summaries of top 4 ranks per community (`summaries_of_communities.csv`)

4. **Visualization**

python scripts/graphics.py

- Generates plots:
  - 3D community scatter plots
  - Degree distributions (histograms)
  - Betweenness vs. Clustering coefficient .html scatterplots
  - Global and per-community taxonomic rank distributions
  - Stacked bar plot and boxplot of rank composition and dispersion in the entire graph
- Saves figures in `results/`

---

## Output

- `project_data/staging_data/`
  - `labels_no_fungi.json` → filtered NCBI and habitat labels
  - `similarity_coo_mat22_no_fungi.npz` → COO sparse similarity matrix
  - `graph_layout_data_no_fungi.json` → node positions (3D) and community assignments
  - `ncbi_bc_cc(2000-70%)_no_fungi.txt` → betweenness centrality and clustering coefficient per node
  - Degree CSV files: `degrees_i.csv`, `degrees_ii.csv`
  - Community rank summaries and full distributions:
    - `summaries_of_communities.csv`
    - `ncbi_node_info_with_communities.csv`
    - `community_rank_distribution.csv`

- `results/` → all plots

---

## Notes

- The pipeline assumes a dense **sparse graph**; This very large datasets (>6000 nodes) requires at least 12GB of RAM.  
- `dsmz_processing.py` is **only imported** and not meant to be executed standalone (even if it works by itself).
- All scripts are written in Python 3.12.3 , and the virtual environment should be recreated using `requirements.txt`.  

