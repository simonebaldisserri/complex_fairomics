from collections import defaultdict
from collections import Counter
import re
import scipy.sparse as sp
import json
import numpy as np
import matplotlib.pyplot as plt

data_files_path = "/home/selvaggio/Desktop/complex networks/Project/"

#Putting all lines (splitted) in array (of arrays)
with open(data_files_path + "DSMZ_Habitat.txt", 'r') as text:
    info = [line.strip().split("\t") for line in text]

bacteria = {}

#Filling bacteria dict as: Bacterium[BD_CODE] = [[NCBI/PATH/TO/BACTERIUM], [LAST/FIVE/HABITATS/OF/PATHS]]
for i in range(len(info)):
    bd_code = info[i][8]
    ncbi_path = [t for t in re.split(r'/ncbi:|/', info[i][3])[1:] if t and not t.startswith("bd:")]
    habitat_paths = info[i][7].split(',')
    habitat_path = set()
    for hp in habitat_paths:
        habitat_path.update([h for h in re.split(r'/OBT:|/', hp)[-5:] if h])
    if bd_code in bacteria:                                                                         #if bacterium already exists it adds new habitats
        bacteria[bd_code][1].update(habitat_path)
    else:
        bacteria[bd_code] = [ncbi_path, habitat_path]
for bd_code in bacteria:
    bacteria[bd_code][1] = list(bacteria[bd_code][1])

#Habitats filter >3
habitat_counter = Counter()
for _, habitats in bacteria.values():
    habitat_counter.update(habitats)

excluded_habitats = {"000001"}
min_bacteria_per_habitat = 3
filtered_habitats = {habitat for habitat, count in habitat_counter.items() if count >= min_bacteria_per_habitat and habitat not in excluded_habitats}

#Taxonomy ID filter >3
ncbi_counter = Counter()
for ncbi, _ in bacteria.values():
    ncbi_counter.update(ncbi)

excluded_ncbi = {"1"}
min_bacteria_per_ncbi = 3
filtered_ncbi = {ncbi for ncbi, count in ncbi_counter.items() if count >= min_bacteria_per_ncbi and ncbi not in excluded_ncbi}

#Creating matrix
matrix = defaultdict(lambda: defaultdict(set))          #this structure automatically creates the dict whenever we call a key that doesn't exist
row_labels = list(filtered_ncbi)
column_labels = list(filtered_habitats)
r_label_to_index = {label: idx for idx, label in enumerate(row_labels)}
c_label_to_index = {label: idx for idx, label in enumerate(column_labels)}

#Storing configuration of matrix in json
with open("labels.json", "w") as f:
    json.dump({
        "row_labels": row_labels,
        "column_labels": column_labels
    }, f, indent=2)

#Filling matrix: Element ij of matrix is the list of bacdive codes (bacteria)
    #that share same ncbi and live same habitat
for bd_code, (taxids, habitats) in bacteria.items():
    bac_taxids_filtered = [t for t in taxids if t in r_label_to_index]
    bac_habitats_filtered = [h for h in habitats if h in c_label_to_index]

    for taxid in bac_taxids_filtered:
        for habitat in bac_habitats_filtered:
            i = r_label_to_index[taxid]
            j = c_label_to_index[habitat]
            matrix[i][j].add(bd_code)


#Creating sparse matrix  (LIst of List, easy to fill)
sparse_mat = sp.lil_matrix((len(row_labels), len(column_labels)), dtype=float)

#Filling
for ncbi, ncbi_idx in r_label_to_index.items():
    bacteria_of_ncbi = set(bd_code for sublist in matrix[ncbi_idx].values() for bd_code in sublist)    #Takes all of the bacteria that shares same "r"-ncbi
    for habitat_idx, bd_list in matrix[ncbi_idx].items():    #For takes the bacdive codes list for each column
        sparse_mat[ncbi_idx, habitat_idx] = float(len(bd_list)) / float(len(bacteria_of_ncbi))   #Element ij of sparse matrix is the number of bacteria that is_a ncbi and lives_in habitat over the total number of bacteria that is_a ncbi

#Conversion of matrix (Compressed Sparse Row, easy for computations)
sparse_mat = sparse_mat.tocsr()

#Matrix times its transpose gives us a similarity matrix ~ adj_matrix [n_ncbiXn_ncbi]
similarity_mat = sparse_mat @ sparse_mat.T

#Filtering on the highest k arcs' values of each node
k=4227
rows, cols, data = [], [], []

for i in range(similarity_mat.shape[0]):
    start = similarity_mat.indptr[i]                    #start of row
    end = similarity_mat.indptr[i+1]                    #end of row
    row_cols = similarity_mat.indices[start:end]        #indices of columns
    row_data = similarity_mat.data[start:end]           #data ordered following previous rules

    if row_data.size > k:
        top_idx = np.argpartition(row_data, -k)[-k:]    #Indices of k highest values
        row_cols = row_cols[top_idx]                    #Filtering column indices with top_indx
        row_data = row_data[top_idx]                    #Filtering related data with top_indx

    rows.extend([i] * len(row_cols))                    #For each element of i row assignes "i”
    cols.extend(row_cols)                               #Adds top-k filtered columns (indices)
    data.extend(row_data)                               #Adds top-k filtered effective values to data array

#Producing final (COO, COOrdinates) matrix
similarity_mat = sp.coo_matrix((data, (rows, cols)), shape=similarity_mat.shape)

#Saving matrix
sp.save_npz(data_files_path + "similarity_coo_mat_1.npz", similarity_mat)

with open(data_files_path + "BD_ncbi_hab_id.txt", 'w') as text:
    for name, values in list(bacteria.items()):
        ncbi_clean = ['N' + elem for elem in values[0] if elem]
        habitat_clean = ['O' + elem for elem in values[1] if elem]
        text.write(f"{name}\t{'\t'.join(ncbi_clean)}\t{'\t'.join(habitat_clean)}\n")
