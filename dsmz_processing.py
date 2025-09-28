from collections import defaultdict
import re
import scipy.sparse as sp
import json

data_files_path = "/home/selvaggio/Desktop/complex networks/Project/"

#Putting all lines (splitted) in array (of arrays)
with open(data_files_path + "DSMZ_Habitat.txt", 'r') as text:
    info = [line.strip().split("\t") for line in text]

bacteria = {}
all_ncbi = set()
all_habitats = set()

#Filling bacteria dict
#and set() of taxID and habitats
for i in range(len(info)):
    if bacteria.get(info[i][8]) and bacteria[info[i][8]][1]!=re.split(r'/OBT:|/', info[i][7])[-1]: #if bacterium already exists
        bacteria[info[i][8]][1] += re.split(r'/OBT:|/', info[i][7])[-1]
        all_habitats.add(re.split(r'/OBT:|/', info[i][7])[-1])    
    else:
        if any(t.startswith("bd:") for t in re.split(r'/ncbi:|/', info[i][3])[1:]):                                               #Deleting bd: specifications in NCBI path
            bacteria[info[i][8]] = [re.split(r'/ncbi:|/', info[i][3])[1:], [re.split(r'/OBT:|/', info[i][7])[-1]]]
            all_ncbi.update(re.split(r'/ncbi:|/', info[i][3])[1:-1])
            all_habitats.add(re.split(r'/OBT:|/', info[i][7])[-1])
        else:
            bacteria[info[i][8]]=[re.split(r'/ncbi:|/', info[i][3])[1:], re.split(r'/OBT:|/', info[i][7])[-1]]
            all_ncbi.update(re.split(r'/ncbi:|/', info[i][3])[1:])
            all_habitats.add(re.split(r'/OBT:|/', info[i][7])[-1])
#Bacterium[BD_CODE] = [[NCBI/PATH/TO/BACTERIUM], MOST_SPECIFIC_HABITAT]

#Creating matrix
matrix = defaultdict(lambda: defaultdict(list))
row_labels = list(all_ncbi)
column_labels = list(all_habitats)
r_label_to_index = {label: idx for idx, label in enumerate(row_labels)}
c_label_to_index = {label: idx for idx, label in enumerate(column_labels)}
index_to_r_label = dict(enumerate(row_labels))

#Storing configuration of matrix in json
with open("labels.json", "w") as f:
    json.dump({
        "row_labels": row_labels,
        "column_labels": column_labels
    }, f, indent=2)

#Creating storing info matrix
for bd_code, (taxids, habitats) in bacteria.items():
    for taxid in filter(None, taxids):
        for habitat in filter(None, habitats):
            i = r_label_to_index.get(taxid)
            j = c_label_to_index.get(habitat)
            if i is not None and j is not None:
                matrix[i][j].append(bd_code)    #Element ij of sparse matrix is the list of bacdive codes (bacteria) that share same ncbi and live same habitat

#Creating sparse matrix  (LIst of List, easy to fill)
sparse_mat = sp.lil_matrix((len(row_labels), len(column_labels)), dtype=float)

#Filling
for ncbi, ncbi_idx in r_label_to_index.items():
    bacteria_of_ncbi = set(bd_code for sublist in matrix[ncbi_idx].values() for bd_code in sublist)    #Takes all of the bacteria that shares same "r"-ncbi
    for habitat_idx, bd_list in matrix[ncbi_idx].items():    #For takes the bacdive codes list for each column
        sparse_mat[ncbi_idx, habitat_idx] = float(len(bd_list)) / float(len(bacteria_of_ncbi))   #Element ij of sparse matrix is the number of bacteria that is_a ncbi and lives_in habitat over the total number of bacteria that is_a ncbi

#Conversion of matrix (Compressed Sparse Row, easy to compute with)
sparse_mat = sparse_mat.tocsr()

similarity_mat = sparse_mat @ sparse_mat.T

#Conversion to final Matrix (COOrdinate list, better to store)
similarity_coo_mat = similarity_mat.tocoo()
sp.save_npz(data_files_path + "similarity_coo_mat.npz", similarity_coo_mat)

# TO LOAD IT USE: similarity_coo_mat = load_npz(data_files_path + "similarity_coo_mat.npz")

with open(data_files_path + "BD_ncbi_hab_id.txt", 'w') as text:
    for name, values in list(bacteria.items()):
        ncbi_clean = [' N' + elem for elem in values[0] if elem is not None ] 
        text.write(f"{name}\t{'\t'.join(ncbi_clean)}\tO{values[1]}\n")