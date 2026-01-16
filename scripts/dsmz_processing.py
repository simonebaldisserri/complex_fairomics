from collections import defaultdict
from collections import Counter
import re
import json
import numpy as np
import os

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DATA_DIR = os.path.join(PROJECT_DIR, "project_data/initial_data")
STAGING_DATA_DIR = os.path.join(PROJECT_DIR, "project_data/staging_data")

#Putting all lines (splitted) in array (of arrays)
with open(os.path.join(RAW_DATA_DIR, "DSMZ_Habitat.txt"), 'r') as text:
    info = [line.strip().split("\t") for line in text]

bacteria = {}

#Filling bacteria dict as: Bacterium[BD_CODE] = [[NCBI/PATH/TO/BACTERIUM], [LAST/FIVE/HABITATS/OF/PATHS]]
for i in range(len(info)):
    bd_code = info[i][8]
    ncbi_path = [t for t in re.split(r'/ncbi:|/', info[i][3])[1:] if t and not t.startswith("bd:")]
    habitat_paths = info[i][7].split(',')
    habitat_path = set()
    for hp in habitat_paths:
        habitat_path.update([h for h in re.split(r'/OBT:|/', hp)[2:] if h])
    if bd_code in bacteria:                                                                         #if bacterium already exists it adds new habitats
        bacteria[bd_code][1].update(habitat_path)
    else:
        bacteria[bd_code] = [ncbi_path, habitat_path]
for bd_code in bacteria:
    bacteria[bd_code][1] = list(bacteria[bd_code][1])

bacteria = {                            #removing fungi bacteria
    bd_code: value
    for bd_code, value in bacteria.items()
    if "4751" not in value[0]
}


#Habitats filter >3
habitat_counter = Counter()
for _, habitats in bacteria.values():
    habitat_counter.update(habitats)

excluded_habitats = {"000001", "000006", "000009", "000010", "000013", "000014", "000039", "000047", "000089", "000158", "000193", "000490"}       
min_bacteria_per_habitat = 2
filtered_habitats = {habitat for habitat, count in habitat_counter.items() if count >= min_bacteria_per_habitat and habitat not in excluded_habitats}

habitat_counter = Counter() #Re-initializing of counter
for _, habitats in bacteria.values():
    for h in habitats:
        if h in filtered_habitats:
            habitat_counter[h] += 1

#Taxonomy ID filter >3
ncbi_counter = Counter()
for ncbi, _ in bacteria.values():
    ncbi_counter.update(ncbi)

excluded_ncbi = {"1"}
min_bacteria_per_ncbi = 2
filtered_ncbi = {ncbi for ncbi, count in ncbi_counter.items() if count >= min_bacteria_per_ncbi and ncbi not in excluded_ncbi}

ncbi_counter = Counter() #Re-initializing of counter
for ncbi, _ in bacteria.values():
    for n in ncbi:
        if n in filtered_ncbi:
            ncbi_counter[n] += 1

#Creating matrix
matrix = defaultdict(lambda: defaultdict(set))          #this structure automatically creates the dict whenever we call a key that doesn't exist
row_labels = list(filtered_ncbi)
column_labels = list(filtered_habitats)
r_label_to_index = {label: idx for idx, label in enumerate(row_labels)}
c_label_to_index = {label: idx for idx, label in enumerate(column_labels)}

#Storing configuration of matrix in json
with open(os.path.join(STAGING_DATA_DIR, "labels_no_fungi.json"), "w") as f:
    json.dump({
        "row_labels": row_labels,
        "column_labels": column_labels,
        "r_label_to_index": r_label_to_index,
        "c_label_to_index": c_label_to_index
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