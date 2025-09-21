import numpy as np
import re

relation_file_path = "/home/selvaggio/Desktop/complex networks/Project/DSMZ_Habitat.txt"

#Putting all lines (splitted) in array (of arrays)
with open(relation_file_path, 'r') as text:
    info = [line.strip().split("\t") for line in text]

bacteria = {}

#Filling bacteria dict

for i in range(len(info)):
    if bacteria.get(info[i][8]):                                                                                    #Verifying if same BD already exists
        bacteria[info[i][8]][1] += re.split(r'/OBT:|/', info[i][7])[-1]       
    else:
        bacteria[info[i][8]]=[re.split(r'/ncbi:|/', info[i][3])[1:], re.split(r'/OBT:|/', info[i][7])[-1]]
#       Bacterium[BD_CODE] = [[NCBI/PATH/TO/BACTERIUM], MOST_SPECIFIC_HABITAT]
        if any(t.startswith("bd:") for t in bacteria[info[i][8]][0]):                                               #Deleting bd: specifications in NCBI path
            bacteria[info[i][8]][0] = [taxid for taxid in bacteria[info[i][8]][0] if not taxid.startswith("bd:")]

with open("/home/selvaggio/Desktop/complex networks/Project/BD_ncbi_hab_id.txt", 'w') as text:
    for name, values in list(bacteria.items()):
        ncbi_clean = [' N' + elem for elem in values[0] if elem is not None ] 
        text.write(f"{name}\t{'\t'.join(ncbi_clean)}\tO{values[1]}\n")