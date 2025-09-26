import numpy as np
import re

data_files_path = "/home/selvaggio/Desktop/complex networks/Project/"

#Putting all lines (splitted) in array (of arrays)
with open(data_files_path + "DSMZ_Habitat.txt", 'r') as text:
    info = [line.strip().split("\t") for line in text]

bacteria = {}
all_ncbi = set()
all_habitats = set()

#Filling bacteria dict

for i in range(len(info)):
    if bacteria.get(info[i][8]):                                                                                    #Verifying if same BD already exists
        bacteria[info[i][8]][1] += re.split(r'/OBT:|/', info[i][7])[-1]
        all_habitats.add(re.split(r'/OBT:|/', info[i][7])[-1])    
    else:
        if any(t.startswith("bd:") for t in re.split(r'/ncbi:|/', info[i][3])[1:]):                                               #Deleting bd: specifications in NCBI path
            bacteria[info[i][8]] = [re.split(r'/ncbi:|/', info[i][3])[1:-1], re.split(r'/OBT:|/', info[i][7])[-1]]
            all_ncbi.update(re.split(r'/ncbi:|/', info[i][3])[1:-1])
            all_habitats.add(re.split(r'/OBT:|/', info[i][7])[-1])
        else:
            bacteria[info[i][8]]=[re.split(r'/ncbi:|/', info[i][3])[1:], re.split(r'/OBT:|/', info[i][7])[-1]]
            all_ncbi.update(re.split(r'/ncbi:|/', info[i][3])[1:])
            all_habitats.add(re.split(r'/OBT:|/', info[i][7])[-1])
#       Bacterium[BD_CODE] = [[NCBI/PATH/TO/BACTERIUM], MOST_SPECIFIC_HABITAT]

with open(data_files_path + "BD_ncbi_hab_id.txt", 'w') as text:
    for name, values in list(bacteria.items()):
        ncbi_clean = [' N' + elem for elem in values[0] if elem is not None ] 
        text.write(f"{name}\t{'\t'.join(ncbi_clean)}\tO{values[1]}\n")


print(len(all_habitats))
print(len(all_ncbi))