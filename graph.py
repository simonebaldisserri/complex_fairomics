import networkx as nx
import json
import scipy.sparse as sp
from tqdm import tqdm 

data_files_path = "/home/selvaggio/Desktop/complex networks/Project/"

similarity_coo_mat = sp.load_npz(data_files_path + "similarity_coo_mat.npz")

with open(data_files_path + "labels.json", "r") as f:
    labels = json.load(f)

row_labels = labels["row_labels"]
column_labels = labels["column_labels"]

print(type(similarity_coo_mat))
print("Numero elementi non zero:", similarity_coo_mat.nnz)

similarity_graph = nx.Graph()
for i, j, w in tqdm(zip(similarity_coo_mat.row, similarity_coo_mat.col, similarity_coo_mat.data), desc = "Costruzione Graph"):
    if w > 0:
        similarity_graph.add_edge(row_labels[i], row_labels[j], weight=w)

