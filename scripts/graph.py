import networkx as nx
import json
import scipy.sparse as sp
from tqdm import tqdm
import community.community_louvain as community_louvain
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

start_time = time.time()

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAGING_DATA_DIR = os.path.join(PROJECT_DIR, "project_data/staging_data")

similarity_coo_mat = sp.load_npz(os.path.join(STAGING_DATA_DIR, "similarity_coo_mat22_no_fungi.npz"))

with open(os.path.join(STAGING_DATA_DIR,  "labels_no_fungi.json"), "r") as f:
    labels = json.load(f)

row_labels = labels["row_labels"]
column_labels = labels["column_labels"]

similarity_graph = nx.Graph()
for i, j, w in tqdm(zip(similarity_coo_mat.row, similarity_coo_mat.col, similarity_coo_mat.data), desc = "Costruzione Graph"):
    if w > 0:
        similarity_graph.add_edge(row_labels[i], row_labels[j], weight=w)
print("Number of non-zero elements of the graph:", similarity_coo_mat.nnz)

end_time = time.time() 
elapsed = end_time - start_time
print(f"Execution time after nx.graph construction: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")

# Degree computation
degrees = [degree for node, degree in similarity_graph.degree(weight="weight")]

df_degrees = pd.DataFrame({"weighted_degree": degrees})
df_degrees.to_csv(os.path.join(STAGING_DATA_DIR, "degrees_i.csv"), index=False)

print("\nSaved degrees_i.csv")
end_time = time.time() 
elapsed = end_time - start_time
print(f"Execution time after saving degrees_i: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")


# 3D spring layout
pos3 = nx.spring_layout(similarity_graph, dim=3, weight="weight")

labels = list(similarity_graph.nodes())
x = [pos3[n][0] for n in labels]
y = [pos3[n][1] for n in labels]
z = [pos3[n][2] for n in labels]

print(f"\nSpring layout completed")
end_time = time.time() 
elapsed = end_time - start_time
print(f"Execution time after 3D nx.spring_layout: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")

    # Louvain community detection
partition = community_louvain.best_partition(similarity_graph, weight="weight")
communities = [partition[n] for n in labels]
    # Number of communities
unique_comms = sorted(set(communities))
N = len(unique_comms)

print(f"\nFound {N} communities")
pos3_list = {node: pos3[node].tolist() for node in labels}

output = {"labels": labels, "pos3": pos3_list, "communities": communities, "unique_comms": list(unique_comms), "N": N}

with open(os.path.join(STAGING_DATA_DIR, "graph_layout_data_no_fungi.json"), "w") as f:
    json.dump(output, f)
    
print("Saved graph_layout_data.json")
end_time = time.time() 
elapsed = end_time - start_time
print(f"Execution time after community detection: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")


#Betweenness Centrality
    # Number of sample nodes
k = 800  


b_c_approx = nx.betweenness_centrality(similarity_graph, k=k, weight="weight", seed=42)

sorted_bc_nodes = sorted(b_c_approx.items(), key = lambda x: x[1], reverse=True)

end_time = time.time() 
elapsed = end_time - start_time
print(f"\nExecution time after Betweenness Centrality: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")


#Further Sparsification for Clustering Coefficient
    # Percentage of edges to keep for each node
PERCENT = 0.70         
    # Max number of neighbors for each node
M_MAX = 2000           

keep_edges = {}

for node in tqdm(list(similarity_graph.nodes()), desc="Top-percent selection"):

    neighbors = similarity_graph[node]
    weights = [(nbr, data["weight"]) for nbr, data in neighbors.items()]
    n_of_edges = len(weights)

    if n_of_edges == 0:
        keep_edges[node] = set()
        continue

    # Number of neighbors to keep
    M_node = int(PERCENT * n_of_edges)
    if M_node > M_MAX:
        M_node = M_MAX
    if M_node < 1:
        M_node = 1

    # Keep all neighbors if less than M_Max
    if n_of_edges <= M_node:
        keep_edges[node] = set(n for n, w in weights)
        continue

    # ordering by decreasing weight
    weights_sorted = sorted(weights, key=lambda x: x[1], reverse=True)

    # keep only top M_node
    keep_edges[node] = set(n for n, w in weights_sorted[:M_node])

for u, v in tqdm(list(similarity_graph.edges()), desc="Cutting arcs"):
    # keep arc only if both nodes want the neighbor
    if not (v in keep_edges[u] and u in keep_edges[v]):
        similarity_graph.remove_edge(u, v)

end_time = time.time() 
elapsed = end_time - start_time
print(f"\nExecution time after sparsification: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
print(f"Number of non-zero elements after sparsification:", similarity_graph.number_of_edges())


# New Degree computation
degrees = [degree for node, degree in similarity_graph.degree(weight="weight")]

df_degrees = pd.DataFrame({"weighted_degree": degrees})
df_degrees.to_csv(os.path.join(STAGING_DATA_DIR, "degrees_ii.csv"), index=False)

print("\nSaved degrees_ii.csv")
end_time = time.time() 
elapsed = end_time - start_time
print(f"Execution time after saving degrees_ii: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")


# Clustering Coefficient
clustering_coeffs = nx.clustering(similarity_graph, weight='weight')

    # ordering by decreasing clustering coefficient values
sorted_cc_nodes = sorted(clustering_coeffs.items(), key=lambda x: x[1], reverse=True)

end_time = time.time()
elapsed = end_time - start_time
print(f"\nExecution time after Clustering Coefficient: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")

# Saving BC-CC info for each node
bc_dict = dict(sorted_bc_nodes)
cc_dict = dict(sorted_cc_nodes)

all_ncbi = set(bc_dict.keys()) | set(cc_dict.keys())

with open(os.path.join(STAGING_DATA_DIR, "ncbi_bc_cc(2000-70%)_no_fungi.txt"), "w") as f:
    f.write("ncbi\tbc\tcc\n")
    for ncbi in sorted(all_ncbi):
        bc_val = bc_dict.get(ncbi, 0.0)
        cc_val = cc_dict.get(ncbi, 0.0)
        f.write(f"{ncbi}\t{bc_val:.25f}\t{cc_val:.25f}\n")

print(f'\nSaved ncbi_bc_cc(2000-70%)_no_fungi.txt')
end_time = time.time()
elapsed = end_time - start_time
print(f"\nTotal execution time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")