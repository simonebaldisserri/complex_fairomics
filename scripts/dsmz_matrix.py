from collections import defaultdict
from collections import Counter
import re
import scipy.sparse as sp
import json
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
from dsmz_processing import matrix
import os

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAGING_DATA_DIR = os.path.join(PROJECT_DIR, "project_data/staging_data")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")


with open(os.path.join(STAGING_DATA_DIR, "labels_no_fungi.json"), "r") as f:
    labels = json.load(f)

row_labels = labels["row_labels"]
column_labels = labels["column_labels"]
r_label_to_index = labels["r_label_to_index"]
c_label_to_index = labels["c_label_to_index"]

#Creating sparse matrix  (LIst of List, easy to fill)
sparse_mat = sp.lil_matrix((len(row_labels), len(column_labels)), dtype=float)

#Filling
for ncbi, ncbi_idx in r_label_to_index.items():
    for habitat_idx, bd_list in matrix[ncbi_idx].items():    #For takes the bacdive codes list for each column
        sparse_mat[ncbi_idx, habitat_idx] = len(bd_list)   #Element ij of sparse matrix is the number of bacteria that is_a ncbi and lives_in habitat


sparse_mat = sparse_mat.tocsr() #Conversion of matrix to CSR_matrix (Compressed Sparse Row, easy for computations)
sparse_mat_csc = sparse_mat.tocsc() #Usage of CSC_matrix (Compressed Sparse Column, easy for computations)
n_rows = sparse_mat.shape[0]
row_sums = np.asarray(sparse_mat.sum(axis=1)).ravel()

rows, cols, data = [], [], []

for i in range(n_rows):                         #Dissimilarity_mat for loop
    row_vals = row_sums[i] + row_sums

    # Subtract 2 * sum_k min(A[i,k], A[j,k]) for rows that share columns
    si, ei = sparse_mat.indptr[i], sparse_mat.indptr[i+1]
    cols_i = sparse_mat.indices[si:ei]
    vals_i = sparse_mat.data[si:ei]

    for kcol, vik in zip(cols_i, vals_i):
        cs, ce = sparse_mat_csc.indptr[kcol], sparse_mat_csc.indptr[kcol+1]
        c_rows = sparse_mat_csc.indices[cs:ce]
        c_vals = sparse_mat_csc.data[cs:ce]
        row_vals[c_rows] -= 2.0 * np.minimum(vik, c_vals)

    
    row_vals[i] = 0.0                               #self distance must be 0
    row_vals[np.abs(row_vals) < 1e-12] = 0.0        #tiny float noise to zero
    row_cols = np.arange(n_rows)
    row_data = row_vals

    nz = row_data != 0
    rows.extend([i]*int(nz.sum()))
    cols.extend(row_cols[nz])
    data.extend(row_data[nz])

#Creating in final (COO, COOrdinates) matrix
dissimilarity_mat = sp.coo_matrix((data, (rows, cols)), shape=(n_rows, n_rows))

similarity_mat = dissimilarity_mat.toarray()

# For Dissimilarity distribution plot
vals = similarity_mat.ravel()

plt.figure(figsize=(8,5))
plt.hist(vals, bins=np.linspace(0, 2000, 200))
plt.xlabel("Dissimilarity value")
plt.ylabel("Frequency")
plt.title("Distribution of pairwise dissimilarities")
plt.savefig(os.path.join(RESULTS_DIR, "distribuzione_dissimilarity_mat22_no_fungi.png"), dpi=300)
print("\n=== Dissimilarity statistics ===")
print("Min:", vals.min(),
      "Max:", vals.max(),
      "Mean:", vals.mean(),
      "Median:", np.median(vals))

#Parameters
lam= np.log(2)/np.median(vals)
threshold = 0.10

np.exp(-lam * similarity_mat, out=similarity_mat)
similarity_mat[similarity_mat < threshold] = 0.0
vals = similarity_mat.ravel()

#Heatmap construction
heatmap = px.imshow(similarity_mat, x=row_labels, y=row_labels, color_continuous_scale="Inferno", aspect="auto",)
heatmap.update_traces(hovertemplate="NCBI X: %{x}<br>" + "NCBI Y: %{y}<br>" + "Similarity: %{z:.5f}<extra></extra>")
heatmap.update_layout(width=900, height=900)

heatmap.write_html(os.path.join(RESULTS_DIR, "similarity_22_heatmap_no_fungi.html"))

#Sparse comeback!
similarity_mat = sp.coo_matrix(similarity_mat)

print("\n=== Similarity statistics ===")
print(f"Number of non zero values: {similarity_mat.nnz}")
print(f"Number of zero values: {(similarity_mat.shape[0] * similarity_mat.shape[0])-similarity_mat.nnz}")
print("Min:", vals.min(),
      "Max:", vals.max(),
      "Mean:", vals.mean(),
      "Median:", np.median(vals))

plt.figure(figsize=(8,5))
plt.hist(vals, bins=np.linspace(0, 1, 100))
plt.xlabel("Similarity value")
plt.ylabel("Frequency")
plt.title("Distribution of pairwise similarities (S = exp(-λD))")
plt.savefig(os.path.join(RESULTS_DIR, "distribuzione_similarity_mat22_no_fungi.png"), dpi=300)

#Saving matrix
sp.save_npz(os.path.join(STAGING_DATA_DIR, "similarity_coo_mat22_no_fungi.npz"), similarity_mat)