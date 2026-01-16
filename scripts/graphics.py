import json
import time
import colorsys
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

start_time = time.time()

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAGING_DATA_DIR = os.path.join(PROJECT_DIR, "project_data/staging_data")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")

# Load data for plots
with open(os.path.join(STAGING_DATA_DIR, "graph_layout_data_no_fungi.json"), "r") as f:
    data = json.load(f)

nodes_df = pd.read_csv(os.path.join(STAGING_DATA_DIR, "ncbi_node_info_with_communities.csv"))
rank_dist_df = pd.read_csv(os.path.join(STAGING_DATA_DIR, "community_rank_distribution.csv"))
summary_df = pd.read_csv(os.path.join(STAGING_DATA_DIR, "summaries_of_communities.csv"))
bc_cc_df = pd.read_csv(os.path.join(STAGING_DATA_DIR, "ncbi_bc_cc(2000-70%)_no_fungi.txt"), sep="\t")
degrees_i_df = pd.read_csv(os.path.join(STAGING_DATA_DIR, "degrees_i.csv"))
degrees_ii_df = pd.read_csv(os.path.join(STAGING_DATA_DIR, "degrees_ii.csv"))

# Degree histogram plot function
def plot_degree_histogram(degrees, bins, filename, xlim=None, ylim=None):
    plt.figure()
    plt.hist(degrees, bins=bins)
    plt.xlabel("Weighted Degree")
    plt.ylabel("Number of nodes")
    plt.title("Weighted Degree Distribution")

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    plt.savefig(os.path.join(RESULTS_DIR, filename), dpi=300, bbox_inches="tight")
    plt.close()

# Generate distinct colors
def distinct_colors(n):
    out = []
    step = 209
    for i in range(n):
        hue = ((i * step) % n) / n
        r, g, b = colorsys.hsv_to_rgb(hue, 1.0, 1.0)
        out.append(f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})")
    return out

labels = data["labels"]
pos3 = data["pos3"]
communities = data["communities"]
unique_comms = data["unique_comms"]
N = data["N"]
x = [pos3[n][0] for n in labels]
y = [pos3[n][1] for n in labels]
z = [pos3[n][2] for n in labels]
palette = distinct_colors(N)
color_map = {comm: palette[i] for i, comm in enumerate(unique_comms)}

# Assign RGB colors to each node
node_colors = [color_map[c] for c in communities]

# 3D SCatterPlot generation
fig = go.Figure(data=[go.Scatter3d(
    x=x, y=y, z=z,
    mode='markers',
    marker=dict(
        size=2,
        color=node_colors,
    ),
    text=labels,
    customdata=communities,
    hovertemplate=
        "<b>NCBI:</b> %{text}<br>" +
        "Community: %{customdata}<br>" +
        "<extra></extra>"
)])

fig.write_html(os.path.join(RESULTS_DIR, "louvain_3d_graph_22_no_fungi.html"))
print("Saved louvain_3d_graph_22.html")

elapsed = time.time() - start_time
print(f"\nExecution time after 3D Scatter Plot generation: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")

bc_cc_df = bc_cc_df.merge(nodes_df[["ncbi", "rank"]], on="ncbi", how="left")

# Adding unknown to missing ranks
bc_cc_df["rank"] = bc_cc_df["rank"].fillna("Unknown")

# BC vs CC scatter plot generation
fig = px.scatter(
    bc_cc_df,
    x="cc",
    y="bc",
    color="rank",
    hover_name="ncbi",
    title="Taxa ecological roles: BC vs CC",
    labels={
        "bc": "Betweenness Centrality",
        "cc": "Clustering Coefficient",
    }
)

fig.write_html(os.path.join(RESULTS_DIR, "bc_cc_scatterplot_22_no_fungi.html"))
print("Saved bc_cc_scatterplot_22_no_fungi.html")

elapsed = time.time() - start_time
print(f"\nExecution time after Betweenness centrality-Clustering coefficient Scatter Plot generation: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")

# Degree distributions plots generation
degrees_i = degrees_i_df["weighted_degree"].values
print(f"Min degree: {degrees_i.min():.2f}")
print(f"Max degree: {degrees_i.max():.2f}")
print(f"Mean degree: {degrees_i.mean():.2f}")

degrees_ii = degrees_ii_df["weighted_degree"].values
print(f"Min degree: {degrees_ii.min():.2f}")
print(f"Max degree: {degrees_ii.max():.2f}")
print(f"Mean degree: {degrees_ii.mean():.2f}")

    # Initial graph
        # 0–100
plot_degree_histogram(degrees_i, bins=np.linspace(0, 100, 100), filename="degree_histogram_i_100.png")

        # 2000–4000
plot_degree_histogram(degrees_i, bins=np.linspace(2000, 4000, 900), filename="degree_histogram_i_2000.png", ylim=(0, 130))

        # 0–4000
plot_degree_histogram(degrees_i, bins=np.linspace(0, 4000, 1100), filename="degree_histogram_i.png", ylim=(0, 130))

    # Graph after further Sparsification
        # 0–350
plot_degree_histogram(degrees_ii, bins=np.linspace(0, 250, 100), filename="degree_histogram_ii_250.png")

        # 1500–1810
plot_degree_histogram(degrees_ii, bins=np.linspace(1500, 1810, 310), filename="degree_histogram_ii_1500.png", ylim=(0, 130))

        # 0–1810
plot_degree_histogram(degrees_ii, bins=np.linspace(0, 1810, 880), filename="degree_histogram_ii.png", ylim=(0, 130))


# Filtering target ranks for rank distribution
KEEP_RANKS = ["species", "genus", "family", "order", "class", "phylum", "subspecies", "Unknown", "no rank", "kingdom"]

nodes_df["rank_norm"] = nodes_df["rank"].apply(lambda r: r if r in KEEP_RANKS else "Unknown")

# Global deepest rank distribution plot generation
global_rank_dist = (nodes_df["rank_norm"].value_counts())
plt.figure(figsize=(10, 5))
global_rank_dist.plot(kind="bar", edgecolor="k")
plt.ylabel("Number of taxa")
plt.xlabel("Rank")
plt.title("Global distribution of deepest taxonomic ranks")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "global_deepest_rank_distribution.png"), dpi=300)
plt.close()
print("Saved global_deepest_rank_distribution.png")

elapsed = time.time() - start_time
print(f"\nExecution time after Global deepest rank distribution plot generation: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")


# Filtering target ranks
RANKS = ["species", "genus", "family", "order", "class", "phylum"]

rank_dist_df = rank_dist_df[rank_dist_df["rank"].isin(RANKS)]

# Size vs #Rank Scatter Plot
for rank in RANKS:
    df = rank_dist_df[rank_dist_df["rank"] == rank]
    plt.figure(figsize=(6, 4))
    plt.scatter(
        df["community_size"],
        df["count"],
        color="red",
        alpha=0.3
    )
    plt.xlabel("Community size")
    plt.ylabel(f"{rank.capitalize()} count")
    plt.title(f"Community size vs {rank.capitalize()} count")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, f"Comm_size_vs_{rank.capitalize()}.png"), dpi=300, bbox_inches="tight")
    plt.close()
print("Saved Comm_size_vs_\"rank\" plots")

elapsed = time.time() - start_time
print(f"\nExecution time after Comm_size_vs_\"rank\" plots generation: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")

# Taxonomic Dispersion (BOXPLOT)
dispersion_data = []

for rank in RANKS:
    df = rank_dist_df[rank_dist_df["rank"] == rank].copy()
    df["dispersion"] = df["count"] / df["community_size"]
    df["rank"] = rank
    dispersion_data.append(df[["rank", "dispersion"]])

disp_df = pd.concat(dispersion_data)

plt.figure(figsize=(7, 5))
disp_df.boxplot(
    column="dispersion",
    by="rank",
    grid=False
)
plt.suptitle("")
plt.title("Taxonomic dispersion per rank")
plt.xlabel("Rank")
plt.ylabel("Dispersion (#taxa / community size)")

plt.savefig(os.path.join(RESULTS_DIR, "taxonomic_dispersion_boxplot_no_fungi.png"), dpi=300)
plt.close()
print("Saved taxonomic_dispersion_boxplot_no_fungi.png")

elapsed = time.time() - start_time
print(f"\nExecution time after Taxonomic Dispersion Boxplot generation: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")

# Rank composition Stacked Bar Plot (top 10 communities)
top10 = (
    rank_dist_df[["community_id", "community_size"]]
    .drop_duplicates()
    .sort_values("community_size", ascending=False)
    .head(10)["community_id"]
    .tolist()
)

stack_df = rank_dist_df[
    rank_dist_df["community_id"].isin(top10)
]

pivot = stack_df.pivot_table(
    index="community_id",
    columns="rank",
    values="percentage",
    fill_value=0
)
    # Consistent order
pivot = pivot[RANKS] 

pivot.plot(
    kind="bar",
    stacked=True,
    figsize=(10, 6),
    colormap="tab20"
)

plt.xlabel("Community ID")
plt.ylabel("Percentage of taxa")
plt.title("Rank composition of top 10 communities")
plt.legend(title="Rank", bbox_to_anchor=(1.02, 1), loc="upper left")

plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "stacked_rank_composition_top10_no_fungi.png"), dpi=300)
plt.close()
print("Saved stacked_rank_composition_top10_no_fungi.png")

elapsed = time.time() - start_time
print(f"\nTotal execution time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")