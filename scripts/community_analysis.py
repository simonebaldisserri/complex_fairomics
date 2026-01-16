from ete3 import NCBITaxa
from collections import defaultdict, Counter
import json
import numpy as np
import time
import csv
import os

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAGING_DATA_DIR = os.path.join(PROJECT_DIR, "project_data/staging_data")

with open(os.path.join(STAGING_DATA_DIR, "labels_no_fungi.json"), "r") as f:
    labels_data = json.load(f)

row_labels = labels_data["row_labels"]  # NCBI taxids as strings

with open(os.path.join(STAGING_DATA_DIR, "graph_layout_data_no_fungi.json"), "r") as f:
    graph_data = json.load(f)

labels = graph_data["labels"]           # nodes used in louvain/layout script
communities = graph_data["communities"] # parallel list of community ids

ncbi = NCBITaxa()

# row_labels = list of NCBI taxids as strings
taxids_int = [int(t) for t in row_labels]

# taxid -> rank (species, genus, family, ...)
rank_dict = ncbi.get_rank(taxids_int)

# taxid (int) -> scientific name
name_dict = ncbi.get_taxid_translator(taxids_int)

# taxid (str) -> scientific name
ncbi_to_name = {str(t): name_dict.get(t, "Unknown") for t in taxids_int}


rank_counter = Counter(rank_dict.values())

print("=== Number of taxa per rank ===")
for rank, count in rank_counter.most_common():
    print(f"{rank} : {count}")


comm_to_nodes = defaultdict(list)
for node, comm_id in zip(labels, communities):
    comm_to_nodes[comm_id].append(node)

communities_sorted = sorted(comm_to_nodes.items(), key=lambda x: len(x[1]), reverse=True)
ncbi_to_communities = defaultdict(list)
community_rank_summary = {}
community_rank_full = {}


MIN_COMMUNITY_SIZE = 2

for comm_id, nodes in communities_sorted:
    if len(nodes) < MIN_COMMUNITY_SIZE:
        continue

    ranks_of_community = []
    for n in nodes:
        tid = int(n)
        ranks_of_community.append(rank_dict.get(tid, "Unknown"))
        ncbi_to_communities[str(n)].append(comm_id)
    
    counter = Counter(ranks_of_community)

    community_rank_full[comm_id] = {
        "community_size": len(nodes),
        "rank_distribution": counter
    }

    top_ranks = counter.most_common(4)

    community_rank_summary[comm_id] = {
        "size": len(nodes),
        "top_ranks": [
            {
                "rank": r,
                "count": c,
                "fraction": c / len(nodes)
            }
            for r, c in top_ranks
        ]
    }


output_file = "community_rank_distribution.csv"

with open(os.path.join(STAGING_DATA_DIR, output_file), "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "community_id",
        "community_size",
        "rank",
        "count",
        "percentage"
    ])

    for comm_id, info in community_rank_full.items():
        size = info["community_size"]
        counter = info["rank_distribution"]

        for rank, count in counter.items():
            percentage = 100.0 * count / size

            writer.writerow([
                comm_id,
                size,
                rank,
                count,
                round(percentage, 3)
            ])

print(f"\nSaved full rank distributions (size>1) to: {output_file}")


output_file = "ncbi_node_info_with_communities.csv"

with open(os.path.join(STAGING_DATA_DIR, output_file), "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["ncbi", "scientific_name", "rank", "communities"])

    for ncbi in sorted(ncbi_to_communities.keys(), key=lambda x: int(x)):
        tid_int = int(ncbi)

        name = ncbi_to_name.get(ncbi, "Unknown")
        rank = rank_dict.get(tid_int, "Unknown_rank")


        comm_list = sorted(set(ncbi_to_communities[ncbi]))
        comm_str = ";".join(map(str, comm_list))

        writer.writerow([ncbi, name, rank, comm_str])

print(f"\nSaved node info table to: {output_file}")

output_file = "summaries_of_communities.csv"

with open(os.path.join(STAGING_DATA_DIR, output_file), "w", newline="") as f:
    writer = csv.writer(f)
    # header
    writer.writerow([
        "community", "size",
        "rank_1", "count_1", "fraction_1",
        "rank_2", "count_2", "fraction_2",
        "rank_3", "count_3", "fraction_3",
        "rank_4", "count_4", "fraction_4"
    ])

    for comm_id, info in community_rank_summary.items():
        row = [comm_id, info["size"]]

        top_ranks = info["top_ranks"]

        # riempi fino a 4 (se meno di 4 rank, pad con vuoti)
        for i in range(4):
            if i < len(top_ranks):
                row.extend([
                    top_ranks[i]["rank"],
                    top_ranks[i]["count"],
                    round(top_ranks[i]["fraction"], 4)
                ])
            else:
                row.extend(["", "", ""])

        writer.writerow(row)

print(f"\nSaved communities' summeries to: {output_file}")