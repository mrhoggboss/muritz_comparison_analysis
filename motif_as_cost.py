import networkx as nx
import pickle
import pandas as pd
import os 
import numpy as np
import ot
import csv
from scipy.stats import pearsonr
webs_pkl_path = "./data/loc2webs.pkl"

with open(webs_pkl_path, "rb") as f: 
    loc2webs = pickle.load(f)

MOTIF_PATH = "./data/motif_vectors"

def compute_optimal_transport(loc1, loc2, motif_dist = 'euclidean', marginal = "uniform"):
    G1 = nx.convert_node_labels_to_integers(loc2webs[loc1], label_attribute="name")
    species_labels1 = list(nx.get_node_attributes(G1, "name").values())
    G2 = nx.convert_node_labels_to_integers(loc2webs[loc2], label_attribute="name")
    species_labels2 = list(nx.get_node_attributes(G2, "name").values())
    motif_G1 = pd.read_csv(os.path.join(MOTIF_PATH, f"{loc1}_motifs.csv"))
    motif_G2 = pd.read_csv(os.path.join(MOTIF_PATH, f"{loc2}_motifs.csv"))
    
    # compute cost matrix
    cost_matrix = np.zeros((len(species_labels1), len(species_labels2)))
    for row_idx, species1 in enumerate(species_labels1):
        motif1 = motif_G1.loc[motif_G1["node"] == species1].to_numpy()[0][1:]
        # print(motif1)
        for col_idx, species2 in enumerate(species_labels2):
            motif2 = motif_G2.loc[motif_G2["node"] == species2].to_numpy()[0][1:]
            # print(motif2)

            if motif_dist == "euclidean":
                cost = np.linalg.norm(motif1 - motif2)
            elif motif_dist == "Hamming":
                cost = np.count_nonzero(motif1!=motif2)
            elif motif_dist == "PearsonR":
                # print(np.array(motif1, dtype=float)), print(motif2)
                cost, _ = pearsonr(np.array(motif1, dtype=float), np.array(motif2, dtype=float))
                cost = 1 - cost
                print(f"Cost between {species1} and {species2} is {cost}")
            else:
                print("motif distance measure not recognized. Euclidean is used.")
                cost = np.linalg.norm(motif1 - motif2)
            
            cost_matrix[row_idx][col_idx] = cost
    
    # compute marginal distrbution
    if marginal == "uniform":
        marg1 = np.ones(len(species_labels1)) / len(species_labels1)
        marg2 = np.ones(len(species_labels2)) / len(species_labels2)
    else:
        print("marginal calculation method not recognized. Uniform is used.")
        marg1 = np.ones(len(species_labels1)) / len(species_labels1)
        marg2 = np.ones(len(species_labels2)) / len(species_labels2)
    
    # compute transport plan and wasserstein distance
    transport_plan = ot.emd(marg1, marg2, cost_matrix)
    wass_distance = ot.emd2(marg1, marg2, cost_matrix)

    return transport_plan, wass_distance, species_labels1, species_labels2

# file io
# hyperparameters
MOTIF_DIST = "PearsonR"
TRANSPORT_PLAN_DIR = f"./data/motif_as_cost_matrix_transport_plans_{MOTIF_DIST}"
TRANSPORT_COST_DIR = "./data/alignment_scores"
MARGINAL = "uniform"

done_pairs = set()  # Use a set for efficient lookups
transport_cost_matrix = np.zeros((len(loc2webs), len(loc2webs)))

if not os.path.exists(TRANSPORT_PLAN_DIR):
    os.makedirs(TRANSPORT_PLAN_DIR)

for idx1, loc1 in enumerate(loc2webs):
    for idx2, loc2 in enumerate(loc2webs):
        if idx1 < idx2 and (loc1, loc2) not in done_pairs and (loc2, loc1) not in done_pairs:
            # Compute transport plan and cost
            transport_plan, transport_cost, species_labels1, species_labels2 = compute_optimal_transport(loc1, loc2, MOTIF_DIST, MARGINAL)

            # Save transport plan
            transport_plan_filename = f"{loc1}_vs_{loc2}_transport_plan.csv"
            with open(os.path.join(TRANSPORT_PLAN_DIR, transport_plan_filename), mode='w', newline='') as file:
                writer = csv.writer(file)
                # Write column headers (with an empty space for row header column)
                writer.writerow([''] + species_labels2)
                # Write each row with its corresponding row header
                for row_header, row in zip(species_labels1, transport_plan):
                    writer.writerow([row_header] + list(row))

            print(f"Transport plan saved to {transport_plan_filename}")

            # Update transport cost matrix
            transport_cost_matrix[idx1][idx2] = transport_cost
            transport_cost_matrix[idx2][idx1] = transport_cost  # Assume symmetry if applicable
            print(f"Transport cost: {transport_cost}")

            # Record done pairs
            done_pairs.add((loc1, loc2))

# Save transport cost matrix
transport_cost_filename = f"motif_as_cost_matrix_transport_cost_{MOTIF_DIST}.csv"
with open(os.path.join(TRANSPORT_COST_DIR, transport_cost_filename), mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write column headers (with an empty space for row header column)
    writer.writerow([''] + list(loc2webs))
    # Write each row with its corresponding row header
    for row_header, row in zip(list(loc2webs), transport_cost_matrix):
        writer.writerow([row_header] + list(row))



