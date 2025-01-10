
import numpy as np
import os
import pickle
import networkx as nx
import pandas as pd

webs_pkl_path = "./data/loc2webs.pkl"

with open(webs_pkl_path, "rb") as f: 
    loc2webs = pickle.load(f)

locs = list(loc2webs)

def compute_ot_cost(G1, G2, alignment_matrix, nodelist1, nodelist2):
    # Retrieve cost matrices and node distributions
    C1 = nx.floyd_warshall_numpy(G1, nodelist = nodelist1)
    C1[C1 == np.inf] = 1000
    C2 = nx.floyd_warshall_numpy(G2, nodelist = nodelist2)
    C2[C2 == np.inf] = 1000

    # Assume T is the alignment matrix (optimal transport plan)
    T = alignment_matrix

    # Calculate the optimal transport cost
    def compute_gw_cost(C1, C2, T):
        """
        Compute the Gromov-Wasserstein cost.

        Parameters:
        -----------
        C1 : ndarray
            Cost matrix for graph G1 (size n1 x n1).
        C2 : ndarray
            Cost matrix for graph G2 (size n2 x n2).
        T : ndarray
            Transport plan between nodes of G1 and G2 (size n1 x n2).

        Returns:
        --------
        float
            The Gromov-Wasserstein cost.
        """
        n1, n2 = T.shape
        gw_cost = 0.0

        # Double loop over all pairs (i, j) in G1 and (k, l) in G2
        for i in range(n1):
            for j in range(n1):
                for k in range(n2):
                    for l in range(n2):
                        # Compute the contribution for this combination
                        diff = (C1[i, j] - C2[k, l]) ** 2
                        gw_cost += diff * T[i, k] * T[j, l]
        
        return gw_cost

    return compute_gw_cost(C1, C2, alignment_matrix)

def compute_ot_cost_wrapper(loc1, loc2, input_filepath):
    # Handle special case for corrupted location names
    if loc1 == "Qui\uf03f\uf03fma":
        loc1 = "Quiçãma"
    if loc2 == "Qui\uf03f\uf03fma":
        loc2 = "Quiçãma"
    row_idx, col_idx = locs.index(loc1), locs.index(loc2)
    # if row_idx <= col_idx:  # Only fill lower triangle
    #     row_idx, col_idx = col_idx, row_idx
    
    # read alignment matrix and nodelists
    # Load the CSV into a Pandas DataFrame
    df = pd.read_csv(input_filepath, index_col=0)  # Use the first column as row headers

    # Extract the node lists
    nodelist1 = df.index.tolist()  # Column headers (graph 1 nodes)
    # print(nodelist1)
    # print(loc2webs[loc1].nodes())
    nodelist2 = df.columns.tolist()  # Row headers (graph 2 nodes)
    # print(nodelist2)
    # print(loc2webs[loc2].nodes())
    # Extract the alignment matrix as a NumPy array
    alignment_matrix = df.to_numpy()
    return compute_ot_cost(loc2webs[loc1], loc2webs[loc2], alignment_matrix, nodelist1, nodelist2), row_idx, col_idx

# file i/o 
input_dir = "./data/kai's_transport_plans_fgw"
output_dir = "./data/alignment_scores"
output_filename = "kai_alignment_quality_fgw_[checking].csv"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

ot_cost_matrix = np.zeros((len(loc2webs), len(loc2webs)))
locs = list(loc2webs)
# Process each file in the input directory 
for filename in os.listdir(input_dir):
    # # - muritz
    # if filename.lower().endswith('_alignment.csv'):
    #     input_filepath = os.path.join(input_dir, filename)

    #     core_name = filename.replace("_alignment.csv", "")
    #     locations = core_name.split("_vs_")
    #     loc1, loc2 = locations[0].replace("_", " "), locations[1].replace("_", " ")

    #     ot_cost, row_idx, col_idx = compute_ot_cost_wrapper(loc1, loc2, input_filepath)
    #     if row_idx <= col_idx:  # Only fill lower triangle
    #         row_idx, col_idx = col_idx, row_idx
    #     ot_cost_matrix[row_idx, col_idx] = ot_cost
    #     print(f"Calculated OT cost between {loc1} and {loc2}: {ot_cost}")

    # - kai
    if filename.lower().endswith('.csv'):
        input_filepath = os.path.join(input_dir, filename)

        core_name = filename.replace(".csv", "")
        locations = core_name.split("_")
        loc1, loc2 = locations[0], locations[1]

        ot_cost, row_idx, col_idx = compute_ot_cost_wrapper(loc1, loc2, input_filepath)
        if row_idx <= col_idx:  # Only fill lower triangle
            row_idx, col_idx = col_idx, row_idx
        ot_cost_matrix[row_idx, col_idx] = ot_cost
        print(f"Calculated OT cost between {loc1} and {loc2}: {ot_cost}")
ot_cost_matrix = ot_cost_matrix + ot_cost_matrix.T
output_filepath = os.path.join(output_dir, output_filename)
np.savetxt(output_filepath, ot_cost_matrix, delimiter=",", fmt="%.10e")

# # testing code
# print(compute_ot_cost_wrapper("Aberdare", "Ai Ais", "./data/muritz_alignment_matrices/Aberdare_vs_Ai_Ais_alignment.csv"))
# print(compute_ot_cost_wrapper("Aberdare", "Ai Ais", "./data/kai's_transport_plans_gw/Aberdare_Ai Ais.csv"))
# print(compute_ot_cost_wrapper("Aberdare", "Akagera", "./data/kai's_transport_plans_gw/Aberdare_Akagera.csv"))
# print(compute_ot_cost_wrapper("Aberdare", "Amboseli", "./data/kai's_transport_plans_gw/Aberdare_Amboseli.csv"))

