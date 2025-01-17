import csv
import numpy as np
from scipy.stats import pearsonr
import os
import pickle

webs_pkl_path = "./data/undir_loc2webs.pkl"

with open(webs_pkl_path, "rb") as f: 
    loc2webs = pickle.load(f)

def load_csv(filepath):
    with open(filepath, 'r') as file:
        reader = csv.reader(file)
        data = list(reader)
    return data

def parse_alignment_matrix(filepath):
    data = load_csv(filepath)
    species_col = data[0][1:]  # Column headers (species from network 2)
    species_row = [row[0] for row in data[1:]]  # Row headers (species from network 1)
    matrix = np.array([[float(value) for value in row[1:]] for row in data[1:]])
    return matrix, species_row, species_col

def parse_motif_vectors(filepath):
    data = load_csv(filepath)
    species = [row[0] for row in data[1:]]
    vectors = {row[0]: np.array([float(value) for value in row[1:]]) for row in data[1:]}
    return species, vectors

def calculate_zeroth_soft_alignment_quality(alignment_filepath, motif_filepath_1, motif_filepath_2):
    """
    Calculate zeroth-degree alignment quality for soft alignments.
        
    Returns:
        float: zeroth-degree alignment quality.
    """
    weighted_cost = 0
    total_weight = 0

    alignment_matrix, species_row, species_col = parse_alignment_matrix(alignment_filepath)
    # G1species2index = {species:species_row.index(species) for species in species_row}
    # G2species2index = {species:species_col.index(species) for species in species_col}
    _, motif_vectors_1 = parse_motif_vectors(motif_filepath_1)
    _, motif_vectors_2 = parse_motif_vectors(motif_filepath_2)

    # Iterate through all pairs in the alignment matrix
    # print(f"{len(motif_vectors_1)} motifs for G1")
    # print(motif_vectors_1)
    # print(f"{len(motif_vectors_2)} motifs for G2")
    # print(motif_vectors_2)
    # print(f"{alignment_matrix.shape} sized alignment matrix")
    for species_i in motif_vectors_1.keys():
        for species_j in motif_vectors_2.keys():
            weight = alignment_matrix[species_row.index(species_i), species_col.index(species_j)]
            if weight > 0:  # Consider non-zero weights
                # Retrieve motif vectors
                vector_i = motif_vectors_1[species_i]
                vector_j = motif_vectors_2[species_j]
                # Compute Pearson correlation
                rho, _ = pearsonr(vector_i, vector_j)
                # print(1 - rho)
                # Accumulate weighted cost
                weighted_cost += weight * (1 - rho)
                total_weight += weight

    # Normalize the cost
    if total_weight == 0:  # Avoid division by zero
        print("zero alignment matrix")
        return 0
    print("weighted sum is: " + str(weighted_cost))
    print("total weight is:" + str(total_weight))
    e_norm = weighted_cost / total_weight
    return e_norm

# def calculate_first_soft_alignment_quality(loc1, loc2, alignment_filepath, motif_filepath_1, motif_filepath_2):
#     """
#     Calculate first-degree alignment quality for soft alignments.
        
#     Returns:
#         float: first-degree alignment quality.
#     """
#     G1 = loc2webs[loc1]
#     G2 = loc2webs[loc2]
#     weighted_cost = 0
#     total_weight = 0

#     alignment_matrix, species_row, species_col = parse_alignment_matrix(alignment_filepath)
#     _, motif_vectors_1 = parse_motif_vectors(motif_filepath_1)
#     _, motif_vectors_2 = parse_motif_vectors(motif_filepath_2)

#     # Iterate through all pairs in the alignment matrix
#     for i, species_i in enumerate(motif_vectors_1.keys()):
#         for j, species_j in enumerate(motif_vectors_2.keys()):
#             weight = alignment_matrix[i, j]

#             weighted_cost = 0
#             num_of_neighbors = 0
#             alpha_paired = {species: False for species in G1.neighbors(species_i)}
#             beta_paired = {species: False for species in G2.neighbors(species_j)}

#             for species_alpha in G1.neighbors(species_i):
#                 index_alpha = index2species[species_alpha]
                
#                 for species_beta in G2.neighbors(species_j):
#                     index_beta = index2species[species_beta]
#                     neighbor_specific_weight = alignment_matrix[index_alpha, index_beta]
#                     if neighbor_specific_weight > 0:
#                         num_of_neighbors += 1
#                         # Retrieve motif vectors
#                         vector_i = motif_vectors_1[species_i]
#                         vector_j = motif_vectors_2[species_j]
#                         # Compute Pearson correlation
#                         rho, _ = pearsonr(vector_i, vector_j)
#                         # Accumulate weighted cost
#                         weighted_cost += neighbor_specific_weight * (1 - rho)
#                         total_weight += neighbor_specific_weight
#                     else:
                



#     # Normalize the cost
#     if total_weight == 0:  # Avoid division by zero
#         print("zero alignment matrix")
#         return 0
        
#     e_norm = weighted_cost / total_weight
#     return e_norm

def parse_location_names(filename):
    """
    Parse the names of two locations from the given filename.

    Parameters:
        filename (str): The input filename in the format "location_1_name_vs_location_2_name_alignment.txt".

    Returns:
        tuple: A tuple containing the names of the two locations with underscores replaced by spaces.
    """
    # Remove the "_alignment.txt" suffix
    core_name = filename.replace("_alignment.csv", "")
    # Split the string into two parts using "_vs_" as the delimiter
    locations = core_name.split("_vs_")
    # Replace underscores with spaces in both locations
    location_1 = locations[0].replace("_", " ")
    location_2 = locations[1].replace("_", " ")
    return location_1, location_2
            
def process_directory(input_dir, output_filename, output_dir, name):
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    alignment_score_matrix = np.zeros((len(loc2webs), len(loc2webs)))
    locs = list(loc2webs)
    # Process each file in the input directory
    for filename in os.listdir(input_dir):
        if filename.lower().endswith('.csv'):
            input_filepath = os.path.join(input_dir, filename)
            
            if name == "muritz":
            # for muritz
                loc1, loc2 = parse_location_names(filename)
            elif name == "kai":
            # for kai's
                core_name = filename.replace(".csv", "")
                locations = core_name.split("_")
                loc1, loc2 = locations[0], locations[1]
            elif name == "wass_motif":
            # for Intermediate Wass Motif
                core_name = filename.replace("_transport_plan.csv", "")
                locations = core_name.split("_vs_")
                loc1, loc2 = locations[0], locations[1]
            else:
                print("Error: Unrecognized name")
                return

            # Handle special case for corrupted location names
            if loc1 == "Qui\uf03f\uf03fma":
                loc1 = "Quiçãma"
            if loc2 == "Qui\uf03f\uf03fma":
                loc2 = "Quiçãma"
            row_idx, col_idx = locs.index(loc1), locs.index(loc2)
            if row_idx <= col_idx:  # Only fill lower triangle
                row_idx, col_idx = col_idx, row_idx
            
            motif_filepath_1 = f"./data/motif_vectors/{loc1}_motifs.csv"
            motif_filepath_2 = f"./data/motif_vectors/{loc2}_motifs.csv"
            
            quality = calculate_zeroth_soft_alignment_quality(input_filepath, motif_filepath_1, motif_filepath_2)
            alignment_score_matrix[row_idx, col_idx] = quality
            print(f"Calculated quality score for {name} between {loc1} and {loc2}: {quality}")

    # Save only the lower triangle to the output file
    lower_triangle_matrix = np.tril(alignment_score_matrix)
    output_filepath = os.path.join(output_dir, output_filename)
    np.savetxt(output_filepath, lower_triangle_matrix, delimiter=",", fmt="%.10e")

process_directory("./data/muritz_alignment_matrices", "muritz_alignment_quality_soft.csv", "./data/alignment_scores", "muritz")
# process_directory("./data/kai's_transport_plans_gw", "kai_gw_alignment_quality_soft.csv", "./data/alignment_scores", "kai")
# process_directory("./data/motif_as_cost_matrix_transport_plans_euclidean", "euclidean_wass_motif_alignment_quality_soft.csv", "./data/alignment_scores", "wass_motif")
# process_directory("./data/motif_as_cost_matrix_transport_plans_PearsonR", "pearsonR_wass_motif_alignment_quality_soft.csv", "./data/alignment_scores", "wass_motif")