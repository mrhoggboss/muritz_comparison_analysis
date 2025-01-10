import csv
import numpy as np
from scipy.stats import pearsonr
import os
import pickle

webs_pkl_path = "./data/loc2webs.pkl"

with open(webs_pkl_path, "rb") as f: 
    loc2webs = pickle.load(f)

def load_csv(filepath):
    with open(filepath, 'r') as file:
        reader = csv.reader(file)
        data = list(reader)
    return data

def parse_alignment_matrix(filepath):
    data = load_csv(filepath)
    species_row = data[0][1:]  # Column headers (species from network 2)
    species_col = [row[0] for row in data[1:]]  # Row headers (species from network 1)
    matrix = np.array([[float(value) for value in row[1:]] for row in data[1:]])
    return matrix, species_row, species_col

def parse_motif_vectors(filepath):
    data = load_csv(filepath)
    species = [row[0] for row in data[1:]]
    vectors = {row[0]: np.array([float(value) for value in row[1:]]) for row in data[1:]}
    return species, vectors

import numpy as np
from scipy.stats import pearsonr

def calculate_soft_alignment_quality(alignment_filepath, motif_filepath_1, motif_filepath_2):
    """
    Calculate alignment quality for soft alignments.
        
    Returns:
        float: Generalized alignment quality.
    """
    weighted_cost = 0
    total_weight = 0

    alignment_matrix, species_row, species_col = parse_alignment_matrix(alignment_filepath)
    _, motif_vectors_1 = parse_motif_vectors(motif_filepath_1)
    _, motif_vectors_2 = parse_motif_vectors(motif_filepath_2)

    # Iterate through all pairs in the alignment matrix
    for i, species_i in enumerate(motif_vectors_1.keys()):  # Rows (network 2)
        for j, species_j in enumerate(motif_vectors_2.keys()):  # Columns (network 1)
            weight = alignment_matrix[i, j]
            if weight > 0:  # Consider non-zero weights
                # Retrieve motif vectors
                vector_i = motif_vectors_1[species_i]
                vector_j = motif_vectors_2[species_j]
                # Compute Pearson correlation
                rho, _ = pearsonr(vector_i, vector_j)
                # Accumulate weighted cost
                weighted_cost += weight * (1 - rho)
                total_weight += weight

    # Normalize the cost
    if total_weight == 0:  # Avoid division by zero
        return 0
    e_norm = weighted_cost / total_weight
    return e_norm

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

# converting soft alignments to hard alignments
def naive_convert(soft_alignment):
    hard_alignment = np.zeros(soft_alignment.shape)
    for row_idx in range(len(hard_alignment)):
        strongest_link = max(soft_alignment[row_idx])
        for col_idx in range(len(hard_alignment[row_idx])):
            if (soft_alignment[row_idx, col_idx] == strongest_link):
                hard_alignment[row_idx, col_idx] = 1
            else:
                hard_alignment[row_idx, col_idx] = 0
    return hard_alignment
            
def process_directory(input_dir, output_filename, output_dir):
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    alignment_score_matrix = np.zeros((len(loc2webs), len(loc2webs)))
    locs = list(loc2webs)
    counter = 0
    pc= 0
    sc = 0
    # Process each file in the input directory
    for filename in os.listdir(input_dir):
        sc += 1
        if filename.lower().endswith('.csv'):
            input_filepath = os.path.join(input_dir, filename)
            
            # for muritz
            # loc1, loc2 = parse_location_names(filename)

            # for kai's
            # core_name = filename.replace(".csv", "")
            # locations = core_name.split("_")
            # loc1, loc2 = locations[0], locations[1]

            # for Intermediate Wass Motif
            core_name = filename.replace("_transport_plan.csv", "")
            locations = core_name.split("_vs_")
            loc1, loc2 = locations[0], locations[1]

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
            
            quality = calculate_soft_alignment_quality(input_filepath, motif_filepath_1, motif_filepath_2)
            alignment_score_matrix[row_idx, col_idx] = quality
            if quality == 0:
                pc += 1
            else:
                counter += 1
            print(f"Calculated quality score between {loc1} and {loc2}: {quality}")

            
    print(f"number of nonzero dissimilarity: {counter}")
    print(f"number of 0 dissimilarity: {pc}")
    print(sc)

    # Save only the lower triangle to the output file
    lower_triangle_matrix = np.tril(alignment_score_matrix)
    output_filepath = os.path.join(output_dir, output_filename)
    np.savetxt(output_filepath, lower_triangle_matrix, delimiter=",", fmt="%.10e")

# process_directory("./data/muritz_alignment_matrices", "muritz_alignment_quality_soft.csv", "./data/alignment_scores")
# process_directory("./data/kai's_transport_plans", "kai_alignment_quality_soft.csv", "./data/alignment_scores")
# process_directory("./data/motif_as_cost_matrix_transport_plans_euclidean", "motif_as_cost_matrix_euclidean_soft.csv", "./data/alignment_scores")
process_directory("./data/motif_as_cost_matrix_transport_plans_PearsonR", "motif_as_cost_matrix_PearsonR_soft.csv", "./data/alignment_scores")