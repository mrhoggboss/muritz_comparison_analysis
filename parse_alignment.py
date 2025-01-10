import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd
import pickle
import os 

webs_pkl_path = "./data/loc2webs.pkl"

with open(webs_pkl_path, "rb") as f: 
    loc2webs = pickle.load(f)

# Step 1: Read the alignment data
def parse_alignment(file_path, filename):
    with open(file_path, 'r') as f:
        data = f.read().strip()
    # Extract node pairs using regex
    pairs = re.findall(r'\(([^,]+),([^)]+)\)', data)
    print(file_path)
    loc1, loc2 = parse_location_names(filename)
    return pairs, loc1, loc2

# Step 1.5: Sort the species names according to what we want
def order_species(csv_path, loc):
    """
    Orders the provided list of species based on the following criteria:
    1. Partition into species that eat mammals and species that don't (based on 'Mammal' column).
    2. Sort each partition in descending order of body mass ('Mass.g').
    3. Combine the two partitions: [eats mammals, sorted by mass] + [doesn't eat mammals, sorted by mass].

    Args:
        csv_path (str): Path to the CSV file containing 'Species', 'Mass.g', and 'Mammal' columns.
        species_list (list): List of species names to be sorted.

    Returns:
        list: A sorted list of species names.
    """

    # Step 1: Read the CSV file
    df = pd.read_csv(csv_path)

    # Step 2: Filter the DataFrame to only include the species in species_list
    df_filtered = df[df['Species'].isin(species_list)]

    # # Step 3: Partition species into two groups based on 'Mammal' column
    # eats_mammals = df_filtered[df_filtered['Mammal'] != 0]  # Species that eat mammals
    # doesnt_eat_mammals = df_filtered[df_filtered['Mammal'] == 0]  # Species that don't eat mammals

    # Step 3: Partition species into two groups base on whether they eat another mammal in this web
    eats_mammals = df_filtered[
        df_filtered['Species'].apply(lambda species: loc2webs[loc].out_degree(species) != 0)
    ]

    doesnt_eat_mammals = df_filtered[
        df_filtered['Species'].apply(lambda species: loc2webs[loc].out_degree(species) == 0)
    ]
    # Step 4: Sort each group by descending body mass ('Mass.g')
    eats_mammals_sorted = eats_mammals.sort_values(by='Mass.g', ascending=False)
    doesnt_eat_mammals_sorted = doesnt_eat_mammals.sort_values(by='Mass.g', ascending=False)

    # Step 5: Combine the two sorted partitions
    ordered_species = list(eats_mammals_sorted['Species']) + list(doesnt_eat_mammals_sorted['Species'])

    return ordered_species

# Step 2: Create heatmap matrix
def create_heatmap_matrix(pairs):
    nodes_1 = order_species('./data/functional_traits.csv', 'Luengue Luiana')
    nodes_2 = order_species('./data/functional_traits.csv', 'West Lunga')
    print(nodes_1)
    print(nodes_2)
    # Map nodes to indices
    node_1_to_idx = {node: idx for idx, node in enumerate(nodes_1)}
    node_2_to_idx = {node: idx for idx, node in enumerate(nodes_2)}
    
    # Create an empty matrix
    heatmap = np.zeros((len(nodes_1), len(nodes_2)))

    shared_species = set(nodes_1).intersection(set(nodes_2))

    # Fill the matrix based on alignments and record colors
    species_colors_1 = dict()
    species_colors_2 = dict()
    
    for node_1 in nodes_1:
        if node_1 in shared_species:
            species_colors_1[node_1] = 'blue'
        else:
            species_colors_1[node_1] = 'black'
    for node_2 in nodes_2:
        if node_2 in shared_species:
            species_colors_2[node_2] = 'blue'
        else:
            species_colors_2[node_2] = 'black'

    for node_1, node_2 in pairs:
        if node_1 != 'NULL' and node_2 != 'NULL':
            if node_1 == node_2:
                heatmap[node_1_to_idx[node_1], node_2_to_idx[node_2]] = 2
                species_colors_1[node_1] = 'green'
                species_colors_2[node_2] = 'green'
            else:
                heatmap[node_1_to_idx[node_1], node_2_to_idx[node_2]] = 1 

    return heatmap, nodes_1, nodes_2, species_colors_1, species_colors_2

def parse_location_names(filename):
    """
    Parse the names of two locations from the given filename.

    Parameters:
        filename (str): The input filename in the format "location_1_name_vs_location_2_name_alignment.txt".

    Returns:
        tuple: A tuple containing the names of the two locations with underscores replaced by spaces.
    """
    # Remove the "_alignment.txt" suffix
    core_name = filename.replace("_alignment.txt", "")
    # Split the string into two parts using "_vs_" as the delimiter
    locations = core_name.split("_vs_")
    # Replace underscores with spaces in both locations
    location_1 = locations[0].replace("_", " ")
    location_2 = locations[1].replace("_", " ")
    return location_1, location_2

def create_alignment_matrix(pairs, loc1, loc2):
    if loc1 == "Qui\uf03f\uf03fma":
        loc1 = "Quiçãma"
    if loc2 == "Qui\uf03f\uf03fma":
        loc2 = "Quiçãma"
    # nodes_1 = order_species('./data/functional_traits.csv', loc1)
    # nodes_2 = order_species('./data/functional_traits.csv', loc2)
    nodes_1 = list(loc2webs[loc1].nodes())
    nodes_2 = list(loc2webs[loc2].nodes())
    # print(nodes_1)
    # print(nodes_2)
    # Map nodes to indices
    node_1_to_idx = {node: idx for idx, node in enumerate(nodes_1)}
    node_2_to_idx = {node: idx for idx, node in enumerate(nodes_2)}
    
    # Create an empty matrix
    heatmap = np.zeros((len(nodes_1), len(nodes_2)))

    for node_1, node_2 in pairs:
        if node_1 != 'NULL' and node_2 != 'NULL':
            heatmap[node_1_to_idx[node_1], node_2_to_idx[node_2]] = 1 

    return heatmap, nodes_1, nodes_2

# Step 3: Update dataset with common names
def species_to_common(csv_path):
    # Read the CSV file
    name_mapping = pd.read_csv(csv_path)
    species_to_common = dict(zip(name_mapping['Species'], name_mapping['Common_Name']))
    ret_dict = dict()

    for node_1, node_2 in pairs:
        if node_1 != 'NULL':
            ret_dict[node_1] = species_to_common.get(node_1)
        if node_2 != 'NULL':
            ret_dict[node_2] = species_to_common.get(node_2)
    
    return ret_dict

# Step 4: Plot heatmap
def plot_heatmap(heatmap, nodes_1, nodes_2, name1, name2, species_to_common, species_colors_1, species_colors_2):
    """
    Plots a heatmap with custom coloring for alignment and optional custom coloring for species names.

    Args:
        heatmap (np.ndarray): 2D array representing the heatmap values.
        nodes_1 (list): List of species on the y-axis.
        nodes_2 (list): List of species on the x-axis.
        name1 (str): Name of the first graph/species set.
        name2 (str): Name of the second graph/species set.
        species_colors (dict): Optional mapping of species names to colors (e.g., {'species_name': 'red'}).
    """
    from matplotlib.colors import ListedColormap

    # Define a custom colormap for the heatmap
    cmap = ListedColormap(['white', 'blue', 'blue'])

    # Normalize data to match colormap
    normalized_heatmap = np.zeros_like(heatmap)
    normalized_heatmap[heatmap == 1] = 1  
    normalized_heatmap[heatmap == 2] = 2 

    plt.figure(figsize=(10, 8))
    plt.imshow(normalized_heatmap, cmap=cmap, interpolation='nearest', aspect='auto')

    # Set axis labels and titles
    plt.xlabel(f"{name2} species")
    plt.ylabel(f"{name1} species")
    plt.title("Muritz Node Alignment Heatmap between Luengue Luiana and West Lunga")

    # Set axis tick labels
    ylabels = [species_to_common[node_2] for node_2 in nodes_2]
    xlabels = [species_to_common[node_1] for node_1 in nodes_1]
    plt.xticks(ticks=np.arange(len(nodes_2)), labels=ylabels, rotation=90)
    plt.yticks(ticks=np.arange(len(nodes_1)), labels=xlabels)
    # print(xlabels)
    # print(ylabels)

    # print(nodes_1)
    # print(nodes_2)

    ax = plt.gca()
    # Color x-axis tick labels
    for tick_label, node in zip(ax.get_xticklabels(), nodes_2):
        tick_label.set_color(species_colors_2.get(node, 'black'))  # Default to black if no color specified
    # Color y-axis tick labels
    for tick_label, node in zip(ax.get_yticklabels(), nodes_1):
        tick_label.set_color(species_colors_1.get(node, 'black'))  # Default to black if no color specified
    plt.tight_layout()
    plt.show()

import csv
import numpy as np

def save_alignment_matrix(matrix, row_species, column_species, output_filepath):
    """
    Save an alignment matrix into a CSV file with row and column headers.

    Parameters:
        matrix (numpy.ndarray): 2D NumPy array representing the alignment matrix.
        row_species (list): List of species names corresponding to the rows.
        column_species (list): List of species names corresponding to the columns.
        output_filepath (str): Path to save the output CSV file.
    """
    # Ensure the matrix dimensions match the species lists
    if len(row_species) != matrix.shape[0] or len(column_species) != matrix.shape[1]:
        raise ValueError("Matrix dimensions do not match the provided species lists.")

    # Write the alignment matrix to a CSV file
    with open(output_filepath, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write the header row
        csvwriter.writerow([""] + column_species)
        # Write each row with its corresponding species name
        for row_name, row_data in zip(row_species, matrix):
            csvwriter.writerow([row_name] + list(row_data))

def process_directory(input_dir, output_dir):
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Process each file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith('alignment.txt'):
            input_filepath = os.path.join(input_dir, filename)
            output_filename = os.path.splitext(filename)[0] + '.csv'
            output_filepath = os.path.join(output_dir, output_filename)
            pairs, loc1, loc2 = parse_alignment(input_filepath, filename)
            heatmap, nodes_1, nodes_2 = create_alignment_matrix(pairs, loc1, loc2)
            save_alignment_matrix(heatmap,nodes_1,nodes_2, output_filepath)
            print (f"Processed {input_filepath} to {output_filepath}")

# Example usage
# file_path = './alignment/Luengue_Luiana_vs_West_Lunga_alignment.txt'
# csv_path = './data/Common_Names_Mammals.csv'

# pairs = parse_alignment(file_path)
# heatmap, nodes_1, nodes_2 = create_alignment_matrix(pairs)

# save_alignment_matrix(heatmap, nodes_1, nodes_2, "./data/testing.csv")

process_directory("./data/muritz_results_raw", "./data/muritz_alignment_matrices")