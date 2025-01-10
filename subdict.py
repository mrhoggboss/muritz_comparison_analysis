import csv

import pandas as pd

# Load the pairwise distance matrix from the CSV file
input_file = './data/kai_pairwise.csv'
matrix = pd.read_csv(input_file, index_col=0)

# List of keys (row and column labels) to extract
keys_to_extract = ['Amboseli', 'Bangweulu', 'Boumba Bek', 'Gile', 'Gombe', 'Gonarezhou', 'Kafue', 'Kafue Flats', 'Kibale', 'Kidepo Valley', 'Kilimanjaro', 'Kizigo', 'Kourtiagou', 'Kruger', 'Lake Nakuru', 'Luengue Luiana', 'Mago', 'Mahango', 'Manovo Gounda Saint Floris', 'Masai Mara', 'Ngotto Forest', 'Nkhotakota', 'Nouabale Ndoki', 'Nyungwe', 'Samburu', 'Simien Mountains', 'Siniaka Minia', 'Tsavo', 'Volcans', 'Zakouma', 'Zinave']


# Extract the submatrix
submatrix = matrix.loc[keys_to_extract, keys_to_extract]

# Save the submatrix to a new CSV file
output_file = './data/kai_pairwise_sample31.csv'
submatrix.to_csv(output_file)

print(f"Submatrix saved to '{output_file}'")
