import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import numpy as np
# Load the data
METRIC = 'PearsonR'
file_path = f'./data/alignment_scores/motif_as_cost_matrix_transport_cost_{METRIC}.csv'
data = pd.read_csv(file_path, encoding='latin1')

# Extract labels and data for the heatmap
row_labels = data['Unnamed: 0']
heatmap_data = data.drop(columns='Unnamed: 0').set_index(row_labels)

# Create a mask for the upper triangular part, including the diagonal
mask = np.triu(np.ones_like(heatmap_data, dtype=bool))

# Plot the heatmap with the mask applied
plt.figure(figsize=(18, 15))
sns.heatmap(heatmap_data, mask=mask, cmap="viridis", xticklabels=True, yticklabels=True) # Normalize the data using Min-Max Scaling
plt.title(f"Intermediate Wasserstein Motif Transport Cost - {METRIC}")
plt.xlabel("Locations")
plt.ylabel("Locations")
plt.show()

# # Shift values to make them positive and magnify differences
# shifted_positive_data = (heatmap_data - heatmap_data.min().min()) * 1

# # Plot the heatmap with the transformed data
# plt.figure(figsize=(15, 12))
# sns.heatmap(shifted_positive_data, cmap="viridis", xticklabels=True, yticklabels=True)
# plt.title("Shifted and Magnified Heatmap of Pearson Correlation Matrix")
# plt.xlabel("Labels")
# plt.ylabel("Labels")
# plt.show()
