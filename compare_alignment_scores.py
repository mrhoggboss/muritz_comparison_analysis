import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pickle

webs_pkl_path = "./data/loc2webs.pkl"

with open(webs_pkl_path, "rb") as f: 
    loc2webs = pickle.load(f)

locs = list(loc2webs)
NAME1 = "kai"
NAME2 = "muritz"
indices_to_exclude = [115, 89] # corresponds to Skeleton Coast and Namaqua
# indices_to_exclude = []
locations_to_exclude = ['Skeleton Coast', 'Namaqua']
# locations_to_exclude = []
def plot_comparison_matrix(mode):
    if mode == 'gw':
        # Load matrices with row and column headers
        matrix1 = pd.read_csv(f"./data/alignment_scores/{NAME1}_alignment_quality_{mode}.csv", index_col=0).drop(index=locations_to_exclude, columns=locations_to_exclude)
        matrix2 = pd.read_csv(f"./data/alignment_scores/{NAME2}_alignment_quality_{mode}.csv", header=None).drop(index=locations_to_exclude, columns=locations_to_exclude)

        # Assign row and column labels from matrix1 to matrix2
        matrix2.index = matrix1.index  # Row labels
        matrix2.columns = matrix1.columns  # Column labels

        # Force matrix1 and matrix2 to be lower triangular
        matrix1_lower = pd.DataFrame(np.tril(matrix1.values), index=matrix1.index, columns=matrix1.columns)
        matrix2_lower = pd.DataFrame(np.tril(matrix2.values), index=matrix2.index, columns=matrix2.columns)

        # Calculate difference
        difference = matrix1_lower - matrix2_lower

        # Plot heatmaps
        fig, axes = plt.subplots(1, 3, figsize = (24, 7))
        sns.heatmap(matrix1_lower, ax=axes[0], cmap="viridis", annot=False)
        sns.heatmap(matrix2_lower, ax=axes[1], cmap="viridis", annot=False)
        sns.heatmap(difference, ax=axes[2], cmap="coolwarm", center=0, annot=False)

        # Add titles
        axes[0].set_title("Kai")
        axes[1].set_title("Muritz")
        axes[2].set_title(f"{mode} Transport Cost Difference (Kai - Muritz)")

        # Add row and column labels
        for ax in axes:
            ax.axis('off')
            # ax.set_xticklabels(locs, rotation=45, ha="right")
            # ax.set_yticklabels(locs, rotation=0)

        plt.tight_layout()
        plt.show()

    else:
        # Load matrices
        matrix1 = pd.read_csv(f"./data/alignment_scores/{NAME1}_alignment_quality_{mode}.csv", header=None).drop(index=indices_to_exclude, axis=0).drop(columns=indices_to_exclude, axis=1)
        matrix2 = pd.read_csv(f"./data/alignment_scores/{NAME2}_alignment_quality_{mode}.csv", header=None).drop(index=indices_to_exclude, axis=0).drop(columns=indices_to_exclude, axis=1)
        # Assign row and column labels from matrix1 to matrix2
        matrix2.index = matrix1.index  # Row labels
        matrix2.columns = matrix1.columns  # Column labels

        # Force matrix1 and matrix2 to be lower triangular
        matrix1_lower = pd.DataFrame(np.tril(matrix1.values), index=matrix1.index, columns=matrix1.columns)
        matrix2_lower = pd.DataFrame(np.tril(matrix2.values), index=matrix2.index, columns=matrix2.columns)
        
        difference = matrix1 - matrix2
        
        # Plot heatmaps
        fig, axes = plt.subplots(1, 3, figsize = (24, 7))
        sns.heatmap(matrix1, ax=axes[0], cmap="viridis", annot=False)
        sns.heatmap(matrix2, ax=axes[1], cmap="viridis", annot=False)
        sns.heatmap(difference, ax=axes[2], cmap="coolwarm", center=0, annot=False)

        # Add titles
        axes[0].set_title(NAME1)
        axes[1].set_title(NAME2)
        axes[2].set_title(f"{mode} Alignment Quality Difference ({NAME1} - {NAME2})")

        # Add row and column labels
        for ax in axes:
            ax.axis('off')
            # ax.set_xticklabels(locs, rotation=45, ha="right")
            # ax.set_yticklabels(locs, rotation=0)
        print(np.sum(difference.to_numpy()))
        print(np.sum(matrix1.to_numpy()))
        plt.tight_layout()
        plt.show()

# plot_comparison_matrix("hard")
plot_comparison_matrix("soft")
# plot_comparison_matrix("gw")

print(locs[115])
print(locs.index('Namaqua'))
# print(locs[54])
# print(locs[96])
# print(loc2webs[locs[115]])

