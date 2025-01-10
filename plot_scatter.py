import csv
import matplotlib.pyplot as plt

def plot_network_distances(csv1_path, csv2_path, output_path="scatter_plot.png"):
    """
    Plots a scatter plot comparing two notions of distance between networks.

    Args:
        csv1_path (str): Path to the first CSV file.
        csv2_path (str): Path to the second CSV file.
        output_path (str): Path to save the scatter plot (default is "scatter_plot.png").

    Returns:
        None
    """
    def read_csv_to_dict(csv_path):
        """
        Reads a CSV file and converts it to a dictionary of pairwise distances.

        Args:
            csv_path (str): Path to the CSV file.

        Returns:
            dict: Dictionary with pairwise distances.
        """
        distances = {}
        with open(csv_path, 'r') as csvfile:
            reader = csv.reader(csvfile)
            headers = next(reader)[1:]  # Skip the first column (row headers)
            for row in reader:
                row_name = row[0]
                for col_name, value in zip(headers, row[1:]):
                    if value:  # Skip empty cells
                        distances[(row_name, col_name)] = float(value)
        return distances

    # Read both CSV files into dictionaries
    data1 = read_csv_to_dict(csv1_path)
    data2 = read_csv_to_dict(csv2_path)

    # Ensure both CSVs have the same pairs
    common_pairs = set(data1.keys()) & set(data2.keys())
    x_values = [data1[pair] for pair in common_pairs]
    y_values = [data2[pair] for pair in common_pairs]
    # print(len(x_values))
    # Plot the scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(x_values, y_values, alpha=0.7)
    plt.xlabel("Muritz Alignment Score")
    plt.ylabel("Entropic GGW distance (undirected graph)")
    # plt.ylabel("Kai's FGW distance")
    # plt.yscale("log")
    plt.title("Comparison of Network Distances between Food Webs")
    plt.grid(True)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_path)
    print(f"Scatter plot saved to {output_path}")

# plot_network_distances("./muritz_comparison/data/muritz_pairwise.csv", "./muritz_comparison/data/kai_pairwise.csv", "./muritz_comparison/plots/muritz_kai_scatter.png")
# plot_network_distances("./muritz_comparison/data/muritz_pairwise.csv", "./muritz_comparison/data/kai_pairwise.csv", "./muritz_comparison/plots/muritz_kai_scatter_logscale.png")
# plot_network_distances("./muritz_comparison/data/muritz_pairwise.csv", "./muritz_comparison/data/ent_pairwise_distance_1e-06.csv", "./muritz_comparison/plots/muritz_GGW_scatter.png")
# plot_network_distances("./muritz_comparison/data/muritz_pairwise_sample.csv", "./muritz_comparison/data/ent_pairwise_distance_1e-06_sample31.csv", "./muritz_comparison/plots/muritz_GGW_scatter_sample.png")
# plot_network_distances("./muritz_comparison/data/muritz_pairwise_sample.csv", "./muritz_comparison/data/kai_pairwise_sample31.csv", "./muritz_comparison/plots/muritz_kai_scatter_sample.png")
