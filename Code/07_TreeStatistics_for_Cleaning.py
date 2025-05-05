from Bio import Phylo
import numpy as np
import matplotlib.pyplot as plt

def calculate_branch_lengths(tree_file):
    # Load the tree
    tree = Phylo.read(tree_file, 'newick')
    
    # Extract all branch lengths
    branch_lengths = []
    for clade in tree.find_clades(order='preorder'):
        if clade.branch_length:
            branch_lengths.append(clade.branch_length)
    
    return branch_lengths

def calculate_bootstrap_values(tree_file):
    # Load the tree
    tree = Phylo.read(tree_file, 'newick')
    
    # Extract all bootstrap values
    bootstrap_values = []
    for clade in tree.find_clades(order='preorder'):
        if clade.confidence is not None:
            bootstrap_values.append(clade.confidence)
    
    return bootstrap_values

def suggest_cutoff_std_dev(values, factor):
    # Calculate the cutoff using the standard deviation method
    mean_value = np.mean(values)
    std_dev = np.std(values)
    cutoff_std_dev = mean_value + factor * std_dev
    print(f"Suggested cutoff (mean + {factor} * std_dev): {cutoff_std_dev}")
    
    return cutoff_std_dev

def analyze_values(values, cutoff_std_dev, title):
    # Calculate statistics
    mean_value = np.mean(values)
    median_value = np.median(values)
    std_dev = np.std(values)
    max_value = np.max(values)
    min_value = np.min(values)
    
    # Print statistics
    print(f"Mean {title}: {mean_value}")
    print(f"Median {title}: {median_value}")
    print(f"Standard Deviation: {std_dev}")
    print(f"Max {title}: {max_value}")
    print(f"Min {title}: {min_value}")
    
    # Plot a histogram of the values with cutoff line
    plt.hist(values, bins=30, edgecolor='black')
    plt.axvline(cutoff_std_dev, color='red', linestyle='dashed', linewidth=2, label=f'Std Dev Cutoff ({cutoff_std_dev:.2f})')
    plt.title(f'Distribution of {title}')
    plt.xlabel(title)
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()

def main():
    # Modify the following line to adjust the parameters
    tree_file = '/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_aligned_protein_sequences_with_isoforms.fasta.treefile'
    factor = 7  # Adjust this for std_dev method
    
    # Analyze branch lengths
    branch_lengths = calculate_branch_lengths(tree_file)
    cutoff_std_dev_bl = suggest_cutoff_std_dev(branch_lengths, factor)
    analyze_values(branch_lengths, cutoff_std_dev_bl, "Branch Lengths")
    
    # Analyze bootstrap values (report statistics, no pruning)
    bootstrap_values = calculate_bootstrap_values(tree_file)
    analyze_values(bootstrap_values, min(bootstrap_values), "Bootstrap Values")

if __name__ == "__main__":
    main()
