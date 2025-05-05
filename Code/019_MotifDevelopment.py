from Bio import SeqIO
from ete3 import NCBITaxa
import os
from collections import defaultdict, Counter
import pandas as pd
import numpy as np

# Initialize NCBITaxa
ncbi = NCBITaxa()

# File paths
motif_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_Motif.fasta"
tree_files = {
    "ITPKA": "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA.treefile",
    "ITPKB": "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKB/ITPKB.treefile",
    "ITPKC": "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKC/ITPKC.treefile",
    "All": "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_Isoform_aligned_protein_sequences.fasta.treefile"
}
isoform_files = {
    "ITPKA": "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_ownAnnotation.fasta",
    "ITPKB": "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKB/ITPKB_ownAnnotation.fasta",
    "ITPKC": "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKC/ITPKC_ownAnnotation.fasta",
    "All": "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_Motif.fasta",
}

# List of all amino acids
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

# Function to check if a motif_id is in a tree file
def is_motif_in_treefile(treefile, motif_id):
    with open(treefile, "r") as file:
        for line in file:
            if motif_id in line:
                return True
    return False

# Function to check if a motif_id is in a FASTA file
def is_motif_in_fastafile(fastafile, motif_id):
    with open(fastafile, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            if record.id == motif_id:
                return True
    return False

# Function to get the class of a taxon given its taxid
def get_class_of_taxon(taxid):
    try:
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        
        # Define the rank order we're interested in
        rank_order = ['class', 'order']
        
        # Look for the first available rank in our preferred order
        for rank in rank_order:
            for tid in reversed(lineage):  # Start from the most specific rank
                if ranks[tid] == rank:
                    return names[tid]
        
        print(f"Warning: No suitable rank found for taxid {taxid}. Lineage: {names}")
    except Exception as e:
        print(f"Error processing taxid {taxid}: {str(e)}")
    return None

# Create dictionaries to store the counts and details of motifs found in each isoform
isoform_counts = {isoform: 0 for isoform in tree_files.keys()}
isoform_details = {isoform: [] for isoform in tree_files.keys()}
isoform_class_counts = {isoform: {} for isoform in tree_files.keys()}
class_sequences = {isoform: defaultdict(list) for isoform in tree_files.keys()}
class_matrices = {isoform: {} for isoform in tree_files.keys()}

# Read the motif file and determine which isoform each motif belongs to
with open(motif_file, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        motif_id = record.id.strip()  # Ensure there are no leading or trailing spaces
        sequence = str(record.seq)
        taxid = motif_id.split('_')[0]
        taxon_class = get_class_of_taxon(taxid)
        found = False
        for isoform, treefile in tree_files.items():
            if is_motif_in_treefile(treefile, motif_id):
                isoform_counts[isoform] += 1
                isoform_details[isoform].append((motif_id, sequence, taxid, taxon_class))
                class_sequences[isoform][taxon_class].append(sequence)
                if taxon_class:
                    if taxon_class in isoform_class_counts[isoform]:
                        isoform_class_counts[isoform][taxon_class] += 1
                    else:
                        isoform_class_counts[isoform][taxon_class] = 1
                found = True
                break
        if not found:
            for isoform, fastafile in isoform_files.items():
                if is_motif_in_fastafile(fastafile, motif_id):
                    isoform_counts[isoform] += 1
                    isoform_details[isoform].append((motif_id, sequence, taxid, taxon_class))
                    class_sequences[isoform][taxon_class].append(sequence)
                    if taxon_class:
                        if taxon_class in isoform_class_counts[isoform]:
                            isoform_class_counts[isoform][taxon_class] += 1
                        else:
                            isoform_class_counts[isoform][taxon_class] = 1
                    found = True
                    break
        if not found:
            print(f"Motif ID: {motif_id} not found in any isoform file")

# Function to generate amino acid frequency matrix
def generate_frequency_matrix(sequences):
    matrix = []
    for i in range(len(sequences[0])):
        column = [seq[i] for seq in sequences]
        total = len(column)
        freq_counter = Counter(column)
        freq_dict = {aa: freq_counter.get(aa, 0) / total for aa in amino_acids}
        matrix.append(freq_dict)
    return matrix

# Generate the frequency matrices for each class in each isoform
for isoform, class_dict in class_sequences.items():
    for class_name, sequences in class_dict.items():
        class_matrices[isoform][class_name] = generate_frequency_matrix(sequences)
        # Save each matrix as a CSV file
        df = pd.DataFrame(class_matrices[isoform][class_name])
        df = df[amino_acids]  # Ensure all amino acids are in the columns
        df.to_csv(f'{isoform}_{class_name}_matrix.csv', index=False)

# Function to compare two matrices using Pearson correlation coefficient
def compare_matrices(matrix1, matrix2):
    correlations = []
    for col1, col2 in zip(matrix1, matrix2):
        col1_values = [col1[aa] for aa in amino_acids]
        col2_values = [col2[aa] for aa in amino_acids]
        if np.sum(col1_values) == 0 or np.sum(col2_values) == 0:
            correlations.append(0)  # Handle cases where columns have no values
        else:
            correlation = np.corrcoef(col1_values, col2_values)[0, 1]
            correlations.append(correlation)
    return np.mean(correlations)

# Compare each matrix against each other matrix (within and across isoforms)
comparison_results = []
for isoform1, class_dict1 in class_matrices.items():
    for class1, matrix1 in class_dict1.items():
        # Compare within the same isoform
        for class2, matrix2 in class_dict1.items():
            if class1 != class2:
                similarity = compare_matrices(matrix1, matrix2)
                comparison_results.append((isoform1, class1, len(class_sequences[isoform1][class1]), 
                                           isoform1, class2, len(class_sequences[isoform1][class2]), similarity))
        # Compare across different isoforms for identical classes
        for isoform2, class_dict2 in class_matrices.items():
            if isoform1 != isoform2 and class1 in class_dict2:
                matrix2 = class_dict2[class1]
                similarity = compare_matrices(matrix1, matrix2)
                comparison_results.append((isoform1, class1, len(class_sequences[isoform1][class1]), 
                                           isoform2, class1, len(class_sequences[isoform2][class1]), similarity))

# Save comparison results to a text file
with open('matrix_comparisons.txt', 'w') as file:
    for result in comparison_results:
        isoform1, class1, count1, isoform2, class2, count2, similarity = result
        file.write(f"Similarity between {class1} ({count1} sequences) in {isoform1} and {class2} ({count2} sequences) in {isoform2}: {similarity:.4f}\n")

# Additionally, write the number of sequences for each class in each isoform
with open('class_sequence_counts.txt', 'w') as file:
    for isoform, class_dict in class_sequences.items():
        file.write(f"\nIsoform: {isoform}\n")
        for class_name, sequences in class_dict.items():
            file.write(f"  {class_name}: {len(sequences)} sequences\n")
