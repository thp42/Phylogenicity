import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# Define the input and output file paths
filtered_sequences_path = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/PSSM_B1.fasta"
pssm_output_path = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/pssm_PAM30_B2.csv"
pam_matrix_path = "/home/pythagoras/Documents/PhD/Evolution/PSSM Matrix/PAM30.txt"  # Update this path as needed

# Standard set of amino acids
standard_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

# Function to read sequences from a FASTA file
def read_fasta_sequences(fasta_file):
    """
    Reads protein sequences from a FASTA file and converts them to uppercase.
    """
    sequences = [str(record.seq).upper() for record in SeqIO.parse(fasta_file, "fasta")]
    return sequences

# Function to calculate the frequency matrix
def calculate_frequency_matrix(sequences):
    """
    Calculates the frequency of each standard amino acid at each position across all sequences.
    """
    if not sequences:
        raise ValueError("No sequences found in the input file.")
        
    sequence_length = len(sequences[0])
    for seq in sequences:
        if len(seq) != sequence_length:
            raise ValueError("All sequences must be of the same length.")
    
    amino_acids = sorted(set(standard_amino_acids))
    
    # Initialize the frequency matrix with zeros
    frequency_matrix = np.zeros((len(amino_acids), sequence_length))
    aa_to_index = {aa: idx for idx, aa in enumerate(amino_acids)}
    
    # Count the frequency of each amino acid at each position
    for seq in sequences:
        for i, aa in enumerate(seq):
            if aa in aa_to_index:
                frequency_matrix[aa_to_index[aa], i] += 1
    
    # Normalize the frequency matrix by the number of sequences
    frequency_matrix /= len(sequences)
    
    return frequency_matrix, amino_acids

# Function to load PAM matrix from file
def load_pam_matrix(file_path):
    """
    Loads the PAM matrix from a file and converts it to a nested dictionary.
    """
    pam = defaultdict(dict)
    with open(file_path, 'r') as f:
        lines = f.readlines()
        # Skip header lines starting with '#'
        lines = [line.strip() for line in lines if not line.startswith('#') and line.strip()]
        header = lines[0].split()
        for line in lines[1:]:
            parts = line.split()
            aa1 = parts[0]
            scores = parts[1:]
            for aa2, score in zip(header, scores):
                pam[aa1][aa2] = int(score)
                pam[aa2][aa1] = int(score)  # Assuming symmetry
    return pam

# Function to generate the PSSM matrix from the frequency matrix using PAM250 scores
def generate_pssm_matrix(frequency_matrix, amino_acids, pam_matrix):
    """
    Generates a Position-Specific Scoring Matrix (PSSM) by applying the PAM250 substitution scores
    to the frequency matrix.
    """
    pssm_matrix = np.zeros(frequency_matrix.shape)
    
    for i in range(frequency_matrix.shape[1]):  # Iterate over each position
        for j, aa in enumerate(amino_acids):    # Iterate over each amino acid
            score = 0
            for k, aa2 in enumerate(amino_acids):  # Iterate over substitution pairs
                substitution_score = pam_matrix[aa].get(aa2, 0)
                score += frequency_matrix[k, i] * substitution_score
            pssm_matrix[j, i] = score
    
    return pssm_matrix

# Read the PAM250 matrix from the file
pam_matrix = load_pam_matrix(pam_matrix_path)

# Read the filtered sequences from the FASTA file
sequences = read_fasta_sequences(filtered_sequences_path)

# Calculate the frequency matrix
frequency_matrix, amino_acids = calculate_frequency_matrix(sequences)

# Generate the PSSM matrix using PAM250
pssm_matrix = generate_pssm_matrix(frequency_matrix, amino_acids, pam_matrix)

# Convert PSSM matrix to DataFrame for easier handling and export
pssm_df = pd.DataFrame(pssm_matrix, index=list(amino_acids))

# Ensure the DataFrame includes all standard amino acids
for aa in standard_amino_acids:
    if aa not in pssm_df.index:
        pssm_df.loc[aa] = [0] * pssm_df.shape[1]

# Reorder the DataFrame to match the standard amino acids order
pssm_df = pssm_df.loc[list(standard_amino_acids)]

# Output the PSSM matrix to a CSV file
pssm_df.to_csv(pssm_output_path)

print(f"PSSM matrix saved to: {pssm_output_path}")
