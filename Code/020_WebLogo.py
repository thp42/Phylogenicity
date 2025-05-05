import pandas as pd
import numpy as np
import random
import subprocess
import glob
import os

# Directory containing the CSV files
directory_path = '/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/Motif/Development'

# Get all CSV files in the directory
csv_files = glob.glob(os.path.join(directory_path, '*.csv'))

# Amino acids (assuming these are the columns in the original PSSM matrix)
# We will get these from the first file and assume they are the same for all files
if csv_files:
    first_file = csv_files[0]
    pssm_matrix = pd.read_csv(first_file)
    amino_acids = list(pssm_matrix.columns)
else:
    amino_acids = []

# Number of sequences to generate
num_sequences = 100

# Function to generate a single sequence based on probabilities in the PSSM
def generate_sequence(pssm):
    sequence = ""
    for i in range(pssm.shape[1]):
        probs = pssm[:, i]
        sequence += random.choices(amino_acids, probs)[0]
    return sequence

# Process each CSV file
for file_path in csv_files:
    # Load the PSSM matrix from the CSV file
    pssm_matrix = pd.read_csv(file_path)

    # Transpose the matrix to match the expected input format for sequence generation
    transposed_matrix = pssm_matrix.transpose()

    # Generate sequences
    sequences = [generate_sequence(transposed_matrix.to_numpy()) for _ in range(num_sequences)]

    # Build the output file paths based on the input file name
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    fasta_file_path = os.path.join(directory_path, f'{base_name}_generated_sequences.fasta')
    output_file_path = os.path.join(directory_path, f'{base_name}_sequence_logo.png')

    # Save sequences to a FASTA file
    with open(fasta_file_path, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f'>seq{i+1}\n')
            f.write(f'{seq}\n')

    # Generate the sequence logo using WebLogo
    result = subprocess.run([
        'weblogo', '-f', fasta_file_path, '-o', output_file_path, 
        '--format', 'svg', '--sequence-type', 'protein',
        '--errorbars', 'no',  # Remove error bars
        '--color-scheme', 'chemistry'  # Apply the chemistry color scheme
    ], capture_output=True, text=True)

    if result.returncode == 0:
        print(f"Sequence logo for {base_name} has been saved as {output_file_path}")
    else:
        print(f"Failed to generate sequence logo for {base_name}:")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
