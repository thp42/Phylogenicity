from Bio import SeqIO
import re

def filter_sequences(input_fasta, output_fasta):
    # Define a pattern to match the descriptions of interest (case-sensitive)
    pattern = re.compile(r"(inositol-trisphosphate 3-kinase C|itpkc|ITPKC)")
    
    # Read the sequences from the input FASTA file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Filter sequences
    filtered_sequences = []
    for seq_record in sequences:
        header = seq_record.description
        if pattern.search(header):
            filtered_sequences.append(seq_record)
    
    # Write the filtered sequences to the output FASTA file
    SeqIO.write(filtered_sequences, output_fasta, "fasta")

# Define the input and output FASTA file paths
input_fasta = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_OrthoDB.fasta"
output_fasta = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKC/ITPKC_Isoform.fasta"

# Run the filter function
filter_sequences(input_fasta, output_fasta)
