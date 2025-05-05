# Import necessary packages
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import multiprocessing
import re

# Function to modify headers
def modify_headers(seqs):
    new_seqs = []
    for seq in seqs:
        # Replace spaces and special characters with underscores
        new_name = re.sub(r'[^A-Za-z0-9]', '_', seq.id)
        new_seq = SeqRecord(seq.seq, id=new_name, description="")
        new_seqs.append(new_seq)
    return new_seqs

# Load your amino acid sequences
fasta_file = "//home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/PSSM_02.fasta"
seqs = list(SeqIO.parse(fasta_file, "fasta"))

# Modify headers
seqs = modify_headers(seqs)

# Save the aligned sequences to the same directory
output_file = os.path.join(os.path.dirname(fasta_file), "ITPK_Motif.fasta")
SeqIO.write(seqs, output_file, "fasta")

# Print the file path for confirmation
print(f"Modified sequences saved to: {output_file}")

# Register a parallel backend (if needed)
num_workers = multiprocessing.cpu_count()
print(f"Number of CPU cores available: {num_workers}")

# Assuming further parallel processing tasks, here's an example using multiprocessing
def process_sequence(seq):
    # Dummy function for example purposes
    return seq

if __name__ == "__main__":
    with multiprocessing.Pool(num_workers) as pool:
        processed_seqs = pool.map(process_sequence, seqs)
    
    # Do something with processed_seqs if needed
    print("Parallel processing complete.")
