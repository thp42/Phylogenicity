from Bio import SeqIO

# File paths
fasta_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/PSSM_02.fasta"
annotation_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKC/ITPKC_ownAnnotation.fasta"
output_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKC/Filtered_Flank_All.fasta"

# Step 1: Read and adjust identifiers from ITPKA_Annotation.fasta
identifiers = []
with open(annotation_file, 'r') as ann_file:
    for line in ann_file:
        if line.startswith('>'):
            # Example: >7868_0_00219e___ becomes 7868_0:00219e
            parts = line[1:].strip().split('_')
            adjusted_id = f"{parts[0]}_{parts[1]}:{parts[2]}"
            identifiers.append(adjusted_id)

# Step 2: Loop through identifiers and match with sequences in the fasta_file
filtered_sequences = []
count_matched = 0
count_not_matched = 0

# Create a dictionary from the fasta file for fast lookup
fasta_dict = {record.id.split(' ')[0]: record for record in SeqIO.parse(fasta_file, "fasta")}

for identifier in identifiers:
    if identifier in fasta_dict:
        filtered_sequences.append(fasta_dict[identifier])
        count_matched += 1
    else:
        count_not_matched += 1

# Step 3: Write filtered sequences to output file
with open(output_file, 'w') as out_file:
    SeqIO.write(filtered_sequences, out_file, "fasta")

# Print the number of matching and non-matching sequences
print(f"Number of matching sequences found: {count_matched}")
print(f"Number of non-matching sequences: {count_not_matched}")
print(f"Filtered sequences saved to {output_file}")
