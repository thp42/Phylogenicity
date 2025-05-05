from Bio import SeqIO

def filter_common_sequences(fl_file, flank_file, output_fl_file, output_flank_file):
    # Read the headers from the FL FASTA file
    fl_records = {record.id: record for record in SeqIO.parse(fl_file, "fasta")}
    flank_records = {record.id: record for record in SeqIO.parse(flank_file, "fasta")}
    
    # Find the common headers
    common_headers = fl_records.keys() & flank_records.keys()
    
    # Print the number of sequences before filtering
    print(f"Number of sequences in FL file before filtering: {len(fl_records)}")
    print(f"Number of sequences in Flank file before filtering: {len(flank_records)}")
    print(f"Number of sequences in common: {len(common_headers)}")
    
    # Filter and write the common sequences to new files
    with open(output_fl_file, 'w') as output_fl_handle, open(output_flank_file, 'w') as output_flank_handle:
        for header in common_headers:
            SeqIO.write(fl_records[header], output_fl_handle, "fasta-2line")
            SeqIO.write(flank_records[header], output_flank_handle, "fasta-2line")
    
    # Print the number of sequences after filtering
    print(f"\nNumber of sequences in both files after filtering: {len(common_headers)}")

# Example usage
fl_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_FL_DNA_Sequences.filtered.fasta"
flank_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_Flank_DNA_Sequences.filtered.fasta"
output_fl_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_FL_DNA_Sequences_Common.fasta"
output_flank_file = "//home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_Flank_DNA_Sequences_Common.fasta"

filter_common_sequences(fl_file, flank_file, output_fl_file, output_flank_file)
