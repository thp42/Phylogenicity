from Bio import SeqIO
from Bio.Seq import Seq

def filter_and_shorten_sequences(input_fasta, output_fasta):
    ambiguous_bases = set('NRYWSKMBDHV')
    stop_codons = ['TAA', 'TAG', 'TGA']
    total_sequences = 0
    ambiguous_filtered = 0
    length_filtered = 0
    stop_codons_removed = 0
    duplicate_filtered = 0
    kept_sequences = 0

    seen_ids = set()

    def is_ambiguous(sequence):
        """Check if the sequence contains any ambiguous nucleotides."""
        return any(base in ambiguous_bases for base in sequence)

    def remove_stop_codon_from_end(sequence):
        """Remove a stop codon only if it is at the end of the sequence."""
        seq_str = str(sequence)
        for codon in stop_codons:
            if seq_str.endswith(codon):
                seq_str = seq_str[:-3]  # Remove the last 3 nucleotides (the stop codon)
                break
        return Seq(seq_str)

    def shorten_id(seq_id):
        """Shorten the sequence ID to 10 characters."""
        return seq_id.split()[0][:10]

    with open(output_fasta, 'w') as output_handle:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            total_sequences += 1

            # Shorten the sequence ID
            record.id = shorten_id(record.id)
            record.description = ""  # Remove description to avoid conflicts

            # Check if the shortened sequence ID is already seen (duplicate)
            if record.id in seen_ids:
                duplicate_filtered += 1
                continue

            # Mark the sequence ID as seen
            seen_ids.add(record.id)

            # Filter by ambiguous bases
            if is_ambiguous(record.seq):
                ambiguous_filtered += 1
                continue

            # Remove stop codon only if it's at the end of the sequence
            original_seq = record.seq
            record.seq = remove_stop_codon_from_end(record.seq)

            # Check if stop codons were removed
            if len(record.seq) < len(original_seq):
                stop_codons_removed += 1

            # Filter by sequence length not divisible by 3
            if len(record.seq) % 3 != 0:
                length_filtered += 1
                continue

            # Write the sequence to the output if it passes all filters
            SeqIO.write(record, output_handle, 'fasta')
            kept_sequences += 1

    print(f"Total sequences processed: {total_sequences}")
    print(f"Sequences filtered due to ambiguous bases: {ambiguous_filtered}")
    print(f"Sequences filtered due to length not divisible by 3: {length_filtered}")
    print(f"Sequences with stop codons removed from the end: {stop_codons_removed}")
    print(f"Sequences filtered due to duplicate names: {duplicate_filtered}")
    print(f"Sequences retained: {kept_sequences}")
    print(f"Filtered sequences written to {output_fasta}")

if __name__ == "__main__":
    input_fasta = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_Flank_DNA_Sequences.fasta"  # Replace with your input file path
    output_fasta = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_Flank_DNA_Sequences.filtered.fasta"  # Replace with your desired output file path

    filter_and_shorten_sequences(input_fasta, output_fasta)
