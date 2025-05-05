import numpy as np
import pandas as pd
import subprocess
import tempfile
import os
import multiprocessing
from Bio import SeqIO
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
import re

# Define the regex pattern
regex_pattern = r'.V..[K][V]..F[E].'

# Define paths
fasta_file = '/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_OrthoDB.fasta'
output_fasta_file = '/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/RegEx_B.fasta'
iupred_path = "/home/pythagoras/Programs/iupred3"
psipred_path = "/home/pythagoras/Programs/psipred"

# Function to search for sequences using a regular expression
def search_with_regex(sequence, regex_pattern):
    matches = []
    for match in re.finditer(regex_pattern, sequence):
        start = match.start()
        end = match.end()
        fragment = sequence[start:end]
        matches.append((start, end, fragment))
    return matches

# Function to read sequences from a FASTA file
def read_fasta_sequences(fasta_file):
    fasta_sequences = list(SeqIO.parse(fasta_file, "fasta"))
    return fasta_sequences

# Function to cut sequence into fragments
def cut_sequence_into_fragments(sequence, fragment_size):
    sequence_length = len(sequence)
    fragments = []
    for i in range(sequence_length - fragment_size + 1):
        fragment = sequence[i:i + fragment_size]
        fragments.append(fragment)
    return fragments

# Function to get IUPRED and ANCHOR scores
def get_iupred_anchor_scores(sequence, seq_name, iupred_path):
    temp_fasta = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta')
    temp_fasta.write(f">{seq_name}\n{sequence}\n")
    temp_fasta.close()
    temp_fasta_name = temp_fasta.name
    
    try:
        # Run IUPred3 for disorder and ANCHOR scores
        iupred_cmd = [
            "python3",
            os.path.join(iupred_path, "iupred3.py"),
            temp_fasta_name,
            "long",
            "--anchor"
        ]
        result = subprocess.run(iupred_cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Error running IUPred3 for {seq_name}: {result.stderr}")
            return None, None

        iupred_scores, anchor_scores = [], []
        for line in result.stdout.splitlines():
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) == 4:
                iupred_scores.append(float(parts[2]))  # IUPRED score
                anchor_scores.append(float(parts[3]))  # ANCHOR score
        return iupred_scores, anchor_scores
    finally:
        os.unlink(temp_fasta_name)

# Function to get secondary structure probabilities from PSIPRED
def get_psipred_probabilities(seq_id, sequence, psipred_path):
    with tempfile.TemporaryDirectory() as temp_dir:
        original_dir = os.getcwd()
        try:
            os.chdir(temp_dir)
            temp_fasta_name = 'input.fasta'
            with open(temp_fasta_name, 'w') as temp_fasta:
                temp_fasta.write(f">{seq_id}\n{sequence}\n")
            
            env = os.environ.copy()
            env['PSIPRED_DATA'] = os.path.join(psipred_path, 'data')
            env['PATH'] += ':/usr/bin'  # Ensure BLAST is available in PATH
            
            runpsipred_cmd = [
                os.path.join(psipred_path, "runpsipred_single"),
                temp_fasta_name
            ]
            result = subprocess.run(runpsipred_cmd, capture_output=True, text=True, env=env)
            
            if result.returncode != 0:
                print(f"Error running PSIPRED for {seq_id}: {result.stderr}")
                return None

            ss2_file = 'input.ss2'
            if not os.path.exists(ss2_file):
                return None

            probabilities = []
            with open(ss2_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.strip() == '':
                        continue
                    parts = line.strip().split()
                    if len(parts) != 6:
                        continue
                    prob_c = float(parts[3])  # Coil
                    prob_h = float(parts[4])  # Helix
                    prob_e = float(parts[5])  # Strand
                    probabilities.append({'Prob_C': prob_c, 'Prob_H': prob_h, 'Prob_E': prob_e})
            return probabilities
        finally:
            os.chdir(original_dir)

# Function to process a single sequence record using regex instead of PSSM
def process_sequence_record_with_regex(args):
    seq_record, iupred_path, psipred_path, iupred_cutoff, anchor_cutoff, ss_cutoff, flank_size, regex_pattern = args
    sequence_name = seq_record.id
    full_header = seq_record.description
    sequence = str(seq_record.seq)

    # Get IUPRED, ANCHOR, and PSIPRED scores
    iupred_scores, anchor_scores = get_iupred_anchor_scores(sequence, sequence_name, iupred_path)
    if iupred_scores is None or anchor_scores is None:
        return []

    psipred_probs = get_psipred_probabilities(sequence_name, sequence, psipred_path)
    if psipred_probs is None:
        return []

    results = []

    # Search for fragments matching the regular expression pattern
    matches = search_with_regex(sequence, regex_pattern)

    for match in matches:
        start, end, fragment = match
        start_flank = max(0, start - flank_size)
        end_flank = min(len(sequence), end + flank_size)
        
        # Exclude the fragment itself from the IUPRED score calculation
        left_flank = iupred_scores[start_flank:start]  # Flanking region before the fragment
        right_flank = iupred_scores[end:end_flank]  # Flanking region after the fragment
        flanking_region_scores = left_flank + right_flank

        if len(flanking_region_scores) > 0:
            mean_iupred = np.mean(flanking_region_scores)
        else:
            mean_iupred = 1  # Handle cases where no flanking region exists

        mean_anchor = np.mean(anchor_scores[start:end])

        mean_coil = np.mean([p['Prob_C'] for p in psipred_probs[start:end]])

        mean_helix = np.mean([p['Prob_H'] for p in psipred_probs[start:end]])

        mean_strand = np.mean([p['Prob_E'] for p in psipred_probs[start:end]])

        if (mean_iupred >= iupred_cutoff and mean_anchor >= anchor_cutoff and
            mean_coil < ss_cutoff['coil'] and mean_helix >= ss_cutoff['helix'] and mean_strand < ss_cutoff['strand']):
            fasta_header = (f">{full_header} pos:{start+1} "
                            f"IUPRED_Score:{mean_iupred:.2f} ANCHOR_Score:{mean_anchor:.2f} "
                            f"Coil:{mean_coil:.2f} Helix:{mean_helix:.2f} Strand:{mean_strand:.2f}")
            
            results.append((fasta_header, fragment))
    
    return results

# Parallel processing setup
sequences_to_score = read_fasta_sequences(fasta_file)
pssm_cutoff, iupred_cutoff, anchor_cutoff = 10, 0.43686, 0.41004
ss_cutoff = {'coil': 0.50391, 'helix': 0.0, 'strand': 1.0}
flank_size = 60

args = [(seq_record, iupred_path, psipred_path, iupred_cutoff, anchor_cutoff, ss_cutoff, flank_size, regex_pattern) for seq_record in sequences_to_score]

num_cores = cpu_count() - 1
with Pool(num_cores) as pool:
    all_results = pool.map(process_sequence_record_with_regex, args)

# Write results to output file
with open(output_fasta_file, "w") as output_conn:
    for result in all_results:
        for header, fragment in result:
            output_conn.write(f"{header}\n{fragment}\n")

print(f"Results written to {output_fasta_file}")
