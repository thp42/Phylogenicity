import json
import time
import concurrent.futures
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from urllib.error import HTTPError, URLError
import socket
from http.client import IncompleteRead  # Import IncompleteRead exception

# Your NCBI API key
API_KEY = "21994fff4aa18a624e4f2a5ccb4441f24e09"

def extract_motif_region(cds_seq, protein_seq, gene_motif, num_codons=60):
    motif_start = protein_seq.find(gene_motif)
    motif_end = motif_start + len(gene_motif)

    cds_motif_start = motif_start * 3
    cds_motif_end = motif_end * 3

    start = max(0, cds_motif_start - (num_codons * 3))
    end = min(len(cds_seq), cds_motif_end + (num_codons * 3))

    return cds_seq[start:end]

def search_ncbi_by_organism_and_gene(organism_name, gene_name, retries=3, delay=5):
    Entrez.email = "your.email@example.com"  # Replace with your actual email address

    # Add a delay before starting each request cycle to avoid overloading NCBI
    time.sleep(0.1)
    
    search_term = f"{organism_name}[Organism] AND {gene_name}[Gene]"
    gene_id = None

    for attempt in range(retries):
        try:
            search_handle = Entrez.esearch(db="gene", term=search_term, api_key=API_KEY)
            search_results = Entrez.read(search_handle)
            search_handle.close()

            if not search_results["IdList"]:
                print(f"No results found for organism '{organism_name}' and gene '{gene_name}'.")
                return None

            gene_id = search_results["IdList"][0]
            print(f"Gene ID for '{gene_name}' in '{organism_name}': {gene_id}")
            break
        except (HTTPError, URLError, socket.timeout, IncompleteRead) as e:
            print(f"Error during search: {e}. Retrying ({attempt + 1}/{retries})...")
            time.sleep(delay)
    
    if not gene_id:
        print("Failed to retrieve gene ID after multiple attempts.")
        return None

    accession_numbers = []
    for attempt in range(retries):
        try:
            fetch_handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="xml", api_key=API_KEY, post=True)
            gene_record = Entrez.read(fetch_handle)
            fetch_handle.close()

            for item in gene_record[0]["Entrezgene_locus"]:
                for sub_item in item["Gene-commentary_products"]:
                    accession_number = sub_item["Gene-commentary_accession"]
                    accession_numbers.append(accession_number)
            break
        except (HTTPError, URLError, socket.timeout, IncompleteRead) as e:
            print(f"Error fetching gene data: {e}. Retrying ({attempt + 1}/{retries})...")
            time.sleep(delay)

    if not accession_numbers:
        print("No mRNA or protein accession numbers found for this gene.")
        return None

    print(f"Accession numbers: {accession_numbers}")
    return accession_numbers

def fetch_fasta_and_cds_and_search_motif(accession_number, gene_motif, output_handle, motif_output_handle, organism_name, num_codons=60):
    Entrez.email = "your.email@example.com"
    
    time.sleep(0.1)

    try:
        fetch_handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text", api_key=API_KEY, post=True)
        record = SeqIO.read(fetch_handle, "genbank")
        fetch_handle.close()
    except HTTPError as e:
        print(f"HTTPError: {e.code} - {e.reason} while fetching nucleotide data for accession number {accession_number}")
        return False

    for feature in record.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(record.seq)
            protein_seq = cds_seq.translate(to_stop=True)

            if gene_motif in protein_seq:
                start_position = protein_seq.find(gene_motif)
                print(f"Gene motif '{gene_motif}' found in translated protein sequence for {accession_number} at position {start_position}!\n")

                # Write the full CDS sequence to the output FASTA
                record_id = f"{accession_number}_{gene_motif}_full_CDS"
                description = f"Full CDS sequence for {accession_number} with gene motif '{gene_motif}'"
                fasta_record = SeqRecord(cds_seq, id=record_id, description=description)
                SeqIO.write(fasta_record, output_handle, "fasta-2line")

                # Extract and write the motif region to the second output FASTA
                motif_region = extract_motif_region(cds_seq, protein_seq, gene_motif, num_codons)
                motif_record_id = f"{accession_number}_{gene_motif}_motif_region"
                motif_description = f"Motif region for {accession_number} with gene motif '{gene_motif}' and {num_codons} codons on each side"
                motif_fasta_record = SeqRecord(motif_region, id=motif_record_id, description=motif_description)
                SeqIO.write(motif_fasta_record, motif_output_handle, "fasta-2line")

                return True
            else:
                print(f"Gene motif '{gene_motif}' NOT found in translated protein sequence for {accession_number}.\n")
    
    return False

def process_sequence(header, sequence, output_fasta, motif_output_fasta, hardcoded_gene_name):
    organism_name = None
    gene_motif = None
    pub_gene_id = None

    annotation_start_index = header.find('{')
    annotation_end_index = header.find('}')
    
    if annotation_start_index != -1 and annotation_end_index != -1:
        annotation = header[annotation_start_index:annotation_end_index + 1]

        try:
            annotation_dict = json.loads(annotation)
            organism_name = annotation_dict.get("organism_name", None)
            pub_gene_id = annotation_dict.get("pub_gene_id", None)
        except json.JSONDecodeError:
            print(f"Error decoding JSON in header: {header.strip()}")

    if organism_name and sequence:
        gene_motif = sequence.strip()
        
        with open(output_fasta, 'a') as output_handle, open(motif_output_fasta, 'a') as motif_output_handle:
            accession_numbers = search_ncbi_by_organism_and_gene(organism_name, hardcoded_gene_name)
            if not accession_numbers and pub_gene_id:
                print(f"Falling back to using pub_gene_id: {pub_gene_id}")
                accession_numbers = search_ncbi_by_organism_and_gene(organism_name, pub_gene_id)

            if accession_numbers:
                for acc_num in accession_numbers:
                    success = fetch_fasta_and_cds_and_search_motif(acc_num, gene_motif, output_handle, motif_output_handle, organism_name)
                    if success:
                        return True
            else:
                print(f"No gene record found for '{hardcoded_gene_name}' or '{pub_gene_id}' in '{organism_name}'.")
    
    return False

def extract_info_and_search_ncbi(fasta_file, output_fasta, motif_output_fasta):
    successful_retrievals = 0
    failed_retrievals = 0
    hardcoded_gene_name = "ITPK"  # Hardcoded gene name to be used in all searches

    sequences_to_process = []

    with open(fasta_file, 'r') as file:
        header = None
        sequence = None

        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header and sequence:
                    sequences_to_process.append((header, sequence))
                header = line
                sequence = None
            else:
                sequence = line

        if header and sequence:
            sequences_to_process.append((header, sequence))

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(process_sequence, header, sequence, output_fasta, motif_output_fasta, hardcoded_gene_name): (header, sequence) for header, sequence in sequences_to_process}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                successful_retrievals += 1
            else:
                failed_retrievals += 1

    print(f"\nTotal successful retrievals: {successful_retrievals}")
    print(f"Total failed retrievals: {failed_retrievals}")
    
# Example usage of the function
fasta_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/Filtered_Flank_ITPKA_All.fasta"
output_fasta = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_FL_DNA_Sequences.fasta"
motif_output_fasta = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_Flank_DNA_Sequences.fasta"
extract_info_and_search_ncbi(fasta_file, output_fasta, motif_output_fasta)