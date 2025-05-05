from ete3 import Tree, TreeStyle, NodeStyle, NCBITaxa
from Bio import SeqIO
import re

# Initialize NCBITaxa
ncbi = NCBITaxa()

# Update the local NCBI taxonomy database
print("Updating the local NCBI taxonomy database...")
ncbi.update_taxonomy_database()

# Function to parse FASTA file and extract taxonomic IDs and headers
def extract_taxids_and_headers_from_fasta(file_path):
    taxids = set()
    headers = {}
    for record in SeqIO.parse(file_path, "fasta"):
        header = record.description
        match = re.search(r'organism_taxid___(\d+)', header)
        if match:
            taxid = int(match.group(1))
            taxids.add(taxid)
            if taxid not in headers:
                headers[taxid] = []
            headers[taxid].append(header)
    return list(taxids), headers

# Function to get the taxonomic lineage from local NCBI database
def get_lineage(taxid):
    try:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        return [(str(tid), names[tid]) for tid in lineage] + [(str(taxid), names[taxid])]
    except Exception as e:
        print(f"Error fetching data for taxid {taxid}: {e}")
        return None

# Path to the FASTA file
fasta_file_path = "/home/pythagoras/Documents/PhD/Evolution/Proteins/CEFIP/CEFIP_aligned_protein_sequences.fasta"

# Extract taxonomic IDs and headers from the FASTA file
taxids, headers = extract_taxids_and_headers_from_fasta(fasta_file_path)

# Fetch taxonomic lineage for each taxid
lineage_dict = {}
for taxid in taxids:
    lineage = get_lineage(taxid)
    if lineage:
        lineage_dict[taxid] = lineage

# Function to get or create a node for a taxid
def get_or_create_node(taxid, name, taxid_to_node):
    if taxid in taxid_to_node:
        return taxid_to_node[taxid]
    else:
        node = Tree(name=f"{name} (TaxID {taxid})")
        taxid_to_node[taxid] = node
        return node

# Create the root of the tree
root = Tree(name="root")
taxid_to_node = {"1": root}

# Build the tree based on the taxonomic hierarchy
for taxid in taxids:
    if taxid in lineage_dict:
        lineage = lineage_dict[taxid]
        current_node = root
        for ancestor, name in lineage:
            node = get_or_create_node(ancestor, name, taxid_to_node)
            if node.up is None and node != current_node:
                current_node.add_child(node)
            current_node = node
        # Extract and add all sequences for the current taxid
        for header in headers[taxid]:
            unique_id = header.split('___')[0].replace('>', '')
            organism_name = re.search(r'organism_name___([^_]+)', header).group(1)
            leaf_node = Tree(name=f"{unique_id} {organism_name} (TaxID {taxid})")
            current_node.add_child(leaf_node)
    else:
        print(f"No lineage found for TaxID: {taxid}")

# Customize the tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = 'c'  # Circular mode

# Customize node style
nstyle = NodeStyle()
nstyle["fgcolor"] = "darkred"
nstyle["size"] = 10

for node in root.traverse():
    node.set_style(nstyle)

# Save the tree to a file in Newick format with node names
newick_file_path = "/home/pythagoras/Documents/PhD/Evolution/Proteins/CEFIP/CEFIP_tax_tree.newick"
root.write(outfile=newick_file_path, format=1)  # Ensure to use format=1 for extended Newick format with node names

# Print the path for confirmation
print(f"Tree saved as Newick file: {newick_file_path}")

# Plot the tree
print("Plotting the tree...")
root.show(tree_style=ts)
