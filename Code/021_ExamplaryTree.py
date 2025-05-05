from ete3 import Tree, NCBITaxa

# Initialize NCBITaxa
ncbi = NCBITaxa()

# Update the local NCBI taxonomy database
print("Updating the local NCBI taxonomy database...")
ncbi.update_taxonomy_database()

# Function to get the taxonomic lineage from local NCBI database
def get_lineage(taxid):
    try:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        return [(str(tid), names[tid]) for tid in lineage] + [(str(taxid), names[taxid])]
    except Exception as e:
        print(f"Error fetching data for taxid {taxid}: {e}")
        return None

# Hardcoded tips (Aves, Mammalia, Cladistia)
tip_names = [
    'Chondrichthyes',  # 1 sequence
    'Cladistia',       # 2 sequences
    'Coelacanthiformes',# 1 sequence
#    'Ceratodontiformes',# 1 sequence
    'Amphibia',        # 10 sequences
    'Actinopteri',     # 145 sequences
    'Mammalia',        # 173 sequences
    'Lepidosauria',    # 16 sequences
    'Crocodylia',      # 1 sequence
    'Testudines',      # 11 sequences
    'Aves',             # 59 sequences
    'Hyperoartia',
    'Saccharomyces cerevisiae',
    'Pseudomonas aeruginosa'
]

# Fetch taxids for the hardcoded tips
tip_taxids = [ncbi.get_name_translator([name])[name][0] for name in tip_names]

# Fetch taxonomic lineage for each taxid
lineage_dict = {}
for taxid in tip_taxids:
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
for taxid in tip_taxids:
    if taxid in lineage_dict:
        lineage = lineage_dict[taxid]
        current_node = root
        for ancestor, name in lineage:
            node = get_or_create_node(ancestor, name, taxid_to_node)
            if node.up is None and node != current_node:
                current_node.add_child(node)
            current_node = node
    else:
        print(f"No lineage found for TaxID: {taxid}")

# Save the tree to a Newick format file with node names
output_file = "phylogenetic_tree_with_names.nw"
root.write(outfile=output_file, format=1)  # format=1 includes internal node names

print(f"Tree saved to {output_file} with node names")
