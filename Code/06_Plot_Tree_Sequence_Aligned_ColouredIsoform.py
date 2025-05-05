from ete3 import Tree, TreeStyle, TextFace, NodeStyle

# Load the Newick tree file
tree_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_Isoform_aligned_protein_sequences.fasta.treefile"
tree = Tree(tree_file)

# Function to create a mapping from FASTA files
def create_id_isoform_mapping(fasta_file, isoform_label):
    id_isoform_map = {}
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                unique_id = line.split('___')[0][1:]  # Extract unique ID
                id_isoform_map[unique_id] = isoform_label
    return id_isoform_map

# Create mappings for each isoform
itpka_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA-Isoform.fasta", "ITPKA")
itpkb_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKB/ITPKB-Isoform.fasta", "ITPKB")
itpkc_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKC/ITPKC-Isoform.fasta", "ITPKC")
none_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/None/ITPK-NoneIsoform.fasta", "None")

#cefip_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/CEFIP/CEFIP_Isoform.fasta", "CEFIP")
#cefipnon_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/CEFIP/CEFIP_Nonisoform.fasta", "nonCEF")
#usp54_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/FL/USP/USP54/USP54.fasta", "USP54")

# Combine all mappings
#id_isoform_map = {**usp_mapping, **usp53_mapping, **usp54_mapping}
id_isoform_map = {**itpka_mapping, **itpkb_mapping, **itpkc_mapping, **none_mapping}

# Annotate the tree and rewrite the leaf names while preserving bootstrap values
for leaf in tree.iter_leaves():
    # Extract the unique ID part
    unique_id = leaf.name.split('___')[0]
    
    isoform = id_isoform_map.get(unique_id, None)
    leaf.add_features(isoform=isoform)
    if isoform:
        leaf.name = f"{unique_id} {isoform}"
    else:
        leaf.name = unique_id

# Define node style function
def my_layout(node):
    if node.is_leaf():
        label = node.name
        label_face = TextFace(label, fsize=14)
        if node.isoform:
            if node.isoform == "ITPKA":
                label_face.fgcolor = "red"
            elif node.isoform == "ITPKB":
                label_face.fgcolor = "green"
            elif node.isoform == "ITPKC":
                label_face.fgcolor = "blue"
        else:
            label_face.fgcolor = "black"
        node.add_face(label_face, column=0)

# Create a TreeStyle object
ts = TreeStyle()
ts.layout_fn = my_layout
ts.show_leaf_name = False
ts.show_branch_support = True  # Ensure bootstrap values are shown
ts.mode = "c"  # Set the tree to circular mode

# Show the tree
tree.show(tree_style=ts)

# Save the modified tree while preserving bootstrap values
modified_tree_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_aligned_protein_sequences_with_isoforms.fasta.treefile"
tree.write(outfile=modified_tree_file)

print(f"Modified tree with isoforms saved as treefile: {modified_tree_file}")
