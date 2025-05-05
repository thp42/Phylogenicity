from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import re

# Load the Newick tree file
tree_file = "/home/pythagoras/Documents/PhD/Evolution/FL/ITPKA/All/ITPK_tax_tree.newick"
tree = Tree(tree_file, format=1)

# Function to create a mapping from FASTA files
def create_id_isoform_mapping(fasta_file, isoform_label):
    id_isoform_map = {}
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                unique_id = line.split('___')[0][1:]  # Extract unique ID
                organism_name = line.split('___organism_name___')[-1].split('___')[0]  # Extract organism name
                if unique_id not in id_isoform_map:
                    id_isoform_map[unique_id] = {'isoforms': set(), 'organism_name': organism_name}
                id_isoform_map[unique_id]['isoforms'].add(isoform_label)
    return id_isoform_map

# Create mappings for each isoform
itpka_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/FL/ITPKA/ITPKA/ITPKA.fasta", "ITPKA")
itpkb_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/FL/ITPKA/ITPKB/ITPKB.fasta", "ITPKB")
itpkc_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/FL/ITPKA/ITPKC/ITPKC.fasta", "ITPKC")

# Combine all mappings
id_isoform_map = {}
for d in [itpka_mapping, itpkb_mapping, itpkc_mapping]:
    for k, v in d.items():
        if k not in id_isoform_map:
            id_isoform_map[k] = {'isoforms': set(), 'organism_name': v['organism_name']}
        id_isoform_map[k]['isoforms'].update(v['isoforms'])

# Annotate the tree and rewrite the leaf names
for leaf in tree.iter_leaves():
    # Extract the unique ID part
    unique_id = leaf.name.split(' ')[0]
    
    data = id_isoform_map.get(unique_id, None)
    if data:
        isoforms = data['isoforms']
        organism_name = data['organism_name']
        leaf.add_features(isoforms=isoforms, organism_name=organism_name)
        leaf.name = f"{unique_id} {','.join(isoforms)} {organism_name}"
    else:
        # Extract organism name directly from the leaf name if not in isoform map
        name_parts = leaf.name.split(' ')
        if len(name_parts) > 1:
            organism_name = ' '.join(name_parts[1:])  # Everything after the unique ID
        else:
            organism_name = ""
        leaf.add_features(isoforms=None, organism_name=organism_name)
        leaf.name = f"{unique_id} {organism_name}"

# Define node style function
def my_layout(node):
    if node.is_leaf():
        label = node.name
        label_face = TextFace(label, fsize=14)
        if node.isoforms:
            if len(node.isoforms) > 1:
                label_face.fgcolor = "purple"  # Color for multiple isoforms
            else:
                isoform = list(node.isoforms)[0]
                if isoform == "ITPKA":
                    label_face.fgcolor = "red"
                elif isoform == "ITPKB":
                    label_face.fgcolor = "green"
                elif isoform == "ITPKC":
                    label_face.fgcolor = "blue"
        else:
            label_face.fgcolor = "black"
        node.add_face(label_face, column=0)

# Create a TreeStyle object
ts = TreeStyle()
ts.layout_fn = my_layout
ts.show_leaf_name = False
ts.mode = "c"  # Set the tree to circular mode

# Show the tree
tree.show(tree_style=ts)

# Save the modified tree as a tree file
modified_tree_file = "/home/pythagoras/Documents/PhD/Evolution/FL/ITPKA/All/ITPK_TaxAligned_tree_plot_with_isoforms.treefile"
tree.write(outfile=modified_tree_file)

# Print the file paths for confirmation
print(f"Modified tree with isoforms saved as treefile: {modified_tree_file}")
