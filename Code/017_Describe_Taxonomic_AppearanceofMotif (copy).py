from ete3 import Tree, TreeStyle, TextFace, NodeStyle, NCBITaxa
from collections import Counter, defaultdict

# Initialize NCBITaxa
ncbi = NCBITaxa()

# Function to read motifs from a file
def read_motifs(fasta_file):
    motifs = {}
    with open(fasta_file, 'r') as file:
        motif_name = None
        motif_seq = None
        for line in file:
            if line.startswith('>'):
                if motif_name and motif_seq:
                    motifs[motif_name] = {'name': motif_name, 'sequence': motif_seq}
                motif_name = line.strip()[1:]  # Extract unique ID
                motif_seq = ''
            else:
                motif_seq += line.strip()
        if motif_name and motif_seq:
            motifs[motif_name] = {'name': motif_name, 'sequence': motif_seq}
    return motifs

# Function to create a mapping from FASTA files
def create_id_isoform_mapping(fasta_file, isoform_label):
    id_isoform_map = {}
    with open(fasta_file, 'r') as file:
        unique_id = None
        for line in file:
            if line.startswith('>'):
                if unique_id:
                    if unique_id not in id_isoform_map:
                        id_isoform_map[unique_id] = {'isoforms': set(), 'motifs': set()}
                    id_isoform_map[unique_id]['isoforms'].add(isoform_label)
                unique_id = line.split('___')[0][1:]  # Extract unique ID
        if unique_id:
            if unique_id not in id_isoform_map:
                id_isoform_map[unique_id] = {'isoforms': set(), 'motifs': set()}
            id_isoform_map[unique_id]['isoforms'].add(isoform_label)
    return id_isoform_map

# Function to get taxonomy name
def get_taxonomy_name(taxid):
    if taxid is None:
        return "Unknown"
    names = ncbi.get_taxid_translator([taxid])
    return names.get(taxid, "Unknown")

def get_lineage(taxid):
    try:
        lineage = ncbi.get_lineage(int(taxid))
        names = ncbi.get_taxid_translator(lineage)
        return [names[taxid] for taxid in lineage]
    except:
        return []
    
def find_common_ancestor_taxonomy(lineages):
    if not lineages:
        return "Unknown"
    common_ancestor = set(lineages[0])
    for lineage in lineages[1:]:
        common_ancestor.intersection_update(lineage)
    if not common_ancestor:
        return "Unknown"
    valid_ancestors = [x for x in common_ancestor if len(x) > 1]  # Filter out single-letter ancestors
    if not valid_ancestors:
        return "Unknown"
    return max(valid_ancestors, key=lambda x: max(lineage.index(x) if x in lineage else -1 for lineage in lineages))

def find_common_ancestor(lineages):
    if not lineages:
        return None
    common_ancestors = set(lineages[0])
    for lineage in lineages[1:]:
        common_ancestors.intersection_update(lineage)
    return max(common_ancestors, key=lambda x: max(lineage.index(x) if x in lineage else -1 for lineage in lineages)) if common_ancestors else None

# Modify the analyze_motif_origins function to return more detailed information
def analyze_motif_origins(tree, id_isoform_map):
    isoform_lineages = defaultdict(list)
    motif_isoforms = defaultdict(set)
    
    for leaf in tree.iter_leaves():
        unique_id = leaf.name.split('___')[0]
        taxid = unique_id.split('_')[0]
        data = id_isoform_map.get(unique_id, {'isoforms': set(), 'motifs': set()})
        
        lineage = get_lineage(taxid)
        
        for isoform in data['isoforms']:
            if data['motifs']:  # Only consider this lineage if the sequence has any motifs
                isoform_lineages[isoform].append(lineage)
                for motif_name, _ in data['motifs']:
                    motif_isoforms[motif_name].add(isoform)
    
    isoform_origins = {}
    for isoform, lineages in isoform_lineages.items():
        common_ancestor = find_common_ancestor_taxonomy(lineages)
        isoform_origins[isoform] = common_ancestor
    
    return isoform_origins, motif_isoforms

# Function to get common taxonomy
def get_common_taxonomy(node, lineage_info):
    descendant_lineages = [lineage_info[leaf.name.split('_')[0]] for leaf in node.iter_leaves() if leaf.name.split('_')[0] in lineage_info]
    if not descendant_lineages:
        return "Unknown"
    common_ancestor_taxid = find_common_ancestor(descendant_lineages)
    return get_taxonomy_name(common_ancestor_taxid)

# Function to annotate tree with taxonomy
def annotate_tree_with_taxonomy(tree, lineage_info, tax_info):
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            tax_id = node.name.split('_')[0]
            node.add_feature("taxonomy", tax_info.get(tax_id, "Unknown"))
        else:
            node.add_feature("taxonomy", get_common_taxonomy(node, lineage_info))
    return tree

# Function to propagate isoform and motif information
def propagate_isoform_motif_info(tree, id_isoform_map):
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            unique_id = node.name.split('_')[0]
            data = id_isoform_map.get(unique_id, {'isoforms': set(), 'motifs': set()})
            node.add_feature("isoforms", data['isoforms'])
            node.add_feature("motifs", data['motifs'])
            node.add_feature("isoform_counts", Counter(data['isoforms']))
            node.add_feature("motif_counts", Counter(motif for motif, _ in data['motifs']))
        else:
            node.add_feature("isoform_counts", sum((child.isoform_counts for child in node.children), Counter()))
            node.add_feature("motif_counts", sum((child.motif_counts for child in node.children), Counter()))
    return tree

# Function to find first appearance of motifs
def find_motif_first_appearance(tree):
    motif_appearances = {}
    for node in tree.traverse("postorder"):
        for motif in node.motif_counts:
            if motif not in motif_appearances:
                motif_appearances[motif] = node
    return motif_appearances

# Updated layout function for the tree
def isoform_motif_layout(node, id_isoform_map):
    if node.is_leaf():
        # Extract the unique ID part
        unique_id = node.name.split('___')[0]
        
        # Get isoform and motif information
        data = id_isoform_map.get(unique_id, {'isoforms': set(), 'motifs': set()})
        isoforms = data['isoforms']
        motifs = data['motifs']
        
        # Collect the motif sequences
        motif_sequences = [seq for _, seq in motifs]
        
        # Create label
        label = f"{unique_id} {','.join(isoforms)} {'; '.join(['Motif: ' + seq for seq in motif_sequences])}"
        
        # Create and add label face
        label_face = TextFace(label, fsize=8)
        if isoforms:
            if "Shroom2" in isoforms:
                label_face.fgcolor = "red"
            elif "Shroom3" in isoforms:
                label_face.fgcolor = "green"
            elif "Shroom4" in isoforms:
                label_face.fgcolor = "blue"
            else:
                label_face.fgcolor = "purple"  # Color for multiple isoforms
        else:
            label_face.fgcolor = "black"
        node.add_face(label_face, column=0)
        
        # Add motif information
        if motifs:
            motif_names = [name for name, _ in motifs]
            motif_face = TextFace(" | ".join(motif_names), fsize=6)
            motif_face.fgcolor = "gray"
            node.add_face(motif_face, column=1)
    else:
        if hasattr(node, "taxonomy"):
            node.add_face(TextFace(node.taxonomy, fsize=6, fgcolor="green"), column=0, position="branch-top")
        
        isoform_info = ", ".join([f"{iso}: {count}" for iso, count in node.isoform_counts.items()])
        if isoform_info:
            node.add_face(TextFace(isoform_info, fsize=6, fgcolor="purple"), column=1, position="branch-top")

# New function to check if a taxon is in a lineage
def is_taxon_in_lineage(taxon_id, lineage):
    taxon_lineage = get_lineage(taxon_id)
    return set(lineage).intersection(set(taxon_lineage)) != set()

# New function to analyze unmapped motifs
def analyze_unmapped_motifs(unmapped_motifs, isoform_origins):
    print("\nAnalyzing unmapped motifs...")
    for isoform, ancestor in isoform_origins.items():
        try:
            # Translate common ancestor name to taxid
            ancestor_taxid = ncbi.get_name_translator([ancestor])[ancestor][0]
            print(f"\nAnalyzing unmapped motifs for {isoform} (Common ancestor: {ancestor}, Taxid: {ancestor_taxid})...")

            motifs_outside_lineage = []
            for motif in unmapped_motifs:
                taxid = int(motif.split('_')[0])
                # Get the lineage of the motif taxon
                motif_lineage = get_lineage(taxid)
                
                #print(f"Motif: {motif}, Taxid: {taxid}")
                #print(f"Motif Lineage (taxids): {motif_lineage}")
                
                # Get lineage names with error handling
                lineage_names = []
                for tid in motif_lineage:
                    if str(tid).isdigit():
                        try:
                            name = get_taxonomy_name(tid)
                            lineage_names.append(name)
                        except Exception as e:
                            print(f"Error getting name for taxid {tid}: {str(e)}")
                    else:
                        lineage_names.append(str(tid))

                if ancestor not in lineage_names:
                    motifs_outside_lineage.append(motif)

            if motifs_outside_lineage:
                print(f"The following unmapped motifs belong to taxa outside the lineage of {isoform}'s common ancestor:")
                for motif in motifs_outside_lineage:
                    taxid = int(motif.split('_')[0])
                    try:
                        taxon_name = get_taxonomy_name(taxid)
                    except Exception as e:
                        taxon_name = f"Unknown (Error: {str(e)})"
                    print(f"  - {motif} ({taxon_name})")
            else:
                print(f"All unmapped motifs belong to taxa within or including the lineage of {isoform}'s common ancestor.")
        except Exception as e:
            print(f"Error during unmapped motif analysis for {isoform}: {str(e)}")
            import traceback
            traceback.print_exc()


# Main execution
if __name__ == "__main__":
    # Load the tree
    tree_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_Isoform_aligned_protein_sequences.fasta.treefile"
    tree = Tree(tree_file, format=1)

    # Read motifs
    motif_file = '/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_Motif.fasta'
    motifs = read_motifs(motif_file)

    # Create mappings for each isoform
    itpka_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKA/ITPKA_ownAnnotation.fasta", "ITPKA")
    itpkb_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKB/ITPKB_ownAnnotation.fasta", "ITPKB")
    itpkc_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKC/ITPKC_ownAnnotation.fasta", "ITPKC")
    #cefip_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/Proteins/CEFIP/CEFIP_aligned_protein_sequences.fasta", "CEFIP")
    #shroom3_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/FL/Shroom3H2/Shroom3/Shroom3_ownAnnotation.fasta", "Shroom3")
    #shroom4_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/FL/Shroom3H2/Shroom4/Shroom4_ownAnnotation.fasta", "Shroom4")
    #shroom_mapping = create_id_isoform_mapping("/home/pythagoras/Documents/PhD/Evolution/FL/Shroom3H2/Shroom3_aligned_protein_sequences.fasta", "ShroomAll")
    # Combine all mappings and add motifs
    id_isoform_map = {}
    mapped_motifs = set()
    for d in [itpka_mapping, itpkb_mapping, itpkc_mapping]:
        for k, v in d.items():
            if k not in id_isoform_map:
                id_isoform_map[k] = {'isoforms': set(), 'motifs': set()}
            id_isoform_map[k]['isoforms'].update(v['isoforms'])
            if k in motifs:
                id_isoform_map[k]['motifs'].add((motifs[k]['name'], motifs[k]['sequence']))
                mapped_motifs.add(motifs[k]['name'])

    # Check for unmapped motifs and print warning
    unmapped_motifs = set(motifs.keys()) - mapped_motifs
    if unmapped_motifs:
        print("Warning: The following motifs were not found in any sequence:")
        for motif in unmapped_motifs:
            print(f"  - {motif}")

    # Get taxonomy information
    tax_info = {}
    lineage_info = {}
    for leaf in tree.iter_leaves():
        parts = leaf.name.split('_')
        tax_id = parts[0]
        if tax_id.isdigit():
            lineage = ncbi.get_lineage(int(tax_id))
            if lineage:
                names = ncbi.get_taxid_translator(lineage)
                lineage_names = [names[taxid] for taxid in lineage]
                tax_info[tax_id] = lineage_names[-1]
                lineage_info[tax_id] = lineage
            else:
                tax_info[tax_id] = "Unknown"
                lineage_info[tax_id] = []

    # Annotate tree
    tree = annotate_tree_with_taxonomy(tree, lineage_info, tax_info)
    tree = propagate_isoform_motif_info(tree, id_isoform_map)

    # Analyze motif origins
    isoform_origins, motif_isoforms = analyze_motif_origins(tree, id_isoform_map)

    # Print motif origins for each isoform
    print("Common ancestors for sequences containing motifs in each isoform:")
    for isoform, common_ancestor in isoform_origins.items():
        print(f"{isoform}: Common ancestor: {common_ancestor}")

    # Analyze unmapped motifs
    analyze_unmapped_motifs(unmapped_motifs, isoform_origins)

    # Create TreeStyle
    ts = TreeStyle()
    ts.layout_fn = lambda node: isoform_motif_layout(node, id_isoform_map)
    ts.show_leaf_name = False
    ts.mode = "c"  # Set to circular mode
    ts.arc_start = -180  # Start from the top
    ts.arc_span = 360  # Full circle
    ts.title.add_face(TextFace("ITPK Isoform and Motif Phylogenetic Tree", fsize=20), column=0)

    # Show the tree
    tree.show(tree_style=ts)

    # Save the tree as an image
    tree.render("/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_Isoform_Motif_annotated_tree.png", tree_style=ts, w=1000, units="px")
    print("Annotated tree image saved as: /home/pythagoras/Documents/PhD/Evolution/FL/Shroom3H2/Shroom3_Isoform_Motif_annotated_tree.png")

    # Save the tree as a treefile
    tree.write(format=1, outfile="/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_Isoform_Motif_annotated_tree.treefile")
    print("Annotated tree saved as: /home/pythagoras/Documents/PhD/Evolution/FL/Shroom3H2/Shroom3H2_Isoform_Motif_annotated_tree.treefile")