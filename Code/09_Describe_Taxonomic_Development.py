from ete3 import NCBITaxa, Tree, TreeStyle, TextFace, NodeStyle, RectFace
import logging
from collections import Counter

# Initialize logging
#logging.basicConfig(filename='/home/pythagoras/Documents/PhD/Evolution/ITPKA/All/debug_log.txt', level=logging.DEBUG, filemode='w')

# Initialize NCBITaxa instance
ncbi = NCBITaxa()

# Load the pruned tree
seq_tree_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_cleaned_tree.treefile"
seq_tree = Tree(seq_tree_file)

# Extract taxonomic IDs from the leaves and fetch taxonomic information
tax_info = {}
lineage_info = {}
for leaf in seq_tree.iter_leaves():
    parts = leaf.name.split('_')
    tax_id = parts[0]
    
    if tax_id.isdigit():
        try:
            lineage = ncbi.get_lineage(int(tax_id))
            if lineage:
                names = ncbi.get_taxid_translator(lineage)
                lineage_names = [names[taxid] for taxid in lineage]
                tax_info[tax_id] = lineage_names[-1]
                lineage_info[tax_id] = lineage
            else:
                tax_info[tax_id] = "Unknown"
                lineage_info[tax_id] = []
        except Exception as e:
            tax_info[tax_id] = "Unknown"
            lineage_info[tax_id] = []
            logging.debug(f"Error fetching lineage for Tax ID {tax_id}: {e}")

def find_common_ancestor(lineages):
    if not lineages:
        return None
    common_ancestors = set(lineages[0])
    for lineage in lineages[1:]:
        common_ancestors.intersection_update(lineage)
    return max(common_ancestors, key=lambda x: max(lineage.index(x) if x in lineage else -1 for lineage in lineages)) if common_ancestors else None

def get_taxonomy_name(taxid):
    if taxid is None:
        return "Unknown"
    names = ncbi.get_taxid_translator([taxid])
    return names.get(taxid, "Unknown")

def get_common_taxonomy(node, lineage_info):
    descendant_lineages = [lineage_info[leaf.name.split('_')[0]] for leaf in node.iter_leaves() if leaf.name.split('_')[0] in lineage_info]
    if not descendant_lineages:
        return "Unknown"
    common_ancestor_taxid = find_common_ancestor(descendant_lineages)
    return get_taxonomy_name(common_ancestor_taxid)

def annotate_tree_with_taxonomy(tree, lineage_info, tax_info):
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            tax_id = node.name.split('_')[0]
            node.add_feature("taxonomy", tax_info.get(tax_id, "Unknown"))
        else:
            node.add_feature("taxonomy", get_common_taxonomy(node, lineage_info))
    return tree

def get_isoform(node_name):
    parts = node_name.split()
    return parts[-1] if parts and parts[-1] in ['ITPKA', 'ITPKB', 'ITPKC', 'None'] else None

def propagate_isoform_counts(tree):
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            isoform = get_isoform(node.name)
            node.add_feature("isoform_counts", Counter({isoform: 1}) if isoform else Counter({"Unknown": 1}))
            node.add_feature("total_leaves", 1)
        else:
            node.add_feature("isoform_counts", sum((child.isoform_counts for child in node.children), Counter()))
            node.add_feature("total_leaves", sum(child.total_leaves for child in node.children))
    return tree

def find_isoform_splits(tree):
    return [node for node in tree.traverse() if not node.is_leaf() and len(node.children) > 1 and 
            len(set.union(*[set(child.isoform_counts.keys()) - {"Unknown"} for child in node.children])) > 1]

def cluster_and_taxonomy_layout(node):
    if node.is_leaf():
        label = node.name
        label_face = TextFace(label, fsize=8)
        isoform = get_isoform(node.name)
        label_face.fgcolor = {"ITPKA": "red", "ITPKB": "green", "ITPKC": "blue"}.get(isoform, "black")
        node.add_face(label_face, column=0)
    else:
        if hasattr(node, "taxonomy"):
            node.add_face(TextFace(node.taxonomy, fsize=6, fgcolor="green"), column=0, position="branch-top")
        if hasattr(node, "isoform_counts"):
            isoform_info = ", ".join([f"{iso}: {count}" for iso, count in node.isoform_counts.items() if iso != "Unknown"])
            if isoform_info:
                node.add_face(TextFace(isoform_info, fsize=6, fgcolor="purple"), column=1, position="branch-top")
        if len(set(node.isoform_counts.keys()) - {"Unknown"}) > 1:
            node.img_style["bgcolor"] = "LightYellow"

def summarize_highly_significant_isoform_divergence(splits, min_proportion=0.01, max_events=10, min_isoform_count=50, min_separation_ratio=0.7, min_child_isoform_count=100):
    def is_highly_significant_divergence(node):
        parent_isoforms = {k: v for k, v in node.isoform_counts.items() if k != 'Unknown' and v >= min_isoform_count}
        if len(parent_isoforms) < 2:
            return False
        child_isoforms = [{k: v for k, v in child.isoform_counts.items() if k != 'Unknown'} for child in node.children]
        if any(sum(child.values()) < min_child_isoform_count for child in child_isoforms):
            return False
        significant_separations = sum(1 for isoform in parent_isoforms if max(child.get(isoform, 0) for child in child_isoforms) / parent_isoforms[isoform] >= min_separation_ratio)
        unique_isoform_transitions = len([isoform for isoform in parent_isoforms if any(child.get(isoform, 0) > 0 for child in child_isoforms)])
        return significant_separations >= 2 and unique_isoform_transitions > 1

    total_leaves = splits[0].total_leaves
    divergence_events = [node for node in splits if is_highly_significant_divergence(node) and node.total_leaves / total_leaves >= min_proportion]
    divergence_events.sort(key=lambda x: x.total_leaves, reverse=True)

    print("\nHighly Significant ITPK Isoform Divergence Events:")
    for i, node in enumerate(divergence_events[:max_events], 1):
        print(f"Divergence Event {i}:")
        print(f"  Taxonomy: {getattr(node, 'taxonomy', 'Unknown')}")
        print(f"  Parent Isoforms: {dict(node.isoform_counts)}")
        print(f"  Child Distributions:")
        for j, child in enumerate(node.children, 1):
            print(f"    Child {j} ({getattr(child, 'taxonomy', 'Unknown')}): {dict(child.isoform_counts)}")
        print(f"  Total leaves: {node.total_leaves}")
        print(f"  Proportion of total tree: {node.total_leaves / total_leaves:.2%}")
        print("  Evolutionary Significance:")
        interpret_highly_significant_divergence(node, total_leaves, min_separation_ratio)
        print()

def interpret_highly_significant_divergence(node, total_leaves, min_separation_ratio):
    parent_isoforms = {k: v for k, v in node.isoform_counts.items() if k != 'Unknown'}
    child_isoforms = [{k: v for k, v in child.isoform_counts.items() if k != 'Unknown'} for child in node.children]
    
    print(f"    - This event represents a major divergence in ITPK evolution, affecting {node.total_leaves / total_leaves:.2%} of all analyzed sequences.")
    print(f"    - The ancestral lineage at the {getattr(node, 'taxonomy', 'Unknown')} level had multiple ITPK isoforms: {', '.join(parent_isoforms.keys())}.")
    
    for isoform in parent_isoforms:
        parent_count = parent_isoforms[isoform]
        child_counts = [child.get(isoform, 0) for child in child_isoforms]
        max_child_count = max(child_counts)
        if max_child_count / parent_count >= min_separation_ratio:
            child_index = child_counts.index(max_child_count)
            child_taxonomy = getattr(node.children[child_index], 'taxonomy', 'Unknown')
            print(f"    - A significant portion of {isoform} ({max_child_count}/{parent_count}, {max_child_count/parent_count:.2%}) transitioned to the {child_taxonomy} lineage.")
    
    if any('Unknown' in child for child in child_isoforms):
        print("    - The presence of unclassified sequences suggests potential novel ITPK variants or intermediates in this evolutionary transition.")

# Main execution
seq_tree = annotate_tree_with_taxonomy(seq_tree, lineage_info, tax_info)
seq_tree = propagate_isoform_counts(seq_tree)

total_leaves = seq_tree.total_leaves
isoform_counts = seq_tree.isoform_counts

print("Overall Isoform Distribution:")
for isoform, count in isoform_counts.items():
    print(f"{isoform}: {count} ({count/total_leaves*100:.2f}%)")
print(f"Non-classified: {total_leaves - sum(isoform_counts.values())} ({(total_leaves - sum(isoform_counts.values()))/total_leaves*100:.2f}%)")

splits = find_isoform_splits(seq_tree)
summarize_highly_significant_isoform_divergence(splits, min_proportion=0.01, max_events=10, min_isoform_count=200, min_separation_ratio=0.4)
#summarize_highly_significant_isoform_divergence(splits, min_proportion=0.01, max_events=10, min_isoform_count=20, min_separation_ratio=0.1)

# Create and configure TreeStyle
seq_ts = TreeStyle()
seq_ts.layout_fn = cluster_and_taxonomy_layout
seq_ts.show_leaf_name = False
seq_ts.mode = "c"
seq_ts.title.add_face(TextFace("USP Phylogenetic Tree", fsize=20), column=0)

# Add legend
legend_faces = [
    TextFace("ITPKA", fsize=10, fgcolor="red"),
    TextFace("ITPKB", fsize=10, fgcolor="green"),
    TextFace("ITPKC", fsize=10, fgcolor="blue"),
    TextFace("Isoform Split", fsize=10),
    TextFace("Taxonomy", fsize=10, fgcolor="green"),
    TextFace("Isoform Distribution", fsize=10, fgcolor="purple")
]

for i, face in enumerate(legend_faces):
    seq_ts.legend.add_face(face, column=0)
    if face.text == "Isoform Split":
        bg_face = RectFace(60, 20, "LightYellow", "LightYellow")
        seq_ts.legend.add_face(bg_face, column=0)

# Show the annotated tree
seq_tree.show(tree_style=seq_ts)