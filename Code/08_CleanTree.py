from ete3 import Tree, TreeStyle, TextFace
import logging


def prune_branches(tree_file, threshold_distance, output_file):
    try:
        # Read the tree with standard Newick format
        tree = Tree(tree_file, format=0)  # Standard Newick format
    except Exception as e:
        print(f"Error reading tree: {e}")
        return None

    eliminated_tips = 0
    
    # Function to prune branches and count eliminated tips
    def prune_and_count(tree):
        eliminated_tips = 0
        nodes_to_prune = []
        for node in tree.traverse():
            if not node.is_root() and node.dist > threshold_distance:
                nodes_to_prune.append(node)
        
        for node in nodes_to_prune:
            eliminated_tips += len(node.get_leaves())
            node.detach()
        
        return eliminated_tips
    
    # Prune the tree and count the eliminated tips
    eliminated_tips = prune_and_count(tree)
    
    # Further clean up empty nodes
    def clean_empty_nodes(tree):
        for node in tree.traverse("postorder"):
            if not node.is_leaf() and len(node.children) == 0:
                node.delete(prevent_nondicotomic=False)
            elif not node.is_leaf() and len(node.children) == 1 and node.support is None and node.name == "":
                child = node.children[0]
                child.dist += node.dist
                node.delete(prevent_nondicotomic=False)
                
        # Additional check for improperly structured nodes
        for node in tree.traverse("postorder"):
            if node.name == "" and len(node.children) == 0:
                node.delete(prevent_nondicotomic=False)
    
    clean_empty_nodes(tree)
    
    # Save the pruned tree
    tree.write(outfile=output_file)
    
    # Print the number of eliminated tips
    print(f"Number of tips eliminated: {eliminated_tips}")
    
    return tree

# Parameters
tree_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_aligned_protein_sequences_with_isoforms.fasta.treefile"
threshold_distance = 0.8432119834893999  # Adjust this threshold as needed
output_file = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/All/ITPK_cleaned_tree.treefile"

# Prune branches and get the pruned tree
pruned_tree = prune_branches(tree_file, threshold_distance, output_file)

if pruned_tree:
    def visualize_tree(tree):
        # Define node style function
        def my_layout(node):
            if node.is_leaf():
                label = node.name
                label_face = TextFace(label, fsize=14)
                if "ITPKA" in node.name:
                    label_face.fgcolor = "red"
                elif "ITPKB" in node.name:
                    label_face.fgcolor = "green"
                elif "ITPKC" in node.name:
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

    # Visualize the pruned tree
    visualize_tree(pruned_tree)