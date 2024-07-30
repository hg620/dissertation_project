from Bio import Phylo

tree_path = "****/tree_double_precision.nwk"
tree = Phylo.read(tree_path, "newick")

# Function to remove support values
def remove_support_values(tree):
    for clade in tree.find_clades():
        clade.confidence = None

# Function to prune non-binary nodes only at the tips
def prune_large_tip_nodes(tree):
    for clade in tree.find_clades(order="postorder"):
        if not clade.is_terminal() and all(child.is_terminal() for child in clade.clades) and len(clade.clades) > 2:
            print(f"Pruning clade with {len(clade.clades)} tip descendants: {clade.clades}")
            # Keep only the first two clades (or modify the logic as needed)
            clade.clades = clade.clades[:2]
            print(f"Remaining clades: {clade.clades}")

print("Original Tree:")
Phylo.draw_ascii(tree)

#modifying the tree
remove_support_values(tree)
prune_large_tip_nodes(tree)

# Printing tree after modifications
print("Modified Tree:")
Phylo.draw_ascii(tree)

# Saving the modified tree back to a new file with nine decimal places
output_path = "*****/tree_double_precision_binary.nwk"
with open(output_path, "w") as outfile:
    Phylo.write(tree, outfile, "newick", format_branch_length="%.9f")
