#!/usr/bin/env python3
"""
Plot phylogenetic tree from Newick format to PNG
"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
import io

def plot_tree(newick_file, output_png):
    """
    Read Newick tree and plot to PNG
    """
    # Read the tree
    tree = Phylo.read(newick_file, "newick")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot the tree
    Phylo.draw(tree, axes=ax, do_show=False)
    
    # Customize
    ax.set_xlabel("Branch Length", fontsize=12)
    ax.set_title("Phylogenetic Tree (Maximum Likelihood)", fontsize=14, fontweight='bold')
    
    # Save
    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Tree plot saved to {output_png}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_tree.py <newick_file> <output_png>")
        sys.exit(1)
    
    newick_file = sys.argv[1]
    output_png = sys.argv[2]
    
    plot_tree(newick_file, output_png)
