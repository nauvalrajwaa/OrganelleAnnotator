#!/usr/bin/env python3
"""
plot_tree.py – Plot a phylogenetic tree from a Newick file.

Usage:
    python plot_tree.py <newick_file> <output_png>
"""

import sys
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo


def plot_tree(newick_file: str, output_png: str):
    if not os.path.exists(newick_file) or os.path.getsize(newick_file) == 0:
        print(f"WARNING: Tree file not found or empty: {newick_file}")
        # Create placeholder image
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No phylogenetic tree available",
                ha="center", va="center", fontsize=14, color="grey")
        ax.axis("off")
        plt.savefig(output_png, dpi=150, bbox_inches="tight")
        plt.close()
        return

    tree = Phylo.read(newick_file, "newick")
    n_terminals = tree.count_terminals()

    fig_height = max(6, n_terminals * 0.35)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    Phylo.draw(tree, axes=ax, do_show=False)

    ax.set_xlabel("Branch Length", fontsize=12)
    ax.set_title("Phylogenetic Tree (Maximum Likelihood)", fontsize=14, fontweight="bold")

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_png) or ".", exist_ok=True)
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Tree plot saved to {output_png}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_tree.py <newick_file> <output_png>")
        sys.exit(1)

    plot_tree(sys.argv[1], sys.argv[2])
