#!/usr/bin/env python3
"""
calculate_rscu.py – Relative Synonymous Codon Usage analysis.

Extracts CDS features from a FASTA file, computes RSCU values, and generates
a heatmap plot.  Gene name resolution uses the shared gene_utils module.

Changes vs. original:
  - Removed embedded Ka/Ks logic (now in run_kaks_analysis.py)
  - Uses shared gene_utils.py for gene name normalisation
  - Uses MAFFT instead of MUSCLE (if alignment is needed)
  - Cleaner output formats

Usage:
    python calculate_rscu.py <input_fasta> <output_dir> [--genetic_code 2]
"""

import sys
import os
import argparse
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id

# Add scripts directory to path for local imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gene_utils import normalise_gene_name


# ── Codon tables ─────────────────────────────────────────────────────────────

def get_synonymous_codons(genetic_code: int) -> dict[str, list[str]]:
    """Return {amino_acid: [codons]} for the given genetic code table."""
    table = unambiguous_dna_by_id[genetic_code]
    aa_to_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in table.forward_table.items():
        aa_to_codons[aa].append(codon)
    return dict(aa_to_codons)


# ── RSCU calculation ─────────────────────────────────────────────────────────

def count_codons(sequence: str) -> dict[str, int]:
    """Count triplet codons in a nucleotide sequence."""
    counts: dict[str, int] = defaultdict(int)
    seq = sequence.upper().replace("-", "").replace("N", "")
    # Trim to multiple of 3
    length = (len(seq) // 3) * 3
    for i in range(0, length, 3):
        codon = seq[i : i + 3]
        if len(codon) == 3:
            counts[codon] += 1
    return dict(counts)


def compute_rscu(codon_counts: dict[str, int], genetic_code: int) -> dict[str, float]:
    """Compute RSCU values from codon counts.

    RSCU_ij = (X_ij / sum_j(X_ij)) * n_i
    where n_i = number of synonymous codons for amino acid i.
    """
    syn_codons = get_synonymous_codons(genetic_code)
    rscu: dict[str, float] = {}

    for aa, codons in syn_codons.items():
        total = sum(codon_counts.get(c, 0) for c in codons)
        n = len(codons)
        for c in codons:
            if total > 0:
                rscu[c] = (codon_counts.get(c, 0) / total) * n
            else:
                rscu[c] = 0.0

    return rscu


# ── Main ─────────────────────────────────────────────────────────────────────

def calculate_rscu(input_fasta: str, output_dir: str, genetic_code: int = 2):
    os.makedirs(output_dir, exist_ok=True)

    # Aggregate codon counts across all CDS records
    total_counts: dict[str, int] = defaultdict(int)
    per_gene: dict[str, dict[str, int]] = {}

    for record in SeqIO.parse(input_fasta, "fasta"):
        gene_name = normalise_gene_name(record.id)
        seq_str = str(record.seq).upper().replace("-", "")

        # Skip very short or non-CDS-like sequences
        if len(seq_str) < 30:
            continue

        codons = count_codons(seq_str)
        per_gene[gene_name] = codons
        for c, n in codons.items():
            total_counts[c] += n

    if not total_counts:
        print("WARNING: No valid CDS sequences found for RSCU analysis.")
        # Write empty outputs so Snakemake succeeds
        pd.DataFrame().to_csv(os.path.join(output_dir, "rscu.tsv"), sep="\t", index=False)
        return

    # Global RSCU
    rscu = compute_rscu(dict(total_counts), genetic_code)

    # Build output table
    syn_codons = get_synonymous_codons(genetic_code)
    rows = []
    for aa in sorted(syn_codons.keys()):
        for codon in sorted(syn_codons[aa]):
            rows.append({
                "AminoAcid": aa,
                "Codon": codon,
                "Count": total_counts.get(codon, 0),
                "RSCU": round(rscu.get(codon, 0.0), 4),
            })

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(output_dir, "rscu.tsv"), sep="\t", index=False)
    print(f"RSCU table saved ({len(df)} codons)")

    # ── Per-gene RSCU matrix for heatmap ──
    gene_rscu_data = {}
    for gene_name, gene_codons in per_gene.items():
        gene_rscu_data[gene_name] = compute_rscu(gene_codons, genetic_code)

    if gene_rscu_data:
        rscu_matrix = pd.DataFrame(gene_rscu_data).T.fillna(0)
        # Only keep codons that have some variation
        rscu_matrix = rscu_matrix.loc[:, rscu_matrix.std() > 0.01]

        if not rscu_matrix.empty:
            plt.figure(figsize=(max(14, rscu_matrix.shape[1] * 0.3), max(6, rscu_matrix.shape[0] * 0.4)))
            sns.heatmap(
                rscu_matrix, cmap="YlOrRd", linewidths=0.5, linecolor="white",
                xticklabels=True, yticklabels=True,
                cbar_kws={"label": "RSCU"},
            )
            plt.title("Relative Synonymous Codon Usage (RSCU)", fontsize=14, fontweight="bold")
            plt.xlabel("Codon")
            plt.ylabel("Gene")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "rscu_heatmap.png"), dpi=150, bbox_inches="tight")
            plt.close()
            print("RSCU heatmap saved")

    # ── Bar plot of global RSCU ──
    plt.figure(figsize=(16, 6))
    bar_df = df[df["Count"] > 0].copy()
    if not bar_df.empty:
        colors = sns.color_palette("husl", n_colors=len(bar_df["AminoAcid"].unique()))
        aa_colors = {aa: colors[i] for i, aa in enumerate(sorted(bar_df["AminoAcid"].unique()))}
        bar_colors = [aa_colors[aa] for aa in bar_df["AminoAcid"]]

        plt.bar(range(len(bar_df)), bar_df["RSCU"], color=bar_colors, edgecolor="white", linewidth=0.5)
        plt.xticks(range(len(bar_df)), bar_df["Codon"], rotation=90, fontsize=8)
        plt.axhline(y=1.0, color="red", linestyle="--", alpha=0.7, label="RSCU = 1.0")
        plt.ylabel("RSCU")
        plt.title("Global RSCU Values by Codon")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "rscu_barplot.png"), dpi=150, bbox_inches="tight")
        plt.close()
        print("RSCU bar plot saved")


def main():
    parser = argparse.ArgumentParser(description="Calculate RSCU from CDS FASTA")
    parser.add_argument("input_fasta", help="Input CDS FASTA file")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--genetic_code", type=int, default=2, help="NCBI genetic code (default: 2)")
    args = parser.parse_args()
    calculate_rscu(args.input_fasta, args.output_dir, args.genetic_code)


if __name__ == "__main__":
    main()
