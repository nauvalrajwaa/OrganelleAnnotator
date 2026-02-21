#!/usr/bin/env python3
"""
analyze_composition.py – GC content and amino acid composition analysis.

Reads a CDS FASTA, computes per-gene GC content and aggregate amino acid
frequencies, and generates plots.

Usage:
    python analyze_composition.py <input_fasta> <output_dir> [--genetic_code 2]
"""

import sys
import os
import argparse

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gene_utils import normalise_gene_name


def analyze_composition(input_fasta: str, output_dir: str, genetic_code: int = 2):
    os.makedirs(output_dir, exist_ok=True)

    records = []
    aa_counts: dict[str, int] = {}

    for record in SeqIO.parse(input_fasta, "fasta"):
        gene_name = normalise_gene_name(record.id)
        seq_str = str(record.seq).upper().replace("-", "")

        if len(seq_str) < 30:
            continue

        # GC content
        gc = gc_fraction(seq_str) * 100
        records.append({
            "Gene": gene_name,
            "Length": len(seq_str),
            "GC_Content": gc,
        })

        # Translate to protein
        try:
            # Trim to codon-complete
            trim_len = (len(seq_str) // 3) * 3
            if trim_len < 3:
                continue
            protein = record.seq[:trim_len].translate(table=genetic_code, to_stop=False)
            for aa in str(protein):
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
        except Exception as exc:
            print(f"Translation warning for {gene_name}: {exc}")

    # ── GC Content ──
    if records:
        df_gc = pd.DataFrame(records)
        df_gc.to_csv(os.path.join(output_dir, "gc_content.tsv"), sep="\t", index=False)

        plt.figure(figsize=(max(10, len(records) * 0.5), 6))
        sns.barplot(data=df_gc, x="Gene", y="GC_Content", hue="Gene",
                    palette="viridis", legend=False)
        plt.axhline(y=df_gc["GC_Content"].mean(), color="red", linestyle="--",
                     alpha=0.7, label=f"Mean GC: {df_gc['GC_Content'].mean():.1f}%")
        plt.title("GC Content per Gene")
        plt.ylabel("GC Content (%)")
        plt.xticks(rotation=45, ha="right")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "gc_content_plot.png"), dpi=150, bbox_inches="tight")
        plt.close()
        print("Generated GC content plot.")
    else:
        print("No CDS records found for GC analysis.")

    # ── Amino Acid Composition ──
    if aa_counts:
        total_aa = sum(aa_counts.values())
        aa_data = []
        for aa in sorted(aa_counts.keys()):
            if aa == "*":
                continue
            pct = (aa_counts[aa] / total_aa) * 100
            aa_data.append({"AminoAcid": aa, "Percentage": pct, "Count": aa_counts[aa]})

        df_aa = pd.DataFrame(aa_data)
        df_aa.to_csv(os.path.join(output_dir, "aa_composition.tsv"), sep="\t", index=False)

        plt.figure(figsize=(12, 6))
        sns.barplot(data=df_aa, x="AminoAcid", y="Percentage", hue="AminoAcid",
                    palette="magma", legend=False)
        plt.title(f"Amino Acid Composition (Total: {total_aa:,} residues)")
        plt.ylabel("Percentage (%)")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "aa_composition_plot.png"), dpi=150, bbox_inches="tight")
        plt.close()
        print("Generated AA composition plot.")
    else:
        print("No amino acid data generated.")


def main():
    parser = argparse.ArgumentParser(description="GC content & AA composition analysis")
    parser.add_argument("input_fasta", help="Input CDS FASTA")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--genetic_code", type=int, default=2)
    args = parser.parse_args()
    analyze_composition(args.input_fasta, args.output_dir, args.genetic_code)


if __name__ == "__main__":
    main()
