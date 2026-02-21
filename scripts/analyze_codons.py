#!/usr/bin/env python3
"""
analyze_codons.py – Start / stop codon analysis for CDS features.

Reads a CDS FASTA file and reports the start and stop codon usage frequency.

Usage:
    python analyze_codons.py <input_fasta> <output_file> [--genetic_code 2]
"""

import sys
import os
import argparse
from collections import Counter

from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gene_utils import normalise_gene_name


def analyze_codons(input_fasta: str, output_file: str, genetic_code: int = 2):
    start_codons: Counter = Counter()
    stop_codons: Counter = Counter()
    gene_details: list[dict] = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq).upper().replace("-", "")
        if len(seq) < 6:
            continue

        gene_name = normalise_gene_name(record.id)
        start = seq[:3]
        stop = seq[-3:]

        start_codons[start] += 1
        stop_codons[stop] += 1
        gene_details.append({
            "gene": gene_name,
            "start_codon": start,
            "stop_codon": stop,
            "length_bp": len(seq),
        })

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)

    with open(output_file, "w") as f:
        f.write(f"=== Start/Stop Codon Analysis (genetic code: {genetic_code}) ===\n\n")
        f.write(f"Total CDS analysed: {len(gene_details)}\n\n")

        f.write("Start Codon Frequencies:\n")
        for codon, count in start_codons.most_common():
            pct = count / len(gene_details) * 100 if gene_details else 0
            f.write(f"  {codon}: {count:>4} ({pct:5.1f}%)\n")

        f.write("\nStop Codon Frequencies:\n")
        for codon, count in stop_codons.most_common():
            pct = count / len(gene_details) * 100 if gene_details else 0
            f.write(f"  {codon}: {count:>4} ({pct:5.1f}%)\n")

        f.write("\nPer-gene Details:\n")
        f.write(f"{'Gene':<20} {'Start':<8} {'Stop':<8} {'Length (bp)':<12}\n")
        f.write("-" * 50 + "\n")
        for g in sorted(gene_details, key=lambda x: x["gene"]):
            f.write(f"{g['gene']:<20} {g['start_codon']:<8} {g['stop_codon']:<8} {g['length_bp']:<12}\n")

    print(f"Codon analysis written to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Analyse start/stop codons in CDS FASTA")
    parser.add_argument("input_fasta", help="Input CDS FASTA file")
    parser.add_argument("output_file", help="Output text file")
    parser.add_argument("--genetic_code", type=int, default=2)
    args = parser.parse_args()
    analyze_codons(args.input_fasta, args.output_file, args.genetic_code)


if __name__ == "__main__":
    main()
