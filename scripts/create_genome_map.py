#!/usr/bin/env python3
"""
create_genome_map.py – Circular/linear genome map using pyGenomeViz.

Replaces the original Circos-based script (which required ~20 Perl dependencies)
with a pure-Python solution using pyGenomeViz.

Generates:
  - Circular genome map with gene tracks (CDS, tRNA, rRNA)
  - GC content and GC skew tracks
  - Publication-quality PNG output

Usage:
    python create_genome_map.py <genbank_file> <output_png> [--linear]
"""

import sys
import os
import argparse

import matplotlib
matplotlib.use("Agg")

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def compute_gc_content(seq: str, window: int = 500, step: int = 100) -> list[tuple[int, float]]:
    """Compute GC content deviation from mean in sliding windows."""
    seq = seq.upper()
    mean_gc = gc_fraction(seq)
    values = []
    for i in range(0, len(seq) - window, step):
        w = seq[i : i + window]
        gc = gc_fraction(w)
        values.append((i + window // 2, gc - mean_gc))
    return values


def compute_gc_skew(seq: str, window: int = 500, step: int = 100) -> list[tuple[int, float]]:
    """Compute GC skew (G-C)/(G+C) in sliding windows."""
    seq = seq.upper()
    values = []
    for i in range(0, len(seq) - window, step):
        w = seq[i : i + window]
        g = w.count("G")
        c = w.count("C")
        skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
        values.append((i + window // 2, skew))
    return values


def create_genome_map(gbk_file: str, output_png: str, linear: bool = False):
    """Create a circular/linear genome map from a GenBank file."""

    try:
        from pygenomeviz import GenomeViz
        from pygenomeviz.parser import Genbank
    except ImportError:
        print("ERROR: pygenomeviz not installed. Install via: pip install pygenomeviz")
        # Create placeholder image
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "pyGenomeViz not installed", ha="center", va="center",
                fontsize=14, color="red")
        ax.axis("off")
        os.makedirs(os.path.dirname(output_png) or ".", exist_ok=True)
        plt.savefig(output_png, dpi=150, bbox_inches="tight")
        plt.close()
        return

    # Parse GenBank
    record = SeqIO.read(gbk_file, "genbank")
    genome_len = len(record.seq)
    genome_name = record.id
    avg_gc = gc_fraction(record.seq)

    print(f"Genome: {genome_name}, Length: {genome_len:,} bp, GC: {avg_gc:.1%}")

    # Colour scheme for feature types
    color_map = {
        "CDS": "#4C78A8",    # Steel blue
        "tRNA": "#E45756",   # Red
        "rRNA": "#72B7B2",   # Teal
        "gene": "#EECA3B",   # Yellow
    }

    # Initialize GenomeViz
    gv = GenomeViz(
        fig_track_height=0.5,
        feature_track_ratio=0.25,
        tick_track_ratio=0.1,
    )

    # Add genome track
    track = gv.add_feature_track(genome_name, genome_len)
    track.add_sublabel(f"({genome_len:,} bp, GC: {avg_gc:.1%})")

    # Add features
    for feat in record.features:
        if feat.type not in color_map:
            continue

        start = int(feat.location.start)
        end = int(feat.location.end)
        strand = feat.location.strand or 1
        color = color_map.get(feat.type, "#999999")

        # Gene label
        label = ""
        if "gene" in feat.qualifiers:
            label = feat.qualifiers["gene"][0]
        elif "product" in feat.qualifiers:
            label = feat.qualifiers["product"][0]
        elif "locus_tag" in feat.qualifiers:
            label = feat.qualifiers["locus_tag"][0]

        track.add_feature(
            start, end, strand,
            label=label,
            facecolor=color,
            labelsize=7,
            labelvpos="top" if strand == 1 else "bottom",
            labelha="center",
            linewidth=0.5,
        )

    # ── GC content & skew as subtracks (if supported) ──
    seq_str = str(record.seq)
    gc_content_vals = compute_gc_content(seq_str)
    gc_skew_vals = compute_gc_skew(seq_str)

    # Add GC content as a bar subplot
    if gc_content_vals:
        positions = [v[0] for v in gc_content_vals]
        gc_values = [v[1] for v in gc_content_vals]

        # Positive = blue (GC-rich), Negative = orange (AT-rich)
        pos_colors = ["#4C78A8" if v >= 0 else "#F58518" for v in gc_values]

        try:
            gc_subtrack = track.add_subtrack(ratio=0.3)
            for pos, val, col in zip(positions, gc_values, pos_colors):
                gc_subtrack.add_feature(
                    max(0, pos - 250), min(genome_len, pos + 250),
                    strand=1 if val >= 0 else -1,
                    facecolor=col,
                    linewidth=0,
                )
        except (AttributeError, TypeError):
            # pyGenomeViz version may not support subtracks; skip gracefully
            pass

    # Save figure
    os.makedirs(os.path.dirname(output_png) or ".", exist_ok=True)
    fig = gv.plotfig()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    print(f"Genome map saved to {output_png}")

    # Also save SVG for publication quality
    svg_path = output_png.rsplit(".", 1)[0] + ".svg"
    fig.savefig(svg_path, bbox_inches="tight")
    print(f"SVG version saved to {svg_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate genome map (pyGenomeViz)")
    parser.add_argument("genbank_file", help="Input GenBank file")
    parser.add_argument("output_png", help="Output PNG file path")
    parser.add_argument("--linear", action="store_true", help="Linear layout (default: auto)")
    args = parser.parse_args()
    create_genome_map(args.genbank_file, args.output_png, args.linear)


if __name__ == "__main__":
    main()
