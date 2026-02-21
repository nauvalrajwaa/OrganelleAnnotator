#!/usr/bin/env python3
"""
run_synteny_analysis.py – Synteny analysis using MUMmer4 (nucmer).

Aligns sample vs. reference genome and generates a Bezier ribbon plot
showing collinear and inverted synteny blocks.

Features:
  - FASTA sanitisation (fixes MUMmer error 400 from complex headers)
  - Robust MUMmer coord parsing
  - Fail-safe outputs (prevents Snakemake MissingOutputException)
  - Publication-quality Bezier ribbon visualisation

Usage:
    python run_synteny_analysis.py <sample.fasta> <ref.fasta> <output.png> <stats.tsv>
"""

import sys
import os
import subprocess
import tempfile
import shutil

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from Bio import SeqIO


# ── Utilities ────────────────────────────────────────────────────────────────

def sanitize_fasta(input_fasta: str, output_fasta: str, generic_name: str) -> int:
    """Rewrite FASTA with clean headers for MUMmer compatibility.

    Returns total genome length.
    """
    total_len = 0
    clean_records = []
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        return 0

    for i, record in enumerate(records):
        record.id = f"{generic_name}_{i + 1}"
        record.description = ""
        clean_records.append(record)
        total_len += len(record.seq)

    SeqIO.write(clean_records, output_fasta, "fasta")
    return total_len


def create_failure_outputs(out_plot: str, out_stats: str, message: str = "Analysis Failed"):
    """Create placeholder outputs so Snakemake does not crash."""
    os.makedirs(os.path.dirname(out_plot) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(out_stats) or ".", exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 2))
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=14, color="red")
    ax.axis("off")
    plt.savefig(out_plot, dpi=100, bbox_inches="tight")
    plt.close()

    pd.DataFrame({"status": ["failed"], "error": [message]}).to_csv(
        out_stats, sep="\t", index=False
    )


# ── MUMmer workflow ──────────────────────────────────────────────────────────

def run_mummer_workflow(sample_fasta: str, ref_fasta: str, output_prefix: str) -> str | None:
    """Run nucmer -> delta-filter -> show-coords pipeline."""
    nucmer = shutil.which("nucmer")
    delta_filter = shutil.which("delta-filter")
    show_coords = shutil.which("show-coords")

    if not (nucmer and delta_filter and show_coords):
        print("ERROR: MUMmer4 tools (nucmer, delta-filter, show-coords) not found in PATH.")
        return None

    tmpdir = tempfile.mkdtemp(prefix="synteny_")

    try:
        tmp_s = os.path.join(tmpdir, "sample.fasta")
        tmp_r = os.path.join(tmpdir, "ref.fasta")
        sanitize_fasta(sample_fasta, tmp_s, "Query")
        sanitize_fasta(ref_fasta, tmp_r, "Ref")

        prefix = os.path.join(tmpdir, "align")
        delta = f"{prefix}.delta"
        filtered = f"{prefix}.filter.delta"
        coords = f"{prefix}.coords"

        # 1. nucmer (--maxmatch for organelles which have repeats/IRs)
        subprocess.run(
            [nucmer, "--maxmatch", f"--prefix={prefix}", tmp_r, tmp_s],
            check=True, capture_output=True,
        )

        # 2. delta-filter (-m many-to-many, -i 85% identity)
        with open(filtered, "w") as f:
            subprocess.run([delta_filter, "-m", "-i", "85", delta], stdout=f, check=True)

        # 3. show-coords (-r ref-sorted, -c coverage, -l length, -T tab)
        with open(coords, "w") as f:
            subprocess.run([show_coords, "-r", "-c", "-l", "-T", filtered], stdout=f, check=True)

        return coords

    except subprocess.CalledProcessError as exc:
        print(f"MUMmer process error: {exc}")
        return None
    except Exception as exc:
        print(f"Error: {exc}")
        return None


def parse_mummer_coords(coords_file: str) -> list[dict]:
    """Parse show-coords tab-delimited output into synteny blocks."""
    blocks = []
    if not os.path.exists(coords_file) or os.path.getsize(coords_file) == 0:
        return blocks

    try:
        # show-coords -T has 4 header lines, then tab-delimited data
        with open(coords_file) as f:
            lines = f.readlines()

        # Find the data start (skip header lines)
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith("[S1]") or "\t" in line and line[0].isdigit():
                data_start = i + 1 if line.startswith("[") else i
                break
        if data_start == 0:
            data_start = 4  # Default skip

        for idx, line in enumerate(lines[data_start:]):
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            try:
                ref_start, ref_end = int(parts[0]), int(parts[1])
                q_start, q_end = int(parts[2]), int(parts[3])
                q_strand = "+" if q_start < q_end else "-"

                blocks.append({
                    "id": idx,
                    "sequences": [
                        {"seq_id": "ref", "start": ref_start, "end": ref_end, "strand": "+"},
                        {"seq_id": "query", "start": min(q_start, q_end),
                         "end": max(q_start, q_end), "strand": q_strand},
                    ],
                })
            except (ValueError, IndexError):
                continue

    except Exception as exc:
        print(f"Parse error: {exc}")

    return blocks


# ── Visualisation ────────────────────────────────────────────────────────────

def get_bezier_path(st_top, et_top, st_bot, et_bot, y_top, y_bot):
    """Generate a smooth Bezier ribbon path between two genomic intervals."""
    mid_y = (y_top + y_bot) / 2
    verts = [
        (st_top, y_top), (st_top, mid_y), (st_bot, mid_y), (st_bot, y_bot),
        (et_bot, y_bot), (et_bot, mid_y), (et_top, mid_y), (et_top, y_top),
        (st_top, y_top),
    ]
    codes = [
        Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    return Path(verts, codes)


def draw_genome_ruler(ax, length, y_pos, label, color="#333333"):
    """Draw a genome ruler line with tick marks."""
    ax.plot([0, length], [y_pos, y_pos], color=color, linewidth=2, zorder=5)
    ax.text(-length * 0.02, y_pos, label, va="center", ha="right",
            fontsize=12, fontweight="bold", color=color)

    interval = 50000
    if length < 50000:
        interval = 5000
    elif length < 200000:
        interval = 20000

    for i in range(0, length + 1, interval):
        ax.plot([i, i], [y_pos - 0.05, y_pos + 0.05], color=color, linewidth=1, zorder=5)
        if i % (interval * 2) == 0:
            ax.text(i, y_pos + 0.07, f"{i / 1000:.0f}k", ha="center", va="bottom",
                    fontsize=9, color=color)


def plot_synteny(blocks, sample_len, ref_len, output_plot):
    """Draw Bezier ribbon synteny plot."""
    fig, ax = plt.subplots(figsize=(15, 7))
    y_sample, y_ref = 0.8, 0.2

    for block in blocks:
        if len(block["sequences"]) < 2:
            continue
        ref = block["sequences"][0]
        query = block["sequences"][1]
        is_fwd = ref["strand"] == query["strand"]

        color = "#2E86AB" if is_fwd else "#D64045"
        alpha = 0.6 if is_fwd else 0.5

        if is_fwd:
            path = get_bezier_path(
                query["start"], query["end"], ref["start"], ref["end"],
                y_sample - 0.02, y_ref + 0.02,
            )
        else:
            path = get_bezier_path(
                query["start"], query["end"], ref["end"], ref["start"],
                y_sample - 0.02, y_ref + 0.02,
            )

        ax.add_patch(patches.PathPatch(path, facecolor=color, edgecolor="none", alpha=alpha))

    max_len = max(sample_len, ref_len)
    draw_genome_ruler(ax, sample_len, y_sample, "Sample")
    draw_genome_ruler(ax, ref_len, y_ref, "Reference")

    legend_els = [
        patches.Patch(facecolor="#2E86AB", alpha=0.6, label="Collinear"),
        patches.Patch(facecolor="#D64045", alpha=0.5, label="Inverted"),
    ]
    ax.legend(handles=legend_els, loc="upper center", bbox_to_anchor=(0.5, 1.1),
              ncol=2, frameon=False)

    ax.set_xlim(-max_len * 0.05, max_len * 1.05)
    ax.set_ylim(0, 1.2)
    ax.axis("off")
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_plot) or ".", exist_ok=True)
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 5:
        print("Usage: python run_synteny_analysis.py <sample.fasta> <ref.fasta> <plot.png> <stats.tsv>")
        sys.exit(1)

    s_fasta, r_fasta, out_plot, out_stats = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    os.makedirs(os.path.dirname(out_plot) or ".", exist_ok=True)

    try:
        print("--- Starting Synteny Analysis ---")

        coords = run_mummer_workflow(s_fasta, r_fasta, out_plot)
        blocks = parse_mummer_coords(coords) if coords else []

        if blocks:
            s_len = sum(len(r) for r in SeqIO.parse(s_fasta, "fasta"))
            r_len = sum(len(r) for r in SeqIO.parse(r_fasta, "fasta"))

            print(f"Plotting {len(blocks)} synteny blocks...")
            plot_synteny(blocks, s_len, r_len, out_plot)

            stats = {"blocks": len(blocks), "sample_len": s_len, "ref_len": r_len, "status": "success"}
            pd.DataFrame([stats]).to_csv(out_stats, sep="\t", index=False)
            print("Synteny analysis complete.")
        else:
            print("WARNING: No synteny blocks found or alignment failed.")
            create_failure_outputs(out_plot, out_stats, "No Synteny Detected")

    except Exception as exc:
        print(f"CRITICAL ERROR: {exc}")
        create_failure_outputs(out_plot, out_stats, f"Error: {str(exc)[:40]}")


if __name__ == "__main__":
    main()
