#!/usr/bin/env python3
"""
run_gbdraw.py – Generate BOTH circular and linear genome diagrams using gbdraw CLI.

Called by Snakemake via `script:` directive.
Searches the source annotator's folder for GenBank files and produces
publication-quality maps strictly in PDF and PNG formats for both modes.
"""

import os
import sys
import glob
import subprocess

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
src_dir     = snakemake.params.src_dir
out_dir     = snakemake.params.out_dir
circ_prefix = snakemake.params.circ_prefix
lin_prefix  = snakemake.params.lin_prefix
log_file    = str(snakemake.log[0])

# Paksa output hanya PDF dan PNG
formats     = ["pdf", "png"]

os.makedirs(out_dir, exist_ok=True)
os.makedirs(os.path.dirname(log_file), exist_ok=True)

# Redirect stdout/stderr to log
log_fh = open(log_file, "w")
sys.stdout = log_fh
sys.stderr = log_fh

print("Starting gbdraw for BOTH Circular and Linear maps...")

# ---------------------------------------------------------------------------
# Find GenBank file in annotator's output folder
# ---------------------------------------------------------------------------
gb_files = sorted(
    glob.glob(os.path.join(src_dir, "*.gb"))
    + glob.glob(os.path.join(src_dir, "*.gbk"))
)
# Filter out empty files
gb_files = [f for f in gb_files if os.path.getsize(f) > 0]

if not gb_files:
    print(f"ERROR: No valid GenBank files found in {src_dir}!")
    log_fh.close()
    sys.exit(1) # Keluar dengan error merah (bukan 0)

gb_file = gb_files[0]
print(f"Found GenBank file: {gb_file}\n")

# ---------------------------------------------------------------------------
# Build CLI Commands (Tasks for Circular and Linear)
# ---------------------------------------------------------------------------
tasks = []

# 1. Setup perintah CIRCULAR
circ_cmd = [
    "gbdraw", "circular",
    "--gbk", gb_file,
    "--separate_strands",
    "-k", "CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin",
    "--block_stroke_width", "1",
    "--block_stroke_color", "black",
    "--axis_stroke_width", "3",
    "--line_stroke_width", "2",
    "--suppress_gc",
    "--suppress_skew",
    "-p", "default",
    "--track_type", "tuckin",
    "--show_labels",
    "--allow_inner_labels",
    "--outer_label_x_radius_offset", "0.90",
    "--outer_label_y_radius_offset", "0.90",
    "--inner_label_x_radius_offset", "0.975",
    "--inner_label_y_radius_offset", "0.975",
    "--definition_font_size", "28",
    "--legend", "upper_left"
]

# Tambahkan file TSV sirkular jika ada di root directory
if os.path.exists("2025-09-19_chloroplast.tsv"):
    circ_cmd.extend(["-t", "2025-09-19_chloroplast.tsv"])
if os.path.exists("qualifier_priority.tsv"):
    circ_cmd.extend(["--qualifier_priority", "qualifier_priority.tsv"])

tasks.append({"name": "CIRCULAR", "prefix": circ_prefix, "cmd": circ_cmd})


# 2. Setup perintah LINEAR
lin_cmd = [
    "gbdraw", "linear",
    "--gbk", gb_file,
    "--show_labels",
    "--separate_strands",
    "--legend", "left",
    "--block_stroke_width", "2",
    "--axis_stroke_width", "5",
    "--definition_font_size", "24"
]

# Tambahkan file TSV linear jika ada di root directory
if os.path.exists("cds_white.tsv"):
    lin_cmd.extend(["-d", "cds_white.tsv"])
if os.path.exists("lambda_specific_table.tsv"):
    lin_cmd.extend(["-t", "lambda_specific_table.tsv"])

tasks.append({"name": "LINEAR", "prefix": lin_prefix, "cmd": lin_cmd})

# ---------------------------------------------------------------------------
# Execute gbdraw for each task and each format
# ---------------------------------------------------------------------------
try:
    for task in tasks:
        print(f"--- Generating {task['name']} maps ---")
        for fmt in formats:
            final_cmd = task["cmd"] + ["-o", task["prefix"], "-f", fmt]
            print(f"Executing: {' '.join(final_cmd)}")
            
            # subprocess.run akan melempar error jika gbdraw gagal
            subprocess.run(final_cmd, check=True)
            print(f"Successfully created {fmt.upper()} map for {task['name']}.\n")

except subprocess.CalledProcessError as e:
    print(f"\nERROR: gbdraw command failed with exit code {e.returncode}")
    log_fh.close()
    sys.exit(1) # Lapor ke Snakemake bahwa proses gagal

print("All diagrams (Circular & Linear) generated successfully.")
log_fh.close()