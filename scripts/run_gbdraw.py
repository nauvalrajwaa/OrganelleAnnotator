#!/usr/bin/env python3
"""
run_gbdraw.py – Generate comprehensive genome diagrams using gbdraw CLI.

Called by Snakemake via `script:` directive.
Searches the source annotator's folder for GenBank files and produces
publication-quality maps strictly in PDF and PNG formats.
"""

import os
import sys
import glob
import subprocess

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
src_dir      = snakemake.params.src_dir
out_prefix   = snakemake.params.out_prefix
out_dir      = snakemake.params.out_dir
source_tool  = snakemake.wildcards.source_tool
sample       = snakemake.wildcards.sample
log_file     = str(snakemake.log[0])

# Paksa output hanya PDF dan PNG (sesuai permintaan)
formats      = ["pdf", "png"]

os.makedirs(out_dir, exist_ok=True)
os.makedirs(os.path.dirname(log_file), exist_ok=True)

# Redirect stdout/stderr to log
log_fh = open(log_file, "w")
sys.stdout = log_fh
sys.stderr = log_fh

print(f"Starting gbdraw strictly for formats: {formats}")

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
print(f"Found GenBank file: {gb_file}")

# ---------------------------------------------------------------------------
# Build CLI Command based on user's advanced parameters
# ---------------------------------------------------------------------------
# Anda bisa mengubah ini menjadi "linear" dari file .smk melalui params
draw_mode = snakemake.params.get("draw_mode", "circular")

if draw_mode == "circular":
    cmd_base = [
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
    
    # Tambahkan file TSV jika ada di folder kerja (mencegah error jika file tidak ada)
    if os.path.exists("2025-09-19_chloroplast.tsv"):
        cmd_base.extend(["-t", "2025-09-19_chloroplast.tsv"])
    if os.path.exists("qualifier_priority.tsv"):
        cmd_base.extend(["--qualifier_priority", "qualifier_priority.tsv"])

else:
    # Mode Linear
    cmd_base = [
        "gbdraw", "linear",
        "--gbk", gb_file,
        "--show_labels",
        "--separate_strands",
        "--legend", "left",
        "--block_stroke_width", "2",
        "--axis_stroke_width", "5",
        "--definition_font_size", "24"
    ]
    
    if os.path.exists("cds_white.tsv"):
        cmd_base.extend(["-d", "cds_white.tsv"])
    if os.path.exists("lambda_specific_table.tsv"):
        cmd_base.extend(["-t", "lambda_specific_table.tsv"])

# ---------------------------------------------------------------------------
# Execute gbdraw for each format
# ---------------------------------------------------------------------------
try:
    for fmt in formats:
        cmd = cmd_base + ["-o", out_prefix, "-f", fmt]
        print(f"\nExecuting: {' '.join(cmd)}")
        
        # subprocess.run akan melempar error jika gbdraw gagal
        subprocess.run(cmd, check=True)
        print(f"Successfully created {fmt.upper()} map.")

except subprocess.CalledProcessError as e:
    print(f"\nERROR: gbdraw command failed with exit code {e.returncode}")
    log_fh.close()
    sys.exit(1) # Lapor ke Snakemake bahwa proses gagal

print("\nDone.")
log_fh.close()