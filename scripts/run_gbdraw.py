#!/usr/bin/env python3
"""
run_gbdraw.py – Generate comprehensive circular genome diagrams using gbdraw Python API.

Called by Snakemake via `script:` directive.
Searches the source annotator's folder for GenBank files and produces
publication-quality circular genome maps with GC content, labels, and
multiple output formats (SVG, PNG, PDF).
"""

import os
import sys
import glob

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
src_dir      = snakemake.params.src_dir
out_prefix   = snakemake.params.out_prefix
out_dir      = snakemake.params.out_dir
formats      = snakemake.params.formats
extra_config = snakemake.params.extra_config
source_tool  = snakemake.wildcards.source_tool
sample       = snakemake.wildcards.sample
log_file     = str(snakemake.log[0])

os.makedirs(out_dir, exist_ok=True)
os.makedirs(os.path.dirname(log_file), exist_ok=True)

# Redirect stdout/stderr to log
log_fh = open(log_file, "w")
sys.stdout = log_fh
sys.stderr = log_fh

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
    print(f"No GenBank files found in {src_dir}, skipping gbdraw for {source_tool}")
    # Touch SVG output so Snakemake is satisfied
    open(snakemake.output.svg, "w").close()
    log_fh.close()
    sys.exit(0)

gb_file = gb_files[0]
print(f"Found GenBank file: {gb_file}")

# ---------------------------------------------------------------------------
# Load and validate the GenBank record
# ---------------------------------------------------------------------------
try:
    record = next(SeqIO.parse(gb_file, "genbank"))
except StopIteration:
    print(f"Could not parse GenBank file: {gb_file}")
    open(snakemake.output.svg, "w").close()
    log_fh.close()
    sys.exit(0)

# Check for gene annotations
feature_types = {f.type for f in record.features}
annotation_types = {"gene", "CDS", "tRNA", "rRNA", "mRNA"}
if not feature_types & annotation_types:
    print(f"GenBank file has no gene annotations (features: {feature_types}), skipping")
    open(snakemake.output.svg, "w").close()
    log_fh.close()
    sys.exit(0)

print(f"Record: {record.id}, length: {len(record)} bp, "
      f"features: {len(record.features)}")

# ---------------------------------------------------------------------------
# Generate diagram using gbdraw Python API
# ---------------------------------------------------------------------------
try:
    from gbdraw.api import assemble_circular_diagram_from_record
    from gbdraw.api.render import save_figure

    # Feature types to display
    selected_features = [
        "CDS", "rRNA", "tRNA", "tmRNA", "ncRNA",
        "misc_RNA", "repeat_region",
    ]

    # Config overrides for comprehensive output
    config_overrides = {
        "strandedness": True,     # Separate + and - strand features
        "show_labels": True,      # Show gene labels
    }
    # Merge any extra config from Snakemake params
    if extra_config and isinstance(extra_config, dict):
        config_overrides.update(extra_config)

    print(f"Assembling circular diagram with config: {config_overrides}")

    # Create the diagram
    canvas = assemble_circular_diagram_from_record(
        record,
        selected_features_set=selected_features,
        output_prefix=out_prefix,
        legend="right",
        config_overrides=config_overrides,
    )

    # Save in all requested formats
    print(f"Saving in formats: {formats}")
    save_figure(canvas, formats)

    print(f"Diagram saved successfully: {out_prefix}.*")

except ImportError as e:
    print(f"gbdraw Python API not available ({e}), falling back to CLI")

    # Fallback: use CLI
    import subprocess
    for fmt in formats:
        cmd = [
            "gbdraw", "circular",
            "--gbk", gb_file,
            "-o", out_prefix,
            "-f", fmt,
            "--separate_strands",
        ]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=False)

except Exception as e:
    print(f"Error generating diagram: {e}")
    import traceback
    traceback.print_exc()

# ---------------------------------------------------------------------------
# Ensure SVG output exists (Snakemake requirement)
# ---------------------------------------------------------------------------
svg_path = snakemake.output.svg
if not os.path.exists(svg_path):
    open(svg_path, "w").close()

print("Done.")
log_fh.close()
