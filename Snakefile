# =============================================================================
# Organelle Annotation Pipeline – Snakemake Workflow
# =============================================================================
# Tools:  Chloe.jl | PGA
#         MFannot (Docker) | fpma | MITOS2 (Docker) | MitoZ (Docker)
#         tRNAscan-SE | Aragorn | Liftoff
# Viz:    OGDraw (Docker) | pyGenomeViz (genome map)
# QC:     BUSCO + custom gene-completeness summary
# Report: Aggregated indexed HTML
# Downstream: RSCU | Codon Usage | Ka/Ks (MAFFT + KaKs_Calculator) |
#             Phylogeny (IQ-TREE) | GC/AA Composition | Genome Map |
#             Synteny (MUMmer4) | NCBI Reference Fetch
# =============================================================================
#
# Directory layout (per-sample):
#   results/{sample}/{tool}/          – annotation outputs
#   results/{sample}/logs/{tool}.log  – per-tool logs
#   results/{sample}/qc/             – QC results
#   results/{sample}/downstream/     – downstream analyses
# =============================================================================

import os
import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Load samples
# ---------------------------------------------------------------------------
samples_df = pd.read_csv(config["samples"], sep="\t", comment="#", dtype=str)
samples_df = samples_df.fillna("")
samples_df = samples_df.set_index("sample", drop=False)

SAMPLES = samples_df["sample"].tolist()
OUTDIR = config["outdir"]

# ---------------------------------------------------------------------------
# Helper: which tools to run per sample
# ---------------------------------------------------------------------------
PLASTID_TOOLS = ["chloe", "pga"]
MITO_TOOLS = ["mfannot", "fpma", "mitos", "mitoz"]
BOTH_TOOLS = ["trnascan", "aragorn", "liftoff"]
ALL_TOOLS = PLASTID_TOOLS + MITO_TOOLS + BOTH_TOOLS

def tools_for_sample(sample):
    """Return list of tools applicable to a sample based on config mode.
    Liftoff is automatically excluded if the sample has no reference."""
    organelle = samples_df.loc[sample, "organelle"]
    mode = config["mode"]

    if mode == "all":
        if organelle == "plastid":
            tools = PLASTID_TOOLS + BOTH_TOOLS
        elif organelle == "mito":
            tools = MITO_TOOLS + BOTH_TOOLS
        else:
            tools = ALL_TOOLS
    elif mode == "plastid":
        tools = PLASTID_TOOLS + BOTH_TOOLS
    elif mode == "mito":
        tools = MITO_TOOLS + BOTH_TOOLS
    elif mode == "select":
        selected = config.get("tools_select", ALL_TOOLS)
        if organelle == "plastid":
            tools = [t for t in selected if t in PLASTID_TOOLS + BOTH_TOOLS]
        elif organelle == "mito":
            tools = [t for t in selected if t in MITO_TOOLS + BOTH_TOOLS]
        else:
            tools = list(selected)
    else:
        tools = ALL_TOOLS

    # Skip liftoff when reference columns are empty
    ref_fasta = get_reference_fasta(sample)
    ref_gff = get_reference_gff(sample)
    if not ref_fasta or not ref_gff:
        tools = [t for t in tools if t != "liftoff"]

    return tools

def get_fasta(sample):
    return samples_df.loc[sample, "fasta"]

def get_genetic_code(sample):
    return samples_df.loc[sample, "genetic_code"]

def get_organelle(sample):
    return samples_df.loc[sample, "organelle"]

def get_reference_fasta(sample):
    """Get per-sample reference FASTA path (empty string if not set)."""
    return samples_df.loc[sample, "reference_fasta"].strip()

def get_reference_gff(sample):
    """Get per-sample reference GFF path (empty string if not set)."""
    return samples_df.loc[sample, "reference_gff"].strip()

# ---------------------------------------------------------------------------
# Collect all expected outputs
# ---------------------------------------------------------------------------
def all_outputs():
    outputs = []
    for s in SAMPLES:
        for tool in tools_for_sample(s):
            outputs.append(f"{OUTDIR}/{s}/{tool}/{s}.done")
        # OGDraw maps — one per source tool that produces a GenBank file
        for src in ogdraw_source_tools(s):
            outputs.append(f"{OUTDIR}/{s}/ogdraw/{src}/{s}.done")
        # QC outputs
        if config["qc"]["enabled"]:
            outputs.append(f"{OUTDIR}/{s}/qc/qc_summary.tsv")
            if get_organelle(s) in ("plastid", "mito"):
                outputs.append(f"{OUTDIR}/{s}/qc/busco/short_summary.txt")
        # Downstream analysis outputs
        if config.get("downstream", {}).get("enabled", False):
            outputs.append(f"{OUTDIR}/{s}/downstream/downstream_report.html")
    # Final report
    outputs.append(f"{OUTDIR}/report/index.html")
    return outputs

rule all:
    input:
        all_outputs()

# ---------------------------------------------------------------------------
# Include per-tool rule files
# ---------------------------------------------------------------------------
include: "rules/chloe.smk"
include: "rules/pga.smk"
include: "rules/mfannot.smk"
include: "rules/fpma.smk"
include: "rules/mitos.smk"
include: "rules/mitoz.smk"
include: "rules/trnascan.smk"
include: "rules/aragorn.smk"
include: "rules/liftoff.smk"
include: "rules/ogdraw.smk"
include: "rules/qc.smk"
include: "rules/report.smk"
include: "rules/downstream.smk"
