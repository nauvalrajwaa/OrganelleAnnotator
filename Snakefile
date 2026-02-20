# =============================================================================
# Organelle Annotation Pipeline – Snakemake Workflow
# =============================================================================
# Tools:  Chloe.jl | PGA | MFannot (Docker) | fpma | MITOS2 (Docker)
#         MitoZ (Docker) | tRNAscan-SE | Aragorn
# QC:     BUSCO + custom gene-completeness summary
# Report: Aggregated indexed HTML
# =============================================================================

import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Load samples
# ---------------------------------------------------------------------------
samples_df = pd.read_csv(config["samples"], sep="\t", comment="#", dtype=str)
samples_df = samples_df.set_index("sample", drop=False)

SAMPLES = samples_df["sample"].tolist()
OUTDIR = config["outdir"]

# ---------------------------------------------------------------------------
# Helper: which tools to run per sample
# ---------------------------------------------------------------------------
PLASTID_TOOLS = ["chloe", "pga"]
MITO_TOOLS = ["mfannot", "fpma", "mitos", "mitoz"]
BOTH_TOOLS = ["trnascan", "aragorn"]
ALL_TOOLS = PLASTID_TOOLS + MITO_TOOLS + BOTH_TOOLS

def tools_for_sample(sample):
    """Return list of tools applicable to a sample based on config mode."""
    organelle = samples_df.loc[sample, "organelle"]
    mode = config["mode"]

    if mode == "all":
        if organelle == "plastid":
            return PLASTID_TOOLS + BOTH_TOOLS
        elif organelle == "mito":
            return MITO_TOOLS + BOTH_TOOLS
        else:
            return ALL_TOOLS
    elif mode == "plastid":
        return PLASTID_TOOLS + BOTH_TOOLS
    elif mode == "mito":
        return MITO_TOOLS + BOTH_TOOLS
    elif mode == "select":
        selected = config.get("tools_select", ALL_TOOLS)
        if organelle == "plastid":
            return [t for t in selected if t in PLASTID_TOOLS + BOTH_TOOLS]
        elif organelle == "mito":
            return [t for t in selected if t in MITO_TOOLS + BOTH_TOOLS]
        else:
            return selected
    return ALL_TOOLS

def get_fasta(sample):
    return samples_df.loc[sample, "fasta"]

def get_genetic_code(sample):
    return samples_df.loc[sample, "genetic_code"]

def get_organelle(sample):
    return samples_df.loc[sample, "organelle"]

# ---------------------------------------------------------------------------
# Collect all expected outputs
# ---------------------------------------------------------------------------
def all_outputs():
    outputs = []
    for s in SAMPLES:
        for tool in tools_for_sample(s):
            outputs.append(f"{OUTDIR}/{tool}/{s}/{s}.done")
        # QC outputs
        if config["qc"]["enabled"]:
            outputs.append(f"{OUTDIR}/qc/summary/{s}.qc_summary.tsv")
            if get_organelle(s) in ("plastid", "mito"):
                outputs.append(f"{OUTDIR}/qc/busco/{s}/short_summary.txt")
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
include: "rules/qc.smk"
include: "rules/report.smk"
