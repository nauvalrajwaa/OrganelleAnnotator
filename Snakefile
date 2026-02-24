# =============================================================================
# Organelle Annotation Pipeline – Snakemake Workflow
# =============================================================================
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
samples_df = pd.read_csv(
    config["samples"], sep="\t", comment="#", dtype=str,
    skipinitialspace=True,
)
samples_df.columns = samples_df.columns.str.strip()
samples_df = samples_df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)
samples_df = samples_df.fillna("")
samples_df = samples_df.set_index("sample", drop=False)

SAMPLES = samples_df["sample"].tolist()
OUTDIR  = "results"

# -- Debug: print parsed samples for troubleshooting -------------------------
print(f"\n  📋 Loaded {len(SAMPLES)} sample(s) from {config['samples']}")
print(f"     Columns: {list(samples_df.columns)}")
for _s in SAMPLES:
    print(f"     ✔ {_s}  →  fasta={samples_df.loc[_s, 'fasta']}  "
          f"organelle={samples_df.loc[_s, 'organelle']}  "
          f"genetic_code={samples_df.loc[_s, 'genetic_code']}")
print(f"  🔧 Mode: {config['mode']}  |  tools_select: {config.get('tools_select', 'ALL')}")
print()

# -- Validation --------------------------------------------------------------
if len(samples_df.columns) < 3:
    raise ValueError(
        "\n\nERROR: samples.tsv has fewer than 3 columns.\n"
        "Make sure the file is TAB-separated (not spaces).\n"
    )
for _s in SAMPLES:
    if "/" in _s or "\\" in _s:
        raise ValueError(
            f"\n\nERROR: Invalid sample name '{_s}'.\n"
            "Sample names must be short IDs (e.g. 'sample1'), not file paths.\n"
            "Ensure samples.tsv is TAB-separated: sample<TAB>fasta<TAB>organelle<TAB>...\n"
        )

# ---------------------------------------------------------------------------
# Tool sets
# ---------------------------------------------------------------------------
PLASTID_TOOLS = ["chloe", "pga"]
MITO_TOOLS    = ["mfannot", "fpma", "mitos", "mitoz"]
BOTH_TOOLS    = ["trnascan", "aragorn", "liftoff"]
ALL_TOOLS     = PLASTID_TOOLS + MITO_TOOLS + BOTH_TOOLS

# tools_select is the AUTHORITATIVE list — only tools listed here will run.
# In any mode, the final tool list is intersected with tools_select.
SELECTED_TOOLS = config.get("tools_select", ALL_TOOLS)

def tools_for_sample(sample):
    """Return list of tools applicable to a sample.

    Steps:
      1. Determine candidate tools based on mode + organelle type.
      2. Intersect with tools_select so only listed tools are processed.
      3. Drop liftoff if the sample has no reference columns.
    """
    organelle = samples_df.loc[sample, "organelle"]
    mode = config["mode"]

    # Step 1 — candidate tools by mode
    if mode == "plastid":
        candidates = PLASTID_TOOLS + BOTH_TOOLS
    elif mode == "mito":
        candidates = MITO_TOOLS + BOTH_TOOLS
    else:  # "all", "select", or anything else
        if organelle == "plastid":
            candidates = PLASTID_TOOLS + BOTH_TOOLS
        elif organelle == "mito":
            candidates = MITO_TOOLS + BOTH_TOOLS
        else:
            candidates = ALL_TOOLS

    # Step 2 — intersect with tools_select (only run what's listed)
    tools = [t for t in candidates if t in SELECTED_TOOLS]

    # Step 3 — skip liftoff when reference columns are empty
    ref_fa  = samples_df.loc[sample, "reference_fasta"].strip()
    ref_gff = samples_df.loc[sample, "reference_gff"].strip()
    if not ref_fa or not ref_gff:
        tools = [t for t in tools if t != "liftoff"]

    return tools

# ---------------------------------------------------------------------------
# gbdraw helpers
# ---------------------------------------------------------------------------
# Annotation tools that MAY produce GenBank (.gb/.gbk) files.
# gbdraw will search each tool's output folder at runtime for *.gb / *.gbk.
GB_PRODUCING_TOOLS = {
    "plastid": ["chloe", "pga", "liftoff"],
    "mito":    ["mitoz", "liftoff"],
}

def gbdraw_source_tools(sample):
    """Return list of tools that produce GenBank files for this sample."""
    organelle = samples_df.loc[sample, "organelle"]
    possible  = GB_PRODUCING_TOOLS.get(organelle, GB_PRODUCING_TOOLS.get("plastid", []))
    active    = tools_for_sample(sample)
    return [t for t in possible if t in active]

# ---------------------------------------------------------------------------
# Collect all expected outputs
# ---------------------------------------------------------------------------
def all_outputs():
    outputs = []
    for s in SAMPLES:
        for tool in tools_for_sample(s):
            outputs.append(f"{OUTDIR}/{s}/{tool}/{s}.done")
        for src in gbdraw_source_tools(s):
            outputs.append(f"{OUTDIR}/{s}/gbdraw/{src}/{s}.done")
        if config["qc"]["enabled"]:
            outputs.append(f"{OUTDIR}/{s}/qc/qc_summary.tsv")
            if samples_df.loc[s, "organelle"] in ("plastid", "mito"):
                outputs.append(f"{OUTDIR}/{s}/qc/busco/short_summary.txt")
        if config.get("downstream", {}).get("enabled", False):
            outputs.append(f"{OUTDIR}/{s}/downstream/downstream_report.html")
        # Per-sample report
        outputs.append(f"{OUTDIR}/{s}/report/index.html")
    return outputs

rule all:
    input:
        all_outputs(),
        expand(OUTDIR + "/{sample}/chloe/{sample}.gb", sample=samples_df.index)

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
include: "rules/gbdraw.smk"
include: "rules/qc.smk"
include: "rules/report.smk"
include: "rules/downstream.smk"
