#!/usr/bin/env python3
"""
generate_report.py – Build an indexed HTML report for the Organelle Annotation Pipeline.

Called by Snakemake via `script:` directive; uses `snakemake` object for I/O.
"""

import csv
import os
import html as html_mod
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
outdir = snakemake.params.outdir
samples = snakemake.params.samples
output_html = snakemake.output.html

TOOL_LABELS = {
    "chloe": "Chloë (Chloe.jl) — Chloroplast Annotator",
    "pga": "PGA — Plastid Genome Annotator",
    "plann": "Plann — Reference-based Plastid Annotator",
    "cpgavas2": "CPGAVAS2 — Chloroplast Annotation & Visualisation",
    "mfannot": "MFannot — Mitochondrial/Plastid Annotator",
    "fpma": "fpma — Fast Plant Mitochondria Annotator",
    "mitos": "MITOS2 — Mitochondrial Genome Annotator",
    "mitoz": "MitoZ — Animal Mitochondrial Genome Annotator",
    "trnascan": "tRNAscan-SE — tRNA Gene Detection",
    "aragorn": "Aragorn — tRNA/tmRNA Detection",
    "liftoff": "Liftoff — Reference-based Annotation Lift-over",
    "ogdraw": "OGDraw — Circular Genome Map Visualisation",
}

TOOL_DESCRIPTIONS = {
    "chloe": "Julia-based chloroplast genome annotator using XGBoost models and suffix-array alignment.",
    "pga": "Perl/BLAST pipeline for rapid batch annotation of plastid genomes against GenBank references.",
    "plann": "Perl tool that transfers annotations from a reference plastid GenBank file via BLAST alignments.",
    "cpgavas2": "Docker-based chloroplast annotator using BLAST+HMMER against curated plant cp protein DB. Detects IRs and produces circular maps.",
    "mfannot": "Comprehensive mitochondrial/plastid annotator (Docker: nbeck/mfannot) using BLAST, HMMER, Exonerate, Erpin.",
    "fpma": "Rust-based fast HMM scanner for presence/absence of mitochondrial genes using HMMER3 nhmmer.",
    "mitos": "Reference-based mitochondrial genome annotator (Docker: quay.io/biocontainers/mitos) for protein-coding genes, tRNAs, and rRNAs.",
    "mitoz": "Docker-based animal mitochondrial genome annotator with circular visualisation (Docker: guanliangmeng/mitoz).",
    "trnascan": "Gold-standard tRNA detection tool using covariance models. Supports organellar/mitochondrial mode (-O).",
    "aragorn": "Lightweight tRNA and tmRNA detection using homology search. Fast and suitable for organelle genomes.",
    "liftoff": "Minimap2-based annotation lift-over from a reference organelle genome. Works for both cp and mt genomes.",
    "ogdraw": "OrganellarGenomeDRAW — generates publication-quality circular and linear genome maps from GenBank files.",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_busco_summary(path):
    """Extract key BUSCO metrics from short_summary.txt."""
    metrics = {}
    if not os.path.exists(path):
        return metrics
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "Complete BUSCOs" in line:
                metrics["complete"] = line.split()[0] if line[0].isdigit() else line
            elif "Complete and single-copy" in line:
                metrics["single_copy"] = line.split()[0] if line[0].isdigit() else line
            elif "Complete and duplicated" in line:
                metrics["duplicated"] = line.split()[0] if line[0].isdigit() else line
            elif "Fragmented" in line:
                metrics["fragmented"] = line.split()[0] if line[0].isdigit() else line
            elif "Missing" in line:
                metrics["missing"] = line.split()[0] if line[0].isdigit() else line
            elif "Total" in line and "BUSCO" in line:
                metrics["total"] = line.split()[0] if line[0].isdigit() else line
    return metrics


def read_qc_summary(path):
    """Read per-tool gene completeness TSV."""
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def list_tool_outputs(sample, tool, outdir):
    """List result files produced by a tool for a sample."""
    tool_dir = os.path.join(outdir, tool, sample)
    if not os.path.isdir(tool_dir):
        return []
    files = []
    for fn in sorted(os.listdir(tool_dir)):
        if fn.endswith(".done"):
            continue
        fp = os.path.join(tool_dir, fn)
        if os.path.isfile(fp):
            size_kb = os.path.getsize(fp) / 1024
            files.append((fn, f"{size_kb:.1f} KB", os.path.relpath(fp, os.path.dirname(output_html))))
    return files


# ---------------------------------------------------------------------------
# Build HTML
# ---------------------------------------------------------------------------
now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

sections = []  # (id, title, html_content)

# -- Overview section -------------------------------------------------------
overview_rows = ""
for s in samples:
    tools_run = []
    for tool in ["chloe", "pga", "plann", "cpgavas2", "mfannot", "fpma", "mitos", "mitoz", "trnascan", "aragorn", "liftoff", "ogdraw"]:
        done = os.path.join(outdir, tool, s, f"{s}.done")
        if os.path.exists(done):
            tools_run.append(tool)
    overview_rows += f"<tr><td>{html_mod.escape(s)}</td><td>{', '.join(tools_run) or 'none'}</td></tr>\n"

sections.append((
    "overview",
    "Pipeline Overview",
    f"""
    <p>Report generated: <strong>{now}</strong></p>
    <table>
      <thead><tr><th>Sample</th><th>Tools Run</th></tr></thead>
      <tbody>{overview_rows}</tbody>
    </table>
    """
))

# -- Per-tool sections -------------------------------------------------------
for tool_id, tool_title in TOOL_LABELS.items():
    tool_html = f"<p class='desc'>{TOOL_DESCRIPTIONS[tool_id]}</p>"
    has_data = False

    for s in samples:
        files = list_tool_outputs(s, tool_id, outdir)
        if not files:
            continue
        has_data = True
        file_rows = "".join(
            f"<tr><td><a href='{rel}'>{html_mod.escape(fn)}</a></td><td>{sz}</td></tr>"
            for fn, sz, rel in files
        )
        tool_html += f"""
        <h4>{html_mod.escape(s)}</h4>
        <table class="files">
          <thead><tr><th>File</th><th>Size</th></tr></thead>
          <tbody>{file_rows}</tbody>
        </table>
        """

    if not has_data:
        tool_html += "<p><em>No samples were processed with this tool.</em></p>"

    sections.append((f"tool-{tool_id}", tool_title, tool_html))

# -- QC Gene Completeness section -------------------------------------------
qc_html = ""
for s in samples:
    qc_path = os.path.join(outdir, "qc", "summary", f"{s}.qc_summary.tsv")
    rows = read_qc_summary(qc_path)
    if not rows:
        continue
    table_rows = ""
    for r in rows:
        gene_list = r.get("gene_names", "")
        short = gene_list[:120] + ("…" if len(gene_list) > 120 else "")
        table_rows += (
            f"<tr>"
            f"<td>{html_mod.escape(r.get('tool',''))}</td>"
            f"<td>{r.get('gene_count','0')}</td>"
            f"<td>{r.get('trna_count','0')}</td>"
            f"<td>{r.get('rrna_count','0')}</td>"
            f"<td title='{html_mod.escape(gene_list)}'>{html_mod.escape(short)}</td>"
            f"</tr>\n"
        )
    qc_html += f"""
    <h4>{html_mod.escape(s)}</h4>
    <table>
      <thead><tr><th>Tool</th><th>Genes</th><th>tRNAs</th><th>rRNAs</th><th>Gene Names</th></tr></thead>
      <tbody>{table_rows}</tbody>
    </table>
    """

if not qc_html:
    qc_html = "<p><em>No QC data available.</em></p>"

sections.append(("qc-genes", "Gene Completeness Summary", qc_html))

# -- BUSCO section -----------------------------------------------------------
busco_html = ""
for s in samples:
    bp = os.path.join(outdir, "qc", "busco", s, "short_summary.txt")
    metrics = read_busco_summary(bp)
    if not metrics:
        busco_html += f"<h4>{html_mod.escape(s)}</h4><p><em>No BUSCO results.</em></p>"
        continue
    metric_rows = "".join(
        f"<tr><td>{html_mod.escape(k)}</td><td>{html_mod.escape(str(v))}</td></tr>"
        for k, v in metrics.items()
    )
    busco_html += f"""
    <h4>{html_mod.escape(s)}</h4>
    <table>
      <thead><tr><th>Metric</th><th>Value</th></tr></thead>
      <tbody>{metric_rows}</tbody>
    </table>
    """

sections.append(("qc-busco", "BUSCO Assessment", busco_html))

# ---------------------------------------------------------------------------
# Assemble final HTML
# ---------------------------------------------------------------------------
nav_items = "".join(
    f'<li><a href="#{sid}">{html_mod.escape(stitle)}</a></li>'
    for sid, stitle, _ in sections
)

body_sections = "".join(
    f'<section id="{sid}"><h2>{html_mod.escape(stitle)}</h2>{scontent}</section>\n'
    for sid, stitle, scontent in sections
)

html_content = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Organelle Annotation Pipeline Report</title>
<style>
  :root {{
    --bg: #f8f9fa; --fg: #212529; --accent: #0d6efd;
    --card-bg: #fff; --border: #dee2e6; --hover: #e9ecef;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    background: var(--bg); color: var(--fg);
    display: flex; min-height: 100vh;
  }}
  nav {{
    position: fixed; top: 0; left: 0; width: 260px; height: 100vh;
    background: #1e293b; color: #e2e8f0; padding: 1.5rem 1rem;
    overflow-y: auto;
  }}
  nav h1 {{ font-size: 1.1rem; margin-bottom: 1rem; color: #93c5fd; }}
  nav ul {{ list-style: none; }}
  nav li {{ margin-bottom: .4rem; }}
  nav a {{
    color: #cbd5e1; text-decoration: none; font-size: .9rem;
    display: block; padding: .3rem .5rem; border-radius: 4px;
  }}
  nav a:hover {{ background: #334155; color: #fff; }}
  main {{
    margin-left: 260px; padding: 2rem 2.5rem; max-width: 1100px; width: 100%;
  }}
  h2 {{
    border-bottom: 2px solid var(--accent); padding-bottom: .4rem;
    margin: 2rem 0 1rem;
  }}
  h4 {{ margin: 1rem 0 .5rem; color: #475569; }}
  .desc {{ color: #64748b; margin-bottom: 1rem; font-style: italic; }}
  table {{
    width: 100%; border-collapse: collapse; margin-bottom: 1rem;
    background: var(--card-bg); border-radius: 6px; overflow: hidden;
    box-shadow: 0 1px 3px rgba(0,0,0,.08);
  }}
  th, td {{
    text-align: left; padding: .55rem .75rem; border-bottom: 1px solid var(--border);
  }}
  th {{ background: #f1f5f9; font-weight: 600; font-size: .85rem; color: #334155; }}
  tr:hover {{ background: var(--hover); }}
  table.files td:first-child {{ font-family: monospace; font-size: .85rem; }}
  a {{ color: var(--accent); }}
  section {{ margin-bottom: 2rem; }}
  footer {{ margin-top: 3rem; color: #94a3b8; font-size: .8rem; text-align: center; }}
  @media (max-width: 768px) {{
    nav {{ display: none; }}
    main {{ margin-left: 0; padding: 1rem; }}
  }}
</style>
</head>
<body>
<nav>
  <h1>&#x1F9EC; Organelle Pipeline</h1>
  <ul>{nav_items}</ul>
</nav>
<main>
  <h1>Organelle Annotation Pipeline Report</h1>
  {body_sections}
  <footer>
    Generated by Organelle Annotation Pipeline &middot; {now}
  </footer>
</main>
</body>
</html>
"""

os.makedirs(os.path.dirname(output_html), exist_ok=True)
with open(output_html, "w") as f:
    f.write(html_content)
