#!/usr/bin/env python3
"""
generate_downstream_report.py – Build HTML report for downstream analysis.

Called by Snakemake via ``script:`` directive; uses ``snakemake`` object.
"""

import os
import base64
from datetime import datetime

import pandas as pd

# ── Snakemake interface ──────────────────────────────────────────────────────

sample = snakemake.params.sample
species = snakemake.params.species
output_html = snakemake.output.html

rscu_tsv = snakemake.input.rscu_tsv
rscu_barplot = snakemake.input.rscu_barplot
codon_stats = snakemake.input.codon_stats
kaks_tsv = snakemake.input.kaks_tsv
gc_plot = snakemake.input.gc_plot
aa_plot = snakemake.input.aa_plot
tree_png = snakemake.input.tree_png
genome_map_png = snakemake.input.genome_map_png
synteny_plot = snakemake.input.synteny_plot
synteny_stats = snakemake.input.synteny_stats


# ── Helpers ──────────────────────────────────────────────────────────────────

def img_b64(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return ""
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()


def read_text(path):
    if not os.path.exists(path):
        return "File not found."
    with open(path) as f:
        return f.read()


def df_to_html(path, max_rows=50):
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return "<p><em>No data available.</em></p>"
        return df.head(max_rows).to_html(classes="table table-striped", index=False)
    except Exception as e:
        return f"<p>Error loading data: {e}</p>"


# ── Build sections ───────────────────────────────────────────────────────────

now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# Composition
gc_b64 = img_b64(gc_plot)
aa_b64 = img_b64(aa_plot)
composition_html = ""
if gc_b64 or aa_b64:
    gc_section = f'<div class="col-md-6 mb-3"><h4>GC Content</h4><img src="data:image/png;base64,{gc_b64}" class="plot-img"></div>' if gc_b64 else ""
    aa_section = f'<div class="col-md-6 mb-3"><h4>Amino Acid Composition</h4><img src="data:image/png;base64,{aa_b64}" class="plot-img"></div>' if aa_b64 else ""
    composition_html = f'<div class="section"><h2>Gene Composition Analysis</h2><div class="row">{gc_section}{aa_section}</div></div>'

# RSCU
rscu_b64 = img_b64(rscu_barplot)
rscu_table = df_to_html(rscu_tsv, 30)
rscu_html = f"""
<div class="section">
  <h2>Relative Synonymous Codon Usage (RSCU)</h2>
  <div class="row">
    <div class="col-md-12 mb-3">
      <img src="data:image/png;base64,{rscu_b64}" class="plot-img">
    </div>
    <div class="col-md-12">
      <h4>Data Table</h4>
      <div class="table-responsive">{rscu_table}</div>
    </div>
  </div>
</div>
"""

# Codon stats
codon_text = read_text(codon_stats)
codon_html = f"""
<div class="section">
  <h2>Start/Stop Codon Analysis</h2>
  <pre>{codon_text}</pre>
</div>
"""

# Ka/Ks
kaks_table = df_to_html(kaks_tsv)
kaks_html = f"""
<div class="section">
  <h2>Pairwise Ka/Ks Analysis</h2>
  <p>Ka/Ks estimation using MAFFT alignment + KaKs_Calculator.</p>
  <div class="table-responsive">{kaks_table}</div>
</div>
"""

# Phylogeny
tree_b64 = img_b64(tree_png)
tree_html = ""
if tree_b64:
    tree_html = f"""
    <div class="section">
      <h2>Phylogeny</h2>
      <p>Maximum Likelihood tree (IQ-TREE, GTR+G model, 1000 ultrafast bootstrap).</p>
      <img src="data:image/png;base64,{tree_b64}" class="plot-img">
    </div>
    """

# Genome map
gmap_b64 = img_b64(genome_map_png)
gmap_html = ""
if gmap_b64:
    gmap_html = f"""
    <div class="section">
      <h2>Genome Map</h2>
      <p>Circular genome visualisation generated with pyGenomeViz.</p>
      <div class="text-center">
        <img src="data:image/png;base64,{gmap_b64}" class="plot-img" style="max-width:90%">
      </div>
    </div>
    """

# Synteny
syn_b64 = img_b64(synteny_plot)
syn_html = ""
if syn_b64:
    syn_table = df_to_html(synteny_stats)
    syn_html = f"""
    <div class="section">
      <h2>Synteny Analysis</h2>
      <p>Genome structure comparison via MUMmer4/nucmer.</p>
      <div class="text-center mb-3">
        <img src="data:image/png;base64,{syn_b64}" class="plot-img" style="max-width:100%">
      </div>
      {syn_table}
    </div>
    """


# ── Assemble HTML ────────────────────────────────────────────────────────────

html_content = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Downstream Analysis – {sample}</title>
<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
<style>
  body {{ padding: 20px; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }}
  .section {{ margin-bottom: 40px; padding: 20px; background: #f8f9fa; border-radius: 8px; }}
  h1, h2 {{ color: #2c3e50; }}
  .plot-img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; padding: 5px; }}
  pre {{ background: #eaeaea; padding: 10px; border-radius: 5px; font-size: 0.85rem; max-height: 400px; overflow-y: auto; }}
  .table {{ font-size: 0.85rem; }}
</style>
</head>
<body>
<div class="container">
  <h1 class="text-center mb-2">Downstream Analysis Report</h1>
  <p class="text-center text-muted">
    Sample: <strong>{sample}</strong> &middot;
    Species: <strong>{species}</strong> &middot;
    Generated: {now}
  </p>

  {composition_html}
  {rscu_html}
  {codon_html}
  {kaks_html}
  {tree_html}
  {gmap_html}
  {syn_html}

  <footer class="text-center text-muted mt-5" style="font-size: 0.8rem;">
    Generated by Organelle Annotation Pipeline – Downstream Module &middot; {now}
  </footer>
</div>
</body>
</html>
"""

os.makedirs(os.path.dirname(output_html), exist_ok=True)
with open(output_html, "w") as f:
    f.write(html_content)

print(f"Downstream report: {output_html}")
