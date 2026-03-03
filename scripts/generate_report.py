#!/usr/bin/env python3
"""
generate_report.py – Build a comprehensive, interactive HTML report for the
Organelle Annotation Pipeline.

Merges:
  - Per-tool annotation results (filtered by tools_select from config.yaml)
  - gbdraw genome diagrams (embedded PNG images, circular & linear)
  - QC gene-completeness tables & BUSCO
  - Downstream analyses (RSCU, codon usage, Ka/Ks, composition, phylogeny,
    genome map, synteny) – previously in a separate downstream_report.html

Called by Snakemake via ``script:`` directive; uses ``snakemake`` object for I/O.
"""

import csv
import os
import base64
import html as html_mod
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
outdir = snakemake.params.outdir
samples = snakemake.params.samples
output_html = snakemake.output.html
tools_select = list(snakemake.params.tools_select)
downstream_enabled = snakemake.params.downstream_enabled
species_name = getattr(snakemake.params, "species_name", "")

# Tool metadata
TOOL_LABELS = {
    "chloe":    "Chloë (Chloe.jl) — Chloroplast Annotator",
    "pga":      "PGA — Plastid Genome Annotator",
    "mfannot":  "MFannot — Mitochondrial/Plastid Annotator",
    "fpma":     "fpma — Fast Plant Mitochondria Annotator",
    "mitos":    "MITOS2 — Mitochondrial Genome Annotator",
    "mitoz":    "MitoZ — Animal Mitochondrial Genome Annotator",
    "trnascan": "tRNAscan-SE — tRNA Gene Detection",
    "aragorn":  "Aragorn — tRNA/tmRNA Detection",
    "liftoff":  "Liftoff — Reference-based Annotation Lift-over",
    "gbdraw":   "gbdraw — Genome Map Visualisation",
}

TOOL_DESCRIPTIONS = {
    "chloe":    "Julia-based chloroplast genome annotator using XGBoost models and suffix-array alignment.",
    "pga":      "Perl/BLAST pipeline for rapid batch annotation of plastid genomes against GenBank references.",
    "mfannot":  "Comprehensive mitochondrial/plastid annotator (Docker: nbeck/mfannot) using BLAST, HMMER, Exonerate, Erpin.",
    "fpma":     "Rust-based fast HMM scanner for presence/absence of mitochondrial genes using HMMER3 nhmmer.",
    "mitos":    "Reference-based mitochondrial genome annotator for protein-coding genes, tRNAs, and rRNAs.",
    "mitoz":    "Docker-based animal mitochondrial genome annotator with circular visualisation.",
    "trnascan": "Gold-standard tRNA detection tool using covariance models. Supports organellar/mitochondrial mode.",
    "aragorn":  "Lightweight tRNA and tmRNA detection using homology search. Fast and suitable for organelle genomes.",
    "liftoff":  "Minimap2-based annotation lift-over from a reference organelle genome.",
    "gbdraw":   "Publication-quality circular and linear genome diagrams from GenBank files (bioconda::gbdraw).",
}

TOOL_ICONS = {
    "chloe": "&#x1F9EC;", "pga": "&#x1F33F;", "mfannot": "&#x1F52C;",
    "fpma": "&#x26A1;", "mitos": "&#x1F9EB;", "mitoz": "&#x1F41F;",
    "trnascan": "&#x1F50E;", "aragorn": "&#x1F3F9;",
    "liftoff": "&#x1F4D0;", "gbdraw": "&#x1F3A8;",
}

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def img_b64(path):
    """Return base64-encoded data URI for an image, or empty string."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return ""
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode()
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg",
            "svg": "image/svg+xml", "gif": "image/gif"}.get(ext, "image/png")
    return f"data:{mime};base64,{data}"


def read_busco_summary(path):
    metrics = {}
    if not os.path.exists(path):
        return metrics
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "Complete BUSCOs" in line:
                metrics["Complete BUSCOs"] = line.split()[0] if line[0].isdigit() else line
            elif "Complete and single-copy" in line:
                metrics["Single-copy"] = line.split()[0] if line[0].isdigit() else line
            elif "Complete and duplicated" in line:
                metrics["Duplicated"] = line.split()[0] if line[0].isdigit() else line
            elif "Fragmented" in line:
                metrics["Fragmented"] = line.split()[0] if line[0].isdigit() else line
            elif "Missing" in line:
                metrics["Missing"] = line.split()[0] if line[0].isdigit() else line
            elif "Total" in line and "BUSCO" in line:
                metrics["Total"] = line.split()[0] if line[0].isdigit() else line
    return metrics


def read_qc_summary(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def list_tool_outputs(sample, tool, outdir_path):
    """List result files produced by a tool for a sample."""
    tool_dir = os.path.join(outdir_path, sample, tool)
    if not os.path.isdir(tool_dir):
        return []
    files = []
    for fn in sorted(os.listdir(tool_dir)):
        if fn.endswith(".done"):
            continue
        fp = os.path.join(tool_dir, fn)
        if os.path.isfile(fp):
            size_kb = os.path.getsize(fp) / 1024
            files.append((fn, f"{size_kb:.1f} KB",
                          os.path.relpath(fp, os.path.dirname(output_html))))
    return files


def read_tsv_to_html_table(path, max_rows=50, table_id=""):
    """Read a TSV file and return an HTML table string."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return "<p class='text-muted'><em>No data available.</em></p>"
    try:
        rows_data = []
        with open(path) as f:
            reader = csv.reader(f, delimiter="\t")
            headers = next(reader, None)
            if not headers:
                return "<p class='text-muted'><em>No data available.</em></p>"
            for row in reader:
                rows_data.append(row)
                if len(rows_data) >= max_rows:
                    break

        tid = f' id="{table_id}"' if table_id else ""
        header_html = "".join(f"<th>{html_mod.escape(str(c))}</th>" for c in headers)
        rows_html = ""
        for row in rows_data:
            cells = "".join(f"<td>{html_mod.escape(str(v))}</td>" for v in row)
            rows_html += f"<tr>{cells}</tr>\n"
        total = len(rows_data)
        total_note = f"<p class='text-muted small'>Showing {total} rows (max {max_rows}).</p>" if total >= max_rows else ""
        return f"""
        <div class="table-responsive">
          <table class="data-table sortable"{tid}>
            <thead><tr>{header_html}</tr></thead>
            <tbody>{rows_html}</tbody>
          </table>
        </div>
        {total_note}
        """
    except Exception as e:
        return f"<p class='text-muted'>Error loading data: {html_mod.escape(str(e))}</p>"


def read_text_file(path):
    if not os.path.exists(path):
        return ""
    with open(path) as f:
        return f.read()


# ---------------------------------------------------------------------------
# Build report sections
# ---------------------------------------------------------------------------
now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
sections = []  # list of (id, icon, title, html_content)

# ---- SECTION: Pipeline Overview -------------------------------------------
overview_rows = ""
for s in samples:
    tools_run = []
    for tool in tools_select:
        if tool == "gbdraw":
            continue
        done = os.path.join(outdir, s, tool, f"{s}.done")
        if os.path.exists(done):
            tools_run.append(tool)
    # gbdraw
    gbdraw_base = os.path.join(outdir, s, "gbdraw")
    if os.path.isdir(gbdraw_base):
        gbdraw_srcs = [d for d in os.listdir(gbdraw_base)
                       if os.path.isdir(os.path.join(gbdraw_base, d))
                       and os.path.exists(os.path.join(gbdraw_base, d, f"{s}.done"))]
        if gbdraw_srcs:
            tools_run.append(f"gbdraw ({', '.join(sorted(gbdraw_srcs))})")

    badges = "".join(f'<span class="badge">{html_mod.escape(t)}</span>' for t in tools_run)
    overview_rows += f"""
    <tr>
      <td><strong>{html_mod.escape(s)}</strong></td>
      <td>{badges or '<em>none</em>'}</td>
      <td><span class="badge badge-info">{len(tools_run)}</span></td>
    </tr>
    """

species_line = f"<p><strong>Species:</strong> <em>{html_mod.escape(species_name)}</em></p>" if species_name else ""

sections.append((
    "overview", "&#x1F4CB;", "Pipeline Overview",
    f"""
    <div class="card">
      <p><strong>Report generated:</strong> {now}</p>
      <p><strong>Active tools:</strong> {', '.join(tools_select)}</p>
      {species_line}
      <table class="data-table">
        <thead><tr><th>Sample</th><th>Tools Completed</th><th>Count</th></tr></thead>
        <tbody>{overview_rows}</tbody>
      </table>
    </div>
    """
))

# ---- SECTIONS: Per-Tool Results (only tools in tools_select) ---------------
for tool_id in tools_select:
    if tool_id == "gbdraw":
        continue  # handled separately

    tool_title = TOOL_LABELS.get(tool_id, tool_id)
    tool_desc = TOOL_DESCRIPTIONS.get(tool_id, "")
    tool_icon = TOOL_ICONS.get(tool_id, "&#x1F527;")
    tool_html = f"<p class='tool-desc'>{html_mod.escape(tool_desc)}</p>"
    has_data = False

    for s in samples:
        files = list_tool_outputs(s, tool_id, outdir)
        if not files:
            continue
        has_data = True
        file_rows = "".join(
            f"<tr><td><code><a href='{html_mod.escape(rel)}'>{html_mod.escape(fn)}</a></code></td>"
            f"<td>{sz}</td></tr>"
            for fn, sz, rel in files
        )
        tool_html += f"""
        <div class="card">
          <h4>{html_mod.escape(s)}</h4>
          <table class="data-table files-table">
            <thead><tr><th>File</th><th>Size</th></tr></thead>
            <tbody>{file_rows}</tbody>
          </table>
        </div>
        """

    if not has_data:
        tool_html += "<div class='card'><p class='text-muted'><em>No samples were processed with this tool.</em></p></div>"

    sections.append((f"tool-{tool_id}", tool_icon, tool_title, tool_html))


# ---- SECTION: gbdraw Genome Visualisation (embedded PNGs) [4] ---------------
if "gbdraw" in tools_select:
    gbdraw_html = f"<p class='tool-desc'>{html_mod.escape(TOOL_DESCRIPTIONS['gbdraw'])}</p>"
    has_gbdraw = False

    for s in samples:
        gbdraw_base = os.path.join(outdir, s, "gbdraw")
        if not os.path.isdir(gbdraw_base):
            continue

        for src_tool in sorted(os.listdir(gbdraw_base)):
            src_dir = os.path.join(gbdraw_base, src_tool)
            if not os.path.isdir(src_dir):
                continue

            # Find PNG files
            png_files = [fn for fn in sorted(os.listdir(src_dir))
                         if fn.lower().endswith(".png")]

            if png_files:
                has_gbdraw = True
                gbdraw_html += f"""
                <div class="card">
                  <h4>{html_mod.escape(s)} &mdash; from <em>{html_mod.escape(src_tool)}</em></h4>
                  <div class="gbdraw-gallery">
                """
                for png_fn in png_files:
                    png_path = os.path.join(src_dir, png_fn)
                    b64 = img_b64(png_path)
                    if b64:
                        label = png_fn.replace(".png", "").replace("_", " ").title()
                        gbdraw_html += f"""
                        <div class="gbdraw-item">
                          <h5>{html_mod.escape(label)}</h5>
                          <img src="{b64}" alt="{html_mod.escape(png_fn)}"
                               class="genome-map-img zoomable"
                               onclick="openLightbox(this)">
                          <p class="img-caption">Click to enlarge &middot; {html_mod.escape(png_fn)}</p>
                        </div>
                        """
                gbdraw_html += "</div></div>"

            # List all files for download
            all_files = []
            for fn in sorted(os.listdir(src_dir)):
                fp = os.path.join(src_dir, fn)
                if not fn.endswith(".done") and os.path.isfile(fp):
                    sz = f"{os.path.getsize(fp)/1024:.1f} KB"
                    rel = os.path.relpath(fp, os.path.dirname(output_html))
                    all_files.append((fn, sz, rel))
            if all_files:
                file_rows = "".join(
                    f"<tr><td><code><a href='{html_mod.escape(rel)}'>{html_mod.escape(fn)}</a></code></td>"
                    f"<td>{sz}</td></tr>"
                    for fn, sz, rel in all_files
                )
                gbdraw_html += f"""
                <details class="file-details">
                  <summary>&#x1F4C1; All output files ({len(all_files)})</summary>
                  <table class="data-table files-table">
                    <thead><tr><th>File</th><th>Size</th></tr></thead>
                    <tbody>{file_rows}</tbody>
                  </table>
                </details>
                """

    if not has_gbdraw:
        gbdraw_html += "<div class='card'><p class='text-muted'><em>No genome diagrams available.</em></p></div>"

    sections.append(("tool-gbdraw", "&#x1F3A8;", TOOL_LABELS["gbdraw"], gbdraw_html))


# ---- SECTION: QC Gene Completeness ----------------------------------------
qc_html = ""
for s in samples:
    qc_path = os.path.join(outdir, s, "qc", "qc_summary.tsv")
    rows = read_qc_summary(qc_path)
    if not rows:
        continue

    # Filter to only show tools in tools_select
    rows = [r for r in rows if r.get("tool", "") in tools_select]
    if not rows:
        continue

    table_rows = ""
    for r in rows:
        gene_list = r.get("gene_names", "")
        genes = [g.strip() for g in gene_list.split(";") if g.strip() and g.strip() != "unknown"]
        gene_count = r.get("gene_count", "0")
        trna_count = r.get("trna_count", "0")
        rrna_count = r.get("rrna_count", "0")

        if genes:
            shown = genes[:15]
            gene_badges = " ".join(f'<span class="gene-badge">{html_mod.escape(g)}</span>' for g in shown)
            if len(genes) > 15:
                gene_badges += f' <span class="text-muted">+{len(genes)-15} more</span>'
        else:
            gene_badges = '<em class="text-muted">No named genes detected</em>'

        table_rows += f"""
        <tr>
          <td><span class="badge">{html_mod.escape(r.get('tool',''))}</span></td>
          <td class="text-center"><strong>{gene_count}</strong></td>
          <td class="text-center">{trna_count}</td>
          <td class="text-center">{rrna_count}</td>
          <td>{gene_badges}</td>
        </tr>
        """

    qc_html += f"""
    <div class="card">
      <h4>{html_mod.escape(s)}</h4>
      <table class="data-table">
        <thead><tr>
          <th>Tool</th><th class="text-center">Total Genes</th>
          <th class="text-center">tRNAs</th><th class="text-center">rRNAs</th>
          <th>Gene Names</th>
        </tr></thead>
        <tbody>{table_rows}</tbody>
      </table>
    </div>
    """

if not qc_html:
    qc_html = "<div class='card'><p class='text-muted'><em>No QC data available.</em></p></div>"

sections.append(("qc-genes", "&#x2705;", "Gene Completeness Summary", qc_html))


# ---- SECTION: BUSCO -------------------------------------------------------
busco_html = ""
for s in samples:
    bp = os.path.join(outdir, s, "qc", "busco", "short_summary.txt")
    metrics = read_busco_summary(bp)
    if not metrics:
        busco_html += f"""
        <div class="card">
          <h4>{html_mod.escape(s)}</h4>
          <p class="text-muted"><em>No BUSCO results available.</em></p>
        </div>
        """
        continue

    metric_bars = ""
    total = int(metrics.get("Total", "0") or "0")
    for k, v in metrics.items():
        if k == "Total":
            continue
        val = int(v) if str(v).isdigit() else 0
        pct = (val / total * 100) if total > 0 else 0
        color_map = {"Complete BUSCOs": "#22c55e", "Single-copy": "#3b82f6",
                     "Duplicated": "#f59e0b", "Fragmented": "#f97316",
                     "Missing": "#ef4444"}
        color = color_map.get(k, "#6b7280")
        metric_bars += f"""
        <div class="busco-bar-row">
          <span class="busco-label">{html_mod.escape(k)}</span>
          <div class="busco-bar-bg">
            <div class="busco-bar-fill" style="width:{pct:.1f}%;background:{color}"></div>
          </div>
          <span class="busco-val">{v} ({pct:.1f}%)</span>
        </div>
        """
    busco_html += f"""
    <div class="card">
      <h4>{html_mod.escape(s)}</h4>
      <p class="text-muted">Total BUSCOs: <strong>{metrics.get('Total', 'N/A')}</strong></p>
      {metric_bars}
    </div>
    """

sections.append(("qc-busco", "&#x1F4CA;", "BUSCO Assessment", busco_html))


# ===========================================================================
# DOWNSTREAM ANALYSIS SECTIONS (merged from downstream_report.html) [1]
# ===========================================================================
if downstream_enabled:
    for s in samples:
        ds_base = os.path.join(outdir, s, "downstream")
        if not os.path.isdir(ds_base):
            continue

        # -- Composition -------------------------------------------------------
        gc_plot_path = os.path.join(ds_base, "composition", "gc_content_plot.png")
        aa_plot_path = os.path.join(ds_base, "composition", "aa_composition_plot.png")
        gc_b64 = img_b64(gc_plot_path)
        aa_b64 = img_b64(aa_plot_path)

        comp_html = ""
        if gc_b64 or aa_b64:
            gc_img = f'<div class="plot-half"><h5>GC Content</h5><img src="{gc_b64}" class="plot-img zoomable" onclick="openLightbox(this)"></div>' if gc_b64 else ""
            aa_img = f'<div class="plot-half"><h5>Amino Acid Composition</h5><img src="{aa_b64}" class="plot-img zoomable" onclick="openLightbox(this)"></div>' if aa_b64 else ""
            comp_html = f'<div class="card"><div class="plot-row">{gc_img}{aa_img}</div></div>'
        else:
            comp_html = "<div class='card'><p class='text-muted'><em>No composition data available.</em></p></div>"

        sections.append(("ds-composition", "&#x1F9EA;", "Gene Composition Analysis", comp_html))

        # -- RSCU --------------------------------------------------------------
        rscu_barplot = os.path.join(ds_base, "rscu", "rscu_barplot.png")
        rscu_heatmap = os.path.join(ds_base, "rscu", "rscu_heatmap.png")
        rscu_tsv_path = os.path.join(ds_base, "rscu", "rscu.tsv")

        rscu_bar_b64 = img_b64(rscu_barplot)
        rscu_heat_b64 = img_b64(rscu_heatmap)
        rscu_table = read_tsv_to_html_table(rscu_tsv_path, 30, "rscu-table")

        rscu_imgs = ""
        if rscu_bar_b64:
            rscu_imgs += f'<div class="mb-3"><h5>RSCU Bar Plot</h5><img src="{rscu_bar_b64}" class="plot-img zoomable" onclick="openLightbox(this)"></div>'
        if rscu_heat_b64:
            rscu_imgs += f'<div class="mb-3"><h5>RSCU Heatmap</h5><img src="{rscu_heat_b64}" class="plot-img zoomable" onclick="openLightbox(this)"></div>'

        rscu_html_section = f"""
        <div class="card">
          {rscu_imgs}
          <details open>
            <summary>&#x1F4CB; RSCU Data Table</summary>
            {rscu_table}
          </details>
        </div>
        """
        sections.append(("ds-rscu", "&#x1F4C8;", "Relative Synonymous Codon Usage (RSCU)", rscu_html_section))

        # -- Codon Analysis ----------------------------------------------------
        codon_stats_path = os.path.join(ds_base, "codons", "codon_stats.txt")
        codon_text = read_text_file(codon_stats_path)
        codon_html = f"""
        <div class="card">
          <pre class="code-block">{html_mod.escape(codon_text) if codon_text else '<em>No codon data available.</em>'}</pre>
        </div>
        """
        sections.append(("ds-codons", "&#x1F524;", "Start / Stop Codon Analysis", codon_html))

        # -- Ka/Ks -------------------------------------------------------------
        kaks_tsv_path = os.path.join(ds_base, "kaks", "kaks_summary.tsv")
        kaks_table = read_tsv_to_html_table(kaks_tsv_path, 50, "kaks-table")
        kaks_html = f"""
        <div class="card">
          <p class="tool-desc">Ka/Ks estimation using MAFFT alignment + KaKs_Calculator.</p>
          {kaks_table}
        </div>
        """
        sections.append(("ds-kaks", "&#x1F52C;", "Pairwise Ka/Ks Analysis", kaks_html))

        # -- Phylogeny ---------------------------------------------------------
        tree_png_path = os.path.join(ds_base, "phylogeny", "tree_plot.png")
        tree_b64 = img_b64(tree_png_path)
        if tree_b64:
            tree_html = f"""
            <div class="card">
              <p class="tool-desc">Maximum Likelihood tree (IQ-TREE, GTR+G model, 1000 ultrafast bootstrap).</p>
              <div class="text-center">
                <img src="{tree_b64}" class="plot-img zoomable" onclick="openLightbox(this)" style="max-width:90%">
              </div>
            </div>
            """
        else:
            tree_html = "<div class='card'><p class='text-muted'><em>No phylogeny data available.</em></p></div>"
        sections.append(("ds-phylogeny", "&#x1F333;", "Phylogenetic Tree", tree_html))

        # -- Genome Map --------------------------------------------------------
        gmap_png_path = os.path.join(ds_base, "genome_map", "genome_map.png")
        gmap_b64 = img_b64(gmap_png_path)
        if gmap_b64:
            gmap_html = f"""
            <div class="card">
              <p class="tool-desc">Circular genome visualisation generated with pyGenomeViz.</p>
              <div class="text-center">
                <img src="{gmap_b64}" class="plot-img zoomable" onclick="openLightbox(this)" style="max-width:85%">
              </div>
            </div>
            """
        else:
            gmap_html = "<div class='card'><p class='text-muted'><em>No genome map available.</em></p></div>"
        sections.append(("ds-genomemap", "&#x1F5FA;", "Genome Map", gmap_html))

        # -- Synteny -----------------------------------------------------------
        syn_plot_path = os.path.join(ds_base, "synteny", "synteny_plot.png")
        syn_stats_path = os.path.join(ds_base, "synteny", "synteny_stats.tsv")
        syn_b64 = img_b64(syn_plot_path)
        syn_table = read_tsv_to_html_table(syn_stats_path, 50, "synteny-table")
        if syn_b64:
            syn_html = f"""
            <div class="card">
              <p class="tool-desc">Genome structure comparison via MUMmer4/nucmer.</p>
              <div class="text-center mb-3">
                <img src="{syn_b64}" class="plot-img zoomable" onclick="openLightbox(this)" style="max-width:100%">
              </div>
              <details>
                <summary>&#x1F4CB; Synteny Statistics</summary>
                {syn_table}
              </details>
            </div>
            """
        else:
            syn_html = "<div class='card'><p class='text-muted'><em>No synteny data available.</em></p></div>"
        sections.append(("ds-synteny", "&#x1F517;", "Synteny Analysis", syn_html))


# ---------------------------------------------------------------------------
# Assemble final HTML
# ---------------------------------------------------------------------------

# Navigation items
nav_items = ""
ds_nav_started = False
for sid, icon, stitle, _ in sections:
    if sid.startswith("ds-") and not ds_nav_started:
        nav_items += '<li class="nav-divider">Downstream Analysis</li>'
        ds_nav_started = True
    active_cls = ' class="active"' if sid == "overview" else ""
    nav_items += f'<li><a href="#{sid}"{active_cls}>{icon} {html_mod.escape(stitle)}</a></li>\n'

# Body sections
body_sections = ""
for sid, icon, stitle, scontent in sections:
    body_sections += f"""
    <section id="{sid}" class="report-section">
      <div class="section-header" onclick="toggleSection(this)">
        <h2>{icon} {html_mod.escape(stitle)}</h2>
        <span class="collapse-icon">&#x25BC;</span>
      </div>
      <div class="section-body">
        {scontent}
      </div>
    </section>
    """

# ---------------------------------------------------------------------------
# Full HTML page
# ---------------------------------------------------------------------------
html_content = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Organelle Annotation Report</title>
<style>
  :root {
    --bg: #f0f2f5;
    --fg: #1a1a2e;
    --accent: #0d6efd;
    --accent-hover: #0b5ed7;
    --card-bg: #ffffff;
    --card-shadow: 0 2px 8px rgba(0,0,0,.06);
    --border: #e2e8f0;
    --nav-bg: #0f172a;
    --nav-fg: #cbd5e1;
    --nav-active: #3b82f6;
    --nav-hover: #1e293b;
    --badge-bg: #e0e7ff;
    --badge-fg: #3730a3;
    --badge-info-bg: #dbeafe;
    --badge-info-fg: #1d4ed8;
    --code-bg: #f1f5f9;
    --th-bg: #f8fafc;
    --hover-row: #f1f5f9;
    --muted: #64748b;
    --success: #22c55e;
    --warning: #f59e0b;
    --danger: #ef4444;
    --transition: 0.2s ease;
  }

  [data-theme="dark"] {
    --bg: #0f172a;
    --fg: #e2e8f0;
    --card-bg: #1e293b;
    --card-shadow: 0 2px 8px rgba(0,0,0,.3);
    --border: #334155;
    --nav-bg: #020617;
    --nav-fg: #94a3b8;
    --nav-hover: #1e293b;
    --badge-bg: #312e81;
    --badge-fg: #a5b4fc;
    --badge-info-bg: #1e3a5f;
    --badge-info-fg: #93c5fd;
    --code-bg: #1e293b;
    --th-bg: #1e293b;
    --hover-row: #334155;
    --muted: #94a3b8;
  }

  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  html { scroll-behavior: smooth; }

  body {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 "Helvetica Neue", Arial, sans-serif;
    background: var(--bg); color: var(--fg);
    display: flex; min-height: 100vh; line-height: 1.6;
    transition: background var(--transition), color var(--transition);
  }

  /* Sidebar */
  nav {
    position: fixed; top: 0; left: 0;
    width: 280px; height: 100vh;
    background: var(--nav-bg); color: var(--nav-fg);
    padding: 0; overflow-y: auto; z-index: 100;
    transition: transform 0.3s ease, background var(--transition);
    display: flex; flex-direction: column;
  }
  .nav-header {
    padding: 1.5rem 1.25rem 1rem;
    border-bottom: 1px solid rgba(255,255,255,.08); flex-shrink: 0;
  }
  .nav-header h1 { font-size: 1.05rem; font-weight: 700; color: #93c5fd; letter-spacing: -0.02em; }
  .nav-header p { font-size: 0.75rem; color: #64748b; margin-top: 0.25rem; }
  .nav-controls {
    padding: 0.75rem 1.25rem;
    border-bottom: 1px solid rgba(255,255,255,.08); flex-shrink: 0;
  }
  .nav-search {
    width: 100%; padding: 0.45rem 0.75rem;
    border: 1px solid #334155; border-radius: 6px;
    background: #1e293b; color: #e2e8f0; font-size: 0.8rem; outline: none;
    transition: border-color var(--transition);
  }
  .nav-search:focus { border-color: #3b82f6; }
  .nav-search::placeholder { color: #475569; }
  nav ul { list-style: none; padding: 0.5rem 0.75rem; flex: 1; overflow-y: auto; }
  nav li { margin-bottom: 2px; }
  nav li.nav-divider {
    font-size: 0.7rem; font-weight: 700; text-transform: uppercase;
    letter-spacing: 0.08em; color: #475569;
    padding: 1rem 0.5rem 0.3rem; margin-top: 0.5rem;
  }
  nav a {
    color: var(--nav-fg); text-decoration: none; font-size: 0.82rem;
    display: block; padding: 0.4rem 0.6rem; border-radius: 6px;
    transition: background var(--transition), color var(--transition);
    white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
  }
  nav a:hover, nav a.active { background: var(--nav-hover); color: #fff; }
  nav a.active { background: var(--nav-active); color: #fff; font-weight: 600; }

  /* Main */
  main {
    margin-left: 280px; padding: 2rem 2.5rem;
    max-width: 1200px; width: 100%;
    transition: margin-left 0.3s ease;
  }
  .main-header {
    display: flex; justify-content: space-between; align-items: center;
    margin-bottom: 2rem; flex-wrap: wrap; gap: 1rem;
  }
  .main-header h1 { font-size: 1.6rem; font-weight: 700; color: var(--fg); }
  .header-actions { display: flex; gap: 0.5rem; align-items: center; }

  /* Buttons */
  .btn {
    padding: 0.4rem 0.9rem; border: 1px solid var(--border); border-radius: 6px;
    background: var(--card-bg); color: var(--fg); font-size: 0.8rem;
    cursor: pointer; transition: all var(--transition);
    display: inline-flex; align-items: center; gap: 0.3rem;
  }
  .btn:hover { background: var(--hover-row); border-color: var(--accent); }
  .btn-primary { background: var(--accent); color: #fff; border-color: var(--accent); }
  .btn-primary:hover { background: var(--accent-hover); }

  /* Sections */
  .report-section {
    margin-bottom: 1.5rem; background: var(--card-bg);
    border-radius: 10px; box-shadow: var(--card-shadow); overflow: hidden;
    transition: background var(--transition), box-shadow var(--transition);
  }
  .section-header {
    display: flex; justify-content: space-between; align-items: center;
    padding: 1rem 1.5rem; cursor: pointer; user-select: none;
    border-bottom: 1px solid var(--border);
    transition: background var(--transition);
  }
  .section-header:hover { background: var(--hover-row); }
  .section-header h2 { font-size: 1.15rem; font-weight: 600; margin: 0; border: none; padding: 0; }
  .collapse-icon { font-size: 0.8rem; color: var(--muted); transition: transform 0.3s ease; }
  .section-header.collapsed .collapse-icon { transform: rotate(-90deg); }
  .section-body { padding: 1.25rem 1.5rem; transition: max-height 0.4s ease, padding 0.3s ease; overflow: hidden; }
  .section-body.collapsed { max-height: 0 !important; padding-top: 0; padding-bottom: 0; }

  /* Cards */
  .card {
    padding: 1rem 1.25rem; margin-bottom: 1rem;
    background: var(--bg); border-radius: 8px; border: 1px solid var(--border);
    transition: background var(--transition), border-color var(--transition);
  }
  h4 { margin: 0 0 0.75rem; color: var(--fg); font-size: 1rem; font-weight: 600; }
  h5 { margin: 0 0 0.5rem; color: var(--muted); font-size: 0.9rem; font-weight: 600; }

  /* Tables */
  .table-responsive { overflow-x: auto; }
  .data-table {
    width: 100%; border-collapse: collapse; font-size: 0.85rem;
    background: var(--card-bg); border-radius: 8px; overflow: hidden;
  }
  .data-table th, .data-table td {
    text-align: left; padding: 0.6rem 0.85rem; border-bottom: 1px solid var(--border);
  }
  .data-table th {
    background: var(--th-bg); font-weight: 600; font-size: 0.8rem;
    color: var(--muted); text-transform: uppercase; letter-spacing: 0.04em;
    position: sticky; top: 0;
  }
  .data-table tbody tr:hover { background: var(--hover-row); }
  .data-table tbody tr:last-child td { border-bottom: none; }
  .sortable th { cursor: pointer; position: relative; }
  .sortable th:hover { color: var(--accent); }
  .sortable th::after { content: " \\21C5"; font-size: 0.7rem; color: var(--muted); }
  .files-table td:first-child { font-family: 'JetBrains Mono', 'Fira Code', monospace; font-size: 0.8rem; }
  .text-center { text-align: center; }
  .mb-3 { margin-bottom: 1rem; }

  /* Badges */
  .badge {
    display: inline-block; padding: 0.15rem 0.55rem; border-radius: 12px;
    font-size: 0.75rem; font-weight: 600;
    background: var(--badge-bg); color: var(--badge-fg); margin: 0.1rem;
  }
  .badge-info { background: var(--badge-info-bg); color: var(--badge-info-fg); }
  .gene-badge {
    display: inline-block; padding: 0.1rem 0.4rem; border-radius: 4px;
    font-size: 0.72rem; font-family: monospace;
    background: var(--code-bg); color: var(--fg); margin: 1px; border: 1px solid var(--border);
  }

  /* Typography */
  .tool-desc { color: var(--muted); font-style: italic; margin-bottom: 1rem; font-size: 0.88rem; }
  .text-muted { color: var(--muted); }
  .small { font-size: 0.8rem; }
  a { color: var(--accent); text-decoration: none; }
  a:hover { text-decoration: underline; }
  code { font-family: 'JetBrains Mono', 'Fira Code', monospace; font-size: 0.85em; background: var(--code-bg); padding: 0.1rem 0.35rem; border-radius: 4px; }
  pre, .code-block {
    background: var(--code-bg); padding: 1rem; border-radius: 8px;
    font-size: 0.82rem; font-family: 'JetBrains Mono', 'Fira Code', monospace;
    max-height: 400px; overflow-y: auto; white-space: pre-wrap;
    word-break: break-all; border: 1px solid var(--border);
  }

  /* BUSCO bars */
  .busco-bar-row { display: flex; align-items: center; gap: 0.75rem; margin-bottom: 0.4rem; }
  .busco-label { font-size: 0.82rem; min-width: 110px; text-align: right; }
  .busco-bar-bg { flex: 1; height: 20px; background: var(--border); border-radius: 10px; overflow: hidden; }
  .busco-bar-fill { height: 100%; border-radius: 10px; transition: width 0.6s ease; }
  .busco-val { font-size: 0.8rem; min-width: 90px; color: var(--muted); }

  /* Images */
  .plot-img {
    max-width: 100%; height: auto; border: 1px solid var(--border);
    border-radius: 8px; padding: 4px; background: #fff;
    transition: transform 0.2s ease, box-shadow 0.2s ease;
  }
  .zoomable { cursor: zoom-in; }
  .zoomable:hover { transform: scale(1.01); box-shadow: 0 4px 16px rgba(0,0,0,.12); }
  .plot-row { display: flex; gap: 1.5rem; flex-wrap: wrap; }
  .plot-half { flex: 1; min-width: 300px; }
  .gbdraw-gallery { display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 1.5rem; margin: 1rem 0; }
  .gbdraw-item { text-align: center; }
  .genome-map-img { max-width: 100%; border-radius: 8px; }
  .img-caption { font-size: 0.75rem; color: var(--muted); margin-top: 0.3rem; }

  /* Lightbox */
  .lightbox-overlay {
    display: none; position: fixed; top: 0; left: 0;
    width: 100vw; height: 100vh; background: rgba(0,0,0,.85);
    z-index: 9999; justify-content: center; align-items: center;
    cursor: zoom-out; backdrop-filter: blur(4px);
  }
  .lightbox-overlay.active { display: flex; }
  .lightbox-overlay img { max-width: 92vw; max-height: 92vh; border-radius: 8px; box-shadow: 0 8px 32px rgba(0,0,0,.4); }
  .lightbox-close {
    position: fixed; top: 1rem; right: 1.5rem;
    color: #fff; font-size: 2rem; cursor: pointer; z-index: 10000;
    opacity: 0.7; transition: opacity var(--transition);
  }
  .lightbox-close:hover { opacity: 1; }

  /* Details */
  details { margin: 0.75rem 0; border: 1px solid var(--border); border-radius: 6px; overflow: hidden; }
  summary {
    padding: 0.6rem 1rem; cursor: pointer; font-weight: 600; font-size: 0.88rem;
    background: var(--th-bg); user-select: none;
  }
  summary:hover { background: var(--hover-row); }
  details[open] > summary { border-bottom: 1px solid var(--border); }
  details > *:not(summary) { padding: 0.75rem 1rem; }

  /* Scroll top */
  .scroll-top {
    position: fixed; bottom: 2rem; right: 2rem;
    width: 42px; height: 42px; border-radius: 50%;
    background: var(--accent); color: #fff; border: none;
    font-size: 1.2rem; cursor: pointer; display: none;
    align-items: center; justify-content: center;
    box-shadow: 0 2px 12px rgba(0,0,0,.2); z-index: 50;
    transition: opacity var(--transition), transform var(--transition);
  }
  .scroll-top:hover { transform: translateY(-2px); }
  .scroll-top.visible { display: flex; }

  /* Footer */
  footer {
    margin-top: 3rem; padding: 1.5rem 0; border-top: 1px solid var(--border);
    color: var(--muted); font-size: 0.78rem; text-align: center;
  }

  /* Mobile menu */
  .menu-toggle {
    display: none; position: fixed; top: 1rem; left: 1rem; z-index: 200;
    background: var(--accent); color: #fff; border: none; border-radius: 8px;
    padding: 0.5rem 0.75rem; font-size: 1.1rem; cursor: pointer;
  }
  @media (max-width: 900px) {
    .menu-toggle { display: block; }
    nav { transform: translateX(-100%); }
    nav.open { transform: translateX(0); }
    main { margin-left: 0; padding: 1rem; padding-top: 3.5rem; }
    .gbdraw-gallery { grid-template-columns: 1fr; }
    .plot-row { flex-direction: column; }
  }
  @media print {
    nav, .menu-toggle, .scroll-top, .header-actions { display: none !important; }
    main { margin-left: 0; max-width: 100%; }
    .report-section { break-inside: avoid; box-shadow: none; border: 1px solid #ccc; }
    .section-body.collapsed { max-height: none !important; padding: 1rem !important; }
    .lightbox-overlay { display: none !important; }
  }
</style>
</head>
<body>

<button class="menu-toggle" onclick="document.querySelector('nav').classList.toggle('open')" aria-label="Toggle menu">&#x2630;</button>

<nav>
  <div class="nav-header">
    <h1>&#x1F9EC; Organelle Pipeline</h1>
    <p>Annotation Report</p>
  </div>
  <div class="nav-controls">
    <input type="text" class="nav-search" placeholder="&#x1F50D; Search sections..." oninput="filterNav(this.value)">
  </div>
  <ul id="nav-list">
    """ + nav_items + """
  </ul>
  <div style="padding:0.75rem 1.25rem;border-top:1px solid rgba(255,255,255,.08);flex-shrink:0;">
    <button class="btn" onclick="toggleTheme()" style="width:100%;justify-content:center;background:#1e293b;color:#e2e8f0;border-color:#334155">
      &#x1F319; Toggle Dark Mode
    </button>
  </div>
</nav>

<main>
  <div class="main-header">
    <h1>Organelle Annotation Report</h1>
    <div class="header-actions">
      <button class="btn" onclick="expandAll()">&#x1F4C2; Expand All</button>
      <button class="btn" onclick="collapseAll()">&#x1F4C1; Collapse All</button>
      <button class="btn btn-primary" onclick="window.print()">&#x1F5A8; Print</button>
    </div>
  </div>

  """ + body_sections + """

  <footer>
    Generated by <strong>Organelle Annotation Pipeline</strong> &middot; """ + now + """
    <br>Tools: """ + ", ".join(tools_select) + (
    f" &middot; Species: <em>{html_mod.escape(species_name)}</em>" if species_name else ""
) + """
  </footer>
</main>

<div class="lightbox-overlay" id="lightbox" onclick="closeLightbox()">
  <span class="lightbox-close" onclick="closeLightbox()">&times;</span>
  <img id="lightbox-img" src="" alt="Enlarged view">
</div>

<button class="scroll-top" id="scrollTopBtn" onclick="window.scrollTo({top:0,behavior:'smooth'})">&#x2191;</button>

<script>
function toggleTheme() {
  var html = document.documentElement;
  var current = html.getAttribute('data-theme');
  html.setAttribute('data-theme', current === 'dark' ? 'light' : 'dark');
  localStorage.setItem('theme', html.getAttribute('data-theme'));
}
(function() {
  var saved = localStorage.getItem('theme');
  if (saved) document.documentElement.setAttribute('data-theme', saved);
})();

function toggleSection(header) {
  header.classList.toggle('collapsed');
  var body = header.nextElementSibling;
  body.classList.toggle('collapsed');
}

function expandAll() {
  document.querySelectorAll('.section-header').forEach(function(h) { h.classList.remove('collapsed'); });
  document.querySelectorAll('.section-body').forEach(function(b) { b.classList.remove('collapsed'); });
}

function collapseAll() {
  document.querySelectorAll('.section-header').forEach(function(h) { h.classList.add('collapsed'); });
  document.querySelectorAll('.section-body').forEach(function(b) { b.classList.add('collapsed'); });
}

function filterNav(query) {
  var q = query.toLowerCase();
  document.querySelectorAll('#nav-list li').forEach(function(li) {
    if (li.classList.contains('nav-divider')) { li.style.display = ''; return; }
    var text = li.textContent.toLowerCase();
    li.style.display = text.indexOf(q) >= 0 ? '' : 'none';
  });
}

var secs = document.querySelectorAll('.report-section');
var navLinks = document.querySelectorAll('#nav-list a');
function updateActiveNav() {
  var current = '';
  secs.forEach(function(s) {
    if (s.getBoundingClientRect().top <= 120) current = s.id;
  });
  navLinks.forEach(function(a) {
    if (a.getAttribute('href') === '#' + current) a.classList.add('active');
    else a.classList.remove('active');
  });
}
window.addEventListener('scroll', function() {
  updateActiveNav();
  var btn = document.getElementById('scrollTopBtn');
  if (window.scrollY > 400) btn.classList.add('visible');
  else btn.classList.remove('visible');
});

function openLightbox(img) {
  var overlay = document.getElementById('lightbox');
  document.getElementById('lightbox-img').src = img.src;
  overlay.classList.add('active');
  document.body.style.overflow = 'hidden';
}
function closeLightbox() {
  document.getElementById('lightbox').classList.remove('active');
  document.body.style.overflow = '';
}
document.addEventListener('keydown', function(e) {
  if (e.key === 'Escape') closeLightbox();
});

document.querySelectorAll('.sortable th').forEach(function(th, idx) {
  th.addEventListener('click', function() {
    var table = th.closest('table');
    var tbody = table.querySelector('tbody');
    var rows = Array.from(tbody.querySelectorAll('tr'));
    var dir = th.dataset.sortDir === 'asc' ? 'desc' : 'asc';
    th.dataset.sortDir = dir;
    rows.sort(function(a, b) {
      var aVal = (a.cells[idx] ? a.cells[idx].textContent.trim() : '');
      var bVal = (b.cells[idx] ? b.cells[idx].textContent.trim() : '');
      var aNum = parseFloat(aVal);
      var bNum = parseFloat(bVal);
      if (!isNaN(aNum) && !isNaN(bNum)) return dir === 'asc' ? aNum - bNum : bNum - aNum;
      return dir === 'asc' ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
    });
    rows.forEach(function(r) { tbody.appendChild(r); });
  });
});

document.querySelectorAll('#nav-list a').forEach(function(a) {
  a.addEventListener('click', function() {
    document.querySelector('nav').classList.remove('open');
  });
});
</script>
</body>
</html>
"""

# Write output
os.makedirs(os.path.dirname(output_html), exist_ok=True)
with open(output_html, "w") as f:
    f.write(html_content)

print(f"Report written to: {output_html}")
