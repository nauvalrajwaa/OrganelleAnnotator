# rules/report.smk – Comprehensive merged HTML report (per sample)
# Merges annotation results + downstream analysis into a single interactive page.

rule generate_report:
    """
    Produce a comprehensive, interactive HTML report that merges:
      - Per-tool annotation results (filtered by tools_select)
      - gbdraw genome diagrams (embedded PNGs)
      - QC gene-completeness tables
      - Downstream analyses (RSCU, codon, Ka/Ks, composition, phylogeny, etc.)
    One report is generated per sample inside results/{sample}/report/.
    """
    input:
        qc_summary = lambda wc: (
            OUTDIR + "/" + wc.sample + "/qc/qc_summary.tsv"
            if config["qc"]["enabled"] else []
        ),
        busco_summary = lambda wc: (
            OUTDIR + "/" + wc.sample + "/qc/busco/short_summary.txt"
            if config["qc"]["enabled"]
               and samples_df.loc[wc.sample, "organelle"] in ("plastid", "mito")
            else []
        ),
        done_markers = lambda wc: [
            OUTDIR + "/" + wc.sample + "/" + tool + "/" + wc.sample + ".done"
            for tool in tools_for_sample(wc.sample)
        ],
        gbdraw_markers = lambda wc: [
            OUTDIR + "/" + wc.sample + "/gbdraw/" + src + "/" + wc.sample + ".done"
            for src in gbdraw_source_tools(wc.sample)
        ],
        # Downstream inputs (optional - may not exist)
        downstream_report = lambda wc: (
            OUTDIR + "/" + wc.sample + "/downstream/downstream_report.html"
            if config.get("downstream", {}).get("enabled", False)
            else []
        ),
    output:
        html = OUTDIR + "/{sample}/report/index.html",
    params:
        samples      = lambda wc: [wc.sample],
        outdir       = OUTDIR,
        tools_select = config.get("tools_select", []),
        downstream_enabled = config.get("downstream", {}).get("enabled", False),
        species_name = config.get("downstream", {}).get("species_name", ""),
    log:
        OUTDIR + "/{sample}/logs/report.log",
    script:
        "../scripts/generate_report.py"
