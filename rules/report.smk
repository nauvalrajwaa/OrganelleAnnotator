# rules/report.smk – Aggregated HTML report (per sample)

rule generate_report:
    """
    Produce an indexed HTML report with per-tool sections, QC summaries,
    gene-completeness tables, and links to all output files.
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
    output:
        html = OUTDIR + "/{sample}/report/index.html",
    params:
        samples = lambda wc: [wc.sample],
        outdir  = OUTDIR,
    log:
        OUTDIR + "/{sample}/logs/report.log",
    script:
        "../scripts/generate_report.py"
