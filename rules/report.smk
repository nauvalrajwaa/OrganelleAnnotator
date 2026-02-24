# rules/report.smk – Aggregated HTML report

rule generate_report:
    """
    Produce an indexed HTML report with per-tool sections, QC summaries,
    gene-completeness tables, and links to all output files.
    """
    input:
        qc_summaries = lambda wc: [
            OUTDIR + "/" + s + "/qc/qc_summary.tsv"
            for s in SAMPLES
            if config["qc"]["enabled"]
        ],
        busco_summaries = lambda wc: [
            OUTDIR + "/" + s + "/qc/busco/short_summary.txt"
            for s in SAMPLES
            if config["qc"]["enabled"]
        ],
        done_markers = lambda wc: [
            OUTDIR + "/" + s + "/" + tool + "/" + s + ".done"
            for s in SAMPLES
            for tool in tools_for_sample(s)
        ],
        gbdraw_markers = lambda wc: [
            OUTDIR + "/" + s + "/gbdraw/" + src + "/" + s + ".done"
            for s in SAMPLES
            for src in gbdraw_source_tools(s)
        ],
    output:
        html = OUTDIR + "/report/index.html",
    params:
        samples = SAMPLES,
        outdir  = OUTDIR,
    log:
        OUTDIR + "/report/report.log",
    script:
        "../scripts/generate_report.py"
