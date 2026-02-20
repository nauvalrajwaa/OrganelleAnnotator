# =============================================================================
# rules/report.smk – Aggregated HTML report
# =============================================================================

rule generate_report:
    """
    Produce an indexed HTML report with per-tool sections, QC summaries,
    gene-completeness tables, and links to all output files.
    """
    input:
        qc_summaries=lambda wc: [
            f"{OUTDIR}/qc/summary/{s}.qc_summary.tsv"
            for s in SAMPLES
            if config["qc"]["enabled"]
        ],
        busco_summaries=lambda wc: [
            f"{OUTDIR}/qc/busco/{s}/short_summary.txt"
            for s in SAMPLES
            if config["qc"]["enabled"]
        ],
        done_markers=lambda wc: [
            f"{OUTDIR}/{tool}/{s}/{s}.done"
            for s in SAMPLES
            for tool in tools_for_sample(s)
        ],
    output:
        html=f"{OUTDIR}/report/index.html",
    params:
        samples=SAMPLES,
        outdir=OUTDIR,
    log:
        f"{OUTDIR}/logs/report/report.log",
    script:
        "../scripts/generate_report.py"
