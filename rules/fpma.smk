# rules/fpma.smk – fpma (Fast Plant Mitochondria Annotator)

rule fpma_annotate:
    """
    Quick presence/absence scan of mitochondrial genes using HMM profiles.
    Produces GFF3, TSV summary table, and optional HTML SVG plot.
    """
    input:
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        done = touch(OUTDIR + "/{sample}/fpma/{sample}.done"),
        gff  = OUTDIR + "/{sample}/fpma/{sample}.gff",
        tsv  = OUTDIR + "/{sample}/fpma/{sample}.presence.tsv",
    params:
        fpma_dir    = config["fpma"]["path"],
        nhmmer_path = config["fpma"]["nhmmer_path"],
        hmms_subdir = config["fpma"]["hmms_subdir"],
        e_value     = config["fpma"]["e_value"],
        plot_flag   = lambda wc: (
            "--plot " + OUTDIR + "/" + wc.sample + "/fpma/" + wc.sample + ".html"
            if config["fpma"]["plot"] else ""
        ),
        out_dir     = OUTDIR + "/{sample}/fpma",
    log:
        OUTDIR + "/{sample}/logs/fpma.log",
    threads:
        config["resources"]["fpma"]["threads"]
    resources:
        mem_mb  = config["resources"]["fpma"]["mem_mb"],
        runtime = config["resources"]["fpma"]["runtime"],
    conda:
        "../envs/fpma.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        FPMA_BIN={params.fpma_dir}/target/release/fpma
        if [ ! -x "$FPMA_BIN" ]; then
            FPMA_BIN=$(command -v fpma 2>/dev/null || echo "")
        fi
        if [ -z "$FPMA_BIN" ]; then
            echo "ERROR: fpma binary not found." >&2
            exit 1
        fi

        "$FPMA_BIN" \
            --plant-mito {input.fasta} \
            --nhmmer-path {params.nhmmer_path} \
            --hmms-path {params.fpma_dir}/hmms/{params.hmms_subdir}/ \
            --e-value {params.e_value} \
            --gff {output.gff} \
            {params.plot_flag} \
            > {output.tsv} \
            2> {log}
        """
