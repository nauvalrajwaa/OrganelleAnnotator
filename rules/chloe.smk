# rules/chloe.smk – Chloe.jl chloroplast genome annotator

rule chloe_annotate:
    """
    Annotate a chloroplast genome using Chloe.jl.
    Chloe's -o flag sets the output DIRECTORY; the actual GFF file is
    written as <basename>.gff inside that directory.
    """
    input:
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        done = touch(OUTDIR + "/{sample}/chloe/{sample}.done"),
        gff  = OUTDIR + "/{sample}/chloe/{sample}.gff",
    params:
        out_dir    = OUTDIR + "/{sample}/chloe",
        chloe_dir  = config["chloe"]["path"],
        extra      = config["chloe"].get("extra", ""),
    log:
        OUTDIR + "/{sample}/logs/chloe.log",
    threads:
        config["resources"]["chloe"]["threads"]
    resources:
        mem_mb  = config["resources"]["chloe"]["mem_mb"],
        runtime = config["resources"]["chloe"]["runtime"],
    conda:
        "../envs/chloe.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        julia --project={params.chloe_dir} \
            {params.chloe_dir}/chloe.jl \
            annotate \
            {input.fasta} \
            -o {params.out_dir} \
            {params.extra} \
            2>&1 | tee {log}

        # Chloe writes <input_basename>.gff inside the output dir.
        # Find and rename to the expected output name.
        CHLOE_GFF=$(find {params.out_dir} -maxdepth 1 -name "*.gff" -type f 2>/dev/null | head -1)
        if [ -n "$CHLOE_GFF" ] && [ "$CHLOE_GFF" != "{output.gff}" ]; then
            mv "$CHLOE_GFF" {output.gff}
        fi
        touch {output.gff}
        """
