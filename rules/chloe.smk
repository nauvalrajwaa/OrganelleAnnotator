# rules/chloe.smk – Chloe.jl chloroplast genome annotator

rule chloe_annotate:
    """
    Annotate a chloroplast genome using Chloe.jl.
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
            {params.chloe_dir}/bin/chloe.jl \
            annotate \
            {input.fasta} \
            -o {params.out_dir}/{wildcards.sample}.gff \
            {params.extra} \
            2>&1 | tee {log}
        """
