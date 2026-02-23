# =============================================================================
# rules/chloe.smk – Chloe.jl chloroplast annotation
# =============================================================================

rule chloe_annotate:
    """
    Annotate a chloroplast genome using Chloe.jl (Julia).
    Produces GFF (default), optionally GenBank, EMBL, SFF.
    """
    input:
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        done=touch(f"{OUTDIR}/chloe/{{sample}}/{{sample}}.done"),
        gff=f"{OUTDIR}/chloe/{{sample}}/{{sample}}.gff",
    params:
        chloe_dir=os.path.join(workflow.basedir, config["chloe"]["path"]),
        out_dir=lambda wc: f"{OUTDIR}/chloe/{wc.sample}",
        reference=config["chloe"]["reference"],
        references_dir=config["chloe"].get("references_dir", ""),
        sensitivity=config["chloe"]["sensitivity"],
        fmt_flags=lambda wc: " ".join(
            [f for f, v in {
                "--sff": config["chloe"]["sff"],
                "--gbk": config["chloe"]["gbk"],
                "--embl": config["chloe"]["embl"],
            }.items() if v]
        ),
        extra_flags=lambda wc: " ".join(
            [f for f, v in {
                "--no-filter": config["chloe"]["no_filter"],
                "--no-transform": config["chloe"]["no_transform"],
            }.items() if v]
        ),
    log:
        f"{OUTDIR}/logs/chloe/{{sample}}.log",
    threads:
        config["resources"]["chloe"]["threads"]
    resources:
        mem_mb=config["resources"]["chloe"]["mem_mb"],
        runtime=config["resources"]["chloe"]["runtime"],
    conda:
        "../envs/chloe.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        # Copy input fasta to output dir with sample name for Chloe
        cp {input.fasta} {params.out_dir}/{wildcards.sample}.fa

        REF_FLAG=""
        if [ -n "{params.references_dir}" ]; then
            REF_FLAG="-r {params.references_dir}"
        fi

        julia --project={params.chloe_dir} \
            -t {threads} \
            {params.chloe_dir}/chloe.jl annotate \
            $REF_FLAG \
            --sensitivity {params.sensitivity} \
            {params.fmt_flags} \
            {params.extra_flags} \
            -o {params.out_dir} \
            {params.out_dir}/{wildcards.sample}.fa \
            2>&1 | tee {log}

        # Rename outputs to standardised names
        for ext in gff sff gbk embl; do
            found=$(find {params.out_dir} -maxdepth 1 -name "*.${{ext}}" -o -name "*.chloe.${{ext}}" 2>/dev/null | head -1)
            if [ -n "$found" ] && [ "$found" != "{params.out_dir}/{wildcards.sample}.${{ext}}" ]; then
                cp "$found" "{params.out_dir}/{wildcards.sample}.${{ext}}" 2>/dev/null || true
            fi
        done

        # Ensure GFF exists even if empty
        touch {output.gff}
        """
