# =============================================================================
# rules/mitoz.smk – MitoZ: animal mitochondrial genome annotator (Docker)
# =============================================================================
# Docker image: guanliangmeng/mitoz (https://github.com/linzhi2013/MitoZ)
# MitoZ annotates protein-coding genes, tRNAs, and rRNAs in animal
# mitochondrial genomes. It also generates a circular visualisation.
# =============================================================================

rule mitoz_annotate:
    """
    Annotate an animal mitochondrial genome using MitoZ via Docker.
    Produces GFF, GenBank, and circular visualisation.
    """
    input:
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        done=touch(f"{OUTDIR}/{{sample}}/mitoz/{{sample}}.done"),
        gff=f"{OUTDIR}/{{sample}}/mitoz/{{sample}}.gff",
    params:
        docker_image=config["mitoz"]["docker_image"],
        genetic_code=lambda wc: get_genetic_code(wc.sample),
        clade=config["mitoz"]["clade"],
        extra=config["mitoz"].get("extra", ""),
        out_dir=lambda wc: f"{OUTDIR}/{wc.sample}/mitoz",
        use_singularity=config["mitoz"].get("use_singularity", False),
    log:
        f"{OUTDIR}/{{sample}}/logs/mitoz.log",
    threads:
        config["resources"]["mitoz"]["threads"]
    resources:
        mem_mb=config["resources"]["mitoz"]["mem_mb"],
        runtime=config["resources"]["mitoz"]["runtime"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        # Resolve absolute paths for container volume mounts
        OUTPUT_ABS=$(mkdir -p {params.out_dir} && cd {params.out_dir} && pwd)

        # Copy input to output dir
        cp {input.fasta} {params.out_dir}/{wildcards.sample}.fasta

        if [ "{params.use_singularity}" = "True" ]; then
            singularity exec \
                --bind "${{OUTPUT_ABS}}":/output \
                docker://{params.docker_image} \
                mitoz annotate \
                    --fastafile /output/{wildcards.sample}.fasta \
                    --outprefix {wildcards.sample} \
                    --thread_number {threads} \
                    --clade {params.clade} \
                    --genetic_code {params.genetic_code} \
                    {params.extra} \
                    2>&1 | tee {log}
        else
            docker run --rm \
                -v "${{OUTPUT_ABS}}":/output \
                -w /output \
                {params.docker_image} \
                mitoz annotate \
                    --fastafile /output/{wildcards.sample}.fasta \
                    --outprefix {wildcards.sample} \
                    --thread_number {threads} \
                    --clade {params.clade} \
                    --genetic_code {params.genetic_code} \
                    {params.extra} \
                    2>&1 | tee {log}
        fi

        # MitoZ outputs into a subdirectory; locate and copy results
        RESULT_DIR=$(find {params.out_dir} -maxdepth 2 -name "*.gff" -exec dirname {{}} \; 2>/dev/null | head -1)
        if [ -n "$RESULT_DIR" ] && [ "$RESULT_DIR" != "{params.out_dir}" ]; then
            cp "$RESULT_DIR"/*.gff {params.out_dir}/{wildcards.sample}.gff 2>/dev/null || true
            cp "$RESULT_DIR"/*.gbk {params.out_dir}/{wildcards.sample}.gbk 2>/dev/null || true
            cp "$RESULT_DIR"/*.png {params.out_dir}/ 2>/dev/null || true
        fi

        # Ensure GFF exists even if empty
        touch {output.gff}
        """
