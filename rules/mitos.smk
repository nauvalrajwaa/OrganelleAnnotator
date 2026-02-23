# =============================================================================
# rules/mitos.smk – MITOS2 mitochondrial genome annotator (Docker)
# =============================================================================
# Docker image: quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0
# MITOS2 annotates protein-coding genes, tRNAs, and rRNAs in metazoan
# and plant mitochondrial genomes using a reference-based approach.
# =============================================================================

rule mitos_annotate:
    """
    Annotate a mitochondrial genome using MITOS2 via Docker/Singularity.
    Produces BED, GFF, FASTA (protein, nucleotide), and detailed results.
    """
    input:
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        done=touch(f"{OUTDIR}/{{sample}}/mitos/{{sample}}.done"),
        gff=f"{OUTDIR}/{{sample}}/mitos/result.gff",
        bed=f"{OUTDIR}/{{sample}}/mitos/result.bed",
    params:
        docker_image=config["mitos"]["docker_image"],
        genetic_code=lambda wc: get_genetic_code(wc.sample),
        ref_db=config["mitos"]["ref_db"],
        ref_dir=config["mitos"]["ref_dir"],
        extra=config["mitos"].get("extra", ""),
        out_dir=lambda wc: f"{OUTDIR}/{wc.sample}/mitos",
        use_singularity=config["mitos"].get("use_singularity", False),
    log:
        f"{OUTDIR}/{{sample}}/logs/mitos.log",
    threads:
        config["resources"]["mitos"]["threads"]
    resources:
        mem_mb=config["resources"]["mitos"]["mem_mb"],
        runtime=config["resources"]["mitos"]["runtime"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        # Resolve absolute paths for container volume mounts
        INPUT_ABS=$(cd "$(dirname {input.fasta})" && pwd)/$(basename {input.fasta})
        OUTPUT_ABS=$(mkdir -p {params.out_dir} && cd {params.out_dir} && pwd)

        # Copy input to output dir
        cp {input.fasta} {params.out_dir}/{wildcards.sample}.fasta

        REF_MOUNT=""
        REF_FLAG=""
        if [ -n "{params.ref_dir}" ] && [ -d "{params.ref_dir}" ]; then
            REF_ABS=$(cd "{params.ref_dir}" && pwd)
            REF_MOUNT="-v ${{REF_ABS}}:/ref_db:ro"
            REF_FLAG="--refdir /ref_db"
        fi

        if [ "{params.use_singularity}" = "True" ]; then
            # Singularity/Apptainer execution
            SINGULARITY_BIND="${{OUTPUT_ABS}}:/output"
            if [ -n "{params.ref_dir}" ] && [ -d "{params.ref_dir}" ]; then
                SINGULARITY_BIND="${{SINGULARITY_BIND}},$(cd {params.ref_dir} && pwd):/ref_db:ro"
                REF_FLAG="--refdir /ref_db"
            fi
            singularity exec \
                --bind "${{SINGULARITY_BIND}}" \
                docker://{params.docker_image} \
                runmitos.py \
                    --input /output/{wildcards.sample}.fasta \
                    --code {params.genetic_code} \
                    --outdir /output \
                    --refseqver {params.ref_db} \
                    $REF_FLAG \
                    {params.extra} \
                    2>&1 | tee {log}
        else
            # Docker execution
            docker run --rm \
                -v "${{OUTPUT_ABS}}":/output \
                $REF_MOUNT \
                {params.docker_image} \
                runmitos.py \
                    --input /output/{wildcards.sample}.fasta \
                    --code {params.genetic_code} \
                    --outdir /output \
                    --refseqver {params.ref_db} \
                    $REF_FLAG \
                    {params.extra} \
                    2>&1 | tee {log}
        fi

        # Ensure expected outputs exist
        touch {output.gff} {output.bed}
        """


rule mitos_extract_proteins:
    """
    Extract protein sequences from MITOS output for downstream QC.
    """
    input:
        done=f"{OUTDIR}/{{sample}}/mitos/{{sample}}.done",
    output:
        proteins=f"{OUTDIR}/{{sample}}/mitos/{{sample}}.proteins.fasta",
    params:
        out_dir=lambda wc: f"{OUTDIR}/{wc.sample}/mitos",
    run:
        import os, glob

        # MITOS outputs protein FASTA files in the result directory
        prot_files = glob.glob(os.path.join(params.out_dir, "*.fas")) + \
                     glob.glob(os.path.join(params.out_dir, "*protein*.fa*"))

        with open(output.proteins, "w") as out:
            for pf in prot_files:
                if "protein" in os.path.basename(pf).lower() or pf.endswith(".faa"):
                    with open(pf) as inf:
                        out.write(inf.read())

        # Also check for the standard MITOS protein output
        std_prot = os.path.join(params.out_dir, "result.faa")
        if os.path.exists(std_prot) and os.path.getsize(output.proteins) == 0:
            import shutil
            shutil.copy(std_prot, output.proteins)

        # Touch to ensure file exists even if empty
        if not os.path.exists(output.proteins):
            open(output.proteins, "w").close()
