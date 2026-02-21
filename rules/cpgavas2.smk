# =============================================================================
# rules/cpgavas2.smk – CPGAVAS2: Chloroplast Genome Annotation (Docker)
# =============================================================================
# CPGAVAS2 annotates chloroplast genomes using BLAST + HMMER against a curated
# plant chloroplast protein database. Detects IRs, produces GenBank output
# and circular genome maps.
# Docker image: equipped with all dependencies.
# Reference: Shi et al. (2019) doi:10.1093/nar/gkz345
# =============================================================================

rule cpgavas2_annotate:
    """
    Annotate a chloroplast genome using CPGAVAS2 via Docker.
    Produces GenBank, GFF, and optional circular map image.
    """
    input:
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        done=touch(f"{OUTDIR}/cpgavas2/{{sample}}/{{sample}}.done"),
        gb=f"{OUTDIR}/cpgavas2/{{sample}}/{{sample}}.gb",
        gff=f"{OUTDIR}/cpgavas2/{{sample}}/{{sample}}.gff",
    params:
        docker_image=config["cpgavas2"]["docker_image"],
        genetic_code=lambda wc: get_genetic_code(wc.sample),
        extra=config["cpgavas2"].get("extra", ""),
        out_dir=lambda wc: f"{OUTDIR}/cpgavas2/{wc.sample}",
        use_singularity=config["cpgavas2"].get("use_singularity", False),
    log:
        f"{OUTDIR}/logs/cpgavas2/{{sample}}.log",
    threads:
        config["resources"]["cpgavas2"]["threads"]
    resources:
        mem_mb=config["resources"]["cpgavas2"]["mem_mb"],
        runtime=config["resources"]["cpgavas2"]["runtime"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        OUTPUT_ABS=$(mkdir -p {params.out_dir} && cd {params.out_dir} && pwd)

        # Copy input to output dir
        cp {input.fasta} {params.out_dir}/{wildcards.sample}.fasta

        if [ "{params.use_singularity}" = "True" ]; then
            singularity exec \
                --bind "${{OUTPUT_ABS}}":/data \
                docker://{params.docker_image} \
                cpgavas2.pl \
                    -i /data/{wildcards.sample}.fasta \
                    -o /data/output \
                    -c {params.genetic_code} \
                    -t {threads} \
                    {params.extra} \
                    2>&1 | tee {log}
        else
            docker run --rm \
                -v "${{OUTPUT_ABS}}":/data \
                -w /data \
                {params.docker_image} \
                cpgavas2.pl \
                    -i /data/{wildcards.sample}.fasta \
                    -o /data/output \
                    -c {params.genetic_code} \
                    -t {threads} \
                    {params.extra} \
                    2>&1 | tee {log}
        fi

        # Locate output files from CPGAVAS2 result directory
        RESULT_DIR="{params.out_dir}/output"
        if [ -d "$RESULT_DIR" ]; then
            # GenBank
            found_gb=$(find "$RESULT_DIR" -maxdepth 2 -name "*.gb" -o -name "*.gbk" 2>/dev/null | head -1)
            if [ -n "$found_gb" ]; then
                cp "$found_gb" {output.gb}
            fi
            # GFF
            found_gff=$(find "$RESULT_DIR" -maxdepth 2 -name "*.gff" -o -name "*.gff3" 2>/dev/null | head -1)
            if [ -n "$found_gff" ]; then
                cp "$found_gff" {output.gff}
            fi
            # Copy any map images
            find "$RESULT_DIR" -maxdepth 2 \( -name "*.png" -o -name "*.svg" \) \
                -exec cp {{}} {params.out_dir}/ \; 2>/dev/null || true
        fi

        # Ensure outputs exist
        touch {output.gb} {output.gff}
        """
