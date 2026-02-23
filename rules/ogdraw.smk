# =============================================================================
# rules/ogdraw.smk – OGDraw: OrganellarGenomeDRAW (circular map visualisation)
# =============================================================================
# OGDraw generates publication-quality circular and linear genome maps from
# GenBank-format annotation files. Works for both chloroplast and mitochondrial
# genomes. Docker image includes all Perl/BioPerl dependencies.
# Reference: Greiner et al. (2019) doi:10.1093/nar/gkz238
# =============================================================================

rule ogdraw_map:
    """
    Generate a circular genome map from a GenBank annotation file using OGDraw.
    Runs after a primary annotator has produced a .gb file.
    Outputs SVG and PNG circular maps.
    """
    input:
        gb=lambda wc: _ogdraw_input_gb(wc.sample),
    output:
        done=touch(f"{OUTDIR}/ogdraw/{{sample}}/{{sample}}.done"),
        svg=f"{OUTDIR}/ogdraw/{{sample}}/{{sample}}_map.svg",
    params:
        docker_image=config["ogdraw"]["docker_image"],
        out_dir=lambda wc: f"{OUTDIR}/ogdraw/{wc.sample}",
        extra=config["ogdraw"].get("extra", ""),
        use_singularity=config["ogdraw"].get("use_singularity", False),
    log:
        f"{OUTDIR}/logs/ogdraw/{{sample}}.log",
    threads: 1
    resources:
        mem_mb=config["resources"]["ogdraw"]["mem_mb"],
        runtime=config["resources"]["ogdraw"]["runtime"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        OUTPUT_ABS=$(cd {params.out_dir} && pwd)
        INPUT_ABS=$(cd "$(dirname {input.gb})" && pwd)/$(basename {input.gb})

        # Copy GenBank to output dir
        cp {input.gb} {params.out_dir}/{wildcards.sample}.gb

        if [ "{params.use_singularity}" = "True" ]; then
            singularity exec \
                --bind "${{OUTPUT_ABS}}":/data \
                docker://{params.docker_image} \
                drawgenemap \
                    --infile /data/{wildcards.sample}.gb \
                    --outfile /data/{wildcards.sample}_map \
                    --format svg \
                    {params.extra} \
                    2>&1 | tee {log}
        else
            docker run --rm \
                -v "${{OUTPUT_ABS}}":/data \
                -w /data \
                {params.docker_image} \
                drawgenemap \
                    --infile /data/{wildcards.sample}.gb \
                    --outfile /data/{wildcards.sample}_map \
                    --format svg \
                    {params.extra} \
                    2>&1 | tee {log}
        fi

        # Ensure output exists
        touch {output.svg}
        """


def _ogdraw_input_gb(sample):
    """
    Determine which GenBank file to use as OGDraw input.
    Prefers: chloe > pga > liftoff (for plastid)
             mfannot > mitos > mitoz > liftoff            (for mito)
    Falls back to whichever .gb exists.
    """
    organelle = samples_df.loc[sample, "organelle"]
    if organelle == "plastid":
        priority = ["chloe", "pga", "liftoff"]
    else:
        priority = ["mfannot", "mitos", "mitoz", "liftoff"]

    for tool in priority:
        gb = f"{OUTDIR}/{tool}/{sample}/{sample}.gb"
        gbk = f"{OUTDIR}/{tool}/{sample}/{sample}.gbk"
        # Return the path; Snakemake will resolve whether the file is produced
        if tool in tools_for_sample(sample):
            return gb
    # Fallback: first annotator's gb
    tools = tools_for_sample(sample)
    if tools:
        return f"{OUTDIR}/{tools[0]}/{sample}/{sample}.gb"
    return f"{OUTDIR}/chloe/{sample}/{sample}.gb"
