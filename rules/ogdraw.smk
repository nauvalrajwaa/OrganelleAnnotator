# =============================================================================
# rules/ogdraw.smk – OGDraw: OrganellarGenomeDRAW (circular map visualisation)
# =============================================================================
# OGDraw generates publication-quality circular and linear genome maps from
# GenBank-format annotation files. Works for both chloroplast and mitochondrial
# genomes. Docker image includes all Perl/BioPerl dependencies.
# Reference: Greiner et al. (2019) doi:10.1093/nar/gkz238
#
# OGDraw runs ONCE PER ANNOTATION TOOL that produces a GenBank (.gb/.gbk) file,
# creating multiple maps per sample in {OUTDIR}/{sample}/ogdraw/{source_tool}/.
# =============================================================================

# Tools that can produce GenBank files for OGDraw input
GB_PRODUCING_TOOLS = {
    "plastid": ["chloe", "pga", "liftoff"],
    "mito":    ["liftoff"],
}


def ogdraw_source_tools(sample):
    """Return list of tools that produce GenBank files for this sample."""
    organelle = samples_df.loc[sample, "organelle"]
    possible = GB_PRODUCING_TOOLS.get(organelle, GB_PRODUCING_TOOLS.get("plastid", []))
    # Only include tools that are actually enabled for this sample
    sample_tools = tools_for_sample(sample)
    return [t for t in possible if t in sample_tools]


def ogdraw_gb_path(sample, source_tool):
    """Return the GenBank file path produced by a given tool for a sample."""
    return f"{OUTDIR}/{sample}/{source_tool}/{sample}.gb"


rule ogdraw_map:
    """
    Generate a circular genome map from a GenBank annotation file using OGDraw.
    Runs for EACH annotation tool that produces a .gb file.
    Wildcard {source_tool} determines which tool's GenBank file is used.
    """
    input:
        gb=lambda wc: ogdraw_gb_path(wc.sample, wc.source_tool),
    output:
        done=touch(f"{OUTDIR}/{{sample}}/ogdraw/{{source_tool}}/{{sample}}.done"),
        svg=f"{OUTDIR}/{{sample}}/ogdraw/{{source_tool}}/{{sample}}_map.svg",
    params:
        docker_image=config["ogdraw"]["docker_image"],
        out_dir=lambda wc: f"{OUTDIR}/{wc.sample}/ogdraw/{wc.source_tool}",
        extra=config["ogdraw"].get("extra", ""),
        use_singularity=config["ogdraw"].get("use_singularity", False),
    log:
        f"{OUTDIR}/{{sample}}/logs/ogdraw_{{source_tool}}.log",
    threads: 1
    resources:
        mem_mb=config["resources"]["ogdraw"]["mem_mb"],
        runtime=config["resources"]["ogdraw"]["runtime"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        # Skip if GenBank file is empty
        if [ ! -s {input.gb} ]; then
            echo "GenBank file is empty, skipping OGDraw for {wildcards.source_tool}" > {log}
            touch {output.svg}
            exit 0
        fi

        OUTPUT_ABS=$(cd {params.out_dir} && pwd)

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
