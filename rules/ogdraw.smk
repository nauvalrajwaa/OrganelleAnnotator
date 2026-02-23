# rules/ogdraw.smk – OGDraw: OrganellarGenomeDRAW (circular map visualisation)
#
# OGDraw runs ONCE PER ANNOTATION TOOL that produces a GenBank (.gb) file,
# creating multiple maps per sample in {OUTDIR}/{sample}/ogdraw/{source_tool}/.
#
# NOTE: GB_PRODUCING_TOOLS and ogdraw_source_tools() are defined in the
# main Snakefile (needed at parse time before rule all).


def ogdraw_gb_path(sample, source_tool):
    """Return the GenBank file path produced by a given tool for a sample."""
    return OUTDIR + "/" + sample + "/" + source_tool + "/" + sample + ".gb"


rule ogdraw_map:
    """
    Generate a circular genome map from a GenBank annotation file using OGDraw.
    Wildcard {source_tool} determines which tool's GenBank file is used.
    """
    input:
        gb = lambda wc: ogdraw_gb_path(wc.sample, wc.source_tool),
    output:
        done = touch(OUTDIR + "/{sample}/ogdraw/{source_tool}/{sample}.done"),
        svg  = OUTDIR + "/{sample}/ogdraw/{source_tool}/{sample}_map.svg",
    params:
        docker_image    = config["ogdraw"]["docker_image"],
        out_dir         = OUTDIR + "/{sample}/ogdraw/{source_tool}",
        extra           = config["ogdraw"].get("extra", ""),
        use_singularity = config["ogdraw"].get("use_singularity", False),
    log:
        OUTDIR + "/{sample}/logs/ogdraw_{source_tool}.log",
    threads: 1
    resources:
        mem_mb  = config["resources"]["ogdraw"]["mem_mb"],
        runtime = config["resources"]["ogdraw"]["runtime"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        if [ ! -s {input.gb} ]; then
            echo "GenBank file is empty, skipping OGDraw for {wildcards.source_tool}" > {log}
            touch {output.svg}
            exit 0
        fi

        OUTPUT_ABS=$(cd {params.out_dir} && pwd)

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

        touch {output.svg}
        """
