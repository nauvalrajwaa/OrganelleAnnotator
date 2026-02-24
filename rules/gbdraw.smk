# rules/gbdraw.smk – gbdraw: Genome diagram generator for organelles (circular map)
#
# gbdraw runs ONCE PER ANNOTATION TOOL that produces a GenBank (.gb/.gbk) file,
# creating multiple maps per sample in {OUTDIR}/{sample}/gbdraw/{source_tool}/.
#
# NOTE: GB_PRODUCING_TOOLS, GB_FILE_EXTENSION, and gbdraw_source_tools() are
# defined in the main Snakefile (needed at parse time before rule all).


def gbdraw_gb_path(sample, source_tool):
    """Return the GenBank file path produced by a given tool for a sample."""
    ext = GB_FILE_EXTENSION.get(source_tool, ".gb")
    return OUTDIR + "/" + sample + "/" + source_tool + "/" + sample + ext


rule gbdraw_map:
    """
    Generate comprehensive circular genome diagrams from an annotated GenBank
    file using gbdraw.  Produces multiple output formats (SVG, PNG, PDF, etc.)
    as configured in config.yaml → gbdraw.formats.
    Wildcard {source_tool} determines which tool's GenBank file is used.
    """
    input:
        gb = lambda wc: gbdraw_gb_path(wc.sample, wc.source_tool),
    output:
        done = touch(OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}.done"),
        svg  = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_map.svg",
    params:
        out_prefix = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_map",
        out_dir    = OUTDIR + "/{sample}/gbdraw/{source_tool}",
        formats    = config["gbdraw"].get("formats", ["svg"]),
        extra      = config["gbdraw"].get("extra", ""),
    log:
        OUTDIR + "/{sample}/logs/gbdraw_{source_tool}.log",
    threads: 1
    resources:
        mem_mb  = config["resources"]["gbdraw"]["mem_mb"],
        runtime = config["resources"]["gbdraw"]["runtime"],
    conda:
        "../envs/gbdraw.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        if [ ! -s {input.gb} ]; then
            echo "GenBank file is empty, skipping gbdraw for {wildcards.source_tool}" > {log}
            touch {output.svg}
            exit 0
        fi

        # Check that the GenBank file actually contains gene/CDS features
        if ! grep -qE '^\s+(gene|CDS|tRNA|rRNA)\s' {input.gb}; then
            echo "GenBank file has no gene annotations, skipping gbdraw for {wildcards.source_tool}" > {log}
            touch {output.svg}
            exit 0
        fi

        # Run gbdraw for each requested output format
        for fmt in {params.formats}; do
            echo "=== Generating $fmt format ===" >> {log}
            gbdraw circular \
                --gbk {input.gb} \
                -o {params.out_prefix} \
                -f "$fmt" \
                {params.extra} \
                2>&1 | tee -a {log}
        done

        # Ensure SVG output exists (required by Snakemake)
        touch {output.svg}
        """
