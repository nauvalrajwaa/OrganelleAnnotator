# rules/gbdraw.smk – gbdraw: Genome diagram generator for organelles (circular map)
#
# gbdraw runs ONCE PER ANNOTATION TOOL that may produce a GenBank (.gb/.gbk)
# file, creating maps per sample in {OUTDIR}/{sample}/gbdraw/{source_tool}/.
#
# Uses the gbdraw Python API (see usage_gbdraw.md) for comprehensive output:
#   - Separated strands (forward/reverse on different tracks)
#   - Gene labels
#   - GC content / GC skew rings
#   - Legend
#   - Multiple output formats (SVG, PNG)
#
# The script searches the annotator's output folder for *.gb / *.gbk files
# at runtime, so it works with any annotator regardless of naming convention.
#
# NOTE: GB_PRODUCING_TOOLS and gbdraw_source_tools() are defined in the
# main Snakefile (needed at parse time before rule all).


rule gbdraw_map:
    """
    Generate comprehensive circular genome diagrams from annotated GenBank
    files using gbdraw Python API.  Searches the source tool's output folder
    for *.gb / *.gbk files and draws each one found.

    Output includes:
      - SVG (vector, editable in Inkscape/Illustrator)
      - PNG (raster, for reports and presentations)
      - Separated strands, gene labels, GC content, legend
    """
    input:
        # Depend on the annotator's .done marker so we run AFTER annotation
        done = OUTDIR + "/{sample}/{source_tool}/{sample}.done",
    output:
        done = touch(OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}.done"),
        svg  = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_map.svg",
    params:
        src_dir      = OUTDIR + "/{sample}/{source_tool}",
        out_prefix   = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_map",
        out_dir      = OUTDIR + "/{sample}/gbdraw/{source_tool}",
        formats      = config["gbdraw"].get("formats", ["svg"]),
        extra_config = config["gbdraw"].get("config_overrides", {}),
    log:
        OUTDIR + "/{sample}/logs/gbdraw_{source_tool}.log",
    threads: 1
    resources:
        mem_mb  = config["resources"]["gbdraw"]["mem_mb"],
        runtime = config["resources"]["gbdraw"]["runtime"],
    conda:
        "../envs/gbdraw.yaml"
    script:
        "../scripts/run_gbdraw.py"
