# rules/gbdraw.smk – gbdraw: Genome diagram generator for organelles
#
# gbdraw runs ONCE PER ANNOTATION TOOL that may produce a GenBank (.gb/.gbk)
# file, creating maps per sample in {OUTDIR}/{sample}/gbdraw/{source_tool}/.
#
# Uses the gbdraw CLI for highly customized, publication-quality output:
#   - Separated strands (forward/reverse on different tracks)
#   - Specific feature filtering (CDS, rRNA, tRNA, etc.)
#   - Custom stroke widths, colors, and fonts
#   - Outputs strictly in PDF and PNG formats
#
# The script searches the annotator's output folder for *.gb / *.gbk files
# at runtime, so it works with any annotator regardless of naming convention.
#
# NOTE: GB_PRODUCING_TOOLS and gbdraw_source_tools() are defined in the
# main Snakefile (needed at parse time before rule all).

rule gbdraw_map:
    """
    Generate comprehensive genome diagrams from annotated GenBank
    files using the gbdraw CLI. Searches the source tool's output folder
    for *.gb / *.gbk files and draws each one found.

    Output includes:
      - PDF (vector, scalable and ideal for publications)
      - PNG (raster, easy to view for reports and presentations)
    """
    input:
        # Depend on the annotator's .done marker so we run AFTER annotation
        done = OUTDIR + "/{sample}/{source_tool}/{sample}.done",
    output:
        done = touch(OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}.done"),
        # Hapus .svg dan wajibkan .png serta .pdf sesuai script baru
        png  = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_map.png",
        pdf  = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_map.pdf",
    params:
        src_dir    = OUTDIR + "/{sample}/{source_tool}",
        out_prefix = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_map",
        out_dir    = OUTDIR + "/{sample}/gbdraw/{source_tool}",
        # Menambahkan parameter draw_mode (circular/linear) agar sejalan dengan Python script
        draw_mode  = config.get("gbdraw", {}).get("draw_mode", "circular"),
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