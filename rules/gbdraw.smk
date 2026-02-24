# rules/gbdraw.smk – gbdraw: Genome diagram generator for organelles
#
# gbdraw runs ONCE PER ANNOTATION TOOL that may produce a GenBank (.gb/.gbk)
# file, creating maps per sample in {OUTDIR}/{sample}/gbdraw/{source_tool}/.
#
# Uses the gbdraw CLI for highly customized, publication-quality output.
# Generates BOTH Circular and Linear maps automatically.

rule gbdraw_map:
    """
    Generate comprehensive genome diagrams from annotated GenBank
    files using the gbdraw CLI. Searches the source tool's output folder
    for *.gb / *.gbk files and draws both circular and linear maps.

    Output includes for BOTH Circular and Linear:
      - PDF (vector, scalable and ideal for publications)
      - PNG (raster, easy to view for reports and presentations)
    """
    input:
        done = OUTDIR + "/{sample}/{source_tool}/{sample}.done",
    output:
        done = touch(OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}.done"),
        # Output Sirkular
        circ_png = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_circular.png",
        circ_pdf = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_circular.pdf",
        # Output Linear
        lin_png  = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_linear.png",
        lin_pdf  = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_linear.pdf",
    params:
        src_dir     = OUTDIR + "/{sample}/{source_tool}",
        out_dir     = OUTDIR + "/{sample}/gbdraw/{source_tool}",
        circ_prefix = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_circular",
        lin_prefix  = OUTDIR + "/{sample}/gbdraw/{source_tool}/{sample}_linear",
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