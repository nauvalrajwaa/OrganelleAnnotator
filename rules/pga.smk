# =============================================================================
# rules/pga.smk – PGA (Plastid Genome Annotator)
# =============================================================================

rule pga_annotate:
    """
    Annotate a plastid genome using PGA (Perl + BLAST).
    Input:  single FASTA → copied into a per-sample target directory.
    Output: GenBank (.gb) file produced by PGA.
    """
    input:
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        done=touch(f"{OUTDIR}/pga/{{sample}}/{{sample}}.done"),
        gb=f"{OUTDIR}/pga/{{sample}}/{{sample}}.gb",
    params:
        pga_dir=config["pga"]["path"],
        ref_dir=config["pga"]["reference_dir"],
        form=config["pga"]["form"],
        ir_min=config["pga"]["ir_min"],
        pidentity=config["pga"]["pidentity"],
        qcoverage=config["pga"]["qcoverage"],
        out_dir=lambda wc: f"{OUTDIR}/pga/{wc.sample}",
        target_dir=lambda wc: f"{OUTDIR}/pga/{wc.sample}/target",
        gb_dir=lambda wc: f"{OUTDIR}/pga/{wc.sample}/gb",
    log:
        f"{OUTDIR}/logs/pga/{{sample}}.log",
    threads:
        config["resources"]["pga"]["threads"]
    resources:
        mem_mb=config["resources"]["pga"]["mem_mb"],
        runtime=config["resources"]["pga"]["runtime"],
    conda:
        "../envs/pga.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.target_dir} {params.gb_dir} $(dirname {log})

        # PGA expects a directory of FASTA files
        cp {input.fasta} {params.target_dir}/{wildcards.sample}.fasta

        perl {params.pga_dir}/PGA.pl \
            -r {params.ref_dir} \
            -t {params.target_dir} \
            -o {params.gb_dir} \
            -f {params.form} \
            -i {params.ir_min} \
            -p {params.pidentity} \
            -q {params.qcoverage} \
            -l {params.out_dir}/warning \
            2>&1 | tee {log}

        # Move generated GenBank file to standard name
        found=$(find {params.gb_dir} -maxdepth 1 -name "*.gb" 2>/dev/null | head -1)
        if [ -n "$found" ]; then
            cp "$found" {output.gb}
        else
            touch {output.gb}
        fi
        """
