# rules/pga.smk – PGA (Plastid Genome Annotator)

rule pga_annotate:
    """
    Annotate a plastid genome using PGA (Perl + BLAST).
    Input:  single FASTA → copied into a per-sample target directory.
    Output: GenBank (.gb) file produced by PGA.
    """
    input:
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        done = touch(OUTDIR + "/{sample}/pga/{sample}.done"),
        gb   = OUTDIR + "/{sample}/pga/{sample}.gb",
    params:
        pga_dir    = os.path.join(workflow.basedir, config["pga"]["path"]),
        ref_dir    = os.path.join(workflow.basedir, config["pga"]["reference_dir"]),
        form       = config["pga"]["form"],
        ir_min     = config["pga"]["ir_min"],
        pidentity  = config["pga"]["pidentity"],
        qcoverage  = config["pga"]["qcoverage"],
        out_dir    = OUTDIR + "/{sample}/pga",
        target_dir = OUTDIR + "/{sample}/pga/target",
        gb_dir     = OUTDIR + "/{sample}/pga/gb",
    log:
        OUTDIR + "/{sample}/logs/pga.log",
    threads:
        config["resources"]["pga"]["threads"]
    resources:
        mem_mb  = config["resources"]["pga"]["mem_mb"],
        runtime = config["resources"]["pga"]["runtime"],
    conda:
        "../envs/pga.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.target_dir} {params.gb_dir} $(dirname {log})

        # PGA writes temp files (_temp1, _temp2, _reference1-4, BLAST DBs)
        # directly into the reference directory.  If multiple samples run in
        # parallel they clobber each other.  Solution: give each sample its
        # OWN copy of the reference directory inside its output tree.
        SAMPLE_REF="{params.out_dir}/ref"
        rm -rf "$SAMPLE_REF"
        cp -r {params.ref_dir} "$SAMPLE_REF"

        # PGA expects a directory of FASTA files.
        # Clean the FASTA header to a simple name (PGA uses it for filenames
        # and BLAST databases; special characters break path handling).
        awk 'NR==1 && /^>/ {{print ">{wildcards.sample}"; next}} {{print}}' \
            {input.fasta} > {params.target_dir}/{wildcards.sample}.fasta

        # PGA needs ABSOLUTE paths — its internal path manipulation creates
        # broken paths (double slashes, garbled DB names) with relative paths.
        ABS_REF=$(cd "$SAMPLE_REF" && pwd)
        ABS_TARGET=$(cd {params.target_dir} && pwd)
        ABS_GB=$(cd {params.gb_dir} && pwd)
        ABS_LOG=$(cd {params.out_dir} && pwd)

        perl {params.pga_dir}/PGA.pl \
            -r "$ABS_REF" \
            -t "$ABS_TARGET" \
            -o "$ABS_GB" \
            -f {params.form} \
            -i {params.ir_min} \
            -p {params.pidentity} \
            -q {params.qcoverage} \
            -l "$ABS_LOG/warning" \
            2>&1 | tee {log}

        # Clean up per-sample reference copy (temp files inside)
        rm -rf "$SAMPLE_REF"

        # Move generated GenBank file to standard name
        found=$(find {params.gb_dir} -maxdepth 1 -name "*.gb" 2>/dev/null | head -1)
        if [ -n "$found" ]; then
            cp "$found" {output.gb}
        else
            touch {output.gb}
        fi
        """
