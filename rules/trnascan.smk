# =============================================================================
# rules/trnascan.smk – tRNAscan-SE: tRNA gene detection
# =============================================================================
# tRNAscan-SE is the gold-standard tool for tRNA detection. It works with
# both organelle types (plastid and mitochondrial) and can use organelle-
# specific models via the -O (organelle/mitochondrial) flag.
# Conda: bioconda::trnascan-se
# =============================================================================

rule trnascan_annotate:
    """
    Detect tRNA genes using tRNAscan-SE.
    Produces tabular output, secondary structure file, and GFF3.
    """
    input:
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        done=touch(f"{OUTDIR}/{{sample}}/trnascan/{{sample}}.done"),
        tsv=f"{OUTDIR}/{{sample}}/trnascan/{{sample}}.trnascan.tsv",
        gff=f"{OUTDIR}/{{sample}}/trnascan/{{sample}}.gff",
        ss=f"{OUTDIR}/{{sample}}/trnascan/{{sample}}.ss",
    params:
        out_dir=lambda wc: f"{OUTDIR}/{wc.sample}/trnascan",
        search_mode=lambda wc: config["trnascan"]["search_mode"],
        score_cutoff=config["trnascan"]["score_cutoff"],
        extra=config["trnascan"].get("extra", ""),
    log:
        f"{OUTDIR}/{{sample}}/logs/trnascan.log",
    threads:
        config["resources"]["trnascan"]["threads"]
    resources:
        mem_mb=config["resources"]["trnascan"]["mem_mb"],
        runtime=config["resources"]["trnascan"]["runtime"],
    conda:
        "../envs/trnascan.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        # Determine search mode flag
        # -O = organellar/mitochondrial model
        # -M mam = mammalian mito model
        # --general = general (default)
        MODE_FLAG=""
        case "{params.search_mode}" in
            organelle|mito)
                MODE_FLAG="-O"
                ;;
            mammalian_mito)
                MODE_FLAG="-M mam"
                ;;
            general)
                MODE_FLAG=""
                ;;
        esac

        tRNAscan-SE \
            $MODE_FLAG \
            --thread {threads} \
            -X {params.score_cutoff} \
            -o {output.tsv} \
            -f {output.ss} \
            --gff {output.gff} \
            {params.extra} \
            {input.fasta} \
            2>&1 | tee {log}

        # Ensure outputs exist even if empty
        touch {output.tsv} {output.gff} {output.ss}
        """
