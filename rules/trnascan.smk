# rules/trnascan.smk – tRNAscan-SE: tRNA gene detection

rule trnascan_annotate:
    """
    Detect tRNA genes using tRNAscan-SE.
    Produces tabular output, secondary structure file, and GFF3.
    """
    input:
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        done = touch(OUTDIR + "/{sample}/trnascan/{sample}.done"),
        tsv  = OUTDIR + "/{sample}/trnascan/{sample}.trnascan.tsv",
        gff  = OUTDIR + "/{sample}/trnascan/{sample}.gff",
        ss   = OUTDIR + "/{sample}/trnascan/{sample}.ss",
    params:
        out_dir      = OUTDIR + "/{sample}/trnascan",
        search_mode  = config["trnascan"]["search_mode"],
        score_cutoff = config["trnascan"]["score_cutoff"],
        extra        = config["trnascan"].get("extra", ""),
    log:
        OUTDIR + "/{sample}/logs/trnascan.log",
    threads:
        config["resources"]["trnascan"]["threads"]
    resources:
        mem_mb  = config["resources"]["trnascan"]["mem_mb"],
        runtime = config["resources"]["trnascan"]["runtime"],
    conda:
        "../envs/trnascan.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

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

        touch {output.tsv} {output.gff} {output.ss}
        """
