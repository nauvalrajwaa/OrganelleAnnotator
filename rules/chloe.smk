# rules/chloe.smk – Chloe.jl chloroplast genome annotator

rule chloe_annotate:
    """
    Annotate a chloroplast genome using Chloe.jl.
    Chloe's -o flag sets the output DIRECTORY; the actual GFF file is
    written as <basename>.gff inside that directory.
    """
    input:
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        # HANYA output GFF. File .done dihapus dari sini!
        gff  = OUTDIR + "/{sample}/chloe/{sample}.gff",
    params:
        out_dir    = OUTDIR + "/{sample}/chloe",
        chloe_dir  = config["chloe"]["path"],
        extra      = config["chloe"].get("extra", ""),
    log:
        OUTDIR + "/{sample}/logs/chloe.log",
    threads:
        config["resources"]["chloe"]["threads"]
    resources:
        mem_mb  = config["resources"]["chloe"]["mem_mb"],
        runtime = config["resources"]["chloe"]["runtime"],
    conda:
        "../envs/chloe.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        # 1. Jalankan Chloe
        julia --project={params.chloe_dir} \
            {params.chloe_dir}/chloe.jl \
            annotate \
            {input.fasta} \
            -o {params.out_dir} \
            {params.extra} \
            2>&1 | tee {log}

        # 2. Cari GFF yang dihasilkan dan rename sesuai standar kita
        CHLOE_GFF=$(find {params.out_dir} -maxdepth 1 -name "*.gff" -type f ! -name "$(basename {output.gff})" 2>/dev/null | head -1)
        
        if [ -n "$CHLOE_GFF" ]; then
            mv "$CHLOE_GFF" {output.gff}
        fi
        
        # 3. Validasi ketat: Gagalkan jika file GFF tidak ada atau 0 byte (mencegah error diam-diam)
        if [ ! -s "{output.gff}" ]; then
            echo "Error: Chloe gagal menghasilkan file GFF yang valid (kosong/tidak ada)!" >&2
            exit 1
        fi
        """

rule chloe_gff_to_genbank:
    """
    Convert Chloe's GFF3 annotation + original FASTA → GenBank (.gb) file.
    Menghasilkan file .done HANYA setelah file .gb berhasil dibuat.
    """
    input:
        gff   = OUTDIR + "/{sample}/chloe/{sample}.gff",
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        gb    = OUTDIR + "/{sample}/chloe/{sample}.gb",
        # PENTING: File .done dipindah ke sini agar gbdraw menunggu konversi ini selesai!
        done  = touch(OUTDIR + "/{sample}/chloe/{sample}.done"),
    log:
        OUTDIR + "/{sample}/logs/chloe_to_genbank.log",
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/gff_to_genbank.py"