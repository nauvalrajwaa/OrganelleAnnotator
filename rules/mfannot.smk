# rules/mfannot.smk – MFannot via Docker (nbeck/mfannot)

rule mfannot_annotate:
    """
    Annotate a mitochondrial (or plastid) genome using MFannot through Docker.
    Docker image: nbeck/mfannot (https://hub.docker.com/r/nbeck/mfannot/)
    """
    input:
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        done       = touch(OUTDIR + "/{sample}/mfannot/{sample}.done"),
        masterfile = OUTDIR + "/{sample}/mfannot/{sample}.new",
    params:
        docker_image    = config["mfannot"]["docker_image"],
        genetic_code    = lambda wc: samples_df.loc[wc.sample, "genetic_code"],
        blast_evalue    = config["mfannot"]["blast_evalue"],
        min_orf_len     = config["mfannot"]["min_orf_len"],
        max_intron_size = config["mfannot"]["max_intron_size"],
        extra           = config["mfannot"]["extra"],
        sqn_flag        = lambda wc: "--sqnformat" if config["mfannot"]["sqn_format"] else "",
        out_dir         = OUTDIR + "/{sample}/mfannot",
        use_singularity = config["mfannot"].get("use_singularity", False),
    log:
        OUTDIR + "/{sample}/logs/mfannot.log",
    threads:
        config["resources"]["mfannot"]["threads"]
    resources:
        mem_mb  = config["resources"]["mfannot"]["mem_mb"],
        runtime = config["resources"]["mfannot"]["runtime"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        INPUT_DIR=$(cd "$(dirname {input.fasta})" && pwd)
        INPUT_BASENAME=$(basename {input.fasta})
        OUTPUT_DIR=$(cd {params.out_dir} 2>/dev/null || mkdir -p {params.out_dir} && cd {params.out_dir} && pwd)

        cp {input.fasta} {params.out_dir}/{wildcards.sample}.fasta

        if [ "{params.use_singularity}" = "True" ]; then
            singularity exec \
                --bind "${{OUTPUT_DIR}}":/data \
                docker://{params.docker_image} \
                mfannot \
                -g {params.genetic_code} \
                --blast {params.blast_evalue} \
                --minorflen {params.min_orf_len} \
                --maxintronsize {params.max_intron_size} \
                {params.sqn_flag} \
                {params.extra} \
                -o /data/{wildcards.sample}.new \
                -l /data/{wildcards.sample}.log \
                /data/{wildcards.sample}.fasta \
                2>&1 | tee {log}
        else
            docker run --rm \
                -v "${{OUTPUT_DIR}}":/data \
                {params.docker_image} \
                mfannot \
                -g {params.genetic_code} \
                --blast {params.blast_evalue} \
                --minorflen {params.min_orf_len} \
                --maxintronsize {params.max_intron_size} \
                {params.sqn_flag} \
                {params.extra} \
                -o /data/{wildcards.sample}.new \
                -l /data/{wildcards.sample}.log \
                /data/{wildcards.sample}.fasta \
                2>&1 | tee {log}
        fi

        touch {output.masterfile}
        """


rule mfannot_to_gff:
    """
    Convert MFannot masterfile output to a simplified GFF3 for downstream QC.
    """
    input:
        masterfile = OUTDIR + "/{sample}/mfannot/{sample}.new",
    output:
        gff = OUTDIR + "/{sample}/mfannot/{sample}.gff",
    log:
        OUTDIR + "/{sample}/logs/mfannot_to_gff.log",
    run:
        import re

        gff_lines = ["##gff-version 3"]
        current_seq = None

        with open(input.masterfile) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    parts = line.split()
                    current_seq = parts[0].lstrip(">").rstrip(";")
                elif line.startswith(";;") and "gene" in line.lower():
                    pass
                m = re.match(
                    r";\s+(\w+)\s+(\d+)-(\d+)\s+(\w+)",
                    line
                )
                if m:
                    feature, start, end, name = m.groups()
                    strand = "+"
                    gff_lines.append(
                        f"{current_seq or wildcards.sample}\tMFannot\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={name};Name={name}"
                    )

        with open(output.gff, "w") as out:
            out.write("\n".join(gff_lines) + "\n")
