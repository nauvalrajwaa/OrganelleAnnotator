# rules/qc.smk – Quality control: BUSCO + custom gene completeness

# ---------------------------------------------------------------------------
# BUSCO genome completeness assessment
# ---------------------------------------------------------------------------
rule busco:
    """
    Run BUSCO on the input FASTA to assess genome completeness.
    """
    input:
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        summary = OUTDIR + "/{sample}/qc/busco/short_summary.txt",
        out_dir = directory(OUTDIR + "/{sample}/qc/busco"),
    params:
        lineage = config["qc"]["busco_lineage"],
        mode    = config["qc"]["busco_mode"],
    log:
        OUTDIR + "/{sample}/logs/busco.log",
    threads:
        config["resources"]["busco"]["threads"]
    resources:
        mem_mb  = config["resources"]["busco"]["mem_mb"],
        runtime = config["resources"]["busco"]["runtime"],
    conda:
        "../envs/busco.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})

        busco \
            -i {input.fasta} \
            -l {params.lineage} \
            -m {params.mode} \
            -o {wildcards.sample} \
            --out_path {output.out_dir}/ \
            -c {threads} \
            --offline \
            2>&1 | tee {log} || true

        BUSCO_DIR="{output.out_dir}/{wildcards.sample}"
        SUMMARY=$(find "$BUSCO_DIR" -name "short_summary*.txt" 2>/dev/null | head -1)
        if [ -z "$SUMMARY" ]; then
            SUMMARY=$(find "{output.out_dir}" -name "short_summary*.txt" 2>/dev/null | head -1)
        fi
        if [ -n "$SUMMARY" ]; then
            cp "$SUMMARY" {output.summary}
        else
            echo "BUSCO did not produce a summary" > {output.summary}
        fi
        """


# ---------------------------------------------------------------------------
# Custom gene completeness summary per sample
# ---------------------------------------------------------------------------
rule gene_completeness_summary:
    """
    Parse annotation outputs from all tools run for a sample and produce
    a unified gene-completeness QC summary TSV.
    """
    input:
        done = lambda wc: [
            OUTDIR + "/" + wc.sample + "/" + tool + "/" + wc.sample + ".done"
            for tool in tools_for_sample(wc.sample)
        ],
    output:
        summary = OUTDIR + "/{sample}/qc/qc_summary.tsv",
    params:
        tools    = lambda wc: tools_for_sample(wc.sample),
        out_base = OUTDIR,
    log:
        OUTDIR + "/{sample}/logs/gene_summary.log",
    run:
        import re, os, csv

        rows = []
        sample = wildcards.sample

        for tool in params.tools:
            gene_names = []
            trna_count = 0
            rrna_count = 0
            gene_count = 0

            tool_dir = os.path.join(params.out_base, sample, tool)

            if tool == "chloe":
                gff = os.path.join(tool_dir, f"{sample}.gff")
                if os.path.exists(gff):
                    with open(gff) as f:
                        for line in f:
                            if line.startswith("#"):
                                continue
                            cols = line.strip().split("\t")
                            if len(cols) >= 9:
                                ftype = cols[2]
                                attrs = cols[8]
                                name_m = re.search(r"Name=([^;]+)", attrs)
                                name = name_m.group(1) if name_m else "unknown"
                                if ftype == "gene":
                                    gene_count += 1
                                    gene_names.append(name)
                                if "trn" in name.lower():
                                    trna_count += 1
                                if "rrn" in name.lower():
                                    rrna_count += 1

            elif tool == "pga":
                gb = os.path.join(tool_dir, f"{sample}.gb")
                if os.path.exists(gb):
                    with open(gb) as f:
                        for line in f:
                            m = re.match(r'\s+/gene="([^"]+)"', line)
                            if m:
                                name = m.group(1)
                                if name not in gene_names:
                                    gene_names.append(name)
                                    gene_count += 1
                                    if "trn" in name.lower():
                                        trna_count += 1
                                    if "rrn" in name.lower():
                                        rrna_count += 1

            elif tool == "mfannot":
                masterfile = os.path.join(tool_dir, f"{sample}.new")
                if os.path.exists(masterfile):
                    with open(masterfile) as f:
                        for line in f:
                            m = re.search(r"gene\s*=\s*(\S+)", line, re.IGNORECASE)
                            if not m:
                                m = re.match(r";\s+(\w+)\s+\d+-\d+", line)
                            if m:
                                name = m.group(1)
                                if name not in gene_names:
                                    gene_names.append(name)
                                    gene_count += 1
                                    if "trn" in name.lower():
                                        trna_count += 1
                                    if "rrn" in name.lower():
                                        rrna_count += 1

            elif tool in ("mitos", "mitoz", "trnascan", "aragorn", "liftoff"):
                gff = os.path.join(tool_dir, f"{sample}.gff")
                if tool == "mitos":
                    gff = os.path.join(tool_dir, "result.gff")
                if os.path.exists(gff) and os.path.getsize(gff) > 0:
                    with open(gff) as f:
                        for line in f:
                            if line.startswith("#"):
                                continue
                            cols = line.strip().split("\t")
                            if len(cols) >= 9:
                                ftype = cols[2]
                                attrs = cols[8]
                                name_m = re.search(r"Name=([^;]+)", attrs)
                                if not name_m:
                                    name_m = re.search(r"gene=([^;]+)", attrs)
                                if not name_m:
                                    name_m = re.search(r"ID=([^;]+)", attrs)
                                name = name_m.group(1) if name_m else "unknown"
                                if ftype in ("gene", "CDS", "tRNA", "rRNA", "mRNA"):
                                    if name not in gene_names:
                                        gene_count += 1
                                        gene_names.append(name)
                                    if ftype == "tRNA" or "trn" in name.lower():
                                        trna_count += 1
                                    if ftype == "rRNA" or "rrn" in name.lower():
                                        rrna_count += 1

            elif tool == "fpma":
                tsv = os.path.join(tool_dir, f"{sample}.presence.tsv")
                if os.path.exists(tsv):
                    with open(tsv) as f:
                        for line in f:
                            cols = line.strip().split("\t")
                            if len(cols) >= 2:
                                name = cols[0]
                                present = cols[1] if len(cols) > 1 else ""
                                if present.strip().lower() in ("true", "1", "yes", "+"):
                                    gene_count += 1
                                    gene_names.append(name)
                                    if "trn" in name.lower():
                                        trna_count += 1
                                    if "rrn" in name.lower():
                                        rrna_count += 1

            elif tool == "gbdraw":
                pass  # Visualisation only

            rows.append({
                "sample": sample,
                "tool": tool,
                "gene_count": gene_count,
                "trna_count": trna_count,
                "rrna_count": rrna_count,
                "gene_names": ";".join(gene_names),
            })

        os.makedirs(os.path.dirname(output.summary), exist_ok=True)
        with open(output.summary, "w", newline="") as out:
            writer = csv.DictWriter(
                out,
                fieldnames=["sample", "tool", "gene_count", "trna_count", "rrna_count", "gene_names"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)
