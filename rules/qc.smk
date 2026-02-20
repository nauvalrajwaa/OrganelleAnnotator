# =============================================================================
# rules/qc.smk – Quality control: BUSCO + custom gene completeness
# =============================================================================

# ---------------------------------------------------------------------------
# BUSCO genome completeness assessment
# ---------------------------------------------------------------------------
rule busco:
    """
    Run BUSCO on the input FASTA to assess genome completeness.
    Works as a general QC step for organelle genomes.
    """
    input:
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        summary=f"{OUTDIR}/qc/busco/{{sample}}/short_summary.txt",
        out_dir=directory(f"{OUTDIR}/qc/busco/{{sample}}"),
    params:
        lineage=config["qc"]["busco_lineage"],
        mode=config["qc"]["busco_mode"],
    log:
        f"{OUTDIR}/logs/qc/busco_{{sample}}.log",
    threads:
        config["resources"]["busco"]["threads"]
    resources:
        mem_mb=config["resources"]["busco"]["mem_mb"],
        runtime=config["resources"]["busco"]["runtime"],
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
            --out_path {OUTDIR}/qc/busco/ \
            -c {threads} \
            --offline \
            2>&1 | tee {log} || true

        # BUSCO writes into a nested dir; copy summary to expected location
        BUSCO_DIR="{OUTDIR}/qc/busco/{wildcards.sample}"
        SUMMARY=$(find "$BUSCO_DIR" -name "short_summary*.txt" 2>/dev/null | head -1)
        if [ -n "$SUMMARY" ]; then
            cp "$SUMMARY" {output.summary}
        else
            echo "BUSCO did not produce a summary (may be expected for small organelle genomes)" > {output.summary}
        fi
        """


# ---------------------------------------------------------------------------
# Custom gene completeness summary per sample
# ---------------------------------------------------------------------------
rule gene_completeness_summary:
    """
    Parse annotation outputs from all tools run for a sample and produce
    a unified gene-completeness QC summary TSV.
    Reports: tool, gene_count, gene_names, tRNA_count, rRNA_count.
    """
    input:
        done=lambda wc: [
            f"{OUTDIR}/{tool}/{wc.sample}/{wc.sample}.done"
            for tool in tools_for_sample(wc.sample)
        ],
    output:
        summary=f"{OUTDIR}/qc/summary/{{sample}}.qc_summary.tsv",
    params:
        tools=lambda wc: tools_for_sample(wc.sample),
        out_base=OUTDIR,
    log:
        f"{OUTDIR}/logs/qc/gene_summary_{{sample}}.log",
    run:
        import re, os, csv

        rows = []
        sample = wildcards.sample

        for tool in params.tools:
            gene_names = []
            trna_count = 0
            rrna_count = 0
            gene_count = 0

            tool_dir = os.path.join(params.out_base, tool, sample)

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
                            # Parse gene annotations from masterfile
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

            elif tool == "mitos":
                # Parse MITOS BED or GFF output
                bed = os.path.join(tool_dir, "result.bed")
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
                                    name_m = re.search(r"gene_id=([^;]+)", attrs)
                                name = name_m.group(1) if name_m else "unknown"
                                if ftype in ("gene", "CDS", "tRNA", "rRNA"):
                                    if name not in gene_names:
                                        gene_count += 1
                                        gene_names.append(name)
                                    if ftype == "tRNA" or "trn" in name.lower():
                                        trna_count += 1
                                    if ftype == "rRNA" or "rrn" in name.lower():
                                        rrna_count += 1
                elif os.path.exists(bed) and os.path.getsize(bed) > 0:
                    with open(bed) as f:
                        for line in f:
                            cols = line.strip().split("\t")
                            if len(cols) >= 4:
                                name = cols[3]
                                if name not in gene_names:
                                    gene_count += 1
                                    gene_names.append(name)
                                if "trn" in name.lower():
                                    trna_count += 1
                                if "rrn" in name.lower():
                                    rrna_count += 1

            elif tool == "mitoz":
                # Parse MitoZ GFF output
                gff = os.path.join(tool_dir, f"{sample}.gff")
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
                                name = name_m.group(1) if name_m else "unknown"
                                if ftype in ("gene", "CDS", "tRNA", "rRNA"):
                                    if name not in gene_names:
                                        gene_count += 1
                                        gene_names.append(name)
                                    if ftype == "tRNA" or "trn" in name.lower():
                                        trna_count += 1
                                    if ftype == "rRNA" or "rrn" in name.lower():
                                        rrna_count += 1

            elif tool == "trnascan":
                # Parse tRNAscan-SE GFF output
                gff = os.path.join(tool_dir, f"{sample}.gff")
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
                                name = name_m.group(1) if name_m else "unknown"
                                if ftype in ("tRNA", "gene"):
                                    if name not in gene_names:
                                        gene_count += 1
                                        gene_names.append(name)
                                    trna_count += 1

            elif tool == "aragorn":
                # Parse Aragorn GFF output
                gff = os.path.join(tool_dir, f"{sample}.gff")
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
                                name = name_m.group(1) if name_m else "unknown"
                                if name not in gene_names:
                                    gene_count += 1
                                    gene_names.append(name)
                                if "trn" in name.lower() or ftype == "tRNA":
                                    trna_count += 1
                                if "tmrna" in name.lower() or ftype == "tmRNA":
                                    pass  # tmRNA counted as gene, not tRNA/rRNA

            rows.append({
                "sample": sample,
                "tool": tool,
                "gene_count": gene_count,
                "trna_count": trna_count,
                "rrna_count": rrna_count,
                "gene_names": ";".join(gene_names),
            })

        with open(output.summary, "w", newline="") as out:
            writer = csv.DictWriter(
                out,
                fieldnames=["sample", "tool", "gene_count", "trna_count", "rrna_count", "gene_names"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)
