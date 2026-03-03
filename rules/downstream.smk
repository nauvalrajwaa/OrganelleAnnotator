# =============================================================================
# rules/downstream.smk – Post-annotation downstream analysis
# =============================================================================

# ── Helper: get a GenBank file for a sample ────────────────────────────────
def get_sample_gbk(sample):
    """Find the best GenBank annotation file for genome map / downstream."""
    organelle = samples_df.loc[sample, "organelle"]
    if organelle == "plastid":
        priority = ["chloe", "pga", "liftoff"]
    else:
        priority = ["mfannot", "mitos", "mitoz", "liftoff"]

    # Coba cari dengan ekstensi .gbk maupun .gb
    extensions = [".gbk", ".gb"]

    for tool in priority:
        for ext in extensions:
            gbk_path = os.path.join(OUTDIR, sample, tool, sample + ext)
            if os.path.exists(gbk_path):
                return gbk_path
                
    for tool in ALL_TOOLS:
        for ext in extensions:
            gbk_path = os.path.join(OUTDIR, sample, tool, sample + ext)
            if os.path.exists(gbk_path):
                return gbk_path
                
    return ""


# ── RULE BARU: Ekstrak CDS langsung dari file GBK hasil anotasi ──────
rule extract_best_cds:
    """Extract CDS sequences from the best annotation GenBank file."""
    input:
        done = lambda wc: [
            OUTDIR + "/" + wc.sample + "/" + tool + "/" + wc.sample + ".done"
            for tool in tools_for_sample(wc.sample)
        ],
    output:
        best_cds = OUTDIR + "/{sample}/downstream/sample_best_cds.fasta"
    params:
        sample = lambda wc: wc.sample,
        fallback_fasta = lambda wc: samples_df.loc[wc.sample, "fasta"]
    run:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        import shutil
        import os

        # Ambil path file GBK terbaik menggunakan fungsi helper kita
        gbk_path = get_sample_gbk(params.sample)
        cds_records = []

        if gbk_path and os.path.exists(gbk_path):
            try:
                # SeqIO bisa membaca file .gb sama seperti .gbk dengan format "genbank"
                record = SeqIO.read(gbk_path, "genbank")
                for feature in record.features:
                    if feature.type == "CDS":
                        # Ambil nama gen
                        gene_name = feature.qualifiers.get("gene", feature.qualifiers.get("locus_tag", ["unknown_gene"]))[0]
                        try:
                            # Ekstrak sekuens dari GBK
                            seq = feature.extract(record.seq)
                            cds_rec = SeqRecord(seq, id=gene_name, description="")
                            cds_records.append(cds_rec)
                        except Exception as e:
                            pass
            except Exception as e:
                print(f"Error parsing {gbk_path}: {e}")

        # Simpan file CDS
        if cds_records:
            SeqIO.write(cds_records, output.best_cds, "fasta")
            print(f"Berhasil mengekstrak {len(cds_records)} gen CDS dari {gbk_path}")
        else:
            print(f"WARNING: Tidak ada gen CDS di {gbk_path}. Memakai genom utuh.")
            shutil.copy(params.fallback_fasta, output.best_cds)


# ── 1. Fetch reference genomes from NCBI ─────────────────────────────────

rule fetch_reference:
    """Download reference organelle genomes for comparative analysis."""
    output:
        ref_dir = directory(OUTDIR + "/{sample}/downstream/references"),
    params:
        species    = lambda wc: config["downstream"].get("species_name", ""),
        email      = lambda wc: config["downstream"].get("email", ""),
        organelle  = lambda wc: samples_df.loc[wc.sample, "organelle"],
        max_genomes = lambda wc: config["downstream"].get("max_ref_genomes", 10),
        min_len    = lambda wc: config["downstream"].get("min_genome_length", 10000),
        max_len    = lambda wc: config["downstream"].get("max_genome_length", 300000),
    conda:
        "../envs/downstream.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_fetch_reference.log",
    shell:
        """
        python {workflow.basedir}/scripts/fetch_organelle_ref.py \
            "{params.species}" "{params.email}" {output.ref_dir} \
            --organelle {params.organelle} \
            --max_genomes {params.max_genomes} \
            --min_len {params.min_len} \
            --max_len {params.max_len} \
        2>&1 | tee {log}
        """


# ── 2. RSCU analysis ─────────────────────────────────────────────────────

rule rscu_analysis:
    """Compute Relative Synonymous Codon Usage from CDS sequences."""
    input:
        fasta = OUTDIR + "/{sample}/downstream/sample_best_cds.fasta",
    output:
        rscu_tsv     = OUTDIR + "/{sample}/downstream/rscu/rscu.tsv",
        rscu_barplot = OUTDIR + "/{sample}/downstream/rscu/rscu_barplot.png",
    params:
        output_dir   = OUTDIR + "/{sample}/downstream/rscu",
        genetic_code = lambda wc: samples_df.loc[wc.sample, "genetic_code"],
    conda:
        "../envs/downstream.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_rscu.log",
    shell:
        """
        python {workflow.basedir}/scripts/calculate_rscu.py \
            {input.fasta} {params.output_dir} \
            --genetic_code {params.genetic_code} \
        2>&1 | tee {log}
        """


# ── 3. Start/stop codon analysis ─────────────────────────────────────────

rule codon_analysis:
    """Analyse start and stop codon usage across CDS features."""
    input:
        fasta = OUTDIR + "/{sample}/downstream/sample_best_cds.fasta",
    output:
        codon_stats = OUTDIR + "/{sample}/downstream/codons/codon_stats.txt",
    params:
        genetic_code = lambda wc: samples_df.loc[wc.sample, "genetic_code"],
    conda:
        "../envs/downstream.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_codons.log",
    shell:
        """
        python {workflow.basedir}/scripts/analyze_codons.py \
            {input.fasta} {output.codon_stats} \
            --genetic_code {params.genetic_code} \
        2>&1 | tee {log}
        """


# ── 4. Ka/Ks analysis ────────────────────────────────────────────────────

rule kaks_analysis:
    """Pairwise Ka/Ks estimation using MAFFT + KaKs_Calculator."""
    input:
        sample_fasta = OUTDIR + "/{sample}/downstream/sample_best_cds.fasta",
        ref_dir      = OUTDIR + "/{sample}/downstream/references",
    output:
        kaks_tsv = OUTDIR + "/{sample}/downstream/kaks/kaks_summary.tsv",
    params:
        output_dir   = OUTDIR + "/{sample}/downstream/kaks",
        genetic_code = lambda wc: samples_df.loc[wc.sample, "genetic_code"],
        method       = lambda wc: config["downstream"].get("kaks_method", "NG"),
    conda:
        "../envs/downstream.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_kaks.log",
    shell:
        """
        REF_FASTA=$(find {input.ref_dir} -name '*_cds.fasta' -type f | head -1)
        if [ -z "$REF_FASTA" ]; then
            echo "No reference CDS FASTA found" > {log}
            echo -e "Gene\\tKa\\tKs\\tKa/Ks\\tMethod" > {output.kaks_tsv}
        else
            python {workflow.basedir}/scripts/run_kaks_analysis.py \
                {input.sample_fasta} "$REF_FASTA" {params.output_dir} \
                --genetic_code {params.genetic_code} \
                --method {params.method} \
            2>&1 | tee {log}
        fi
        """


# ── 5. GC content & amino acid composition ───────────────────────────────

rule composition_analysis:
    """Analyse GC content and amino acid composition per CDS gene."""
    input:
        fasta = OUTDIR + "/{sample}/downstream/sample_best_cds.fasta",
    output:
        gc_plot = OUTDIR + "/{sample}/downstream/composition/gc_content_plot.png",
        aa_plot = OUTDIR + "/{sample}/downstream/composition/aa_composition_plot.png",
    params:
        output_dir   = OUTDIR + "/{sample}/downstream/composition",
        genetic_code = lambda wc: samples_df.loc[wc.sample, "genetic_code"],
    conda:
        "../envs/downstream.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_composition.log",
    shell:
        """
        python {workflow.basedir}/scripts/analyze_composition.py \
            {input.fasta} {params.output_dir} \
            --genetic_code {params.genetic_code} \
        2>&1 | tee {log}
        """


# ── 6. Prepare phylogenetic supermatrix ──────────────────────────────────

rule prepare_phylo:
    """Build supermatrix from shared genes using MAFFT alignment."""
    input:
        sample_fasta = OUTDIR + "/{sample}/downstream/sample_best_cds.fasta",
        ref_dir      = OUTDIR + "/{sample}/downstream/references",
    output:
        supermatrix = OUTDIR + "/{sample}/downstream/phylogeny/supermatrix.fasta",
        partitions  = OUTDIR + "/{sample}/downstream/phylogeny/partitions.nex",
    params:
        output_dir = OUTDIR + "/{sample}/downstream/phylogeny",
        min_genes  = lambda wc: config["downstream"].get("phylo_min_genes", 4),
    conda:
        "../envs/phylo.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_prepare_phylo.log",
    shell:
        """
        python {workflow.basedir}/scripts/prepare_phylo.py \
            {input.sample_fasta} {input.ref_dir} {params.output_dir} \
            --min_genes {params.min_genes} \
        2>&1 | tee {log}
        """


# ── 7. IQ-TREE phylogeny ─────────────────────────────────────────────────

rule phylogeny_tree:
    """Run IQ-TREE maximum-likelihood phylogeny on supermatrix."""
    input:
        supermatrix = OUTDIR + "/{sample}/downstream/phylogeny/supermatrix.fasta",
        partitions  = OUTDIR + "/{sample}/downstream/phylogeny/partitions.nex",
    output:
        treefile = OUTDIR + "/{sample}/downstream/phylogeny/phylogeny.treefile",
    params:
        prefix    = OUTDIR + "/{sample}/downstream/phylogeny/phylogeny",
        model     = lambda wc: config["downstream"].get("phylo_model", "GTR+G"),
        bootstrap = lambda wc: config["downstream"].get("phylo_bootstrap", 1000),
    conda:
        "../envs/phylo.yaml"
    threads:
        lambda wc: config["resources"].get("downstream", {}).get("threads", 4)
    log:
        OUTDIR + "/{sample}/logs/downstream_iqtree.log",
    shell:
        """
        if [ ! -s {input.supermatrix} ]; then
            echo "Empty supermatrix; skipping IQ-TREE" > {log}
            touch {output.treefile}
            exit 0
        fi
        iqtree2 -s {input.supermatrix} \
            -p {input.partitions} \
            -m {params.model} \
            -B {params.bootstrap} \
            --prefix {params.prefix} \
            -T {threads} \
            --redo \
        2>&1 | tee {log} || \
        iqtree -s {input.supermatrix} \
            -p {input.partitions} \
            -m {params.model} \
            -B {params.bootstrap} \
            --prefix {params.prefix} \
            -T {threads} \
            --redo \
        2>&1 | tee -a {log}
        """


# ── 8. Plot phylogenetic tree ────────────────────────────────────────────

rule plot_tree:
    """Render phylogenetic tree as a PNG image."""
    input:
        treefile = OUTDIR + "/{sample}/downstream/phylogeny/phylogeny.treefile",
    output:
        tree_png = OUTDIR + "/{sample}/downstream/phylogeny/tree_plot.png",
    conda:
        "../envs/downstream.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_plot_tree.log",
    shell:
        """
        python {workflow.basedir}/scripts/plot_tree.py \
            {input.treefile} {output.tree_png} \
        2>&1 | tee {log}
        """


# ── 9. Genome map (pyGenomeViz) ──────────────────────────────────────────

rule genome_map:
    """Generate circular genome map using pyGenomeViz."""
    input:
        done = lambda wc: [
            OUTDIR + "/" + wc.sample + "/" + tool + "/" + wc.sample + ".done"
            for tool in tools_for_sample(wc.sample)
        ],
    output:
        genome_map_png = OUTDIR + "/{sample}/downstream/genome_map/genome_map.png",
    params:
        sample = lambda wc: wc.sample,
        outdir = OUTDIR,
    conda:
        "../envs/downstream.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_genome_map.log",
    run:
        gbk = get_sample_gbk(params.sample)
        if gbk and os.path.exists(gbk):
            shell(
                "python {workflow.basedir}/scripts/create_genome_map.py "
                f"'{gbk}' {{output.genome_map_png}} "
                "2>&1 | tee {log}"
            )
        else:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(8, 3))
            ax.text(0.5, 0.5, "No GenBank file available for genome map",
                    ha="center", va="center", fontsize=12, color="grey")
            ax.axis("off")
            os.makedirs(os.path.dirname(output.genome_map_png), exist_ok=True)
            plt.savefig(output.genome_map_png, dpi=150)
            plt.close()


# ── 10. Synteny analysis (MUMmer4) ───────────────────────────────────────

rule synteny_analysis:
    """Compare genome structure against a reference using nucmer."""
    input:
        sample_fasta = lambda wc: samples_df.loc[wc.sample, "fasta"], # Tetap pakai genom utuh
        ref_dir      = OUTDIR + "/{sample}/downstream/references",
    output:
        synteny_plot  = OUTDIR + "/{sample}/downstream/synteny/synteny_plot.png",
        synteny_stats = OUTDIR + "/{sample}/downstream/synteny/synteny_stats.tsv",
    conda:
        "../envs/synteny.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_synteny.log",
    shell:
        """
        # Ambil file genom referensi utuh, BUKAN file CDS
        REF_FASTA=$(find {input.ref_dir} -name '*.fasta' ! -name '*_cds.fasta' -type f | head -1)
        if [ -z "$REF_FASTA" ]; then
            echo "No reference FASTA found" > {log}
            python {workflow.basedir}/scripts/run_synteny_analysis.py \
                {input.sample_fasta} {input.sample_fasta} \
                {output.synteny_plot} {output.synteny_stats}
        else
            python {workflow.basedir}/scripts/run_synteny_analysis.py \
                {input.sample_fasta} "$REF_FASTA" \
                {output.synteny_plot} {output.synteny_stats} \
            2>&1 | tee {log}
        fi
        """


# ── 11. Downstream HTML report ───────────────────────────────────────────

rule downstream_report:
    """Generate a comprehensive downstream analysis HTML report per sample."""
    input:
        rscu_tsv       = OUTDIR + "/{sample}/downstream/rscu/rscu.tsv",
        rscu_barplot   = OUTDIR + "/{sample}/downstream/rscu/rscu_barplot.png",
        codon_stats    = OUTDIR + "/{sample}/downstream/codons/codon_stats.txt",
        kaks_tsv       = OUTDIR + "/{sample}/downstream/kaks/kaks_summary.tsv",
        gc_plot        = OUTDIR + "/{sample}/downstream/composition/gc_content_plot.png",
        aa_plot        = OUTDIR + "/{sample}/downstream/composition/aa_composition_plot.png",
        tree_png       = OUTDIR + "/{sample}/downstream/phylogeny/tree_plot.png",
        genome_map_png = OUTDIR + "/{sample}/downstream/genome_map/genome_map.png",
        synteny_plot   = OUTDIR + "/{sample}/downstream/synteny/synteny_plot.png",
        synteny_stats  = OUTDIR + "/{sample}/downstream/synteny/synteny_stats.tsv",
    output:
        html = OUTDIR + "/{sample}/downstream/downstream_report.html",
    params:
        sample  = lambda wc: wc.sample,
        species = lambda wc: config["downstream"].get("species_name", wc.sample),
    conda:
        "../envs/downstream.yaml"
    log:
        OUTDIR + "/{sample}/logs/downstream_report.log",
    script:
        "../scripts/generate_downstream_report.py"
