# =============================================================================
# rules/downstream.smk – Post-annotation downstream analysis
# =============================================================================
# Modules:
#   1. fetch_reference   – Fetch reference genomes from NCBI
#   2. rscu_analysis     – Relative Synonymous Codon Usage
#   3. codon_analysis    – Start/stop codon frequencies
#   4. kaks_analysis     – Pairwise Ka/Ks (MAFFT + KaKs_Calculator)
#   5. composition       – GC content & amino acid composition
#   6. prepare_phylo     – Supermatrix construction (MAFFT)
#   7. phylogeny_tree    – Maximum-likelihood tree (IQ-TREE)
#   8. plot_tree         – Tree visualisation
#   9. genome_map        – Circular genome map (pyGenomeViz)
#  10. synteny_analysis  – Synteny comparison (MUMmer4/nucmer)
#  11. downstream_report – Aggregated downstream HTML report
# =============================================================================


# ── Helper: get a GenBank file for a sample ────────────────────────────────
# Picks the first available GenBank from annotation tools (preference order)
def get_sample_gbk(sample):
    """Find the best GenBank annotation file for genome map / downstream."""
    organelle = get_organelle(sample)
    if organelle == "plastid":
        priority = ["chloe", "pga", "plann", "cpgavas2", "liftoff"]
    else:
        priority = ["mfannot", "mitos", "mitoz", "liftoff"]

    for tool in priority:
        gbk_path = f"{OUTDIR}/{tool}/{sample}/{sample}.gbk"
        if os.path.exists(gbk_path):
            return gbk_path
    # Fallback: any .gbk file
    for tool in ALL_TOOLS:
        gbk_path = f"{OUTDIR}/{tool}/{sample}/{sample}.gbk"
        if os.path.exists(gbk_path):
            return gbk_path
    return ""


def get_sample_cds_fasta(sample):
    """Find or derive a CDS FASTA for downstream analysis."""
    # Check common output paths from annotation tools
    candidates = [
        f"{OUTDIR}/liftoff/{sample}/{sample}_cds.fasta",
        f"{OUTDIR}/mitos/{sample}/{sample}_cds.fasta",
        f"{OUTDIR}/chloe/{sample}/{sample}_cds.fasta",
        f"{OUTDIR}/pga/{sample}/{sample}_cds.fasta",
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    # Fall back to the original FASTA (some analyses can still use it)
    return get_fasta(sample)


# ── 1. Fetch reference genomes from NCBI ─────────────────────────────────

rule fetch_reference:
    """Download reference organelle genomes for comparative analysis."""
    output:
        ref_dir=directory(f"{OUTDIR}/downstream/{{sample}}/references"),
    params:
        species=lambda wc: config["downstream"].get("species_name", ""),
        email=lambda wc: config["downstream"].get("email", ""),
        organelle=lambda wc: get_organelle(wc.sample),
        max_genomes=lambda wc: config["downstream"].get("max_ref_genomes", 10),
        min_len=lambda wc: config["downstream"].get("min_genome_length", 10000),
        max_len=lambda wc: config["downstream"].get("max_genome_length", 300000),
    conda:
        "../envs/downstream.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/fetch_reference.log",
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
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        rscu_tsv=f"{OUTDIR}/downstream/{{sample}}/rscu/rscu.tsv",
        rscu_barplot=f"{OUTDIR}/downstream/{{sample}}/rscu/rscu_barplot.png",
    params:
        output_dir=f"{OUTDIR}/downstream/{{sample}}/rscu",
        genetic_code=lambda wc: get_genetic_code(wc.sample),
    conda:
        "../envs/downstream.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/rscu.log",
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
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        codon_stats=f"{OUTDIR}/downstream/{{sample}}/codons/codon_stats.txt",
    params:
        genetic_code=lambda wc: get_genetic_code(wc.sample),
    conda:
        "../envs/downstream.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/codons.log",
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
        sample_fasta=lambda wc: get_fasta(wc.sample),
        ref_dir=f"{OUTDIR}/downstream/{{sample}}/references",
    output:
        kaks_tsv=f"{OUTDIR}/downstream/{{sample}}/kaks/kaks_summary.tsv",
    params:
        output_dir=f"{OUTDIR}/downstream/{{sample}}/kaks",
        genetic_code=lambda wc: get_genetic_code(wc.sample),
        method=lambda wc: config["downstream"].get("kaks_method", "NG"),
    conda:
        "../envs/downstream.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/kaks.log",
    shell:
        """
        # Use first reference genome for pairwise comparison
        REF_FASTA=$(find {input.ref_dir} -name '*.fasta' -type f | head -1)
        if [ -z "$REF_FASTA" ]; then
            echo "No reference FASTA found" > {log}
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
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        gc_plot=f"{OUTDIR}/downstream/{{sample}}/composition/gc_content_plot.png",
        aa_plot=f"{OUTDIR}/downstream/{{sample}}/composition/aa_composition_plot.png",
    params:
        output_dir=f"{OUTDIR}/downstream/{{sample}}/composition",
        genetic_code=lambda wc: get_genetic_code(wc.sample),
    conda:
        "../envs/downstream.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/composition.log",
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
        sample_fasta=lambda wc: get_fasta(wc.sample),
        ref_dir=f"{OUTDIR}/downstream/{{sample}}/references",
    output:
        supermatrix=f"{OUTDIR}/downstream/{{sample}}/phylogeny/supermatrix.fasta",
        partitions=f"{OUTDIR}/downstream/{{sample}}/phylogeny/partitions.nex",
    params:
        output_dir=f"{OUTDIR}/downstream/{{sample}}/phylogeny",
        min_genes=lambda wc: config["downstream"].get("phylo_min_genes", 4),
    conda:
        "../envs/phylo.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/prepare_phylo.log",
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
        supermatrix=f"{OUTDIR}/downstream/{{sample}}/phylogeny/supermatrix.fasta",
        partitions=f"{OUTDIR}/downstream/{{sample}}/phylogeny/partitions.nex",
    output:
        treefile=f"{OUTDIR}/downstream/{{sample}}/phylogeny/phylogeny.treefile",
    params:
        prefix=f"{OUTDIR}/downstream/{{sample}}/phylogeny/phylogeny",
        model=lambda wc: config["downstream"].get("phylo_model", "GTR+G"),
        bootstrap=lambda wc: config["downstream"].get("phylo_bootstrap", 1000),
    conda:
        "../envs/phylo.yaml"
    threads:
        lambda wc: config["resources"].get("downstream", {}).get("threads", 4)
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/iqtree.log",
    shell:
        """
        # Skip if supermatrix is empty
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
        treefile=f"{OUTDIR}/downstream/{{sample}}/phylogeny/phylogeny.treefile",
    output:
        tree_png=f"{OUTDIR}/downstream/{{sample}}/phylogeny/tree_plot.png",
    conda:
        "../envs/downstream.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/plot_tree.log",
    shell:
        """
        python {workflow.basedir}/scripts/plot_tree.py \
            {input.treefile} {output.tree_png} \
        2>&1 | tee {log}
        """


# ── 9. Genome map (pyGenomeViz, replaces Circos) ─────────────────────────

rule genome_map:
    """Generate circular genome map using pyGenomeViz (no Perl/Circos needed)."""
    input:
        # Use the checkpoint .done marker to ensure annotation has finished
        done=lambda wc: [
            f"{OUTDIR}/{tool}/{wc.sample}/{wc.sample}.done"
            for tool in tools_for_sample(wc.sample)
        ],
    output:
        genome_map_png=f"{OUTDIR}/downstream/{{sample}}/genome_map/genome_map.png",
    params:
        sample=lambda wc: wc.sample,
        outdir=OUTDIR,
    conda:
        "../envs/downstream.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/genome_map.log",
    run:
        # Find GenBank file dynamically
        gbk = get_sample_gbk(params.sample)
        if gbk and os.path.exists(gbk):
            shell(
                "python {workflow.basedir}/scripts/create_genome_map.py "
                f"'{gbk}' {{output.genome_map_png}} "
                "2>&1 | tee {log}"
            )
        else:
            # Create placeholder
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
        sample_fasta=lambda wc: get_fasta(wc.sample),
        ref_dir=f"{OUTDIR}/downstream/{{sample}}/references",
    output:
        synteny_plot=f"{OUTDIR}/downstream/{{sample}}/synteny/synteny_plot.png",
        synteny_stats=f"{OUTDIR}/downstream/{{sample}}/synteny/synteny_stats.tsv",
    conda:
        "../envs/synteny.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/synteny.log",
    shell:
        """
        # Use first reference for synteny
        REF_FASTA=$(find {input.ref_dir} -name '*.fasta' -type f | head -1)
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
        rscu_tsv=f"{OUTDIR}/downstream/{{sample}}/rscu/rscu.tsv",
        rscu_barplot=f"{OUTDIR}/downstream/{{sample}}/rscu/rscu_barplot.png",
        codon_stats=f"{OUTDIR}/downstream/{{sample}}/codons/codon_stats.txt",
        kaks_tsv=f"{OUTDIR}/downstream/{{sample}}/kaks/kaks_summary.tsv",
        gc_plot=f"{OUTDIR}/downstream/{{sample}}/composition/gc_content_plot.png",
        aa_plot=f"{OUTDIR}/downstream/{{sample}}/composition/aa_composition_plot.png",
        tree_png=f"{OUTDIR}/downstream/{{sample}}/phylogeny/tree_plot.png",
        genome_map_png=f"{OUTDIR}/downstream/{{sample}}/genome_map/genome_map.png",
        synteny_plot=f"{OUTDIR}/downstream/{{sample}}/synteny/synteny_plot.png",
        synteny_stats=f"{OUTDIR}/downstream/{{sample}}/synteny/synteny_stats.tsv",
    output:
        html=f"{OUTDIR}/downstream/{{sample}}/downstream_report.html",
    params:
        sample=lambda wc: wc.sample,
        species=lambda wc: config["downstream"].get("species_name", wc.sample),
    conda:
        "../envs/downstream.yaml"
    log:
        f"{OUTDIR}/logs/downstream/{{sample}}/report.log",
    script:
        "../scripts/generate_downstream_report.py"
