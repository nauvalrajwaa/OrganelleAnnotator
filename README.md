# Organelle Annotation Pipeline

A **Snakemake** workflow for comprehensive organelle genome annotation using **12 tools**, unified QC, an aggregated HTML report, and **integrated downstream analysis** (RSCU, Ka/Ks, phylogeny, genome maps, synteny). Supports both **chloroplast/plastid** and **mitochondrial** genomes with automatic tool selection based on organelle type.

---

## Table of Contents

- [Features](#features)
- [Tool Overview](#tool-overview)
- [Downstream Analysis](#downstream-analysis)
- [Quick Start](#quick-start)
- [Directory Structure](#directory-structure)
- [Configuration](#configuration)
  - [Run Modes](#run-modes)
  - [Sample Sheet](#sample-sheet)
  - [Tool-specific Settings](#tool-specific-settings)
  - [Downstream Settings](#downstream-settings)
- [Prerequisites](#prerequisites)
  - [First-time Setup](#first-time-setup)
- [Output & Report](#output--report)
- [Advanced Usage](#advanced-usage)

---

## Features

- **12 annotation tools** covering plastid, mitochondrial, tRNA, and reference-based annotation
- **Automatic tool routing** — plastid samples get plastid tools, mito samples get mito tools
- **Flexible run modes** — run all, plastid-only, mito-only, or hand-pick tools
- **Mixed execution backends** — Conda, Docker, and Singularity/Apptainer supported
- **Unified QC layer** — BUSCO + cross-tool gene completeness comparison
- **Integrated downstream analysis** — RSCU, codon usage, Ka/Ks, phylogeny, composition, genome maps, synteny
- **Self-contained HTML reports** with per-tool sections, gene tables, BUSCO metrics, and downstream results
- **Cluster-ready** — SLURM, SGE, PBS via Snakemake's `--cluster` interface
- **Circular genome maps** via OGDraw (Docker) or pyGenomeViz (Python-native)

---

## Tool Overview

### Plastid / Chloroplast Annotation

| Tool | Type | Description | Reference |
|------|------|-------------|-----------|
| **Chloë** (Chloe.jl) | Conda (Julia) | XGBoost + suffix-array annotator optimised for angiosperm chloroplasts. Produces GFF, GenBank, EMBL, SFF. | [chloe.plastid.org](https://chloe.plastid.org) |
| **PGA** | Conda (Perl + BLAST) | Batch plastid genome annotator using TBLASTN against a curated reference GenBank collection. Detects IRs. | [PGA GitHub](https://github.com/quxiaojian/PGA) |
| **CPGAVAS2** | Docker | Comprehensive chloroplast annotator (BLAST + HMMER) with IR detection and circular map generation. | Shi et al. (2019) |

### Mitochondrial Annotation

| Tool | Type | Description | Reference |
|------|------|-------------|-----------|
| **MFannot** | Docker | Comprehensive mito/plastid annotator using BLAST, HMMER, Exonerate, Erpin. Excellent for intron-rich genomes. | [MFannot Docker](https://hub.docker.com/r/nbeck/mfannot) |
| **fpma** | Conda (Rust + HMMER) | Fast presence/absence scan of 43 core + 31 tRNA genes in angiosperm mitochondrial genomes via HMM profiles. | [fpma GitHub](https://github.com/liftoff/fpma) |
| **MITOS2** | Docker | Reference-based mitochondrial annotator for protein-coding genes, tRNAs, rRNAs. Supports metazoan and fungal refs. | Donath et al. (2019) |
| **MitoZ** | Docker | Animal mitochondrial genome annotator with automatic circular visualisation. Multiple clade-specific models. | [MitoZ GitHub](https://github.com/linzhi2013/MitoZ) |

### Both Organelles (tRNA / Reference-based)

| Tool | Type | Description | Reference |
|------|------|-------------|-----------|
| **tRNAscan-SE** | Conda | Gold-standard tRNA detection using covariance models. Organellar mode (`-O`) for mito/plastid tRNAs. | Lowe & Chan (2016) |
| **Aragorn** | Conda | Lightweight tRNA and tmRNA detection. Very fast; suitable as a second-opinion tRNA caller. | Laslett & Canback (2004) |
| **Liftoff** | Conda (minimap2) | Reference-based annotation lift-over. Maps features from a reference GFF+FASTA to a target genome. | Shumate & Salzberg (2021) |

### Visualisation

| Tool | Type | Description | Reference |
|------|------|-------------|-----------|
| **OGDraw** | Docker | Generates publication-quality circular and linear genome maps from GenBank annotation files. | Greiner et al. (2019) |

### Quality Control

| Tool | Type | Description |
|------|------|-------------|
| **BUSCO** | Conda | Genome completeness assessment against lineage-specific orthologue databases. |
| **Gene Completeness Summary** | Built-in (Python) | Cross-tool comparison of detected genes, tRNAs, and rRNAs, parsed from each annotator's output. |

---

## Downstream Analysis

Integrated post-annotation analyses that run automatically after annotation completes (enable via `downstream.enabled: true` in config). All downstream analyses produce per-sample results and a unified HTML report.

| Analysis | Tool / Method | Description |
|----------|:---:|-------------|
| **RSCU** | BioPython | Relative Synonymous Codon Usage heatmap and bar plots |
| **Codon Usage** | BioPython | Start/stop codon frequency analysis across all CDS |
| **Ka/Ks** | MAFFT + KaKs_Calculator2 | Pairwise synonymous/non-synonymous substitution rates vs. reference |
| **Phylogeny** | MAFFT + IQ-TREE | Supermatrix from shared genes → ML tree with ultrafast bootstrap |
| **GC/AA Composition** | BioPython | Per-gene GC content and aggregate amino acid frequency plots |
| **Genome Map** | pyGenomeViz | Circular genome visualisation (pure Python — no Perl/Circos needed) |
| **Synteny** | MUMmer4 / nucmer | Genome structure comparison with Bezier ribbon visualisation |
| **Reference Fetch** | NCBI Entrez | Automatic download of related reference genomes for comparative analysis |

**Key improvements over the original post-assembly pipeline:**

- **Circos replaced** with pyGenomeViz (Python-native, no 20+ Perl dependencies)
- **MUSCLE replaced** with MAFFT (faster, better for divergent sequences)
- **Duplicate code eliminated** — shared `gene_utils.py` for gene name normalisation
- **Duplicate Ka/Ks removed** — single implementation using proper KaKs_Calculator2
- **Hardcoded species removed** — fully configurable via `config.yaml`
- **Removed** `mitos_to_genbank.py` (redundant — main pipeline tools already produce GenBank)

---

## Quick Start

```bash
# 1. Edit the sample sheet
#    Columns: sample, fasta, organelle (plastid|mito), genetic_code
vim config/samples.tsv

# 2. Review / edit configuration
vim config/config.yaml

# 3. Dry-run to see what will be executed
snakemake -n --configfile config/config.yaml

# 4. Run the pipeline (all tools for each sample's organelle type)
snakemake --cores 8 --use-conda --configfile config/config.yaml

# 5. Run only plastid annotation tools
snakemake --cores 8 --use-conda --config mode=plastid

# 6. Run only mitochondrial annotation tools
snakemake --cores 8 --use-conda --config mode=mito

# 7. Run a hand-picked subset of tools
snakemake --cores 8 --use-conda --config mode=select

# Docker must be available for: MFannot, MITOS2, MitoZ, CPGAVAS2, OGDraw
```

---

## Directory Structure

```
Organelle_annotation/
├── Snakefile                          # Main workflow (tool routing + includes)
├── config/
│   ├── config.yaml                    # Pipeline configuration (all tools)
│   └── samples.tsv                    # Sample sheet (TSV)
├── rules/
│   ├── chloe.smk                      # Chloë (Julia)
│   ├── pga.smk                        # PGA (Perl + BLAST)
│   ├── cpgavas2.smk                   # CPGAVAS2 (Docker)
│   ├── mfannot.smk                    # MFannot (Docker)
│   ├── fpma.smk                       # fpma (Rust + HMMER)
│   ├── mitos.smk                      # MITOS2 (Docker)
│   ├── mitoz.smk                      # MitoZ (Docker)
│   ├── trnascan.smk                   # tRNAscan-SE
│   ├── aragorn.smk                    # Aragorn
│   ├── liftoff.smk                    # Liftoff (minimap2 lift-over)
│   ├── ogdraw.smk                     # OGDraw (Docker, visualisation)
│   ├── qc.smk                         # BUSCO + gene completeness
│   ├── report.smk                     # HTML report generation
│   └── downstream.smk                 # Post-annotation downstream analyses
├── scripts/
│   ├── gene_utils.py                  # Shared gene name maps & utilities
│   ├── generate_report.py             # Main report builder
│   ├── generate_downstream_report.py  # Downstream HTML report
│   ├── fetch_organelle_ref.py         # NCBI reference genome fetcher
│   ├── calculate_rscu.py              # RSCU analysis
│   ├── analyze_codons.py              # Start/stop codon analysis
│   ├── run_kaks_analysis.py           # Ka/Ks (MAFFT + KaKs_Calculator)
│   ├── prepare_phylo.py               # Supermatrix builder (MAFFT)
│   ├── plot_tree.py                   # Phylogenetic tree plotter
│   ├── analyze_composition.py         # GC/AA composition analysis
│   ├── create_genome_map.py           # Genome map (pyGenomeViz)
│   └── run_synteny_analysis.py        # Synteny analysis (MUMmer4)
├── envs/
│   ├── aragorn.yaml                   # Conda: aragorn
│   ├── busco.yaml                     # Conda: BUSCO
│   ├── chloe.yaml                     # Conda: Julia
│   ├── downstream.yaml                # Conda: downstream analysis tools
│   ├── fpma.yaml                      # Conda: HMMER + Rust
│   ├── liftoff.yaml                   # Conda: liftoff + minimap2
│   ├── pga.yaml                       # Conda: Perl + BLAST
│   ├── phylo.yaml                     # Conda: IQ-TREE + MAFFT
│   ├── synteny.yaml                   # Conda: MUMmer4
│   └── trnascan.yaml                  # Conda: tRNAscan-SE
├── repo/
│   ├── Chloe.jl/                      # Chloë source
│   ├── PGA/                           # PGA source
│   ├── Mfannot/                       # MFannot reference (Docker used)
│   └── fpma/                          # fpma source (Rust)
└── results/                           # Output (created by pipeline)
    ├── <tool>/<sample>/               # Per-tool annotation outputs
    ├── qc/
    │   ├── busco/<sample>/            # BUSCO results
    │   └── summary/<sample>.tsv       # Cross-tool gene completeness
    ├── downstream/<sample>/           # Downstream analysis results
    │   ├── rscu/                      # RSCU heatmaps and tables
    │   ├── codons/                    # Codon usage stats
    │   ├── kaks/                      # Ka/Ks summary tables
    │   ├── composition/               # GC content and AA plots
    │   ├── phylogeny/                 # Supermatrix, tree, partition files
    │   ├── genome_map/                # Circular genome map (PNG + SVG)
    │   ├── synteny/                   # Synteny plot and stats
    │   ├── references/                # Downloaded reference genomes
    │   └── downstream_report.html     # Per-sample downstream HTML report
    ├── logs/                          # Per-rule log files
    └── report/index.html              # Final HTML report
```

---

## Configuration

### Run Modes

Set `mode` in `config/config.yaml`:

| Mode | Plastid tools | Mito tools | Both/QC tools | Description |
|------|:---:|:---:|:---:|-------------|
| `all` | Chloë, PGA, CPGAVAS2 | MFannot, fpma, MITOS2, MitoZ | tRNAscan-SE, Aragorn, Liftoff | All compatible tools per sample |
| `plastid` | ✓ | — | ✓ | Plastid tools + both-organelle tools |
| `mito` | — | ✓ | ✓ | Mito tools + both-organelle tools |
| `select` | (user picks) | (user picks) | (user picks) | Only tools listed in `tools_select` |

### Sample Sheet

`config/samples.tsv` — tab-separated:

| Column | Description | Example |
|--------|-------------|---------|
| `sample` | Unique sample identifier | `arabidopsis_cp` |
| `fasta` | Path to input FASTA file | `/data/genomes/ath_cp.fasta` |
| `organelle` | `plastid` or `mito` | `plastid` |
| `genetic_code` | NCBI genetic code number | `11` (Plant Plastid) |

Common genetic codes for organelles:

| Code | Name | Use case |
|------|------|----------|
| 1 | Standard | Default |
| 2 | Vertebrate Mitochondrial | Vertebrate mt |
| 4 | Mold/Protozoan/Coelenterate Mito | Fungal/protist mt |
| 5 | Invertebrate Mitochondrial | Invertebrate mt |
| 11 | Bacterial/Plant Plastid | Plant cp and mt |

### Tool-specific Settings

All tool parameters are in `config/config.yaml` under their respective keys. Key settings requiring user input:

| Tool | Required Config | Notes |
|------|----------------|-------|
| **Liftoff** | `liftoff.reference_fasta`, `liftoff.reference_gff` | Reference FASTA + GFF3 from a related species |
| **CPGAVAS2** | — | Docker image auto-pulls |
| **MFannot** | — | Docker image: `nbeck/mfannot` |
| **MITOS2** | `mitos.ref_db` | Choose reference DB version (e.g. `refseq63m`) |
| **MitoZ** | `mitoz.clade` | Choose clade model (e.g. `Chordata`) |
| **fpma** | `fpma.hmms_subdir` | Choose HMM set (e.g. `angiosperm_hmms`, `fern_hmms`) |
| **OGDraw** | — | Docker image: `chlorobox/ogdraw:1.3.1` |

Docker images (`mfannot`, `mitos`, `mitoz`, `cpgavas2`, `ogdraw`) also support **Singularity/Apptainer** — set `use_singularity: true` in the respective config section.

### Downstream Settings

Enable integrated downstream analysis in `config/config.yaml`:

```yaml
downstream:
  enabled: true
  species_name: "Arabidopsis thaliana"   # For NCBI reference fetching
  email: "user@example.com"              # Required by NCBI Entrez
  max_ref_genomes: 10                    # Max references to download
  kaks_method: "NG"                      # NG, LWL, YN, MYN, GY
  phylo_min_genes: 4                     # Min shared genes for supermatrix
  phylo_model: "GTR+G"                   # IQ-TREE substitution model
  phylo_bootstrap: 1000                  # Ultrafast bootstrap replicates
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `enabled` | `true` | Enable/disable all downstream analyses |
| `species_name` | `""` | Species name for NCBI reference search |
| `email` | `""` | Email for NCBI Entrez API (required) |
| `max_ref_genomes` | `10` | Maximum reference genomes to download |
| `min_genome_length` | `10000` | Minimum genome length filter (bp) |
| `max_genome_length` | `300000` | Maximum genome length filter (bp) |
| `kaks_method` | `"NG"` | Ka/Ks calculation method |
| `phylo_min_genes` | `4` | Minimum genes for supermatrix construction |
| `phylo_model` | `"GTR+G"` | IQ-TREE substitution model |
| `phylo_bootstrap` | `1000` | Bootstrap replicates for ML tree |

---

## Prerequisites

### Software Requirements

| Tool | Requirement | Install |
|------|-------------|---------|
| **Snakemake** | ≥ 7.0 | `conda install -c bioconda snakemake` |
| **Conda** | Miniconda or Mamba | [docs.conda.io](https://docs.conda.io) |
| **Docker** | For MFannot, MITOS2, MitoZ, CPGAVAS2, OGDraw | [docker.com](https://docker.com) |
| **Chloë** | Julia ≥ 1.9 | Via conda or [`juliaup`](https://julialang.org/downloads/) |
| **PGA** | Perl ≥ 5.26, BLAST+ ≥ 2.8.1 | Via conda env |
| **fpma** | Rust toolchain, HMMER3 | `cargo build --release` |
| **Liftoff** | Python ≥ 3.8, minimap2 | Via conda env |
| **tRNAscan-SE** | ≥ 2.0.12 | Via conda env |
| **Aragorn** | ≥ 1.2.41 | Via conda env |
| **BUSCO** | ≥ 5.4 | Via conda env |

### First-time Setup

```bash
# Pull Docker images
docker pull nbeck/mfannot
docker pull quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0
docker pull guanliangmeng/mitoz:3.6
docker pull lipme/cpgavas2:latest
docker pull chlorobox/ogdraw:1.3.1

# Build fpma binary
cd repo/fpma && cargo build --release && cd ../..

# Install Chloë Julia dependencies
cd repo/Chloe.jl && julia --project=. -e 'using Pkg; Pkg.instantiate()' && cd ../..

# (Optional) Clone Chloë full reference set
# git clone https://github.com/ian-small/chloe_references
```

---

## Output & Report

### Per-tool Outputs

| Tool | Key output files |
|------|-----------------|
| **Chloë** | `{sample}.gff`, `{sample}.gbk`, `{sample}.sff`, `{sample}.embl` |
| **PGA** | `{sample}.gb` (GenBank) |
| **CPGAVAS2** | `{sample}.gb`, `{sample}.gff`, circular map images |
| **MFannot** | `{sample}.new` (masterfile), `{sample}.gff` |
| **fpma** | `{sample}.gff`, `{sample}.presence.tsv`, `{sample}.html` (SVG plot) |
| **MITOS2** | `result.gff`, `result.bed`, protein FASTA |
| **MitoZ** | `{sample}.gff`, `{sample}.gbk`, circular PNG |
| **tRNAscan-SE** | `{sample}.trnascan.tsv`, `{sample}.gff`, `{sample}.ss` (secondary structure) |
| **Aragorn** | `{sample}.aragorn.txt`, `{sample}.gff` |
| **Liftoff** | `{sample}.gff`, `{sample}.unmapped.txt`, `{sample}.gb` |
| **OGDraw** | `{sample}_map.svg` (circular genome map) |

### QC Outputs

| File | Content |
|------|---------|
| `qc/busco/{sample}/short_summary.txt` | BUSCO completeness metrics (Complete, Fragmented, Missing) |
| `qc/summary/{sample}.qc_summary.tsv` | Cross-tool comparison: gene count, tRNA count, rRNA count, gene names per tool |

### HTML Report

`results/report/index.html` — self-contained, navigable report containing:

1. **Pipeline Overview** — samples processed and which tools were run
2. **Per-Tool Sections** — output file listings with sizes and download links
3. **Gene Completeness Summary** — side-by-side table of genes/tRNAs/rRNAs detected by each tool per sample
4. **BUSCO Assessment** — completeness metrics (Complete, Single-copy, Duplicated, Fragmented, Missing)
5. **Downstream Analysis** — links to per-sample downstream reports with status indicators

### Downstream Outputs

Each sample gets a dedicated `downstream/<sample>/` directory containing:

| Directory | Key Files | Description |
|-----------|-----------|-------------|
| `rscu/` | `rscu.tsv`, `rscu_barplot.png`, `rscu_heatmap.png` | RSCU values and visualisations |
| `codons/` | `codon_stats.txt` | Start/stop codon frequencies |
| `kaks/` | `kaks_summary.tsv` | Per-gene Ka, Ks, Ka/Ks values |
| `composition/` | `gc_content_plot.png`, `aa_composition_plot.png` | GC & amino acid plots |
| `phylogeny/` | `supermatrix.fasta`, `phylogeny.treefile`, `tree_plot.png` | Alignment, ML tree, tree plot |
| `genome_map/` | `genome_map.png`, `genome_map.svg` | Circular genome map (pyGenomeViz) |
| `synteny/` | `synteny_plot.png`, `synteny_stats.tsv` | Synteny ribbon plot & statistics |
| `references/` | `*.fasta`, `*.gbk` | Downloaded NCBI reference genomes |
| — | `downstream_report.html` | Self-contained per-sample downstream report |

---

## Advanced Usage

```bash
# Cluster execution (SLURM)
snakemake --cores 100 --use-conda \
  --cluster "sbatch -p normal -c {threads} --mem={resources.mem_mb}M -t {resources.runtime}" \
  --configfile config/config.yaml

# Singularity instead of Docker
snakemake --cores 8 --use-conda --use-singularity \
  --configfile config/config.yaml

# Only regenerate the report from existing outputs
snakemake --cores 1 --use-conda results/report/index.html

# Force re-run a specific tool for a sample
snakemake --cores 4 --use-conda -f results/chloe/sample1/sample1.done

# Run selected tools only
snakemake --cores 8 --use-conda \
  --config mode=select tools_select="[chloe,pga,trnascan,aragorn]"

# Override output directory
snakemake --cores 8 --use-conda --config outdir=my_results

# Run annotation only (disable downstream analysis)
snakemake --cores 8 --use-conda --config downstream="{enabled: false}"

# Run downstream analysis for a specific sample
snakemake --cores 4 --use-conda results/downstream/sample1/downstream_report.html

# Run phylogeny only for a sample
snakemake --cores 4 --use-conda results/downstream/sample1/phylogeny/phylogeny.treefile
```

---

## License

See individual tool repositories for their respective licenses.
