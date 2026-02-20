# ========================================================
# Organelle Annotation Pipeline
# ========================================================
#
# A Snakemake workflow for organelle genome annotation using multiple tools:
#
#   Plastid/Chloroplast:
#     - Chloë (Chloe.jl) – Julia-based chloroplast annotator
#     - PGA              – Perl/BLAST plastid genome annotator
#
#   Mitochondrial:
#     - MFannot          – Docker-based comprehensive annotator
#     - fpma             – Fast HMM-based gene presence/absence scanner
#
#   QC:
#     - BUSCO            – Genome completeness assessment
#     - Gene completeness summary (custom)
#
#   Report:
#     - Indexed HTML report with per-tool sections
#
# =============================================================================

## Quick Start

```bash
# 1. Edit the sample sheet
#    Columns: sample, fasta, organelle (plastid|mito), genetic_code
vim config/samples.tsv

# 2. Review / edit configuration
vim config/config.yaml

# 3. Dry-run to see what will be executed
snakemake -n --configfile config/config.yaml

# 4. Run the pipeline
snakemake --cores 8 --use-conda --configfile config/config.yaml

# 5. Run only specific tools via mode override
snakemake --cores 8 --use-conda --config mode=plastid

# 6. Run with Docker support (needed for MFannot)
snakemake --cores 8 --use-conda --configfile config/config.yaml
# (Docker must be available; MFannot rule calls `docker run` directly)
```

## Directory Structure

```
Organelle_annotation/
├── Snakefile                     # Main workflow
├── config/
│   ├── config.yaml               # Pipeline configuration
│   └── samples.tsv               # Sample sheet (TSV)
├── rules/
│   ├── chloe.smk                 # Chloë rules
│   ├── pga.smk                   # PGA rules
│   ├── mfannot.smk               # MFannot rules (Docker)
│   ├── fpma.smk                  # fpma rules
│   ├── qc.smk                    # BUSCO + gene completeness
│   └── report.smk                # HTML report generation
├── scripts/
│   └── generate_report.py        # Report builder
├── envs/
│   ├── busco.yaml                # Conda env: BUSCO
│   ├── chloe.yaml                # Conda env: Julia
│   ├── fpma.yaml                 # Conda env: HMMER + Rust
│   └── pga.yaml                  # Conda env: Perl + BLAST
├── Chloe.jl/                     # Chloë source
├── PGA/                          # PGA source
├── Mfannot/                      # MFannot reference (Docker used)
├── fpma/                         # fpma source
└── results/                      # Output (created by pipeline)
    ├── chloe/<sample>/           # Chloë outputs (GFF, GBK, etc.)
    ├── pga/<sample>/             # PGA outputs (GenBank)
    ├── mfannot/<sample>/         # MFannot outputs (masterfile)
    ├── fpma/<sample>/            # fpma outputs (GFF, TSV, HTML)
    ├── qc/
    │   ├── busco/<sample>/       # BUSCO results
    │   └── summary/<sample>.tsv  # Gene completeness summaries
    ├── logs/                     # Per-rule log files
    └── report/index.html         # Final HTML report
```

## Configuration

### Run Modes (`config.yaml → mode`)

| Mode      | Description                                           |
|-----------|-------------------------------------------------------|
| `all`     | Run all tools compatible with each sample's organelle |
| `plastid` | Run only plastid tools (Chloë, PGA)                  |
| `mito`    | Run only mitochondrial tools (MFannot, fpma)          |
| `select`  | Run only tools listed in `tools_select`               |

### Sample Sheet (`config/samples.tsv`)

Tab-separated with columns:

| Column        | Description                                       |
|---------------|---------------------------------------------------|
| `sample`      | Unique sample identifier                          |
| `fasta`       | Absolute or relative path to input FASTA file     |
| `organelle`   | `plastid` or `mito`                               |
| `genetic_code`| NCBI genetic code (e.g. 11=Plant Plastid, 4=Mold) |

## Prerequisites

| Tool    | Requirement                                    |
|---------|------------------------------------------------|
| Chloë   | Julia ≥ 1.9 (installed via conda or `juliaup`) |
| PGA     | Perl ≥ 5.26, BLAST+ ≥ 2.8.1                   |
| MFannot | Docker (image: `nbeck/mfannot`)                |
| fpma    | Rust toolchain, HMMER3 (`nhmmer`)              |
| BUSCO   | BUSCO ≥ 5.4 (via conda)                        |

### First-time setup

```bash
# Pull MFannot Docker image
docker pull nbeck/mfannot

# Build fpma binary
cd fpma && cargo build --release && cd ..

# Install Chloë Julia dependencies
cd Chloe.jl && julia --project=. -e 'using Pkg; Pkg.instantiate()' && cd ..

# (Optional) Clone Chloë references for full reference set
# git clone https://github.com/ian-small/chloe_references
```

## Output Report

The pipeline generates `results/report/index.html` — a self-contained HTML
report with:

1. **Pipeline Overview** — sample list and which tools were run
2. **Per-Tool Sections** — output files with download links for each tool
3. **Gene Completeness** — cross-tool comparison of detected genes, tRNAs, rRNAs
4. **BUSCO Assessment** — genome completeness metrics

## Advanced Usage

```bash
# Cluster execution (SLURM example)
snakemake --cores 100 --use-conda \
  --cluster "sbatch -p normal -c {threads} --mem={resources.mem_mb}M -t {resources.runtime}" \
  --configfile config/config.yaml

# Only run QC for already-annotated samples
snakemake --cores 4 --use-conda results/report/index.html

# Force re-run a specific tool for a sample
snakemake --cores 4 --use-conda -f results/chloe/sample1/sample1.done
```
# OrganelleAnnotator
