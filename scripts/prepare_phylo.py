#!/usr/bin/env python3
"""
prepare_phylo.py – Build a supermatrix from shared genes and run IQ-TREE.

For each gene in common across sample + references, aligns with MAFFT,
concatenates into a supermatrix, writes a NEXUS partition file, and
(optionally) invokes IQ-TREE for maximum-likelihood phylogeny.

Changes vs. original:
  - Replaced MUSCLE with MAFFT (faster, handles divergent sequences better)
  - Uses shared gene_utils.py for gene name resolution
  - Cleaner partition file generation

Usage:
    python prepare_phylo.py <sample_fasta> <ref_dir> <output_dir> \\
        [--min_genes 4] [--run_iqtree]
"""

import sys
import os
import argparse
import subprocess
import logging
from collections import OrderedDict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gene_utils import normalise_gene_name

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_genes_from_fasta(fasta_path: str) -> dict[str, str]:
    """Load gene sequences from a CDS FASTA. Returns {normalised_name: seq}."""
    genes = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        name = normalise_gene_name(rec.id)
        seq = str(rec.seq).upper().replace("-", "")
        if name not in genes or len(seq) > len(genes[name]):
            genes[name] = seq
    return genes


def load_all_taxa(sample_fasta: str, ref_dir: str) -> dict[str, dict[str, str]]:
    """Load gene sequences for all taxa (sample + references).

    Returns {taxon_name: {gene_name: sequence}}.
    """
    taxa = {}

    # Sample
    sample_name = os.path.splitext(os.path.basename(sample_fasta))[0]
    taxa[sample_name] = load_genes_from_fasta(sample_fasta)
    log.info(f"Sample '{sample_name}': {len(taxa[sample_name])} genes")

    # References from ref_dir (any .fasta/.fa/.fna files)
    if ref_dir and os.path.isdir(ref_dir):
        for fn in sorted(os.listdir(ref_dir)):
            if fn.endswith((".fasta", ".fa", ".fna")):
                path = os.path.join(ref_dir, fn)
                taxon = os.path.splitext(fn)[0]
                taxa[taxon] = load_genes_from_fasta(path)
                log.info(f"Reference '{taxon}': {len(taxa[taxon])} genes")

    return taxa


def find_common_genes(taxa: dict[str, dict[str, str]], min_taxa: int = 2) -> list[str]:
    """Find genes present in at least min_taxa taxa."""
    from collections import Counter
    gene_count = Counter()
    for taxon_genes in taxa.values():
        for gene in taxon_genes:
            gene_count[gene] += 1

    common = [g for g, c in gene_count.items() if c >= min_taxa]
    return sorted(common)


def align_gene_mafft(sequences: dict[str, str]) -> dict[str, str] | None:
    """Align sequences for a single gene using MAFFT. Returns {taxon: aligned_seq}."""
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        for taxon, seq in sequences.items():
            tmp.write(f">{taxon}\n{seq}\n")
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            ["mafft", "--auto", "--quiet", "--preservecase", tmp_path],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            return None

        aligned = {}
        from io import StringIO
        for rec in SeqIO.parse(StringIO(result.stdout), "fasta"):
            aligned[rec.id] = str(rec.seq)
        return aligned

    except Exception as exc:
        log.warning(f"MAFFT failed: {exc}")
        return None
    finally:
        os.unlink(tmp_path)


# ── Supermatrix construction ─────────────────────────────────────────────────

def build_supermatrix(
    taxa: dict[str, dict[str, str]],
    common_genes: list[str],
) -> tuple[dict[str, str], list[tuple[str, int, int]]]:
    """Build concatenated supermatrix from aligned genes.

    Returns:
        supermatrix: {taxon: concatenated_aligned_sequence}
        partitions: [(gene_name, start_pos, end_pos), ...]
    """
    taxon_names = sorted(taxa.keys())
    supermatrix = {t: "" for t in taxon_names}
    partitions = []
    current_pos = 1

    for gene in common_genes:
        # Collect sequences for this gene from all taxa that have it
        gene_seqs = {}
        for taxon in taxon_names:
            if gene in taxa[taxon]:
                gene_seqs[taxon] = taxa[taxon][gene]

        if len(gene_seqs) < 2:
            continue

        # Align
        aligned = align_gene_mafft(gene_seqs)
        if not aligned:
            continue

        # Determine alignment length
        aln_len = max(len(s) for s in aligned.values())

        # Add to supermatrix (pad missing taxa with gaps)
        for taxon in taxon_names:
            if taxon in aligned:
                seq = aligned[taxon].ljust(aln_len, "-")
            else:
                seq = "-" * aln_len
            supermatrix[taxon] += seq

        partitions.append((gene, current_pos, current_pos + aln_len - 1))
        current_pos += aln_len

    return supermatrix, partitions


def write_supermatrix(supermatrix: dict[str, str], output_path: str):
    """Write supermatrix as aligned FASTA."""
    records = []
    for taxon, seq in sorted(supermatrix.items()):
        records.append(SeqRecord(Seq(seq), id=taxon, description=""))
    SeqIO.write(records, output_path, "fasta")


def write_partition_file(partitions: list[tuple[str, int, int]], output_path: str):
    """Write RAxML/IQ-TREE partition file."""
    with open(output_path, "w") as f:
        f.write("#nexus\nbegin sets;\n")
        for gene, start, end in partitions:
            f.write(f"  charset {gene} = {start}-{end};\n")
        f.write("end;\n")


# ── Main ─────────────────────────────────────────────────────────────────────

def prepare_phylo(sample_fasta: str, ref_dir: str, output_dir: str,
                  min_genes: int = 4, run_iqtree: bool = False):
    os.makedirs(output_dir, exist_ok=True)

    # Load all taxa
    taxa = load_all_taxa(sample_fasta, ref_dir)
    if len(taxa) < 2:
        log.error("Need at least 2 taxa (sample + 1 reference) for phylogeny.")
        return

    # Find common genes
    min_taxa = max(2, len(taxa) // 2)  # Present in at least half of taxa
    common = find_common_genes(taxa, min_taxa=min_taxa)
    log.info(f"Found {len(common)} genes present in >= {min_taxa} taxa")

    if len(common) < min_genes:
        log.warning(f"Too few common genes ({len(common)} < {min_genes}); skipping phylogeny.")
        # Write empty outputs so Snakemake doesn't fail on missing files
        open(os.path.join(output_dir, "supermatrix.fasta"), "w").close()
        open(os.path.join(output_dir, "partitions.nex"), "w").close()
        return

    # Build supermatrix
    supermatrix, partitions = build_supermatrix(taxa, common)

    # Check matrix is not empty
    total_len = max(len(s) for s in supermatrix.values()) if supermatrix else 0
    if total_len == 0:
        log.warning("Supermatrix is empty after alignment; skipping.")
        return

    log.info(f"Supermatrix: {len(supermatrix)} taxa x {total_len} bp ({len(partitions)} gene partitions)")

    # Write outputs
    matrix_path = os.path.join(output_dir, "supermatrix.fasta")
    partition_path = os.path.join(output_dir, "partitions.nex")
    write_supermatrix(supermatrix, matrix_path)
    write_partition_file(partitions, partition_path)
    log.info(f"Written: {matrix_path}, {partition_path}")

    # Optionally run IQ-TREE
    if run_iqtree:
        run_iqtree_ml(matrix_path, partition_path, output_dir)


def run_iqtree_ml(matrix_path: str, partition_path: str, output_dir: str):
    """Run IQ-TREE maximum-likelihood phylogeny."""
    prefix = os.path.join(output_dir, "phylogeny")
    cmd = [
        "iqtree2", "-s", matrix_path,
        "-p", partition_path,
        "-m", "GTR+G",
        "-B", "1000",          # Ultrafast bootstrap
        "--prefix", prefix,
        "-T", "AUTO",
        "--redo",
    ]
    log.info(f"Running IQ-TREE: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            log.error(f"IQ-TREE failed:\n{result.stderr[:500]}")
        else:
            log.info(f"IQ-TREE completed. Tree: {prefix}.treefile")
    except FileNotFoundError:
        # Try iqtree (v1 name) as fallback
        cmd[0] = "iqtree"
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if result.returncode != 0:
                log.error(f"IQ-TREE failed:\n{result.stderr[:500]}")
        except Exception as exc:
            log.error(f"IQ-TREE not found: {exc}")
    except Exception as exc:
        log.error(f"IQ-TREE error: {exc}")


def main():
    parser = argparse.ArgumentParser(description="Build supermatrix & run phylogeny")
    parser.add_argument("sample_fasta", help="Sample CDS FASTA")
    parser.add_argument("ref_dir", help="Directory with reference CDS FASTA files")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--min_genes", type=int, default=4)
    parser.add_argument("--run_iqtree", action="store_true", help="Run IQ-TREE after building supermatrix")
    args = parser.parse_args()
    prepare_phylo(args.sample_fasta, args.ref_dir, args.output_dir,
                  args.min_genes, args.run_iqtree)


if __name__ == "__main__":
    main()
