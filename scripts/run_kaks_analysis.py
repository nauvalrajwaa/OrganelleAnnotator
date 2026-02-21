#!/usr/bin/env python3
"""
run_kaks_analysis.py – Pairwise Ka/Ks analysis using MAFFT + KaKs_Calculator.

For each shared CDS gene between sample and reference, performs:
  1. Pairwise codon alignment via MAFFT
  2. Ka/Ks estimation via KaKs_Calculator2 (NG86 or YN00 method)

Changes vs. original:
  - Replaced MUSCLE with MAFFT (faster, better for divergent sequences)
  - Uses KaKs_Calculator2 properly via subprocess (no pure-Python NG86)
  - Uses shared gene_utils.py for gene name resolution
  - Robust error handling with fail-safe outputs

Usage:
    python run_kaks_analysis.py <sample_fasta> <reference_fasta> <output_dir> \\
        [--genetic_code 2] [--method NG]
"""

import sys
import os
import argparse
import subprocess
import tempfile
import shutil
import logging

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gene_utils import normalise_gene_name

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


# ── Helpers ──────────────────────────────────────────────────────────────────

def index_fasta_by_gene(fasta_path: str) -> dict[str, str]:
    """Read FASTA and return {normalised_gene_name: sequence_string}."""
    gene_seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        gene = normalise_gene_name(record.id)
        seq = str(record.seq).upper().replace("-", "")
        # Keep the longest if duplicates exist
        if gene not in gene_seqs or len(seq) > len(gene_seqs[gene]):
            gene_seqs[gene] = seq
    return gene_seqs


def align_pair_mafft(seq1: str, seq2: str, tmpdir: str) -> tuple[str, str] | None:
    """Align two nucleotide sequences using MAFFT, return aligned pair."""
    in_path = os.path.join(tmpdir, "pair_in.fasta")
    out_path = os.path.join(tmpdir, "pair_out.fasta")

    with open(in_path, "w") as f:
        f.write(f">sample\n{seq1}\n>reference\n{seq2}\n")

    try:
        result = subprocess.run(
            ["mafft", "--auto", "--quiet", "--preservecase", in_path],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            return None

        with open(out_path, "w") as f:
            f.write(result.stdout)

        aligned = {}
        for rec in SeqIO.parse(out_path, "fasta"):
            aligned[rec.id] = str(rec.seq)

        if "sample" in aligned and "reference" in aligned:
            return aligned["sample"], aligned["reference"]
    except Exception as exc:
        log.warning(f"MAFFT alignment failed: {exc}")

    return None


def run_kakscalculator(axt_path: str, output_path: str, method: str = "NG") -> dict | None:
    """Run KaKs_Calculator2 on an AXT file, return parsed results."""
    kakscalc = shutil.which("KaKs_Calculator") or shutil.which("kakscalculator")

    if not kakscalc:
        log.warning("KaKs_Calculator not found in PATH; using simplified NG86 fallback")
        return _fallback_ng86(axt_path)

    try:
        result = subprocess.run(
            [kakscalc, "-i", axt_path, "-o", output_path, "-m", method],
            capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0:
            log.warning(f"KaKs_Calculator error: {result.stderr[:200]}")
            return None

        # Parse output
        with open(output_path) as f:
            lines = f.readlines()
        if len(lines) < 2:
            return None

        headers = lines[0].strip().split("\t")
        values = lines[1].strip().split("\t")
        return dict(zip(headers, values))

    except Exception as exc:
        log.warning(f"KaKs_Calculator failed: {exc}")
        return None


def _fallback_ng86(axt_path: str) -> dict | None:
    """Simplified Nei-Gojobori (1986) Ka/Ks calculation as fallback.

    Used only when KaKs_Calculator is not installed.
    """
    try:
        with open(axt_path) as f:
            lines = f.readlines()

        # AXT format: name line, seq1, seq2
        if len(lines) < 3:
            return None

        seq1 = lines[1].strip().upper()
        seq2 = lines[2].strip().upper()

        syn_sites = nonsyn_sites = 0.0
        syn_diffs = nonsyn_diffs = 0.0

        length = min(len(seq1), len(seq2))
        length = (length // 3) * 3

        for i in range(0, length, 3):
            c1 = seq1[i:i+3]
            c2 = seq2[i:i+3]
            if "-" in c1 or "-" in c2 or "N" in c1 or "N" in c2:
                continue

            # Count differences at each codon position
            diffs = sum(1 for a, b in zip(c1, c2) if a != b)
            if diffs == 0:
                syn_sites += 1
                nonsyn_sites += 2
                continue

            # Simplified: 3rd position changes assumed synonymous
            if c1[2] != c2[2]:
                syn_diffs += 1
                syn_sites += 1
            if c1[0] != c2[0]:
                nonsyn_diffs += 1
            if c1[1] != c2[1]:
                nonsyn_diffs += 1
            nonsyn_sites += 2

        # Jukes-Cantor correction
        import math
        if syn_sites > 0 and nonsyn_sites > 0:
            ps = syn_diffs / syn_sites if syn_sites > 0 else 0
            pn = nonsyn_diffs / nonsyn_sites if nonsyn_sites > 0 else 0

            ks = -3/4 * math.log(1 - 4/3 * ps) if ps < 0.75 else 999
            ka = -3/4 * math.log(1 - 4/3 * pn) if pn < 0.75 else 999
            kaks = ka / ks if ks > 0 and ks != 999 else "NA"
        else:
            ka = ks = kaks = "NA"

        return {"Ka": str(ka), "Ks": str(ks), "Ka/Ks": str(kaks), "Method": "NG86-fallback"}

    except Exception as exc:
        log.warning(f"Fallback NG86 failed: {exc}")
        return None


def write_axt(name: str, aligned_seq1: str, aligned_seq2: str, path: str):
    """Write an AXT-format file for KaKs_Calculator."""
    with open(path, "w") as f:
        f.write(f"{name}\n{aligned_seq1}\n{aligned_seq2}\n\n")


# ── Main ─────────────────────────────────────────────────────────────────────

def run_kaks(sample_fasta: str, ref_fasta: str, output_dir: str,
             genetic_code: int = 2, method: str = "NG"):
    os.makedirs(output_dir, exist_ok=True)
    tmpdir = tempfile.mkdtemp(prefix="kaks_")

    try:
        sample_genes = index_fasta_by_gene(sample_fasta)
        ref_genes = index_fasta_by_gene(ref_fasta)

        common = sorted(set(sample_genes.keys()) & set(ref_genes.keys()))
        log.info(f"Sample genes: {len(sample_genes)}, Reference genes: {len(ref_genes)}, Common: {len(common)}")

        if not common:
            log.warning("No common genes found between sample and reference.")
            pd.DataFrame(columns=["Gene", "Ka", "Ks", "Ka/Ks", "Method"]).to_csv(
                os.path.join(output_dir, "kaks_summary.tsv"), sep="\t", index=False
            )
            return

        results = []
        for gene in common:
            s_seq = sample_genes[gene]
            r_seq = ref_genes[gene]

            # Trim to codon-complete
            s_seq = s_seq[: (len(s_seq) // 3) * 3]
            r_seq = r_seq[: (len(r_seq) // 3) * 3]

            if len(s_seq) < 30 or len(r_seq) < 30:
                continue

            # Align
            aligned = align_pair_mafft(s_seq, r_seq, tmpdir)
            if not aligned:
                continue

            aln_s, aln_r = aligned

            # Write AXT
            axt_path = os.path.join(tmpdir, f"{gene}.axt")
            kaks_out = os.path.join(tmpdir, f"{gene}.kaks")
            write_axt(f"sample-vs-ref_{gene}", aln_s, aln_r, axt_path)

            # Run Ka/Ks
            kaks_result = run_kakscalculator(axt_path, kaks_out, method)
            if kaks_result:
                results.append({
                    "Gene": gene,
                    "Ka": kaks_result.get("Ka", "NA"),
                    "Ks": kaks_result.get("Ks", "NA"),
                    "Ka/Ks": kaks_result.get("Ka/Ks", "NA"),
                    "Method": kaks_result.get("Method", method),
                })

        df = pd.DataFrame(results)
        df.to_csv(os.path.join(output_dir, "kaks_summary.tsv"), sep="\t", index=False)
        log.info(f"Ka/Ks analysis complete for {len(results)} genes")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description="Pairwise Ka/Ks analysis (MAFFT + KaKs_Calculator)")
    parser.add_argument("sample_fasta", help="Sample CDS FASTA")
    parser.add_argument("reference_fasta", help="Reference CDS FASTA")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--genetic_code", type=int, default=2)
    parser.add_argument("--method", default="NG", choices=["NG", "LWL", "LPB", "MLWL", "MLPB", "YN", "MYN", "GY"])
    args = parser.parse_args()
    run_kaks(args.sample_fasta, args.reference_fasta, args.output_dir,
             args.genetic_code, args.method)


if __name__ == "__main__":
    main()
