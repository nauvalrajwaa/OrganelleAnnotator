#!/usr/bin/env python3
"""
fetch_organelle_ref.py – Fetch reference organelle genomes from NCBI.

Downloads one or more reference genomes (FASTA + GFF) for a given species
or taxonomic group from NCBI Entrez.  Used to provide reference sequences
for Ka/Ks analysis, phylogeny construction, and synteny comparison.

Usage (CLI):
    python fetch_organelle_ref.py <species_name> <email> <output_dir> \\
        [--organelle mito|plastid] [--genetic_code 2] \\
        [--max_genomes 10] [--min_len 10000] [--max_len 200000]

Usage (Snakemake):
    Called via `script:` directive.
"""

import sys
import os
import time
import argparse
import logging

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

try:
    from Bio import Entrez, SeqIO
except ImportError:
    log.error("BioPython is required: pip install biopython")
    sys.exit(1)


# ── NCBI Search ──────────────────────────────────────────────────────────────

ORGANELLE_FILTERS = {
    "mito": "mitochondrion[filter]",
    "plastid": "chloroplast[filter] OR plastid[filter]",
}


def build_query(species: str, organelle: str, min_len: int, max_len: int) -> str:
    """Build an NCBI Entrez search query."""
    organelle_filter = ORGANELLE_FILTERS.get(organelle, "mitochondrion[filter]")
    return (
        f'("{species}"[Organism]) AND ({organelle_filter}) '
        f'AND ("{min_len}"[SLEN]:"{max_len}"[SLEN]) '
        f'AND (refseq[filter] OR "complete genome")'
    )


def search_ncbi(query: str, email: str, max_results: int = 10) -> list[str]:
    """Search NCBI Nucleotide database, return list of accession IDs."""
    Entrez.email = email
    log.info(f"NCBI query: {query}")

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results, usehistory="y")
    results = Entrez.read(handle)
    handle.close()

    ids = results.get("IdList", [])
    total = results.get("Count", "0")
    log.info(f"Found {total} records, retrieving up to {len(ids)} IDs")
    return ids


def fetch_sequences(ids: list[str], email: str, output_dir: str) -> list[str]:
    """Fetch GenBank records and write FASTA + GFF files.

    Returns list of accession strings that were successfully saved.
    """
    Entrez.email = email
    os.makedirs(output_dir, exist_ok=True)
    saved = []

    for gid in ids:
        try:
            # Fetch GenBank format
            handle = Entrez.efetch(db="nucleotide", id=gid, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            accession = record.id.replace(".", "_")
            fasta_path = os.path.join(output_dir, f"{accession}.fasta")
            gbk_path = os.path.join(output_dir, f"{accession}.gbk")

            # Write FASTA
            SeqIO.write(record, fasta_path, "fasta")
            # Write GenBank (useful for downstream)
            SeqIO.write(record, gbk_path, "genbank")

            log.info(f"Saved {accession}: {record.description[:80]}... ({len(record.seq):,} bp)")
            saved.append(accession)

            time.sleep(0.4)  # Respect NCBI rate limit

        except Exception as exc:
            log.warning(f"Failed to fetch ID {gid}: {exc}")
            continue

    return saved


def expand_search(species: str, organelle: str, email: str,
                  min_len: int, max_len: int, max_results: int,
                  output_dir: str) -> list[str]:
    """Broaden search to genus level if species-level returns too few hits."""
    query = build_query(species, organelle, min_len, max_len)
    ids = search_ncbi(query, email, max_results)

    if len(ids) < 3:
        genus = species.split()[0] if " " in species else species
        log.info(f"Few results for '{species}'; expanding to genus '{genus}'")
        query = build_query(genus, organelle, min_len, max_len)
        ids = search_ncbi(query, email, max_results)

    if not ids:
        log.warning("No reference genomes found on NCBI for given criteria.")
        return []

    return fetch_sequences(ids, email, output_dir)


# ── CLI / Snakemake entry points ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Fetch reference organelle genomes from NCBI Entrez."
    )
    parser.add_argument("species", help="Species or genus name (e.g. 'Arabidopsis thaliana')")
    parser.add_argument("email", help="Email for NCBI Entrez (required by NCBI)")
    parser.add_argument("output_dir", help="Output directory for downloaded files")
    parser.add_argument("--organelle", default="mito", choices=["mito", "plastid"])
    parser.add_argument("--max_genomes", type=int, default=10)
    parser.add_argument("--min_len", type=int, default=10000)
    parser.add_argument("--max_len", type=int, default=300000)

    args = parser.parse_args()

    saved = expand_search(
        species=args.species,
        organelle=args.organelle,
        email=args.email,
        min_len=args.min_len,
        max_len=args.max_len,
        max_results=args.max_genomes,
        output_dir=args.output_dir,
    )
    log.info(f"Downloaded {len(saved)} reference genome(s) to {args.output_dir}")


if __name__ == "__main__":
    main()
