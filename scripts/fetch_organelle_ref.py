#!/usr/bin/env python3
"""
fetch_organelle_ref.py – Fetch reference organelle genomes from NCBI.

Downloads one or more reference genomes (FASTA + GFF) for a given species
or taxonomic group from NCBI Entrez. Used to provide reference sequences
for Ka/Ks analysis, phylogeny construction, and synteny comparison.
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
    from Bio.SeqRecord import SeqRecord
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
    organelle_filter = ORGANELLE_FILTERS.get(organelle, "chloroplast[filter] OR plastid[filter]")
    return (
        f'("{species}"[Organism]) AND ({organelle_filter}) '
        f'AND ("{min_len}"[SLEN]:"{max_len}"[SLEN]) '
        f'AND (refseq[filter] OR "complete genome")'
    )


def search_ncbi(query: str, email: str, max_results: int = 10) -> list[str]:
    """Search NCBI Nucleotide database, return list of accession IDs."""
    Entrez.email = email
    log.info(f"NCBI query: {query}")

    try:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results, usehistory="y")
        results = Entrez.read(handle)
        handle.close()
    except Exception as e:
        log.error(f"Gagal terhubung ke NCBI: {e}")
        return []

    ids = results.get("IdList", [])
    total = results.get("Count", "0")
    log.info(f"Found {total} records, retrieving up to {len(ids)} IDs")
    return ids


def fetch_sequences(ids: list[str], email: str, output_dir: str) -> list[str]:
    """Fetch GenBank records, write whole genome FASTA, GBK, and extracted CDS FASTA.

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

            # Bersihkan titik pada ID (misal: NC_012345.1 -> NC_012345_1)
            accession = record.id.replace(".", "_")
            fasta_path = os.path.join(output_dir, f"{accession}.fasta")
            cds_path   = os.path.join(output_dir, f"{accession}_cds.fasta")
            gbk_path   = os.path.join(output_dir, f"{accession}.gbk")

            # 1. Write whole genome FASTA
            SeqIO.write(record, fasta_path, "fasta")
            
            # 2. Write GenBank
            SeqIO.write(record, gbk_path, "genbank")

            # 3. Extract and write CDS
            cds_records = []
            for feature in record.features:
                if feature.type == "CDS":
                    # Ambil nama gen atau locus tag sebagai ID
                    gene_names = feature.qualifiers.get("gene", feature.qualifiers.get("locus_tag", ["unknown"]))
                    gene_name = gene_names[0]
                    protein_id = feature.qualifiers.get("protein_id", ["unknown_id"])[0]
                    
                    # Ekstrak sekuens nukleotida untuk CDS ini
                    try:
                        seq = feature.extract(record.seq)
                        cds_rec = SeqRecord(
                            seq,
                            # Format ini sangat penting agar terbaca oleh prepare_phylo.py
                            id=f"{accession}_{gene_name}",
                            description=f"protein_id={protein_id}"
                        )
                        cds_records.append(cds_rec)
                    except Exception as e:
                        log.debug(f"Could not extract CDS {gene_name} from {accession}: {e}")

            if cds_records:
                SeqIO.write(cds_records, cds_path, "fasta")
                log.info(f"Saved {accession}: Genome ({len(record.seq):,} bp) and {len(cds_records)} CDS sequences.")
            else:
                log.warning(f"Saved {accession}, but found 0 CDS features. Skipping CDS fasta.")

            saved.append(accession)
            time.sleep(0.5)  # Respect NCBI rate limit (max 3 requests per second)

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
    
    # Argumen bisa ditangkap dari eksekusi Snakemake
    parser.add_argument("species", nargs='?', help="Species or genus name (e.g. 'Arabidopsis thaliana')")
    parser.add_argument("email", nargs='?', help="Email for NCBI Entrez (required by NCBI)")
    parser.add_argument("output_dir", nargs='?', help="Output directory for downloaded files")
    
    # Fitur opsional jika Anda ingin script membaca langsung dari config.yaml
    parser.add_argument("--config", help="Path to config.yaml (optional)")
    parser.add_argument("--organelle", default="plastid", choices=["mito", "plastid"])
    parser.add_argument("--max_genomes", type=int, default=10)
    parser.add_argument("--min_len", type=int, default=10000)
    parser.add_argument("--max_len", type=int, default=300000)

    args = parser.parse_args()

    species = args.species
    email = args.email
    out_dir = args.output_dir

    # Jika argument species kosong, mari kita coba baca dari config.yaml (jika dideklarasikan)
    if args.config and os.path.exists(args.config):
        try:
            import yaml
            with open(args.config, 'r') as f:
                config_data = yaml.safe_load(f)
                downstream = config_data.get('downstream', {})
                if not species: species = downstream.get('species_name')
                if not email: email = downstream.get('email')
        except ImportError:
            log.warning("PyYAML tidak terinstall. Mengabaikan --config.")

    # Validasi akhir
    if not species or not email or not out_dir:
        parser.error("Species, email, dan output_dir harus diisi! (Bisa via argumen posisi atau Snakemake params)")

    log.info(f"Memulai pencarian NCBI untuk spesies: {species}")

    saved = expand_search(
        species=species,
        organelle=args.organelle,
        email=email,
        min_len=args.min_len,
        max_len=args.max_len,
        max_results=args.max_genomes,
        output_dir=out_dir,
    )
    log.info(f"Downloaded {len(saved)} reference genome(s) to {out_dir}")


if __name__ == "__main__":
    main()
