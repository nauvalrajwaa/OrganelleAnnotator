"""
Convert a GFF3 file + FASTA → GenBank format using BioPython.
Called by Snakemake; uses snakemake.input / snakemake.output / snakemake.log.
"""
import sys
import logging
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def setup_logging(log_path):
    logging.basicConfig(
        filename=log_path,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    return logging.getLogger(__name__)


def load_fasta(fasta_path: str) -> dict[str, SeqRecord]:
    """Return {seq_id: SeqRecord} from a FASTA file."""
    return SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))


def annotate_records(gff_path: str, seq_dict: dict) -> list[SeqRecord]:
    """
    Parse GFF3 and merge features onto the matching SeqRecord objects.
    Returns a list of annotated SeqRecord objects.
    """
    records = []
    with open(gff_path) as gff_handle:
        for rec in GFF.parse(gff_handle, base_dict=seq_dict):
            rec.annotations["molecule_type"] = "DNA"
            records.append(rec)
    return records


def write_genbank(records: list[SeqRecord], out_path: str):
    """Write annotated records to a GenBank file."""
    with open(out_path, "w") as out_handle:
        SeqIO.write(records, out_handle, "genbank")


def main():
    log = setup_logging(snakemake.log[0])  # noqa: F821

    gff_path   = snakemake.input.gff        # noqa: F821
    fasta_path = snakemake.input.fasta      # noqa: F821
    out_path   = snakemake.output.gb        # noqa: F821

    log.info(f"Loading FASTA: {fasta_path}")
    seq_dict = load_fasta(fasta_path)

    log.info(f"Parsing GFF3: {gff_path}")
    records = annotate_records(gff_path, seq_dict)

    if not records:
        raise ValueError(
            f"No records produced — check that sequence IDs in {gff_path} "
            f"match those in {fasta_path}."
        )

    log.info(f"Writing GenBank: {out_path} ({len(records)} record(s))")
    write_genbank(records, out_path)
    log.info("Done.")


if __name__ == "__main__":
    main()