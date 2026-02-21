#!/usr/bin/env python3
"""
gene_utils.py – Shared gene name maps and utilities for downstream analysis.

Centralises gene synonym resolution, standard gene sets, and common helpers
so they are not duplicated across multiple analysis scripts.
"""

from collections import OrderedDict

# =============================================================================
# Standard mitochondrial protein-coding genes (vertebrate + plant names)
# =============================================================================
MITO_GENES_VERTEBRATE = [
    "atp6", "atp8", "cox1", "cox2", "cox3",
    "cob", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6",
]

MITO_GENES_PLANT = MITO_GENES_VERTEBRATE + [
    "atp1", "atp4", "atp9",
    "ccmB", "ccmC", "ccmFc", "ccmFn",
    "matR", "mttB",
    "rpl2", "rpl5", "rpl10", "rpl16",
    "rps1", "rps2", "rps3", "rps4", "rps7", "rps10", "rps11", "rps12",
    "rps13", "rps14", "rps19",
    "sdh3", "sdh4",
]

# =============================================================================
# Standard chloroplast protein-coding genes
# =============================================================================
PLASTID_GENES = [
    "accD", "atpA", "atpB", "atpE", "atpF", "atpH", "atpI",
    "ccsA", "cemA", "clpP",
    "infA",
    "matK",
    "ndhA", "ndhB", "ndhC", "ndhD", "ndhE", "ndhF", "ndhG", "ndhH", "ndhI",
    "ndhJ", "ndhK",
    "petA", "petB", "petD", "petG", "petL", "petN",
    "psaA", "psaB", "psaC", "psaI", "psaJ",
    "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI",
    "psbJ", "psbK", "psbL", "psbM", "psbN", "psbT", "psbZ",
    "rbcL", "rpoA", "rpoB", "rpoC1", "rpoC2",
    "rpl2", "rpl14", "rpl16", "rpl20", "rpl22", "rpl23", "rpl32",
    "rpl33", "rpl36",
    "rps2", "rps3", "rps4", "rps7", "rps8", "rps11", "rps12", "rps14",
    "rps15", "rps16", "rps18", "rps19",
    "ycf1", "ycf2", "ycf3", "ycf4",
]

# =============================================================================
# Gene synonym map – normalises alternative names to a canonical form
# =============================================================================
GENE_SYNONYMS: dict[str, str] = {
    # Cytochrome oxidase
    "COX1": "cox1", "COI": "cox1", "co1": "cox1", "COX I": "cox1",
    "COX2": "cox2", "COII": "cox2", "co2": "cox2", "COX II": "cox2",
    "COX3": "cox3", "COIII": "cox3", "co3": "cox3", "COX III": "cox3",
    # Cytochrome b
    "CYTB": "cob", "cytb": "cob", "Cyt b": "cob", "CYB": "cob",
    "COB": "cob", "Cytb": "cob",
    # ATPase
    "ATP6": "atp6", "ATPase6": "atp6", "ATPase 6": "atp6",
    "ATP8": "atp8", "ATPase8": "atp8", "ATPase 8": "atp8",
    # NADH dehydrogenase
    "ND1": "nad1", "NAD1": "nad1", "NADH1": "nad1",
    "ND2": "nad2", "NAD2": "nad2", "NADH2": "nad2",
    "ND3": "nad3", "NAD3": "nad3", "NADH3": "nad3",
    "ND4": "nad4", "NAD4": "nad4", "NADH4": "nad4",
    "ND4L": "nad4l", "NAD4L": "nad4l", "NADH4L": "nad4l",
    "ND5": "nad5", "NAD5": "nad5", "NADH5": "nad5",
    "ND6": "nad6", "NAD6": "nad6", "NADH6": "nad6",
    # rRNA
    "12S": "rrnS", "s-rRNA": "rrnS", "s_rRNA": "rrnS",
    "rns": "rrnS", "rrn12": "rrnS", "12S rRNA": "rrnS",
    "16S": "rrnL", "l-rRNA": "rrnL", "l_rRNA": "rrnL",
    "rnl": "rrnL", "rrn16": "rrnL", "16S rRNA": "rrnL",
}


def normalise_gene_name(raw: str) -> str:
    """Resolve a gene name to its canonical lowercase form.

    Strips common suffixes/prefixes, deduplicates case variants, and checks
    the synonym table.  Returns the canonical form or the lowercased input
    if no synonym is found.
    """
    name = raw.strip()
    # Direct match in synonyms table
    if name in GENE_SYNONYMS:
        return GENE_SYNONYMS[name]
    # Case-insensitive match
    for key, val in GENE_SYNONYMS.items():
        if name.lower() == key.lower():
            return val
    # Strip common annotations like "(gene)" or numbers after underscore
    base = name.split("(")[0].strip().split("_exon")[0].strip()
    if base in GENE_SYNONYMS:
        return GENE_SYNONYMS[base]
    return name.lower()


def get_gene_set(organelle: str) -> list[str]:
    """Return the standard gene set for the given organelle type."""
    if organelle in ("plastid", "chloroplast", "cp"):
        return PLASTID_GENES
    elif organelle in ("mito", "mitochondria", "mt"):
        return MITO_GENES_VERTEBRATE
    elif organelle in ("plant_mito", "plant_mt"):
        return MITO_GENES_PLANT
    return MITO_GENES_VERTEBRATE + PLASTID_GENES


# =============================================================================
# Genetic code tables
# =============================================================================
GENETIC_CODE_NAMES = {
    1: "Standard",
    2: "Vertebrate Mitochondrial",
    3: "Yeast Mitochondrial",
    4: "Mold, Protozoan, Coelenterate Mitochondrial",
    5: "Invertebrate Mitochondrial",
    11: "Bacterial / Plant Plastid",
}


def validate_genetic_code(code: int) -> int:
    """Validate the genetic code number; raise ValueError if unsupported."""
    valid = {1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25}
    if code not in valid:
        raise ValueError(f"Unsupported genetic code: {code}. Valid: {sorted(valid)}")
    return code
