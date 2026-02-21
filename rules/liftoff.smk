# =============================================================================
# rules/liftoff.smk – Liftoff: reference-based annotation lift-over
# =============================================================================
# Liftoff uses minimap2 to map annotations from a reference genome (GFF + FASTA)
# to a target genome. Works for both chloroplast and mitochondrial genomes.
# Especially useful when a well-annotated relative is available.
# Conda: bioconda::liftoff
# Reference: Shumate & Salzberg (2021) doi:10.1093/bioinformatics/btaa1016
# =============================================================================

rule liftoff_annotate:
    """
    Transfer annotations from a reference organelle genome to the target
    assembly using Liftoff (minimap2-based lift-over).
    Requires reference FASTA + GFF3 from a closely related species.
    """
    input:
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        done=touch(f"{OUTDIR}/liftoff/{{sample}}/{{sample}}.done"),
        gff=f"{OUTDIR}/liftoff/{{sample}}/{{sample}}.gff",
        unmapped=f"{OUTDIR}/liftoff/{{sample}}/{{sample}}.unmapped.txt",
    params:
        out_dir=lambda wc: f"{OUTDIR}/liftoff/{wc.sample}",
        ref_fasta=config["liftoff"]["reference_fasta"],
        ref_gff=config["liftoff"]["reference_gff"],
        min_coverage=config["liftoff"]["min_coverage"],
        min_identity=config["liftoff"]["min_identity"],
        extra=config["liftoff"].get("extra", ""),
    log:
        f"{OUTDIR}/logs/liftoff/{{sample}}.log",
    threads:
        config["resources"]["liftoff"]["threads"]
    resources:
        mem_mb=config["resources"]["liftoff"]["mem_mb"],
        runtime=config["resources"]["liftoff"]["runtime"],
    conda:
        "../envs/liftoff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        if [ -z "{params.ref_fasta}" ] || [ ! -f "{params.ref_fasta}" ]; then
            echo "ERROR: Liftoff requires a reference FASTA (liftoff.reference_fasta in config)." >&2
            touch {output.gff} {output.unmapped}
            exit 0
        fi
        if [ -z "{params.ref_gff}" ] || [ ! -f "{params.ref_gff}" ]; then
            echo "ERROR: Liftoff requires a reference GFF (liftoff.reference_gff in config)." >&2
            touch {output.gff} {output.unmapped}
            exit 0
        fi

        liftoff \
            -g {params.ref_gff} \
            -o {output.gff} \
            -u {output.unmapped} \
            -p {threads} \
            -s {params.min_coverage} \
            -a {params.min_identity} \
            -dir {params.out_dir}/intermediate \
            {params.extra} \
            {input.fasta} \
            {params.ref_fasta} \
            2>&1 | tee {log}

        # Ensure outputs exist
        touch {output.gff} {output.unmapped}
        """


rule liftoff_to_gb:
    """
    Convert Liftoff GFF3 output to GenBank format for downstream tools (OGDraw).
    Uses a lightweight Python conversion.
    """
    input:
        gff=f"{OUTDIR}/liftoff/{{sample}}/{{sample}}.gff",
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        gb=f"{OUTDIR}/liftoff/{{sample}}/{{sample}}.gb",
    log:
        f"{OUTDIR}/logs/liftoff/{{sample}}_to_gb.log",
    conda:
        "../envs/liftoff.yaml"
    shell:
        r"""
        set -euo pipefail

        # Use BioPython if available for proper GFF→GenBank conversion
        python3 -c "
import sys, os

gff_path = '{input.gff}'
fasta_path = '{input.fasta}'
gb_path = '{output.gb}'

if os.path.getsize(gff_path) == 0:
    open(gb_path, 'w').close()
    sys.exit(0)

try:
    from BCBio import GFF
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    # Parse FASTA
    records = dict()
    for rec in SeqIO.parse(fasta_path, 'fasta'):
        records[rec.id] = rec

    # Parse GFF and attach to records
    with open(gff_path) as gff_fh:
        for rec in GFF.parse(gff_fh, base_dict=records):
            pass  # features are attached to records

    with open(gb_path, 'w') as out:
        SeqIO.write(records.values(), out, 'genbank')
except ImportError:
    # Fallback: just touch the file
    print('WARNING: bcbio-gff not available; GenBank conversion skipped.', file=sys.stderr)
    open(gb_path, 'w').close()
" 2>{log}

        touch {output.gb}
        """
