# rules/liftoff.smk – Liftoff: reference-based annotation lift-over

rule liftoff_annotate:
    """
    Lift annotation from a reference organelle genome onto the target assembly.
    Requires reference FASTA + GFF3 from a closely related species.
    """
    input:
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        done     = touch(OUTDIR + "/{sample}/liftoff/{sample}.done"),
        gff      = OUTDIR + "/{sample}/liftoff/{sample}.gff",
        unmapped = OUTDIR + "/{sample}/liftoff/{sample}.unmapped.txt",
    params:
        out_dir      = OUTDIR + "/{sample}/liftoff",
        ref_fasta    = lambda wc: samples_df.loc[wc.sample, "reference_fasta"].strip(),
        ref_gff      = lambda wc: samples_df.loc[wc.sample, "reference_gff"].strip(),
        min_coverage = config["liftoff"]["min_coverage"],
        min_identity = config["liftoff"]["min_identity"],
        extra        = config["liftoff"].get("extra", ""),
    log:
        OUTDIR + "/{sample}/logs/liftoff.log",
    threads:
        config["resources"]["liftoff"]["threads"]
    resources:
        mem_mb  = config["resources"]["liftoff"]["mem_mb"],
        runtime = config["resources"]["liftoff"]["runtime"],
    conda:
        "../envs/liftoff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        if [ -z "{params.ref_fasta}" ] || [ ! -f "{params.ref_fasta}" ]; then
            echo "ERROR: Liftoff requires a reference FASTA (reference_fasta in samples.tsv)." >&2
            touch {output.gff} {output.unmapped}
            exit 0
        fi
        if [ -z "{params.ref_gff}" ] || [ ! -f "{params.ref_gff}" ]; then
            echo "ERROR: Liftoff requires a reference GFF (reference_gff in samples.tsv)." >&2
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

        touch {output.gff} {output.unmapped}
        """


rule liftoff_to_gb:
    """
    Convert Liftoff GFF3 output to GenBank format for downstream tools (OGDraw).
    """
    input:
        gff   = OUTDIR + "/{sample}/liftoff/{sample}.gff",
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        gb = OUTDIR + "/{sample}/liftoff/{sample}.gb",
    log:
        OUTDIR + "/{sample}/logs/liftoff_to_gb.log",
    conda:
        "../envs/liftoff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})

        python3 -c "
from Bio import SeqIO
from BCBio import GFF
import sys

gff_path = '{input.gff}'
fasta_path = '{input.fasta}'
gb_path = '{output.gb}'

if __import__('os').path.getsize(gff_path) == 0:
    open(gb_path, 'w').close()
    sys.exit(0)

# Read FASTA sequences
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))

# Parse GFF and add features to sequences
with open(gff_path) as gff_fh:
    for rec in GFF.parse(gff_fh, base_dict=seq_dict):
        rec.annotations['molecule_type'] = 'DNA'
        if not rec.annotations.get('topology'):
            rec.annotations['topology'] = 'circular'

# Write GenBank
records = list(seq_dict.values())
for r in records:
    r.annotations.setdefault('molecule_type', 'DNA')
SeqIO.write(records, gb_path, 'genbank')
" 2>&1 | tee {log}

        touch {output.gb}
        """
