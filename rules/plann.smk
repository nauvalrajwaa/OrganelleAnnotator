# =============================================================================
# rules/plann.smk – Plann: reference-based plastid genome annotator
# =============================================================================
# Plann transfers annotations from a reference GenBank plastid genome to a
# new assembly using BLAST alignments. Ideal for closely related species.
# Conda: bioconda::plann
# Reference: Huang & Cronk (2015) doi:10.1186/s12859-015-0472-1
# =============================================================================

rule plann_annotate:
    """
    Annotate a plastid genome via reference-based transfer using Plann.
    Requires a reference GenBank file from a closely related species.
    Produces GenBank and tbl outputs.
    """
    input:
        fasta=lambda wc: get_fasta(wc.sample),
    output:
        done=touch(f"{OUTDIR}/plann/{{sample}}/{{sample}}.done"),
        gb=f"{OUTDIR}/plann/{{sample}}/{{sample}}.gb",
        gff=f"{OUTDIR}/plann/{{sample}}/{{sample}}.gff",
    params:
        out_dir=lambda wc: f"{OUTDIR}/plann/{wc.sample}",
        reference_gb=config["plann"]["reference_gb"],
        extra=config["plann"].get("extra", ""),
    log:
        f"{OUTDIR}/logs/plann/{{sample}}.log",
    threads:
        config["resources"]["plann"]["threads"]
    resources:
        mem_mb=config["resources"]["plann"]["mem_mb"],
        runtime=config["resources"]["plann"]["runtime"],
    conda:
        "../envs/plann.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        if [ -z "{params.reference_gb}" ] || [ ! -f "{params.reference_gb}" ]; then
            echo "ERROR: Plann requires a reference GenBank file (plann.reference_gb in config)." >&2
            echo "Provide a GenBank file from a closely related plastid genome." >&2
            touch {output.gb} {output.gff}
            exit 0
        fi

        # Run Plann
        plann \
            -reference {params.reference_gb} \
            -genome {input.fasta} \
            -out {params.out_dir}/{wildcards.sample} \
            {params.extra} \
            2>&1 | tee {log}

        # Plann outputs: <prefix>.gb, <prefix>.tbl, <prefix>.fsa
        # Locate and standardise the GenBank output
        found_gb=$(find {params.out_dir} -maxdepth 1 -name "*.gb" -o -name "*.gbk" 2>/dev/null | head -1)
        if [ -n "$found_gb" ] && [ "$found_gb" != "{output.gb}" ]; then
            cp "$found_gb" {output.gb}
        fi

        # Convert GenBank to GFF3 using a simple parser
        python3 -c "
import re, sys, os

gff_lines = ['##gff-version 3']
gb_path = '{output.gb}'
if not os.path.exists(gb_path) or os.path.getsize(gb_path) == 0:
    with open('{output.gff}', 'w') as out:
        out.write('\n'.join(gff_lines) + '\n')
    sys.exit(0)

seq_id = '{wildcards.sample}'
with open(gb_path) as fh:
    for line in fh:
        # Extract LOCUS name
        if line.startswith('LOCUS'):
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[1]
        # Extract gene annotations
        m = re.match(r'\s+(gene|CDS|tRNA|rRNA)\s+(?:complement\()?(\d+)\.\.(\d+)', line)
        if m:
            ftype, start, end = m.group(1), m.group(2), m.group(3)
            strand = '-' if 'complement' in line else '+'
            gff_lines.append(seq_id+'\tplann\t'+ftype+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tID='+ftype+'_'+start)

with open('{output.gff}', 'w') as out:
    out.write('\n'.join(gff_lines) + '\n')
" 2>>{log}

        # Ensure outputs exist
        touch {output.gb} {output.gff}
        """
