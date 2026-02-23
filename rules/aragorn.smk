# rules/aragorn.smk – Aragorn: tRNA and tmRNA detection

rule aragorn_annotate:
    """
    Detect tRNA and tmRNA genes using Aragorn.
    Produces tabular/text output and a GFF3 conversion.
    """
    input:
        fasta = lambda wc: samples_df.loc[wc.sample, "fasta"],
    output:
        done = touch(OUTDIR + "/{sample}/aragorn/{sample}.done"),
        txt  = OUTDIR + "/{sample}/aragorn/{sample}.aragorn.txt",
        gff  = OUTDIR + "/{sample}/aragorn/{sample}.gff",
    params:
        out_dir              = OUTDIR + "/{sample}/aragorn",
        genetic_code         = lambda wc: samples_df.loc[wc.sample, "genetic_code"],
        topology             = config["aragorn"]["topology"],
        search_both_strands  = config["aragorn"]["search_both_strands"],
        extra                = config["aragorn"].get("extra", ""),
    log:
        OUTDIR + "/{sample}/logs/aragorn.log",
    threads:
        config["resources"]["aragorn"]["threads"]
    resources:
        mem_mb  = config["resources"]["aragorn"]["mem_mb"],
        runtime = config["resources"]["aragorn"]["runtime"],
    conda:
        "../envs/aragorn.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log})

        TOPO_FLAG=""
        case "{params.topology}" in
            circular) TOPO_FLAG="-c" ;;
            linear)   TOPO_FLAG="-l" ;;
        esac

        STRAND_FLAG=""
        if [ "{params.search_both_strands}" = "True" ]; then
            STRAND_FLAG="-d"
        fi

        aragorn \
            $TOPO_FLAG \
            $STRAND_FLAG \
            -gc{params.genetic_code} \
            -w \
            {params.extra} \
            {input.fasta} \
            -o {output.txt} \
            2>&1 | tee {log}

        # Convert Aragorn text output to GFF3
        python3 -c "
import re, sys

gff = ['##gff-version 3']
seq_id = None
with open('{output.txt}') as fh:
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            seq_id = line.split()[0].lstrip('>')
            continue
        m = re.match(r'\s*\d+\s+(tRNA-\S+|tmRNA)\s+(\\[?c?\\[?)(\\d+),(\\d+)\\]?\\)?\\s*(\\(\\w+\\))?', line)
        if not m:
            m = re.match(r'\s*\d+\s+(tRNA-\S+|tmRNA)\s+c\\[(\\d+),(\\d+)\\]\s*(\\(\\w+\\))?', line)
            if m:
                name, start, end = m.group(1), m.group(2), m.group(3)
                anticodon = m.group(4) or ''
                strand = '-'
                gff.append(seq_id+'\taragorn\ttRNA\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tID='+name+';Name='+name+';anticodon='+anticodon)
                continue
            m = re.match(r'\s*\d+\s+(tRNA-\S+|tmRNA)\s+\\[?(\\d+),(\\d+)\\]?\s*(\\(\\w+\\))?', line)
            if m:
                name, start, end = m.group(1), m.group(2), m.group(3)
                anticodon = m.group(4) or ''
                strand = '+'
                gff.append(seq_id+'\taragorn\ttRNA\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tID='+name+';Name='+name+';anticodon='+anticodon)
                continue

with open('{output.gff}', 'w') as out:
    out.write('\n'.join(gff) + '\n')
" 2>>{log}

        touch {output.txt} {output.gff}
        """
