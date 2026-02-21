#!/usr/bin/env python3
"""
Generate native Circos (circos.ca) configuration and data files from a GenBank
file, then run circos to produce a publication-quality circular genome map.

Tracks (from outermost to innermost):
  1. Karyotype backbone with tick marks
  2. Forward-strand genes (CDS / tRNA / rRNA) – outer ring
  3. Reverse-strand genes – inner ring
  4. GC content histogram (deviation from mean)
  5. GC skew histogram
  6. Gene labels

Usage:
    python create_circos_plot.py <gbk_file> <output_png>
"""

import sys
import os
import subprocess
import shutil
import tempfile
import textwrap
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _clean(name):
    """Sanitise a gene name for Circos labels (remove parentheses etc.)."""
    return name.replace("(", "_").replace(")", "").replace(" ", "_")


# ---------------------------------------------------------------------------
# Data generation
# ---------------------------------------------------------------------------

def generate_karyotype(record, outdir):
    """Write karyotype file.  One chromosome = the full genome."""
    genome_len = len(record.seq)
    chrom_id = "chr1"
    path = os.path.join(outdir, "karyotype.txt")
    with open(path, "w") as fh:
        fh.write(f"chr - {chrom_id} {record.id} 0 {genome_len} chr1\n")
    return chrom_id, genome_len


def generate_gene_tracks(record, chrom_id, outdir):
    """Write highlight (tile) data for genes, split by strand."""
    fp_fwd = os.path.join(outdir, "genes_forward.txt")
    fp_rev = os.path.join(outdir, "genes_reverse.txt")
    fp_labels = os.path.join(outdir, "gene_labels.txt")

    color_map = {
        "CDS":  "blue",
        "tRNA": "red",
        "rRNA": "green",
    }

    with open(fp_fwd, "w") as fwd, \
         open(fp_rev, "w") as rev, \
         open(fp_labels, "w") as lab:
        for feat in record.features:
            if feat.type not in color_map:
                continue
            start = int(feat.location.start)
            end   = int(feat.location.end)
            strand = feat.location.strand
            color  = color_map.get(feat.type, "grey")

            gene_name = ""
            if "gene" in feat.qualifiers:
                gene_name = feat.qualifiers["gene"][0]
            elif "product" in feat.qualifiers:
                gene_name = feat.qualifiers["product"][0]

            line = f"{chrom_id} {start} {end} fill_color={color}\n"
            if strand == 1:
                fwd.write(line)
            else:
                rev.write(line)

            # Labels for CDS and rRNA (tRNA names are short; label all)
            safe = _clean(gene_name) if gene_name else feat.type
            lab.write(f"{chrom_id} {start} {end} {safe}\n")


def generate_gc_content(record, chrom_id, outdir, window=500, step=100):
    """Write GC-content histogram data (deviation from mean)."""
    seq = str(record.seq).upper()
    genome_len = len(seq)
    avg_gc = gc_fraction(record.seq)

    path = os.path.join(outdir, "gc_content.txt")
    with open(path, "w") as fh:
        for i in range(0, genome_len - window, step):
            w = seq[i : i + window]
            gc = gc_fraction(w)
            val = gc - avg_gc           # deviation
            fh.write(f"{chrom_id} {i} {i + window} {val:.6f}\n")


def generate_gc_skew(record, chrom_id, outdir, window=500, step=100):
    """Write GC-skew histogram data  (G-C)/(G+C)."""
    seq = str(record.seq).upper()
    genome_len = len(seq)

    path = os.path.join(outdir, "gc_skew.txt")
    with open(path, "w") as fh:
        for i in range(0, genome_len - window, step):
            w = seq[i : i + window]
            g = w.count("G")
            c = w.count("C")
            skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
            fh.write(f"{chrom_id} {i} {i + window} {skew:.6f}\n")


# ---------------------------------------------------------------------------
# Circos configuration
# ---------------------------------------------------------------------------

def write_circos_conf(outdir, genome_len):
    """Write the main circos.conf and supporting configuration blocks.
    All file paths are relative so Circos can run from outdir (avoids spaces)."""

    # ---- ticks.conf --------------------------------------------------------
    ticks_conf = textwrap.dedent("""\
    show_ticks       = yes
    show_tick_labels = yes

    <ticks>
    radius    = 1r
    color     = black
    thickness = 2p

    <tick>
    spacing        = 1000u
    size           = 10p
    show_label     = yes
    label_size     = 20p
    label_offset   = 10p
    format         = %d
    suffix         = kb
    multiplier     = 1e-3
    </tick>

    <tick>
    spacing        = 500u
    size           = 6p
    show_label     = no
    </tick>
    </ticks>
    """)

    with open(os.path.join(outdir, "ticks.conf"), "w") as fh:
        fh.write(ticks_conf)

    # ---- ideogram.conf -----------------------------------------------------
    ideogram_conf = textwrap.dedent("""\
    <ideogram>
    <spacing>
    default = 0.005r
    break   = 0.5r
    </spacing>

    radius           = 0.90r
    thickness         = 20p
    fill              = yes
    fill_color        = grey
    stroke_color      = dgrey
    stroke_thickness  = 2p

    show_label        = yes
    label_font        = default
    label_radius      = 1r + 60p
    label_size        = 30p
    label_parallel    = yes
    </ideogram>
    """)

    with open(os.path.join(outdir, "ideogram.conf"), "w") as fh:
        fh.write(ideogram_conf)

    # ---- main circos.conf --------------------------------------------------
    # All paths are RELATIVE — circos will be invoked with cwd = outdir
    circos_conf = textwrap.dedent(f"""\
    # Circos configuration – mitochondrial genome map
    karyotype = karyotype.txt

    <<include ideogram.conf>>
    <<include ticks.conf>>

    chromosomes_units = 1
    chromosomes_display_default = yes

    # ---- image settings ---
    <image>
    <<include etc/image.conf>>
    dir*    = .
    file*   = circos.png
    radius* = 1500p
    </image>

    # ================================================================
    # HIGHLIGHTS – gene tracks (coloured tiles)
    # ================================================================
    <highlights>

    # Forward-strand genes (outer ring)
    <highlight>
    file   = genes_forward.txt
    r0     = 0.88r
    r1     = 0.95r
    </highlight>

    # Reverse-strand genes (next ring inwards)
    <highlight>
    file   = genes_reverse.txt
    r0     = 0.80r
    r1     = 0.87r
    </highlight>

    </highlights>

    # ================================================================
    # PLOTS – histograms for GC content & GC skew, plus gene labels
    # ================================================================
    <plots>

    # -- Gene labels --
    <plot>
    type       = text
    file       = gene_labels.txt
    color      = black
    r0         = 0.95r
    r1         = 1.15r
    label_size = 18p
    label_font = condensed
    label_parallel = yes

    <rules>
    <rule>
    condition  = 1
    label_size = 16p
    </rule>
    </rules>
    </plot>

    # -- GC content (deviation from mean) --
    <plot>
    type       = histogram
    file       = gc_content.txt
    r0         = 0.60r
    r1         = 0.78r
    min        = -0.15
    max        = 0.15
    color      = black
    fill_under = yes
    fill_color = blue
    thickness  = 1

    <backgrounds>
    <background>
    color = vvlgrey
    </background>
    </backgrounds>

    <axes>
    <axis>
    color     = lgrey
    thickness = 1
    spacing   = 0.05r
    </axis>
    </axes>

    <rules>
    <rule>
    condition  = var(value) < 0
    fill_color = orange
    </rule>
    </rules>
    </plot>

    # -- GC skew (G-C)/(G+C) --
    <plot>
    type       = histogram
    file       = gc_skew.txt
    r0         = 0.40r
    r1         = 0.58r
    min        = -0.30
    max        = 0.30
    color      = black
    fill_under = yes
    fill_color = dgreen
    thickness  = 1

    <backgrounds>
    <background>
    color = vvlgrey
    </background>
    </backgrounds>

    <axes>
    <axis>
    color     = lgrey
    thickness = 1
    spacing   = 0.05r
    </axis>
    </axes>

    <rules>
    <rule>
    condition  = var(value) < 0
    fill_color = purple
    </rule>
    </rules>
    </plot>

    </plots>

    <<include etc/colors_fonts_patterns.conf>>
    <<include etc/housekeeping.conf>>
    """)

    conf_path = os.path.join(outdir, "circos.conf")
    with open(conf_path, "w") as fh:
        fh.write(circos_conf)

    return conf_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(gbk_file, output_png):
    # Parse GenBank
    gbk_file = os.path.abspath(gbk_file)
    output_png = os.path.abspath(output_png)

    print(f"Parsing GenBank: {gbk_file}")
    record = SeqIO.read(gbk_file, "genbank")
    genome_len = len(record.seq)
    avg_gc = gc_fraction(record.seq)
    print(f"Genome length: {genome_len:,} bp  |  GC: {avg_gc:.1%}")

    # Use a temp directory WITHOUT spaces (circos cannot handle spaces in paths)
    tmpdir = tempfile.mkdtemp(prefix="circos_")
    print(f"Working directory: {tmpdir}")

    try:
        # Generate data files
        chrom_id, glen = generate_karyotype(record, tmpdir)
        generate_gene_tracks(record, chrom_id, tmpdir)
        generate_gc_content(record, chrom_id, tmpdir)
        generate_gc_skew(record, chrom_id, tmpdir)

        # Write Circos configuration (all paths relative)
        conf_path = write_circos_conf(tmpdir, glen)
        print(f"Circos config written to {conf_path}")

        # Run Circos from inside tmpdir so relative paths resolve
        circos_exe = shutil.which("circos")
        if not circos_exe:
            print("ERROR: circos executable not found in PATH.")
            print("Install via: conda install -c bioconda circos")
            sys.exit(1)

        print("Running circos ...")
        
        # Build environment with PERL5LIB pointing to the conda env's
        # perl lib dirs.  This is needed when the conda env path contains
        # spaces (circos uses #!/usr/bin/env perl which may resolve to
        # the system perl instead of the conda perl).
        env = os.environ.copy()
        conda_prefix = env.get("CONDA_PREFIX", "")
        if conda_prefix:
            import glob
            perl_lib_dirs = glob.glob(os.path.join(conda_prefix, "lib", "perl5", "**"), recursive=False)
            perl_lib_dirs += glob.glob(os.path.join(conda_prefix, "lib", "perl5", "site_perl", "**"), recursive=False)
            # Also add the base perl5 dir itself
            perl_lib_dirs.insert(0, os.path.join(conda_prefix, "lib", "perl5", "site_perl"))
            perl_lib_dirs.insert(0, os.path.join(conda_prefix, "lib", "perl5"))
            env["PERL5LIB"] = ":".join(perl_lib_dirs)
            print(f"CONDA_PREFIX: {conda_prefix}")
            print(f"PERL5LIB set with {len(perl_lib_dirs)} paths")
        
        # Use the conda perl explicitly if available
        conda_perl = os.path.join(conda_prefix, "bin", "perl") if conda_prefix else None
        if conda_perl and os.path.isfile(conda_perl):
            perl_exe = conda_perl
        else:
            perl_exe = shutil.which("perl")
        
        # Get circos script path (not the wrapper, the actual perl script)
        circos_bin = os.path.join(conda_prefix, "bin", "circos") if conda_prefix else circos_exe

        result = subprocess.run(
            [perl_exe, circos_bin, "-conf", "circos.conf", "-nosvg"],
            capture_output=True, text=True, timeout=300,
            cwd=tmpdir,
            env=env,
        )

        # Print circos debug output
        if result.stdout:
            print(result.stdout[-2000:])
        if result.stderr:
            print(result.stderr[-2000:])

        if result.returncode != 0:
            print(f"Circos exited with code {result.returncode}")
            sys.exit(1)

        # Copy output PNG to final destination
        circos_png = os.path.join(tmpdir, "circos.png")
        if not os.path.exists(circos_png):
            print(f"ERROR: circos did not produce {circos_png}")
            print("Files in tmpdir:", os.listdir(tmpdir))
            sys.exit(1)

        os.makedirs(os.path.dirname(output_png), exist_ok=True)
        shutil.copy2(circos_png, output_png)
        print(f"Circos plot saved: {output_png}")

        # Also keep data files next to the output for reference
        data_copy = os.path.join(os.path.dirname(output_png), "circos_data")
        if os.path.exists(data_copy):
            shutil.rmtree(data_copy)
        shutil.copytree(tmpdir, data_copy)

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_circos_plot.py <gbk_file> <output_png>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
