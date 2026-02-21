#!/usr/bin/env python3
"""
Synteny Analysis Pro (Fixed for Snakemake & MUMmer)
Features:
- FASTA Sanitization (Fixes delta-filter error 400)
- Robust MUMmer Parsing
- Fail-safe Outputs (Prevents Snakemake MissingOutputException)
- Bezier Ribbons Visualization
"""
import sys
import os
import subprocess
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from Bio import SeqIO
import shutil
import tempfile

# --- 1. UTILITIES ---

def sanitize_fasta(input_fasta, output_fasta, generic_name):
    """
    MUMmer crashes with complex headers. This rewrites input to a temp file
    with a clean header (e.g., 'sample', 'ref') to prevent 'error no: 400'.
    Returns the total length of the genome.
    """
    total_len = 0
    clean_records = []
    
    # Read original
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        return 0
    
    # If multiple contigs, we concatenate them virtually or keep distinct indices
    # For MUMmer plotting simplicity, we usually take the longest or treat as one concat
    # Here we just rename them sequentially to avoid special char issues
    
    for i, record in enumerate(records):
        # Rename header to simple alphanumeric
        record.id = f"{generic_name}_{i+1}"
        record.description = ""
        clean_records.append(record)
        total_len += len(record.seq)
        
    SeqIO.write(clean_records, output_fasta, "fasta")
    return total_len

def create_failure_outputs(out_plot, out_stats, message="Analysis Failed"):
    """Creates placeholder files so Snakemake doesn't crash."""
    print(f"Generating failure outputs: {message}")
    
    # 1. Create Image with Error Message
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=14, color='red')
    ax.axis('off')
    plt.savefig(out_plot, dpi=100, bbox_inches='tight')
    plt.close()
    
    # 2. Create Empty CSV
    df = pd.DataFrame({'status': ['failed'], 'error': [message]})
    df.to_csv(out_stats, sep='\t', index=False)

# --- 2. ALIGNMENT WORKFLOW ---

def run_mummer_workflow(sample_fasta, reference_fasta, output_prefix):
    """Runs robust MUMmer pipeline with sanitized inputs"""
    nucmer = shutil.which("nucmer")
    delta_filter = shutil.which("delta-filter")
    show_coords = shutil.which("show-coords")
    
    if not (nucmer and delta_filter and show_coords):
        print("Error: MUMmer tools missing.")
        return None
    
    # Create Temp files for sanitized inputs
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp_s, \
         tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tmp_r:
        
        # Sanitize Inputs
        s_len = sanitize_fasta(sample_fasta, tmp_s.name, "Query")
        r_len = sanitize_fasta(reference_fasta, tmp_r.name, "Ref")
        
        tmp_s_path = tmp_s.name
        tmp_r_path = tmp_r.name

    prefix = output_prefix.replace('.png', '').replace('.xmfa', '')
    delta = f"{prefix}.delta"
    filtered = f"{prefix}.filter.delta"
    coords = f"{prefix}.coords"
    
    try:
        # 1. Align (nucmer)
        # --maxmatch allows finding all matches regardless of uniqueness (good for organelles)
        subprocess.run([nucmer, "--maxmatch", "--prefix=" + prefix, tmp_r_path, tmp_s_path], 
                       check=True, capture_output=True)
        
        # 2. Filter (delta-filter)
        # -m: many-to-many, -i 85: min identity
        with open(filtered, 'w') as f:
            subprocess.run([delta_filter, "-m", "-i", "85", delta], stdout=f, check=True)
            
        # 3. Coords (show-coords)
        with open(coords, 'w') as f:
            subprocess.run([show_coords, "-r", "-c", "-l", "-T", filtered], stdout=f, check=True)
            
    except subprocess.CalledProcessError as e:
        print(f"MUMmer process error: {e}")
        # Cleanup temps
        os.remove(tmp_s_path)
        os.remove(tmp_r_path)
        return None
    except Exception as e:
        print(f"General error: {e}")
        os.remove(tmp_s_path)
        os.remove(tmp_r_path)
        return None

    # Cleanup temps
    try:
        os.remove(tmp_s_path)
        os.remove(tmp_r_path)
    except: pass
    
    return coords

# --- 3. PARSING ---

def parse_mummer_coords(coords_file):
    blocks = []
    if not os.path.exists(coords_file) or os.path.getsize(coords_file) == 0: return blocks
    
    try:
        # Read flexibly
        try: df = pd.read_csv(coords_file, sep='\t')
        except: df = pd.read_csv(coords_file, sep='\t', skiprows=4, header=None)
        
        if df.empty: return blocks

        for idx, row in df.iterrows():
            vals = row.values
            # show-coords -T output structure: [S1] [E1] [S2] [E2] ...
            # Ensure we cast to int
            ref_start, ref_end = int(vals[0]), int(vals[1])
            q_start, q_end = int(vals[2]), int(vals[3])
            
            # Determine Strand based on Query coordinates
            q_strand = '+' if q_start < q_end else '-'
            
            # Normalize Query Coords (min, max) for plotting logic
            final_q_start, final_q_end = min(q_start, q_end), max(q_start, q_end)

            blocks.append({
                'id': idx,
                'sequences': [
                    {'seq_id': 'ref', 'start': ref_start, 'end': ref_end, 'strand': '+'},
                    {'seq_id': 'query', 'start': final_q_start, 'end': final_q_end, 'strand': q_strand}
                ]
            })
    except Exception as e:
        print(f"Parse Error: {e}")
        
    return blocks

# --- 4. VISUALIZATION ---

def get_bezier_path(start_top, end_top, start_bot, end_bot, y_top, y_bot):
    """Generates a smooth ribbon path"""
    mid_y = (y_top + y_bot) / 2
    verts = [
        (start_top, y_top), (start_top, mid_y), (start_bot, mid_y), (start_bot, y_bot),
        (end_bot, y_bot), (end_bot, mid_y), (end_top, mid_y), (end_top, y_top),
        (start_top, y_top)
    ]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, 
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    return Path(verts, codes)

def draw_genome_ruler(ax, length, y_pos, label, color='#333333'):
    ax.plot([0, length], [y_pos, y_pos], color=color, linewidth=2, zorder=5)
    ax.text(-length*0.02, y_pos, label, va='center', ha='right', fontsize=12, fontweight='bold', color=color)
    
    interval = 50000
    if length < 50000: interval = 5000
    elif length < 200000: interval = 20000
    
    for i in range(0, length + 1, interval):
        tick_h = 0.05
        ax.plot([i, i], [y_pos - tick_h, y_pos + tick_h], color=color, linewidth=1, zorder=5)
        if i % (interval * 2) == 0:
            ax.text(i, y_pos + tick_h + 0.02, f"{i/1000:.0f}k", ha='center', va='bottom', fontsize=9, color=color)

def plot_synteny_ribbon(blocks, sample_len, ref_len, output_plot):
    fig, ax = plt.subplots(figsize=(15, 7))
    y_sample, y_ref = 0.8, 0.2
    
    # Draw Ribbons
    for block in blocks:
        if len(block['sequences']) < 2: continue
        ref, query = block['sequences'][0], block['sequences'][1]
        is_fwd = (ref['strand'] == query['strand'])
        
        color = '#2E86AB' if is_fwd else '#D64045' # Teal vs Red
        alpha = 0.6 if is_fwd else 0.5
        
        if is_fwd:
            path = get_bezier_path(query['start'], query['end'], ref['start'], ref['end'], y_sample - 0.02, y_ref + 0.02)
        else:
            path = get_bezier_path(query['start'], query['end'], ref['end'], ref['start'], y_sample - 0.02, y_ref + 0.02)
            
        ax.add_patch(patches.PathPatch(path, facecolor=color, edgecolor='none', lw=0, alpha=alpha))

    # Draw Rulers
    max_len = max(sample_len, ref_len)
    draw_genome_ruler(ax, sample_len, y_sample, "Sample", '#333333')
    draw_genome_ruler(ax, ref_len, y_ref, "Reference", '#333333')
    
    # Legend
    legend_elements = [patches.Patch(facecolor='#2E86AB', alpha=0.6, label='Collinear'),
                       patches.Patch(facecolor='#D64045', alpha=0.5, label='Inverted')]
    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, frameon=False)
    
    ax.set_xlim(-max_len * 0.05, max_len * 1.05)
    ax.set_ylim(0, 1.2)
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()

# --- 5. MAIN ---

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python synteny_viz.py <sample> <ref> <plot.png> <stats.txt>")
        sys.exit(1)
        
    s_fasta, r_fasta, out_plot, out_stats = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    
    # Setup Dirs
    os.makedirs(os.path.dirname(out_plot), exist_ok=True)
    
    try:
        print("--- Starting Synteny Analysis ---")
        
        # 1. Run Alignment (with sanitized inputs inside)
        coords = run_mummer_workflow(s_fasta, r_fasta, out_plot)
        
        blocks = []
        if coords:
            blocks = parse_mummer_coords(coords)
        
        if blocks:
            # 2. Calc Lengths (Original Files)
            s_len = sum(len(r) for r in SeqIO.parse(s_fasta, "fasta"))
            r_len = sum(len(r) for r in SeqIO.parse(r_fasta, "fasta"))
            
            # 3. Plot
            print(f"Generating plot for {len(blocks)} blocks...")
            plot_synteny_ribbon(blocks, s_len, r_len, out_plot)
            
            # 4. Stats
            stats = {'blocks': len(blocks), 'status': 'success'}
            pd.DataFrame([stats]).to_csv(out_stats, sep='\t', index=False)
            print("Analysis Complete.")
            
        else:
            print("WARNING: No synteny blocks found or alignment failed.")
            create_failure_outputs(out_plot, out_stats, "No Synteny Detected")

    except Exception as e:
        print(f"CRITICAL ERROR: {e}")
        # Ensure Snakemake doesn't fail due to missing outputs
        create_failure_outputs(out_plot, out_stats, f"Script Error: {str(e)[:30]}")
        # Optional: sys.exit(1) if you WANT the pipeline to stop on error
        # But for 'MissingOutputException' prevention, we exit 0 usually