import sys
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def analyze_composition(input_fasta, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    records = []
    aa_counts = {}
    
    # Valid protein coding genes in mitochondria (usually)
    valid_genes = ['atp6', 'atp8', 'cox1', 'cox2', 'cox3', 'cob', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    
    print(f"Reading {input_fasta}...")
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        name = record.id.lower()
        if any(g in name for g in valid_genes) and not ("trn" in name or "rrn" in name):
            seq_str = str(record.seq).upper().replace("-", "")
            
            # GC Content
            gc = gc_fraction(seq_str) * 100
            
            records.append({
                "Gene": record.id,
                "Length": len(seq_str),
                "GC_Content": gc
            })
            
            # Translation (Simple vertebrate mitochondrial)
            # We assume sequence is nucletoide. If it's already protein (rare from this pipeline flow?), this needs check.
            # Based on previous scripts, input seems to be nucleotides.
            # Using table 2 (Vertebrate Mitochondrial)
            try:
                # Pad if needed
                pad = len(seq_str) % 3
                if pad > 0:
                     seq_str_t = seq_str[:-pad]
                else:
                     seq_str_t = seq_str
                     
                # Translate
                protein = record.seq.translate(table=2, to_stop=False)
                
                for aa in str(protein):
                    if aa not in aa_counts: aa_counts[aa] = 0
                    aa_counts[aa] += 1
            except Exception as e:
                print(f"Translation warning for {name}: {e}")

    # 1. GC Content Plot
    if records:
        df_gc = pd.DataFrame(records)
        df_gc.to_csv(os.path.join(output_dir, "gc_content.tsv"), sep="\t", index=False)
        
        plt.figure(figsize=(10, 6))
        sns.barplot(data=df_gc, x="Gene", y="GC_Content", hue="Gene", palette="viridis", legend=False)
        plt.title("GC Content per Gene")
        plt.ylabel("GC Content (%)")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "gc_content_plot.png"))
        plt.close()
        print("Generated GC content plot.")
    else:
        print("No CDS records found for GC analysis.")

    # 2. AA Composition Plot
    if aa_counts:
        # Sort by AA code
        sorted_aa = sorted(aa_counts.keys())
        total_aa = sum(aa_counts.values())
        
        aa_data = []
        for aa in sorted_aa:
            if aa == "*": continue # Skip stop
            perc = (aa_counts[aa] / total_aa) * 100
            aa_data.append({"AminoAcid": aa, "Percentage": perc, "Count": aa_counts[aa]})
            
        df_aa = pd.DataFrame(aa_data)
        df_aa.to_csv(os.path.join(output_dir, "aa_composition.tsv"), sep="\t", index=False)
        
        plt.figure(figsize=(12, 6))
        sns.barplot(data=df_aa, x="AminoAcid", y="Percentage", hue="AminoAcid", palette="magma", legend=False)
        plt.title(f"Amino Acid Composition (Total AA: {total_aa})")
        plt.ylabel("Percentage (%)")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "aa_composition_plot.png"))
        plt.close()
        print("Generated AA composition plot.")
    else:
        print("No Amino Acid data generated.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python analyze_composition.py input.fasta output_dir")
        sys.exit(1)
        
    analyze_composition(sys.argv[1], sys.argv[2])
