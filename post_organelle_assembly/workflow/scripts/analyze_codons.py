from Bio import SeqIO
import sys
import pandas as pd

def analyze_codons(fasta_file, output_file, genetic_code=2):
    valid_genes = ['atp6', 'atp8', 'cox1', 'cox2', 'cox3', 'cob', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    
    records = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        name = record.id
        # Simple filter
        is_cds = any(g in name.lower() for g in valid_genes) and not ("trn" in name.lower() or "rrn" in name.lower())
        
        if is_cds:
            seq = str(record.seq).upper().replace("-", "")
            if len(seq) < 3: continue
            
            start_codon = seq[:3]
            stop_codon = seq[-3:]
            
            # Check length check
            remainder = len(seq) % 3
            
            records.append({
                'Gene': name,
                'Length': len(seq),
                'Start_Codon': start_codon,
                'Stop_Codon': stop_codon,
                'DivisibleBy3': remainder == 0
            })
            
    df = pd.DataFrame(records)
    
    # Calculate frequencies
    start_counts = df['Start_Codon'].value_counts()
    stop_counts = df['Stop_Codon'].value_counts()
    
    with open(output_file, 'w') as f:
        f.write(f"Analyzed {len(df)} CDS sequences.\n\n")
        f.write("--- Start Codon Usage ---\n")
        f.write(start_counts.to_string())
        f.write("\n\n--- Stop Codon Usage ---\n")
        f.write(stop_counts.to_string())
        f.write("\n\n--- Detailed Table ---\n")
        f.write(df.to_string())
        
if __name__ == "__main__":
    analyze_codons(sys.argv[1], sys.argv[2])
