#!/usr/bin/env python3
"""
Robust Ka/Ks Analysis Script (Universal MUSCLE Support)
Features:
- Supports both MUSCLE v3 (-in -out) and v5 (-align -output)
- Extracts CDS from Reference (.gbk) & Sample (.fas)
- Fail-safe: Always generates output TSV to prevent Snakemake crash
"""
import sys
import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil

# --- CONFIGURATION ---
GENE_SYNONYMS = {
    'cox1': ['coi', 'coxi'], 'cox2': ['coii'], 'cox3': ['coiii'],
    'cob': ['cytb'], 'nad1': ['nd1'], 'nad2': ['nd2'], 'nad3': ['nd3'],
    'nad4': ['nd4'], 'nad4l': ['nd4l'], 'nad5': ['nd5'], 'nad6': ['nd6'],
    'atp6': ['atpase6'], 'atp8': ['atpase8']
}

def get_tool_path(tool_name):
    return shutil.which(tool_name)

def clean_gene_name(name):
    """Normalize gene names to lowercase standard"""
    n = name.lower().replace('_', '').replace('-', '')
    for std, syns in GENE_SYNONYMS.items():
        if n == std or n in syns:
            return std
    return n

def parse_reference_gbk(gbk_file):
    """Extract CDS sequences from GenBank file"""
    cds_dict = {}
    if not os.path.exists(gbk_file): return cds_dict
    
    for rec in SeqIO.parse(gbk_file, "genbank"):
        for feature in rec.features:
            if feature.type == "CDS":
                gene = feature.qualifiers.get("gene", ["unknown"])[0]
                try:
                    seq = feature.extract(rec.seq)
                    clean_name = clean_gene_name(gene)
                    if clean_name not in ["unknown", "rrna", "trna"]:
                        cds_dict[clean_name] = str(seq)
                except:
                    pass 
    return cds_dict

def find_sample_genes(fasta_file, ref_genes):
    """Attempt to find genes in sample fasta via Header Matching"""
    found_genes = {}
    if not os.path.exists(fasta_file): return found_genes
    
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    for rec in records:
        header = rec.description.lower()
        # Clean ID just in case
        clean_id = clean_gene_name(rec.id)
        
        for gene in ref_genes:
            # Check 1: Exact ID match
            if clean_id == gene:
                found_genes[gene] = str(rec.seq)
                break
            # Check 2: Substring in description (e.g. "tarsius gene=cox1")
            # Boundary check to avoid finding 'nad1' inside 'nad11'
            if f" {gene} " in f" {header} " or f"_{gene}" in header or f"={gene}" in header:
                found_genes[gene] = str(rec.seq)
                break
                
    return found_genes

def run_muscle_universal(seq1, seq2, gene_name, out_dir):
    """Aligns sequences handling both MUSCLE v3 and v5 syntax"""
    muscle = get_tool_path("muscle")
    if not muscle: return None
    
    in_file = os.path.join(out_dir, f"{gene_name}.unaligned.fas")
    out_file = os.path.join(out_dir, f"{gene_name}.aligned.fas")
    
    # Write temp unaligned
    records = [
        SeqRecord(Seq(seq1), id="Reference", description=""),
        SeqRecord(Seq(seq2), id="Sample", description="")
    ]
    SeqIO.write(records, in_file, "fasta")
    
    # Strategy: Try v3 syntax first, if error, try v5
    
    # Command v3: muscle -in in.fa -out out.fa -quiet
    cmd_v3 = [muscle, "-in", in_file, "-out", out_file, "-quiet"]
    
    # Command v5: muscle -align in.fa -output out.fa
    cmd_v5 = [muscle, "-align", in_file, "-output", out_file]
    
    success = False
    
    # Try v3
    try:
        subprocess.run(cmd_v3, check=True, stderr=subprocess.PIPE, timeout=60)
        success = True
    except subprocess.CalledProcessError:
        # If v3 fails (likely "Unknown option"), try v5
        try:
            # v5 is often noisier, we capture output
            subprocess.run(cmd_v5, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=60)
            success = True
        except Exception as e:
            print(f"  > Alignment failed for {gene_name}: {e}")
    except Exception as e:
        print(f"  > General error for {gene_name}: {e}")
        
    # Clean input
    if os.path.exists(in_file): os.remove(in_file)
    
    return out_file if success and os.path.exists(out_file) else None

def convert_to_axt(aligned_fasta, axt_file):
    """Convert aligned FASTA to AXT format for KaKs_Calculator"""
    try:
        recs = list(SeqIO.parse(aligned_fasta, "fasta"))
        if len(recs) < 2: return False
        
        s1 = str(recs[0].seq)
        s2 = str(recs[1].seq)
        
        with open(axt_file, 'w') as f:
            f.write(f"{recs[0].id}-{recs[1].id}\n")
            f.write(f"{s1}\n{s2}\n")
        return True
    except:
        return False

def run_kaks_calculator(axt_file, out_file):
    kaks_exe = get_tool_path("KaKs_Calculator")
    if not kaks_exe: return False
    
    cmd = [kaks_exe, "-i", axt_file, "-o", out_file, "-m", "NG"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except:
        return False

# --- MAIN ---

if __name__ == "__main__":
    # Args: <sample_fas> <ref_gbk> <out_align_dir> <out_summary_tsv>
    if len(sys.argv) != 5:
        # Create dummy if args missing just to handle strange edge cases
        sys.exit(1)
        
    sample_fas = sys.argv[1]
    ref_gbk = sys.argv[2]
    out_dir = sys.argv[3]
    out_tsv = sys.argv[4]
    
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    
    print("--- Starting Robust Ka/Ks Analysis (Universal MUSCLE) ---")
    
    ref_genes = parse_reference_gbk(ref_gbk)
    print(f"Ref genes: {len(ref_genes)}")
    
    sample_genes = find_sample_genes(sample_fas, list(ref_genes.keys()))
    print(f"Sample genes matched: {len(sample_genes)}")
    
    common = sorted(list(set(ref_genes.keys()) & set(sample_genes.keys())))
    
    results = []
    kaks_available = get_tool_path("KaKs_Calculator") is not None
    
    if not common:
        print("WARNING: No common genes found.")
        results.append({"Gene": "None", "Status": "No_Common_Genes"})
    else:
        for gene in common:
            print(f"Processing {gene}...")
            
            # 1. Align
            aligned_file = run_muscle_universal(ref_genes[gene], sample_genes[gene], gene, out_dir)
            
            gene_result = {
                "Gene": gene, 
                "Ka": "NA", "Ks": "NA", "Ka/Ks": "NA", "P-Value": "NA", 
                "Status": "Aligned"
            }
            
            if aligned_file:
                # 2. KaKs
                if kaks_available:
                    axt_file = aligned_file.replace(".aligned.fas", ".axt")
                    kaks_out = aligned_file + ".kaks"
                    
                    if convert_to_axt(aligned_file, axt_file):
                        if run_kaks_calculator(axt_file, kaks_out):
                            try:
                                df = pd.read_csv(kaks_out, sep='\t')
                                res = df.iloc[0]
                                gene_result.update({
                                    "Ka": res.get("Ka", "NA"),
                                    "Ks": res.get("Ks", "NA"),
                                    "Ka/Ks": res.get("Ka/Ks", "NA"),
                                    "P-Value": res.get("P-Value(Fisher)", "NA"),
                                    "Status": "Success"
                                })
                            except:
                                gene_result["Status"] = "Parse_Error"
                        else:
                            gene_result["Status"] = "Calc_Failed"
                else:
                    gene_result["Status"] = "Tool_Missing"
            else:
                gene_result["Status"] = "Align_Failed"
            
            results.append(gene_result)

    # ALWAYS write the output file to satisfy Snakemake
    if not results:
        results.append({"Gene": "Error", "Status": "Script_Failed"})
        
    df_out = pd.DataFrame(results)
    df_out.to_csv(out_tsv, sep='\t', index=False)
    print(f"Analysis complete. Results written to {out_tsv}")
    
    sys.exit(0)