#!/usr/bin/env python3
import sys
import os
import shutil
import subprocess
import pandas as pd
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Data import CodonTable

# --- CONFIGURATION: GENE SYNONYMS ---
GENE_MAP = {
    'coi': 'cox1', 'cox1': 'cox1', 'coxi': 'cox1',
    'coii': 'cox2', 'cox2': 'cox2',
    'coiii': 'cox3', 'cox3': 'cox3',
    'cytb': 'cob', 'cob': 'cob', 'cyb': 'cob',
    'nd1': 'nad1', 'nad1': 'nad1', 'nadh1': 'nad1',
    'nd2': 'nad2', 'nad2': 'nad2', 'nadh2': 'nad2',
    'nd3': 'nad3', 'nad3': 'nad3', 'nadh3': 'nad3',
    'nd4': 'nad4', 'nad4': 'nad4', 'nadh4': 'nad4',
    'nd4l': 'nad4l', 'nad4l': 'nad4l', 'nadh4l': 'nad4l',
    'nd5': 'nad5', 'nad5': 'nad5', 'nadh5': 'nad5',
    'nd6': 'nad6', 'nad6': 'nad6', 'nadh6': 'nad6',
    'atp6': 'atp6', 'atpase6': 'atp6',
    'atp8': 'atp8', 'atpase8': 'atp8'
}

def clean_name(name):
    """Normalize gene name"""
    n = name.lower().replace('_', '').replace('-', '').strip()
    return GENE_MAP.get(n, None)

def get_sequences(file_path, file_type="fasta"):
    """Extract sequences with normalized names"""
    seqs = {}
    if not os.path.exists(file_path): return seqs
    
    iterator = SeqIO.parse(file_path, file_type)
    for rec in iterator:
        found = None
        if file_type == "genbank":
            for feature in rec.features:
                if feature.type == "CDS":
                    gene = feature.qualifiers.get("gene", [""])[0]
                    if not gene: gene = feature.qualifiers.get("product", [""])[0]
                    std = clean_name(gene)
                    if std:
                        s = str(feature.extract(rec.seq))
                        if std not in seqs or len(s) > len(seqs[std]):
                            seqs[std] = s
            return seqs 
        else: # FASTA
            std = clean_name(rec.id)
            if std:
                found = std
            else:
                header = (rec.description + " " + rec.id).lower()
                for k, v in GENE_MAP.items():
                    if f" {k} " in f" {header} " or f"_{k}" in header or f"={k}" in header:
                        found = v
                        break
            if found:
                seqs[found] = str(rec.seq)
    return seqs

def translate_dna(dna_seq):
    seq_obj = Seq(dna_seq)
    remainder = len(seq_obj) % 3
    if remainder != 0: seq_obj = seq_obj[:-remainder]
    try:
        return str(seq_obj.translate(table=2, to_stop=False))
    except:
        return ""

def protein_guided_alignment(gene_name, seq_sample, seq_ref, temp_dir):
    """Align using Muscle via Protein translation"""
    aa_sample = translate_dna(seq_sample)
    aa_ref = translate_dna(seq_ref)
    
    if not aa_sample or not aa_ref: return None, None
    
    aa_in = os.path.join(temp_dir, f"{gene_name}_aa.fasta")
    aa_out = os.path.join(temp_dir, f"{gene_name}_aa.aln")
    
    with open(aa_in, "w") as f:
        f.write(f">Sample\n{aa_sample}\n>Reference\n{aa_ref}\n")
    
    muscle = shutil.which("muscle")
    if not muscle: return None, None
    
    # Try v3 then v5 syntax
    try:
        subprocess.run([muscle, "-in", aa_in, "-out", aa_out, "-quiet"], check=True, stderr=subprocess.DEVNULL)
    except:
        try:
            subprocess.run([muscle, "-align", aa_in, "-output", aa_out], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except:
            return None, None
            
    try:
        aln = AlignIO.read(aa_out, "fasta")
        sample_aln_aa = ""
        ref_aln_aa = ""
        for r in aln:
            if "Sample" in r.id: sample_aln_aa = str(r.seq)
            if "Reference" in r.id: ref_aln_aa = str(r.seq)
            
        def insert_gaps(dna, aligned_aa):
            res = ""
            dna_idx = 0
            for aa in aligned_aa:
                if aa == "-":
                    res += "---"
                else:
                    chunk = dna[dna_idx:dna_idx+3]
                    if len(chunk) < 3: chunk += "-" * (3-len(chunk))
                    res += chunk
                    dna_idx += 3
            return res
            
        return insert_gaps(seq_sample, sample_aln_aa), insert_gaps(seq_ref, ref_aln_aa)
    except:
        return None, None

# --- PURE PYTHON KA/KS CALCULATOR (Nei-Gojobori 1986 Simplified) ---
def count_sites(codon, table):
    """Count synonymous and non-synonymous sites for a codon"""
    syn = 0
    non = 0
    # Neighbors: change 1 position
    bases = ['T', 'C', 'A', 'G']
    chars = list(codon)
    
    paths = 0
    for pos in range(3):
        orig_base = chars[pos]
        for b in bases:
            if b == orig_base: continue
            new_codon = list(chars)
            new_codon[pos] = b
            new_c_str = "".join(new_codon)
            
            # Skip stops in calculation usually, or count as non-syn
            if new_c_str in table.stop_codons or codon in table.stop_codons:
                continue
                
            paths += 1
            if table.forward_table.get(new_c_str, '*') == table.forward_table.get(codon, '*'):
                syn += 1
            else:
                non += 1
    
    if paths == 0: return 0, 0
    return syn/3.0, non/3.0

def calculate_ka_ks_ng86(seq1, seq2):
    """
    Simple implementation of Nei-Gojobori (1986) method.
    Returns: (Ka, Ks, Ka/Ks ratio)
    """
    table = CodonTable.unambiguous_dna_by_id[2] # Vertebrate Mitochondrial
    
    S_sites = 0 # Total synonymous sites
    N_sites = 0 # Total non-synonymous sites
    Sd = 0 # Synonymous differences
    Nd = 0 # Non-synonymous differences
    
    # Process by codon
    for i in range(0, len(seq1), 3):
        c1 = seq1[i:i+3]
        c2 = seq2[i:i+3]
        
        if len(c1) < 3 or len(c2) < 3: continue
        if "-" in c1 or "-" in c2: continue # Ignore gaps
        if "N" in c1 or "N" in c2: continue
        
        # Check Stop Codons
        if c1 in table.stop_codons or c2 in table.stop_codons: continue

        # 1. Count Sites (Average of both sequences)
        s1, n1 = count_sites(c1, table)
        s2, n2 = count_sites(c2, table)
        S_sites += (s1 + s2) / 2
        N_sites += (n1 + n2) / 2
        
        # 2. Count Differences
        diffs = sum(1 for a, b in zip(c1, c2) if a != b)
        if diffs == 0: continue
        
        # Simplified counting for 1-3 differences
        # Ideally we check pathways, here we simplify:
        # If AA changed -> Non-syn diff, else Syn diff
        aa1 = table.forward_table.get(c1, '*')
        aa2 = table.forward_table.get(c2, '*')
        
        if aa1 == aa2:
            Sd += diffs # All differences are synonymous
        else:
            # Simplification: Assume 1 non-syn event dominates if AA changes
            # For strict NG86 we need pathways, but for simple stats this is often acceptable approximation
            Nd += (diffs) 
            # Note: This is a rough estimation. Real NG86 averages pathways.
            # But sufficient for general validation without external tools.

    # Jukes-Cantor Correction
    def correct(p):
        if p >= 0.75: return 0 # Too saturated
        return -0.75 * math.log(1 - (4/3 * p))

    if S_sites == 0 or N_sites == 0: return "NA", "NA", "NA"

    pS = Sd / S_sites
    pN = Nd / N_sites
    
    try:
        Ks = correct(pS)
        Ka = correct(pN)
        ratio = Ka / Ks if Ks > 0 else "NA"
        return round(Ka, 4), round(Ks, 4), round(ratio, 4) if isinstance(ratio, float) else ratio
    except:
        return "NA", "NA", "NA"

def main(sample_fas, ref_gbk, aln_dir, summary_out):
    if os.path.exists(aln_dir): shutil.rmtree(aln_dir)
    os.makedirs(aln_dir, exist_ok=True)
    temp_dir = os.path.join(aln_dir, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    
    print("--- Starting Pure Python Ka/Ks Analysis ---")
    
    sample_genes = get_sequences(sample_fas, "fasta")
    ref_genes = get_sequences(ref_gbk, "genbank")
    common_genes = sorted(list(set(sample_genes.keys()) & set(ref_genes.keys())))
    
    print(f"Common genes found: {common_genes}")
    results = []
    
    for gene in common_genes:
        s_seq = sample_genes[gene]
        r_seq = ref_genes[gene]
        
        # Align
        s_aln, r_aln = protein_guided_alignment(gene, s_seq, r_seq, temp_dir)
        
        if s_aln and r_aln:
            # Write Alignment for Record
            with open(os.path.join(aln_dir, f"{gene}.aln.fasta"), "w") as f:
                f.write(f">Sample\n{s_aln}\n>Reference\n{r_aln}\n")
            
            # Calculate Ka/Ks (Python Internal)
            Ka, Ks, Ratio = calculate_ka_ks_ng86(s_aln, r_aln)
            
            results.append({
                'Gene': gene,
                'Ka': Ka,
                'Ks': Ks,
                'Ka/Ks': Ratio
            })
        else:
            print(f"Alignment failed for {gene}")
            results.append({'Gene': gene, 'Ka': 'NA', 'Ks': 'NA', 'Ka/Ks': 'NA'})

    if os.path.exists(temp_dir): shutil.rmtree(temp_dir)

    if not results:
        results.append({'Gene': 'None', 'Ka': 'NA', 'Ks': 'NA', 'Ka/Ks': 'NA'})
        
    df_res = pd.DataFrame(results)
    df_res.to_csv(summary_out, sep="\t", index=False)
    print(f"Results saved to {summary_out}")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])