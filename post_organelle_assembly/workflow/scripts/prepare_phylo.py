#!/usr/bin/env python3
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess
import os
import sys
import shutil
import glob

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
    """Normalize gene name using the map"""
    n = name.lower().replace('_', '').replace('-', '').strip()
    if n in GENE_MAP: return GENE_MAP[n]
    return GENE_MAP.get(n, None)

def extract_cds_from_gbk(gbk_input):
    """Returns a dict of dicts: {accession: {gene: seq_string}}"""
    species_cds = {}
    files = []
    
    # Handle single file or directory input
    if os.path.isfile(gbk_input):
        files = [gbk_input]
        ref_dir = os.path.dirname(gbk_input)
        if not ref_dir: ref_dir = "resources/references"
        if os.path.exists(ref_dir):
            files += glob.glob(os.path.join(ref_dir, "*.gbk"))
            files = list(set(files)) 
    elif os.path.isdir(gbk_input):
        files = glob.glob(os.path.join(gbk_input, "*.gbk"))

    print(f"Loading references from {len(files)} files...")

    for fpath in files:
        for record in SeqIO.parse(fpath, "genbank"):
            acc = record.id
            organism = record.annotations.get("organism", acc).replace(" ", "_")
            taxon_id = organism
            
            genes = {}
            for feature in record.features:
                if feature.type == "CDS":
                    gene_raw = feature.qualifiers.get("gene", [""])[0]
                    if not gene_raw:
                        gene_raw = feature.qualifiers.get("product", [""])[0]
                    
                    std_name = clean_name(gene_raw)
                    if not std_name: continue
                    
                    try:
                        seq = feature.extract(record.seq)
                        if std_name not in genes or len(seq) > len(genes[std_name]):
                            genes[std_name] = str(seq)
                    except: pass
            
            if genes:
                species_cds[taxon_id] = genes
        
    return species_cds

def run_muscle_robust(input_file, output_file):
    """Runs Muscle handling both v3 (-in -out) and v5 (-align -output)"""
    muscle_exe = shutil.which("muscle")
    if not muscle_exe: return False

    # v3 syntax
    cmd_v3 = [muscle_exe, "-in", input_file, "-out", output_file, "-quiet"]
    # v5 syntax
    cmd_v5 = [muscle_exe, "-align", input_file, "-output", output_file]

    try:
        subprocess.run(cmd_v3, check=True, stderr=subprocess.PIPE, timeout=120)
        return True
    except:
        try:
            subprocess.run(cmd_v5, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
            return True
        except:
            return False

def main(sample_fasta, ref_gbk, output_aln, output_partition):
    temp_dir = os.path.join(os.path.dirname(output_aln), "tmp_phylo")
    os.makedirs(temp_dir, exist_ok=True)

    # 1. Load Sample
    print(f"Parsing Sample: {sample_fasta}")
    sample_genes = {}
    for rec in SeqIO.parse(sample_fasta, "fasta"):
        header = (rec.description + " " + rec.id).lower()
        found_gene = None
        for key, val in GENE_MAP.items():
            if f" {key} " in f" {header} " or f"_{key}" in header or f"={key}" in header:
                found_gene = val
                break
        if found_gene:
            sample_genes[found_gene] = str(rec.seq)

    # 2. Load References
    refs = extract_cds_from_gbk(ref_gbk)
    
    # 3. Process Genes
    candidate_genes = sorted(list(set(sample_genes.keys())))
    print(f"Genes found in sample: {candidate_genes}")
    
    all_taxa = ["Sample"] + sorted(list(refs.keys()))
    concatenated_seqs = {taxa: "" for taxa in all_taxa}
    partitions = []
    current_len = 0
    valid_genes_count = 0
    
    for gene in candidate_genes:
        gene_records = []
        
        # Add Sample
        gene_records.append(SeqRecord(Seq(sample_genes[gene]), id="Sample", description=""))
        
        # Add Refs
        refs_with_gene = 0
        for ref_id, ref_genes in refs.items():
            if gene in ref_genes:
                gene_records.append(SeqRecord(Seq(ref_genes[gene]), id=ref_id, description=""))
                refs_with_gene += 1
        
        if refs_with_gene == 0:
            print(f"Skipping {gene}: Not found in references.")
            continue
            
        # Write Temp & Align
        tmp_in = os.path.join(temp_dir, f"{gene}.in.fas")
        tmp_out = os.path.join(temp_dir, f"{gene}.aln.fas")
        SeqIO.write(gene_records, tmp_in, "fasta")
        
        print(f"Aligning {gene} ({len(gene_records)} taxa)...")
        if run_muscle_robust(tmp_in, tmp_out) and os.path.exists(tmp_out):
            try:
                aln = AlignIO.read(tmp_out, "fasta")
                aln_len = aln.get_alignment_length()
                seq_map = {rec.id: str(rec.seq) for rec in aln}
                
                for taxa in all_taxa:
                    if taxa in seq_map:
                        concatenated_seqs[taxa] += seq_map[taxa]
                    else:
                        concatenated_seqs[taxa] += "-" * aln_len
                
                # --- FIX IS HERE ---
                # Changed from "DNA, gene = ..." to "charset gene = ..."
                partitions.append(f"charset {gene} = {current_len + 1}-{current_len + aln_len}")
                current_len += aln_len
                valid_genes_count += 1
            except Exception as e:
                print(f"Error reading alignment for {gene}: {e}")
        else:
            print(f"Alignment failed for {gene}")

    # 4. Write Output
    if valid_genes_count == 0:
        print("ERROR: No genes aligned. Creating dummy output.")
        with open(output_aln, 'w') as f: f.write(">Dummy\nN")
        with open(output_partition, 'w') as f: f.write("#NEXUS")
    else:
        # Write FASTA Supermatrix
        final_recs = []
        for taxa, seq in concatenated_seqs.items():
            if len(seq) > 0 and set(seq) != {'-'}:
                final_recs.append(SeqRecord(Seq(seq), id=taxa, description=""))
        SeqIO.write(final_recs, output_aln, "fasta")
        
        # Write NEXUS Partition (Corrected)
        with open(output_partition, "w") as f:
            f.write("#NEXUS\n")
            f.write("BEGIN SETS;\n")
            for p in partitions:
                f.write(f"\t{p};\n")
            f.write("END;\n")

    print(f"Done. Supermatrix: {current_len} bp, Taxa: {len(final_recs)}")
    if os.path.exists(temp_dir): shutil.rmtree(temp_dir)

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python prepare_phylo.py <sample> <ref> <out_aln> <out_part>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])