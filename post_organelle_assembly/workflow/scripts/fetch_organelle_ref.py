import time
import sys
import os
import argparse
from Bio import Entrez
from urllib.error import HTTPError, URLError

def robust_entrez_search(term, retmax):
    """Wrapper search yang sabar menghadapi server NCBI."""
    attempt = 0
    max_attempts = 5
    
    while attempt < max_attempts:
        try:
            handle = Entrez.esearch(db="nucleotide", term=term, retmax=retmax, idtype="acc")
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"]
        
        except (HTTPError, URLError, RuntimeError) as e:
            print(f"   [!] Connection warning (attempt {attempt+1}/{max_attempts}): {e}")
            print("       -> Sleeping for 20 seconds...")
            time.sleep(20)
            attempt += 1
            
    print(f"   [X] Failed to search term after {max_attempts} attempts.")
    return []

def get_parent_taxa(organism):
    """Mencari Genus dan Family dari spesies target."""
    try:
        # 1. Cari TaxID
        handle = Entrez.esearch(db="taxonomy", term=organism)
        record = Entrez.read(handle)
        handle.close()
        
        if not record['IdList']: return []
        tax_id = record['IdList'][0]
        
        # 2. Fetch Lineage
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        if "Lineage" in records[0]:
            # Lineage biasanya string: "cellular organisms; Eukaryota; Viridiplantae; ..."
            # Kita split jadi list
            full_lineage = records[0]["Lineage"].split("; ")
            # Kita return urutan terbalik (dari Genus ke atas)
            return full_lineage[::-1]
            
    except Exception as e:
        print(f"   [!] Taxonomy lookup warning: {e}")
        return []
    return []

def build_search_list(target_species, expand_lineage):
    """Membuat daftar organisme yang akan dicari berdasarkan config."""
    search_targets = []
    
    # 1. Target Utama (Selalu dicari)
    search_targets.append((target_species, 20)) 
    
    # 2. Jika Mode Expanded aktif
    if expand_lineage:
        print(f"   [INFO] Mode Expanded aktif. Mencari kerabat dari: {target_species}...")
        parents = get_parent_taxa(target_species)
        
        taxa_to_add = parents[:2] 
        
        for taxa in taxa_to_add:
            search_targets.append((taxa, 5))
            print(f"   [INFO] Menambahkan kerabat ke pencarian: {taxa}")
            
    return search_targets

def search_and_fetch(args):
    target_species = args.species
    organelle_query = args.organelle_type  # Configurable organelle type
    
    mode_str = "Expanded (Kerabat)" if args.expand else "Specific (Target)"
    print(f"--- [NCBI FETCH] Mode: {mode_str} | Target: {target_species} | Organelle: {organelle_query} ---")
    
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
        
    unique_ids = set()
    search_list = build_search_list(target_species, args.expand)
    
    # MULAI PENCARIAN
    for organism, retmax in search_list:
        print(f"   > Searching for: {organism} (Limit: {retmax})...")
        
        # Query yang presisi
        term = (f'"{organism}"[Organism] AND ({organelle_query}) '
                f'AND (complete genome[Title] OR complete sequence[Title]) '
                f'AND {args.min_len}:{args.max_len}[Sequence Length]')
        
        ids = robust_entrez_search(term, retmax)
        
        if ids:
            new_ids = [x for x in ids if x not in unique_ids]
            if new_ids:
                print(f"     [+] Found {len(new_ids)} new sequence(s).")
                unique_ids.update(new_ids)
            else:
                print(f"     [.] Found IDs but already in list.")
        else:
            print(f"     [-] No match found for {organism}.")
        
        time.sleep(2) 

    # FALLBACK
    if not unique_ids:
        print("\n[WARNING] Tidak ada ID ditemukan sama sekali. Menggunakan FALLBACK NCBI Standard.")
        fallback_id = "NC_006084.1" # Default fallback
        print(f"   -> Fallback ID: {fallback_id}")
        unique_ids.add(fallback_id)

    # DOWNLOADING
    num_ids = len(unique_ids)
    print(f"\n--- Downloading {num_ids} sequences ---")
    
    if num_ids == 0:
        print("[FATAL] Tidak ada ID untuk didownload.")
        sys.exit(1)

    id_list = list(unique_ids)
    
    try:
        # Download FASTA
        print("   > Downloading FASTA...")
        net_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        fasta_data = net_handle.read()
        net_handle.close()
        
        with open(args.output_fasta, "w") as out_f:
            out_f.write(fasta_data)
            
        # Download GBK
        print("   > Downloading GBK...")
        net_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
        gbk_data = net_handle.read()
        net_handle.close()
        
        with open(args.output_gbk, "w") as out_g:
            out_g.write(gbk_data)
            
        print(f"   > Success! Output saved to: {args.output_fasta}")
        
    except Exception as e:
        print(f"[FATAL ERROR] Download failed: {e}")
        if os.path.exists(args.output_fasta): os.remove(args.output_fasta)
        if os.path.exists(args.output_gbk): os.remove(args.output_gbk)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch Organelle Reference Genomes")
    parser.add_argument("--species", required=True, help="Target species name")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file path")
    parser.add_argument("--output_gbk", required=True, help="Output GenBank file path")
    parser.add_argument("--organelle_type", default="mitochondrion", help="Organelle type: mitochondrion, plastid, or chloroplast")
    parser.add_argument("--min_len", type=int, default=10000, help="Minimum sequence length")
    parser.add_argument("--max_len", type=int, default=200000, help="Maximum sequence length")
    parser.add_argument("--expand", action="store_true", help="Expand search to lineage")
    parser.add_argument("--api_key", help="NCBI API Key")
    
    args = parser.parse_args()
    
    # If the user asks for a tree, we should probably force expand if it's not set?
    # But let's leave it to Snakemake params.
    
    search_and_fetch(args)