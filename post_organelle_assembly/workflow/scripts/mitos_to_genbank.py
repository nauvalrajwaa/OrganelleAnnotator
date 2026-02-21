#!/usr/bin/env python3
"""
Convert GFF output (from MITOS/MITFI) to GenBank format
Reads GFF file and creates proper GenBank file
"""
import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

def parse_gff_attributes(attr_str):
    """Parse GFF attribute string into dictionary"""
    attributes = {}
    for attr in attr_str.strip().split(';'):
        if not attr.strip():
            continue
        if '=' in attr:
            key, val = attr.strip().split('=', 1)
            attributes[key.strip()] = val.strip()
    return attributes

def parse_mitos_to_genbank(gff_file, fasta_file, output_gbk):
    """Convert MITOS GFF output to GenBank format"""
    
    # Read sequence
    # Handle single or multiple records (take first if multiple is common)
    try:
        seq_record = SeqIO.read(fasta_file, "fasta")
    except ValueError:
        with open(fasta_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_record = record
                break
    
    # Parse GFF file
    features = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            source = parts[1]
            feature_type = parts[2]
            start = int(parts[3]) - 1  # Convert to 0-based
            end = int(parts[4])
            strand_val = parts[6]
            strand = 1 if strand_val == '+' else -1
            attributes = parse_gff_attributes(parts[8])
            
            gene_name = attributes.get('Name', '')
            
            # Map feature types
            feat_type = None
            
            # MITOS uses 'gene' for CDS
            if feature_type == 'gene' and source == 'mitos':
                feat_type = 'CDS'
            # MITFI uses 'tRNA' and 'rRNA'
            elif feature_type == 'tRNA':
                feat_type = 'tRNA'
            elif feature_type == 'rRNA':
                feat_type = 'rRNA'
            
            if feat_type:
                # Create feature
                location = FeatureLocation(start, end, strand=strand)
                feature = SeqFeature(location, type=feat_type)
                
                # Add qualifiers
                if gene_name:
                    feature.qualifiers['gene'] = [gene_name]
                    feature.qualifiers['product'] = [gene_name]
                
                # Add gene_id if available
                if 'gene_id' in attributes:
                    feature.qualifiers['gene_id'] = [attributes['gene_id']]
                
                features.append(feature)
        
        # Sort features by start position
        features.sort(key=lambda x: x.location.start)
    
    # Add features to record
    seq_record.features = features
    seq_record.annotations['molecule_type'] = 'DNA'
    seq_record.annotations['topology'] = 'circular'
    
    # Write GenBank
    SeqIO.write(seq_record, output_gbk, "genbank")
    print(f"GenBank file created: {output_gbk}")
    print(f"Total features: {len(features)}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python mitos_to_genbank.py <gff_file> <fasta_file> <output_gbk>")
        sys.exit(1)
    
    parse_mitos_to_genbank(sys.argv[1], sys.argv[2], sys.argv[3])
