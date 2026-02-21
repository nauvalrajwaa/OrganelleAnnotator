import pandas as pd
import sys
import base64
import os
from datetime import datetime

def get_base64_image(image_path):
    if not os.path.exists(image_path):
        return ""
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')

def read_text_file(file_path):
    if not os.path.exists(file_path):
        return "File not found."
    with open(file_path, "r") as f:
        return f.read()

def generate_html(rscu_plot, rscu_table, codon_stats, kaks_summary, output_html, species_name, tree_file=None, gc_plot=None, aa_plot=None, tree_plot=None, circos_plot=None, synteny_plot=None, synteny_summary=None):
    
    # Process RSCU Data
    rscu_df = pd.read_csv(rscu_table, sep="\t")
    rscu_html_table = rscu_df.head(20).to_html(classes="table table-striped", index=False) # Show top 20
    
    # Process Ka/Ks Data
    try:
        kaks_df = pd.read_csv(kaks_summary, sep="\t")
        if kaks_df.empty:
            kaks_html_table = "<p>No common genes found between sample and reference for Ka/Ks analysis.</p>"
        else:
            kaks_html_table = kaks_df.to_html(classes="table table-striped", index=False)
    except Exception as e:
        kaks_html_table = f"<p>Error loading Ka/Ks data: {e}. File might be empty.</p>"
    
    # Process Codon Stats
    codon_text = read_text_file(codon_stats)
    
    # Process RSCU Plot
    plot_b64 = get_base64_image(rscu_plot)
    
    # Process Tree (optional)
    tree_section = ""
    if tree_plot and os.path.exists(tree_plot):
        tree_b64 = get_base64_image(tree_plot)
        tree_section = f"""
        <div class="section">
            <h2>Phylogeny</h2>
            <p>Maximum Likelihood Tree (IQ-TREE with GTR+G model, 1000 ultrafast bootstrap replicates).</p>
            <img src="data:image/png;base64,{tree_b64}" class="plot-img" alt="Phylogenetic Tree">
            <p><small>Tree file: {tree_file if tree_file else 'N/A'}</small></p>
        </div>
        """
    
    # Process GC Plot
    gc_section = ""
    if gc_plot and os.path.exists(gc_plot):
         gc_b64 = get_base64_image(gc_plot)
         gc_section = f"""
         <div class="col-md-6 mb-3">
             <h4>GC Content</h4>
             <img src="data:image/png;base64,{gc_b64}" class="plot-img" alt="GC Content Plot">
         </div>
         """

    # Process AA Plot
    aa_section = ""
    if aa_plot and os.path.exists(aa_plot):
         aa_b64 = get_base64_image(aa_plot)
         aa_section = f"""
         <div class="col-md-6 mb-3">
             <h4>Amino Acid Composition</h4>
             <img src="data:image/png;base64,{aa_b64}" class="plot-img" alt="AA Composition Plot">
         </div>
         """

    composition_html = ""
    if gc_section or aa_section:
        composition_html = f"""
        <div class="section">
            <h2>Gene Composition Analysis</h2>
            <div class="row">
                {gc_section}
                {aa_section}
            </div>
        </div>
        """
    
    # Process Circos Plot
    circos_html = ""
    if circos_plot and os.path.exists(circos_plot):
        circos_b64 = get_base64_image(circos_plot)
        circos_html = f"""
        <div class="section">
            <h2>Genome Visualization (Circos)</h2>
            <p>Circular genome map showing gene positions, orientations, and GC content.</p>
            <div class="text-center">
                <img src="data:image/png;base64,{circos_b64}" class="plot-img" alt="Circos Genome Plot" style="max-width: 90%;">
            </div>
        </div>
        """
    
    # Process Synteny Analysis
    synteny_html = ""
    if synteny_plot and os.path.exists(synteny_plot):
        synteny_b64 = get_base64_image(synteny_plot)
        synteny_table_html = ""
        
        if synteny_summary and os.path.exists(synteny_summary):
            synteny_df = pd.read_csv(synteny_summary, sep="\t")
            synteny_table_html = synteny_df.to_html(index=False, classes="table table-sm table-striped")
        
        synteny_html = f"""
        <div class="section">
            <h2>Synteny Analysis</h2>
            <p>Genome structure comparison showing conserved blocks, inversions, and rearrangements.</p>
            <div class="text-center mb-3">
                <img src="data:image/png;base64,{synteny_b64}" class="plot-img" alt="Synteny Plot" style="max-width: 100%;">
            </div>
            {synteny_table_html if synteny_table_html else ""}
        </div>
        """

    # HTML Template
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Mitochondrial Genome Analysis Report</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
        <style>
            body {{ padding: 20px; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }}
            .section {{ margin-bottom: 40px; padding: 20px; background-color: #f8f9fa; border-radius: 8px; }}
            h1, h2 {{ color: #2c3e50; }}
            .plot-img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; padding: 5px; }}
            pre {{ background: #eaeaea; padding: 10px; border-radius: 5px; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1 class="text-center mb-4">Mitochondrial Analysis Report</h1>
            <p class="text-center text-muted">Generated on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} for species: <strong>{species_name}</strong></p>
            
            {composition_html}

            <div class="section">
                <h2>Relative Synonymous Codon Usage (RSCU)</h2>
                <div class="row">
                    <div class="col-md-12 mb-3">
                        <h4>Visualization</h4>
                        <img src="data:image/png;base64,{plot_b64}" class="plot-img" alt="RSCU Plot">
                    </div>
                    <div class="col-md-12">
                        <h4>Data Table (Top 20 rows)</h4>
                        <div class="table-responsive">
                            {rscu_html_table}
                        </div>
                        <p><em>Full data available in <code>results/rscu/rscu.tsv</code></em></p>
                    </div>
                </div>
            </div>
            
            <div class="section">
                <h2>Start/Stop Codon Analysis</h2>
                <div class="row">
                    <div class="col-md-12">
                        <pre>{codon_text}</pre>
                    </div>
                </div>
            </div>
            
            <div class="section">
                <h2>Pairwise Analysis (Ka/Ks Proxy)</h2>
                <p>Simple pairwise comparison against reference.</p>
                <div class="table-responsive">
                    {kaks_html_table}
                </div>
            </div>

            {tree_section}
            
            {circos_html}
            
            {synteny_html}
            
        </div>
    </body>
    </html>
    """
    
    with open(output_html, "w") as f:
        f.write(html_content)

if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("Usage: python generate_html_report.py rscu_plot rscu_table codon_stats kaks_summary species_name output_html tree_file gc_plot aa_plot")
        sys.exit(1)
    
    # Args are fixed position based on Snakefile
    # 1: rscu_plot
    # 2: rscu_table
    # 3: codon_stats
    # 4: kaks_summary
    # 5: species_name
    # 6: output_html
    # 7: tree_file
    # 8: gc_plot
    # 9: aa_plot
    
    tree = sys.argv[7] if len(sys.argv) > 7 else None
    gc_plot = sys.argv[8] if len(sys.argv) > 8 else None
    aa_plot = sys.argv[9] if len(sys.argv) > 9 else None
    tree_plot = sys.argv[10] if len(sys.argv) > 10 else None
    circos_plot = sys.argv[11] if len(sys.argv) > 11 else None
    synteny_plot = sys.argv[12] if len(sys.argv) > 12 else None
    synteny_summary = sys.argv[13] if len(sys.argv) > 13 else None
        
    generate_html(
        sys.argv[1], 
        sys.argv[2], 
        sys.argv[3], 
        sys.argv[4], 
        sys.argv[6], 
        sys.argv[5],
        tree,
        gc_plot,
        aa_plot,
        tree_plot,
        circos_plot,
        synteny_plot,
        synteny_summary
    )
