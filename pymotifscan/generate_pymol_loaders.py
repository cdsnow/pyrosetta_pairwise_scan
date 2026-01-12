#!/usr/bin/env python3
"""
Generate PyMOL loader scripts (.pml) for candidates in candidates_v4.
Creates one .pml file per signature group (e.g. F1_ASP_BB_N).
"""

import os
import glob
import argparse

def main():
    parser = argparse.ArgumentParser(description="Generate PyMOL loader scripts for pruned candidates")
    parser.add_argument('--indir', '-i', default='candidates_pruned', help='Input directory (e.g. candidates_pruned)')
    parser.add_argument('--out-prefix', default='load_v4_', help='Prefix for generated .pml files')
    args = parser.parse_args()

    if not os.path.exists(args.indir):
        print(f"Error: Directory {args.indir} not found.")
        return

    # Find all first-level subdirectories (the groups like F1_ASP_BB_N)
    groups = sorted([d for d in os.listdir(args.indir) if os.path.isdir(os.path.join(args.indir, d))])
    
    print(f"Found {len(groups)} groups in {args.indir}.")

    for group in groups:
        group_path = os.path.join(args.indir, group)
        pdbs = sorted(glob.glob(os.path.join(group_path, "*.pdb")))
        
        if not pdbs:
            continue
            
        pml_filename = f"{args.out_prefix}{group}.pml"
        print(f"  Generating {pml_filename} ({len(pdbs)} signatures)...")
        
        with open(pml_filename, 'w') as f:
            f.write(f"# PyMOL Loader for {group} (Pruned v4)\n")
            # Using relative path for the cd command
            f.write(f"cd {group_path}\n\n")
            
            for pdb in pdbs:
                basename = os.path.basename(pdb)
                sig_name = basename.replace(".pdb", "")
                # Sanitize sig_name for PyMOL (replace dots or other special chars if any)
                # But based on original file, it uses dots directly in object names? 
                # Actually original used sig.ARG_NE
                
                f.write(f"delete sig.{sig_name}\n")
                f.write(f"load {basename}, sig.{sig_name}\n")
        
    print("\nDone.")

if __name__ == "__main__":
    main()
