#!/usr/bin/env python3
import os
import json
import numpy as np
from scipy.spatial.distance import cdist
import argparse

def parse_pdb_coords(lines):
    """Extracts coordinates from PDB lines."""
    coords = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except ValueError:
                continue
    return np.array(coords)

def load_target_coords(target_path):
    """Loads target coordinates from PDB file."""
    with open(target_path, 'r') as f:
        lines = f.readlines()
    return parse_pdb_coords(lines)

def get_model_coords(pdb_path, model_num):
    """Extracts coordinates for a specific model from a multi-model PDB."""
    coords = []
    current_model_num = None
    in_model = False
    model_lines = []
    
    # Heuristic: if file has no MODEL records, assume it is model 1
    has_models = False
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('MODEL'):
                has_models = True
                current_model_num = int(line.split()[1])
                if current_model_num == model_num:
                    in_model = True
            elif line.startswith('ENDMDL'):
                if in_model:
                    return parse_pdb_coords(model_lines)
                in_model = False
            elif in_model:
                model_lines.append(line)
            elif not has_models and model_num == 1 and (line.startswith('ATOM') or line.startswith('HETATM')):
                # Case for single-model file without MODEL/ENDMDL tags
                model_lines.append(line)
                
    if not has_models and model_num == 1:
        return parse_pdb_coords(model_lines)
        
    return None

def main():
    parser = argparse.ArgumentParser(description="Validate Rosetta energy scan against steric clashes.")
    parser.add_argument('--results', '-r', default='energy_scan_results_v2.json', help='Path to results JSON')
    parser.add_argument('--target', '-t', default='flagorigin_and_loop.pdb', help='Path to target PDB')
    parser.add_argument('--output', '-o', default='steric_validation.csv', help='Output CSV file')
    args = parser.parse_args()

    print(f"Loading target: {args.target}")
    target_coords = load_target_coords(args.target)
    print(f"Target has {len(target_coords)} atoms.")

    print(f"Loading results: {args.results}")
    with open(args.results, 'r') as f:
        data = json.load(f)
    
    results = data.get('results', {})
    print(f"Found {len(results)} candidate files in results.")

    print("Validating sterics...")
    print(f"{'PDB':<60} | {'Model':<5} | {'Min Dist':<10} | {'Delta fa_rep':<12}")
    print("-" * 100)

    # Prepare for stats
    validation_data = []

    # Iterate through results
    for pdb_path, entry in results.items():
        if not os.path.exists(pdb_path):
            continue
            
        for model in entry.get('models', []):
            if 'error' in model or 'delta_fa_rep' not in model:
                continue
            
            model_num = model['model_num']
            delta_rep = model['delta_fa_rep']
            
            # Get coords for this model
            motif_coords = get_model_coords(pdb_path, model_num)
            
            if motif_coords is None or len(motif_coords) == 0:
                continue
            
            # Compute all-vs-all distances
            dists = cdist(target_coords, motif_coords)
            min_dist = np.min(dists)
            
            validation_data.append({
                'pdb': pdb_path,
                'model': model_num,
                'min_dist': min_dist,
                'delta_fa_rep': delta_rep
            })

    # Sort by min_dist to find clashes
    validation_data.sort(key=lambda x: x['min_dist'])

    # Write to CSV
    with open(args.output, 'w') as f:
        f.write("pdb_file,model_num,min_distance,delta_fa_rep\n")
        for item in validation_data:
            f.write(f"{item['pdb']},{item['model']},{item['min_dist']:.4f},{item['delta_fa_rep']:.4f}\n")

    # Print summary of closest contacts (potential clashes)
    print("\nTop 20 Closest Contacts (Potential Clashes):")
    print(f"{'PDB Filename (truncated)':<40} | {'Mod':<3} | {'Dist':<6} | {'fa_rep':<10}")
    for item in validation_data[:20]:
        short_name = item['pdb'].replace('candidates_v2/', '')
        if len(short_name) > 40: short_name = "..." + short_name[-37:]
        print(f"{short_name:<40} | {item['model']:<3} | {item['min_dist']:.2f} | {item['delta_fa_rep']:.2f}")

    # Analyze correlation for low distance vs low repulsion (The dangerous quadrant)
    # Define thresholds
    clash_threshold = 2.0  # Angstroms (severe clash usually < 2.0 or 2.5)
    rep_threshold = 10.0   # Rosetta units
    
    hidden_clashes = [x for x in validation_data if x['min_dist'] < clash_threshold and x['delta_fa_rep'] < rep_threshold]
    
    print(f"\nAnalysis:")
    print(f"Total models analyzed: {len(validation_data)}")
    print(f"Models with min_dist < {clash_threshold}A: {len([x for x in validation_data if x['min_dist'] < clash_threshold])}")
    print(f"Models with low repulsion (< {rep_threshold}) but min_dist < {clash_threshold}A (Hidden Clashes): {len(hidden_clashes)}")
    
    if hidden_clashes:
        print("\nTop 10 Hidden Clashes (Low Repulsion, Short Distance):")
        hidden_clashes.sort(key=lambda x: x['min_dist'])
        for item in hidden_clashes[:10]:
            short_name = item['pdb'].replace('candidates_v2/', '')
            if len(short_name) > 40: short_name = "..." + short_name[-37:]
            print(f"{short_name:<40} | {item['model']:<3} | {item['min_dist']:.2f} | {item['delta_fa_rep']:.2f}")

if __name__ == "__main__":
    main()
