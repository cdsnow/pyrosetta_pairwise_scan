#!/usr/bin/env python3
import os
import json
import numpy as np
from scipy.spatial.distance import cdist
import argparse

def parse_pdb_coords(lines, heavy_only=False):
    """Extracts coordinates from PDB lines."""
    coords = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if heavy_only:
                element = line[76:78].strip().upper()
                name = line[12:16].strip().upper()
                
                # Known heavy elements in biology
                known_heavy = {'C', 'N', 'O', 'S', 'P', 'F', 'CL', 'BR', 'I', 'FE', 'MG', 'CA', 'ZN', 'MN', 'NA', 'K', 'LI', 'AG', 'AU', 'CU', 'CO', 'NI'}
                
                is_hydrogen = False
                
                if element == 'H' or element == 'D':
                    is_hydrogen = True
                elif element in known_heavy:
                    is_hydrogen = False
                else:
                    # Element is missing or trash (e.g. '1'). Fallback to name.
                    if name.startswith('H'):
                         # Careful with Hg (Mercury). But assuming protein context without Mercury for now.
                         # If we really need to support Hg, we check resname or if name is specifically HG and no HG1/HG2.
                         if name != 'HG': 
                             is_hydrogen = True
                    elif len(name) > 1 and name[0].isdigit() and name[1] == 'H':
                        # e.g. 1HD2, 2HD2
                        is_hydrogen = True
                
                if is_hydrogen:
                    continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except ValueError:
                continue
    return np.array(coords)

def load_target_coords(target_path, heavy_only=False):
    """Loads target coordinates from PDB file."""
    with open(target_path, 'r') as f:
        lines = f.readlines()
    return parse_pdb_coords(lines, heavy_only)

def get_model_coords(pdb_path, model_num, heavy_only=False):
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
                    return parse_pdb_coords(model_lines, heavy_only)
                in_model = False
            elif in_model:
                model_lines.append(line)
            elif not has_models and model_num == 1 and (line.startswith('ATOM') or line.startswith('HETATM')):
                # Case for single-model file without MODEL/ENDMDL tags
                model_lines.append(line)
                
    if not has_models and model_num == 1:
        return parse_pdb_coords(model_lines, heavy_only)
        
    return None

def main():
    parser = argparse.ArgumentParser(description="Validate Rosetta energy scan against steric clashes.")
    parser.add_argument('--results', '-r', default='energy_scan_results_v2.json', help='Path to results JSON')
    parser.add_argument('--target', '-t', default='flagorigin_and_loop.pdb', help='Path to target PDB')
    parser.add_argument('--output', '-o', default='steric_validation.csv', help='Output CSV file')
    parser.add_argument('--heavy-only', action='store_true', help='Only consider heavy atoms for clash detection')
    args = parser.parse_args()

    print(f"Loading target: {args.target}")
    target_coords = load_target_coords(args.target, args.heavy_only)
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
            motif_coords = get_model_coords(pdb_path, model_num, args.heavy_only)
            
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
        
        # We need to re-load coords to identify atoms for these top 10
        # This is inefficient but fine for just 10 items.
        
        for item in hidden_clashes[:20]: # Show top 20
            short_name = item['pdb'].replace('candidates_v2/', '')
            if len(short_name) > 40: short_name = "..." + short_name[-37:]
            
            # Find the atoms
            motif_coords = get_model_coords(item['pdb'], item['model'], args.heavy_only)
            if motif_coords is None: continue
            
            dists = cdist(target_coords, motif_coords)
            t_idx, m_idx = np.unravel_index(np.argmin(dists), dists.shape)
            
            # We need to get the actual atom lines to know the types. 
            # Re-reading file is painful. 
            # For now, just print the stats.
            print(f"{short_name:<40} | {item['model']:<3} | {item['min_dist']:.2f} | {item['delta_fa_rep']:.2f}")

if __name__ == "__main__":
    main()
