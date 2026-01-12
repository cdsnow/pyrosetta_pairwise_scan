#!/usr/bin/env python3
"""
Generate a pruned_results.json by filtering merged_results.json to only include
models that were kept in candidates_v3.
"""

import json
import os
import argparse
import glob

def parse_pdb_models_kept(pdb_file):
    """Parse model numbers from a PDB file."""
    model_nums = set()
    if not os.path.exists(pdb_file):
        return model_nums
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("MODEL"):
                try:
                    num = int(line.split()[1])
                    model_nums.add(num)
                except (ValueError, IndexError):
                    continue
    # If no MODEL lines found, assume it was a single model (model 1)
    if not model_nums:
        model_nums.add(1)
    return model_nums

def main():
    parser = argparse.ArgumentParser(description="Prune results JSON based on kept candidates")
    parser.add_argument('--results', '-r', default='merged_results.json', help='Input merged results JSON')
    parser.add_argument('--input-dir', '-i', default='candidates_transformed', help='Directory used for the scan (transformed)')
    parser.add_argument('--pruned-dir', '-p', default='candidates_pruned', help='Directory with pruned PDBs')
    parser.add_argument('--output', '-o', default='pruned_results.json', help='Output results JSON')
    args = parser.parse_args()

    print(f"Loading {args.results}...")
    with open(args.results) as f:
        data = json.load(f)

    results = data.get('results', {})
    pruned_results = {}
    
    total_in = 0
    total_out = 0
    
    print("Filtering results...")
    for pdb_path, entry in results.items():
        # pdb_path in results is likely relative to the execution dir, e.g. "candidates_v2/..."
        # We need to find the corresponding file in candidates_v3
        rel_path = os.path.relpath(pdb_path, args.input_dir)
        v3_path = os.path.join(args.pruned_dir, rel_path)
        
        if not os.path.exists(v3_path):
            # Signature was entirely pruned
            continue
            
        kept_models = parse_pdb_models_kept(v3_path)
        
        models = entry.get('models', [])
        total_in += len(models)
        
        # Filter models list
        new_models = [m for m in models if m.get('model_num') in kept_models]
        
        if new_models:
            new_entry = entry.copy()
            new_entry['models'] = new_models
            # Update pdb_file path to the pruned one if desired, or keep original?
            # Keeping original reference but could update to v3_path.
            new_entry['pdb_file'] = v3_path 
            pruned_results[v3_path] = new_entry
            total_out += len(new_models)

    # Update metadata
    new_data = data.copy()
    new_data['results'] = pruned_results
    if 'metadata' in new_data:
        new_data['metadata']['candidates_dir'] = args.pruned_dir
        new_data['metadata']['total_processed'] = total_out
        new_data['metadata']['note'] = "Pruned via structural redundancy and fa_rep > 10"

    print(f"Writing {args.output}...")
    with open(args.output, 'w') as f:
        json.dump(new_data, f, indent=2)

    print("\nSummary:")
    print(f"  Input signatures:  {len(results)}")
    print(f"  Output signatures: {len(pruned_results)}")
    print(f"  Input models:      {total_in}")
    print(f"  Output models:     {total_out}")

if __name__ == "__main__":
    main()
