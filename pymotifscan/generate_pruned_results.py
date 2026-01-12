#!/usr/bin/env python3
"""
Generate a pruned_results.json by filtering merged_results.json to only include
models that were kept in candidates_pruned.
Also enriches the JSON with metadata (DOF, waters) from PDB remarks and optimizes structure.
"""

import json
import os
import argparse
import glob

def parse_pdb_metadata(pdb_file):
    """
    Parse model metadata from a PDB file.
    Returns a dict: { model_num: {'spanning_dof': int, 'bridging_waters': int} }
    """
    metadata = {}
    if not os.path.exists(pdb_file):
        return metadata
    
    current_model = None
    current_data = {}
    
    # Defaults if remarks are missing
    default_data = {'spanning_dof': None, 'bridging_waters': None}

    with open(pdb_file) as f:
        for line in f:
            if line.startswith("MODEL"):
                try:
                    current_model = int(line.split()[1])
                    current_data = default_data.copy()
                except (ValueError, IndexError):
                    current_model = None
            elif line.startswith("REMARK 501 Spanning DOF:"):
                if current_model is not None:
                    try:
                        current_data['spanning_dof'] = int(line.split(':')[1].strip())
                    except ValueError:
                        pass
            elif line.startswith("REMARK 503 Bridging Waters:"):
                if current_model is not None:
                    try:
                        current_data['bridging_waters'] = int(line.split(':')[1].strip())
                    except ValueError:
                        pass
            elif line.startswith("ENDMDL"):
                if current_model is not None:
                    metadata[current_model] = current_data
                    current_model = None

    # Handle single-model files without MODEL/ENDMDL tags explicitly if needed,
    # but the pipeline produces multi-model files.
    # If the file had no MODEL tags but had data, we might miss it with this logic.
    # However, the input files are known to be multi-model or at least have MODEL tags.
    if not metadata and current_model is not None: 
         # Case where ENDMDL might be missing at EOF
         metadata[current_model] = current_data

    return metadata

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
    
    # Global target data (extracted from the first valid model encountered)
    target_data = None

    print("Filtering and enriching results...")
    for pdb_path, entry in results.items():
        # pdb_path in results is likely relative to the execution dir, e.g. "candidates_v2/..."
        # We need to find the corresponding file in candidates_pruned
        rel_path = os.path.relpath(pdb_path, args.input_dir)
        v3_path = os.path.join(args.pruned_dir, rel_path)
        
        if not os.path.exists(v3_path):
            # Signature was entirely pruned
            continue
            
        # Parse metadata for all models in this PDB
        pdb_metadata = parse_pdb_metadata(v3_path)
        kept_models = set(pdb_metadata.keys())
        
        models = entry.get('models', [])
        total_in += len(models)
        
        new_models = []
        for m in models:
            m_num = m.get('model_num')
            if m_num in kept_models:
                # Extract target data once
                if target_data is None and 'target' in m:
                    target_data = m['target']
                
                # Create optimized model entry
                new_m = {
                    'model_num': m_num,
                    'original_score': m.get('original_score'),
                    'original_id': m.get('original_id'),
                    # Keep motif scores
                    'motif': m.get('motif'),
                    # Keep deltas
                    'delta_total': m.get('delta_total'),
                    'delta_fa_rep': m.get('delta_fa_rep'),
                    'delta_hbond': m.get('delta_hbond'),
                    # Add new metadata
                    'spanning_dof': pdb_metadata[m_num].get('spanning_dof'),
                    'bridging_waters': pdb_metadata[m_num].get('bridging_waters'),
                    # Keep other useful fields
                    'chains_found': m.get('chains_found'),
                    'n_residues': m.get('n_residues')
                }
                new_models.append(new_m)
        
        if new_models:
            new_entry = {
                'pdb_file': v3_path,
                'models': new_models
            }
            pruned_results[v3_path] = new_entry
            total_out += len(new_models)

    # Update metadata
    new_data = {}
    if 'metadata' in data:
        new_data['metadata'] = data['metadata'].copy()
        new_data['metadata']['candidates_dir'] = args.pruned_dir
        new_data['metadata']['total_processed'] = total_out
        new_data['metadata']['note'] = "Pruned via structural redundancy and fa_rep > 10"
        
        # Add the target info to top-level metadata
        if target_data:
            new_data['metadata']['target_scores'] = target_data

    new_data['results'] = pruned_results

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