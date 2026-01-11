#!/usr/bin/env python3
"""
Prune candidate motifs to create candidates_v3.

Pruning criteria:
1. Remove any candidate with delta_fa_rep > 10 (steric clashes)
2. Remove structurally redundant candidates within the same signature:
   - If two motifs are highly similar (avg min distance < threshold)
   - AND one is worse in BOTH delta_total AND delta_hbond
   - Remove the worse one

Uses spatial similarity metric from measure_motif_redundancy.py
"""

import os
import sys
import json
import argparse
import numpy as np
from collections import defaultdict
import shutil


def is_heavy(name):
    """Check if atom name is a heavy atom (not hydrogen)."""
    n = name.strip()
    if n.startswith('H'):
        return False
    if n[0].isdigit() and len(n) > 1 and n[1] == 'H':
        return False
    return True


def get_element(name):
    """Heuristic element guess from atom name."""
    n = name.strip()
    if n[0].isdigit():
        n = n[1:]
    return n[0]  # C, N, O, S, etc.


def parse_model_coords(pdb_lines):
    """Parse coordinates from PDB lines for a single model."""
    coords = []
    elements = []

    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_name = line[12:16].strip()
            if not is_heavy(atom_name):
                continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            coords.append([x, y, z])
            elements.append(get_element(atom_name))

    return np.array(coords) if coords else np.array([]), elements


def parse_pdb_models(pdb_file):
    """Parse all models from a multi-model PDB file."""
    models = {}
    current_model_num = None
    current_lines = []

    with open(pdb_file) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model_num = int(line.split()[1])
                current_lines = []
            elif line.startswith("ENDMDL"):
                if current_model_num is not None and current_lines:
                    coords, elements = parse_model_coords(current_lines)
                    models[current_model_num] = {
                        'coords': coords,
                        'elements': elements
                    }
                current_model_num = None
                current_lines = []
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                current_lines.append(line)

    # Handle single-model files without MODEL/ENDMDL
    if not models and current_lines:
        coords, elements = parse_model_coords(current_lines)
        models[1] = {'coords': coords, 'elements': elements}

    return models


def compute_spatial_similarity(m1, m2):
    """
    Compute spatial similarity between two models.

    Returns average minimum distance to atom of matching element.
    Lower values = more similar.
    """
    c1 = m1['coords']
    e1 = m1['elements']
    c2 = m2['coords']
    e2 = m2['elements']

    if len(c1) == 0 or len(c2) == 0:
        return 999.0

    dists = []

    # Direction A -> B
    for i, elem in enumerate(e1):
        matches = [j for j, e in enumerate(e2) if e == elem]
        if not matches:
            dists.append(5.0)  # Penalty for missing element
            continue

        diff = c2[matches] - c1[i]
        d = np.min(np.sqrt(np.sum(diff**2, axis=1)))
        dists.append(d)

    # Direction B -> A
    for i, elem in enumerate(e2):
        matches = [j for j, e in enumerate(e1) if e == elem]
        if not matches:
            dists.append(5.0)
            continue

        diff = c1[matches] - c2[i]
        d = np.min(np.sqrt(np.sum(diff**2, axis=1)))
        dists.append(d)

    return np.mean(dists)


def is_dominated(model_a, model_b):
    """
    Check if model_a is dominated by model_b.

    model_a is dominated if model_b is better (more negative) in BOTH
    delta_total AND delta_hbond.
    """
    return (model_b['delta_total'] < model_a['delta_total'] and
            model_b['delta_hbond'] < model_a['delta_hbond'])


def extract_model_from_pdb(pdb_path, model_num):
    """Extract a single model's PDB lines from a multi-model file."""
    lines = []
    in_model = False
    current_model_num = None

    with open(pdb_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model_num = int(line.split()[1])
                if current_model_num == model_num:
                    in_model = True
                    lines.append(line)
            elif line.startswith("ENDMDL"):
                if in_model:
                    lines.append(line)
                    return lines
                in_model = False
            elif in_model:
                lines.append(line)

    return lines


def write_pruned_pdb(output_path, models_to_keep, source_pdb):
    """Write a new PDB file with only the specified models."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, 'w') as out:
        for model_num in sorted(models_to_keep):
            model_lines = extract_model_from_pdb(source_pdb, model_num)
            for line in model_lines:
                out.write(line)


def compute_com_radius(model):
    """Compute Center of Mass and Maximum Radius for a model."""
    coords = model['coords']
    if len(coords) == 0:
        return np.zeros(3), 0.0
    
    com = np.mean(coords, axis=0)
    # Radius = max distance from CoM
    diff = coords - com
    dists = np.sqrt(np.sum(diff**2, axis=1))
    radius = np.max(dists)
    
    return com, radius


def prune_signature(pdb_path, model_results, similarity_threshold=1.0, verbose=False):
    """
    Prune models within a single signature (PDB file) using a strict redundancy filter.

    1. Filter out steric clashes (delta_fa_rep > 10).
    2. Sort remaining models by combined score: (delta_total + delta_hbond) ascending.
    3. Iterate from best to worst:
       - Keep the current model.
       - Remove ANY subsequent models that are structurally similar (redundancy score < threshold).
         (We no longer check for strict Pareto domination; we just pick the best combined scorer).
    
    Returns list of model numbers to keep.
    """
    # Step 1: Filter by fa_rep
    valid_models = []
    for model in model_results:
        if 'error' in model:
            continue
        if model.get('delta_fa_rep', 999) > 10.0:
            continue
        valid_models.append(model)

    if not valid_models:
        return []

    if len(valid_models) == 1:
        return [valid_models[0]['model_num']]

    # Step 2: Load structural data for valid models
    pdb_models = parse_pdb_models(pdb_path)

    # Create list of full model data
    model_data = []
    for m in valid_models:
        model_num = m['model_num']
        if model_num in pdb_models:
            structure = pdb_models[model_num]
            com, radius = compute_com_radius(structure)
            
            # Compute combined score for sorting
            combined_score = m.get('delta_total', 0) + m.get('delta_hbond', 0)
            
            model_data.append({
                'model_num': model_num,
                'energy': m,
                'structure': structure,
                'com': com,
                'radius': radius,
                'score': combined_score
            })

    if not model_data:
        return []

    # Step 3: Sort by combined score (best/lowest first)
    model_data.sort(key=lambda x: x['score'])

    # Step 4: Greedy elimination
    kept_models = []
    
    # We consume the list. model_data is sorted best to worst.
    while model_data:
        # Pick the best remaining model
        current = model_data.pop(0)
        kept_models.append(current)
        
        survivors = []
        c_com = current['com']
        c_rad = current['radius']
        
        for other in model_data:
            # Optimization: CoM Filter
            o_com = other['com']
            dist_sq = np.sum((c_com - o_com)**2)
            max_reach = c_rad + other['radius'] + similarity_threshold + 5.0 
            
            if dist_sq > max_reach**2:
                survivors.append(other)
                continue

            # Check structural similarity
            similarity = compute_spatial_similarity(current['structure'], other['structure'])
            
            if similarity < similarity_threshold:
                # Models are structurally similar (redundant).
                # Since 'current' has a better combined score (due to sorting), we prune 'other'.
                if verbose:
                    print(f"    Pruning model {other['model_num']} (redundant with {current['model_num']}, sim={similarity:.3f})")
                continue
            
            # If not similar, keep it.
            survivors.append(other)
            
        model_data = survivors

    # Return models to keep
    return [x['model_num'] for x in kept_models]


def main():
    parser = argparse.ArgumentParser(
        description="Prune candidate motifs based on steric clashes and structural redundancy"
    )
    parser.add_argument(
        '--results', '-r',
        default='merged_results.json',
        help='Path to energy scan results JSON'
    )
    parser.add_argument(
        '--input-dir', '-i',
        default='candidates_v2',
        help='Input candidates directory'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='candidates_v3',
        help='Output directory for pruned candidates'
    )
    parser.add_argument(
        '--max-fa-rep',
        type=float,
        default=10.0,
        help='Maximum delta_fa_rep threshold (default: 10.0)'
    )
    parser.add_argument(
        '--similarity-threshold',
        type=float,
        default=1.0,
        help='Structural similarity threshold for redundancy (default: 1.0)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print detailed pruning decisions'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Only report what would be pruned, do not write files'
    )

    args = parser.parse_args()

    # Load results
    print(f"Loading results from: {args.results}")
    with open(args.results) as f:
        data = json.load(f)

    results = data.get('results', {})
    print(f"Found {len(results)} candidate files in results")

    # Statistics
    stats = {
        'total_files': 0,
        'total_models_input': 0,
        'filtered_fa_rep': 0,
        'filtered_redundancy': 0,
        'total_models_output': 0,
        'files_with_survivors': 0,
        'files_empty': 0
    }

    # Group results by signature (PDB file)
    print(f"\nPruning candidates...")
    print(f"  Max fa_rep: {args.max_fa_rep}")
    print(f"  Similarity threshold: {args.similarity_threshold}")
    print()

    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)

    for pdb_path, entry in results.items():
        stats['total_files'] += 1

        models = entry.get('models', [])
        stats['total_models_input'] += len(models)

        # Count models filtered by fa_rep
        valid_models = [m for m in models if 'error' not in m and m.get('delta_fa_rep', 999) <= args.max_fa_rep]
        stats['filtered_fa_rep'] += len(models) - len(valid_models)

        if not os.path.exists(pdb_path):
            if args.verbose:
                print(f"  Skipping {pdb_path} (file not found)")
            continue
            
        if args.verbose:
            print(f"  Processing {os.path.basename(pdb_path)} ({len(models)} models)...")

        # Prune this signature
        models_to_keep = prune_signature(
            pdb_path,
            models,
            similarity_threshold=args.similarity_threshold,
            verbose=args.verbose
        )

        # Count redundancy removals
        stats['filtered_redundancy'] += len(valid_models) - len(models_to_keep)
        stats['total_models_output'] += len(models_to_keep)

        if models_to_keep:
            stats['files_with_survivors'] += 1

            # Construct output path
            rel_path = os.path.relpath(pdb_path, args.input_dir)
            output_path = os.path.join(args.output_dir, rel_path)

            if not args.dry_run:
                write_pruned_pdb(output_path, models_to_keep, pdb_path)

            if args.verbose:
                print(f"  {os.path.basename(pdb_path)}: {len(models)} -> {len(models_to_keep)} models")
        else:
            stats['files_empty'] += 1
            if args.verbose:
                print(f"  {os.path.basename(pdb_path)}: {len(models)} -> 0 models (all pruned)")

    # Print summary
    print("\n" + "="*60)
    print("PRUNING SUMMARY")
    print("="*60)
    print(f"Input files:              {stats['total_files']}")
    print(f"Input models:             {stats['total_models_input']}")
    print(f"Filtered (fa_rep > {args.max_fa_rep}): {stats['filtered_fa_rep']}")
    print(f"Filtered (redundancy):    {stats['filtered_redundancy']}")
    print(f"Output models:            {stats['total_models_output']}")
    print(f"Output files (non-empty): {stats['files_with_survivors']}")
    print(f"Files fully pruned:       {stats['files_empty']}")
    print()

    retention_rate = (stats['total_models_output'] / stats['total_models_input'] * 100
                      if stats['total_models_input'] > 0 else 0)
    print(f"Retention rate: {retention_rate:.1f}%")

    if args.dry_run:
        print("\n[DRY RUN - no files written]")
    else:
        print(f"\nPruned candidates written to: {args.output_dir}")


if __name__ == '__main__':
    main()
