import os
import glob
import numpy as np
import collections
import argparse
import itertools

def is_heavy(name):
    n = name.strip()
    if n.startswith('H'): return False
    if n[0].isdigit() and len(n) > 1 and n[1] == 'H': return False
    return True

def get_element(name):
    # Heuristic element guess
    n = name.strip()
    if n[0].isdigit(): n = n[1:]
    return n[0] # C, N, O, S, etc.

def parse_pdb_models(pdb_file):
    models = []
    current_model = None
    
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = {'id': '?', 'atoms': [], 'coords': [], 'elements': []}
            elif line.startswith("REMARK") and "ID:" in line:
                if current_model:
                    parts = line.split("ID:")
                    if len(parts) > 1: current_model['id'] = parts[1].strip()
            # Handle implicit model 1
            elif line.startswith("ATOM") and current_model is None:
                current_model = {'id': '?', 'atoms': [], 'coords': [], 'elements': []}
                
            if line.startswith("ATOM") and current_model:
                atom_name = line[12:16].strip()
                if not is_heavy(atom_name): continue
                
                # Parse Coords
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                current_model['atoms'].append(atom_name)
                current_model['coords'].append([x, y, z])
                current_model['elements'].append(get_element(atom_name))
                
            elif line.startswith("ENDMDL"):
                if current_model:
                    current_model['coords'] = np.array(current_model['coords'])
                    models.append(current_model)
                    current_model = None
                    
    if current_model:
        current_model['coords'] = np.array(current_model['coords'])
        models.append(current_model)
        
    return models

def compute_spatial_similarity(m1, m2):
    # Metric: Average minimum distance to atom of matching element
    
    c1 = m1['coords']
    e1 = m1['elements']
    
    c2 = m2['coords']
    e2 = m2['elements']
    
    if len(c1) == 0 or len(c2) == 0: return 999.0
    
    dists = []
    
    # Direction A -> B
    for i, elem in enumerate(e1):
        # Find matching elements in B
        matches = [j for j, e in enumerate(e2) if e == elem]
        if not matches:
            # Penalty for missing element match? 
            # Or ignore?
            # User said "computing the average minimum distance... to an atom of matching element".
            # If no matching element exists, the distance is undefined.
            # Let's assign a penalty distance (e.g. 5.0 A) or max observed?
            dists.append(5.0) 
            continue
            
        # Dist to all matches
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", default="candidates", help="Directory to scan")
    args = parser.parse_args()
    
    pdb_files = glob.glob(os.path.join(args.dir, "**/*.pdb"), recursive=True)
    pdb_files.sort()
    
    print("--- Measuring Spatial Redundancy (Avg Min Dist) ---")
    
    scores = []
    
    for p in pdb_files:
        if "load_" in p: continue
        sig = os.path.basename(p).replace(".pdb", "")
        
        models = parse_pdb_models(p)
        if len(models) < 2: continue
        
        # Calculate pairwise scores for first few models to gauge spread
        # Compare Model 1 vs Rest
        ref = models[0]
        
        model_scores = []
        for m in models[1:]:
            score = compute_spatial_similarity(ref, m)
            model_scores.append(score)
            
        avg_score = np.mean(model_scores)
        max_score = np.max(model_scores)
        
        scores.append((p, avg_score, max_score, len(models)))
        
        # If score is low (< 0.5?), it suggests high redundancy despite potentially different atom counts
        if avg_score > 2.0: # Arbitrary high threshold for "Interesting/Different"
             pass 
             
    # Top redundant scores (Lowest)
    print("\n--- Top High-Redundancy Groups (Lowest Metric Score) ---")
    scores.sort(key=lambda x: x[1]) # Sort by Avg Score Ascending
    for p, avg, mx, n in scores[:10]:
        print(f"{p} (N={n}) | Metric: {avg:.4f} | Max: {mx:.4f}")
        
    print("\n--- Top Low-Redundancy Groups (Highest Metric Score - Distinct) ---")
    # Sort Descending
    scores.sort(key=lambda x: x[1], reverse=True)
    for p, avg, mx, n in scores[:10]:
        print(f"{p} (N={n}) | Metric: {avg:.4f} | Max: {mx:.4f}")

if __name__ == "__main__":
    main()
