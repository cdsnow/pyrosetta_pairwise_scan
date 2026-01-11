import numpy as np
import sys
import os
import itertools

def is_heavy(name):
    n = name.strip()
    if n.startswith('H'): return False
    if n[0].isdigit() and len(n) > 1 and n[1] == 'H': return False
    return True

def get_element(name):
    n = name.strip()
    if n[0].isdigit(): n = n[1:]
    return n[0]

def parse_pdb_models(pdb_file):
    models = []
    current_model = None
    
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = {'id': int(line.split()[1]), 'atoms': [], 'coords': [], 'elements': []}
            elif line.startswith("ATOM") and current_model is not None:
                atom_name = line[12:16].strip()
                if not is_heavy(atom_name): continue
                
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
                    
    # Handle implicit single model if no ENDMDL/MODEL tags (though grep showed 13 models)
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
        matches = [j for j, e in enumerate(e2) if e == elem]
        if not matches: 
            dists.append(5.0) 
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

def compute_rmsd_no_superposition(m1, m2):
    # Assumes identical atom ordering and count
    c1 = m1['coords']
    c2 = m2['coords']
    
    if c1.shape != c2.shape:
        return -1.0 # Mismatch
        
    diff = c1 - c2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

def main():
    pdb_file = "candidates_v4/F6_ASP_SC/ARG_NE_NH1.pdb"
    if not os.path.exists(pdb_file):
        print(f"File {pdb_file} not found.")
        sys.exit(1)
        
    models = parse_pdb_models(pdb_file)
    print(f"Loaded {len(models)} models from {pdb_file}")
    
    pairs = list(itertools.combinations(models, 2))
    print(f"Analyzing {len(pairs)} pairs...")
    
    sim_scores = []
    rmsd_scores = []
    
    print(f"{ 'Model A':<8} {'Model B':<8} {'SimScore':<10} {'RMSD':<10}")
    print("-" * 40)
    
    for m1, m2 in pairs:
        sim = compute_spatial_similarity(m1, m2)
        rmsd = compute_rmsd_no_superposition(m1, m2)
        
        sim_scores.append(sim)
        if rmsd >= 0:
            rmsd_scores.append(rmsd)
            
        print(f"{m1['id']:<8} {m2['id']:<8} {sim:.4f}     {rmsd:.4f}")

    if sim_scores:
        print("\nSummary Statistics:")
        print(f"Motif Redundancy Score: Min={min(sim_scores):.4f}, Max={max(sim_scores):.4f}, Mean={np.mean(sim_scores):.4f}")
    
    if rmsd_scores:
        print(f"No-Superposition RMSD:  Min={min(rmsd_scores):.4f}, Max={max(rmsd_scores):.4f}, Mean={np.mean(rmsd_scores):.4f}")

if __name__ == "__main__":
    main()
