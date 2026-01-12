import sys
import numpy as np
from scipy.spatial.distance import cdist
import argparse

def parse_pdb_coords_with_names(lines):
    atoms = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                name = line[12:16].strip()
                res = line[17:20].strip()
                resseq = line[22:26].strip()
                chain = line[21].strip()
                element = line[76:78].strip()
                
                # Filter hydrogens (using same logic as validate_sterics)
                known_heavy = {'C', 'N', 'O', 'S', 'P', 'F', 'CL', 'BR', 'I', 'FE', 'MG', 'CA', 'ZN', 'MN', 'NA', 'K', 'LI', 'AG', 'AU', 'CU', 'CO', 'NI'}
                
                is_hydrogen = False
                
                if element == 'H' or element == 'D':
                    is_hydrogen = True
                elif element in known_heavy:
                    is_hydrogen = False
                else:
                    # Element is missing or trash (e.g. '1'). Fallback to name.
                    if name.startswith('H'):
                         if name != 'HG': 
                             is_hydrogen = True
                    elif len(name) > 1 and name[0].isdigit() and name[1] == 'H':
                        is_hydrogen = True
                
                if is_hydrogen:
                    continue

                atoms.append({
                    'coord': [x, y, z],
                    'name': name,
                    'res': res,
                    'resseq': resseq,
                    'chain': chain,
                    'line': line.strip()
                })
            except ValueError:
                continue
    return atoms

def get_model_lines(pdb_path, model_num):
    model_lines = []
    current_model_num = None
    in_model = False
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
                    return model_lines
                in_model = False
            elif in_model:
                model_lines.append(line)
            elif not has_models and model_num == 1 and (line.startswith('ATOM') or line.startswith('HETATM')):
                model_lines.append(line)
                
    if not has_models and model_num == 1:
        return model_lines
    return model_lines

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_path')
    parser.add_argument('model_num', type=int)
    args = parser.parse_args()

    target_path = 'flagorigin_and_loop.pdb'
    candidate_path = args.pdb_path
    model_num = args.model_num

    print(f"Loading target: {target_path}")
    with open(target_path, 'r') as f:
        target_atoms = parse_pdb_coords_with_names(f.readlines())

    print(f"Loading candidate: {candidate_path} model {model_num}")
    model_lines = get_model_lines(candidate_path, model_num)
    candidate_atoms = parse_pdb_coords_with_names(model_lines)

    if not candidate_atoms:
        print("Error: No atoms found for candidate model (or all were filtered)")
        sys.exit(1)

    t_coords = np.array([a['coord'] for a in target_atoms])
    c_coords = np.array([a['coord'] for a in candidate_atoms])

    dists = cdist(t_coords, c_coords)
    min_dist = np.min(dists)
    print(f"Minimum Distance: {min_dist:.4f}")

    # Find indices of closest pair
    idx = np.unravel_index(np.argmin(dists), dists.shape)
    t_idx, c_idx = idx

    t_atom = target_atoms[t_idx]
    c_atom = candidate_atoms[c_idx]

    print("\nClosest Pair:")
    print(f"Target Atom: {t_atom['line']}")
    print(f"Cand   Atom: {c_atom['line']}")
    print(f"Distance: {dists[t_idx, c_idx]:.4f}")
