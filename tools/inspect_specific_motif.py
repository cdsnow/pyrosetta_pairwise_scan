import json
import sys

def inspect_arg_motifs(results_file, target_pdb):
    with open(results_file, 'r') as f:
        data = json.load(f)
    
    results = data.get('results', {})
    
    # Find the entry
    entry = results.get(target_pdb)
    if not entry:
        print(f"Entry {target_pdb} not found in results.")
        # Try to find partial match
        for k in results.keys():
            if target_pdb in k:
                print(f"Found similar key: {k}")
                entry = results[k]
                break
    
    if not entry:
        return

    print(f"Inspecting {target_pdb}")
    print(f"Total models: {len(entry.get('models', []))}")
    
    print(f"{'Model':<5} | {'Total':<10} | {'HBond':<10} | {'fa_rep':<10} | {'Status'}")
    print("-" * 60)
    
    for m in entry.get('models', []):
        if 'error' in m:
            print(f"{m['model_num']:<5} | ERROR: {m['error']}")
            continue
            
        rep = m.get('delta_fa_rep', 999)
        status = "CLASH" if rep > 10 else "OK"
        
        print(f"{m['model_num']:<5} | {m.get('delta_total', 0):<10.2f} | {m.get('delta_hbond', 0):<10.2f} | {rep:<10.2f} | {status}")

if __name__ == "__main__":
    inspect_arg_motifs('merged_results.json', 'candidates_v2/F6_ASP_SC/ARG_N_NE_NH2.pdb')
